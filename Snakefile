import os
import re
import subprocess
import json

with open("config.json", "r") as f:
  config = json.load(f)

GWASDIR = config['gwasdir']
OUTDIR = config['outdir']
OUTDIR2 = "/user/home/ny19205/scratch/mreve_testing/" # tmp local storage for data/ folder
LDREFHOST = config['ldrefhost']
LDREFPATH = OUTDIR + "/reference/" + config['ldrefname']
RFHOST = config['rfhost']

os.makedirs(OUTDIR + '/job_reports', exist_ok=True)
os.makedirs(OUTDIR + '/reference', exist_ok=True)
os.makedirs(OUTDIR + '/resources', exist_ok=True)
os.makedirs(OUTDIR + '/neo4j', exist_ok=True)

IDLIST = OUTDIR + "/resources/ids_all.txt" # this is list of ids of datasets in GWASDIR # this makes sure data for them exists at the location 
INSTRUMENTLIST = OUTDIR + "/resources/instruments.txt"

IDLIST_test = OUTDIR + "/resources/ids.txt"
#IDLIST_test2 = OUTDIR + "/resources/ids_double_test_10.txt"
ID_test = open(IDLIST_test, 'r').read().strip().split('\n') # list of 5 IDs for which I manually added/produced data -- testing rule check_extract_master


NTHREAD = 5
CHUNKS = list(range(1,2)) 
nCHUNKS = len(CHUNKS)


rule install_r_packages: 
    output:
        expand('{OUTDIR}/r_packages_installed.csv', OUTDIR=OUTDIR)
    shell:
        "Rscript scripts/install_r_packages.R"


rule get_genes: # YES
    output:
        expand('{OUTDIR}/resources/genes.rdata', OUTDIR=OUTDIR)
    shell:
        "Rscript scripts/get_genes.R {output}"


rule get_id_info: # YES
    input:
        expand('{IDLIST}', IDLIST=IDLIST)
    output:
        expand('{IDLIST}.rdata', IDLIST=IDLIST)
    shell:
        "Rscript scripts/get_ids.R {input} {GWASDIR} {output}"

rule download_ldref: # YES
    output:
        expand("{LDREFPATH}.bed", LDREFPATH=LDREFPATH),
        expand("{LDREFPATH}.bim", LDREFPATH=LDREFPATH),
        expand("{LDREFPATH}.fam", LDREFPATH=LDREFPATH)
    shell:
        "curl -s {LDREFHOST} | tar xzvf - -C {OUTDIR}/reference"

        
rule create_ldref_sqlite: #YES # new version of the rule # had to run in screen
    input:
        expand("{LDREFPATH}.bed", LDREFPATH=LDREFPATH),
        expand("{LDREFPATH}.bim", LDREFPATH=LDREFPATH),
        expand("{LDREFPATH}.fam", LDREFPATH=LDREFPATH)
    output:
        expand("{LDREFPATH}.sqlite", LDREFPATH=LDREFPATH)
    shell:
        "Rscript -e 'gwasvcf::set_plink()' -e 'gwasvcf::create_ldref_sqlite(\"{LDREFPATH}\", \"{output}\")' "
    

rule download_rfobj: # YES
    output:
        expand("{OUTDIR}/reference/rf.rdata", OUTDIR=OUTDIR)
    shell:
        "wget -O {output} {RFHOST}"


#rule instrument_list: # THIS WORKS but not on BP: don't have access to GWASDIR - so I'm manually adding INSTRUMENTLIST output (created on BC4) to dir and skipping this step
#    input:
# HERE ALSO SHOULD BE CHECKING FOR clump.txt?
#        expand('{IDLIST}', IDLIST=IDLIST)
#    output:
#        expand('{INSTRUMENTLIST}', INSTRUMENTLIST=INSTRUMENTLIST)
#    shell:
#        './scripts/instrument_list.py --dirs {GWASDIR} --idlists {IDLIST} --output {output}'


# this currenntly does not work on BP/BC4 becaus of bcftools problems
# can't run like this anyway - need to add threads to process each id
# copied the output produced locally to output filders
# also curretly gwas files and clup files are in the same place as output - this has to be changed when gwas folder issues are solved

#rule extract_master:
#    input:
#        expand('{INSTRUMENTLIST}', INSTRUMENTLIST=INSTRUMENTLIST),
#        expand('{LDREFPATH}.sqlite', LDREFPATH=LDREFPATH)
#    output:
#        '{OUTDIR}/data/{id}/ml.csv.gz'
#    shell:
#        """
#  mkdir -p {OUTDIR}/data/{wildcards.id}
#  Rscript scripts/extract_masterlist.r \
#  --snplist {INSTRUMENTLIST} \
#  --gwasdir {GWASDIR} \
#  --out {output} \
#  --dbfile {LDREFPATH}.sqlite \
#  --id {wildcards.id} \
#  --instrument-list \
#  --get-proxies yes
#        """


rule check_extract_master:
    input: 
        expand('data/{id}/ml.csv.gz', OUTDIR=OUTDIR, id=ID_test), # later add {OUTDIR}
    output:
        expand('{OUTDIR}/resources/flag_extract_master', OUTDIR=OUTDIR)
    shell:
        "touch {output}"

rule mr:
    input:
        flag = expand('{OUTDIR}/resources/flag_extract_master', OUTDIR=OUTDIR),
        rf = expand('{OUTDIR}/reference/rf.rdata', OUTDIR=OUTDIR),
        idlist = expand('{IDLIST}.rdata', IDLIST=IDLIST),
    output:
        'data/{id}/mr.rdata' # put {OUTPUT} back + and in outdir in R command too
    shell: 
        """
    Rscript scripts/mr.R \
    --idlist {input.idlist} \
    --outdir {OUTDIR2} \
    --id {wildcards.id} \
    --rf {input.rf} \
    --what eve \
    --threads {NTHREAD}
        """

rule check_mr:
    input:
    	expand('data/{id}/mr.rdata', id=ID_test), # later add {OUTDIR}
    output:
    	expand('{OUTDIR}/resources/flag_mr', OUTDIR=OUTDIR)
    shell:
    	"touch {output}"


rule heterogeneity:
    input:
        flag = expand('{OUTDIR}/resources/flag_mr', OUTDIR=OUTDIR),
        mr = 'data/{id}/mr.rdata', # put {OUTPUT} back
        idlist = expand('{IDLIST}.rdata', IDLIST=IDLIST)
    output:
        'data/{id}/heterogeneity.rdata' # put {OUTPUT} back + and in outdir in R command too
    shell:
        """
    Rscript scripts/heterogeneity.R \
    --idlist {input.idlist} \
    --outdir {OUTDIR2} \
    --id {wildcards.id} \
    --what eve \
    --threads {NTHREAD}
        """


rule check_heterogeneity:
    input:
        expand('data/{id}/heterogeneity.rdata', id=ID_test), # later add {OUTDIR}
    output:
    	expand('{OUTDIR}/resources/flag_heterogeneity', OUTDIR=OUTDIR)
    shell:
    	"touch {output}"


rule neo4j_mr:
    input:
        instruments = expand('{INSTRUMENTLIST}', INSTRUMENTLIST=INSTRUMENTLIST),
        idlist = expand('{IDLIST}.rdata', IDLIST=IDLIST),
        mr = 'data/{id}/mr.rdata',
        het = 'data/{id}/heterogeneity.rdata'
    output:
        'data/{id}/neo4j_stage/{id}_mr.csv.gz',
        'data/{id}/neo4j_stage/{id}_moe.csv.gz',
        'data/{id}/neo4j_stage/{id}_int.csv.gz',
        'data/{id}/neo4j_stage/{id}_het.csv.gz',
        'data/{id}/neo4j_stage/{id}_met.csv.gz',
        'data/{id}/neo4j_stage/{id}_vt.csv.gz',
        'data/{id}/neo4j_stage/{id}_inst.csv.gz'
    shell:
        """
        Rscript scripts/prepare_neo4j_mr.R \
            --outdir {OUTDIR2} \
            --snplist {input.instruments} \
            --idlist {input.idlist} \
            --id {wildcards.id} \
            --headerid {wildcards.id} 
        """ # need to do header id part differently


rule check_neo4j_mr:
    input:
          expand('data/{id}/heterogeneity.rdata', id=ID_test), # later add {OUTDIR}
          expand('data/{id}/neo4j_stage/{id}_mr.csv.gz', id=ID_test), # later add {OUTDIR},
          expand('data/{id}/neo4j_stage/{id}_moe.csv.gz', id=ID_test), # later add {OUTDIR},
          expand('data/{id}/neo4j_stage/{id}_int.csv.gz', id=ID_test), # later add {OUTDIR},
          expand('data/{id}/neo4j_stage/{id}_het.csv.gz', id=ID_test), # later add {OUTDIR},
          expand('data/{id}/neo4j_stage/{id}_met.csv.gz', id=ID_test), # later add {OUTDIR},
          expand('data/{id}/neo4j_stage/{id}_vt.csv.gz', id=ID_test), # later add {OUTDIR},
          expand('data/{id}/neo4j_stage/{id}_inst.csv.gz', id=ID_test), # later add {OUTDIR}'
          expand('{OUTDIR}/resources/neo4j_stage/header_mr.csv.gz', OUTDIR=OUTDIR),
          expand('{OUTDIR}/resources/neo4j_stage/header_moe.csv.gz', OUTDIR=OUTDIR),
          expand('{OUTDIR}/resources/neo4j_stage/header_int.csv.gz', OUTDIR=OUTDIR),
          expand('{OUTDIR}/resources/neo4j_stage/header_het.csv.gz', OUTDIR=OUTDIR),
          expand('{OUTDIR}/resources/neo4j_stage/header_met.csv.gz', OUTDIR=OUTDIR),
          expand('{OUTDIR}/resources/neo4j_stage/header_vt.csv.gz', OUTDIR=OUTDIR),
          expand('{OUTDIR}/resources/neo4j_stage/header_inst.csv.gz', OUTDIR=OUTDIR)

    output:
        expand('{OUTDIR}/resources/flag_neo4j_mr', OUTDIR=OUTDIR)
    shell:
        "touch {output}"


# in this rule we don't use data/ folder, do OUTDIR is ok
rule neo4j_others:
    input:
      genes = expand('{OUTDIR}/resources/genes.rdata', OUTDIR=OUTDIR),
      instruments = expand('{INSTRUMENTLIST}', INSTRUMENTLIST=INSTRUMENTLIST),
      idlist = expand('{IDLIST}.rdata', IDLIST=IDLIST)
    output:
        expand('{OUTDIR}/resources/neo4j_stage/genes.csv.gz', OUTDIR=OUTDIR),
        expand('{OUTDIR}/resources/neo4j_stage/gv.csv.gz', OUTDIR=OUTDIR),
        expand('{OUTDIR}/resources/neo4j_stage/traits.csv.gz', OUTDIR=OUTDIR),
        expand('{OUTDIR}/resources/neo4j_stage/variants.csv.gz', OUTDIR=OUTDIR)
    shell:
        """Rscript scripts/prepare_neo4j_other.R \
            --outdir {OUTDIR} \
            --snplist {input.instruments} \
            --idlist {input.idlist} \
            --genes {input.genes}
        """


rule collect_neo4j_mr:
    input:
        flag = expand('{OUTDIR}/resources/flag_neo4j_mr',OUTDIR=OUTDIR),
        idlist = expand('{IDLIST}', IDLIST=IDLIST_test) # 5 ids
    output:
        mr =  '{OUTDIR}/resources/neo4j_stage/{chunk}_mr.csv.gz', 
        moe = '{OUTDIR}/resources/neo4j_stage/{chunk}_moe.csv.gz',
        int = '{OUTDIR}/resources/neo4j_stage/{chunk}_int.csv.gz',
        het = '{OUTDIR}/resources/neo4j_stage/{chunk}_het.csv.gz',
        met = '{OUTDIR}/resources/neo4j_stage/{chunk}_met.csv.gz',
        vt =  '{OUTDIR}/resources/neo4j_stage/{chunk}_vt.csv.gz', 
        inst ='{OUTDIR}/resources/neo4j_stage/{chunk}_inst.csv.gz'
    wildcard_constraints:
        chunk="\d+" # only files that contain digits in the name - otherwise picks up 'header' files
    shell:
        "python scripts/collect_neo4j_mr.py --outdir {OUTDIR} --idlist {input.idlist}  --chunks_all {nCHUNKS} --chunk_current {wildcards.chunk}"



rule check_collect_neo4j_mr:
    input:
      expand('{OUTDIR}/resources/neo4j_stage/{chunk}_mr.csv.gz', OUTDIR=OUTDIR, chunk=CHUNKS),
        expand('{OUTDIR}/resources/neo4j_stage/{chunk}_moe.csv.gz', OUTDIR=OUTDIR, chunk=CHUNKS),
        expand('{OUTDIR}/resources/neo4j_stage/{chunk}_int.csv.gz', OUTDIR=OUTDIR, chunk=CHUNKS),
        expand('{OUTDIR}/resources/neo4j_stage/{chunk}_het.csv.gz', OUTDIR=OUTDIR, chunk=CHUNKS),
        expand('{OUTDIR}/resources/neo4j_stage/{chunk}_met.csv.gz', OUTDIR=OUTDIR, chunk=CHUNKS),
        expand('{OUTDIR}/resources/neo4j_stage/{chunk}_vt.csv.gz', OUTDIR=OUTDIR, chunk=CHUNKS),
        expand('{OUTDIR}/resources/neo4j_stage/{chunk}_inst.csv.gz', OUTDIR=OUTDIR, chunk=CHUNKS),
    output:
        expand('{OUTDIR}/resources/flag_collect_neo4j_mr', OUTDIR=OUTDIR)
    shell:
        "touch {output}"


