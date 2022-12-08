import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--chunks_all', default = 5, required=False)
parser.add_argument('--chunk_current', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--idlist', required=True)
args = parser.parse_args()

idlist = args.idlist
outdir = args.outdir
chunk_current = int(args.chunk_current)
chunks_all= int(args.chunks_all)


assert(chunk_current <= chunks_all)


#local
#chunks_all=list(range(1,6))
#chunk_current = 1
#idlist =  '/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project5/mr-eve/resources/ids.txt'
#outdir = '/Users/ny19205/OneDrive - University of Bristol/Documents - OneDrive/Mini-project5/mr-eve/'

# read id list  
ids = open(idlist, 'r').read().strip().split('\n')

# split list of ids into N of CHUNKS subarrays : 10 ids / 5 chunks = 5 subarrays by 2 ids 
ids_list = np.array_split(np.array(ids), chunks_all) 
print("all ids: {}".format(ids_list))
l = ids_list[chunk_current-1]
print("current chunk: {}".format(l))

def inputfile(what, ids, outdir):
  # construct input file name
  # NB add outdir to first pos! + "/data"
  return ["data/" + x + "/neo4j_stage/" + x + "_" + what + ".csv.gz" for x in ids]

def outputfile(what, outdir):
   # construct output file name
  return outdir + "/resources/neo4j_stage/" + str(chunk_current) + "_" + what +".csv.gz"
  
def cmd(inputs, output):
  # concatenate data in this defined chunk, within thesame type (i.e. 'what')
  return "cat {} | gunzip -c | awk NF | gzip -c > {}".format(inputs, output)


os.system(cmd(" ".join(inputfile("mr", l, outdir)),  outputfile("mr", outdir)))
os.system(cmd(" ".join(inputfile("moe", l, outdir)), outputfile("moe", outdir)))
os.system(cmd(" ".join(inputfile("int", l, outdir)), outputfile("int", outdir)))
os.system(cmd(" ".join(inputfile("het", l, outdir)), outputfile("het", outdir)))
os.system(cmd(" ".join(inputfile("met", l, outdir)), outputfile("met", outdir)))
os.system(cmd(" ".join(inputfile("vt", l, outdir)),  outputfile("vt", outdir)))
os.system(cmd(" ".join(inputfile("inst", l, outdir)),outputfile("inst", outdir)))

