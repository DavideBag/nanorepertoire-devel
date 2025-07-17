#!/usr/bin/env python
#This is a script to go through a CD-HIT clusters output file and for each entry like
#@R7
#MAQVRLVESGGGLAQAGGSLRLSCEASGFTSDDWAIGWFRQAPGKEREGVSCIRHSTQTTAYADSVKGRFNISGDSAKNTVYLQMNSLKPKDTAVYYCAALILFMGTYYDPIDLLGYEYENWGQGIQVTVSS
#Extract the CDR3 sequence and write it out in a file in fasta format
#The CDR3 begins with amino acids which commence after a YYC and terminates with the amino acids which precede a WGQ
#Skip any CDR3 sequence that has already been encountered, and write out total unique CDR3s at the end
#Also make a histogram of lengths of these CDR3s
# UPDATE - the script also outputs a tsv file with all CDR3 even if not unique: this is to be used later

import sys
import re
import Bio
from Bio.Seq import Seq
## from Bio.Alphabet import generic_dna ##not used anymore in biopython
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="cdhit clusters file",
                    action="store", dest = "input")
parser.add_argument("-c", "--cdr", help="output unique CDR3 fasta file",
                    action="store", dest = "cdr")
parser.add_argument("-t", "--tsv", help="output all CDR3 TSV file",
                    action="store", dest = "tsv")
parser.add_argument("-o", "--hist", help="output histogram file",
                    action="store", dest = "hist")
## additional options to enrich the fasta name
parser.add_argument("-s", "--sample", help="sample identifier",
                    action="store", dest = "sample")
parser.add_argument("-m", "--immuno", help="immunisation of the sample",
                    action="store", dest = "immuno")
parser.add_argument("-b", "--boost", help="boost round of the immunisation",
                    action="store", dest = "boost")
args = parser.parse_args()

bigset = set()

preseq = "YYC"
postseq = "WGQ"
# my_rna.find(re.findall(r"A.{2}G", str(my_rna))[0])

filein = open(args.input)
fileout = open(args.cdr, "w")
tsvfile = open(args.tsv, "w")
histo = open(args.hist, "w")

tsvfile.write("ID\tCDR3\tsequence\tunique\n")

linecount = 0
for line in filein:
  linecount += 1
  if linecount % 2 == 1:
    header = line.rstrip()
    fastaheader = header.replace("@",">") #+ "_S" + args.sample + "_I" + args.immuno + "_B" + args.boost
    identifier = header.replace("@", "")
  else:
    sequence = Seq(line.rstrip())
    ## modified pattern search using a motif, i.e. introduced a regex
    ## meaning the regex itself might not yield any results, so additional
    ## conditionals need to be introduced before using .find()
    ## NB: always finds using first regex match
    ### startaa = sequence.find(preseq) + 3
    if len(re.findall(r"T.{2}Y.{1}C", str(sequence))) > 0:
      startaa = sequence.find(re.findall(r"T.{2}Y.{1}C", str(sequence))[0]) + 3
    else:
      startaa = -1
    ### endaa = sequence.find(postseq)
    if len(re.findall(r".{6}TVSS", str(sequence))) > 0:
      endaa = sequence.find(re.findall(r".{6}TVSS", str(sequence))[0])
    else:
      endaa = -1
    if startaa == -1 or endaa == -1:
      tsvfile.write(identifier + "\t" + "NA" + "\t" + line.rstrip() + "\tno-cdr3\n")
      continue
    targetbit = line.rstrip()[startaa:endaa]
    if targetbit in bigset or len(targetbit) > 50 or len(targetbit) < 1:
      tsvfile.write(identifier + "\t" + targetbit + "\t" + line.rstrip() + "\tnon-unique\n")
      continue
    bigset.add(targetbit)
    fileout.write(fastaheader + "\n")
    fileout.write(targetbit + "\n")
    tsvfile.write(identifier + "\t" + targetbit + "\t" + line.rstrip() + "\tunique\n")

filein.close()
fileout.close()

#print("Got total number of CDR3s: "+str(len(bigset)))

bighist = {}
for i in range(51):
  bighist[i] = 0
for sequence in bigset:
  bighist[len(sequence)] += 1
histo.write("Size,Count")
for i in range(51):
  histo.write(str(i) + "," + str(bighist[i]) + "\n")
