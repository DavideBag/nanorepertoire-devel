#!/usr/bin/env python
#This is a python script to read the CD-HIT output to produce summary info on clusters
#Give the name of the cluster representative, the cluster size, and the cluster average identity to representative

import sys

filein = open(sys.argv[1])
fileout = open(sys.argv[2], "w")

clusterrep = None
clustermembers = 0
clusteridentity = 1
clusteridentitysum = 0

fileout.write("Representative,Count,Identity\n")
for line in filein:
    if ">Cluster" in line:
        if clusterrep == None:
            continue
        if clustermembers == 1:
            clusteridentity = 1
        else:
            clusteridentity = float(clusteridentitysum) / (clustermembers - 1)
        fileout.write(clusterrep+","+str(clustermembers)+","+str(clusteridentity)+"\n")
        clustermembers = 0
        clusteridentitysum = 0
    else:
        clustermembers += 1
        #print(line)
        readname = line.rstrip().split(">")[1].split("...")[0]
        if "*" in line:
            clusterrep = readname
        else:
            #print(line)
            thisidentity = float(line.rstrip().split("at ")[1].split("%")[0])
            clusteridentitysum += thisidentity
if clustermembers == 1:
    clusteridentity = 1
else:
    clusteridentity = float(clusteridentitysum) / (clustermembers - 1)
fileout.write(clusterrep+","+str(clustermembers)+","+str(clusteridentity)+"\n")


filein.close()
fileout.close()
