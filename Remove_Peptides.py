#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 12:33:23 2021
@author: joe
"""


import numpy as np
import pandas as pd
import math
from scipy.spatial import KDTree
import os

#USE NVT OR MD FILES:
'''
nvtin = "nvt.gro"
nvtout = "nvt_pr.gro"
intpr = "nvt.tpr"
xtcin = "nvt.xtc"
xtcout = "nvt_pr.xtc"
'''

'''
nvtin = "nvt.gro"
nvtout = "nvt_pr.gro"
intpr = "md.tpr"
xtcin = "md.xtc"
xtcout = "md_pr.xtc"
'''


nvtin = "nvt.gro"

## First open the file and strip the lines to list "lines".

lines = []

with open (nvtin) as f:
    for num, line in enumerate(f, 1):
        line = [line.rstrip('\n')]
        lines.append(line)

## Remove the header and the end.

lines.pop(0);lines.pop(0);lines.pop(-1)

## For each line in lines split it to th seperate parts.

l = []
for n,v in enumerate(lines):
    indices = [0,5,8,12,15,20,22,28,30,36,38,44]
    p = (v[0])
    parts = [p[i:j] for i,j in zip(indices, indices[1:]+[None])]
    #print(parts)
    value = []
    value.append(int(n+1))
    value.append(str(parts[1]))
    value.append(str(parts[3]))
    value.append(float(parts[6]))
    value.append(float(parts[8]))
    value.append(float(parts[10]))
    #print(value)
    l.append(value)

## Creat the dataframe from the split.

df = pd.DataFrame(l)
df.columns = ['Index', 'MolName', 'AtomName', 'x', 'y', 'z']
print(df)
## Make a new DataFrame for each AtomType, add all to Dict.

Atom_Dict = dict(tuple(df.groupby("AtomName")))

#print(Atom_Dict)

## Make new DF called SurfaceAtoms of just surface gold.

AU = Atom_Dict["AUB"]
Pythag = []
SA = []


for index, row in AU.iterrows():
    x = float(row['x'])
    y = float(row['y'])
    z = float(row['z'])
    pyth = (math.sqrt(x**2+y**2+z**2))
    if pyth > 6.5:
        SA.append(row)
        Pythag.append(pyth)
    else:
        continue


SurfaceAtoms = pd.DataFrame(SA)
SurfaceAtoms["Pythag"] = Pythag
#print(SurfaceAtoms)


SAxyz = SurfaceAtoms[["x","y","z"]].to_numpy()

#print(SAxyz)

k = KDTree(SAxyz)

pts = np.array([[7, 7, 7]])
print(k.query(pts[0]))

kout = []
SG = Atom_Dict[" SG"]

with open ("DISTANCES.txt", "w") as f:
    f.write("ATOM, INDEX, DISTANCE, ATOM, INDEX \n")

Diss_S_Index = []


for index, row in SG.iterrows():
    LINE = []
    SGx = float(row['x'])
    SGy = float(row['y'])
    SGz = float(row['z'])

    pt = np.array([[SGx, SGy, SGz]])
    LINE.append("SG")
    LINE.append(index + 1)
    d = ((k.query(pt[0])))
    print(d)
    if d[0] >= 0.5:
        Diss_S_Index.append(index+1)
    else:
        continue

    LINE.append(d[0])
    LINE.append("AU")
    LINE.append(d[1])
    print(LINE)
    line = str(LINE)
    with open ("DISTANCES.txt", "a") as f:
        f.write(line + "\n")


num_diss = (len(Diss_S_Index))
num_total = (len(SG))
Dissociated_Perc = ((num_diss/num_total)*100)

print(Diss_S_Index, "sulfur indexs of dissociated")
print(num_diss, "peptides dissociated")
print(Dissociated_Perc, "% of peptides have moved further than 0.5 nm from the surface ")

First_Diss = [i - 294 for i in Diss_S_Index]
print(First_Diss)

Last_Diss = [i + 297 for i in First_Diss]
print(Last_Diss)

with open('Diss.ndx', 'w') as f:
    f.write("[Dissociated] \n")
    c=1
    for value in First_Diss:
        for index in range (value, value+297):
            if c % 15 == 0:
                f.write(str(index)+ " ")
                f.write("\n")
            else:
                f.write(str(index)+" ")
            c+=1
    f.write("\n")

print("FINISHED LAH")


## Combine the Diss.ndx and original ndx with 0 being diss.
os.system("cat Diss.ndx index.ndx > index2.ndx")
'''
fr = "sed  \'s/inxtc/" +str(xtcin)+ "/g\' index.bash > indexf.bash"
os.system(fr)
fr = "sed -i \'s/tpr/" +str(intpr)+ "/g\' indexf.bash"
os.system(fr)
fr = "sed -i \'s/outxtc/" +str(xtcout)+ "/g\' indexf.bash "
os.system(fr)
'''
os.system("bash index.bash")
