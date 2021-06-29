#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 17:11:00 2021
@author: joe
"""
## Import all packages

import numpy as np
import pandas as pd
import math
from scipy.spatial import KDTree

##### DEFINING THE SYSTEM

FNO = "SurfaceX.gro"
NUMPEPT = 800 # TOTAL FOR WHOLE PARTICLE
SIZE = 14 # AGAIN WHOLE PARTICLE
FileName = "TEST \n" #
FileLen = 288064 # SHOULD AUTOMATE THIS!

#####

## Functions from paper (cite!!).

# Creates the mesh.

def spherical_coordinate(x,y):

    return [math.cos(x) * math.cos(y),

            math.sin(x) * math.cos(y), math.sin(y)]

# Creates the actual N points.
def NX(n, x):
    pts=[]
    start = (-1. + 1. / (n - 1.) )
    increment = (2. - 2. / (n -1.)) / (n -1.)

    for j in range(0, n):
        s = start + j * increment
        pts.append(
        spherical_coordinate(
        s * x, math.pi / 2. *
        math.copysign(1, s) *
        (1. - math.sqrt(1. -abs(s)))
        ))

    return pts

# Generates complete set of points in a list of lists.
def generate_points(n):
    return NX(n, 0.1 + 1.2 * n)

## My Functions

# For the general points created, multiply by radius.
def xyz(col):

    global points
    points = []
    c = 0
    for v in genplot:
        value = ((genplot[c][col])/2*SIZE)
        points.append(value)
        c+=1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def fs(leng,column,pfn):
    global variable
    variable = []
    f = open (pfn, "r")

    for line in f:
        data = f.readlines()

    variable = []
    value = 0
    while value < leng:
        if value>0:
            dx = data[value].split()
            variable.append(dx[column])
            value+=1
        else:
         value+=1
    return(variable)

# Get all of the peptide info.
def split(leng, pfn):
    #pfn = "peptide.gro"
    peptide_gro = {}

    v = 0
    while v<6:
        peptide_gro["var"+str(v)] = fs(leng,v,pfn)
        v+=1

    global Pep_Name
    Pep_Name = (peptide_gro['var0'])
    global Atom_Name
    Atom_Name = (peptide_gro['var1'])
    global Atom_Number
    Atom_Number = (peptide_gro['var2'])
    global X_Coor
    X_Coor = (peptide_gro['var3'])
    global Y_Coor
    Y_Coor = (peptide_gro['var4'])
    global Z_Coor
    Z_Coor = (peptide_gro['var5'])

## Some vector functions; coudl probably use math instead.
def CP(V1, V2):
    global V3
    cx =  (V1[1]*V2[2])-(V1[2]*V2[1])
    cy =  (V1[2]*V2[0])-(V1[0]*V2[2])
    cz =  (V1[0]*V2[1])-(V1[1]*V2[0])
    V3 = [cx, cy, cz]
    return V3

def DP(V1, V2):
    global V3
    V3 = ((V1[0]*V2[0])+(V1[1]*V2[1])+(V1[2]*V2[2]))
    return V3

def NO(V1):
    global V3
    V3 = math.sqrt((V1[0]**2)+(V1[1]**2)+(V1[2]**2))
    return V3

def MULTI(V1, V2):
    global V3
    V3 = [V1[0]*V2,V1[1]*V2,V1[2]*V2]
    return V3

def PROD(V1, V2):
    global V3
    V3 = [V1[0]*V2[0],V1[1]*V2[1],V1[2]*V2[2]]
    return V3

def DIV(V1, V2):
    global V3
    V3 = [V1[0]/V2, V1[1]/V2, V1[2]/V2]
    return V3

def PLUS(V1, V2):
    global V3
    V3 = [V1[0]+V2[0],V1[1]+V2[1],V1[2]+V2[2]]
    return V3

def MINUS(V1,V2):
    global V3
    V3 = [V1[0]-V2[0],V1[1]-V2[1],V1[2]-V2[2]]
    return V3

## Making the nanoparticle.

with open ("au-cubic.gro", "r") as f:
    with open (FNO, "w") as fo:
        fo.write(FileName)
        fo.write("%5s \n" %(FileLen))
        ln = 0
        b = 0
        for line in f:
            if ln < 2:
                ln+=1
            elif ln > 5627880:
                break
            else:
                col = line.split()
                pp = col[0]
                atom = str(col[1])
                x = float(col[-6])
                y = float(col[-5])
                z = float(col[-4])
                r = math.sqrt((x**2)+(y**2)+(z**2))
                if ln % 2 == 0:
                    if r < SIZE/2:
                        if x > 0 and y > 0 and z > 0:  # # # #for surface# # #
                            fo.write(line)
                            ln+=1
                            b+=1
                        else:
                            ln+=1
                    else:
                        ln+=1
                else:
                    if b == 1:
                        fo.write(line)
                        ln+=1
                        b-=1
                    else:
                        ln+=1

## Generating a mesh
genplot = (generate_points(NUMPEPT))
xyz(0)
XPOINTS = points
xyz(1)
YPOINTS = points
xyz(2)
ZPOINTS = points

lines = []

with open (FNO, "r") as f:
    for num, line in enumerate(f, 1):
        line = [line.rstrip('\n')]
        lines.append(line)
## Remove the header
lines.pop(0);lines.pop(0)
## For each line in lines split it to th seperate parts.

l = []
for n,v in enumerate(lines):

    indices = [0,5,8,12,15,20,22,28,30,36,38,44]
    p = (v[0])
    parts = [p[i:j] for i,j in zip(indices, indices[1:]+[None])]
    #print(parts)
    value = []
    #value.append(n)
    value.append(int(n))
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

## Make new DF called SurfaceAtoms of just surface gold.
AU = Atom_Dict["AUB"]
Pythag = []
SA = []

for index, row in AU.iterrows():
    x = float(row['x'])
    y = float(row['y'])
    z = float(row['z'])
    pyth = (math.sqrt(x**2+y**2+z**2))
    if pyth > ((SIZE/2)-0.5):
        SA.append(row)
        Pythag.append(pyth)
    else:
        continue

SurfaceAtoms = pd.DataFrame(SA)
SurfaceAtoms["Pythag"] = Pythag
SAxyz = SurfaceAtoms[["x","y","z"]].to_numpy()
kd = KDTree(SAxyz)
kout = []
#print(kd.data)
c = 0
Ax = []
Ay = []
Az = []

NNDISTANCES = []
for v in XPOINTS:
    LINE = []
    pt = np.array([[XPOINTS[c], YPOINTS[c], ZPOINTS[c]]])

    d = kd.query(pt[0])
    i = (d[1])
    AUxyz = (kd.data[i])
    Ax.append(AUxyz[0])
    Ay.append(AUxyz[1])
    Az.append(AUxyz[2])

    NNDISTANCES.append(d[0])
    line = str(LINE)

    c+=1

################################################################################
sortedx = []
sortedy = []
sortedz = []



c1 = 0
#NNDISTANCES is the sulfur points.
for v in NNDISTANCES:

    D = 0.3
    nn = v
    #The poitns 00 is zero, AU is the NNgold. So is the current Sulfur pos.
    #Sx is an arbitary point along the axis of 00,S1 used to calc theta.
    OO = [0,0,0]
    AU = [Ax[c1],Ay[c1], Az[c1]]
    So = [XPOINTS[c1], YPOINTS[c1], ZPOINTS[c1]]
    t_v = 1
    V = MINUS(So,OO)
    VN = NO(V)
    U = [((V[0])/VN),((V[1])/VN),((V[2])/VN)]
    DU = MULTI(U,t_v)
    Sx = PLUS(So,DU)
    #Calculating angle Sx So Au. Note the vectors of points So - point.
    AuSo = MINUS(So, AU)
    SxSo = MINUS(So, Sx)
    AuSoSxSo = DP(AuSo, SxSo)
    NOAuSo = NO(AuSo)
    NOSxSo = NO(SxSo)
    P = (NOAuSo * NOSxSo)
    BRA = (AuSoSxSo/P)
    AngleSo = math.acos(BRA)
    #Calculate AngleS1
    AngleS1 = math.sin((nn*(math.sin(AngleSo)))/D)
    #Caculate AngleAU
    AngleAU = math.pi - AngleS1 - AngleSo
    d = ((D*(math.sin(AngleAU)))/(math.sin(AngleSo)))
    #Adding d to So along the axis.
    V = MINUS(So,OO)
    VN = NO(V)
    U = [((V[0])/VN),((V[1])/VN),((V[2])/VN)]
    DU = MULTI(U,d)
    XYZ_NEW = PLUS(So,DU)

    #Now rotate each point on the peptide to the appropriate pos.
    split(299, "peptide.gro")
    num = SIZE/2
    OS  = [num,0,0]
    #OSf = [XPOINTS[c1], YPOINTS[c1], ZPOINTS[c1]]

    OSf = (XYZ_NEW)
    CPosf = CP(OS, OSf)
    Nosf = NO(CPosf)
    k = DIV(CPosf,Nosf)
    #print(k)
    thetaprod = (DP(OS, OSf))/ (NO(OS) * NO(OSf))
    theta = math.acos(thetaprod)
    c2 = 0
    #The peptide
    #the peptide points



    peptide_in = []
    PNs = []
    ANas = []
    c2s = []
    xs = []
    ys = []
    zs = []
    
    
    
    
    for v in X_Coor:
        with open (FNO, "a") as f:
            xcoor = float(X_Coor[c2])+(SIZE/2)
            ycoor = float(Y_Coor[c2])
            zcoor = float(Z_Coor[c2])
            OA = [xcoor, ycoor, zcoor]
            SA = MINUS(OA, OS)
            SfAf1 = MULTI(SA,math.cos(theta))

            SfAf2_1 = MULTI(SA,math.sin(theta))
            SfAf2_2 = CP(k,SfAf2_1)

            SfAf3_1 = DP(k,SA)
            SfAf3_2 = MULTI(k, SfAf3_1)
            SfAf3_3 = MULTI(SfAf3_2, (1-math.cos(theta)))

            SfAf_pre = PLUS(SfAf3_3, SfAf2_2)
            SfAf = PLUS(SfAf_pre, SfAf1)

            OAf = PLUS(OSf, SfAf)

            PN = (Pep_Name[c2])
            ANa = (Atom_Name[c2])
            ANu = (Atom_Number[c2])
            x = float(OAf[0])
            y = float(OAf[1])
            z = float(OAf[2])
            
            sortedx.append(x)
            sortedy.append(y)
            sortedz.append(z)

            if x > 0 and y > 0 and z > 0:
                #print(x,y,z)
                PNs.append(PN)
                ANas.append(ANa)
                c2s.append(c2)
                xs.append(x)
                ys.append(y)
                zs.append(z)

            else:
                PNs.clear()
                ANas.clear()
                c2s.clear()
                xs.clear()
                ys.clear()
                zs.clear()
                break
            
        c2+=1
                    
                
    with open(FNO, 'a') as f:
        for c,v in enumerate(PNs):
            #print("%8s%7s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f \n" % (PNs[c], ANas[c], c2s[c],xs[c],ys[c],zs[c],0,0,0))
            f.write("%8s%7s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f \n" % (PNs[c], ANas[c], c2s[c],xs[c],ys[c],zs[c],0,0,0))
    c1+=1
    
with open (FNO, "a") as f:
   f.write("%10.5f%10.5f%10.5f \n" % (30,30,30))   
    
print(peptide_in)    
xyzDF = pd.DataFrame(peptide_in)
print(xyzDF)
print(xyzDF.describe)


xyzDF = pd.DataFrame(
    {'xxx': sortedx,
     'yyy': sortedy,
     'zzz': sortedz
    })

print(xyzDF)
print(xyzDF.describe())

