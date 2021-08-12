#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 13:36:53 2021

@author: joe
"""

import numpy as np
import matplotlib as plt
import pandas as pd




for R in range(0,24):
    #print(R)
    if R < 1:
        FN =  "energyR%i.xvg" %R
        REP = "R%i" %R
        pre_settings = []
        settings = []
        data = []
        legends = []
        
        with open(FN, "r") as f:
            for line in f:
                line = str(line)
                if line.startswith("#"):
                    continue
                elif line.startswith("@"):
                    pre_settings.append(line)
                else:
                    data.append(line.split())
        for line in pre_settings:
            #print(line)
            if '"' in line:
                settings.append(line)      
        for line in settings:
            line = line.split()
            if line[1] == "xaxis":
                xaxis = line[3] + line[4]
                xaxis = xaxis[1:-1]
            elif line[1] == "yaxis":
                yaxis = line[3]
                yaxis = yaxis[1:-1]
            elif line[2] == "legend":
                legends.append(line[3])
            else:
                continue
    

        DATAFRAMES = {}
        t = []
        for i in data:
            t.append(i[0])
        
        
        for n, energy in enumerate(legends,1):
            d = []
            for i in data:
                d.append(i[n])
            #now d is full
            fin = list(zip(t,d))
            df = pd.DataFrame(fin, columns=["Time(ps)", "R%s" %R])

            DATAFRAMES[str(energy[1:-1])] = df
        #print(DATAFRAMES)         
    else:
        FN =  "energyR%i.xvg" %R
        REP = "R%i" %R
        pre_settings = []
        settings = []
        data = []
        legends = []
        
        with open(FN, "r") as f:
            for line in f:
                line = str(line)
                if line.startswith("#"):
                    continue
                elif line.startswith("@"):
                    pre_settings.append(line)
                else:
                    data.append(line.split())
        for line in pre_settings:
            #print(line)
            if '"' in line:
                settings.append(line)      
        for line in settings:
            line = line.split()
            if line[1] == "xaxis":
                xaxis = line[3] + line[4]
                xaxis = xaxis[1:-1]
            elif line[1] == "yaxis":
                yaxis = line[3]
                yaxis = yaxis[1:-1]
            elif line[2] == "legend":
                legends.append(line[3])
            else:
                continue
    

        
        for n, energy in enumerate(legends,1):
            d = []
            for i in data:
                d.append(i[n])
            #now d is full
            NAME  = str(energy[1:-1])
            dicto = (DATAFRAMES[NAME])

            dicto['R%s'%R] = d
            #######df = pd.DataFrame(d, columns=["Time(ps)", "R%s" %R])

            #########DATAFRAMES[str(energy[1:-1])] = df

KEYS = ([*DATAFRAMES])         
print(KEYS[0])
print(DATAFRAMES[KEYS[0]])
POTENTIAL_PLOT = (DATAFRAMES[KEYS[0]])
            
