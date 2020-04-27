#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 18:38:38 2019

@author: rohitagarwal
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn-whitegrid')
sns.set_style("white")

#Import data
reaction = pd.read_csv("example.csv")

#creates a Pandas dataframe
df = pd.DataFrame(reaction)

# Create Fig and gridspec
fig = plt.figure(figsize=(7, 7), dpi= 80)
grid = plt.GridSpec(3, 1, height_ratios=[1.5, 5, 1.5], hspace=0, wspace=0.2)

# Define the axes
ax_top = fig.add_subplot(grid[0, 0], xticklabels=[])
ax_main = fig.add_subplot(grid[1, 0],xticklabels=[])
ax_bottom = fig.add_subplot(grid[2, 0])

#define colors for the plot
cpos = [[0.75, 0.9, 0.75]]
ctor = [[0.9, 0.68, 0.23]]
crea = [[1.0, 0, 0]]
cneg = [[0.75, 0.75, 0.95]]

#group by tag and loop !!! Both If with mulplot[0] and count works
xstart = 1150
xend = 1800
binwidth = 40
count = 0
cutdat = 0
bins = np.arange(xstart, xend, binwidth)
for mulplot in df.groupby('Tag'):

   xpos = mulplot[1]["PosBurstPosition"][df['PosBurstRate']>cutdat]
   ypos = mulplot[1]["PosBurstRate"][df['PosBurstRate']>cutdat]
   xneg = mulplot[1]["NegBurstPosition"][df['NegBurstRate']<-cutdat]
   yneg = mulplot[1]["NegBurstRate"][df['NegBurstRate']<-cutdat]
   if mulplot[0] == "noBreak": 

       ax_main.scatter(xpos,ypos,s = 20,c = cpos)
       ax_main.scatter(xneg,yneg,s = 20,c = cneg)
       ax_top.hist([xpos], bins=bins, histtype='stepfilled', orientation='vertical', color=cpos)  #top histogram
       ax_bottom.hist([xneg], bins=bins, histtype='stepfilled', orientation='vertical', color=cneg) #bottom histogram
       ax_bottom.invert_yaxis()
   if mulplot[0] == "reactionBreak":

#       ax_main.scatter(xpos,ypos,s = 20,c = cpos)
#       ax_main.scatter(xneg,yneg,s = 20,c = cneg)
      ax_main.scatter(xpos,ypos,s = 50,c = crea, edgecolors='black', linewidths=0.5)
      ax_main.scatter(xneg,yneg,s = 50,c = crea, edgecolors='black', linewidths=0.5)
   if mulplot[0] == "torqueBreak":

       ax_main.scatter(xpos,ypos,s = 30,c = ctor, edgecolors='black', linewidths=0.5)
       ax_main.scatter(xneg,yneg,s = 30,c = ctor, edgecolors='black', linewidths=0.5)
       
ax_main.set_xlim((xstart, xend))
ax_main.set_ylim((-1.7, 3))             
ax_top.set_xlim((xstart, xend))
ax_bottom.set_xlim((xstart, xend))
#plt.rcParams.update({'font.size': 20})
plt.rc('font', size=20)
plt.savefig("bubblepos.png")
plt.show()