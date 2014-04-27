#!/usr/bin/env python

import singleHelpers
from helpers import *
from xsec import *
import cPickle
import sys
import os
import os.path
import glob
import re

orderDef = [
    "CombSplitAll",
    "Jets01SplitCatAll",
    "Jet2SplitCutsGFSplit",
    "Jets01PassCatAll" ,
    "Jets01FailCatAll" ,

    "Jets01PassPtG10BB",
    "Jets01PassPtG10BO",
    "Jets01PassPtG10BE",
    "Jets01PassPtG10OO",
    "Jets01PassPtG10OE",
    "Jets01PassPtG10EE",
                          
    "Jets01FailPtG10BB",
    "Jets01FailPtG10BO",
    "Jets01FailPtG10BE",
    "Jets01FailPtG10OO",
    "Jets01FailPtG10OE",
    "Jets01FailPtG10EE",

    "Jet2CutsVBFPass",
    "Jet2CutsGFPass",
    "Jet2CutsFailVBFGF",
  ]

def getData(fileString,matchString=r"_([-\d.]+)\.txt\.out",dontMatchStrings=[],doSort=True):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  result = []
  fNames =  glob.glob(fileString)
  if doSort:
    fNames.sort(key=sortfun)
  #print fileString
  #print fNames
  for fname in fNames: 
    dontMatch = False
    for dont in dontMatchStrings:
      if re.search(dont,fname):
        dontMatch = True
    if dontMatch:
        continue
    tmpF = open(fname)
    match = re.search(matchString,fname)
    obs = -10.0
    median = -10.0
    low2sig = -10.0
    low1sig = -10.0
    high1sig = -10.0
    high2sig = -10.0
    xNum = -10.0
    if match:
      xNum = match.group(1)
    for line in tmpF:
      obsMatch = re.search(r"Observed[\s]Limit:[^.\d]*< ([.\deE]+)",line)
      low2sigMatch = re.search(r"Expected.*2\.5.:[^.\d]*< ([.\deE]+)",line)
      low1sigMatch = re.search(r"Expected.*16.0[^.\d]*< ([.\deE]+)",line)
      medianMatch = re.search(r"Expected.*50\.0.*< ([.\deE]+)",line)
      high1sigMatch = re.search(r"Expected.*84.0.*< ([.\deE]+)",line)
      high2sigMatch = re.search(r"Expected.*97.5.*< ([.\deE]+)",line)
      if obsMatch:
        obs = obsMatch.group(1)
      if low2sigMatch:
        low2sig = low2sigMatch.group(1)
      if low1sigMatch:
        low1sig = low1sigMatch.group(1)
      if medianMatch:
        median = medianMatch.group(1)
      if high1sigMatch:
        high1sig = high1sigMatch.group(1)
      if high2sigMatch:
        high2sig = high2sigMatch.group(1)
    median = float(median)
    obs = float(obs)
    high1sig = float(high1sig)
    low1sig = float(low1sig)
    high2sig = float(high2sig)
    low2sig = float(low2sig)
    high1sig = abs(high1sig-median)
    low1sig = abs(low1sig-median)
    high2sig = abs(high2sig-median)
    low2sig = abs(low2sig-median)
    thisPoint = [xNum,median,obs,high1sig,low1sig,high2sig,low2sig]
    if thisPoint.count("-10.0")>0:
        continue
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
      
  return result

if __name__ == "__main__":
  print

  combData = getData("statsInput/CombSplitAll_7P8TeV_*.txt.out")
  print r"\begin{tabular}{ | c | c | c | c | c |} \hline"
  print r"Mass $[\GeVcc{}]$ & expected & observed & 1 $\sigma$ band &  2 $\sigma$ band \\ \hline"
  for line in combData:
    #print r"{0:.0f} & {1:.1f} & X.X & $^{{+{3:0.1f}}}_{{-{4:0.1f}}}$ & $^{{+{5:0.1f}}}_{{-{6:0.1f}}}$ \\ \hline".format(*[float(i) for i in line])
    print r"{0:.0f} & {1:.1f} & {2:0.1f} & $^{{+{3:0.1f}}}_{{-{4:0.1f}}}$ & $^{{+{5:0.1f}}}_{{-{6:0.1f}}}$ \\ \hline".format(*[float(i) for i in line])
  print r"\end{tabular}"

  
  print

  catsData = getData("statsInput/*_7P8TeV_125.0.txt.out",matchString=r"([0-9a-zA-Z]+)_7P8TeV_125\.0\.txt\.out",doSort=False)
  catsData.sort(key=lambda x: orderDef.index(x[0]))
  for i in range(len(catsData)):
    catsData[i][0] = TITLEMAP[catsData[i][0]]
  print r"\begin{tabular}{ | l | c | c | c | c |} \hline"
  print r"Category(s) & expected & observed & 1 $\sigma$ band &  2 $\sigma$ band \\ \hline"
  for line in catsData:
    #print r"{0} & {1:.1f} & X.X & $^{{+{3:0.1f}}}_{{-{4:0.1f}}}$ & $^{{+{5:0.1f}}}_{{-{6:0.1f}}}$ \\ \hline".format(*line)
    print r"{0} & {1:.1f} & {2:0.1f} & $^{{+{3:0.1f}}}_{{-{4:0.1f}}}$ & $^{{+{5:0.1f}}}_{{-{6:0.1f}}}$ \\ \hline".format(*line)
  print r"\end{tabular}"

  print
