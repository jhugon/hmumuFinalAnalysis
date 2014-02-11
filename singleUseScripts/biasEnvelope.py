#!/usr/bin/env python

import singleHelpers
from helpers import *
import cPickle
import sys
import os
import os.path
import glob
import re
from numpy import mean,median


if __name__ == "__main__":
  fns = glob.glob("biasMaxPkl*.pkl")
  energies = []
  for fn in fns:
    match = re.match("biasMaxPkl_([78P14]+TeV)_signif([.0-9]+).pkl",fn)
    assert(match)
    energyStr = match.group(1)
    if not energyStr in energies:
      energies.append(energyStr)
  energies.sort()
  finalDict = {}
  categorySet = None
  for energy in energies:
    globStr = "biasMaxPkl_"+energy
    dataList = []
    if categorySet == None:
      categorySet = set()
    refNameSet = set()
    for fn in glob.glob(globStr+"*.pkl"):
      match = re.match(globStr+"_signif([.0-9]+).pkl",fn)
      signif = float(match.group(1))
      f = open(fn)
      data = cPickle.load(f)
      f.close()
      dataList.append((data,signif))
      for cat in data:
        if not cat in categorySet:
          categorySet.add(cat)
        for refName in data[cat]:
          if not refName in refNameSet:
            refNameSet.add(refName)
    dataList =  sorted(dataList,key= lambda x: x[1])
    catList = sortCatNames(list(categorySet))
    refNameList = list(refNameSet)
    plainCompareTable = "#########################\n### "+energy+" Max Bias Comparison Table (Nevts)\n\n"
    latexCompareTable = "%%%%%%%%%%%%%%%%%%%%%%%%%\n%%% "+energy+" Max Bias Comparison Table (Nevts)\n\n"
    plainCompareTable += "{0:20}".format("Inject Signif:")
    for dataTup in dataList:
      data, signif = dataTup
      plainCompareTable += "{0:^15.0f}".format(signif)
    plainCompareTable += "\n"
    catMaxDict = {}
    for cat in catList:
      plainCompareTable += "{0:20}".format(cat)
      catMax = 0.
      for dataTup in dataList:
        data, signif = dataTup
        if data.has_key(cat):
          subMax = 0.
          for refName in refNameList:
            if data[cat].has_key(refName):
                tmp = data[cat][refName]
                if tmp > subMax:
                  subMax = tmp
          plainCompareTable += "{0:^15.1f}".format(subMax)
          if subMax > catMax:
            catMax = subMax
        else:
          plainCompareTable += "{0:^15}".format('-')
      catMaxDict[cat] = catMax
      plainCompareTable += "\n"

    plainCompareTable += "\n\n"
    latexCompareTable += "\n\n"

    finalDict[energy] = catMaxDict
        
    print latexCompareTable
    print plainCompareTable

  plainCompareTable = "#########################\n### "+energy+" Max Bias Final Table (Nevts)\n\n"
  latexCompareTable = "%%%%%%%%%%%%%%%%%%%%%%%%%\n%%% "+energy+" Max Bias Final Table (Nevts)\n\n"
  plainCompareTable += "{0:20}".format("Energy:")
  for energy in energies:
    plainCompareTable += "{0:^15}".format(energy)
  plainCompareTable += "\n"
  for cat in sortCatNames(list(categorySet)):
    plainCompareTable += "{0:20}".format(cat)
    for energy in energies:
      plainCompareTable += "{0:^15.1f}".format(finalDict[energy][cat])
    plainCompareTable += "\n"
  plainCompareTable += "\n\n"
  latexCompareTable += "\n\n"
      
  print latexCompareTable
  print plainCompareTable
    
    

