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
from numpy import mean,median

errorsDictFile = open("pklfiles/oneSig.pkl")
ERRORSDICT = cPickle.load(errorsDictFile)
errorsDictFile.close()

def printMaxBias(dataList,sigList,energyStr):
  dataOut = {}
  categories = sortCatNames(dataList[0].keys())
  for data in dataList:
    assert(categories == sortCatNames(data.keys()))

  print "+++++++++++++++++++++++++++++++++++++++++++++"
  print "+++++++++++++++++++++++++++++++++++++++++++++"
  print "+++++++++++++++++++++++++++++++++++++++++++++"
  print

  print "For "+energyStr
  print
  
  lineStr = "{0:24}".format("Sig Injected: ")
  for sig in sigList:
    lineStr+=  "{0:>10}".format("{0}-Sigma".format(sig))
  print lineStr
  
  for cat in categories:
    dataOut[cat] = {}
    refNames = sorted(dataList[0][cat].keys())
    for data in dataList:
      assert(refNames == sorted(data[cat].keys()))
    title = TITLEMAP[cat]
    print ("{0:^"+str(24+len(sigList)*10)+"}").format(title)
    for refName in refNames:
      refTitle = PDFTITLEMAP[refName]
      lineStr = "{0:24}".format(refTitle)
      maxBias = 0.
      for sig,data in zip(sigList,dataList):
        maxBias = max(maxBias,data[cat][refName])
        lineStr += "{0:10.1f}".format(data[cat][refName])
      dataOut[cat][int(refName[0])]=maxBias
      print lineStr

  print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n"*3
  print

  orders = [3,4,5,6]
  lineStr = "{0:24}".format("Orders:")
  for i in orders:
    lineStr += "{0:>10}".format(i)
  print lineStr
  for cat in categories:
    if "PassPtG10BO" in cat or "Jet2CutsVBFPass" in cat or "Jet2CutsGFPass" in cat:
      continue
    lineStr = "{0:24}".format(TITLEMAP[cat])
    for i in orders:
      lineStr += "{0:10.1f}".format(dataOut[cat][i])
    print lineStr

  orders = [5,6,7,8]
  lineStr = "{0:24}".format("Orders:")
  for i in orders:
    lineStr += "{0:>10}".format(i)
  print lineStr
  lineStr = "{0:24}".format("0,1-Jet Tight BO")
  for i in orders:
    lineStr += "{0:10.1f}".format(dataOut["Jets01PassPtG10BO"][i])
  print lineStr

  orders = [2,3,4,5]
  lineStr = "{0:24}".format("Orders:")
  for i in orders:
    lineStr += "{0:>10}".format(i)
  print lineStr
  for cat in ["Jet2CutsVBFPass","Jet2CutsGFPass"]:
    lineStr = "{0:24}".format(TITLEMAP[cat])
    for i in orders:
      lineStr += "{0:10.1f}".format(dataOut[cat][i])
    print lineStr

  print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n"*3
  print
  print r"\begin{tabular}{|l|r|r|r|r|} \hline"

  

  orders = [3,4,5,6]
  lineStr = "{0:24}".format("Bernstein Order")
  for i in orders:
    lineStr += " & {0:>10}".format(i)
  print lineStr+r" \\ \hline"
  for cat in categories:
    if "PassPtG10BO" in cat or "Jet2CutsVBFPass" in cat or "Jet2CutsGFPass" in cat:
      continue
    lineStr = "{0:24}".format(TITLEMAP[cat])
    for i in orders:
      lineStr += " & {0:10.0f}\\%".format(dataOut[cat][i]/ERRORSDICT[energyStr][cat][125]*100.)
    print lineStr+r" \\ \hline"

  orders = [5,6,7,8]
  lineStr = "{0:24}".format("Orders:")
  for i in orders:
    lineStr += " & {0:>10}".format(i)
  print lineStr+r" \\ \hline"
  lineStr = "{0:24}".format("0,1-Jet Tight BO")
  for i in orders:
    lineStr += " & {0:10.0f}\\%".format(dataOut["Jets01PassPtG10BO"][i]/ERRORSDICT[energyStr]["Jets01PassPtG10BO"][125]*100.)
  print lineStr+r" \\ \hline"

  orders = [2,3,4,5]
  lineStr = "{0:24}".format("Orders:")
  for i in orders:
    lineStr += " & {0:>10}".format(i)
  print lineStr+r" \\ \hline"
  for cat in ["Jet2CutsVBFPass","Jet2CutsGFPass"]:
    lineStr = "{0:24}".format(TITLEMAP[cat])
    for i in orders:
      lineStr += " & {0:10.0f}\\%".format(dataOut[cat][i]/ERRORSDICT[energyStr][cat][125]*100.)
    print lineStr+r" \\ \hline"

  print r"\end{tabular}"


if __name__ == "__main__":
  for energyStr in ["7TeV","8TeV"]:
    sigList = [0,3]
    dataList = []
    for signif in sigList:
      fString = "BiasDataVRefFunc_{0}_{1}Sig.pkl".format(energyStr,signif)
      try:
        f = open(fString)
        data = cPickle.load(f)
        f.close()
        dataList.append(data)
      except:
        pass
    if len(dataList) != len(sigList):
      continue
    printMaxBias(dataList,sigList,energyStr)
  
