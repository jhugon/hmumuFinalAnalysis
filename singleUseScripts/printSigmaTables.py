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

# Do match and don't match w/o extensions.  \.expsig and \.sig are added automatically
def getDataSig(fileString,matchString=r"_([-\d.]+)\.txt\.sig",dontMatchStrings=[],doSort=True,getPValue=False,xMax=150,xMin=120):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  #print fileString
  #print matchString
  result = {}
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
    exp = -20.0
    xNum = -10.0
    if match:
      xNum = match.group(1)
      if float(xNum) > xMax:
        continue
      if float(xNum) < xMin:
        continue
    for line in tmpF:
      obsMatch = None
      if getPValue:
        obsMatch = re.search(r"p-value = ([.\deE]+)",line)
      else:
        obsMatch = re.search(r"^Significance:[\s]+([.\deE]+)",line)
      if obsMatch:
        obs = obsMatch.group(1)
    expFname = os.path.splitext(fname)[0]+".expsig"
    """
    try:
      tmpFExp = open(expFname)
      for line in tmpFExp:
        obsMatch = None
        if getPValue:
          obsMatch = re.search(r"p-value = ([.\deE]+)",line)
        else:
          obsMatch = re.search(r"^Significance:[\s]+([.\deE]+)",line)
        if obsMatch:
          exp = obsMatch.group(1)
    except Exception:
      print("Expected Significance Not Found: "+expFname)
    """
    thisPoint = [float(xNum),float(obs),float(exp)]
    #print thisPoint
    if thisPoint.count("-10.0")>0:
        continue
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result[float(xNum)] = thisPoint
  #print result
  return result



if __name__ == "__main__":
  print

  combData = getDataSig("statsInput/CombSplitAll_7P8TeV_*.txt.sig")
  zeroOneJetData = getDataSig("statsInput/Jets01SplitCatAll_7P8TeV_*.txt.sig")
  twoJetData = getDataSig("statsInput/Jet2SplitCutsGFSplit_7P8TeV_*.txt.sig")
  for i in range(120,151):
    i = float(i)
    lineComb = [float('NaN'),float('NaN'),float('NaN')]
    line01 = [float('NaN'),float('NaN'),float('NaN')]
    line2 = [float('NaN'),float('NaN'),float('NaN')]
    if combData.has_key(i):
      lineComb = combData[i]
    if zeroOneJetData.has_key(i):
      line01 = zeroOneJetData[i]
    if twoJetData.has_key(i):
      line2 = twoJetData[i]

    lineArr = [lineComb[1],line01[1],line2[1]]
    lineArr += [scipy.stats.norm.sf(ele) for ele in lineArr]
    lineArr = [int(i)] + lineArr
    print r"{0} & {1:.2f} & {2:.2f} & {3:.2f} & {4:.2f} & {5:.2f} & {6:.2f} \\ \hline".format(*lineArr)
