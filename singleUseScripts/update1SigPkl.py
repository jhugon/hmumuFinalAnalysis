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
from copy import deepcopy
from singleUseScripts.biasPklToMu import getSMSigCounts

newTxtF = open("2JetVBFTight1SigBands.txt")

newData = {}
thisMass = None

for line in newTxtF:
  if line[:4] == "==> ":
    match = re.match(r"^==> [0-9A-Za-z]+_7TeV_([0-9]+)\.0\.txt\.mu <==$",line)
    assert(match)
    thisMass = int(match.group(1))
  if line[:5] == "Best ":
    match = re.search(r"[0-9.]+/\+([0-9.]+)[\s]+\(68% CL\)",line)
    if not match:
      print line
      continue
    thisLimit = float(match.group(1))
    newData[thisMass] = thisLimit*getSMSigCounts("Jet2CutsVBFPass",thisMass,"7TeV")
    print thisMass,newData[thisMass]
newTxtF.close()

newData[150] = 33.*getSMSigCounts("Jet2CutsVBFPass",thisMass,"7TeV")

oldPklF = open("pklfiles/oneSig.pkl")
data = cPickle.load(oldPklF)
oldPklF.close()

for mass in [120,125,130,135,140,150]:
  print "{0} {1:10.2f} {2:10.2f}".format(mass,data['7TeV']['Jet2CutsVBFPass'][mass],newData[mass])
  #print "{0} {1:10.2f} {2:10.2f}".format(mass,data['7TeV']['Jet2CutsVBFPass'][mass]/getSMSigCounts("Jet2CutsVBFPass",mass,"7TeV"),newData[mass]/getSMSigCounts("Jet2CutsVBFPass",mass,"7TeV"))

data['7TeV']['Jet2CutsVBFPass'] = newData
newPklF = open("pklfiles/oneSig.pkl.new",'w')
cPickle.dump(data,newPklF)
newPklF.close()
