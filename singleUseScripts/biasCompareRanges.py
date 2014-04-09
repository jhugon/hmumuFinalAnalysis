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

errorsDictFile = open("pklfiles/oneSig.pkl")
ERRORSDICT = cPickle.load(errorsDictFile)
errorsDictFile.close()

def mergeDictMoreToys(data,newData,multiJob=False):
  newDataKeys = newData.keys()
  dataKeys = data.keys()
  for key in newDataKeys:
    if key == 'meta':
      data[key]['nData'] = newData[key]['nData']
      if multiJob:
        data[key]['nToys'] += newData[key]['nToys']
    elif not (key in dataKeys):
      data[key] = newData[key]
    else:
      masses = sorted(newData[key].keys())
      massesOld = sorted(data[key].keys())
      assert(len(masses)==len(massesOld))
      for old,new in zip(masses,massesOld):
        if old != new:
          print "masses not equal:",old,new
          assert(False)
      for mass in masses:
        subkeys = sorted(newData[key][mass])
        subkeysOld = sorted(data[key][mass])
        assert(len(subkeys) == len(subkeysOld))
        for old,new in zip(subkeys,subkeysOld):
          if old != new:
            print "subkeys not equal:",old,new
            assert(False)
        for subkey in subkeys:
          if subkey == "orderTrue":
            assert(newData[key][mass][subkey] == data[key][mass][subkey])
          elif type(newData[key][mass][subkey])==list:
            assert(type(data[key][mass][subkey])==list)
            data[key][mass][subkey].extend(newData[key][mass][subkey])
          else:
            subsubkeys = sorted(newData[key][mass][subkey])
            subsubkeysOld = sorted(data[key][mass][subkey])
            assert(len(subsubkeys) == len(subsubkeysOld))
            for old,new in zip(subsubkeys,subsubkeysOld):
              if old != new:
                print "subsubkeys not equal:",old,new
                assert(False)
            for subsubkey in subsubkeys:
              assert(type(data[key][mass][subkey][subsubkey])==list)
              assert(type(newData[key][mass][subkey][subsubkey])==list)
              data[key][mass][subkey][subsubkey].extend(newData[key][mass][subkey][subsubkey])

def loadBiasData(pklDir):
  allData = {}
  inputPklFiles = glob.glob(pklDir+"*.pkl")
  foundJobGroupPkl = False
  foundNotJobGroupPkl = False
  for inputPkl in inputPklFiles:
    if "_jobGrp" in inputPkl:
      foundJobGroupPkl = True
    else:
      foundNotJobGroupPkl = True
  if foundNotJobGroupPkl and foundJobGroupPkl:
      print "Error: found .pkl files containing '_jobGrp' and not containing '_jobGrp'."
      print "  Can only process one or the other.  "
      print "  Please delete one or the other and try again. Exiting."
      sys.exit(1)
  if foundNotJobGroupPkl:
    for inputPkl in inputPklFiles:
      print "Running over input pkl file: "+inputPkl
      tmpF = open(inputPkl)
      tmpD = cPickle.load(tmpF)
      allData[tmpD['meta']['catName']] = tmpD
      tmpF.close()
  else:
    # Identify basenames to combine job groups
    basenames = set()
    for inputPklFn in inputPklFiles:
      match = re.match(r"(.+jobGrp)[\d]+\.pkl",inputPklFn)
      assert(match)
      tmpBase = match.group(1)
      if not (tmpBase+"*.pkl") in basenames:
          basenames.add((tmpBase+"*.pkl"))
    for globStr in basenames:
      fns = glob.glob(globStr)
      resultData = None
      for tmpFn in fns:
        tmpF = open(tmpFn)
        tmpD = cPickle.load(tmpF)
        if resultData == None:
          resultData = tmpD 
        else:
          mergeDictMoreToys(resultData,tmpD,True)
        tmpF.close()
      allData[resultData['meta']['catName']] = resultData
  return allData


if __name__ == "__main__":
  from scipy import *
  data = loadBiasData("output/")
  massGroups = {}
  massGroups["120to150"] = [120,125,130,135,140,145,150]
  massGroups["120to140"] = [120,125,130,135,140]
  masses = set()
  for grp in massGroups:
    masses = masses.union(set(massGroups[grp]))
  nBiasResult = {}

  energyStr = None
  signif = None
  for cat in data:
    dataCat = data[cat]
    if energyStr == None:
      energyStr = dataCat['meta']['energyStr']
    elif dataCat['meta']['energyStr'] != energyStr:
      print "Error: energyStr for category '{0}' doesn't match other categories, exiting."
      sys.exit(1)
    if signif == None:
      signif = dataCat['meta']['sigInjectNsigma']
    elif signif != dataCat['meta']['sigInjectNsigma']:
        print "Error: sigInjectNsigma for category '{0}' doesn't match other categories, exiting."
        sys.exit(1)
    nBiasResult[cat] = {}
    for grp in massGroups:
      nBiasResult[cat][grp] = 0.
    refFuncs = dataCat.keys()
    refFuncs.pop(refFuncs.index("meta"))
    for refFunc in refFuncs:
      for mass in masses:
        nBias = median(array(dataCat[refFunc][mass]['MSSM']['n']) - array(dataCat[refFunc][mass]['nTrue']))
        nBias = abs(nBias)
        for grp in massGroups:
          if mass in massGroups[grp]:
            nBiasResult[cat][grp] = max(nBiasResult[cat][grp],nBias)

  print "For {0} {1} sigma signal injected".format(energyStr,signif)

  massGroupList = sorted(massGroups)
  line = "{0:20}".format("")
  for grp in massGroupList:
    line += "{0:>10}".format(grp.replace("to","-"))
  print line
  for cat in sortCatNames(nBiasResult.keys()):
    line = "{0:20}".format(TITLEMAP[cat])
    for grp in massGroupList:
      nBias = nBiasResult[cat][grp]
      statErr = ERRORSDICT[energyStr][cat][125.]
      relBias = nBias / statErr
      line += "{0:>9.0f}%".format(relBias*100.)
    print line 

  print
  print "Nbias:"
  line = "{0:20}".format("")
  for grp in massGroupList:
    line += "{0:>10}".format(grp.replace("to","-"))
  print line
  for cat in sortCatNames(nBiasResult.keys()):
    line = "{0:20}".format(TITLEMAP[cat])
    for grp in massGroupList:
      nBias = nBiasResult[cat][grp]
      line += "{0:>10.3g}".format(nBias)
    print line 
