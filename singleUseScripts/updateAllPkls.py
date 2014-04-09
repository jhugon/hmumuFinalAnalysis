#!/usr/bin/env python

import singleHelpers
from helpers import *
import cPickle
import sys
import os
import os.path
import glob
import re
from copy import deepcopy

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

class PklUpdater(object):
  def __init__(self,oldDir,newDir):
    self.oldDir = oldDir
    self.newDir = newDir
    print "Loading old data from ",oldDir
    self.oldData = loadBiasData(oldDir)
    print "Loading new data from ",newDir
    self.newData = loadBiasData(newDir)
    print "Data Loaded."
    energyStr = None
    for cat in self.oldData:
      if energyStr == None:
        energyStr = self.oldData[cat]['meta']['energyStr']
      else:
        if not energyStr == self.oldData[cat]['meta']['energyStr']:
          print "Error: not all old data cats have the same 'energyStr', exiting."
          sys.exit(1)
    for cat in self.newData:
      if not energyStr == self.newData[cat]['meta']['energyStr']:
        print "Error: a new data cat was found with a differing 'energyStr' from the old data, exiting."
        sys.exit(1)
    self.energyStr = energyStr

  def getUpdatedData(self):
    categories = sorted(self.oldData.keys())
    result = {}
    for cat in categories:
      if not self.newData.has_key(cat):
        print "Warning: new data dir '{0}' doesn't have category '{1}' found in old data dir, skipping category".format(self.newDir,cat)
        continue
      oldDataCat = self.oldData[cat]
      newDataCat = self.newData[cat]

      result[cat] = {}
      for refName in oldDataCat:
        if refName == "Bernstein":
          order = self.getCorrectOrder(cat)
          result[cat][refName] = newDataCat["{0}{1}".format(order,refName)]
        else:
          result[cat][refName] = oldDataCat[refName]
    return result

  def writeUpdatedPkls(self,outDir):
    updatedData = self.getUpdatedData()
    for cat in updatedData:
      outFileName = "{0}/biasData_{1}_{2}_jobGrp12345.pkl".format(outDir,cat,self.energyStr)
      outFile = open(outFileName,'w')
      cPickle.dump(updatedData[cat],outFile)
      outFile.close()
    print "Done writing updated pkls."

  def getCorrectOrder(self,cat):
     bernRefOrderDict = {}
     bernRefOrderDict['Jets01PassPtG10BB'] = 4
     bernRefOrderDict['Jets01PassPtG10BO'] = 6
     bernRefOrderDict['Jets01PassPtG10BE'] = 5
     bernRefOrderDict['Jets01PassPtG10OO'] = 4
     bernRefOrderDict['Jets01PassPtG10OE'] = 3
     bernRefOrderDict['Jets01PassPtG10EE'] = 3
     bernRefOrderDict['Jets01FailPtG10BB'] = 4
     bernRefOrderDict['Jets01FailPtG10BO'] = 4 
     bernRefOrderDict['Jets01FailPtG10BE'] = 4 
     bernRefOrderDict['Jets01FailPtG10OO'] = 4 
     bernRefOrderDict['Jets01FailPtG10OE'] = 3 
     bernRefOrderDict['Jets01FailPtG10EE'] = 4 
     bernRefOrderDict['Jet2CutsVBFPass']   = 2 
     bernRefOrderDict['Jet2CutsGFPass']    = 3 
     bernRefOrderDict['Jet2CutsFailVBFGF'] = 4 
     return bernRefOrderDict[cat]

if __name__ == "__main__":
  pklUpdater = PklUpdater("oldAllFuncPkls/","newBernPkls/")
  pklUpdater.writeUpdatedPkls("allFuncPkls/")

  ## Testing:
  #from scipy import *
  #data = loadBiasData("oldAllFuncPkls/")
  #data = data["Jet2CutsVBFPass"]
  #print "oldPkls 'Old': ",median(array(data['Old'][125]['MSSM']['n']) - array(data['Old'][125]['nTrue']))

  #data = loadBiasData("newBernPkls/")
  #data = data["Jet2CutsVBFPass"]
  #print "newPkls '2Bernstein': ",median(array(data['2Bernstein'][125]['MSSM']['n']) - array(data['2Bernstein'][125]['nTrue']))

  #data = loadBiasData("allFuncPkls/")
  #data = data["Jet2CutsVBFPass"]
  #print "upPkls 'Old': ",median(array(data['Old'][125]['MSSM']['n']) - array(data['Old'][125]['nTrue']))
  #print "upPkls 'Bernstein': ",median(array(data['Bernstein'][125]['MSSM']['n']) - array(data['Bernstein'][125]['nTrue']))
