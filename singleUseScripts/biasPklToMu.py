#!/usr/bin/env python

import cPickle
import sys
import os
import os.path
import glob
import re
from numpy import mean,median

def mergeDicts(data,newData,multiJob=False):
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

##############################################333
##############################################333
##############################################333
##############################################333

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
          mergeDicts(resultData,tmpD,True)
        tmpF.close()
      allData[resultData['meta']['catName']] = resultData
  return allData

def createMuDict(data):
  muDict = {}
  for catName in data:
    muDict[catName] = {}
    pdfAltNamesDict = data[catName]['meta']['pdfAltNamesDict']
    for refPdfName in pdfAltNamesDict:
      muDict[catName][refPdfName] = {'meta':data[catName]['meta']}
      for altPdfName in pdfAltNamesDict[refPdfName]:
        muDict[catName][refPdfName][altPdfName] = {}
        for hmass in data[catName]['meta']['sigMasses']:
          nSigSM = 1.  ## Needs to be found
          nSigRef = data[catName][refPdfName][hmass]['nTrue']
          nSigAlt = data[catName][refPdfName][hmass][altPdfName]['n']
          nSigAltUnc = data[catName][refPdfName][hmass][altPdfName]['err']
          muAlt = [i/nSigSM for i in nSigAlt]
          muR = [(iAlt-iRef)/nSigSM for iAlt,iRef in zip(nSigAlt,nSigRef)]
          muAltUnc = [i/nSigSM for i in nSigAltUnc]
          muDict[catName][refPdfName][altPdfName][hmass] = {}
          muDict[catName][refPdfName][altPdfName][hmass]["muR"]      = muR
          muDict[catName][refPdfName][altPdfName][hmass]["muAlt"]    = muAlt
          muDict[catName][refPdfName][altPdfName][hmass]["muAltUnc"] = muAltUnc
  return muDict

def createSummaryMuDict(data)
  muDict = {}
  for catName in data:
    muDict[catName] = {}
    pdfAltNamesDict = data[catName]['meta']['pdfAltNamesDict']
    for refPdfName in pdfAltNamesDict:
      muDict[catName][refPdfName] = {'meta':data[catName]['meta']}
      for altPdfName in pdfAltNamesDict[refPdfName]:
        muDict[catName][refPdfName][altPdfName] = {}
        for hmass in data[catName]['meta']['sigMasses']:
          muRList = data[catName][refPdfName][altPdfName][hmass]["muR"]
          muAltList = data[catName][refPdfName][altPdfName][hmass]["muAlt"]
          muAltUncList = data[catName][refPdfName][altPdfName][hmass]["muAltUnc"]
          muR = median(muRList)
          muAlt = median(muAltList)
          muAltUnc = median(muAltUncList)
          muDict[catName][refPdfName][altPdfName][hmass] = {}
          muDict[catName][refPdfName][altPdfName][hmass]["muR"]      = muR
          muDict[catName][refPdfName][altPdfName][hmass]["muAlt"]    = muAlt
          muDict[catName][refPdfName][altPdfName][hmass]["muAltUnc"] = muAltUnc
  return muDict

if __name__ == "__main__":

  biasData = loadBiasData("output/")
  muDict = createMuDict(biasData)
  muSummaryDict = createSummaryMuDict(muDict)
  print muSummaryDict
