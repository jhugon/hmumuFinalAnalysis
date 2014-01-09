#!/usr/bin/env python

import singleHelpers
from helpers import *
from xsec import *
from effReaderFromFile import *
import cPickle
import sys
import os
import os.path
import glob
import re
from numpy import mean,median

def getProductionPdf(name):
  if "Jets01PassPtG10BB" in name:
    order = 5
  elif "Jets01PassPtG10BO" in name:
    order = 6
  elif "Jets01PassPtG10BE" in name:
    order = 6
  elif "Jets01PassPtG10OO" in name:
    order = 5
  elif "Jets01PassPtG10OE" in name:
    order = 5
  elif "Jets01PassPtG10EE" in name:
    order = 5
  elif "Jets01FailPtG10BB" in name:
    order = 5
  elif "Jets01FailPtG10BO" in name:
    order = 5
  elif "Jets01FailPtG10BE" in name:
    order = 5
  elif "Jets01FailPtG10OO" in name:
    order = 6
  elif "Jets01FailPtG10OE" in name:
    order = 5
  elif "Jets01FailPtG10EE" in name:
    order = 6
  elif "Jet2CutsVBFPass" in name:
    order = 3
  elif "Jet2CutsGFPass" in name:
    order = 5
  elif "Jet2CutsFailVBFGF" in name:
    order = 5
  else:
    print "Error: makePDFBakBernsteinProd: category '"+str(name)+"' not recognized, exiting."
    sys.exit(1)
  return str(order)+"Bernstein"

def getSMSigCounts(catName,higgsMass,energy="8TeV"):
  lumi = lumiDict[energy]
  lumi *= 1000.
  xsTotal = 0.
  for prodMode in ['gg','vbf','wh','zh']:
    name = prodMode + "Hmumu"+str(higgsMass)+"_"+energy
    effReader = effReaderFromFile('fitresults',
                                prodMode,
                                energy,
                                catName)
    efficiencies = effReader.getEff()
    eff = efficiencies[ '{0:.1f}'.format(higgsMass) ]
    higgsMassStr = "{0:.1f}".format(higgsMass)
    tmpName = name.replace("125",higgsMassStr)
    xs = xsec[tmpName]
    xsTotal += xs*eff
  return xsTotal*lumi

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
    muDict[catName] = {'meta':data[catName]['meta']}
    pdfAltNamesDict = data[catName]['meta']['pdfAltNamesDict']
    pdfRefNames = data[catName]['meta']['refPdfNameList']
    for refPdfName in pdfRefNames:
      muDict[catName][refPdfName] = {}
      for altPdfName in pdfAltNamesDict[refPdfName]:
        muDict[catName][refPdfName][altPdfName] = {}
        for hmass in data[catName]['meta']['sigMasses']:
          nSigSM = getSMSigCounts(catName,hmass)
          nSigRef = data[catName][refPdfName][hmass]['nTrue']
          nSigRefUnc = data[catName][refPdfName][hmass]['errTrue']
          nSigAlt = data[catName][refPdfName][hmass][altPdfName]['n']
          nSigAltUnc = data[catName][refPdfName][hmass][altPdfName]['err']
          muAlt = [i/nSigSM for i in nSigAlt]
          muRef = [i/nSigSM for i in nSigRef]
          muR = [(iAlt-iRef)/nSigSM for iAlt,iRef in zip(nSigAlt,nSigRef)]
          muAltUnc = [i/nSigSM for i in nSigAltUnc]
          muRefUnc = [i/nSigSM for i in nSigRefUnc]
          muDict[catName][refPdfName][altPdfName][hmass] = {}
          muDict[catName][refPdfName][altPdfName][hmass]["muR"]      = muR
          muDict[catName][refPdfName][altPdfName][hmass]["muAlt"]    = muAlt
          muDict[catName][refPdfName][altPdfName][hmass]["muAltUnc"] = muAltUnc
          muDict[catName][refPdfName][altPdfName][hmass]["muRef"]    = muRef
          muDict[catName][refPdfName][altPdfName][hmass]["muRefUnc"] = muRefUnc
          muDict[catName][refPdfName][altPdfName][hmass]["nSigSM"]   = nSigSM
  return muDict

def createSummaryMuDict(data):
  muDict = {}
  for catName in data:
    muDict[catName] = {'meta':data[catName]['meta']}
    pdfRefNames = data[catName]['meta']['refPdfNameList']
    altPdfName = getProductionPdf(catName)
    muDict[catName]['meta']['altPdf'] = altPdfName
    muDict[catName]['meta'].pop('pdfAltNamesDict')
    for refPdfName in pdfRefNames:
      muDict[catName][refPdfName] = {}
      for hmass in data[catName]['meta']['sigMasses']:
        nSigSM = data[catName][refPdfName][altPdfName][hmass]["nSigSM"]
        muRList = data[catName][refPdfName][altPdfName][hmass]["muR"]
        muAltList = data[catName][refPdfName][altPdfName][hmass]["muAlt"]
        muAltUncList = data[catName][refPdfName][altPdfName][hmass]["muAltUnc"]
        muRefList = data[catName][refPdfName][altPdfName][hmass]["muRef"]
        muRefUncList = data[catName][refPdfName][altPdfName][hmass]["muRefUnc"]
        muR = median(muRList)
        muAlt = median(muAltList)
        muAltUnc = median(muAltUncList)
        muRef = median(muRefList)
        muRefUnc = median(muRefUncList)
        muDict[catName][refPdfName][hmass] = {}
        muDict[catName][refPdfName][hmass]["muR"]      = muR
        muDict[catName][refPdfName][hmass]["muAlt"]    = muAlt
        muDict[catName][refPdfName][hmass]["muAltUnc"] = muAltUnc
        muDict[catName][refPdfName][hmass]["muRef"]    = muRef
        muDict[catName][refPdfName][hmass]["muRefUnc"] = muRefUnc
        muDict[catName][refPdfName][hmass]["nSigSM"] = nSigSM
  return muDict

if __name__ == "__main__":

  ## Create detailed mu data dicts
  #biasData = loadBiasData("output/")
  #muDict = createMuDict(biasData)
  #muDictFile = open("biasMuDictDetail.pkl",'w')
  #cPickle.dump(muDict,muDictFile)
  #muDictFile.close()

  # Read in detailed mu data dict
  muDictFile = open("biasMuDictDetail.pkl")
  muDict = cPickle.load(muDictFile)
  muDictFile.close()

  #############################
  # Create Summary Dictionary #
  #############################
  # Summary dictionary takes median values of muR, muAlt, 
  # and muAltUnc over all of the toys.  Also, orders
  # besides the ones chosen in the function getProductionPdf
  # are not included.
  muSummaryDict = createSummaryMuDict(muDict)

  #muSummaryDictFile = open("biasMuDictSummary.pkl",'w')
  #cPickle.dump(muSummaryDict,muSummaryDictFile)
  #muSummaryDictFile.close()
  #print muSummaryDict

  # Look at Summary Data
  for catName in muSummaryDict:
    print catName
    for refName in muSummaryDict[catName]:
      if refName == "meta":
        continue
      print "  ",refName
      print "    {0:6} {1:>7} {2:>7} {3:>7} {4:>7} {5:>7}".format('hmass','muR','muAlt','muAltUnc','muRef','nSigSM')
      for hmass in sorted(muSummaryDict[catName][refName].keys()):
        # muR is (N(alt)-N(ref))/N(SM) similar to what we used before
        muR = muSummaryDict[catName][refName][hmass]['muR']  
        # muAlt is N(alt)/N(SM), just in case you want it
        muAlt = muSummaryDict[catName][refName][hmass]['muAlt']
        # muAltUnc is the uncertainty on N(alt) divided by N(SM)
        muAltUnc = muSummaryDict[catName][refName][hmass]['muAltUnc']
        # muRef is N(ref)/N(SM), just in case you want it
        muRef = muSummaryDict[catName][refName][hmass]['muRef']
        # muRefUnc is the uncertainty on N(ref) divided by N(SM)
        # You shouldn't need it
        muRefUnc = muSummaryDict[catName][refName][hmass]['muRefUnc']
        # nSigSM is the expected number of signal events in this category
        nSigSM = muSummaryDict[catName][refName][hmass]['nSigSM']
        print "    {0:6.2f} {1:7.2f} {2:7.2f} {3:7.2f} {4:7.2f} {5:7.2f}".format(hmass,muR,muAlt,muAltUnc,muRef,nSigSM)
