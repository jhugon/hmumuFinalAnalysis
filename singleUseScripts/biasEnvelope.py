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

errorsDictFile = open("pklfiles/oneSig.pkl")
ERRORSDICT = cPickle.load(errorsDictFile)
errorsDictFile.close()

def biasEnvelope(useRefFuncs,relativeForMH=125):
  print "+++++++++++++++++++++++++++++++++++++++++++++"
  print "+++++++++++++++++++++++++++++++++++++++++++++"
  print "+++++++++++++++++++++++++++++++++++++++++++++"
  print "BIAS ENVELOPE FOR THESE REFERENCE FUNCTIONS:"
  print useRefFuncs
  print
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
    dataList =  sorted(dataList,key= lambda x: x[1])
    catList = sortCatNames(list(categorySet))
    refNameList = useRefFuncs
    plainCompareTable = "#########################\n### "+energy+" Max Bias Comparison Table (Nevts)\n"
    plainCompareTable += "Included Reference Functions:"
    for refName in refNameList:
      plainCompareTable += " "+refName
    plainCompareTable += "\n\n"
    plainCompareTable += "{0:20}".format("Inject Signif:")
    latexCompareTable = "%%%%%%%%%%%%%%%%%%%%%%%%%\n%%% "+energy+" Max Bias Comparison Table (Nevts)\n\n"
    latexCompareTable += r"\begin{tabular}{|l|"+"r|"*len(refNameList)+"} \\hline \n"
    latexCompareTable += "{0}".format("Injected Significance")
    for dataTup in dataList:
      data, signif = dataTup
      plainCompareTable += "{0:^15.0f}".format(signif)
      latexCompareTable += " & ${0:.0f}\sigma$".format(signif)
    plainCompareTable += "\n"
    latexCompareTable += " \\\\ \\hline \\hline \n"
    catMaxDict = {}
    for cat in catList:
      plainCompareTable += "{0:20}".format(cat)
      latexCompareTable += TITLEMAP[cat]
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
            else:
              print "Error: refName: {0} not found for {1} {2} {3}sigma, exiting.".format(refName,cat,energy,signif)
              sys.exit(1)
          plainCompareTable += "{0:^15.1f}".format(subMax)
          latexCompareTable += " & {0:.1f}".format(subMax)
          if subMax > catMax:
            catMax = subMax
        else:
          plainCompareTable += "{0:^15}".format('-')
          latexCompareTable += " & - "
      catMaxDict[cat] = catMax
      plainCompareTable += "\n"
      latexCompareTable += " \\\\ \\hline \n"

    plainCompareTable += "\n\n"
    latexCompareTable += "\\end{tabular}\n"
    latexCompareTable += "\\\\ Included Reference Functions: "
    for refName in refNameList:
      latexCompareTable += PDFTITLEMAP[refName]+", "
    latexCompareTable = latexCompareTable[:-2] + "\n\n"

    finalDict[energy] = catMaxDict
        
    print latexCompareTable
    print plainCompareTable

  plainCompareTable = "#########################\n###  Max Bias Final Table (Nevts)\n"
  plainCompareTable += "Included Reference Functions:"
  for refName in refNameList:
    plainCompareTable += " "+refName
  plainCompareTable += "\n\n"
  latexCompareTable = "%%%%%%%%%%%%%%%%%%%%%%%%%\n%%%  Max Bias Final Table (Nevts)\n\n"
  latexCompareTable += r"\begin{tabular}{|l|"+"r|"*len(energies)+"} \\hline \n"
  relativeCompareTable = "%%%%%%%%%%%%%%%%%%%%%%%%%\n%%%  Max Bias Final Table (Nevts/deltaNevts) deltaNevts for mH={0}\n\n".format(relativeForMH)
  relativeCompareTable += r"\begin{tabular}{|l|"+"r|"*len(energies)+"} \\hline \n"
  plainCompareTable += "{0:20}".format("Energy:")
  latexCompareTable += "{0}".format("Energy")
  relativeCompareTable += "{0}".format("Energy")
  for energy in energies:
    plainCompareTable += "{0:^15}".format(energy)
    latexCompareTable += " & {0:}".format(energy)
    relativeCompareTable += " & {0:}".format(energy)
  plainCompareTable += "\n"
  latexCompareTable += "\\\\ \\hline \\hline \n"
  relativeCompareTable += "\\\\ \\hline \\hline \n"
  for cat in sortCatNames(list(categorySet)):
    plainCompareTable += "{0:20}".format(cat)
    latexCompareTable += "{0:}".format(TITLEMAP[cat])
    relativeCompareTable += "{0:}".format(TITLEMAP[cat])
    for energy in energies:
      plainCompareTable += "{0:^15.1f}".format(finalDict[energy][cat])
      latexCompareTable += " & {0:.1f}".format(finalDict[energy][cat])
      relativeCompareTable += " & {0:.1f}".format(finalDict[energy][cat]/ERRORSDICT[energy][cat][relativeForMH])
    plainCompareTable += "\n"
    latexCompareTable += " \\\\ \\hline \n"
    relativeCompareTable += " \\\\ \\hline \n"
  plainCompareTable += "\n\n"
  latexCompareTable += "\\end{tabular}\n"
  relativeCompareTable += "\\end{tabular}\n"
  latexCompareTable += "\\\\ Included Reference Functions: "
  relativeCompareTable += "\\\\ Included Reference Functions: "
  for refName in refNameList:
    latexCompareTable += PDFTITLEMAP[refName]+", "
    relativeCompareTable += PDFTITLEMAP[refName]+", "
  latexCompareTable = latexCompareTable[:-2] + "\n\n"
  relativeCompareTable = relativeCompareTable[:-2] + "\n\n"
      
  print latexCompareTable
  print plainCompareTable
  print relativeCompareTable
    
  dictName = "BakParameterizationUncDict"
  dictCompareTable = "#########################\n###  Max Bias Final Dict (Nevts)\n\n"
  dictCompareTable += dictName+" = {'7TeV':{},'8TeV':{}}\n"
  for energy in energies:
    for cat in sortCatNames(list(categorySet)):
      dictCompareTable += dictName+"['"+energy+"']"+"{0:<21}".format("['"+cat+"']")+" = {0:.2f}\n".format(finalDict[energy][cat])
    
  print dictCompareTable


if __name__ == "__main__":
  onlyVoitRefs = ["Old","VoigtPMm2","VoigtPExpMm2"]
  voitAndSMRefs = ["Old","ExpMOverSq","VoigtPMm2","VoigtPExpMm2"]
  allRefs = ["Old","ExpMOverSq","SumExp","VoigtPMm2","Bernstein","VoigtPExpMm2"]

  biasEnvelope(onlyVoitRefs)
  biasEnvelope(voitAndSMRefs)
  biasEnvelope(allRefs)
  
