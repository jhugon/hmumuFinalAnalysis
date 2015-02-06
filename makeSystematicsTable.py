#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import sys
from math import ceil

mass = '125'

energies=[
'7TeV',
'8TeV'
]

massCut = "dimuonMass > 110. && dimuonMass < 170."
jet01Cuts = " && !(jetLead_pt>40. && jetSub_pt>30. && ptMiss<40.)"
jet2Cuts = " && jetLead_pt>40. && jetSub_pt>30. && ptMiss<40."

datasets = ["gg","vbf"]
categories = [
  "Jets01PassPtG10",
  "Jets01FailPtG10",
  "Jet2CutsVBFPass",
  "Jet2CutsGFPass",
  "Jet2CutsFailVBFGF"
]

combNames = ["QCDscale","pdf_","_j"]
combNameMap = {
  "QCDscale":"Scale",
  "pdf_":"PDF",
  "_j":"Jet",
}

######################################
######################################
######################################
######################################
######################################

allNuisanceKeys = sorted(nuisanceMap.keys())
#allNuisanceKeys.remove(allNuisanceKeys.find("QCDscale_VH"))
allNuisanceKeys.remove("QCDscale_VH")

# Now for a tables
for iTable in range(1):
  nErrors = 20
  nuisanceKeys = allNuisanceKeys[iTable*nErrors:min((iTable+1)*nErrors,len(allNuisanceKeys))]
  tableStr = ""
  extraCols = len(nuisanceKeys)
  tableStr += r"\begin{tabular}{|c|l|c|"+"c|"*extraCols+r"} \hline" + "\n"
  tableStr += r"\multicolumn{3}{|c|}{} & \multicolumn{"+str(extraCols)+r"}{|c|}{Systematic Error (Relative Systematic Error)} \\ \hline" + "\n"
  #tableStr += r"\begin{tabular}{|c|l|c|"+"cc|"*extraCols+r"} \hline" + "\n"
  #tableStr += r"\multicolumn{3}{|c|}{} & \multicolumn{"+str(2*extraCols)+r"}{|c|}{Systematic Error (Relative Systematic Error)} \\ \hline" + "\n"
  tableStr += r"Sample & Category"
  combNamesInTitleSet = set()
  for error in nuisanceKeys:
    foundCombName = False
    for combName in combNames:
      if combName in error:
        foundCombName = True
        if not combName in combNamesInTitleSet:
          combNamesInTitleSet.add(combName)
          combNameTitle = combNameMap[combName]
          tableStr += r" & %s" % nuisanceMap.getTitle(combNameTitle)
        break
    if not foundCombName:
      tableStr += r" & %s" % nuisanceMap.getTitle(error)
  tableStr += r"\\ \hline"+ "\n"
  for energy in energies:
    for ds in datasets:
      dsLong = ds+"Hmumu{0}_{1}".format(mass,energy)
      dsLabel = "GF"
      if "vbf" in ds:
        dsLabel = "VBF"
      tableStr += "\multirow{"+str(len(categories))+"}{*}"
      tableStr += "{%s %s} \n" % (dsLabel,energy.replace("TeV"," TeV"))
      # now on to normal stuff
      for cat in categories:
        lineStr = ""
        label = TITLEMAP[cat]
        lineStr += " & "+label + " &"
        varianceCombMap = {}
        for error in nuisanceKeys:
          foundCombName = False
          for combName in combNames:
            if combName in error:
              foundCombName = True
              if not varianceCombMap.has_key(combName):
                varianceCombMap[combName] = 0.
                lineStr += "\n    "
                lineStr +=  r" {"+combName+"}\%"
                if error != nuisanceKeys[-1]:
                  lineStr += " &"
                lineStr += "    % "+combName+"  "
              err = nuisanceMap(error,dsLong,cat,mass)
              if err:
                err = abs(abs(err)-1.)
                varianceCombMap[combName] += err**2
              break
          if not foundCombName:
            err = nuisanceMap(error,dsLong,cat,mass)
            lineStr += "\n    "
            if err:
              err = abs(abs(err)-1.)
              lineStr +=  r" {0:.0f}\%".format(err*100)
            else: 
              lineStr +=  r" $<1$\%"
            if error != nuisanceKeys[-1]:
              lineStr += " &"
            lineStr += "    % "+error+"  "
        formatArgs = {}
        for combName in varianceCombMap.keys():
            err = sqrt(varianceCombMap[combName])
            if err >= 0.01:
              err = "{0:.0f}".format(err*100)
            else:
              err = "$<1$"
            formatArgs[combName] = err
        lineStr = lineStr.format(**formatArgs)
        lineStr += "\n"+r"    \\ "+ "\n"
        tableStr += lineStr
      tableStr = tableStr[:-1] + r" \hline "+ "\n"
  tableStr += r"\end{tabular}" + "\n"
  
  print
  print tableStr
  print
