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

combNames = [
 "QCDscale",
  "pdf_",
  #"_j",
]
combNameMap = {
  "QCDscale":"QCD Scale",
  "pdf_":"PDF",
  #"_j":"Jet Energy",
  "PU":"Pileup",
  "PUID":"Pileup Jet Rejection",
}

######################################
######################################
######################################
######################################
######################################

## add back seperate PU and PUID systematics
## They were combined into CMS_eff_j in
## ad0070889faafc6b1ac7b80620e2698ad888784f
nuisanceMap.PU = {
      'gg' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : 1.0222,
          'Jet2CutsGFPass' : 1.0103,
          'Jet2CutsFailVBFGF' : 1.0112,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : 1.0222,
          'Jet2CutsGFPass' : 1.0103,
          'Jet2CutsFailVBFGF' : 1.0112,
          },
        },
      'vbf' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : 1.0207,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : 1.0207,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      'wh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      'zh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      }
nuisanceMap.PUID = {
      'gg' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0109,
          'Jet2CutsGFPass' : 1.0158,
          'Jet2CutsVBFPass' : 1.0368,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0142,
          'Jet2CutsGFPass' : 1.0174,
          'Jet2CutsVBFPass' : 1.0386,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'vbf' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0134,
          'Jet2CutsGFPass' : 1.0160,
          'Jet2CutsVBFPass' : 1.0324,
          'Jets01FailPtG10' : 1.0115,
          'Jets01PassPtG10' : 1.0108,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0155,
          'Jet2CutsGFPass' : 1.0168,
          'Jet2CutsVBFPass' : 1.0328,
          'Jets01FailPtG10' : 1.0129,
          'Jets01PassPtG10' : 1.0128,
          },
        },
      'wh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'zh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      }

######################################
######################################
######################################
######################################
######################################

data = {}

allNuisanceKeys = sorted(nuisanceMap.keys())
allNuisanceKeys.remove("QCDscale_VH")
allNuisanceKeys.remove("CMS_eff_j")
allNuisanceKeys.append("PU")
allNuisanceKeys.append("PUID")

headerStr = r"Sample & Category"
nHeaders = 2
combNamesInTitleSet = set()
for error in allNuisanceKeys:
  foundCombName = False
  for combName in combNames:
    if combName in error:
      foundCombName = True
      if not combName in combNamesInTitleSet:
        data[combName] = {}
        combNamesInTitleSet.add(combName)
        combNameTitle = combNameMap[combName]
        headerStr += r" & %s" % nuisanceMap.getTitle(combNameTitle)
        nHeaders +=1
      break
  if not foundCombName:
    data[error] = {}
    headerStr += r" & %s" % nuisanceMap.getTitle(error)
    nHeaders +=1
headerStr += r"\\ \hline"+ "\n"

tableStr = ""
extraCols = len(allNuisanceKeys)
tableStr += r"\begin{tabular}{ c l "+"c "*(nHeaders-2)+r"} \hline" + "\n"
tableStr += r"\multicolumn{3}{ c }{} & \multicolumn{"+str((nHeaders-3))+r"}{ c }{Systematic Error (Relative Systematic Error)} \\ \hline" + "\n"
tableStr += headerStr

## Finish setting up data
for key in data:
  for energy in energies:
    data[key][energy] = {}
    for ds in datasets:
      data[key][energy][ds] = {}

# Make Table 1
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
      iColumn = 2
      varianceCombMap = {}
      for error in allNuisanceKeys:
        foundCombName = False
        for combName in combNames:
          if combName in error:
            foundCombName = True
            if not varianceCombMap.has_key(combName):
              varianceCombMap[combName] = 0.
              lineStr += "\n    "
              lineStr +=  r" {"+combName+"}"
              if iColumn < nHeaders-1:
                lineStr += " &"
              else:
                lineStr += "  "
              lineStr += "    % "+combName+"  "
              iColumn += 1
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
            lineStr +=  r" {0:.0f}".format(err*100)
            data[error][energy][ds][cat] = err
          else: 
            lineStr +=  r" $<1$"
            data[error][energy][ds][cat] = None
          if iColumn < nHeaders-1:
            lineStr += " &"
          else:
            lineStr += "  "
          lineStr += "    % "+error+"  "
          iColumn += 1
      formatArgs = {}
      for combName in varianceCombMap.keys():
          err = sqrt(varianceCombMap[combName])
          if err >= 0.01:
            data[combName][energy][ds][cat] = err
            err = "{0:.0f}".format(err*100)
          else:
            data[combName][energy][ds][cat] = None
            err = "$<1$"
          formatArgs[combName] = err
      lineStr = lineStr.format(**formatArgs)
      lineStr += "\n"+r"    \\ "+ "\n"
      tableStr += lineStr
    tableStr = tableStr[:-1] + r" \hline "+ "\n"
tableStr += r"\end{tabular}" + "\n"


############################################3

tableStr2 = ""
tableStr2 += r"\begin{tabular}{ l c c c c } \hline" + "\n"
lineStr = " "*15+" &"
for energy in energies:
  for ds in datasets:
    dsLabel = "GF"
    if "vbf" in ds:
      dsLabel = "VBF"
    valStr = "%s %s" % (dsLabel,energy.replace("TeV"," TeV"))
    lineStr += " {0:9} &".format(valStr)

tableStr2 += lineStr[:-1] + " \\\\ \\hline \n"

for error in sorted(data):
  lineStr = ""
  errorTitle = nuisanceMap.getTitle(error)
  if combNameMap.has_key(error):
    errorTitle = combNameMap[error]
  lineStr += "{0:15} &".format(errorTitle)
  for energy in energies:
    for ds in datasets:
      minUnc = 1e9
      maxUnc = 0.
      for cat in categories:
        unc = data[error][energy][ds][cat]
        if unc == None:
          unc = 0.
        if unc < minUnc:
          minUnc = unc
        if unc > maxUnc:
          maxUnc = unc
      valStr = ""
      minUnc *= 100.
      maxUnc *= 100.
      if minUnc == maxUnc or (minUnc < 1. and maxUnc < 1.):
        if minUnc < 1.:
          valStr = r"$<1$"
        else:
          valStr = r"{0:.0f}".format(minUnc)
      else:
        if minUnc < 1.:
          valStr = r"$<1$-{0:.0f}".format(maxUnc)
        else:
          valStr = "{0:.0f}-{1:.0f}".format(minUnc,maxUnc)
      lineStr += " {0:9} &".format(valStr)
  tableStr2 += lineStr[:-1] + " \\\\ \n"
tableStr2 += r"\hline" + "\n"
tableStr2 += r"\end{tabular}" + "\n"

############################################3

tableStr3 = ""
tableStr3 += r"\begin{tabular}{ l c c } \hline" + "\n"
lineStr = " "*15+" &"
for ds in datasets:
  dsLabel = "GF"
  if "vbf" in ds:
    dsLabel = "VBF"
  lineStr += " {0:9} &".format(dsLabel)

tableStr3 += lineStr[:-1] + " \\\\ \\hline \n"

for error in sorted(data):
  lineStr = ""
  errorTitle = nuisanceMap.getTitle(error)
  if combNameMap.has_key(error):
    errorTitle = combNameMap[error]
  lineStr += "{0:15} &".format(errorTitle)
  for ds in datasets:
      minUnc = 1e9
      maxUnc = 0.
      for energy in energies:
        for cat in categories:
          unc = data[error][energy][ds][cat]
          if unc == None:
            unc = 0.
          if unc < minUnc:
            minUnc = unc
          if unc > maxUnc:
            maxUnc = unc
      valStr = ""
      minUnc *= 100.
      maxUnc *= 100.
      if minUnc == maxUnc or (minUnc < 1. and maxUnc < 1.):
        if minUnc < 1.:
          valStr = r"$<1$"
        else:
          valStr = r"{0:.0f}".format(minUnc)
      else:
        if minUnc < 1.:
          valStr = r"$<1$-{0:.0f}".format(maxUnc)
        else:
          valStr = "{0:.0f}-{1:.0f}".format(minUnc,maxUnc)
      lineStr += " {0:9} &".format(valStr)
  tableStr3 += lineStr[:-1] + " \\\\ \n"
tableStr3 += r"\hline" + "\n"
tableStr3 += r"\end{tabular}" + "\n"

##############################3


#print
#print tableStr
#print

print
print tableStr2
print

print
print tableStr3
print

##############################3

testStr = ""
testStr += r"""
\documentclass[12pt]{article}
 
\usepackage{fullpage}
\usepackage{multirow}
\usepackage{tabu}

\newcommand{\hmm}{\ensuremath{H\to\mu^+\mu^-}}
\newcommand{\hee}{\ensuremath{H\to e^+ e^-}}

\newcommand{\BF}{\ensuremath{\mathcal{B}}}
 
\begin{document}
\tiny
 
"""

testStr += tableStr2

testStr += "\n"
testStr += "\n"
testStr += r"\vspace{3em}"
testStr += "\n"
testStr += "\n"

testStr += tableStr3

testStr += "\n"
testStr += "\n"
testStr += r"\vspace{3em}"
testStr += "\n"
testStr += "\n"

testStr += tableStr
 
testStr += r"""

\end{document}

"""

testfile = open("test.tex",'w')
testfile.write(testStr)
testfile.close()

