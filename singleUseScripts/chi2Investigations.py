#! /usr/bin/env python

from ROOT import gSystem

import singleHelpers
import datetime
import sys
import os
import re
import math
import cPickle
import ROOT as root
root.gSystem.Load('libRooFit')
root.gROOT.SetBatch(True)
import scipy.stats

from helpers import *
from makeCards import *
from xsec import *
from fitOrderChooser import makePDFBakSumExp

from numpy import mean, median, corrcoef, percentile
from numpy import std as stddev

root.gErrorIgnoreLevel = root.kWarning
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(2) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT
PRELIMINARYSTRING="CMS Internal"

def rooArgSetToList(a):
  result = []
  iter = a.createIterator()
  item = iter.Next()
  while item:
    result.append(item)
    item = iter.Next()
  return result

def getPDFList(w):
  pdfs = w.allPdfs()
  return rooArgSetToList(pdfs)

def debugFR(fr):
  result = ""
  result += "{0:20}: {1:10.3g}\n".format("Status",fr.status())
  result += "{0:20}: {1:10.3g}\n".format("EDM",fr.edm())
  result += "{0:20}: {1:10.3g}\n".format("NLL",fr.minNll())
  result += "{0:20}: {1:10.3g}\n".format("Cov. Matrix Qual",fr.edm())
  result += "{0:20}: {1:10}\n".format("N Invalid NLL",fr.numInvalidNLL())
  result += "Status History:\n"
  for i in range(fr.numStatusHistory()):
    result += "  {1} ({0})\n".format(fr.statusCodeHistory(i),fr.statusLabelHistory(i))

  rowString = ""
  rowString += "{0:30}".format("Variable Name")
  rowString += "{0:10} +/- {1:<10}  ".format("Fit Value","HESSE Err")
  rowString += " {0:<10} {1:<10}  ".format("High Err","Low Err")

  rowString += "{0:20}".format("[{0}]".format("Variable Limits"))
  result += rowString + "\n"

  fpf_s = fr.floatParsFinal()
  for i in range(fpf_s.getSize()):
    nuis_s = fpf_s.at(i)
    name   = nuis_s.GetName();
    # .getVal() .getError() .getMin() .getMax()

    rowString = ""
    rowString += "{0:30}".format(re.sub(r".*TeV_","",name))
    rowString += "{0:10.3g} +/- {1:<10.3g}  ".format(nuis_s.getVal(),nuis_s.getError())
    rowString += "+{0:<10.3g} {1:<10.3g}  ".format(nuis_s.getErrorHi(),nuis_s.getErrorLo())

    rowString += "{0:20}".format("[{0:.1g},{1:.1g}]".format(nuis_s.getMin(),nuis_s.getMax()))
    result += rowString + "\n"
  constPars = fr.constPars()
  for i in range(constPars.getSize()):
    nuis_s = constPars.at(i)
    name   = nuis_s.GetName();
    # .getVal() .getError() .getMin() .getMax()

    rowString = ""
    rowString += "{0:30}".format(re.sub(r".*TeV_","",name))
    rowString += "{0:10.3g}  FIXED  ".format(nuis_s.getVal())

    rowString += "{0:20}".format("[{0:.3g},{1:.3g}]".format(nuis_s.getMin(),nuis_s.getMax()))
    result += rowString + "\n"
  return result 

def debugChi2(pdf,data):
  result = ""
  observables = rooArgSetToList(pdf.getObservables(data))
  assert(len(observables)==1)
  binning = observables[0].getBinning()
  nBins = binning.numBins()
  pdfParams = rooArgSetToList(pdf.getParameters(data))
  nparams = len(pdfParams)
  nparamsFree = 0
  for p in pdfParams:
    if not p.isConstant():
      nparamsFree += 1
  ndf = nBins - nparamsFree
  errsForChi2 = root.RooAbsData.Poisson
  errsForChi2Label = "Errors Poisson"
  chi2Var = pdf.createChi2(data,root.RooFit.DataError(errsForChi2))
  chi2 = chi2Var.getVal()
  chi2Prob = scipy.stats.chi2.sf(chi2,ndf)
  normChi2 = chi2/ndf
  result += "{0}:\n".format("Chi2 Info Using "+errsForChi2Label+" "+str(errsForChi2))
  result += "{0:20}: {1:<10.3g}\n".format("chi2",chi2)
  result += "{0:20}: {1:<10.3g}\n".format("chi2/ndf",normChi2)
  result += "{0:20}: {1:<10.3g}\n".format("chi2 Prob",chi2Prob)
  result += "{0:20}: {1:<10.3g}\n".format("ndf",ndf)
  result += "{0:20}: {1:<10.3g}\n".format("nBins",nBins)
  result += "{0:20}: {1:<10.3g}\n".format("nParams Free",nparamsFree)
  result += "{0:20}: {1:<10.3g}\n".format("nParams",nparams)
  return result, chi2Prob

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

inputFile = root.TFile("output/debug_RooFit_Bernstein_Jets01PassPtG10BB_8TeV_job0.root")

w = inputFile.Get("debugWorkspace")
#w.Print()

pdfs =  getPDFList(w)
pdf = None
for i in pdfs:
  if "MSSM" in i.GetName():
    pdf = i

dsFromWS = w.data("toyDataJets01PassPtG10BB8TeV0")
nEvents = int(dsFromWS.sumEntries())

dimuonMass = w.var("dimuonMass")

log = open("output/chi2Investigate.log",'w')

nToys = 1000

probHistBinnedGen = root.TH1F("probHistBinnedGen","",5,0,1)
probHistBinnedGen.GetYaxis().SetRangeUser(0.15*nToys,0.3*nToys)
probHistBinnedGen.Sumw2()
setHistTitles(probHistBinnedGen,"#chi^{2} Probability of Fit","N Toys")

probHistUnBinnedGen = root.TH1F("probHistUnBinnedGen","",5,0,1)
probHistUnBinnedGen.GetYaxis().SetRangeUser(0.15*nToys,0.3*nToys)
probHistUnBinnedGen.Sumw2()
setHistTitles(probHistUnBinnedGen,"#chi^{2} Probability of Fit","N Toys")

for iToy in range(nToys):
  toyData = pdf.generateBinned(root.RooArgSet(dimuonMass),nEvents)
  dh = toyData

  fr = pdf.fitTo(dh,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  chi2Str, chi2Prob = debugChi2(pdf,dh)
  probHistBinnedGen.Fill(chi2Prob)
  log.write(chi2Str)

for iToy in range(nToys):
  toyData = pdf.generate(root.RooArgSet(dimuonMass),nEvents)
  dh = toyData.binnedClone()

  fr = pdf.fitTo(dh,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  chi2Str, chi2Prob = debugChi2(pdf,dh)
  probHistUnBinnedGen.Fill(chi2Prob)
  log.write(chi2Str)

canvas = root.TCanvas()
probHistBinnedGen.Draw()
drawStandardCaptions(canvas,"Toys Generated Binned","#chi^{2} Calculated After Binned ML Fit","{0} Toys".format(nToys))
saveAs(canvas,"output/chi2Investigate_chi2ProbGenBinned")

canvas = root.TCanvas()
probHistUnBinnedGen.Draw()
drawStandardCaptions(canvas,"Toys Generated Un-Binned","#chi^{2} Calculated After Binned ML Fit","{0} Toys".format(nToys))
saveAs(canvas,"output/chi2Investigate_chi2ProbGenNonBinned")

log.close()
