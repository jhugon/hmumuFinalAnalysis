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

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

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
  #for errsForChi2, errsForChi2Label in zip([root.RooAbsData.Expected,root.RooAbsData.SumW2,root.RooAbsData.Poisson],["Errors from PDF","Errors Data Weights^2","Errors Poisson"]):
  for errsForChi2, errsForChi2Label in zip([getattr(root.RooAbsData,"None"),root.RooAbsData.Auto,root.RooAbsData.SumW2,root.RooAbsData.Poisson],["Errors None","Errors Auto","Errors Data Weights^2","Errors Poisson"]):
    chi2Var = pdf.createChi2(data,root.RooFit.DataError(errsForChi2))
    chi2 = chi2Var.getVal()
    normChi2 = chi2/ndf
    result += "{0}:\n".format("Chi2 Info Using "+errsForChi2Label+" "+str(errsForChi2))
    result += "{0:20}: {1:<10.3g}\n".format("chi2",chi2)
    result += "{0:20}: {1:<10.3g}\n".format("chi2/ndf",normChi2)
    result += "{0:20}: {1:<10.3g}\n".format("ndf",ndf)
    result += "{0:20}: {1:<10.3g}\n".format("nBins",nBins)
    result += "{0:20}: {1:<10.3g}\n".format("nParams Free",nparamsFree)
    result += "{0:20}: {1:<10.3g}\n".format("nParams",nparams)
  return result

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

dimuonMass = w.var("dimuonMass")

nToys = 10
toyDataSets = []
for i in range(nToys):
  toyDataSets.append(w.data("toyDataJets01PassPtG10BB8TeV"+str(i)))

log = open("output/fitInvestigate.log",'w')

rmpList = []
for iToy,ds in enumerate(toyDataSets):
  if iToy > 0:
    break
  dh = ds.binnedClone()
  log.write("#####################################################\n")
  log.write("Toy: {0}\n".format(iToy))

  log.write("\nDefault Fit:\n")
  fr = pdf.fitTo(ds,
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="Default Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_default")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="Default Fit; New Band",
                        doLinearErrs=False,
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_default_sampleErrBands")
  rmpList.append(rmp)

  log.write("\nMINOS Fit:\n")
  fr = pdf.fitTo(ds,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_minos")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Fit; New Band",
                        doLinearErrs=False,
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_minos_sampleErrBand")
  rmpList.append(rmp)

  log.write("\nOldMinuit MINOS Fit:\n")
  fr = pdf.fitTo(ds,
                 root.RooFit.Minimizer("OldMinuit","migrad"),
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="Old MINUIT & MINOS Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_oldminuitminos")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="Old MINUIT & MINOS Fit; New Band",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_oldminuitminos_sampleErrBand")
  rmpList.append(rmp)


  log.write("\nMinuit2 MINOS Fit:\n")
  fr = pdf.fitTo(ds,
                 root.RooFit.Minimizer("Minuit2","migrad"),
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINUIT2 & MINOS Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_minuit2minos")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdf,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINUIT2 & MINOS Fit; New Band",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_minuit2minos_sampleErrBand")
  rmpList.append(rmp)

  log.write("\nDefault Hist Fit:\n")
  pdf.Print()
  dh.Print()
  fr = pdf.fitTo(dh,
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdf,dh,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="Default Hist Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_defaultHist")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdf,dh,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="Default Hist Fit; New Band",
                        doLinearErrs=False,
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_defaultHist_sampleErrBand")
  rmpList.append(rmp)
  log.write(debugChi2(pdf,dh))

  log.write("\nMINOS Hist Fit:\n")
  fr = pdf.fitTo(dh,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(debugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdf,dh,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Hist Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_minosHist")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdf,dh,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Hist Fit; New Band",
                        doLinearErrs=False,
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate_Toy"+str(iToy)+"_minosHist_sampleErrBand")
  rmpList.append(rmp)
  log.write(debugChi2(pdf,dh))

  log.write("\n")

log.close()
