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

dimuonMass = w.var("dimuonMass")

## Load PDF from file
pdfs =  rooArgSet2List(w.allPdfs())
pdf = None
sigPdf = None
for i in pdfs:
  name = i.GetName()
  if "MSSM" in name:
    pdf = i
  elif "sigPDF" in name and name[-1] == "E" :
    sigPdf = i
nSig = w.var('nSig')

nToys = 10
toyDataSets = []
for i in range(nToys):
  toyDataSets.append(w.data("toyDataJets01PassPtG10BB8TeV"+str(i)))

nBak = root.RooRealVar("nBak","N_{bkg}",toyDataSets[0].sumEntries(),nSig.getMax()*2,toyDataSets[0].sumEntries()*2)
pdfE = root.RooExtendPdf("bakE","Extended Background PDF",pdf,nBak)
pdfSB = root.RooAddPdf("pdfSB","S+B PDF",root.RooArgList(pdfE,sigPdf))

bakParamsFixed = {}

bakParams = rooArgSet2List(pdf.getParameters(toyDataSets[0]))
for i in bakParams:
  name = i.GetName()
  bakParamsFixed[name] = i.isConstant()
print bakParamsFixed

log = open("output/fitInvestigate2.log",'w')

nSigFitVals = []
rmpList = []
for iToy,ds in enumerate(toyDataSets):

  dh = ds.binnedClone()
  log.write("#####################################################\n")
  log.write("Toy: {0}\n".format(iToy))

  # Prefit
  pdfSB.fitTo(dh,
                 PRINTLEVEL,
                )


  log.write("\nS+B Fit:\n")
  fr = pdfSB.fitTo(ds,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  nSigSB = nSig.getVal()
  log.write(rooDebugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdfSB,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate2_Toy"+str(iToy)+"_SBFit")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdfSB,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Fit; New Band",
                        doLinearErrs=False,
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate2_Toy"+str(iToy)+"_SBFit_sampleErrBand")
  rmpList.append(rmp)
  log.write(rooDebugChi2(pdfSB,dh))

  log.write("\nB Fit, fix B params, then S+B Fit:\n")
  frB = pdf.fitTo(ds,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(rooDebugFR(frB))
  for i in rooArgSet2List(pdf.getParameters(ds)):
    i.setConstant(True)
  fr = pdfSB.fitTo(ds,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(rooDebugFR(fr))
  rmp = RooModelPlotter(dimuonMass,pdfSB,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Fit",
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate2_Toy"+str(iToy)+"_SBFitFixedB")
  rmpList.append(rmp)
  rmp = RooModelPlotter(dimuonMass,pdfSB,ds,fr,
                        "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                        caption2="MINOS Fit; New Band",
                        doLinearErrs=False,
                        legEntryData="Toy Data", legEntryModel="MSSM Bkg. Model"
                        )
  rmp.draw("output/fitInvestigate2_Toy"+str(iToy)+"_SBFitFixedB_sampleErrBand")
  rmpList.append(rmp)
  log.write(rooDebugChi2(pdfSB,dh))

  log.write("\n")
  for i in rooArgSet2List(pdf.getParameters(ds)):
    i.setConstant(bakParamsFixed[i.GetName()])
  nSigFitVals.append((nSigSB,nSig.getVal()))

log.write("\n\n##NSig fit values:\n{0:10} {1:10}\n".format("S+B Fit","B Fit->S+B"))
for sb,bfixsb in nSigFitVals:
  log.write("{0:10.1f} {1:10.1f}\n".format(sb,bfixsb))

log.close()
