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

from singleUseScripts.biasPklToMu import getSMSigCounts

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
mssmPdf = None
sigPdf = None
for i in pdfs:
  name = i.GetName()
  if "MSSM" in name:
    mssmPdf = i
  elif "sigPDF" in name and name[-1] == "E" :
    sigPdf = i

mssmPdf.SetName("mssmBak")

data = w.data("realDataJets01PassPtG10BB8TeV")

wNew2 = root.RooWorkspace("wNew2")
wNew2Import = getattr(wNew2,"import")

trashParamList, trashBakNormTup, trashDebug, trashOrder = makePDFBakBernsteinProd("BernPDF_BB",data,dimuonMass,110,160,wNew2Import,order=4)
bernPdf = wNew2.pdf("bak")
bernPdf.SetName("bernBak")

nSig = w.var('nSig')
nBak = root.RooRealVar("nBak","N_{bkg}",data.sumEntries(),nSig.getMax()*2,data.sumEntries()*2)

mssmPdfE = root.RooExtendPdf("mssmBakE","Extended Background PDF",mssmPdf,nBak)
bernPdfE = root.RooExtendPdf("bernBakE","Extended Background PDF",bernPdf,nBak)
mssmPdfSB = root.RooAddPdf("mssmPdfSB","S+B PDF",root.RooArgList(mssmPdfE,sigPdf))
bernPdfSB = root.RooAddPdf("bernPdfSB","S+B PDF",root.RooArgList(bernPdfE,sigPdf))

log = open("output/fitInvestigate3.log",'w')

nSigFitVals = []
rmpList = []
ds = data
dh = ds.binnedClone()
ds = dh

nTimesSMToInject = 5.2
nSigToInject = nTimesSMToInject*getSMSigCounts("Jets01PassPtG10BB",125)
nSigToInject = int(nSigToInject)
assert(nSigToInject*2 < nSig.getMax())
dataSigInject = sigPdf.generate(root.RooArgSet(dimuonMass),nSigToInject)
print "N Events Signal Expected: ", nSigToInject
print "N Events Signal Injected: ", dataSigInject.sumEntries()
dataSigInject.append(data)
print "N Events data+signal: ", dataSigInject.sumEntries()
print "N Events data: ", data.sumEntries()
nTimesSMToInjectStr=("{0:.2g}".format(nTimesSMToInject)).replace(".","p")

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nMSSM S+B Fit:\n")
fr = mssmPdfSB.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,mssmPdfSB,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="S+B Fit",
                      caption3="No Signal Injected",
                      pdfDotLineName="mssmBak",
                      legEntryModel="S+MSSM"
                      )
rmp.draw("output/fitInvestigate3_mssm_SBFit")
rmpList.append(rmp)
log.write(rooDebugChi2(mssmPdfSB,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nMSSM B Fit:\n")
fr = mssmPdf.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,mssmPdf,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="B-Only Fit",
                      caption3="No Signal Injected",
                      legEntryModel="MSSM"
                      )
rmp.draw("output/fitInvestigate3_mssm_BFit")
rmpList.append(rmp)
log.write(rooDebugChi2(mssmPdf,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nBernstein 4 S+B Fit (All Free Coefs):\n")
fr = bernPdfSB.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdfSB,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="S+B Fit",
                      caption3="No Signal Injected",
                      caption4="p0 Free",
                      pdfDotLineName="bernBak",
                      legEntryModel="S+4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFreeP0_SBFit")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdfSB,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nBernstein 4 B Fit (All Free Coefs):\n")
fr = bernPdf.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdf,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="B-Only Fit",
                      caption3="No Signal Injected",
                      caption4="p0 Free",
                      legEntryModel="4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFreeP0_BFit")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdf,dh))

################################################3
for param in rooArgSet2List(bernPdf.getParameters(ds)):
  if param.GetName() == "BernPDF_BB_B0":
    param.setVal(1e-6)
    param.setConstant(True)
log.write("\n"+'#'*30+"\n")
log.write("\nBernstein 4 S+B Fit (Fix p0):\n")
fr = bernPdfSB.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdfSB,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="S+B Fit",
                      caption3="No Signal Injected",
                      caption4="p0 = 10^{-6}",
                      pdfDotLineName="bernBak",
                      legEntryModel="S+4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFixP0_SBFit")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdfSB,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nBernstein 4 B Fit (Fix p0):\n")
fr = bernPdf.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdf,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="B-Only Fit",
                      caption3="No Signal Injected",
                      caption4="p0 = 10^{-6}",
                      legEntryModel="4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFixP0_BFit")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdf,dh))

#######################################################
#######################################################
#######################################################
#######################################################
#######################################################

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\n"+'#'*30+"\n")
log.write("\n"+'#'*30+"\n")
log.write("\n"+'#'*30+"\n")
log.write("Injected {0}*SM Higgs Signal, {1} events\n".format(nTimesSMToInject,nSigToInject))
ds = dataSigInject
dh = ds.binnedClone()
ds = dh

for param in rooArgSet2List(bernPdf.getParameters(ds)):
  param.setConstant(False)

log.write("\nSignal Injected: MSSM S+B Fit:\n")
fr = mssmPdfSB.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,mssmPdfSB,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="S+B Fit",
                      caption3="{0:.2g} #times SM Injected".format(nTimesSMToInject),
                      pdfDotLineName="mssmBak",
                      legEntryModel="S+MSSM"
                      )
rmp.draw("output/fitInvestigate3_mssm_SBFit_"+nTimesSMToInjectStr+"xInjected")
rmpList.append(rmp)
log.write(rooDebugChi2(mssmPdfSB,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nSignal Injected: MSSM B Fit:\n")
fr = mssmPdf.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,mssmPdf,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="B-Only Fit",
                      caption3="{0:.2g} #times SM Injected".format(nTimesSMToInject),
                      legEntryModel="MSSM"
                      )
rmp.draw("output/fitInvestigate3_mssm_BFit_"+nTimesSMToInjectStr+"xInjected")
rmpList.append(rmp)
log.write(rooDebugChi2(mssmPdf,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nSignal Injected: Bernstein 4 S+B Fit (All Free Coefs):\n")
fr = bernPdfSB.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdfSB,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="S+B Fit",
                      caption3="{0:.2g} #times SM Injected".format(nTimesSMToInject),
                      caption4="p0 Free",
                      pdfDotLineName="bernBak",
                      legEntryModel="S+4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFreeP0_SBFit_"+nTimesSMToInjectStr+"xInjected")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdfSB,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nSignal Injected: Bernstein 4 B Fit (All Free Coefs):\n")
fr = bernPdf.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdf,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="B-Only Fit",
                      caption3="{0:.2g} #times SM Injected".format(nTimesSMToInject),
                      caption4="p0 Free",
                      legEntryModel="4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFreeP0_BFit_"+nTimesSMToInjectStr+"xInjected")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdf,dh))

################################################3
for param in rooArgSet2List(bernPdf.getParameters(ds)):
  if param.GetName() == "BernPDF_BB_B0":
    param.setVal(1e-6)
    param.setConstant(True)
log.write("\n"+'#'*30+"\n")
log.write("\nSignal Injected: Bernstein 4 S+B Fit (Fix p0):\n")
fr = bernPdfSB.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdfSB,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="S+B Fit",
                      caption3="{0:.2g} #times SM Injected".format(nTimesSMToInject),
                      caption4="p0 = 10^{-6}",
                      pdfDotLineName="bernBak",
                      legEntryModel="S+4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFixP0_SBFit_"+nTimesSMToInjectStr+"xInjected")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdfSB,dh))

################################################3
log.write("\n"+'#'*30+"\n")
log.write("\nSignal Injected: Bernstein 4 B Fit (Fix p0):\n")
fr = bernPdf.fitTo(ds,
               root.RooFit.Minos(True),
               PRINTLEVEL,
               root.RooFit.Save()
              )
nSigSB = nSig.getVal()
log.write(rooDebugFR(fr))
rmp = RooModelPlotter(dimuonMass,bernPdf,ds,fr,
                      "0,1-Jet Tight BB","8TeV",lumiDict["8TeV"],
                      caption2="B-Only Fit",
                      caption3="{0:.2g} #times SM Injected".format(nTimesSMToInject),
                      caption4="p0 = 10^{-6}",
                      legEntryModel="4Bern"
                      )
rmp.draw("output/fitInvestigate3_bernFixP0_BFit_"+nTimesSMToInjectStr+"xInjected")
rmpList.append(rmp)
log.write(rooDebugChi2(bernPdf,dh))

log.close()
#print "STILL USING BINNED HISTS FOR EVERYTHING!!!\n"*5
