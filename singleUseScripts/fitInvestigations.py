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

class PDFBakMSSM(object):
  def __init__(self,dimuonMass):
    channelName = "MSSM"

    self.dimuonMass = dimuonMass
    self.bwWidth = root.RooRealVar(channelName+"_bwWidth","bwWidth",3.80102314028)
    self.bwmZ = root.RooRealVar(channelName+"_bwmZ","bwmZ",91.055064672)
    self.expParam = root.RooRealVar(channelName+"_expParam","expParam",-1e-03,-1e-01,1e-01)
    self.mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0.5,0,1)

    self.phoExpMmumu = root.RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",root.RooArgList(self.dimuonMass,self.expParam))
    self.bwExpMmumu  = root.RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(self.dimuonMass,self.bwmZ,self.bwWidth,self.expParam))
    self.pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(self.bwExpMmumu,self.phoExpMmumu),root.RooArgList(self.mixParam))
    self.pdf = self.pdfMmumu


class PDFBakMSSMPrime(object):
  def __init__(self,dimuonMass):
    channelName = "MSSMPrime"

    self.dimuonMass = dimuonMass
    self.bwWidth = root.RooRealVar(channelName+"_bwWidth","bwWidth",3.80102314028)
    self.bwmZ = root.RooRealVar(channelName+"_bwmZ","bwmZ",91.055064672)
    self.expParam = root.RooRealVar(channelName+"_expParam","expParam",-1e-03,-1e-01,1e-01)
    #self.phoCoef = root.RooRealVar(channelName+"_phoCoef","Photon Term Coef",0.5,0,1)
    #self.bwCoef = root.RooRealVar(channelName+"_bwCoef","Breit-Wigner Term Coef",0.5,0,1)
    #self.phoCoef = root.RooRealVar(channelName+"_phoCoef","Photon Term Coef",1e-10)
    #self.bwCoef = root.RooRealVar(channelName+"_bwCoef","Breit-Wigner Term Coef",0.5,0,1)
    self.phoCoef = root.RooRealVar(channelName+"_phoCoef","Photon Term Coef",0.5,0,1)
    self.bwCoef = root.RooRealVar(channelName+"_bwCoef","Breit-Wigner Term Coef",1e-10)

    #self.phoExpMmumu = root.RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",root.RooArgList(self.dimuonMass,self.expParam))
    #self.bwExpMmumu  = root.RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(self.dimuonMass,self.bwmZ,self.bwWidth,self.expParam))
    #self.phoExpMmumuE = root.RooExtendPdf("phoExpMmumuE","phoExpMmumuE",self.phoExpMmumu,self.phoCoef)
    #self.bwExpMmumuE = root.RooExtendPdf("bwExpMmumuE","bwExpMmumuE",self.bwExpMmumu,self.bwCoef)
    #self.pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(self.bwExpMmumuE,self.phoExpMmumuE))
    #self.pdf = self.pdfMmumu

    self.pdfMmumu = root.RooGenericPdf("bak","@4*exp(@0*@1)*pow(@0,-2)+@5*exp(@0*@1)/(pow(@0-@2,2)+0.25*pow(@3,2))",root.RooArgList(self.dimuonMass,self.expParam,self.bwmZ,self.bwWidth,self.phoCoef,self.bwCoef))
    self.pdf = self.pdfMmumu

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

### Load PDF from file
#pdfs =  rooArgSet2List(w.allPdfs())
#pdf = None
#for i in pdfs:
#  if "MSSM" in i.GetName():
#    pdf = i

## Make PDF here
#pdfObj = PDFBakMSSM(dimuonMass)
pdfObj = PDFBakMSSMPrime(dimuonMass)
pdf = pdfObj.pdf

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

  # Prefit
  pdf.fitTo(dh,
                 PRINTLEVEL,
                )


  log.write("\nMINOS Fit:\n")
  fr = pdf.fitTo(ds,
                 root.RooFit.Minos(True),
                 PRINTLEVEL,
                 root.RooFit.Save()
                )
  log.write(rooDebugFR(fr))
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
  log.write(rooDebugChi2(pdf,dh))

  log.write("\n")

log.close()
