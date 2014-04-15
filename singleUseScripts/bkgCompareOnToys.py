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

from multiprocessing import Pool
import itertools
from itertools import repeat as itrRepeat

from helpers import *
from makeCards import *
from xsec import *
from fitOrderChooser import makePDFBakSumExp,makePDFBakBernstein

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

class PlotBgkFits:
  def __init__(self,catName,energyStr,dataFileNames,outPrefix,pdfsToTry,nToys):
    catName = catName[0]
    self.catName = catName
    self.energyStr = energyStr
    self.nToys = nToys
    self.outPrefix = outPrefix
    randomGenerator = root.RooRandom.randomGenerator()
    randomGenerator.SetSeed(10001)
    dataTree = root.TChain()
    for i in dataFileNames:
      dataTree.Add(i+"/outtree"+catName)
    dataTree.SetCacheSize(10000000);
    dataTree.AddBranchToCache("*");

    dimuonMass = root.RooRealVar("dimuonMass","M(#mu#mu) [GeV/c^{2}]",110.,160.)
    dimuonMass.setRange("low",110,120) # Silly ranges for old fit functionality
    dimuonMass.setRange("high",130,160)
    dimuonMass.setRange("signal",120,130)
    dimuonMass.setRange("signalfit",110,140)
    dimuonMass.setBins(50)
    self.dimuonMass = dimuonMass

    # Hack to Make makePDFBakOld work
    minMassZ = 88.
    maxMassZ = 94.
    dimuonMassZ = root.RooRealVar("dimuonMass","dimuonMass",minMassZ,maxMassZ)

    ### Load data
    realData = root.RooDataSet("realData"+catName+energyStr,
                                    "realData"+catName+energyStr,
                                        dataTree,root.RooArgSet(dimuonMass)
                                      )
    nData = realData.sumEntries()
    realDataHist = realData.binnedClone()

    realDataZ = root.RooDataSet("realDataZ"+catName+energyStr,
                                    "realDataZ"+catName+energyStr,
                                        dataTree,root.RooArgSet(dimuonMassZ)
                                      )


    ### Make Bak Pdfs
    self.pdfTitles = []
    self.pdfList = []
    self.wList = []
    self.frDebugDict = {}

    self.outStr = "################################################################"
    self.outStr += "\n\n"
    self.outStr += catName + energyStr
    self.outStr += "\n\n"
    self.latexStrDetail = ""
    self.outStrDetail = "##############################################################\n\n"
    self.outStrDetail += catName + energyStr + "\n\n"
    for pdfName in pdfsToTry:
      pdfOrder = None
      pdfBaseName = pdfName
      orderMatch = re.match(r"([\d]+)(.+)",pdfBaseName)
      if orderMatch:
        pdfBaseName = orderMatch.group(2)
        pdfOrder = int(orderMatch.group(1))

      w = root.RooWorkspace("w"+pdfBaseName)
      wImport = getattr(w,"import")
      pdfFunc = globals()["makePDFBak"+pdfBaseName]
      tmpParamList,tmpNormTup,tmpDebug,tmpOrder = pdfFunc(pdfName+catName+energyStr,realData,dimuonMass,110,160,wImport,dimuonMassZ,realDataZ,order=pdfOrder)
      pdf = w.pdf("bak")
      pdf.SetName(pdfName)
      self.pdfList.append(pdf)
      pdfTitle = PDFTITLEMAP[pdfBaseName]
      if tmpOrder != None:
        pdfTitle = str(tmpOrder)+"-"+pdfTitle
      self.pdfTitles.append(pdfTitle)
      self.wList.append(w)

    binWidth = 1
    if "Jet2" in catName or "VBF" in catName:
        binWidth *= 2.5
    elif "BO" in catName:
        binWidth *= 1
    elif "BE" in catName:
        binWidth *= 2.5
    elif "OO" in catName:
        binWidth *= 2.5
    elif "OE" in catName:
        binWidth *= 2.5
    elif "EE" in catName:
        binWidth *= 2.5
    elif "FF" in catName:
        binWidth *= 2.5
    elif "CC" in catName:
        binWidth *= 2.5
    elif "BB" in catName:
        binWidth *= 1

    binning = dimuonMass.getBinning()
    xlow = binning.lowBound()
    xhigh = binning.highBound()
    dimuonMass.setBins(int((xhigh-xlow)/binWidth))

    ## Generate Data time
    pdfToGenFrom = self.pdfList[0]
    nEvents = realData.sumEntries()
    nEventsVar = root.RooFit.RooConst(nEvents)
    pdfEToGenFrom = root.RooExtendPdf(pdfToGenFrom.GetName()+"Extended","",pdfToGenFrom,nEventsVar)
    for iToy in range(nToys):
      toyData = pdfEToGenFrom.generate(root.RooArgSet(dimuonMass),root.RooFit.Extended())
      self.toyFitCompare(toyData,iToy)
    
  def toyFitCompare(self,toyData,iToy):
    dimuonMass = self.dimuonMass
    frList = []
    for pdf in self.pdfList:
      fr = pdf.fitTo(toyData, 
                         root.RooFit.Hesse(True), 
                         root.RooFit.Minos(True), # Doesn't Help, just makes it run longer
                         root.RooFit.Save(True),
                         PRINTLEVEL
                       )

      frList.append(fr)
    
    rcm = RooCompareModels(dimuonMass,toyData,
                              self.pdfList,frList,self.pdfTitles,
                              TITLEMAP[self.catName]+", Toy: {0}".format(iToy),self.energyStr,lumiDict[self.energyStr]
                              )
    rcm.draw(self.outPrefix+"_Comb_"+self.energyStr+"_"+self.catName+"_Toy{0}".format(iToy))

if __name__ == "__main__":
  canvas = root.TCanvas()
  outDir = "output/"

  nToys = 10

  pdfsToTry = ["MSSM","Bernstein","ExpMOverSq","VoigtPMm2","VoigtPExpMm2","Old","SumExp"]
  #pdfsToTry = ["MSSM","1Bernstein","2Bernstein","3Bernstein","4Bernstein","5Bernstein","1SumExp","2SumExp","3SumExp"]

  categories = []

  jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
  jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

  categoriesAll = ["BB","BO","BE","OO","OE","EE"]
  #categories += [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10BO",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10BE",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10OE",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #categories += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  categories += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  #categories += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  #categories += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

  dataDir = getDataStage2Directory()
  #dataDir = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
  #dataDir = "/cms/data/store/user/jhugon/hmumu/stage2/"

  dataFns8TeV = [
    "SingleMuRun2012Av1-22Jan2013",
    "SingleMuRun2012Bv1-22Jan2013",
    "SingleMuRun2012Cv1-22Jan2013",
    "SingleMuRun2012Dv1-22Jan2013",
    ]

  dataFns7TeV = [
    "SingleMuRun2011Av1",
    "SingleMuRun2011Bv1"
    ]
  dataFns7TeV = [dataDir+i+".root" for i in dataFns7TeV]
  dataFns8TeV = [dataDir+i+".root" for i in dataFns8TeV]

  bkgFitList = []
  for energy,dataFns in zip(["7TeV","8TeV"],[dataFns7TeV,dataFns8TeV]):
    for category in categories:
      bkgFits = PlotBgkFits(category,energy,dataFns,outDir+"bkgToyFits",pdfsToTry,nToys)
      bkgFitList.append(bkgFits)
