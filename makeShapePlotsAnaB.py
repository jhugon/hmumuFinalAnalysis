#!/usr/bin/env python


from helpers import *
from xsec import *
import makeCards
import math
import os.path
import glob
import random

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

def makePlot(filename,title,saveName,energyStr,binWidth=1.):
  root.gROOT.SetBatch(True)

  f = root.TFile(filename)
  tree = None
  for key in f.GetListOfKeys():
    key.Print()
    if key.GetClassName() == "TTree":
      tree = key.ReadObj()
      break

  dimuonMass = root.RooRealVar("M","dimuonMass",110,160)
  binning = dimuonMass.getBinning()
  xlow = binning.lowBound()
  xhigh = binning.highBound()
  dimuonMass.setBins(int((xhigh-xlow)/binWidth))
  dimuonMass.SetTitle("M(#mu#mu) [GeV/c^{2}]")


  rooDataset = root.RooDataSet("data","",tree,root.RooArgSet(dimuonMass))

  channelName="blarg"
  InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91.187, 30., 105.)
  ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -2., 2.)

  pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",root.RooArgList(dimuonMass,InvPolMass,ExpMass))

  fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
  fr.SetName("bak"+"_fitResult")
  fr.Print()

  rmp = RooModelPlotter(dimuonMass,pdfMmumu,rooDataset,fr,
                                    title,energyStr,lumiDict[energyStr],
                                    caption2="Analysis B"
                                    )
  rmp.draw(saveName)
  return rmp
  
if __name__ == "__main__":
  rmpList = []
  outdir = "shapes/"

  tmpRmp = makePlot("analysisB/Data_7TeV_2J_cl1_v2.root","2-Jet Non-VBF",outdir+"AnaB_2J1_7TeV","7TeV",binWidth=2.5)
  rmpList.append(tmpRmp)
  tmpRmp = makePlot("analysisB/Data_7TeV_2J_cl2_v2.root","2-Jet VBF",outdir+"AnaB_2J2_7TeV","7TeV",binWidth=2.5)
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_7TeV_1J_v2.root","1-Jet",outdir+"AnaB_1J_7TeV","7TeV")
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_7TeV_0J1_v2.root","0-Jet Loose",outdir+"AnaB_0J1_7TeV","7TeV")
  rmpList.append(tmpRmp)
  
  tmpRmp = makePlot("analysisB/Data_7TeV_0J2_v2.root","0-Jet Tight",outdir+"AnaB_0J2_7TeV","7TeV")
  rmpList.append(tmpRmp)

  #################################
  
  tmpRmp = makePlot("analysisB/Data_2J0_LoosePU_ReRecoMSF.root","2-Jet Non-VBF",outdir+"AnaB_2J0_8TeV","8TeV",binWidth=2.5)
  rmpList.append(tmpRmp)
  tmpRmp = makePlot("analysisB/Data_2J1_LoosePU_ReRecoMSF.root","2-Jet VBF Loose",outdir+"AnaB_2J1_8TeV","8TeV",binWidth=2.5)
  rmpList.append(tmpRmp)
  tmpRmp = makePlot("analysisB/Data_2J2_LoosePU_ReRecoMSF.root","2-Jet VBF Tight",outdir+"AnaB_2J2_8TeV","8TeV",binWidth=2.5)
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_1J_LoosePU_ReRecoMSF.root","1-Jet",outdir+"AnaB_1J_8TeV","8TeV")
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_0J1_LoosePU_ReRecoMSF.root","0-Jet Loose",outdir+"AnaB_0J1_8TeV","8TeV")
  rmpList.append(tmpRmp)
  
  tmpRmp = makePlot("analysisB/Data_0J2_LoosePU_ReRecoMSF.root","0-Jet Tight",outdir+"AnaB_0J2_8TeV","8TeV")
  rmpList.append(tmpRmp)
  
