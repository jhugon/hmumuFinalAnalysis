#!/usr/bin/env python


from helpers import *
from xsec import *
import makeCards
import math
import os.path
import glob
import random
import sys

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

def makePlot(filename,title,saveName,energyStr,binWidth=1.,signal=None):
  print "Making: ",title, " ", energyStr
  print "  infile: ",filename
  assert(signal)

  f = root.TFile(filename)
  tree = None
  for key in f.GetListOfKeys():
  #  key.Print()
    if key.GetClassName() == "TTree":
      tree = key.ReadObj()
      break

  fSig = root.TFile(signal)
  sigHist = None
  for key in fSig.GetListOfKeys():
  #  key.Print()
    if key.GetClassName() == "TH1F":
      sigHist = key.ReadObj()
      break
  if not sigHist:
    print("Error: sigHist not found for signal file name: "+str(signal))
    sys.exit(1)

  dimuonMass = root.RooRealVar("M","dimuonMass",110,160)
  dimuonMass.setRange("sb_low",110,120)
  dimuonMass.setRange("sb_high",130,160)
  binning = dimuonMass.getBinning()
  xlow = binning.lowBound()
  xhigh = binning.highBound()
  dimuonMass.setBins(int((xhigh-xlow)/binWidth))
  dimuonMass.SetTitle("M(#mu#mu) [GeV/c^{2}]")

  rooDataset = root.RooDataSet("data","",tree,root.RooArgSet(dimuonMass))

  channelName="blarg"
  #InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91.187, 30., 105.)
  #ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -2., 2.)

  #pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",root.RooArgList(dimuonMass,InvPolMass,ExpMass))

  InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 90., -20., 100.)
  #InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 90., 80., 90.)
  #InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 90., 50., 95.)
  #InvPolMass.setConstant();
  ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -1., 1.)
  #ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -1., 0.020)

  pdfMmumu = root.RooGenericPdf("bak","exp(@2*@0)*pow(@0-@1,-2)",root.RooArgList(dimuonMass,InvPolMass,ExpMass))

  #InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 90., -20., 150.)
  #ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -5., 5.)

  #pdfMmumu = root.RooGenericPdf("bak","exp(@2*@0)*pow(@0-@1,-2)",root.RooArgList(dimuonMass,InvPolMass,ExpMass))

  rooDataset.Print()
  print rooDataset.mean(dimuonMass)
  print rooDataset.sigma(dimuonMass)
  pdfMmumu.Print()
  fr = pdfMmumu.fitTo(rooDataset,PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("sb_low,sb_high"))
  fr.SetName("bak"+"_fitResult")
  ExpMass.Print()
  InvPolMass.Print()
  fr.Print()

  ### Signal time
  
  sigDataHist = root.RooDataHist("sigDataHist","",root.RooArgList(dimuonMass),sigHist)

  meanG1 = root.RooRealVar("MeanG1","MeanG1",125,125-10,125+3)
  meanG2 = root.RooRealVar("MeanG2","MeanG2",125,125-3,125+3)
  widthG1 = root.RooRealVar("WidthG1","WidthG1",7.049779267, 0.0,50.0)
  widthG2 = root.RooRealVar("WidthG2","WidthG2",1.830513636, 0.0, 4.0)
  mixGG = root.RooRealVar("mixGG","mixGG", 0.1140210709, 0.0,1.0)

  gaus1 = root.RooGaussian("gaus1","gaus1",dimuonMass,meanG1,widthG1)
  gaus2 = root.RooGaussian("gaus2","gaus2",dimuonMass,meanG2,widthG2)

  pdfSignal = root.RooAddPdf("pdfSignal","pdfMmumuGG",gaus1,gaus2,mixGG)
  frSig = pdfSignal.fitTo(sigDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
  frSig.SetName("sig"+"_fitResult")
  #frSig.Print()
  meanG1.setConstant(1)
  meanG2.setConstant(1)
  widthG1.setConstant(1)
  widthG2.setConstant(1)
  mixGG.setConstant(1)

  sigN = sigHist.Integral()
  #print "original sigN: ",sigN
  sigN /= 500.
  if energyStr == "7TeV":
    sigN *= 1.0
  else:
    sigN *= 1.021

  ## Like Ana A
  sigN *= 20.
  sigSF = 20.

  ### Like Ana B
  #sigN *= 100.
  #sigSF = 100.
  #if "2-Jet" in title and not ("Non-VBF" in title):
  #  sigN /= 2.
  #  sigSF /= 2.
  #  # with hack to match PASv24 (incorrect!!!)
  #  # Was norming to 1 GeV bins while plotting w.r.t. 2.5 GeV bins on y-axis
  #  #sigN /= 2.5

  #print "sigN: ",sigN
  #print "sigSF: ",sigSF
  
  legEntrySignal = "SM Higgs#times{0:.0f}".format(sigSF)

  ### Plot time

  title = "Analysis B: "+title

  rmp = RooModelPlotter(dimuonMass,pdfMmumu,rooDataset,fr,
                                    title,energyStr,lumiDict[energyStr],
                                    nSignal=sigN,signalPdf=pdfSignal,
                                    legEntrySignal=legEntrySignal,
                                    )
  #sigDataHist.plotOn(rmp.frame,root.RooFit.MarkerColor(root.kRed),root.RooFit.LineColor(root.kRed))
  rmp.draw(saveName)
  return rmp

def fitLikeMoscow(filename):
  M = root.RooRealVar("M", "M",110,160)
  c61 = root.RooRealVar("c61","c61", 0., -1., 1.);
  c63 = root.RooRealVar("c63","c63", 90, -20, 100);
  bkg = root.RooGenericPdf("bkg","bkg"," exp(c61*M)*pow(M-c63,-2) ",root.RooArgList(M,c61,c63));

  fd = root.TFile(filename)
  tree = fd.Get("BB")
  data = root.RooDataSet("data","",tree,root.RooArgSet(M))
  
  M.setRange("sb_lo",110,120);
  M.setRange("sb_hi",130,160);

  r = bkg.fitTo(data, root.RooFit.Range("sb_lo, sb_hi"), root.RooFit.Save(),root.RooFit.PrintLevel(-1)) ;
  r.Print()
  
if __name__ == "__main__":
  root.gROOT.SetBatch(True)
  rmpList = []
  outdir = "shapes/"

  #fitLikeMoscow("analysisB/Data_7TeV_2J_cl1_v2.root")
  #sys.exit(0)

  tmpRmp = makePlot("analysisB/Data_7TeV_2J_cl1_v2.root","2-Jet Non-VBF",outdir+"AnaB_2J1_7TeV","7TeV",binWidth=2.5,signal="analysisB/7TeV_2J1_Sig125_v3.root")
  rmpList.append(tmpRmp)
  tmpRmp = makePlot("analysisB/Data_7TeV_2J_cl2_v2.root","2-Jet VBF",outdir+"AnaB_2J2_7TeV","7TeV",binWidth=2.5,signal="analysisB/7TeV_2J2_Sig125_v3.root")
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_7TeV_1J_v2.root","1-Jet",outdir+"AnaB_1J_7TeV","7TeV",signal="analysisB/7TeV_1J_Sig125_v3.root")
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_7TeV_0J1_v2.root","0-Jet Loose",outdir+"AnaB_0J1_7TeV","7TeV",signal="analysisB/7TeV_0J1_Sig125_v3.root")
  rmpList.append(tmpRmp)
  
  tmpRmp = makePlot("analysisB/Data_7TeV_0J2_v2.root","0-Jet Tight",outdir+"AnaB_0J2_7TeV","7TeV",signal="analysisB/7TeV_0J2_Sig125_v3.root")
  rmpList.append(tmpRmp)

  #################################
  
  tmpRmp = makePlot("analysisB/Data_2J0_LoosePU_ReRecoMSF.root","2-Jet Non-VBF",outdir+"AnaB_2J0_8TeV","8TeV",binWidth=2.5,signal="analysisB/8TeV_2J1_Sig125.root")
  rmpList.append(tmpRmp)
  tmpRmp = makePlot("analysisB/Data_2J1_LoosePU_ReRecoMSF.root","2-Jet VBF Loose",outdir+"AnaB_2J1_8TeV","8TeV",binWidth=2.5,signal="analysisB/8TeV_2J2_Sig125.root")
  rmpList.append(tmpRmp)
  tmpRmp = makePlot("analysisB/Data_2J2_LoosePU_ReRecoMSF.root","2-Jet VBF Tight",outdir+"AnaB_2J2_8TeV","8TeV",binWidth=2.5,signal="analysisB/8TeV_2J3_Sig125.root")
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_1J_LoosePU_ReRecoMSF.root","1-Jet",outdir+"AnaB_1J_8TeV","8TeV",signal="analysisB/8TeV_1J_Sig125.root")
  rmpList.append(tmpRmp)

  tmpRmp = makePlot("analysisB/Data_0J1_LoosePU_ReRecoMSF.root","0-Jet Loose",outdir+"AnaB_0J1_8TeV","8TeV",signal="analysisB/8TeV_0J1_Sig125.root")
  rmpList.append(tmpRmp)
  
  tmpRmp = makePlot("analysisB/Data_0J2_LoosePU_ReRecoMSF.root","0-Jet Tight",outdir+"AnaB_0J2_8TeV","8TeV",signal="analysisB/8TeV_0J2_Sig125.root")
  rmpList.append(tmpRmp)
  
