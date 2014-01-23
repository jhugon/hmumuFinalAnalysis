#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root
root.gROOT.SetBatch(True)

import singleHelpers
from helpers import *
from xsec import *

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

"""
  const int c_iLowerFitBorder = 100;
  const int c_iUpperFitBorder = 300;

  TTree* p_trData = (TTree*)p_flData->Get("OutputTree");

  RooRealVar My_Minv("My_Minv","",c_iLowerFitBorder,c_iUpperFitBorder);
  
  RooDataSet  p_rdsData("data_obs","data_obs",RooArgSet(My_Minv),Import(*p_trData));
  
  double c_dMZ0;
  double c_dWidthZ0;
 
  if ( p_iCategory == 1 ) {
    c_dMZ0     = 90.95;
    c_dWidthZ0 =  3.97;
  } 
  if ( p_iCategory == 2 ) {
    c_dMZ0     = 90.91;
    c_dWidthZ0 =  4.01;
  } 
  if ( p_iCategory == 3 ) {
    c_dMZ0     = 90.950;
    c_dWidthZ0 =  3.898;
  } 

  RooRealVar p_rrvMZ0("p_rrvMZ0","mass of Z0",c_dMZ0);
  RooRealVar p_rrvWidthZ0("widthZ0","width of Z0",c_dWidthZ0);
  RooRealVar p_rrvLambda("lambda","p.d.f. contribution",-1e-03,-1e-01,-1e-04);
  RooRealVar p_rrvAlpha("alpha","strength of the Z Breit Wigner",0.5,0.,1.);

  RooGenericPdf p_rgpRelBW("relBW","exp(lambda*My_Minv)*widthZ0/(((My_Minv-p_rrvMZ0)*(My_Minv-p_rrvMZ0) + widthZ0*widthZ0*0.25))",RooArgSet( My_Minv, p_rrvLambda, p_rrvAlpha, p_rrvWidthZ0, p_rrvMZ0));

  RooGenericPdf p_rgpPhoton("photon","exp(lambda*My_Minv)/pow(My_Minv,2)",RooArgSet(My_Minv,p_rrvLambda));

  RooAddPdf p_rapBackground("background","relBW + photon",RooArgList(p_rgpRelBW,p_rgpPhoton),p_rrvAlpha);
 
  p_rrvMZ0.setConstant(kTRUE);
  p_rrvWidthZ0.setConstant(kTRUE);

  p_rapSignalBkg.fitTo(p_rdsData);

"""

canvas = root.TCanvas()

massVeryLowRange = [80,95]
massLowRange = [110,120]
massHighRange = [130,160]

maxMass = massHighRange[1]
minMass = massVeryLowRange[0]
mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
mMuMu.setRange("z",88,94)
mMuMu.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
mMuMu.setRange("low",massLowRange[0],massLowRange[1])
mMuMu.setRange("high",massHighRange[0],massHighRange[1])
mMuMu.setRange("signal",massLowRange[1],massHighRange[0])
mMuMu.setRange("lowsighigh",massLowRange[0],massHighRange[1])

bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
bwSig = root.RooRealVar("bwSig","bwSig",0.0,30.0)
expLambda = root.RooRealVar("expLambda","expLambda",-1e-03,-1,1.)
mixParam = root.RooRealVar("mixParam","mixParam",0.99,0,1)

#Just for fitting Z-peak
bwMmumu  = root.RooGenericPdf("bwMmumu","(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(mMuMu,bwmZ,bwSig))
#Real PDFs
phoExpMmumu = root.RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",root.RooArgList(mMuMu,expLambda))
bwExpMmumu  = root.RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(mMuMu,bwmZ,bwSig,expLambda))


pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(bwExpMmumu,phoExpMmumu),root.RooArgList(mixParam))
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(bwMmumu,phoMmumu),root.RooArgList(mixParam))
#pdfMmumu = phoExpMmumu
#pdfMmumu = bwExpMmumu

#####################################################################

filename = getDataStage2Directory()+"/DYJetsToLL_8TeV.root"
#filename = getDataStage2Directory()+"/ttbar_8TeV.root"
print filename
f = root.TFile(filename)
samples = ['ttbar_8TeV','DYJetsToLL_8TeV']
samples += ['ggHmumu125_8TeV','vbfHmumu125_8TeV']
scaleFactors = {}
#scaleFactors['ggHmumu125_8TeV'] = 100.
#scaleFactors['ttbar_8TeV'] = 10.
mDiMu = None
files = []
hstack = root.THStack("hstack","")
for sample in samples:
  f = root.TFile(getDataStage2Directory()+"/"+sample+".root")
  #mDiMuTmp = f.Get("Jets01FailPtG10BE/mDiMu")
  mDiMuTmp = f.Get("Jets01PassPtG10BB/mDiMu")
  #mDiMuTmp = f.Get("Jet2CutsGFPass/mDiMu")
  mDiMuTmp.Rebin(2)
  mDiMuTmp.Scale(xsec[sample]*1000.*lumiDict['8TeV']/nEventsMap[sample])
  mDiMuTmp.SetFillColor(colors[re.sub(r"_.*","",sample)])
  if scaleFactors.has_key(sample):
    mDiMuTmp.Scale(scaleFactors[sample])
  hstack.Add(mDiMuTmp)
  if mDiMu == None:
    mDiMu = mDiMuTmp.Clone("mDiMuHist")
    print "creating hist from ",sample
  else:
    mDiMu.Add(mDiMuTmp)
    print "adding hist from ",sample
  files.append(f)

#####################################################################

mMuMuRooDataHist = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

bwMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
bwmZ.setConstant(True)
bwSig.setConstant(True)

#plotMmumuZ = mMuMu.frame(root.RooFit.Range("z"))
#mMuMuRooDataHist.plotOn(plotMmumuZ)
#bwMmumu.plotOn(plotMmumuZ,root.RooFit.Range("z"))
##axisHistZ = TH2F("axisHistZ","",1,88,94,1,0,2000)
##axisHistZ.Draw()
##plotMmumuZ.Draw("same")
#plotMmumuZ.Draw()
#yay = raw_input("Press Enter to continue...")
#canvas.Clear()

#fr = pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("lowsighigh"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("lowsighigh"),root.RooFit.SumW2Error(False),PRINTLEVEL)
chi2 = pdfMmumu.createChi2(mMuMuRooDataHist,root.RooFit.Range("lowsighigh"))

plotMmumu = mMuMu.frame(root.RooFit.Range("lowsighigh"))

#############################################################################

mMuMu.Print()
plotMmumu.Print()
mMuMuRooDataHist.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("lowsighigh"))
#pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("phoExpMmumu"),root.RooFit.LineStyle(2),root.RooFit.LineColor(root.kGreen+1),root.RooFit.Range("lowsighigh"))
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("bwExpMmumu"),root.RooFit.LineStyle(2),root.RooFit.LineColor(root.kRed+1),root.RooFit.Range("lowsighigh"))

#axisHist = root.TH2F("axisHist","",1,110,160,1,0,1200)
axisHist = root.TH2F("axisHist","",1,massLowRange[0],massHighRange[1],1,0,1200)
axisHist.Draw()
hstack.Draw('histsame')
plotMmumu.Draw('same')
canvas.SaveAs("TestMSSM.png")

###########################################################################

print("chi2/ndf: {}".format(plotMmumu.chiSquare(pdfMmumu.getParameters(mMuMuRooDataHist).getSize())))
for i in [bwmZ,bwSig,expLambda,mixParam]:
  print("{0}: {1:.3g} +/- {2:.3%}".format(i.GetName(),i.getVal(),abs(i.getError()/i.getVal())))
yay = raw_input("Press Enter to continue...")

