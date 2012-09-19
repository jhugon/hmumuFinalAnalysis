#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

canvas = root.TCanvas()

mMuMu = root.RooRealVar("mMuMu","mMuMu",110.0,180.0)
mMuMu.setRange("low",110,120)
mMuMu.setRange("high",130,180)
#mva = root.RooRealVar("mva","mva",-1,1)
mva = root.RooRealVar("mva","mva",-0.8,-0.1)

lmZ = root.RooRealVar("lmZ","lmZ",85,95)
lWidth = root.RooRealVar("lWidth","lWidth",0.0,30.0)

landauMmumu = root.RooLandau("landauMmumu","landauMmumu",mMuMu,lmZ,lWidth)

gWidth = root.RooRealVar("gWidth","gWidth",0.0,10.0)
gmZ = root.RooRealVar("gmZ","gmZ",0.0)
gausMmumu = root.RooGaussian("gausMmumu","gausMmumu",mMuMu,gmZ,gWidth)

bwWidth = root.RooRealVar("bwWidth","bwWidth",0.0,30.0)
bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
bwMmumu = root.RooBreitWigner("bwMmumu","bwMmumu",mMuMu,bwmZ,bwWidth)

voitWidth = root.RooRealVar("voitWidth","voitWidth",0.0,30.0)
voitmZ = root.RooRealVar("voitmZ","voitmZ",85,95)
voitSig = root.RooRealVar("voitSig","voitSig",0.0,30.0)
voitMmumu = root.RooVoigtian("voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

lnA = root.RooRealVar("lnA","lnA",0,100)
lnB = root.RooRealVar("lnB","lnB",0.0,30.0)
lognormalMmumu = root.RooLognormal("lognormalMmumu","lognormalMmumu",mMuMu,lnA,lnB)

#pdfMmumu = root.RooFFTConvPdf("pdfMmumu","pdfMmumu",mMuMu,landauMmumu,gausMmumu)
#pdfMmumu = landauMmumu
#pdfMmumu = lognormalMmumu
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(landauMmumu,bwMmumu))
pdfMmumu = bwMmumu
#pdfMmumu = voitMmumu

f = root.TFile("input/DYJetsToLL.root")

mDiMu = f.Get("mDiMu")
mDiMu.Rebin(2)

templateMmumu = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

pdfMmumu.fitTo(templateMmumu,root.RooFit.Range("low,high"))
#pdfMmumu.fitTo(templateMmumu)

plotMmumu = mMuMu.frame()

templateMmumu.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu)
plotMmumu.Draw()
canvas.SaveAs("learnMmumu.png")

###################

BDTHistMuonOnlyVMass = f.Get("BDTHistMuonOnlyVMass")
tmpHist = BDTHistMuonOnlyVMass.Clone("tmpHist")
tmpAxis = tmpHist.GetXaxis()
BDTHistLowControl = tmpHist.ProjectionY("_LowControl",tmpAxis.FindBin(110.0),tmpAxis.FindBin(115.0))
BDTHistHighControl = tmpHist.ProjectionY("_HighControl",tmpAxis.FindBin(125.0),tmpAxis.FindBin(150.0))

print("N Events in Low Control: {0:.1f}, in High Control: {1:.1f}".format(BDTHistLowControl.Integral(),BDTHistHighControl.Integral()))

BDTHist = BDTHistLowControl.Clone()
BDTHist.Add(BDTHistHighControl)

BDTHist.Rebin(BDTHist.GetNbinsX()/100)

templateMva = root.RooDataHist("templateMva","templateMva",root.RooArgList(mva),BDTHist)

pdfMva = root.RooHistPdf("pdfMva","pdfMva",root.RooArgSet(mva),templateMva)

pdfMva.fitTo(templateMva)

###################

plotMva = mva.frame()

templateMva.plotOn(plotMva)
pdfMva.plotOn(plotMva)
plotMva.Draw()
canvas.SaveAs("learnMva.png")

###################

canvas.Clear()
canvas.Divide(2,2)

multipdf = root.RooProdPdf("multipdf","multipdf",RooArgList(landauMmumu,pdfMva))

mDiMuBinning = root.RooFit.Binning(50,110,150)
yvar = root.RooFit.YVar(mva,root.RooFit.Binning(50,-1,1))
print mDiMu
print mDiMuBinning
print yvar
multipdfHist = multipdf.createHistogram("hh_model_1",mMuMu,mDiMuBinning,yvar)

canvas.cd(3)

multipdfHist.GetXaxis().SetRangeUser(110,150)
multipdfHist.Draw("surf")
BDTHistMuonOnlyVMass.SetLineColor(kRed)
BDTHistMuonOnlyVMass.Rebin2D(BDTHistMuonOnlyVMass.GetNbinsX()/400,BDTHistMuonOnlyVMass.GetNbinsY()/50)
BDTHistMuonOnlyVMass.GetXaxis().SetRangeUser(110,150)
BDTHistMuonOnlyVMass.Draw("surf same")

#canvas.SaveAs("learn2D.png")

BDTHistMuonOnlyVMass.Draw("surf")
multipdfHist.GetXaxis().SetRangeUser(110,150)
multipdfHist.Draw("surf same")

canvas.cd(1)
BDTHistMuonOnlyVMass.Draw("surf")
multipdfHist.Draw("surf same")

canvas.cd(2)
BDTHistMuonOnlyVMass.Draw("colz")

canvas.cd(4)
multipdfHist.Draw("colz")

canvas.SaveAs("learn2D.png")

