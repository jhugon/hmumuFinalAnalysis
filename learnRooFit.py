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

mMuMu = root.RooRealVar("mMuMu","mMuMu",110,140.0)
#mMuMu.setRange("low",110,120)
#mMuMu.setRange("high",130,180)
mMuMu.setRange("low",110,124)
mMuMu.setRange("high",126,140)
#mva = root.RooRealVar("mva","mva",-1,1)
mva = root.RooRealVar("mva","mva",-0.8,-0.1)

lmZ = root.RooRealVar("lmZ","lmZ",85,95)
lWidth = root.RooRealVar("lWidth","lWidth",0.0,30.0)

landauMmumu = root.RooLandau("landauMmumu","landauMmumu",mMuMu,lmZ,lWidth)

gWidth = root.RooRealVar("gWidth","gWidth",0.0,50.0)
gmZ = root.RooRealVar("gmZ","gmZ",91.19)
gausMmumu = root.RooGaussian("gausMmumu","gausMmumu",mMuMu,gmZ,gWidth)

#bwWidth = root.RooRealVar("bwWidth","bwWidth",0.0,30.0)
bwWidth = root.RooRealVar("bwWidth","bwWidth",2.4952)
bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
bwMmumu = root.RooBreitWigner("bwMmumu","bwMmumu",mMuMu,bwmZ,bwWidth)

voitWidth = root.RooRealVar("voitWidth","voitWidth",2.4952)
voitmZ = root.RooRealVar("voitmZ","voitmZ",85,95)
voitSig = root.RooRealVar("voitSig","voitSig",0.0,30.0)
voitMmumu = root.RooVoigtian("voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

"""
lnA = root.RooRealVar("lnA","lnA",0,100)
lnB = root.RooRealVar("lnB","lnB",0.0,30.0)
lognormalMmumu = root.RooLognormal("lognormalMmumu","lognormalMmumu",mMuMu,lnA,lnB)
"""

expParam = root.RooRealVar("expParam","expParam",-1,0)
expMmumu = root.RooExponential("expMmumu","expMmumu",mMuMu,expParam)
"""
expTimesVoig = root.RooProdPdf("expTimesVoig","expTimesVoig",root.RooArgList(expMmumu,voitMmumu))

dyPdf1 = root.RooGenericPdf("dyPdf1","dyPdf1","gmZ*TMath::Power(mMuMu*mMuMu-gmZ*gmZ,2)/(TMath::Power(mMuMu*mMuMu-gmZ*gmZ,2)+mMuMu*mMuMu*gWidth*gWidth/gmZ/gmZ)",RooArgList(mMuMu,gmZ,gWidth))

dyPdf2 = root.RooProdPdf("dyPdf2","dyPdf2",root.RooArgList(expMmumu,dyPdf1))
dyPdf = root.RooAddPdf("dyPdf","dyPdf",root.RooArgList(dyPdf2,expTimesVoig))

alpha = root.RooRealVar("alpha","alpha",0,50)
kappa = root.RooRealVar("kappa","kappa",0,5)
dyHighMassPdf = root.RooGenericPdf("dyHighMassPdf","dyHighMassPdf","TMath::Exp(-alpha*TMath::Power(mMuMu,kappa))",root.RooArgList(mMuMu,alpha,kappa))
"""

a1 = root.RooRealVar("a1","a1",-5,5)
a2 = root.RooRealVar("a2","a2",-5,5)
a3 = root.RooRealVar("a3","a3",-5,5)
a4 = root.RooRealVar("a4","a4",-5,5)
a5 = root.RooRealVar("a5","a5",-5,5)
a6 = root.RooRealVar("a6","a6",-5,5)
a7 = root.RooRealVar("a7","a7",-5,5)
polyMmumu = root.RooChebychev("polyMmumu","polyMmumu",mMuMu,root.RooArgList(a1,a2))
#polyMmumu = root.RooPolynomial("polyMmumu","polyMmumu",mMuMu,root.RooArgList(a1,a2,a3,a4,a5,a6,a7))

bwCoef = root.RooRealVar("bwCoef","bwCoef",0,1)
expCoef = root.RooRealVar("expCoef","expCoef",0,1)
polyPlusExp = root.RooAddPdf("polyPlusExp","bwPlusExp",root.RooArgList(expMmumu,polyMmumu),root.RooArgList(expCoef))
gausPlusExp = root.RooAddPdf("gausPlusExp","bwPlusExp",root.RooArgList(expMmumu,gausMmumu),root.RooArgList(expCoef))

#pdfMmumu = root.RooFFTConvPdf("pdfMmumu","pdfMmumu",mMuMu,landauMmumu,gausMmumu)
#pdfMmumu = landauMmumu
#pdfMmumu = expTimesVoig
#pdfMmumu = dyPdf
#pdfMmumu = dyHighMassPdf
#pdfMmumu = lognormalMmumu
#pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(landauMmumu,bwMmumu))
#pdfMmumu = bwMmumu
#pdfMmumu = voitMmumu
#pdfMmumu = polyMmumu
#pdfMmumu = expMmumu
pdfMmumu = gausPlusExp
#pdfMmumu = polyPlusExp


f = root.TFile("input/DYJetsToLL.root")
#f = root.TFile("input/ttbar.root")

mDiMu = f.Get("mDiMu")
mDiMu.Rebin(2)

templateMmumu = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

fitResult = pdfMmumu.fitTo(templateMmumu,root.RooFit.Range("low,high"),root.RooFit.Save(True))
#pdfMmumu.fitTo(templateMmumu,root.RooFit.Range(130,160))
chi2 = pdfMmumu.createChi2(templateMmumu)

plotMmumu = mMuMu.frame()

templateMmumu.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu)
plotMmumu.Draw()

print fitResult
fitResult.Print()
print("chi2/ndf: {}".format(chi2.getVal()/(mDiMu.GetNbinsX()-1)))
yay = raw_input("Press Enter to continue...")
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

templatePdfMva = root.RooHistPdf("templatePdfMva","templatePdfMva",root.RooArgSet(mva),templateMva)

gMvaWidth = root.RooRealVar("gMvaWidth","gMvaWidth",0.0,10)
gMvaMean = root.RooRealVar("gMvaMean","gMvaMean",-1,1)
gausPdfMva = root.RooGaussian("gausPdfMva","gausPdfMva",mva,gMvaMean,gMvaWidth)

pdfMva = root.RooFFTConvPdf("pdfMva","pdfMva",mva,templatePdfMva,gausPdfMva)
#pdfMva = templatePdfMva

pdfMva.fitTo(templateMva)

###################

plotMva = mva.frame()

templateMva.plotOn(plotMva)
pdfMva.plotOn(plotMva)
templatePdfMva.plotOn(plotMva,root.RooFit.LineStyle(root.kDashed))
plotMva.Draw()
canvas.SaveAs("learnMva.png")

###################

canvas.Clear()
canvas.Divide(2,2)

multipdf = root.RooProdPdf("multipdf","multipdf",RooArgList(pdfMmumu,pdfMva))

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


############################################

canvas.Clear()
plotMmumu = mMuMu.frame()

pdfMmumu.paramOn(plotMmumu)
pdfMva.paramOn(plotMmumu)
multipdf.paramOn(plotMmumu)
plotMmumu.Draw()
canvas.SaveAs("learnParams.png")

