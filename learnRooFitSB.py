#! /usr/bin/env python

from ROOT import gSystem

import sys
import os
import re
import math
from ROOT import *
gSystem.Load('libRooFit')
import ROOT as root

from helpers import *

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

canvas = root.TCanvas()

maxMass = 160.
minMass = 100.
mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
mMuMu.setRange("all",minMass,maxMass)
mMuMu.setRange("sigFit",110,140)

#####################################################################

mean = root.RooRealVar("Mean","Mean",125.,100.,150.)
mean2 = root.RooRealVar("Mean2","Mean2",125.,100.,150.)

width = root.RooRealVar("Width","Width",2.0,0.1,20.0)
width2 = root.RooRealVar("Width2","Width2",0.5,0.1,20.0)

mix = root.RooRealVar("mix","mix",0.5,0.0,1.0)

sigParams = [mean,mean2,width,width2,mix]

gaus1 = root.RooGaussian("gaus1","gaus1",mMuMu,mean,width)
gaus2 = root.RooGaussian("gaus2","gaus2",mMuMu,mean2,width2)

pdfSig = root.RooAddPdf("pdfSig","pdfSig",gaus1,gaus2,mix)

#####################################################################

channelName = ''

InvPolMass = root.RooRealVar("InvPolMass","InvPolMass", 91., 70., 95.)
pdfBak = root.RooGenericPdf("bak","@0/(@0-@1)/(@0-@1)",root.RooArgList(mMuMu,InvPolMass))
bakParams = [InvPolMass]

#InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91., 70., 95.)
#ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -1., 1.)
#pdfBak = root.RooGenericPdf("bak","TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",root.RooArgList(mMuMu,InvPolMass,ExpMass))
#bakParams = [InvPolMass,ExpMass]

#p1 = root.RooRealVar(channelName+"_p1","p1", 0.0, -1., 1.)
#p2 = root.RooRealVar(channelName+"_p2","p2", 0.0, -1., 1.)
#p3 = root.RooRealVar(channelName+"_p3","p3", 0.0, -1., 1.)
#pdfBak = root.RooGenericPdf("bak","TMath::Exp(@0*@0*@1 + @0*@2 + @3*TMath::Log(@0) )",root.RooArgList(mMuMu,p1,p2,p3))
#bakParams = [p1,p2,p3]

#####################################################################

bakStrength = root.RooRealVar("bakStrength","bakStrength",1.0,1e6)
sigStrength = root.RooRealVar("sigStrength","sigStrength",0.0,1e6)
epdfSig = root.RooExtendPdf("epdfSig","epdfSig",pdfSig,sigStrength)
epdfBak = root.RooExtendPdf("epdfBak","epdfBak",pdfBak,bakStrength)
pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(epdfBak,epdfSig))

extParams = [sigStrength,bakStrength]

#####################################################################

fbak = root.TFile("input/preApproveSample/DYJetsToLL_8TeV.root")
fsig = root.TFile("input/preApproveSample/ggHmumu125_8TeV.root")
#f = root.TFile("input/ttbar.root")

mDiMuBak = fbak.Get("IncPreselPtG10BB/mDiMu")
mDiMuSig = fsig.Get("IncPreselPtG10BB/mDiMu")
#mDiMu = f.Get("VBFPresel/mDiMu")
#mDiMu.Rebin(2)

#####################################################################

mMuMuRooDataHistBak = root.RooDataHist("templatebak","template",root.RooArgList(mMuMu),mDiMuBak)
mMuMuRooDataHistSig = root.RooDataHist("templatesig","template",root.RooArgList(mMuMu),mDiMuSig)

epdfSig.fitTo(mMuMuRooDataHistSig,root.RooFit.Range("sigFit"),root.RooFit.SumW2Error(False),PRINTLEVEL)

for i in sigParams:
  i.setConstant(True)

pdfMmumu.fitTo(mMuMuRooDataHistBak,root.RooFit.Range("all"),root.RooFit.SumW2Error(False),PRINTLEVEL)

plotMmumu = mMuMu.frame(minMass,maxMass)

mMuMuRooDataHistBak.plotOn(plotMmumu)

pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("epdfSig"),root.RooFit.LineColor(root.kGreen+1))
pdfMmumu.plotOn(plotMmumu)
pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("epdfBak"),root.RooFit.LineColor(root.kRed+1))

plotMmumu.Draw()

chi2ondf =  plotMmumu.chiSquare()
tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.05)
tlatex.SetTextAlign(22)
tlatex.DrawLatex(0.75,0.85,"#chi^{{2}}/NDF = {0:.2f}".format(chi2ondf))
canvas.SaveAs("testRes.png")

for i in sigParams+bakParams+extParams:
  print("{0}: {1:.3g} +/- {2:.3%}".format(i.GetName(),i.getVal(),abs(i.getError()/i.getVal())))
yay = raw_input("Press Enter to continue...")

