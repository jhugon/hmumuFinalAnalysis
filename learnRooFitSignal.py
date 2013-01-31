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

maxMass = 140.
minMass = 110.
mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)

"""
mean = root.RooRealVar("Mean","Mean",125.,100.,150.)
mean2 = root.RooRealVar("Mean2","Mean2",0.0)

width = root.RooRealVar("Width","Width",2.0,0.1,20.0)
width2 = root.RooRealVar("Width2","Width2",0.5,0.1,20.0)

Alpha = root.RooRealVar("Alpha","Alpha",1.0,0.1,20.0)
n = root.RooRealVar("n","n",2.0,0.2,10.0)
mix = root.RooRealVar("mix","mix",0.5,0.0,1.0)

#gaus1 = root.RooGaussian("gaus1","gaus1",mMuMu,mean,width)
gaus1 = root.RooCBShape("gaus1","gaus1",mMuMu,mean,width,Alpha,n)
gaus2 = root.RooGaussian("gaus2","gaus2",mMuMu,mean2,width2)
#mix = root.RooFormulaVar("mix","mMuMu > Alpha",root.RooArgList(mMuMu,Alpha))

pdfMmumu = root.RooFFTConvPdf("pdf","pdf",mMuMu,gaus1,gaus2)



"""
#Got 1.64 Chi2/ndf
mean = root.RooRealVar("Mean","Mean",125.,100.,150.)

width = root.RooRealVar("Width","Width",2.0,0.1,20.0)
width2 = root.RooRealVar("Width2","Width2",0.5,0.1,20.0)

Alpha = root.RooRealVar("Alpha","Alpha",1.0,0.1,20.0)
n = root.RooRealVar("n","n",2.0,0.2,10.0)
mix = root.RooRealVar("mix","mix",0.5,0.0,1.0)

#gaus1 = root.RooGaussian("gaus1","gaus1",mMuMu,mean,width)
gaus1 = root.RooCBShape("gaus1","gaus1",mMuMu,mean,width,Alpha,n)
gaus2 = root.RooGaussian("gaus2","gaus2",mMuMu,mean,width2)
#mix = root.RooFormulaVar("mix","mMuMu > Alpha",root.RooArgList(mMuMu,Alpha))

pdfMmumu = root.RooAddPdf("pdf","pdf",gaus1,gaus2,mix)

#####################################################################

f = root.TFile("input/freezeSample/ggHmumu125_7TeV.root")
#f = root.TFile("input/ttbar.root")

#mDiMu = f.Get("mDiMu")
mDiMu = f.Get("IncPresel/mDiMu")
#mDiMu = f.Get("VBFPresel/mDiMu")
#mDiMu.Rebin(2)

#####################################################################

mMuMuRooDataHist = root.RooDataHist("template","template",root.RooArgList(mMuMu),mDiMu)

pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL)

plotMmumu = mMuMu.frame(minMass,maxMass)

mMuMuRooDataHist.plotOn(plotMmumu)

pdfMmumu.plotOn(plotMmumu)

plotMmumu.Draw()

chi2ondf =  plotMmumu.chiSquare()
tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.05)
tlatex.SetTextAlign(22)
tlatex.DrawLatex(0.75,0.85,"#chi^{{2}}/NDF = {0:.2f}".format(chi2ondf))
canvas.SaveAs("testRes.png")

for i in [mean,width,width2,mix,Alpha,n]:
  print("{0}: {1:.3g} +/- {2:.3%}".format(i.GetName(),i.getVal(),abs(i.getError()/i.getVal())))
yay = raw_input("Press Enter to continue...")

