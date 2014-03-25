#!/usr/bin/env python

import singleHelpers
from helpers import *
from xsec import *
import cPickle
import sys
import os
import os.path
import glob
import re

root.gROOT.SetBatch(True)

canvas = root.TCanvas()

errorsDictFile = open("pklfiles/oneSig.pkl")
ERRORSDICT = cPickle.load(errorsDictFile)
errorsDictFile.close()

fileNotFix = root.TFile("order_Shape_Jets01PassPtG10BB_Bernstein4_params_pullHist_master.root")
fileFix = root.TFile("order_Shape_Jets01PassPtG10BB_Bernstein4_params_pullHist_bernFixAParam.root")

histNotFix = fileNotFix.Get("pullHist")
histFix = fileFix.Get("pullHist")

histFix.SetLineColor(root.kRed+1)
histNotFix.SetFillStyle(0)
histFix.SetFillStyle(0)

axisHist = root.TH2F("axisHist","",1,110,160,1,-5,5)
setHistTitles(axisHist,"M(#mu#mu) [GeV/c^{2}]","#frac{Data-Fit}{#sqrt{Fit}}")

zeroLine = root.TGraph()
zeroLine.SetLineStyle(2)
zeroLine.SetPoint(0,110,0.)
zeroLine.SetPoint(1,160,0.)

axisHist.Draw()
zeroLine.Draw("L")
histNotFix.Draw("same")
histFix.Draw("same")

drawStandardCaptions(canvas,"0,1-Jet Tight BB","#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}",preliminaryString="CMS Internal")

canvas.SaveAs("PullsCompare.png")
canvas.SaveAs("PullsCompare.pdf")
canvas.SaveAs("PullsCompare.eps")

### Difference

histNotFix.Add(histFix,-1)
histNotFix.SetLineColor(1)
histNotFix.UseCurrentStyle()
setHistTitles(histNotFix,"M(#mu#mu) [GeV/c^{2}]","Free Pull - Fixed Pull")
histNotFix.Draw()

drawStandardCaptions(canvas,"0,1-Jet Tight BB","#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}",preliminaryString="CMS Internal")

canvas.SaveAs("PullsDiff.png")
canvas.SaveAs("PullsDiff.pdf")
canvas.SaveAs("PullsDiff.eps")

###################
# Fits

histNotFix = fileNotFix.Get("fitHist")
histFix = fileFix.Get("fitHist")

histFix.SetLineColor(root.kRed+1)
histNotFix.SetFillStyle(0)
histFix.SetFillStyle(0)

axisHist = root.TH2F("axisHistFit","",1,110,160,1,0,900)
setHistTitles(axisHist,"M(#mu#mu) [GeV/c^{2}]","Fitted Events/Bin")

axisHist.Draw()
#zeroLine.Draw("L")
histNotFix.Draw("same")
histFix.Draw("same")

drawStandardCaptions(canvas,"0,1-Jet Tight BB","#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}",preliminaryString="CMS Internal")

canvas.SaveAs("Pulls_FitsCompare.png")
canvas.SaveAs("Pulls_FitsCompare.pdf")
canvas.SaveAs("Pulls_FitsCompare.eps")

  
