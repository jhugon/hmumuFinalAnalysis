#!/usr/bin/env python

import singleHelpers
import math
import ROOT as root
from helpers import *
import datetime
import sys
import os.path
import copy
import time

from ROOT import gSystem
gSystem.Load('libRooFit')
root.gROOT.SetBatch(True)

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT


x = root.RooRealVar("x","x",0,10)
x.setRange("all",0,10)
x.setRange("fitRange",0,2)
x.setBins(20)
observables = root.RooArgSet(x)

a = root.RooRealVar("a","a",0.1,-5,5)
pdf = root.RooGenericPdf("pdfPol","0.2+@0*@1",root.RooArgList(x,a))

data = pdf.generateBinned(root.RooArgSet(x),1000)
nEvents = data.sumEntries()

fr = pdf.fitTo(data,root.RooFit.Save(),root.RooFit.Minos(),PRINTLEVEL,root.RooFit.Range("fitRange"))
print rooDebugFR(fr)

frame = x.frame()
curveData = data.plotOn(frame)
curveBand = pdf.plotOn(frame,root.RooFit.VisualizeError(fr,1,True),root.RooFit.Range("all"))
curveLine = pdf.plotOn(frame,root.RooFit.Range("all"))
curveData2 = data.plotOn(frame)

canvas = root.TCanvas("canvas")
frame.Draw()
rptip = RooPredictionTolerenceIntervalPlotter(x,pdf,data,fr)
rptip.drawTolerence()
rptip.drawPrediction()
rptip.drawStatOnly()

frame.Draw("same")

saveAs(canvas,"output/PredTolIntervals")
