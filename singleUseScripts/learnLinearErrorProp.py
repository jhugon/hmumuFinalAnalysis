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
x.setBins(20)
observables = root.RooArgSet(x)

a = root.RooRealVar("a","a",0.1,-5,5)
pdf = root.RooGenericPdf("pdfPol","0.2+@0*@1",root.RooArgList(x,a))

data = pdf.generateBinned(root.RooArgSet(x),100)
nEvents = data.sumEntries()

fr = pdf.fitTo(data,root.RooFit.Save(),root.RooFit.Minos(),PRINTLEVEL)
aVal = a.getVal()
print rooDebugFR(fr)

frame = x.frame()
data.plotOn(frame)
pdf.plotOn(frame,root.RooFit.VisualizeError(fr,1,True))
pdf.plotOn(frame)
data.plotOn(frame)
a.setVal(aVal+a.getError())
pdf.plotOn(frame)
a.setVal(aVal-a.getError())
pdf.plotOn(frame)
a.setVal(aVal)

canvas = root.TCanvas("canvas")
frame.Draw()

x.setRange("all",0,10)
pdfSF = nEvents
# I can draw my own PDF line!!!
pdfGraph2 = root.TGraph()
pdfGraph2.SetMarkerColor(root.kRed+1)
pdfGraph2.SetLineColor(root.kRed+1)
for iBin in range(100):
  rangeName = "bin{0}".format(iBin)
  x.setRange(rangeName,iBin*10./100-0.25,iBin*10./100+0.25)
  pdfVal = pdf.createIntegral(observables,observables,rangeName).getVal()
  pdfGraph2.SetPoint(iBin,iBin*10./100,pdfVal*pdfSF)
pdfGraph2.Draw("L")

#######################################################333
#######################################################333
#######################################################333
#######################################################333
#######################################################333

covaa = fr.correlation(a,a)*a.getError()*a.getError()
  
errUpGraph = root.TGraph()
errDownGraph = root.TGraph()
errUpGraph.SetLineColor(root.kRed)
errDownGraph.SetLineColor(root.kRed)
errUpGraph.SetMarkerColor(root.kRed)
errDownGraph.SetMarkerColor(root.kRed)
for iBin in range(20):
  xNow = iBin*0.5
  x.setVal(xNow)
  a.setVal(aVal)
  aDeriv = pdf.derivative(a).getVal()
  a.setVal(aVal)
  variance = 0.
  variance += aDeriv*aDeriv*covaa
  unc = sqrt(variance)
  myVariance =  0.
  myVariance += (xNow*a.getError())**2 ## a unc 
  myUnc = sqrt(myVariance)
  
  rangeName = "binAgain{0}".format(iBin)
  x.setRange(rangeName,xNow-0.25,xNow+0.25)
  a.setVal(aVal)
  pdfVal = pdf.createIntegral(observables,observables,rangeName).getVal()
  a.setVal(aVal)
  funcVal = pdf.getVal()
  relUnc = unc/funcVal

  myFuncVal = 0.2+a.getVal()*xNow
  myRelUnc = myUnc/myFuncVal

  print 
  print xNow
  print "aDeriv: ",aDeriv,"  covaa: ",covaa, "a unc: ",sqrt(covaa)
  print "unc: ",unc,myUnc
  print "relUnc: ",relUnc,myRelUnc
  print "funcVal: ",funcVal,myFuncVal
  print "pdfVal*pdfSF: ",pdfVal*pdfSF

  errUpGraph.SetPoint(iBin,xNow,pdfVal*pdfSF*(1.+relUnc))
  errDownGraph.SetPoint(iBin,xNow,pdfVal*pdfSF*(1.-relUnc))
  #errUpGraph.SetPoint(iBin,xNow,pdfVal*pdfSF+unc*pdfVal/funcVal*pdfSF)
  #errDownGraph.SetPoint(iBin,xNow,pdfVal*pdfSF-unc*pdfVal/funcVal*pdfSF)
errUpGraph.Draw("P")
errDownGraph.Draw("P")

saveAs(canvas,"output/ErrBand")


### Need to scale uncertainty by the scalefactor from the function to the plotted fit
### Can just take the plotted value and divide by pdf.getVal() at the xpoint

### Important to reset all params to best fit values after derivative operations


