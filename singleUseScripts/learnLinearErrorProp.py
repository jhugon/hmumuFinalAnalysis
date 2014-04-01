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
b = root.RooRealVar("b","b",3.,-5,5)
pdf = root.RooGenericPdf("pdfPol","0.2+@0*@1+@0*@0*@2",root.RooArgList(x,a,b))

data = pdf.generateBinned(root.RooArgSet(x),10000)
nEvents = data.sumEntries()

fr = pdf.fitTo(data,root.RooFit.Save(),root.RooFit.Minos(),PRINTLEVEL)
print rooDebugFR(fr)

frame = x.frame()
data.plotOn(frame)
pdf.plotOn(frame,root.RooFit.VisualizeError(fr,1,True))
pdf.plotOn(frame)
data.plotOn(frame)

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
covab = fr.correlation(a,b)*a.getError()*b.getError()
covbb = fr.correlation(b,b)*b.getError()*b.getError()
aVal = a.getVal()
bVal = b.getVal()
  
errUpGraph = root.TGraph()
errDownGraph = root.TGraph()
errUpGraph.SetLineColor(root.kRed)
errDownGraph.SetLineColor(root.kRed)
errUpGraph.SetMarkerColor(root.kRed)
errDownGraph.SetMarkerColor(root.kRed)
for iBin in range(20):
  xNow = iBin*0.5
  x.setVal(xNow)
  a.setVal(aVal); b.setVal(bVal)
  aDeriv = pdf.derivative(a).getVal()
  a.setVal(aVal); b.setVal(bVal)
  bDeriv = pdf.derivative(b).getVal()
  a.setVal(aVal); b.setVal(bVal)
  variance = 0.
  variance += aDeriv*aDeriv*covaa
  variance += aDeriv*bDeriv*covab
  variance += bDeriv*bDeriv*covbb
  unc = sqrt(variance)
  myVariance =  0.
  myVariance += (xNow*a.getError())**2 ## a unc 
  myVariance += (xNow**2*b.getError())**2  ## b unc
  myVariance += (xNow**3*b.getError()*a.getError()*fr.correlation(a,b))  ## ab unc
  myUnc = sqrt(myVariance)
  
  rangeName = "binAgain{0}".format(iBin)
  x.setRange(rangeName,xNow-0.25,xNow+0.25)
  a.setVal(aVal); b.setVal(bVal)
  pdfVal = pdf.createIntegral(observables,observables,rangeName).getVal()
  a.setVal(aVal); b.setVal(bVal)
  funcVal = pdf.getVal()
  relUnc = unc/funcVal
  myFuncVal = 0.2+a.getVal()*xNow+b.getVal()*xNow**2

  myFuncVal = 0.2+a.getVal()*xNow+b.getVal()*xNow**2
  myRelUnc = myUnc/myFuncVal

  #print 
  #print xNow
  #print "unc: ",unc,myUnc
  #print "relUnc: ",relUnc,myRelUnc
  #print "funcVal: ",funcVal,myFuncVal

  errUpGraph.SetPoint(iBin,xNow,pdfVal*pdfSF*(1.+relUnc))
  errDownGraph.SetPoint(iBin,xNow,pdfVal*pdfSF*(1.-relUnc))
  #errUpGraph.SetPoint(iBin,xNow,pdfVal*pdfSF+unc)
  #errDownGraph.SetPoint(iBin,xNow,pdfVal*pdfSF-unc)
errUpGraph.Draw("P")
errDownGraph.Draw("P")

saveAs(canvas,"output/ErrBand")


### Need to scale uncertainty by the scalefactor from the function to the plotted fit
### Can just take the plotted value and divide by pdf.getVal() at the xpoint

### Important to reset all params to best fit values after derivative operations


