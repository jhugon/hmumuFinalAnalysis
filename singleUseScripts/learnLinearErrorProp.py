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
binning = x.getBinning()
binWidth = binning.averageBinWidth()
observables = root.RooArgSet(x)

a = root.RooRealVar("a","a",0.01,-5,5)
b = root.RooRealVar("b","b",-0.0005,-5,5)
pdf = root.RooGenericPdf("pdfPol","0.2+@0*@1+@0*@0*@2",root.RooArgList(x,a,b))

data = pdf.generateBinned(root.RooArgSet(x),300)
nEvents = data.sumEntries()

fr = pdf.fitTo(data,root.RooFit.Save(),root.RooFit.Minos(),PRINTLEVEL)
aVal = a.getVal()
bVal = b.getVal()
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

b.setVal(bVal+b.getError())
pdf.plotOn(frame,root.RooFit.LineColor(root.kGreen+1))
b.setVal(bVal-b.getError())
pdf.plotOn(frame,root.RooFit.LineColor(root.kGreen+1))
b.setVal(bVal)

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

errGraph = root.TGraphAsymmErrors()
errGraph.SetLineColor(root.kRed)
errGraph.SetMarkerColor(root.kRed)
for iBin in range(20):
  xNow = binWidth/2.+iBin*binWidth

  rangeName = "binAgain_{0}".format(iBin)
  x.setRange(rangeName,xNow-binWidth/2.,xNow+binWidth/2.)
  pdfInt =  pdf.createIntegral(observables,observables,rangeName)
  pdfVal = pdfInt.getVal()
  errGraph.SetPoint(iBin,xNow,nEvents * pdfVal)

  fpf = rooArgSet2List(fr.floatParsFinal())
  parList = rooArgSet2List(pdf.getParameters(data))
  pdfDerivParList = []
  for par in parList:
    originalVal = par.getVal()
    paramUnc = par.getError()

    par.setVal(originalVal+paramUnc)
    pdfErrVal = pdfInt.getVal()
    pointErr = abs(pdfErrVal-pdfVal)*nEvents
    par.setVal(originalVal-paramUnc)
    pointErr = abs(pdfErrVal-pdfVal)*nEvents
    pointErr = max(abs(pdfErrVal-pdfVal)*nEvents,pointErr)

    pdfDerivParList.append(pointErr/paramUnc)
    par.setVal(originalVal)

  totalVariance = 0.
  for i in range(len(parList)):
    iParam = parList[i]
    iUnc = iParam.getError()
    iDeriv = pdfDerivParList[i]
    for j in range(i,len(parList)):
      jParam = parList[j]
      jUnc = jParam.getError()
      jDeriv = pdfDerivParList[j]
      totalVariance += iDeriv*jDeriv*iUnc*jUnc*fr.correlation(iParam,jParam)
    
  totalUncertainty = sqrt(totalVariance)
  errGraph.SetPointError(iBin,0,0,totalUncertainty,totalUncertainty)
      
      
      

errGraph.Draw("P")

saveAs(canvas,"output/ErrBand")

### Need to scale uncertainty by the scalefactor from the function to the plotted fit
### Can just take the plotted value and divide by pdf.getVal() at the xpoint

### Important to reset all params to best fit values after derivative operations


