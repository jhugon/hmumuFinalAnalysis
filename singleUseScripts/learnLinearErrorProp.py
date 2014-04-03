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
  print iBin*10/100.,pdfVal*pdfSF
pdfGraph2.Draw("L")

#######################################################333
#######################################################333
#######################################################333
#######################################################333
#######################################################333

def myFunc(x,a):
  return 0.2+a*x
def myFuncIntX(x,a,xmin=0.,xmax=10.):
  return 0.2*xmax - 0.2*xmin + a*xmax**2/2. - a*xmin**2/2.
def myFuncDerivX(x,a):
  return a
def myNormFactorDerivA(x,a,xmin=0.,xmax=10.):
  return -0.5*(xmax**2-xmin**2)*myFuncIntX(x,a,xmin,xmax)**2

covaa = fr.correlation(a,a)*a.getError()*a.getError()
  
errGraph = root.TGraphAsymmErrors()
errGraph.SetLineColor(root.kRed)
errGraph.SetMarkerColor(root.kRed)
for iBin in range(20):
  xNow = binWidth/2.+iBin*binWidth
  intAll = 10**2/2.*a.getVal()
  C = 1./myFuncIntX(xNow,a.getVal())
  myPdfVal = nEvents * C * myFuncIntX(xNow,a.getVal(),xNow-binWidth/2.,xNow+binWidth/2.)
  errGraph.SetPoint(iBin,xNow,myPdfVal)


  deriv = nEvents * C * myFuncDerivX(xNow,a.getVal())
  deriv += nEvents * myFunc(xNow,a.getVal()) * myNormFactorDerivA(xNow,a.getVal())

  error = deriv*sqrt(covaa)

  errGraph.SetPointError(iBin,0,0,error,error)
  print
  print "x: ",xNow
  print "pdf: ",myPdfVal
  print "func deriv: ",myFuncDerivX(xNow,a.getVal())
  print "C: ",C
  print "dC/da: ", myNormFactorDerivA(xNow,a.getVal())
  print "pdf deriv: ",deriv
  print "unc(a) :", sqrt(covaa)
  print "error: ",error
errGraph.Draw("P")

err2Graph = root.TGraphAsymmErrors()
err2Graph.SetLineColor(root.kGreen+1)
err2Graph.SetMarkerColor(root.kGreen+1)
for iBin in range(20):
  xNow = binWidth/2.+iBin*binWidth
  intAll = 10**2/2.*a.getVal()
  C = 1./myFuncIntX(xNow,a.getVal())
  myPdfVal = nEvents * C * myFunc(xNow,a.getVal())*5
  err2Graph.SetPoint(iBin,xNow,myPdfVal)
err2Graph.Draw("P")

saveAs(canvas,"output/ErrBand")

### Need to scale uncertainty by the scalefactor from the function to the plotted fit
### Can just take the plotted value and divide by pdf.getVal() at the xpoint

### Important to reset all params to best fit values after derivative operations


