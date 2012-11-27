#!/usr/bin/env python

from helpers import *
import ROOT as root
import os
import sys

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gROOT.SetBatch(True)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

dataDir = "input/"

root.gErrorIgnoreLevel = root.kWarning

def makePlusMinusPlots(plus,minus):
  #plus.Reset()
  #minus.Reset()
  minus.SetMarkerColor(root.kRed+1)
  minus.SetLineColor(root.kRed+1)
  plus.SetMarkerColor(root.kBlack+1)
  plus.SetLineColor(root.kBlack+1)

  if plus.InheritsFrom("TH2"):
    pass
  else:
    setHistTitles(qOverPt,"|q|/p_{T}","Events/Bin")
  plus.Draw()
  minus.Draw("same")

inFile = root.TFile(dataDir+"EndpointSingleMuRun2012Cv2.root")

canvas = root.TCanvas()

qOverPt = inFile.Get("qOverPt")
qOverPt.UseCurrentStyle()
setHistTitles(qOverPt,"q/p_{T}","Events/Bin")
qOverPt.Draw()
saveAs(canvas,"qOverPtBoth")

histList = ["qOverPt","qOverPtVEta","qOverPtVPhi","qOverPtVPt"]

for i in histList:
  plus = inFile.Get(i+"Plus")
  minus = inFile.Get(i+"Minus")
  makePlusMinusPlots(plus,minus)
  if i == "qOverPt":
    saveAs(canvas,i)

def findEndpoint(tree,):
  ptInvRange = [0.0,0.04]
  ptInv = root.RooRealVar("ptInv","1/p_{T}",ptInvRange[0],ptInvRange[1])
  ptInv.setRange("peak",0.02,0.03)
  ptInv.setRange("peakMore",0.015,0.035)
  pt = root.RooRealVar("pt","p_{T}",20.0,500.0)
  eta = root.RooRealVar("eta","#eta",-2.4,2.4)
  phi = root.RooRealVar("phi","#phi",-3.14,3.14)
  q = root.RooRealVar("q","q",-1.0,1.0)
  k = root.RooRealVar("k","k",0.0,-10000.,10000.)

  norm = root.RooRealVar("norm","norm",0.0,10000000)
  l = root.RooRealVar("l","l",(ptInvRange[0]+ptInvRange[1])/2.)
  h = root.RooRealVar("h","h",(ptInvRange[1]-ptInvRange[0])/2.)
  l.Print()
  h.Print()

  mi = -100.
  ma = 100.
  polyVarList = []
  order = 5
  for polyVarI in range(order):
    polyVarName = "a{0}".format(polyVarI)
    polyVarList.append(root.RooRealVar(polyVarName,polyVarName,0.0))
  for polyVar in polyVarList:
    polyVar.setConstant(True)

  cuts = "abs(eta)>2.098"
  dataVars = root.RooArgSet(pt,eta,phi,ptInv,q)
  data = root.RooDataSet("data","data",tree,dataVars,cuts)
  data.Print()

  corrPtInvVars = root.RooArgList(ptInv,q,k)
  corrPtInv = root.RooFormulaVar("corrPtInv","abs(q*ptInv+k)",corrPtInvVars)

  x = root.RooFormulaVar("x","(ptInv-l)/h",root.RooArgList(ptInv,l,h))

  polyVars = root.RooArgList(*polyVarList)
  #polyPDF = root.RooChebychev("polyPDF","polyPDF",corrPtInv,polyVars)
  #polyPDF = root.RooPolynomial("polyPDF","polyPDF",corrPtInv,polyVars)
  polyPDF = root.RooChebychev("polyPDF","polyPDF",ptInv,polyVars)
  #polyPDF = root.RooPolynomial("polyPDF","polyPDF",x,polyVars)

  #polyExtPDF = root.RooExtendPdf("polyExtPDF","polyExtPDF",polyPDF,norm)
  #polyExtPDF.fitTo(data,PRINTLEVEL,root.RooFit.Range("peak"),root.RooFit.Extended(True))
  #polyExtPDF.fitTo(data,PRINTLEVEL,root.RooFit.Range("peakMore"))

  polyVarValues = []
  for polyVar in polyVarList:
    polyVar.setRange(mi,ma)
    polyVar.setConstant(False)
    polyPDF.fitTo(data)
    polyVarValues.append(polyVar.getVal())
    polyVar.setVal(0.0)
    polyVar.setConstant(True)

  print polyVarValues
  for polyVar,val in zip(polyVarList,polyVarValues):
    polyVar.setRange(mi,ma)
    polyVar.setConstant(False)
    polyVar.setVal(val)


  plotPtInv = ptInv.frame()
  data.plotOn(plotPtInv)
  polyPDF.plotOn(plotPtInv)
  plotPtInv.Draw()

  """
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(norm.GetName(),norm.getVal(),norm.getError()/norm.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a1.GetName(),a1.getVal(),a1.getError()/a1.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a2.GetName(),a2.getVal(),a2.getError()/a2.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a3.GetName(),a3.getVal(),a3.getError()/a3.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a4.GetName(),a4.getVal(),a4.getError()/a4.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a5.GetName(),a5.getVal(),a5.getError()/a5.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a6.GetName(),a6.getVal(),a6.getError()/a6.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a7.GetName(),a7.getVal(),a7.getError()/a7.getVal()))
  print("{0:<5} {1:<9.3g} {1:<8.3%}".format(a8.GetName(),a8.getVal(),a8.getError()/a8.getVal()))
  #print("{0:<5} {1:<9.3g} {1:<8.3%}".format(k.GetName(),k.getVal(),k.getError()/k.getVal()))
  """

  for a1 in polyVarList:
    print("{0:<5} {1:<9.3g}".format(a1.GetName(),a1.getVal()))

def silly():
  #q = root.RooRealVar("q","q",-0.8,0.8)
  q = root.RooRealVar("q","q",-1,1)
  #mi = -1e10
  #ma = 1e10
  mi = 0.001
  ma = 0.1

  mean = root.RooRealVar("mean","mean",0.1)
  sig = root.RooRealVar("sig","sig",0.4)

  a1 = root.RooRealVar("a1","a1",2.0,-10,10)
  a2 = root.RooRealVar("a2","a2",-0.5,-10,10.)
  a3 = root.RooRealVar("a3","a3",1.,-10,10)
  a4 = root.RooRealVar("a4","a4",1.,-10,10)


  x = q

  #polyVars = root.RooArgList(a1,a2,a3,a4,a5,a6,a7,a8)
  #polyVars = root.RooArgList(a1,a2,a3,a4,a5,a6,a7)
  #polyVars = root.RooArgList(a1,a2,a3,a4,a5,a6)
  #polyVars = root.RooArgList(a1,a2,a3,a4,a5)
  #polyVars = root.RooArgList(a1,a2,a3,a4)
  #polyVars = root.RooArgList(a1,a2,a3)
  polyVars = root.RooArgList(a1,a2,a3,a4)
  #polyPDF = root.RooChebychev("polyPDF","polyPDF",corrPtInv,polyVars)
  #polyPDF = root.RooPolynomial("polyPDF","polyPDF",corrPtInv,polyVars)
  polyPDF = root.RooChebychev("polyPDF","polyPDF",x,polyVars)
  #polyPDF = root.RooPolynomial("polyPDF","polyPDF",x,polyVars)

  gausPDF = root.RooGaussian("gausPDF","gausPDF",x,mean,sig)

  #polyExtPDF = root.RooExtendPdf("polyExtPDF","polyExtPDF",polyPDF,norm)
  #polyExtPDF.fitTo(data,PRINTLEVEL,root.RooFit.Range("peak"),root.RooFit.Extended(True))
  #polyExtPDF.fitTo(data,PRINTLEVEL,root.RooFit.Range("peakMore"))

  #data = gausPDF.generate(root.RooArgSet(q),1000)
  data = polyPDF.generate(root.RooArgSet(q),10000)

  polyPDF.fitTo(data)

  plotPtInv = x.frame()
  data.plotOn(plotPtInv)
  polyPDF.plotOn(plotPtInv)
  plotPtInv.Draw()

def silly2(tree):
  ptInvRange = [0.005,0.04]
  cuts = "abs(eta)>2.098"

  hist = root.TH1F("ptInv","ptInv",200,ptInvRange[0],ptInvRange[1])
  tree.Draw("ptInv >> ptInv")

  order = 18
  #params = linearChi2TH1(hist,order,funcName="cheb")
  params = linearChi2TH1(hist,order+1,funcName="poly")

  ptInv = root.RooRealVar("ptInv","1/p_{T}",ptInvRange[0],ptInvRange[1])
  ptInv.setRange("peak",0.02,0.03)
  ptInv.setRange("peakMore",0.015,0.035)
  pt = root.RooRealVar("pt","p_{T}",20.0,500.0)
  eta = root.RooRealVar("eta","#eta",-2.4,2.4)
  phi = root.RooRealVar("phi","#phi",-3.14,3.14)
  q = root.RooRealVar("q","q",-1.0,1.0)
  k = root.RooRealVar("k","k",0.0,-10000.,10000.)

  polyVarList = []
  for polyVarI in range(order):
    polyVarName = "a{0}".format(polyVarI)
    initVal = params[polyVarI+1]
    minVal = initVal/5.
    maxVal = initVal*5.
    if initVal<0.0:
      tmp = maxVal
      maxVal = minVal
      minVal = tmp
    tmpVar = root.RooRealVar(polyVarName,polyVarName,initVal,minVal,maxVal)
    polyVarList.append(tmpVar)

  dataVars = root.RooArgSet(pt,eta,phi,ptInv,q)
  data = root.RooDataSet("data","data",tree,dataVars,cuts)
  data.Print()

  corrPtInvVars = root.RooArgList(ptInv,q,k)
  corrPtInv = root.RooFormulaVar("corrPtInv","abs(q*ptInv+k)",corrPtInvVars)

  polyVars = root.RooArgList(*polyVarList)
  #polyPDF = root.RooChebychev("polyPDF","polyPDF",corrPtInv,polyVars)
  #polyPDF = root.RooPolynomial("polyPDF","polyPDF",corrPtInv,polyVars)
  #polyPDF = root.RooChebychev("polyPDF","polyPDF",ptInv,polyVars)
  polyPDF = root.RooPolynomial("polyPDF","polyPDF",ptInv,polyVars)

  polyPDF.fitTo(data)

  plotPtInv = ptInv.frame()
  data.plotOn(plotPtInv)
  polyPDF.plotOn(plotPtInv)
  plotPtInv.Draw()

  for a1 in polyVarList:
    print("{0:<5} {1:<9.3g} +/- {2:<9.3g}".format(a1.GetName(),a1.getVal(),a1.getError()))

tree = inFile.Get("tree")
#findEndpoint(tree)
#silly()
silly2(tree)
saveAs(canvas,"endPoint")
