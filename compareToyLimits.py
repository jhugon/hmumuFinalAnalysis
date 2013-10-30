#!/usr/bin/env python

from helpers import *
import sys
import os
import glob
import re
import time

from scipy import exp
from scipy.stats import chi2
from scipy.stats import norm

import ROOT as root
from ROOT import gSystem

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

####################################################

HISTCOUNTER=0

def massFromFn(x):
  srch = re.search(r"_([0-9.]+)\.txt",x)
  assert(srch)
  return float(srch.group(1))

class LimitData(object):
  def __init__(self,filename):
    self.filename = filename
    self.mass = massFromFn(filename)
    f = root.TFile(filename)
    t = f.Get("limit")
    self.n = t.GetEntries()/6
    self.median = []
    self.p1sig = []
    self.p2sig = []
    self.m1sig = []
    self.m2sig = []
    self.obs = []
    for i in range(0,t.GetEntries(),6):
      t.GetEntry(i)
      self.m2sig.append(t.limit)
      t.GetEntry(i+1)
      self.m1sig.append(t.limit)
      t.GetEntry(i+2)
      self.median.append(t.limit)
      t.GetEntry(i+3)
      self.p1sig.append(t.limit)
      t.GetEntry(i+4)
      self.p2sig.append(t.limit)
      t.GetEntry(i+5)
      self.obs.append(t.limit)
    del t
    f.Close()

class LimitDataVMass(object):
  def __init__(self,globstr,title,color):
    self.globstr = globstr
    self.title = title
    self.color = color
    self.fns = glob.glob(globstr)
    self.fns = sorted(self.fns,key=massFromFn)
    self.masses = [massFromFn(fn) for fn in self.fns]
    self.limitData = [LimitData(fn) for fn in self.fns]
    nToys = None
    for d in self.limitData:
      if nToys:
        assert(nToys==d.n)
      else:
        nToys=d.n
    self.nToys = nToys

  def fillMedianGraph(self,graph,iToy):
    for i,ld in zip(range(len(self.limitData)),self.limitData):
      mass = ld.mass
      graph.SetPoint(i,mass,ld.median[iToy])
  def fillObsGraph(self,graph,iToy):
    for i,ld in zip(range(len(self.limitData)),self.limitData):
      mass = ld.mass
      graph.SetPoint(i,mass,ld.obs[iToy])
  def fill1sigGraph(self,graph,iToy):
    for i,ld in zip(range(len(self.limitData)),self.limitData):
      mass = ld.mass
      graph.SetPoint(i,mass,ld.median[iToy])
      graph.SetPointError(i,0.0,0.0,ld.median[iToy]-ld.m1sig[iToy],ld.p1sig[iToy]-ld.median[iToy])
  def fill2sigGraph(self,graph,iToy):
    for i,ld in zip(range(len(self.limitData)),self.limitData):
      mass = ld.mass
      graph.SetPoint(i,mass,ld.median[iToy])
      graph.SetPointError(i,0.0,0.0,ld.median[iToy]-ld.m2sig[iToy],ld.p2sig[iToy]-ld.median[iToy])

  def plotToys(self,fnBase,ylimits=[0,100.]):
    global HISTCOUNTER
    canvas = root.TCanvas()
    for iToy in range(self.nToys):
      obs = root.TGraph()
      median = root.TGraph()
      median.SetLineStyle(2)
      oneSig = root.TGraphAsymmErrors()
      oneSig.SetFillColor(root.kGreen)
      oneSig.SetLineStyle(0)
      twoSig = root.TGraphAsymmErrors()
      twoSig.SetFillColor(root.kYellow)
      twoSig.SetLineStyle(0)
      
      self.fillMedianGraph(median,iToy)
      self.fillObsGraph(obs,iToy)
      self.fill1sigGraph(oneSig,iToy)
      self.fill2sigGraph(twoSig,iToy)

      axisHist = root.TH2F("axisHist"+str(HISTCOUNTER),"",1,110,160,1,ylimits[0],ylimits[1])
      HISTCOUNTER+=1
      axisHist.GetXaxis().SetTitle("m_{H} [GeV/c^{2}]")
      axisHist.GetYaxis().SetTitle("95% CL Limit on #sigma/#sigma_{SM} (H#rightarrow#mu#mu)")

      axisHist.Draw()
      twoSig.Draw("3")
      oneSig.Draw("3")
      median.Draw("L")
      obs.Draw("LP")

      canvas.RedrawAxis()

      canvas.SaveAs(fnBase+"_toy"+str(iToy)+".png")
      canvas.SaveAs(fnBase+"_toy"+str(iToy)+".pdf")
      canvas.Clear()

def plotTogether(listOfLimitData,fnBase,showFirstBands=True,ylimits=[0,150],title=""):
    assert(len(listOfLimitData)>0)
    global HISTCOUNTER
    canvas = root.TCanvas()
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(0.04)
    for iToy in range(listOfLimitData[0].nToys):
      oneSig = root.TGraphAsymmErrors()
      oneSig.SetFillColor(root.kGreen)
      oneSig.SetLineStyle(0)
      twoSig = root.TGraphAsymmErrors()
      twoSig.SetFillColor(root.kYellow)
      twoSig.SetLineStyle(0)
      listOfLimitData[0].fill1sigGraph(oneSig,iToy)
      listOfLimitData[0].fill2sigGraph(twoSig,iToy)

      obsList = []
      medianList = []
      for ld in listOfLimitData:
        obs = root.TGraph()
        median = root.TGraph()
        median.SetLineStyle(2)
        median.SetLineColor(ld.color)
        obs.SetLineColor(ld.color)
        median.SetMarkerColor(ld.color)
        obs.SetMarkerColor(ld.color)
        ld.fillMedianGraph(median,iToy)
        ld.fillObsGraph(obs,iToy)
        obsList.append(obs)
        medianList.append(median)

      axisHist = root.TH2F("axisHist"+str(HISTCOUNTER),"",1,110,160,1,ylimits[0],ylimits[1])
      HISTCOUNTER+=1
      axisHist.GetXaxis().SetTitle("m_{H} [GeV/c^{2}]")
      axisHist.GetYaxis().SetTitle("95% CL Limit on #sigma/#sigma_{SM} (H#rightarrow#mu#mu)")

      axisHist.Draw()
      if showFirstBands:
        twoSig.Draw("3")
        oneSig.Draw("3")
      for median,obs in zip(medianList,obsList):
        median.Draw("L")
        obs.Draw("LP")

      canvas.RedrawAxis()

      #leg = root.TLegend(0.172,0.905,0.73,0.585)
      leg = root.TLegend(0.172,0.905,0.73,0.635)
      leg.SetLineColor(0)
      leg.SetFillColor(0)
      for obs,ld in zip(obsList,listOfLimitData):
        leg.AddEntry(obs,ld.title,"lp")
      leg.Draw()

      tlatex.SetTextAlign(12)
      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
      #tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Toy: {0}".format(iToy))
      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.60,"Toy: {0}".format(iToy))
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,title)

      canvas.SaveAs(fnBase+"_toy"+str(iToy)+".png")
      canvas.SaveAs(fnBase+"_toy"+str(iToy)+".pdf")
      canvas.Clear()
  

if __name__ == "__main__":
  root.gROOT.SetBatch(True)

  inDir = "/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/customToys/"
  outDir = "output/"

  limitDataBaselineVBF = LimitDataVMass(inDir+"Baseline/Jet2CutsVBFPass_8TeV_*.txt.outToys.root","Baseline, 110-160 GeV/c^{2}",1)
  limitData20BaselineVBF = LimitDataVMass(inDir+"Baseline20GeVWindow/Jet2CutsVBFPass_8TeV_*.txt.outToys.root","Baseline, 20 GeV/c^{2} Wide",4)
  limitData20BernsteinVBF = LimitDataVMass(inDir+"Bern20GeVWindow/Jet2CutsVBFPass_8TeV_*.txt.outToys.root","2-Bernstein, 20 GeV/c^{2} Wide",2)
  limitDataBaselineBB = LimitDataVMass(inDir+"Baseline/Jets01PassPtG10BB_8TeV_*.txt.outToys.root","Baseline, 110-160 GeV/c^{2}",1)
  limitData20BaselineBB = LimitDataVMass(inDir+"Baseline20GeVWindow/Jets01PassPtG10BB_8TeV_*.txt.outToys.root","Baseline, 20 GeV/c^{2} Wide",4)
  limitData20BernsteinBB = LimitDataVMass(inDir+"Bern20GeVWindow/Jets01PassPtG10BB_8TeV_*.txt.outToys.root","4-Bernstein, 20 GeV/c^{2} Wide",2)

  limitData20BernsteinVBF.plotToys(outDir+"limits_VBFTight_20Bernstein",ylimits=[0,100.])
  limitData20BaselineVBF.plotToys(outDir+"limits_VBFTight_20Baseline",ylimits=[0,100.])
  limitDataBaselineVBF.plotToys(outDir+"limits_VBFTight_Baseline",ylimits=[0,100.])
  limitData20BernsteinBB.plotToys(outDir+"limits_TightBB_20Bernstein",ylimits=[0,100.])
  limitData20BaselineBB.plotToys(outDir+"limits_TightBB_20Baseline",ylimits=[0,100.])
  limitDataBaselineBB.plotToys(outDir+"limits_TightBB_Baseline",ylimits=[0,100.])

  plotTogether([limitDataBaselineBB,limitData20BaselineBB,limitData20BernsteinBB],outDir+"compareLimits_TightBB",title="0,1-Jet Tight BB")
  plotTogether([limitDataBaselineVBF,limitData20BaselineVBF,limitData20BernsteinVBF],outDir+"compareLimits_VBFTight",title="2-Jet VBF Tight")
