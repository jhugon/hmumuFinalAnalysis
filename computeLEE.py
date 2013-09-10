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
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

####################################################

def convertN(sig,sig0,N0):
  # Convert significance to p-value
  p = norm.sf(sig)
  p0 = norm.sf(sig0)
  # Get the test statistic value corresponding to the p-value
  u = chi2.isf(p*2,1)
  u0 = chi2.isf(p0*2,1)
  # The main equation
  N = N0 * exp(-(u-u0)/2.)
  return N

class SigGlobalObj:
  def __init__(self,sigLocal,sig0,N0):
    # Convert significance to p-value
    pLocal = norm.sf(sigLocal)
    p0 = norm.sf(sig0)
    
    # Get the test statistic value corresponding to the p-value
    u = chi2.isf(pLocal*2,1)
    u0 = chi2.isf(p0*2,1)
    
    # The main equations
    N = N0 * exp(-(u-u0)/2.)
    pGlobal = N + chi2.sf(u,1)/2.
    
    # Further info
    sigGlobal = norm.isf(pGlobal)
    trialFactor = pGlobal/pLocal

    self.sigGlobal = sigGlobal
    self.sigLocal = sigLocal
    self.sig0 = sig0
    self.pGlobal = pGlobal
    self.pLocal = pLocal
    self.p0 = p0
    self.N0 = N0
    self.N = N
    self.u0 = u0
    self.u = u
    self.trialFactor = trialFactor

  def __call__(self):
    return self.sigGlobal

  def __str__(self):
    result = ""
    result += "*"*70 + "\n"
    result += "Input: sigLocal: %.2f sig0: %.2f N0: %.2f \n" % (self.sigLocal,self.sig0,self.N0)
    result += "local p-value corresponding to local significance: %.2e\n" % self.pLocal
    result += "p-value corresponding to significance sig0: %.2e\n" % self.p0
    result += "\n"
    result += "test stat values:\n"
    result += "  u: %.2f\n" % self.u
    result += "  u0: %.2f\n" % self.u0
    result += "  N: %.2f\n" % self.N
    result += "\n"
    result += "global p-value:      %.5f\n" % self.pGlobal
    result += "global significance: %.2f\n" % self.sigGlobal
    result += "trial factor         %.2f\n" % self.trialFactor
    result += "*"*70 + "\n"
    return result

class GlobalSigEvaluator:
  def __init__(self,globber):
    minData = 5
    self.filenames = glob.glob(globber)
    self.sigDict = {}
    for fn in self.filenames:
      fData = self.loadToys(fn)
      match = re.match(r"[a-zA-Z]+([\d]+)\.root",fn)
      if not match:
        print("Error: Can't match filename: "+match)
        continue
      if len(fData) < minData:
        print("Warning: Throwing away small nData: %i for %s" % (len(fData),fn))
        continue
      massStr = match.group(1)
      self.sigDict[massStr] = fData
    minLen = 1e6
    maxLen = 0.
    for key in self.sigDict:
      minLen = min(minLen,len(self.sigDict[key]))
      maxLen = max(maxLen,len(self.sigDict[key]))
    print("Have %.0f complete sets out of maxset: %.0f" % (minLen,maxLen))
    self.graphs = []
    self.histMax = root.TH1F("sigMax","",50,0.,5.)
    self.histMax.GetXaxis().SetTitle("Max Significance")
    self.histMax.GetYaxis().SetTitle("Fraction of Events")
    self.listMax = []
    self.hist2D = root.TH2F("sigHist","",50,110,160,50,0.,5.)
    self.hist2D.GetXaxis().SetTitle("m_{H} [GeV/c^{2}]")
    self.hist2D.GetYaxis().SetTitle("Significance")
    dictKeys = sorted(self.sigDict.keys(),key=float)
    self.trialList = []
    for i in range(minLen):
      sigListVMass = []
      graph = root.TGraph()
      graph.SetMarkerColor(1)
      graph.SetMarkerSize(1)
      graph.SetMarkerStyle(8)
      for key,iGraph in zip(dictKeys,range(len(dictKeys))):
        mass = float(key)
        sig = self.sigDict[key][i]
        self.hist2D.Fill(mass,sig)
        graph.SetPoint(iGraph,mass,sig)
        sigListVMass.append(sig)
      self.graphs.append(graph)
      self.trialList.append(sigListVMass)
      self.histMax.Fill(max(sigListVMass))
      self.listMax.append(max(sigListVMass))

    self.histMax.Scale(1./minLen)

    avgN1sig, avgN1sigErr = self.countCrossings(1.)
    sigList = [0.5,1.,1.5,2.,2.5,3.,3.5]
    avgNList = []
    avgNErrList = []
    predAvgNList = []
    predAvgNErrList = []
    for sig in sigList:
      avgN, avgNErr = self.countCrossings(sig)
      avgNList.append(avgN)
      avgNErrList.append(avgNErr)
      predAvgN = convertN(sig,1.0,avgN1sig)
      predAvgNErr = avgN1sigErr*predAvgN/avgN1sig
      predAvgNList.append(predAvgN)
      predAvgNErrList.append(predAvgNErr)

    graphNAvgMeas = root.TGraphErrors()
    graphNAvgPred = root.TGraphErrors()
    graphNAvgPred.SetLineColor(root.kRed)
    graphNAvgPred.SetMarkerColor(root.kRed)
    setHistTitles(graphNAvgPred,"Significance","<N>")
    setHistTitles(graphNAvgMeas,"Significance","<N>")
    self.graphNAvgMeas = graphNAvgMeas
    self.graphNAvgPred = graphNAvgPred
    print
    print("From 1 sigma <N> = %.3f +/- %.3f:" % (avgN1sig,avgN1sigErr))
    iPoint=0
    for sig,pred,meas,predErr,measErr in zip(sigList,predAvgNList,avgNList,predAvgNErrList,avgNErrList):
      print("  %.2f sigma predict <N>: %.3f +/- %.3f, measure <N>: %.3f +/- %.3f" % (sig,pred,predErr,meas,measErr))
      graphNAvgPred.SetPoint(iPoint,sig+0.05,pred)
      graphNAvgPred.SetPointError(iPoint,0.,predErr)
      if meas != 0.:
        graphNAvgMeas.SetPoint(iPoint,sig,meas)
        graphNAvgMeas.SetPointError(iPoint,0.,measErr)
      iPoint += 1

  def loadToys(self,filename):
    f = root.TFile(filename)
    limitTree = f.Get("limit")
    result = []
    if limitTree:
      for iEvent in range(limitTree.GetEntries()):
        limitTree.GetEntry(iEvent)
        result.append(float(limitTree.limit))
    f.Close()
    return result

  def drawHist(self):
    canvas = root.TCanvas("canvas")
    wh = canvas.GetWindowHeight()
    ww = canvas.GetWindowWidth()
    ww *= 1.5
    wh *= 1.5
    canvas.SetWindowSize(int(ww),int(wh))
    canvas.Divide(2,2)
    canvas.cd(1)
    self.hist2D.Draw("colztext")
    canvas.cd(2)
    self.histMax.Draw()
    canvas.cd(3)
    self.graphNAvgPred.Draw("APE")
    self.graphNAvgMeas.Draw("PE")
    canvas.Update()
    raw_input("press enter to continue...")

  def drawGraphs(self):
    canvas = root.TCanvas("canvas3")
    self.hist2D.Reset()
    for g, i in zip(self.graphs,range(len(self.graphs))):
      print "drawing graph..."
      self.hist2D.SetTitle("Trial: %i" % i)
      self.hist2D.Draw()
      g.Draw("PL")
      canvas.Update()
      time.sleep(0.5)

  def drawMassCompare(self):
    canvas = root.TCanvas("canvas4")
    wh = canvas.GetWindowHeight()
    ww = canvas.GetWindowWidth()
    ww *= 2
    canvas.SetWindowSize(ww,wh)
    canvas.Divide(2)
    pad1 = canvas.cd(1)
    pad2 = canvas.cd(2)
    self.hist2D.Reset()
    for g, i in zip(self.graphs,range(len(self.graphs))):
      print "drawing graphs..."
      pad1.cd()
      self.hist2D.SetTitle("Trial: %i" % i)
      self.hist2D.Draw()
      g.Draw("PL")
    
      pad2.cd()
      #rmp = RooModelPlotter(xVar,pdf,data,fr,"Trial: %i" % i,"8TeV",99,canvas=pad2)
      #rmp.draw("/tmp/stupid")

      canvas.Update()
      time.sleep(0.5)


  def countCrossings(self,nSigma):
    nCrossList = []
    nSpentUpList = []
    for trial in self.trialList:
      up = True
      nCrossings = 0
      nSpentUp =0
      for sig in trial:
        if sig >= nSigma and not up:
          up = True
          nCrossings += 1
        elif sig < nSigma and up:
          up = False
        elif up:
          nSpentUp
      nCrossList.append(nCrossings)
      nSpentUpList.append(nSpentUp)
    avgN = float(sum(nCrossList))/len(nCrossList)
    stdDevN = 0.
    for i in nCrossList:
      stdDevN += (i-avgN)**2
    stdDevN /= len(nCrossList)-1.5
    stdDevN = sqrt(stdDevN)
    avgNErr = stdDevN/sqrt(len(nCrossList))
    return avgN, stdDevN

  def getSigGlobal(self,sigLocal,sig0):
    nAvg, nAvgErr = self.countCrossings(sig0)
    sg = SigGlobalObj(sigLocal,sig0,nAvg)
    return sg

  def getSigGlobalTraditional(self,sigLocal):
    result = 0.
    for i in self.listMax:
      if sigLocal <= i:
        result += 1.
    result /= len(self.listMax)
    return result


if __name__ == "__main__":
  os.chdir("playWithLEE/")
  gse = GlobalSigEvaluator("Add*.root")
  gse.drawHist()
  #gse.drawGraphs()
  sigLocal = 2.27435
  print "for local significance: %.3f" % sigLocal
  print "global significance is:"
  for sig0 in [0.25,0.5,1.,1.5,2.]:
    print "  for sig0: %10.3f sigGlobal: %10.3f" % (sig0,gse.getSigGlobal(sigLocal,sig0)())
  print gse.getSigGlobal(sigLocal,1.0)
  print "traditional global p-value (may have large stat errror: ) %10.3f" % gse.getSigGlobalTraditional(sigLocal)
