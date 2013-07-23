#!/usr/bin/env python

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

def computeGlobalSig(sigLocal,sig0,N0):
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

  print ("*"*70)
  print("Input: sigLocal: %.2f sig0: %.2f N0: %.2f " % (sigLocal,sig0,N0))
  print("local p-value corresponding to local significance: %.2e" % pLocal)
  print("p-value corresponding to significance sig0: %.2e" % p0)
  print
  print("test stat values:")
  print("  u: %.2f" % u)
  print("  u0: %.2f" % u0)
  print("  N: %.2f" % N)
  print
  print("global p-value:      %.5f" % pGlobal)
  print("global significance: %.2f" % sigGlobal)
  print("trial factor         %.2f" % trialFactor)
  print ("*"*70)
  return sigGlobal

class GlobalSigEvaluator:
  def __init__(self,globber):
    self.filenames = glob.glob(globber)
    self.sigDict = {}
    for fn in self.filenames:
      fData = self.loadToys(fn)
      match = re.match(r"[a-zA-Z]+([\d]+)\.root",fn)
      if not match:
        print("Error: Can't match filename: "+match)
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
    self.histMax.GetYaxis().SetTitle("Event")
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

    self.histMax.Scale(1./minLen)

    self.countCrossings(1.)

  def loadToys(self,filename):
    f = root.TFile(filename)
    limitTree = f.Get("limit")
    result = []
    for iEvent in range(limitTree.GetEntries()):
      limitTree.GetEntry(iEvent)
      result.append(float(limitTree.limit))
    f.Close()
    return result

  def drawHist(self):
    canvas = root.TCanvas("canvas")
    self.hist2D.Draw("colztext")
    canvas2 = root.TCanvas("canvas2")
    self.histMax.Draw()
    raw_input("press enter to continue...")

  def drawGraphs(self):
    canvas = root.TCanvas("canvas3")
    self.hist2D.Reset()
    latex = root.TLatex()
    latex.SetNDC()
    for g, i in zip(self.graphs,range(len(self.graphs))):
      print "drawing graph..."
      self.hist2D.SetTitle("Trial: %i" % i)
      self.hist2D.Draw()
      g.Draw("PL")
      canvas.Update()
      time.sleep(1)

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
    print("<nCrossings>: %.3f  for sigma %.1f" % (float(sum(nCrossList))/len(nCrossList),nSigma))
    print("<nSpentUp>: %.3f " % (float(sum(nSpentUpList))))
    return nCrossList

if __name__ == "__main__":
  gse = GlobalSigEvaluator("Add*.root")
  #gse.drawHist()
  #gse.drawGraphs()
