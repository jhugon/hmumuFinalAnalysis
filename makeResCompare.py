#!/usr/bin/env python

from helpers import *
import math
import os.path
import glob
import numpy
import matplotlib.pyplot as mpl

root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

#from ROOT import gSystem
#gSystem.Load('libRooFit')

from makeCards import makePDFSig

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

class ResCompare:
  def __init__(self,fileNames,titles,categories):
    assert(len(fileNames) == len(titles))
    self.files = []
    self.histsCatDict = {}
    self.titles = titles
    self.nFiles = len(fileNames)
    self.categories = categories

    self.colors = [1,root.kRed+1,root.kBlue+1,root.kGreen+1,root.kCyan,root.kPink]

    for i in categories:
      self.histsCatDict[i] = []

    for f in fileNames:
      tmpFile = root.TFile(f)
      self.files.append(tmpFile)
      for i in categories:
        strToGet = i + '/mDiMu'
        strToGet = os.path.normpath(strToGet)
        if strToGet[0] == '/':
            strToGet = strToGet[1:]
        tmpHist = tmpFile.Get(strToGet)
        self.histsCatDict[i].append(tmpHist)
    
    minMass = 100
    maxMass = 150

    self.workspace = root.RooWorkspace("w")
    wImport = getattr(self.workspace,"import")
    mMuMu = root.RooRealVar("mMuMu","m_{#mu#mu} [GeV]",minMass,maxMass)
    wImport(mMuMu)
    self.mMuMu = mMuMu
    for i in categories:
     for j in range(self.nFiles):
       makePDFSig(i+str(j),self.histsCatDict[i][j],mMuMu,minMass,maxMass,wImport,i)

    #self.workspace.Print()

    self.canvas = root.TCanvas("canvas")
    
  def plotData(self,saveNameBase):
    datasets = []
    for i in self.categories:
      self.canvas.Clear()
      frame = self.mMuMu.frame()
      frame.SetTitle("")
      frame.SetYTitle("Signal MC Events")
      urLegendPos = [0.70,0.67,0.9,0.9]
      leg = root.TLegend(*urLegendPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      for j in range(self.nFiles):
        tmpDataset = self.workspace.data(i+str(j)+"_Template")
        rooLCol = root.RooFit.LineColor(self.colors[j])
        rooMCol = root.RooFit.MarkerColor(self.colors[j])
        rooNameStr = "Curve_"+i+str(j)+"_Template"
        rooName = root.RooFit.Name(rooNameStr)
        tmpDataset.plotOn(frame,rooLCol,rooMCol,rooName)
        tmpDatasetH = frame.getHist(rooNameStr)
        leg.AddEntry(tmpDatasetH,self.titles[j],"p")
        datasets.append(tmpDataset)
      frame.Draw()
      leg.Draw()
      saveAs(self.canvas,saveNameBase+i)

  def plotPDF(self,saveNameBase):
    curves = []
    for i in self.categories:
      self.canvas.Clear()
      frame = self.mMuMu.frame()
      frame.SetTitle("")
      frame.SetYTitle("Signal MC Events")
      urLegendPos = [0.70,0.67,0.9,0.9]
      leg = root.TLegend(*urLegendPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      for j in range(self.nFiles):
        tmpCurve = self.workspace.pdf(i+str(j))
        rooLCol = root.RooFit.LineColor(self.colors[j])
        rooNameStr = "Curve_"+i+str(j)
        rooName = root.RooFit.Name(rooNameStr)
        tmpCurve.plotOn(frame,rooLCol,rooName)
        tmpCurveH = frame.getCurve(rooNameStr)
        leg.AddEntry(tmpCurveH,self.titles[j],"l")
        curves.append(tmpCurve)
      frame.Draw()
      leg.Draw()
      saveAs(self.canvas,saveNameBase+i)

if __name__ == "__main__":

  categories = ["IncPresel","VBFPresel","IncBDTCut","VBFBDTCut"]

  infiles = []
  titles = []
  infiles.append("input/ggHmumu125_8TeV.root")
  infiles.append("input/smearing/ggHmumu125_8TeV.root")
  infiles.append("input/rochester/ggHmumu125_8TeV.root")
  infiles.append("input/muscle/ggHmumu125_8TeV.root")
  titles.append("Default")
  titles.append("Smearing")
  titles.append("Rochester")
  titles.append("Muscle")

  #infiles.append("input/vbfHmumu125_8TeV.root")
  #infiles.append("input/smearing/vbfHmumu125_8TeV.root")
  #infiles.append("input/rochester/vbfHmumu125_8TeV.root")
  #infiles.append("input/muscle/vbfHmumu125_8TeV.root")

  rs = ResCompare(infiles,titles,categories)
  rs.plotData("resData")
  rs.plotPDF("resShape")
