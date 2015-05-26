#!/usr/bin/env python

import singleHelpers
from helpers import *
from xsec import *

def getHistogram(directory,sample,energy,catname,cutString,binSize=0.5,aroundZ=False):
    fullCutString = treeCut(catname,cutString,eventWeights=True,muonRequirements=True,KDString=None)
    tmpFLoc = directory+sample+"_"+energy+".root"
    tmpF = root.TFile(tmpFLoc)
    treename = "outtree"
    minMass = 110.
    maxMass = 160.
    if aroundZ:
      minMass = 88.
      maxMass = 94.
    tmpTree = tmpF.Get(treename)
    tmpTree.SetCacheSize(10000000);
    tmpTree.AddBranchToCache("*");
    histName = sample+energy+catname
    if aroundZ:
      histName += "_Z_"
    histName +="{0:f}".format(time.time()).replace('.','')
    nBins = int((maxMass-minMass)/binSize)
    tmpHist = root.TH1F(histName,"",nBins,minMass,maxMass)
    drawString = "dimuonMass >> {0}".format(histName)
    tmpTree.Draw(drawString,fullCutString)
    tmpHist.SetDirectory(0)
    return tmpHist

if __name__ == "__main__":
  root.gROOT.SetBatch(True)

  outdir = "output/"
  indir = getDataStage2Directory()
  
  mcSetsDict = {
    "DYJetsToLL":"Drell-Yan",
    "ttbar":"t#bar{t}",
  }
  
  periods = ["7TeV","8TeV"]
  periods = ["8TeV"]
  #periods = ["7TeV"]
  categoriesAll = ["BB","BO","BE","OO","OE","EE"]
  
  jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
  jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

  categories = []
  
  categories += [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts]]
  ####categories += [["Jets01FailPtG10BO",  "dimuonPt>10." +jet01PtCuts]]
  
  #categories += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #categories += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  #categories += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  #categories += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  #categories += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  
  canvas = root.TCanvas("canvas")
  for energy in periods:
    for category in categories:
      binSize = 2.0
      hists = []
      hstack = root.THStack()
      for sample in ["ttbar","DYJetsToLL"]:
        hist = getHistogram(indir,sample,energy,category[0],category[1],binSize=binSize)
        sample_fullname = sample+"_"+energy
        scaleFactor = xsec[sample_fullname]*1000.0*lumiDict[energy]/nEventsMap[sample_fullname]*efficiencyMap[energy]*mcPlotScaleFactorMap[energy]
        hist.Scale(scaleFactor)
        hist.SetFillColor(colors[sample])
        hist.SetLineColor(colors[sample])
        hstack.Add(hist)
        hists.append(hist)
      canvas.Clear()
      ymax = 2000
      axisHist = root.TH2F("axisHist"+energy+category[0],"",1,110,160,1,0,ymax)
      setHistTitles(axisHist,"m_{#mu#mu} [GeV]","Entries/{0:0.1f} GeV".format(binSize))
      axisHist.GetXaxis().CenterTitle(True)
      axisHist.Draw()
      hstack.Draw("same")
      canvas.RedrawAxis()
      energyStr = energy
      if re.search(r"[\d]TeV",energyStr):
        energyStr = energyStr.replace("TeV"," TeV")
      lumiStr = "{0:.1f} fb^{{-1}} ({1})".format(lumiDict[energy],energyStr)
      drawStandardCaptions(canvas,lumiStr,TITLEMAP[category[0]],preliminaryString="")
      canvas.SaveAs("test.png")
