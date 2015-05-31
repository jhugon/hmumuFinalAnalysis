#!/usr/bin/env python

import singleHelpers
from helpers import *
from xsec import *
from makeCards import makePDFBak

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

def getYMaxAuto(plotHist,logy=False):
  yTopOfFrame = 0.
  frameSF = 1.2
  frameRightSideSF = 1.5
  if logy:
    frameSF = 1.5
    frameRightSideSF = 5.

  nBins = plotHist.GetNbinsX()
  for iBin in range(1,nBins+1):
    yMaxThisBin = plotHist.GetBinContent(iBin)+plotHist.GetBinError(iBin)
    if logy:
      yTopOfFrameTry = yMaxThisBin**frameSF
      yTopOfFrame = max(yTopOfFrame,yTopOfFrameTry)
    else:
      yTopOfFrameTry = frameSF*yMaxThisBin
      yTopOfFrame = max(yTopOfFrame,yTopOfFrameTry)
    
  return yTopOfFrame

def setupHistogram(hist,sample,energy):
  sample_fullname = sample+"_"+energy
  scaleFactor = xsec[sample_fullname]*1000.0*lumiDict[energy]/nEventsMap[sample_fullname]*efficiencyMap[energy]*mcPlotScaleFactorMap[energy]
  hist.Scale(scaleFactor)
  hist.SetFillColor(colors[sample])
  hist.SetLineColor(colors[sample])

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
  minMass = 110.
  maxMass = 160.
  dimuonMass = root.RooRealVar("dimuonMass","dimuonMass",minMass,maxMass)
  controlRegionVeryLow=[60,110]
  controlRegionLow=[110,120]
  controlRegionHigh=[130,160]
  SIGNALFIT = [110.,140.]
  dimuonMass.setRange("exprange",120.,controlRegionHigh[1])
  dimuonMass.setRange("whole",controlRegionLow[0],controlRegionHigh[1])
  dimuonMass.setRange("low",controlRegionLow[0],controlRegionLow[1])
  dimuonMass.setRange("high",controlRegionHigh[0],controlRegionHigh[1])
  dimuonMass.setRange("signal",controlRegionLow[1],controlRegionHigh[0])
  dimuonMass.setRange("signalfit",SIGNALFIT[0],SIGNALFIT[1])
  minMassZ = 88.
  maxMassZ = 94.
  dimuonMassZ = root.RooRealVar("dimuonMass","dimuonMass",minMassZ,maxMassZ)
  for energy in periods:
    for category in categories:
      binSize = 2.0
      hists = []
      hstack = root.THStack()
      hsum = None
      zhists = []
      zhsum = None
      for sample in ["ttbar","DYJetsToLL"]:
        hist = getHistogram(indir,sample,energy,category[0],category[1],binSize=binSize)
        setupHistogram(hist,sample,energy)
        hstack.Add(hist)
        hists.append(hist)
        if hsum == None:
          hsum = hist.Clone(energy+category[0]+"sumHist")
        else:
          hsum.Add(hist)
        # now for zhist
        zhist = getHistogram(indir,sample,energy,category[0],category[1],binSize=binSize,aroundZ=True)
        setupHistogram(zhist,sample,energy)
        zhists.append(zhist)
        if zhsum == None:
          zhsum = zhist.Clone(energy+category[0]+"zsumHist")
        else:
          zhsum.Add(hist)

      workspace = root.RooWorkspace(energy+category[0]+"workspace")
      wImport = getattr(workspace,"import")
      rooHistName = "bkg"
      zrooHistName = rooHistName+"Z"
      rooHist = root.RooDataHist(rooHistName,rooHistName,root.RooArgList(dimuonMass),hsum)
      zrooHist = root.RooDataHist(zrooHistName,zrooHistName,root.RooArgList(dimuonMassZ),zhsum)
      makePDFBak(category[0],rooHist,dimuonMass,110,160,wImport,dimuonMassZ,zrooHist)
      pdf = workspace.pdf("bak")

      canvas.Clear()
      ymax = getYMaxAuto(hsum)
      axisHist = root.TH2F("axisHist"+energy+category[0],"",1,110,160,1,0,ymax)
      setHistTitles(axisHist,"m_{#mu#mu} [GeV]","Entries/{0:0.1f} GeV".format(binSize))
      axisHist.GetXaxis().CenterTitle(True)
      axisHist.Draw()
      hstack.Draw("same")

      frame = dimuonMass.frame()
      rooHist.plotOn(frame,root.RooFit.Invisible())
      pdf.plotOn(frame,root.RooFit.LineColor(1))
      frame.Draw("same")

      #legPos = [gStyle.GetPadLeftMargin()+0.03,0.55,0.68,0.93-gStyle.GetPadTopMargin()]
      #leg = root.TLegend(*legPos)
      #leg.SetFillColor(0)
      #leg.SetLineColor(0)
      #leg.AddEntry(expGraph,"Median Expected Limit","lp")
      #leg.AddEntry(oneSigGraph,"#pm1 #sigma Expected Limit","f")
      #leg.AddEntry(twoSigGraph,"#pm2 #sigma Expected Limit","f")
      #self.legPos = legPos
      #self.leg = leg
      #leg.Draw()

      canvas.RedrawAxis()
      energyStr = energy
      if re.search(r"[\d]TeV",energyStr):
        energyStr = energyStr.replace("TeV"," TeV")
      lumiStr = "{0:.1f} fb^{{-1}} ({1})".format(lumiDict[energy],energyStr)
      drawStandardCaptions(canvas,lumiStr,TITLEMAP[category[0]],preliminaryString="Simulated Background Data")
      canvas.SaveAs(outdir+"bkgTTCheck_"+category[0]+energy+".png")

      #####################################################3
      #####################################################3
      origHistTotals = []
      for orighist in hists:
        origHistTotals.append(orighist.Integral(1,orighist.GetNbinsX()))
      total = sum(origHistTotals)
      origHistFractions = [float(x)/total for x in origHistTotals]

      for frac in [0.2,0.5,0.9]:
        #print "Justin: frac:         ",frac
        #print "Justin: total:        ",total
        #print "Justin: ttbar before: ",hists[0].Integral(1,hists[0].GetNbinsX())
        #print "Justin: DY before:    ",hists[1].Integral(1,hists[1].GetNbinsX())
        hstackfrac = root.THStack()
        hsumfrac = None
        for iorig, orighist in enumerate(hists):
          hist = orighist.Clone(energy+category[0]+"frac{0}".format(frac))
          if iorig == 0: # ttbar
            #print "Justin: ttbar before: ",hist.Integral(1,hist.GetNbinsX())
            hist.Scale(float(frac)/origHistFractions[iorig])
            #print "Justin: ttbar after:  ",hist.Integral(1,hist.GetNbinsX())
          else: # others
            #print "Justin: DY before:    ",hist.Integral(1,hist.GetNbinsX())
            hist.Scale(float(1.0-frac)/origHistFractions[iorig])
            #print "Justin: DY after:     ",hist.Integral(1,hist.GetNbinsX())
          hstackfrac.Add(hist)
          if hsumfrac == None:
            hsumfrac = hist.Clone(energy+category[0]+"frac{0}sumHist".format(frac))
          else:
            hsumfrac.Add(hist)
        rooHistFrac = root.RooDataHist(rooHistName+"frac{0}sumHist".format(frac),rooHistName+"frac{0}sumHist".format(frac),root.RooArgList(dimuonMass),hsumfrac)
        pdf.fitTo(rooHistFrac)

        canvas.Clear()
        axisHist.Draw()
        hstackfrac.Draw("same")

        frameFrac = dimuonMass.frame()
        rooHistFrac.plotOn(frameFrac,root.RooFit.Invisible())
        pdf.plotOn(frameFrac,root.RooFit.LineColor(1))
        frameFrac.Draw("same")

        #leg.Draw()

        canvas.RedrawAxis()
        energyStr = energy
        if re.search(r"[\d]TeV",energyStr):
          energyStr = energyStr.replace("TeV"," TeV")
        lumiStr = "{0:.1f} fb^{{-1}} ({1})".format(lumiDict[energy],energyStr)
        drawStandardCaptions(canvas,lumiStr,TITLEMAP[category[0]],preliminaryString="Simulated Background Data")
        canvas.SaveAs(outdir+"bkgTTCheck_"+category[0]+energy+"_TTPerc{0:.0f}".format(frac*100)+".png")
