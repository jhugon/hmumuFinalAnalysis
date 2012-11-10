#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import os.path
import sys
import glob

def getSigLimits(fn, histName, quantile=0.9999):
  f = root.TFile(fn)
  h = f.Get(histName)
  l = getMedianAndQuantileInterval(h,quantile)
  f.Close()
  return l[2],l[0]

def makeSillyHist(inHist,outHist,cutVal):
      cutBin = inHist.GetYaxis().FindBin(cutVal)
      nBinsY = inHist.GetYaxis().GetNbins()
      nBinsX = inHist.GetXaxis().GetNbins()
      #print("max bin: {0}, cutBin: {1}".format(nBinsY,cutBin))
      for iX in range(0,nBinsX+2):
        mySum = 0.
        for iY in range(cutBin,nBinsY+2):
          mySum += inHist.GetBinContent(iX,iY)
        outHist.SetBinContent(iX,mySum)

def makeSillyFiles(inDir,outDir, runperiod,
            rebinY=1,step=0.01
    ):
  incCuts = set()
  vbfCuts = set()
  incLimits = getSigLimits(dataDir+incSigFilename,"BDTHistMuonOnly")
  vbfLimits = getSigLimits(dataDir+vbfSigFilename,"BDTHistVBF")
  #print("Inc Limits: {0} \tVBF Limits: {1}".format(incLimits,vbfLimits))
  incHistname = "BDTHistMuonOnlyVMass"
  vbfHistname = "BDTHistVBFVMass"
  filenames = glob.glob(dataDir+"/*_"+RUNPERIOD+".root")
  for fn in filenames:
    infile = root.TFile(fn)
    inc = infile.Get(incHistname)
    vbf = infile.Get(vbfHistname)
    inc.Rebin2D(1,rebinY)
    vbf.Rebin2D(1,rebinY)
    out = infile.Get("mDiMu")
  
    #outFileName = os.path.splitext(fn)[0]
    outFileName = os.path.basename(fn)
    outFileName = outDir +"/"+ outFileName
    outFileName = os.path.normpath(outFileName)
    print(outFileName)
    outFile = root.TFile(outFileName,"RECREATE")
    outFile.cd()
  
    iCut = 0
    while True:
      incCutVal = incLimits[0]+iCut*step
      if incCutVal > incLimits[1]:
        break
      incCutStr = "{0:.2f}".format(incCutVal)
      if incCutStr not in incCuts:
        incCuts.add(incCutStr)
      outInc = out.Clone("mDiMuBDTIncCut"+incCutStr)
      outInc.Reset()
      makeSillyHist(inc,outInc,incCutVal)
      outInc.Write()
      iCut += 1
    iCut = 0
    while True:
      vbfCutVal = vbfLimits[0]+iCut*step
      if vbfCutVal > vbfLimits[1]:
        break
      vbfCutStr = "{0:.2f}".format(vbfCutVal)
      if vbfCutStr not in vbfCuts:
        vbfCuts.add(vbfCutStr)
      outVBF = out.Clone("mDiMuBDTVBFCut"+vbfCutStr)
      outVBF.Reset()
      makeSillyHist(vbf,outVBF,vbfCutVal)
      outVBF.Write()
      iCut += 1
    outFile.Close()

  return incCuts, vbfCuts

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description="Creates Limit v. BDT Cut Plots For Analysis")
  parser.add_argument("-p","--plotOnly", help="Only Re-does the plotting stage of processing",action="store_true",default=False)
  args = parser.parse_args()

  root.gROOT.SetBatch(True)
  NPROCS = 1

  periods = ["7TeV","8TeV"]
  lumiDict = {}
  lumiDict["7TeV"] = 5.0
  lumiDict["8TeV"] = 20.0
  

  dataDir = "input/"
  outDir = "bdtCutInputs/"
  incCutDict = {}
  vbfCutDict = {}
  if not args.plotOnly:
    for RUNPERIOD in periods:
      incSigFilename = "ggHmumu125_"+RUNPERIOD+".root"
      vbfSigFilename = "vbfHmumu125_"+RUNPERIOD+".root"
      incCuts, vbfCuts = makeSillyFiles(dataDir,outDir,RUNPERIOD)
      incCutDict[RUNPERIOD] = incCuts
      vbfCutDict[RUNPERIOD] = vbfCuts
  
    ############################################################
    print "Starting makeCards Portion..."
    sys.stdout.flush()
  
    from makeCards import *
  
    signalNames=["ggHmumu125","vbfHmumu125","wHmumu125","zHmumu125"]
    backgroundNames= ["DYJetsToLL","ttbar","WW","WZ","ZZ"]
  
    dataDict = {}
    dataDict["8TeV"] = [
      #"SingleMuRun2012Av1",
      #"SingleMuRun2012Bv1",
      #"SingleMuRun2012Cv1",
      #"SingleMuRun2012Cv2"
    ]
    dataDict["7TeV"] = [
      #"SingleMuRun2011Av1",
      #"SingleMuRun2011Bv1"
    ]
    
    MassRebin = 1 # 4 Bins per GeV originally
    controlRegionVeryLow=[80,110]
    controlRegionLow=[110,120]
    controlRegionHigh=[130,160]
    shape=True
    toyData=False
    print("Creating Threads...")
    threads = []
    for p in periods:
      for cutStr in incCutDict[p]:
        threads.append(
          ThreadedCardMaker(
            #__init__ args:
            outDir,["mDiMuBDTIncCut"],
            appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
            rebin=[MassRebin], bakShape=shape,
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=cutStr,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,
            #write args:
            outfilename=outDir+"BDTIncCut"+"_"+p+"_"+cutStr+".txt",lumi=lumiDict[p]
          )
        )
      for cutStr in vbfCutDict[p]:
        threads.append(
          ThreadedCardMaker(
            #__init__ args:
            outDir,["mDiMuBDTVBFCut"],
            appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
            rebin=[MassRebin], bakShape=shape,
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=cutStr,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,
            #write args:
            outfilename=outDir+"BDTVBFCut"+"_"+p+"_"+cutStr+".txt",lumi=lumiDict[p]
          )
        )
    nThreads = len(threads)
    print("nProcs: {0}".format(NPROCS))
    print("nCards: {0}".format(nThreads))
  
    threadsNotStarted = copy.copy(threads)
    threadsRunning = []
    threadsDone = []
    while True:
      iThread = 0
      while iThread < len(threadsRunning):
          alive = threadsRunning[iThread].is_alive()
          if not alive:
            tmp = threadsRunning.pop(iThread)
            threadsDone.append(tmp)
          else:
            iThread += 1
  
      nRunning = len(threadsRunning)
      if nRunning < NPROCS and len(threadsNotStarted) > 0:
          tmp = threadsNotStarted.pop()
          tmp.start()
          threadsRunning.append(tmp)
  
      nRunning = len(threadsRunning)
      nNotStarted = len(threadsNotStarted)
      if nRunning == 0 and nNotStarted == 0:
          break
  
      time.sleep(0.1)
  
    ############################################################
    print "Starting combine tool Portion..."
    sys.stdout.flush()
  
    import os
    import subprocess
  
    os.chdir(outDir)
    print(os.getcwd())
  
    for f in glob.glob("*.txt"):
      outfilename = f + ".out"
      outfile = open(outfilename,"w")
      print("Running combine on {0}".format(f))
      sys.stdout.flush()
      subprocess.call(["combine","-M","Asymptotic",f],stdout=outfile,stderr=outfile)

    os.chdir("..")

  ############################################################
  print "Starting Limits Plots Portion..."
  sys.stdout.flush()

  os.chdir(outDir)

  titleMap={
    "BDTIncCut":"Inclusive",
    "BDTVBFCut":"VBF"
  }

  import makeLimitPlots as limits

  setStyle()
  canvas = root.TCanvas()
  #canvas.SetLogx(1)
  #canvas.SetLogy(1)
  ylimits=[0.0,30.0]

  for period in ["7TeV","8TeV"]:
    allfiles = glob.glob("*_"+period+"_*.txt.out")
    
    ## Limit v. Lumi
    energyStr = ""
    plots = set()
    for fn in allfiles:
      match = re.search(r"(.+)_(.+)_[-\d.]+\.txt\.out",fn)
      badPlot = re.search(r"Silly",fn)
      badPlot2 = re.search(r"Silly",fn)
      if match and not (badPlot or badPlot2):
        plots.add(match.group(1))
        energyStr = match.group(2)

    caption2 = "#sqrt{s}="+energyStr
    caption3 = "L={0:.1f} fb^{{-1}}".format(lumiDict[period])
    legend = root.TLegend(0.58,0.70,0.9,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    for plotName in plots:
      data = limits.getData(plotName+"_"+energyStr+"_*.txt.out")
      if len(data)<=1:
        continue
      incPlot = limits.RelativePlot(data,canvas,legend,titleMap[plotName],caption2=caption2,caption3=caption3,ylimits=ylimits,energyStr=energyStr,xlabel="BDT Cut Value")
      saveAs(canvas,plotName+"_"+energyStr)

  print("Done.")
