#!/usr/bin/env python


import optparse
parser = optparse.OptionParser(description="Makes Shape Diagnostic Plots from Datacards")
parser.add_option("--signalInject", help="Sets a caption saying that signal was injected with strength",type=float,default=20.0)
parser.add_option("-b","--binWidthOverride", help="Overrides the default bin widths and sets all binning to this widht [GeV]",type=float,default=0.0)
#parser.add_option("--plotSignalStrength", help="Plots a signal bump with this strength",type=float,default=5.0)
#parser.add_option("--plotSignalBottom", help="Plots a signal bump on the bottom (bool)",action="store_true",default=True)
#parser.add_option("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
#parser.add_option("-r","--rebinOverride", help="Rebin All plots with this rebinning, overriding all internal configuration",type=int,default=0)
args, fakeargs = parser.parse_args()

import singleHelpers
from helpers import *
from xsec import *
import math
import os.path
import glob
import random

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

PRELIMINARYSTRING="CMS"

def makePDFBakMSSM(name,rooDataset,dimuonMass,workspaceImportFn,mZ=91.19,sigmaZ=5.):
    debug = ""
    debug += "### makePDFBakMSSM: "+name+"\n"
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    bwWidth = root.RooRealVar(channelName+"_bwWidth","bwWidth",sigmaZ)
    bwmZ = root.RooRealVar(channelName+"_bwmZ","bwmZ",mZ)
    expParam = root.RooRealVar(channelName+"_expParam","expParam",-1e-03,-1e-01,1e-01)
    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0.5,0,1)

    phoExpMmumu = root.RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",root.RooArgList(dimuonMass,expParam))
    bwExpMmumu  = root.RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(dimuonMass,bwmZ,bwWidth,expParam))
    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(bwExpMmumu,phoExpMmumu),root.RooArgList(mixParam))
    

    ### Debug Time
    #frameZ = dimuonMassZ.frame()
    #frameZ.SetName("bak_PlotZ")
    #rooDatasetZ.plotOn(frameZ)
    #bwMmumuZ.plotOn(frameZ)
    #canvas = root.TCanvas()
    #frameZ.Draw()
    #canvas.SaveAs("debug_"+name+channelName+"_Z.png")

    # Back to everywhere else
   
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)
    fr.Print()

    rooParamList = [bwmZ,bwWidth,expParam,mixParam]
    paramList = [expParam.GetName(),mixParam.GetName()]

    for param in rooParamList:
      param.setConstant(True)

    ## Let them float!
    expParam.setConstant(False)
    mixParam.setConstant(False)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    ### Debug Time
    #frame = dimuonMass.frame()
    #frame.SetName("bak_Plot")
    #rooDataset.plotOn(frame)
    ##pdfMmumu.plotOn(frame,root.RooFit.Range(110,160))
    #pdfMmumu.plotOn(frame,root.RooFit.Range(minMass,maxMass))
    #canvas = root.TCanvas()
    #frame.Draw()
    #canvas.SaveAs("debug_"+name+channelName+".png")

    #Norm Time
    bakNormTup = None
    if False:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    #debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None

class ShapePlotter:
  def __init__(self,filename,outDir,titleMap,signalInject=20.,binWidthOverride=0,energyStr=None):
    self.signalInject=signalInject
    self.rmpList = []
    self.titleMap = titleMap
    self.filename = filename

    self.energyStr = energyStr
          
    self.lumi = float(lumiDict[self.energyStr])
    self.lumiStr = "L = {0:.1f} fb^{{-1}}".format(self.lumi)
      
    self.data = {}
    self.f = root.TFile(filename)
    channelNameMatch = re.search(r"mass_(cat[0123])_dist\.root",filename)
    if not channelNameMatch:
      print "fn doesn't contain channel name!!"
      sys.exit(1)
    channelName = channelNameMatch.group(1)
    channelTitle = titleMap[channelName]
    data_obs = self.f.Get("data_mass_"+channelName)
    mMuMu = rooArgSet2List(data_obs.get())[0]

    sigHistName = "th1f_sig_ggh_mass_m125_"+channelName
    sigHist = self.f.Get(sigHistName)

    binWidth = 1
    yMax = None
    if channelName == "cat0":
        binWidth *= 1.
        yMax = 4000.
    elif channelName == "cat1":
        binWidth *= 1.
        yMax = 3500.
    elif channelName == "cat2":
        binWidth *= 2.
        yMax = 30.
        sigHist.Rebin(2)
    elif channelName == "cat3":
        binWidth *= 2.
        yMax = 120.
        sigHist.Rebin(2)

    if binWidthOverride > 0:
      binWidth = binWidthOverride

    mMuMu.setRange("plotRange",110,160)
    binning = mMuMu.getBinning()
    xlow = binning.lowBound()
    xhigh = binning.highBound()
    mMuMu.setBins(int((xhigh-xlow)/binWidth))
    mMuMu.SetTitle("m_{ee} [GeV]")

    w = root.RooWorkspace("w_"+channelName)
    wImport = getattr(w,"import")
    makePDFBakMSSM(channelName,data_obs,mMuMu,wImport)

    bakPDF = w.pdf("bak")
    fr = bakPDF.fitTo(data_obs,root.RooFit.Save(),root.RooFit.Range("plotRange"),PRINTLEVEL)

    saveName = outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName
    saveName = re.sub(r"_[0-9P]+TeV_","_"+self.energyStr+"_",saveName)
    saveName = re.sub(r"([\d]+)\.[\d]+",r"\1",saveName)
    saveName = outDir+"/massPlot_"+channelName

    stupidFrame = mMuMu.frame()
    stupidData = data_obs.plotOn(stupidFrame,root.RooFit.Range("plotRange"))
    #print stupidData
    #stupidData.Print()
    #print stupidData.GetMaximum()
    #sys.exit(0)

    legEntrySignal = "SM Higgs #times 10^{6}"

    #Plot Time
    rmp = RooModelPlotter(mMuMu,bakPDF,data_obs,fr,
                          channelTitle,self.energyStr.replace("TeV"," TeV"),self.lumi,
#                          nSignal=nSignal,signalPdf=sigPDF,
                          legEntrySignal=legEntrySignal,
                          RangeName="plotRange",
                          preliminaryString=PRELIMINARYSTRING,
                          #yMax = stupidData.GetMaximum()*1.1,
                          yMax = yMax,
                          sigHist = sigHist,
                          caption1="H #rightarrow e^{+}e^{-}"
                          )
    rmp.draw(saveName)
    #rmp.drawWithParams(saveName+"_params",["mixParam","bwWidth","bwmZ","expParam"])

    self.rmpList.append(rmp)

if __name__ == "__main__":
  root.gROOT.SetBatch(True)

  outDir = "shapes/"

  fns = glob.glob("etc/heeData/mass_cat*_dist.root")

  titleMap = {
    "cat0":"0,1-jet BB",
    "cat1":"0,1-jet Not BB",
    "cat2":"2-jet Tight",
    "cat3":"2-jet Loose",
  }
  

  shapePlotterList = []
  for fn in fns:
    energyStr = "8TeV"
    print fn
    s = ShapePlotter(fn,outDir,titleMap,signalInject=args.signalInject,binWidthOverride=args.binWidthOverride,energyStr=energyStr)
    shapePlotterList.append(s)
