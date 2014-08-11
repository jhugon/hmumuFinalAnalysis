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
    self.channelName = channelName
    channelTitle = titleMap[channelName]
    data_obs = self.f.Get("data_mass_"+channelName)
    mMuMu = rooArgSet2List(data_obs.get())[0]

    sigHistName = "th1f_sig_ggh_mass_m125_"+channelName
    sigHist = self.f.Get(sigHistName)

    w = root.RooWorkspace("w_"+channelName)
    self.w = w
    wImport = getattr(w,"import")

    # Make a signal PDF
    sigPDFName = self.makeSigPdf(mMuMu,sigHist)
    sigPDF = w.pdf(sigPDFName)

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
                          #sigHist = sigHist,
                          nSignal=2*sigHist.Integral(),  # The 2 is a hack to make it work!!
                          signalPdf=sigPDF,
                          caption1="H #rightarrow e^{+}e^{-}"
                          )
    rmp.draw(saveName)
    #rmp.drawWithParams(saveName+"_params",["mixParam","bwWidth","bwmZ","expParam"])

    self.rmpList.append(rmp)

    self.nSignal = sigHist.Integral()
    self.fwhm = calcFWHM(sigPDF,mMuMu,110,160,0.02)
    # Do table work
    obsVarSet = root.RooArgSet(mMuMu)
    fwhmRangeName = "myIntRangeforFWHM_{0}".format(channelName)
    mMuMu.setRange(fwhmRangeName,125-0.5*self.fwhm,125.+0.5*self.fwhm) 
    pdfFrac = bakPDF.createIntegral(obsVarSet,obsVarSet,fwhmRangeName)
    pdfFracNorm = bakPDF.createIntegral(obsVarSet,obsVarSet,"plotRange")
    self.nBkg = pdfFrac.getVal()/pdfFracNorm.getVal() * data_obs.sumEntries("{0} > {1} && {0} < {2}".format(mMuMu.GetName(),110,160))
    self.nData = data_obs.sumEntries("{0} > {1} && {0} < {2}".format(mMuMu.GetName(),125-0.5*self.fwhm,125.+0.5*self.fwhm))

  def makeSigPdf(self,dimuonMass,sigHist):
    #self.wSig = root.RooWorkspace("wSig_"+self.channelName)
    self.wSig = self.w
    sigData = root.RooDataHist("datSig_"+self.channelName,"",root.RooArgList(dimuonMass),sigHist,1.)

    channelName = self.channelName
    name = ""
    meanG1 = root.RooRealVar(channelName+"_"+name+"_MeanG1",
                             channelName+"_"+name+"_MeanG1", 
                             #124.5)
                             124.5,115.,140.)
    meanG2 = root.RooRealVar(channelName+"_"+name+"_MeanG2",
                             channelName+"_"+name+"_MeanG2", 
                             #123.0)
                             123.,115.,140.)
    
    widthG1 = root.RooRealVar(channelName+"_"+name+"_WidthG1",
                              channelName+"_"+name+"_WidthG1", 
                             1.,0.1,15.)
    widthG2 = root.RooRealVar(channelName+"_"+name+"_WidthG2",
                              channelName+"_"+name+"_WidthG2", 
                             2.,0.1,15.)
      
    mixGG = root.RooRealVar(channelName+"_"+name+"_mixGG",
                            channelName+"_"+name+"_mixGG", 
                            0.9,0.0,1.0)
    gaus1 = root.RooGaussian(channelName+"_"+name+"_gaus1",
                             channelName+"_"+name+"_gaus1",
                             dimuonMass,meanG1,widthG1)
    gaus2 = root.RooGaussian(channelName+"_"+name+"_gaus2",
                             channelName+"_"+name+"_gaus2",
                             dimuonMass,meanG2,widthG2)
    sigPdfName = channelName+"__sig_Pdf"
    result = root.RooAddPdf(sigPdfName,
                                name,
                                gaus1,gaus2,mixGG)
    fitRange = root.RooFit.Range(120.,130.)
    if "cat1" in channelName:
      fitRange = root.RooFit.Range(117.,131.)
    if "cat0" in channelName:
      fitRange = root.RooFit.Range(118.,130.)
      #fitRange = root.RooFit.Range(122.,128.)
    if "cat2" in channelName:
      fitRange = root.RooFit.Range(113.,132.)
      #fitRange = root.RooFit.Range(115.,130.)
    fr = result.fitTo(sigData,root.RooFit.Save(),root.RooFit.SumW2Error(True),PRINTLEVEL,fitRange)
    fr.Print()

    frame = dimuonMass.frame(root.RooFit.Range(110.,160.))
    sigData.plotOn(frame)
    result.plotOn(frame,root.RooFit.Range(110.,160.))
    canvas = root.TCanvas()
    frame.Draw()
    canvas.SaveAs("sigHistTestFit_"+channelName+".png")

    getattr(self.wSig,"import")(sigData)
    getattr(self.wSig,"import")(result)
    getattr(self.wSig,"import")(gaus1)
    getattr(self.wSig,"import")(gaus2)
    return sigPdfName

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
    #if not "cat2" in fn:
    #    continue
    energyStr = "8TeV"
    print fn
    s = ShapePlotter(fn,outDir,titleMap,signalInject=args.signalInject,binWidthOverride=args.binWidthOverride,energyStr=energyStr)
    shapePlotterList.append(s)

  print "{0:20} & {1:>6} & {2:>6} & & & & {3:>6} & {4:>6}".format("cat","fwhm","nSig*10^5","nBkg","nData")
  for s in shapePlotterList:
    print r"{0:20} & {1:6.2f} & {2:6.2f} & & & & {3:6.1f} & {4:6.1f} &  \\".format(s.channelName,s.fwhm,s.nSignal/10.,s.nBkg,s.nData)
