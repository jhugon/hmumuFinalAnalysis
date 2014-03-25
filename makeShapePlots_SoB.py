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

class ShapePlotter:
  def __init__(self,filename,outDir,fitDir,titleMap,signalInject=20.,binWidthOverride=0):
    self.signalInject=signalInject
    self.rmpList = []
    self.titleMap = titleMap
    self.filename = filename
    self.textFileName = os.path.splitext(filename)[0]+".txt"
    self.fitFileName = fitDir+"/"+os.path.split(self.textFileName)[1]+".root"
    self.processNameMap, self.params, self.normErrMap = getattr(self,"readCard")(self.textFileName)

    self.lumi = -1
    self.lumiStr = ""
    self.energyStr = ""
    self.massStr = ""

    self.SoB = {}
    self.SoSpB = {}
    self.SosqrtSpB = {}
    self.events = {}
    self.datasets = {}
    self.signalPDF = {}
    self.sigYield = {}

    tmpMatch = re.search(r"([\w]*)_(.+)_([.0-9]+)\.root",filename)

    if tmpMatch:
      self.energyStr = tmpMatch.group(2)
      if   (self.energyStr == "8TeV"):
        self.lumi = float(19.8)
      elif (self.energyStr == "7TeV"):
        self.lumi = float(5.05)
        
      #self.lumi = float(tmpMatch.group(3))
      self.lumiStr = "L = {0:.1f} fb^{{-1}}".format(self.lumi)

      self.massStr = tmpMatch.group(3)
      
    self.data = {}
    self.f = root.TFile(filename)
    for channelKey in self.f.GetListOfKeys():
      if channelKey.GetClassName() != "RooWorkspace":
        continue
      channelNameOrig = channelKey.GetName()
      channelName = re.sub("[\d]+TeV","",channelNameOrig)
      channelTitle = channelName
      if titleMap.has_key(channelName):
        channelTitle = titleMap[channelName]
      channelWS = channelKey.ReadObj()
      mMuMu = channelWS.var("dimuonMass")
      bakPDF = channelWS.pdf("bak")
      sigPDF = channelWS.pdf("ggH")
      data_obs = channelWS.data("data_obs")
      rooDataTitle = data_obs.GetTitle()
      binWidth = 1
      if "Jet2" in channelName or "VBF" in channelName:
          binWidth *= 2.5
      elif "BO" in channelName:
          binWidth *= 1
      elif "BE" in channelName:
          binWidth *= 2.5
      elif "OO" in channelName:
          binWidth *= 2.5
      elif "OE" in channelName:
          binWidth *= 2.5
      elif "EE" in channelName:
          binWidth *= 2.5
      elif "FF" in channelName:
          binWidth *= 2.5
      elif "CC" in channelName:
          binWidth *= 2.5
      elif "BB" in channelName:
          binWidth *= 1
      if binWidthOverride > 0:
        binWidth = binWidthOverride

      binning = mMuMu.getBinning()
      xlow = binning.lowBound()
      xhigh = binning.highBound()
      #mMuMu.setRange("shapePlot",xlow,xhigh)
      mMuMu.setBins(int((xhigh-xlow)/binWidth))
      mMuMu.SetTitle("M(#mu#mu) [GeV/c^{2}]")

      saveName = outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName
      saveName = re.sub(r"([\d]+)\.[\d]+",r"\1",saveName)

      # Get Fit Result
      fr = None
      try:
        rf = root.TFile(self.fitFileName)
        if rf.IsZombie() or not rf.IsOpen():
          raise IOError("ROOT File not open or Zombie")
        self.fitFile = rf
        #fr = rf.Get("fit_s")
        fr = rf.Get("fit_b")
      except Exception as e:
        print("Warning, Couldn't find ML fit file: {0}\n{1}".format(self.fitFileName,e))
        # Backup Fit
        fr = bakPDF.fitTo(data_obs,root.RooFit.Save(True),root.RooFit.PrintLevel(-1))

      # Signal Stuff
      nSignal = 0.
      for key in self.processNameMap[channelNameOrig]:
        if key != "bak":
          nSignal += self.processNameMap[channelNameOrig][key]
      nSignal *= signalInject
      legEntrySignal = "SM Higgs#times{0:.0f}".format(signalInject)

      #Set the PDF pars value from the FitResults
      setPDFfromFR(fr,bakPDF,data_obs)
   
      #Plot Time
      rmp = RooModelPlotter(mMuMu,bakPDF,data_obs,fr,
                            channelTitle,self.energyStr,self.lumi,
                            nSignal=nSignal,signalPdf=sigPDF,
                            legEntrySignal=legEntrySignal
                            )
      rmp.draw(saveName)

      #Pull Distribution Time
      saveNameSplit = os.path.split(saveName)
      saveNamePulls = saveNameSplit[0]+"/"+"pulls_"+saveNameSplit[1]
      rmp.drawPulls(saveNamePulls)

      self.rmpList.append(rmp)      

      #if ("Jet2" in channelNameOrig):
      #  continue
      
      #if ("Jets01Fail" in channelNameOrig):
      #  continue

      #
      # Range
      #
      #
      ##########################################################################
      # Indentify the narrow Gaussian sigma 
      # This part is somewhat hardcoded, in the sense that one needs to know
      # the names of the pdf and which is the narrow Gaussian
      # for the time being (05/08/2013) the name is always _ggHmumu125, with the
      # 125 even though the mass is different
      baseName = channelNameOrig+"_ggHmumu125_"+self.energyStr
      narrowGausName = baseName + "_gaus2"

      # get the narrow Gaussian
      narrowGaus     = sigPDF.pdfList().find(narrowGausName)

      #print narrowGausName 
      #narrowGaus.Print("v")

      paramSet    = narrowGaus.getParameters(data_obs)
      narrowSigma = paramSet.find(baseName + "_WidthG2")
      #print "%s = %s" % (baseName + "_WidthG2", narrowSigma.getVal())
      ##########################################################################

      
      SoBRange = root.RooFit.Range("SoBRange")

      # The range for the S/B calculation is dynamic and this works
      # only for 125 GeV/c2. *** to be fixed ***
      #lowSoBRange  = float(self.massStr)-3*narrowSigma.getVal()
      #highSoBRange = float(self.massStr)+2*narrowSigma.getVal()
      lowSoBRange  = 122.#float(self.massStr)-3*1.6
      highSoBRange = 128.#float(self.massStr)+2*1.6
      mMuMu.setRange("SoBRange", lowSoBRange, highSoBRange)

      mMuMuArgSet = root.RooArgSet(mMuMu)
      rooNormSet  = root.RooFit.NormSet(mMuMuArgSet)
      nDataTotal  = data_obs.sumEntries("dimuonMass > 0.0")

      # bak integral
      bakInt = bakPDF.createIntegral(mMuMuArgSet,SoBRange,rooNormSet)
      nBak   = bakInt.getVal()*nDataTotal

      # sig integral
      sigInt = sigPDF.createIntegral(mMuMuArgSet,SoBRange,rooNormSet)
      nSig   = sigInt.getVal()*(nSignal/signalInject) # this remove the signalInjection

      # for debugging purposes
      print "%s %f %f %f" % (channelNameOrig,
                             nBak,
                             nSig,
                             nSig/nBak)

      self.SoB      [channelNameOrig] = nSig/nBak 
      self.SoSpB    [channelNameOrig] = nSig/(nSig+nBak)
      self.SosqrtSpB[channelNameOrig] = nSig/( (nSig+nBak)**0.5 )
      self.events   [channelNameOrig] = nDataTotal
      self.datasets [channelNameOrig] = data_obs
      self.signalPDF[channelNameOrig] = sigPDF
      self.sigYield [channelNameOrig] = nSig
      

  def readCard(self,fn):
    f = open(fn)
    foundBin = False
    binList = []
    processList = []
    rateList = []
    paramMap = {}
    normErrMap = {}
    for line in f:
      if re.search("^bin",line):
        if foundBin:
          m =  re.findall("[\s]+[\w]+",line)
          binList.extend([i for i in m])
        else:
          foundBin = True
      if re.search("^process[\s]+[a-zA-Z]+",line):
          m =  re.findall("[\s]+[\w]+",line)
          processList.extend([i for i in m])
      if re.search("^rate[\s]+[-+eE.0-9]+",line):
          m =  re.findall("[\s]+[-+eE.0-9]+",line)
          rateList.extend([float(i) for i in m])
      paramMatch = re.search(r"([a-zA-Z0-9_]+)[\s]+param[\s]+([-.+eE0-9]+)[\s]+([-.+eE0-9]+)",line)
      if paramMatch:
        gs = paramMatch.groups()
        paramMap[gs[0]] = [gs[1],gs[2]]
      normMatch = re.search(r"bkN([a-zA-Z0-9_]+)[\s]+gmN[\s]+([-.+eE0-9]+)[\s-]+([.0-9]+)[\s-]+",line)
      if normMatch:
        gs = normMatch.groups()
        normErrMap[gs[0]] = 1.0/sqrt(float(gs[1]))
    binList = [x.replace(r" ","") for x in binList]
    processList = [x.replace(r" ","") for x in processList]
    result = {}
    for i in binList:
      if not result.has_key(i):
        result[i] = {}
    for b,p,r in zip(binList,processList,rateList):
      result[b][p] = r
    return result, paramMap, normErrMap

  # merge two ShapePlotter in one
  def merge(self, anotherSP):

    # quick sanity check
    #print "self.massStr = ", self.massStr
    #print "anotherSP.massStr = ", anotherSP.massStr
    if (self.massStr != anotherSP.massStr):
      print "ERROR merge: different mass values in the ShapePlotters\n"
      sys.exit(0)

    # merging...  
    self.lumi = 0
    self.energyStr = ""

    for key in anotherSP.datasets.keys():

      # here the assumption is that the keys of "anotherSP"
      # are different from the ones in "self"

      # avoid overwriting
      if (key in self.datasets.keys()):
        print "key %s already present -> skipping and not overwriting...\n"
        continue
      
      self.SoB      [key] = anotherSP.SoB      [key]
      self.SoSpB    [key] = anotherSP.SoSpB    [key]
      self.SosqrtSpB[key] = anotherSP.SosqrtSpB[key]
      self.events   [key] = anotherSP.events   [key]
      self.datasets [key] = anotherSP.datasets [key]
      self.signalPDF[key] = anotherSP.signalPDF[key]
      self.sigYield [key] = anotherSP.sigYield [key]

    

  def drawSoB(self,name="sobtest",title="",vetoList=[],folder=""):

    # sum the datasets with weights and plot them
    # define the x-axis variable
    dimuonMass = root.RooRealVar("dimuonMass","M(#mu#mu) [GeV/c^{2}]",110,170)
    datasetSoB = None
    signalPdfList = []
    signalPdfCoeffList = []

    # compute the total number of events in the datasets you want to combine
    # -> needed for the weighting procedure
    nevents = 0
    sumsobdotevents = 0
    sumsospbdotevents = 0

    nsignal = 0
    sumsobdotsigevents = 0

    for dname in self.datasets.keys():

      # if you need to veto some of the dataset:
      #  - only Jets01
      #  - only Jet2
      #  - ...
      if ( dname in vetoList ):
        continue

      dataset = self.datasets[dname]
      nDataTotal = dataset.sumEntries("dimuonMass > 0.0")
      nevents += nDataTotal
      sumsobdotevents += self.SoB[dname] * nDataTotal
      sumsospbdotevents += self.SoSpB[dname] * nDataTotal

      #print "self.sigYield[%s] = %f" % (dname, self.sigYield[dname])
      nsignal += self.sigYield[dname]
      sumsobdotsigevents += self.SoB[dname] * self.sigYield[dname]

    #print "nsignal = ", nsignal
    #id = 0
    # loop again on the datasets  
    for dname in self.datasets.keys():

      if ( dname in vetoList ):
        continue

      #print dname
      dataset = self.datasets[dname]

      # construct the weight (making sure the total num events stays the same)
      SoBWeight = root.RooRealVar("SoBWeight_"+name,"SoB Weight",1)
      SoBWeight.setVal( self.SoB[dname] * nevents / sumsobdotevents )

      # the set containing the observable and the weight
      observablesForRooDatasetWeight = root.RooArgSet(dimuonMass,
                                                      SoBWeight)

      # constructing the dataset with a weight
      dataset.addColumn( SoBWeight )
      dataset_withWeights = root.RooDataSet(dname + "_withWeights",
                                            dname + "_withWeights",
                                            observablesForRooDatasetWeight,
                                            root.RooFit.Import(dataset),
                                            root.RooFit.WeightVar(SoBWeight)
                                            )
      # merging the datasets
      if datasetSoB == None:
        datasetSoB  = dataset_withWeights.Clone("datasetSoB")
      else:
        datasetSoB.append(dataset_withWeights)

      sigWeight = root.RooRealVar("sigWeight_"+dname,"sig Weight",1)
      sigWeight.setVal( self.SoB[dname] * nsignal / sumsobdotsigevents )
      
      signalPdfList.append( self.signalPDF[dname] )
      signalPdfCoeffList.append( sigWeight )

      
    # signal pdf
    coefflist = None
    for coeff in signalPdfCoeffList:
      if coefflist == None:
        coefflist = root.RooArgList(coeff)
      else:
       coefflist.add(coeff)
       
    pdflist = None
    for pdf in signalPdfList:
      if pdflist == None:
        pdflist = root.RooArgList(pdf)
      else:
       pdflist.add(pdf)

    signalPdfSum = root.RooAddPdf("signalPdfSum","signalPdfSum",
                                  pdflist,coefflist)


    # fit with new values from Anna 13 Jun 2013
    InvPolMass = root.RooRealVar("InvPolMass","InvPolMass", 91.187, 30., 105.)
    ExpMass    = root.RooRealVar("ExpMass","ExpMass", 0.0, -2., 2.)

    bakPDF = root.RooGenericPdf("bak",
                                "TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",
                                root.RooArgList(dimuonMass,InvPolMass,ExpMass))

    mMin = 110.
    mMax = 170.
    dimuonMass.setRange("whole",mMin,mMax)
    # 1 GeV bin
    dimuonMass.setBins( int(mMax-mMin) )
    if (name == "Jet2SplitCutsGFSplit"):
      # 2.5 GeV bin
      dimuonMass.setBins( int ((mMax-mMin)/2.5) )
    
    fr = bakPDF.fitTo(datasetSoB,
                      root.RooFit.Range("whole"),
                      root.RooFit.SumW2Error(False),
                      root.RooFit.PrintLevel(-1),
                      root.RooFit.Save(True))
    # end fitting part


    # uncomment for debugging purpouses
    # print "totalweight = %f " % print totalweight
    # print nevents, " vs ", datasetSoB.sumEntries("dimuonMass > 0.0")
    # print "nsignal = %f" % nsignal
    # print "sumsobdotsigevents = %f" % sumsobdotsigevents

    #Plot Time
    rmp = RooModelPlotter(dimuonMass,bakPDF,datasetSoB,fr,
                          title,self.energyStr,self.lumi,
                          None,None,
                          self.signalInject, # sum pdf already norm to num events (you can check it if you put 1)
                          signalPdfSum,
                          "Signal m_{H}=126 GeV #times %d" % ( self.signalInject )
                          )

    savename = outDir + name
    rmp.draw(savename+"_"+self.energyStr+"_"+self.massStr+"_SoB")
    self.rmpList.append(rmp)

    zoom = "zoom"
    dimuonMass.setRange(zoom,120,135)
    rmp = RooModelPlotter(dimuonMass,bakPDF,datasetSoB,fr,
                          title,self.energyStr,self.lumi,
                          None,None,
                          7., # hardcoded
                          signalPdfSum,
                          "Signal m_{H}=126 GeV #times 7", #hardcoded
                          zoom
                          )

    rmp.draw(savename+"_"+self.energyStr+"_"+self.massStr+"_zoomedat126_SoB")
    self.rmpList.append(rmp)

    rmp.drawBkgSub(savename+"_"+self.energyStr+"_"+self.massStr+"_bkgSub_SoB")
    self.rmpList.append(rmp)


#~48 Charactars Max
titleMap = {
  "AllCat":"All Categories Comb.",
  "IncCat":"Non-VBF Categories Comb.",
  "VBFCat":"VBF Categories Comb.",

  "IncPresel":"Non-VBF Preselection",
  "VBFPresel":"VBF Preselection",

  "Pt0to30":"p_{T}^{#mu#mu} #in [0,30]",
  "Pt30to50":"p_{T}^{#mu#mu} #in [30,50]",
  "Pt50to125":"p_{T}^{#mu#mu} #in [50,125]",
  "Pt125to250":"p_{T}^{#mu#mu} #in [125,250]",
  "Pt250":"p_{T}^{#mu#mu}>250",

  "VBFLoose":"VBFL",
  "VBFMedium":"VBFM",
  "VBFTight":"VBFT",
  "VBFVeryTight":"VBFVT",

  "BDTCut":"BDT Cut Combination",
  "IncBDTCut":"Non-VBF BDT Cut",
  "VBFBDTCut":"VBF BDT Cut",

  "BDTCutCat":"BDT Cut Cat. Combination",
  "IncBDTCutCat":"Non-VBF BDT Cut",
  "VBFBDTCutCat":"VBF BDT Cut",

  "IncPreselCat":"Non-VBF Cat. Preselection",
  "VBFPreselCat":"VBF Cat. Preselection",

  "IncBDTCutBB":"Non-VBF BDT Cut BB",
  "IncBDTCutBO":"Non-VBF BDT Cut BO",
  "IncBDTCutBE":"Non-VBF BDT Cut BE",
  "IncBDTCutOO":"Non-VBF BDT Cut OO",
  "IncBDTCutOE":"Non-VBF BDT Cut OE",
  "IncBDTCutEE":"Non-VBF BDT Cut EE",
  "IncBDTCutNotBB":"Non-VBF BDT Cut !BB",
  "VBFBDTCutBB":"VBF BDT Cut BB",
  "VBFBDTCutNotBB":"VBF BDT Cut !BB",
  "IncPreselBB":"Non-VBF Preselection BB",
  "IncPreselBO":"Non-VBF Preselection BO",
  "IncPreselBE":"Non-VBF Preselection BE",
  "IncPreselOO":"Non-VBF Preselection OO",
  "IncPreselOE":"Non-VBF Preselection OE",
  "IncPreselEE":"Non-VBF Preselection EE",
  "IncPreselNotBB":"Non-VBF Preselection !BB",
  "VBFPreselBB":"VBF Preselection BB",
  "VBFPreselNotBB":"VBF Preselection !BB",

  "IncPreselPtG10BB":"Non-VBF BB",
  "IncPreselPtG10BO":"Non-VBF BO",
  "IncPreselPtG10BE":"Non-VBF BE",
  "IncPreselPtG10OO":"Non-VBF OO",
  "IncPreselPtG10OE":"Non-VBF OE",
  "IncPreselPtG10EE":"Non-VBF EE",
  "IncPreselPtG10NotBB":"Non-VBF !BB",

  "IncPreselPtG":"Non-VBF Not Combined",

  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
  "Jets01PassCatAll" : "0,1-Jet Tight ",
                        
  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
  "Jets01FailPtG10EE": "0,1-Jet Loose EE",
  "Jets01FailCatAll" : "0,1-Jet Loose ",
                        
  "Jets01SplitCatAll": "0,1-Jet ",


  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet VBF Loose",

  "Jet2SplitCutsGFSplit" : "2-Jet ",
  "CombSplitAll" : "Combination",
}


        
if __name__ == "__main__":
  #root.gROOT.SetBatch(True)

  dataDir = "statsCards/"
  outDir  = "shapes/"
  fitDir  = "statsInput/"

  plotRange = [110.,170]
  normRange = [110.,120.,130.,170]

  rebin=1

  shapePlotterList = []
  for fn in glob.glob(dataDir+"*.root"):
    if re.search("P[\d.]+TeV",fn):
        continue

    if ("126" not in fn):
        continue
    #if ("8TeV" not in fn):
    #    continue

    if ("CombSplitAll" not in fn):
      continue
    
    print fn
    s = ShapePlotter(fn,outDir,fitDir,titleMap,signalInject=args.signalInject,binWidthOverride=args.binWidthOverride)
    #print s.SoB
    #print s.sumSoB
    #print s.datasets
    shapePlotterList.append(s)


  catVetos = {}
  catVetos["Jets01SplitCatAll"] = [
                                   "Jet2CutsVBFPass",
                                   "Jet2CutsGFPass",
                                   "Jet2CutsFailVBFGF"
                                  ]

  catVetos["Jets01PassCatAll"] = [
                                  "Jets01FailPtG10BB",
                                  "Jets01FailPtG10BO",
                                  "Jets01FailPtG10BE",
                                  "Jets01FailPtG10OO",
                                  "Jets01FailPtG10OE",
                                  "Jets01FailPtG10EE",
                                    
                                  "Jet2CutsVBFPass",
                                  "Jet2CutsGFPass",
                                  "Jet2CutsFailVBFGF"
                                 ]
    

  catVetos["Jets01FailCatAll"] = [
                                  "Jets01PassPtG10BB",
                                  "Jets01PassPtG10BO",
                                  "Jets01PassPtG10BE",
                                  "Jets01PassPtG10OO",
                                  "Jets01PassPtG10OE",
                                  "Jets01PassPtG10EE",
                                  
                                  "Jet2CutsVBFPass",
                                  "Jet2CutsGFPass",
                                  "Jet2CutsFailVBFGF"
                                  ]

  catVetos["Jet2SplitCutsGFSplit"] = [
                                      "Jets01PassPtG10BB",
                                      "Jets01PassPtG10BO",
                                      "Jets01PassPtG10BE",
                                      "Jets01PassPtG10OO",
                                      "Jets01PassPtG10OE",
                                      "Jets01PassPtG10EE",
                                      
                                      "Jets01FailPtG10BB",
                                      "Jets01FailPtG10BO",
                                      "Jets01FailPtG10BE",
                                      "Jets01FailPtG10OO",
                                      "Jets01FailPtG10OE",
                                      "Jets01FailPtG10EE"
                                    ]

  
  catVetos["CombSplitAll"] = []
  
  
  sp_merged = None
  for s in shapePlotterList:

    if (sp_merged == None):
      sp_merged = s
    else:
      sp_merged.merge(s)

#    for cat in catVetos:
#
#      title    = titleMap[cat]
#      vetoList = [ x+s.energyStr for x in catVetos[cat] ]
#      print cat, title, vetoList
#      
#      s.drawSoB(cat, title+" (S/B weighted)", vetoList, outDir)

  # merge
  print "Drawing the full combination"
  sp_merged.drawSoB("CombSplitAll7P8",
                    "Combination (S/B weighted)",
                    [], outDir)

