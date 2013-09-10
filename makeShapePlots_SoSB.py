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


class ModelPlotter:
  def __init__(self,xVar,pdf,data,fr,title,energyStr,lumi,backgroundPDFName=None,signalPDFName=None,nSignal=0,signalPdf=None,signalLegEntry=None,RangeName=""):
    self.xVar = xVar
    self.pdf = pdf
    self.data = data
    self.fr = fr
    self.title = title
    self.energyStr = energyStr
    self.lumi = lumi
    self.backgroundPDFName = backgroundPDFName
    self.signalPDFName = signalPDFName
    self.nSignal = nSignal
    nowStr = str(int(time.time()*1e6))
    self.nowStr = nowStr

    self.lumiStr = "L = {0:.1f} fb^{{-1}}".format(lumi)

    xtitle = xVar.GetTitle()

    binning = xVar.getBinning()
    self.binning = binning
    nBins = self.binning.numBins()
    binWidth = (self.binning.highBound()-self.binning.lowBound())/nBins

    canvas = root.TCanvas("canvas"+nowStr)
    self.canvas = canvas

    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(root.gStyle.GetLabelSize())
    self.tlatex.SetTextAlign(22)

    doLinearErrs = True
    errVisArg = root.RooFit.VisualizeError(fr,1,doLinearErrs)
    errColorArg = root.RooFit.FillColor(root.kCyan)
    lineColorArg = root.RooFit.LineColor(root.kBlue)
    lineWidthArg = root.RooFit.LineWidth(2)
    lineDrawOptArg = root.RooFit.DrawOption("L")
    graphDrawOptArg = root.RooFit.DrawOption("PEZ")
    sigLineColorArg = root.RooFit.LineColor(root.kRed)

    binningArg = root.RooFit.Binning("")
    xVar.setRange("RMPRange",self.binning.lowBound(),self.binning.highBound())
    rangeArg = root.RooFit.Range("RMPRange")

    tmpDataHistName = data.GetName()+nowStr
    tmpBakPDFName = pdf.GetName()+nowStr
    tmpBakPDFErrorName = pdf.GetName()+"Error"+nowStr
    tmpDataHistNameArg = root.RooFit.Name(tmpDataHistName)
    tmpBakPDFNameArg = root.RooFit.Name(tmpBakPDFName)
    tmpBakPDFErrorNameArg = root.RooFit.Name(tmpBakPDFErrorName)

    nData = data.sumEntries()
    normErr = nData**(-0.5)
    self.normErr = normErr

    self.rangename = None
    # Main Frame
    frame       = xVar.frame(root.RooFit.Title(""))

    if (RangeName):
      self.rangename = RangeName
      frame = xVar.frame(root.RooFit.Title(""),root.RooFit.Range(self.rangename))

    self.frame = frame
    data.plotOn(frame,graphDrawOptArg,binningArg)
      
    if backgroundPDFName != None:
      bakCompArg = root.RooFit.Components(backgroundPDFName)
      pdf.plotOn(frame,errVisArg,errColorArg,bakCompArg,rangeArg,tmpBakPDFErrorNameArg)
      pdf.plotOn(frame,lineDrawOptArg,lineColorArg,lineWidthArg,bakCompArg,rangeArg,tmpBakPDFNameArg)
    else:
      pdf.plotOn(frame,errVisArg,errColorArg,rangeArg,tmpBakPDFErrorNameArg)
      pdf.plotOn(frame,lineDrawOptArg,lineColorArg,lineWidthArg,rangeArg,tmpBakPDFNameArg)

    data.plotOn(frame,graphDrawOptArg,binningArg,tmpDataHistNameArg)


    frame.SetTitle("")
    frame.GetXaxis().SetLabelSize(0)
    frame.GetYaxis().SetLabelSize(0.050)
    frame.GetYaxis().SetTitleSize(0.055*1.2)
    frame.GetYaxis().SetTitleOffset(
        0.85*frame.GetYaxis().GetTitleOffset()
        )

    
    unitMatch =  re.search(r"GeV([\s]*/[\s]*c\^\{2\}|[\s]*/[\s]*c)?",xtitle)
    units = ""
    if unitMatch:
      units = " "+unitMatch.group(0)
    frame.SetYTitle("Events/"+str(binWidth)+units)

    # Pulls Frame
    pullsHist = self.makePullPlotHist(frame,tmpDataHistName,tmpBakPDFName)
    self.pullsHist = pullsHist
    pullsHist.SetLineColor(root.kBlue)
    pullsHist.SetLineWidth(2)
    pullsHist.SetFillColor(856)
    pullsHist.SetFillStyle(1001)
    setHistTitles(pullsHist,xtitle,"#frac{Data-Fit}{#sqrt{Fit}}")

    pullsHist.GetXaxis().SetTitle(xtitle)
    pullsHist.GetXaxis().CenterTitle(1)
    self.pullsHist.GetXaxis().SetTitleSize(0.1334)
    self.pullsHist.GetXaxis().SetLabelSize(0.1213)
    self.pullsHist.GetXaxis().SetTitleOffset(
      self.pullsHist.GetXaxis().GetTitleOffset()*0.85
        )

    self.pullsHist.GetYaxis().CenterTitle(1)
    self.pullsHist.GetYaxis().SetTitleSize(0.097*1.2)
    self.pullsHist.GetYaxis().SetLabelSize(0.097)
    self.pullsHist.GetYaxis().SetTitleOffset(0.70*0.9)


    # add the signal PDF if it is there
    if signalPDFName != None and signalPdf == None:
        pdfList = pdf.pdfList()
        for i in range(pdfList.getSize()):
          tmp = pdfList[i]
          if tmp.GetName() == signalPDFName:
            signalPdf = tmp
            break
          elif tmp.InheritsFrom("RooExtendPdf"):
            srvr = tmp.findServer(signalPDFName)
            if srvr:
                signalPdf = srvr
                break
        self.signalPdf = signalPdf
    if signalPdf != None:
        sigPdfToDraw, componentToDraw = self.getProperSigPdfAndComponent(
                            xVar,data,signalPdf,self.nSignal
                            )

        sigPdfToDraw.plotOn(frame,      lineDrawOptArg,sigLineColorArg,lineWidthArg,componentToDraw),
      
        
    # Legend
    self.phonyFitLegHist = root.TH1F("phonyFit"+nowStr,"",1,0,1)
    self.phonyFitLegHist.SetFillColor(root.kCyan)
    self.phonyFitLegHist.SetLineColor(root.kBlue)
    self.phonyFitLegHist.SetLineWidth(2)

    self.phonyDatLegHist = root.TH1F("phonyFitDat"+nowStr,"",1,0,1)
    self.phonyDatLegHist.SetMarkerColor(1)
    self.phonyDatLegHist.SetLineColor(1)
    self.phonyDatLegHist.SetLineWidth(2)

    self.phonySigLegHist = root.TH1F("phonySigDat"+nowStr,"",1,0,1)
    self.phonySigLegHist.SetLineColor(root.kRed)
    self.phonySigLegHist.SetLineWidth(2)
    
    legPos = [0.65,0.65,1.0-gStyle.GetPadRightMargin()-0.01,1.0-gStyle.GetPadTopMargin()-0.01]
    self.legPos = legPos
    self.leg = root.TLegend(*legPos)
    self.leg.SetFillColor(0)
    self.leg.SetLineColor(0)
    self.leg.AddEntry(self.phonyDatLegHist,"Data","lp")
    self.leg.AddEntry(self.phonyFitLegHist,"Background Model","lf")

    if signalPdf != None:
      if signalLegEntry != None:
        self.leg.AddEntry(self.phonySigLegHist,signalLegEntry,"l")
      else:
        self.leg.AddEntry(self.phonySigLegHist,"Signal","l")


  def draw(self,filenameNoExt):
    nowStr = self.nowStr
    self.canvas.SetLogy(0)
    self.canvas.cd()
    pad1 = root.TPad("pad1"+nowStr,"",0.02,0.30,0.98,0.98,0)
    pad2 = root.TPad("pad2"+nowStr,"",0.02,0.01,0.98,0.29,0)
    self.pad1 = pad1
    self.pad2 = pad2
  
    pad1.SetBottomMargin(0.005);
    pad2.SetTopMargin   (0.005);
    pad2.SetBottomMargin(0.33);
  
    pad1.Draw() # Projections pad
    pad2.Draw() # Residuals   pad

    pad1Width = pad1.XtoPixel(pad1.GetX2())
    pad1Height = pad1.YtoPixel(pad1.GetY1())
    pad2Height = pad2.YtoPixel(pad2.GetY1())
    #pad1ToPad2FontScalingFactor = float(pad1Width)/pad2Height
    pad1ToPad2FontScalingFactor = float(pad1Height)/pad2Height
    canvasToPad1FontScalingFactor = float(self.canvas.YtoPixel(self.canvas.GetY1()))/pad1.YtoPixel(pad1.GetY1())
    canvasToPad2FontScalingFactor = float(self.canvas.YtoPixel(self.canvas.GetY1()))/pad2.YtoPixel(pad2.GetY1())
  
    # Main Pad
    pad1.cd();
    self.frame.Draw()
    #self.addPDFNormError(pad1)
    self.leg.Draw()

    # Pulls Pad
    pad2.cd();
    self.pullsHist.Draw("")
  
    # Text
    self.pad1.cd()
    #self.tlatex.SetTextSize(root.gStyle.GetLabelSize())
    self.tlatex.SetTextSize(0.04*canvasToPad1FontScalingFactor)
    self.tlatex.SetTextAlign(12)
    self.tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    self.tlatex.SetTextAlign(32)
    self.tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.title)
    #self.tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,"#sqrt{{s}}={0}, L={1:.1f} fb^{{-1}}".format(self.energyStr,self.lumi))

    self.tlatex.SetTextAlign(32)
    if (self.lumi != 0):
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.820,self.lumiStr)
    if (self.lumi == 0):
      self.tlatex.SetTextSize(0.04)
      self.tlatex.DrawLatex(self.legPos[0]-0.03,0.850,"#sqrt{s}=7 TeV L =  5.0 fb^{-1} ")
      self.tlatex.DrawLatex(self.legPos[0]-0.03,0.770,"#sqrt{s}=8 TeV L = 19.7 fb^{-1}")
            
    energyStr = self.energyStr
    if re.search(r"[\d]TeV",energyStr):
      energyStr = energyStr.replace("TeV"," TeV")
    if (self.energyStr != ""):
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.875,"#sqrt{s} = "+self.energyStr)

    saveAs(self.canvas,filenameNoExt)


  def drawPulls(self,filenameNoExt):
    self.canvas.cd()
    xMin = -5.
    xMax = 5.
    nBins = 20
    self.pullDistHist = root.TH1F("pullDist"+self.nowStr,"",nBins,xMin,xMax)
    setHistTitles(self.pullDistHist,"(Data-Fit)/#sqrt{Fit}","Events/%s" % (getBinWidthStr(self.pullDistHist)))
    self.pullDistHist.Sumw2()
    self.pullDistHist.SetMarkerColor(1)
    self.pullDistHist.SetLineColor(1)

    for i in range(1,self.pullsHist.GetXaxis().GetNbins()+1):
      self.pullDistHist.Fill(self.pullsHist.GetBinContent(i))

    fitFunc = root.TF1("fitFunc"+self.nowStr,"gaus",xMin,xMax)
    fitFunc.SetLineColor(root.kBlue)
    fitResult = self.pullDistHist.Fit(fitFunc,"LEMSQ")
    chi2 = fitFunc.GetChisquare()
    ndf = fitFunc.GetNDF()
    #print("chi2: {0:.2g}/{1}".format(chi2,ndf))
    #nParams =  fitFunc.GetNumberFreeParameters()
    #for i in range(nParams):
    #    parName = fitFunc.GetParName(i)
    #    val = fitFunc.GetParameter(i)
    #    err = fitFunc.GetParError(i)
    #    print("name: {}, value: {}, error: {}".format(parName,val,err))

    mean = fitFunc.GetParameter(1)
    meanErr = fitFunc.GetParError(1)
    sigma = fitFunc.GetParameter(2)
    sigmaErr = fitFunc.GetParError(2)
    
    self.pullDistHist.Draw()

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(root.gStyle.GetLabelSize())
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.title)
    if (self.energyStr != ""):
      tlatex.DrawLatex(0.98-gStyle.GetPadRightMargin(),0.875,"#sqrt{s} = "+self.energyStr)
    if (self.lumi != 0):
      tlatex.DrawLatex(0.98-gStyle.GetPadRightMargin(),0.825,self.lumiStr)

    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.875,"#chi^{{2}}/NDF = {0:.2g}".format(float(chi2)/ndf))
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.825,"#mu = {0:.2f} #pm {1:.2f}".format(mean,meanErr))
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.775,"#sigma = {0:.2f} #pm {1:.2f}".format(sigma,sigmaErr))
    #saveAs(self.canvas,filenameNoExt)


    
  def getProperSigPdfAndComponent(self,var,data,pdf,n):
    nData = data.sumEntries()
    nDummy = nData-n
    nowStr = self.nowStr
    self.dummyUnifPdf = root.RooUniform("dummyUnifPdf"+nowStr,"",root.RooArgSet(var))
    self.nDumbUnif = root.RooRealVar("dummyUnifN"+nowStr,"dummyUnifN",nDummy)
    self.nDumbSig = root.RooRealVar("dummySigN"+nowStr,"dummySigN",n)
    self.dummyUnifPdfE = root.RooExtendPdf("dummyUnifPdfE"+nowStr,"",self.dummyUnifPdf,self.nDumbUnif)
    self.dummySigPdfE = root.RooExtendPdf("dummySigPdfE"+nowStr,"",pdf,self.nDumbSig)
    self.dummyPdf = root.RooAddPdf("dummyPdf"+nowStr,"dummyPdf",root.RooArgList(self.dummySigPdfE,self.dummyUnifPdfE))
    return self.dummyPdf, root.RooFit.Components(pdf.GetName())

  def makePullPlotHist(self,frame,histPlotName,pdfPlotName):
    """
    Makes pulls that are (data-fit)/sqrt(fit) where fit is the average value of the
    PDF within the data histogram bin.
    """

    #print "\n\n\nStarting makePullPlotHist\n==============================================\n"
    hist = frame.findObject(histPlotName)
    curve = frame.findObject(pdfPlotName)
    assert(hist)
    assert(curve)

    nBins = self.binning.numBins()
    lowB  = self.binning.lowBound() 
    highB = self.binning.highBound() 
   
    if (self.rangename):
      binWidth = (highB - lowB) / nBins
      lowB  = self.xVar.getBinning(self.rangename).lowBound()
      highB = self.xVar.getBinning(self.rangename).highBound()
      nBins = (highB - lowB)/binWidth


    pullsHist = root.TH1F("pulls_"+histPlotName+"_"+pdfPlotName,"",
                          int(nBins), lowB, highB
                          )
    #print "pullsHist", pullsHist.GetNbinsX(), pullsHist.GetXaxis().GetXmin(), pullsHist.GetXaxis().GetXmax()
      
    x = root.Double(0.)
    y = root.Double(0.)

    curve.GetPoint(0,x,y) # Get Curve Start X
    xCurveMin = float(x)
    curve.GetPoint(curve.GetN()-1,x,y) # Get Curve End X
    xCurveMax = float(x)

    iBin = 1
    for i in range(1,self.binning.numBins()+1):
      hist.GetPoint(i-1,x,y)
      pull = float(y)
      if (float(x) < lowB or float(x) > highB):
        continue
      #print("hist bin: %10i, x: %10.2f, y: %10.2f" % (iBin,float(x),float(y)))
      if x > xCurveMin and x < xCurveMax:
        curvePoint = curve.interpolate(x)
        pull -= curvePoint
        pull /= sqrt(curvePoint)
        #print(" curve interpolation: %10.2f" % (curvePoint))
      else:
        pull = 0.
        #print(" Warning: x outside of curve range: [ %10.2f %10.2f ]" % (xCurveMin,xCurveMax))
      #print(" pull: %10.2f" % (pull))
      pullsHist.SetBinContent(iBin,pull)
      iBin += 1
      
    return pullsHist




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

    self.SoSB      = {}
    self.events    = {}
    self.datasets  = {}
    self.signalPDF = {}
    self.sigYield  = {}
    self.bakPDF    = {}
    self.bakYield  = {}

    tmpMatch = re.search(r"([\w]*)_(.+)_([.0-9]+)\.root",filename)

    if tmpMatch:
      self.energyStr = tmpMatch.group(2)
      if   (self.energyStr == "8TeV"):
        self.lumi = float(19.7)
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
      #if "8TeV" not in channelNameOrig:
      #if "Jets01PassPtG10BB8TeV" not in channelNameOrig:
      #  continue
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
      xlow  = binning.lowBound()
      xhigh = binning.highBound()
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
      nBackground = 0.
      nSignal = 0.
      for key in self.processNameMap[channelNameOrig]:
        #print key, self.processNameMap[channelNameOrig][key]
        if key != "bak":
          nSignal += self.processNameMap[channelNameOrig][key]
        if key == "bak":
          nBackground += self.processNameMap[channelNameOrig][key]

      nSignal *= signalInject
      signalLegEntry = "SM Higgs#times{0:.0f}".format(signalInject)

      #Set the PDF pars value from the FitResults
      setPDFfromFR(fr,bakPDF,data_obs)

      #Plot Time
      rmp = RooModelPlotter(mMuMu,
                            bakPDF,
                            data_obs,
                            fr,
                            channelTitle,self.energyStr,self.lumi,
                            nSignal=nSignal,signalPdf=sigPDF,
                            signalLegEntry=signalLegEntry
                            )
      rmp.draw(saveName)

      #Pull Distribution Time
      saveNameSplit = os.path.split(saveName)
      saveNamePulls = saveNameSplit[0]+"/"+"pulls_"+saveNameSplit[1]
      #rmp.drawPulls(saveNamePulls)

      self.rmpList.append(rmp)      

      #
      # Range
      #
      #
      ##########################################################################
      
      theRange = root.RooFit.Range("theRange")
      mMuMu.setRange("theRange", 122, 128)

      mMuMuArgSet = root.RooArgSet(mMuMu)
      rooNormSet  = root.RooFit.NormSet(mMuMuArgSet)
      nDataTotal  = data_obs.sumEntries("dimuonMass > 0.0")

      # bak integral
      bakInt = bakPDF.createIntegral(mMuMuArgSet,theRange,rooNormSet)
      nBak   = bakInt.getVal()*nDataTotal

      # sig integral
      sigInt = sigPDF.createIntegral(mMuMuArgSet,theRange,rooNormSet)
      nSig   = sigInt.getVal()*(nSignal/signalInject) # this remove the signalInjection

      # for debugging purposes
      print "%s %f %f %f %f %f %f %f" % (channelNameOrig,
                                         nBak,
                                         nSig,
                                         nSig/nBak,
                                         nSig/(nSig+nBak),
                                         nSignal/signalInject,
                                         nSignal,
                                         nBackground
                                         )

      self.SoSB     [channelNameOrig] = nSig/(nSig+nBak)
      self.events   [channelNameOrig] = nDataTotal

      self.datasets [channelNameOrig] = data_obs

      sigPDF.SetName("sig_"+channelNameOrig)
      self.signalPDF[channelNameOrig] = sigPDF
      self.sigYield [channelNameOrig] = nSignal/signalInject

      bakPDF.SetName("bak_"+channelNameOrig)
      self.bakPDF[channelNameOrig] = bakPDF
      self.bakYield [channelNameOrig] = nBackground
      
      self.fitresult = fr
      

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
      
      self.SoSB [key] = anotherSP.SoSB [key]
      self.events   [key] = anotherSP.events   [key]

      self.datasets [key] = anotherSP.datasets [key]

      self.signalPDF[key] = anotherSP.signalPDF[key]
      self.sigYield [key] = anotherSP.sigYield [key]

      self.bakPDF   [key] = anotherSP.bakPDF   [key]
      self.bakYield [key] = anotherSP.bakYield [key]

    

  def drawSoSB(self,name="sosbtest",title="",vetoList=[],folder=""):

    # sum the datasets with weights and plot them
    # define the x-axis variable
    dimuonMass = root.RooRealVar("dimuonMass","M(#mu#mu) [GeV/c^{2}]",110,160)
    datasetSoSB = None
    signalPdfList = []
    signalPdfCoeffList = []
    bakPdfList = []
    bakPdfCoeffList = []

    # compute the total number of events in the datasets you want to combine
    # -> needed for the weighting procedure
    nevents = 0
    sum_sosb_dot_events = 0

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

      sum_sosb_dot_events += self.SoSB[dname] * nDataTotal


    # loop again on the datasets  
    for dname in self.datasets.keys():

      if ( dname in vetoList ):
        continue

      #print dname
      dataset = self.datasets[dname]

      # construct the weight (making sure the total num events stays the same)
      Weight = root.RooRealVar("Weight_"+dname,"SoSB Weight",1)
      Weight.setVal( self.SoSB[dname] * nevents / sum_sosb_dot_events )

      # the set containing the observable and the weight
      observablesForRooDatasetWeight = root.RooArgSet(dimuonMass,
                                                      Weight)

      # constructing the dataset with a weight
      dataset.addColumn( Weight )
      dataset_withWeights = root.RooDataSet(dname + "_withWeights",
                                            dname + "_withWeights",
                                            observablesForRooDatasetWeight,
                                            root.RooFit.Import(dataset),
                                            root.RooFit.WeightVar(Weight)
                                            )
      # merging the datasets
      if datasetSoSB == None:
        datasetSoSB  = dataset_withWeights.Clone("datasetSoSB")
      else:
        datasetSoSB.append(dataset_withWeights)

      sigWeight = root.RooRealVar("sigWeight_"+dname,"sig Weight",1)
      sigWeight.setVal( self.sigYield[dname] * Weight.getVal() )
      
      signalPdfList.append( self.signalPDF[dname] )
      signalPdfCoeffList.append( sigWeight )
      
      bakWeight = root.RooRealVar("bakWeight_"+dname,"bak Weight",1)
      bakWeight.setVal( self.bakYield[dname] * Weight.getVal() )

      bakPdfList.append( self.bakPDF[dname] )
      bakPdfCoeffList.append( bakWeight )

      #print "%s, Weight.getVal() = %f " % (dname, Weight.getVal())
      #print "sigWeight.getVal() = %f "  % sigWeight.getVal()
      #print "bakWeight.getVal() = %f "  % bakWeight.getVal()


      
    # signal pdf
    coefflistSig = None
    for coeff in signalPdfCoeffList:
      if coefflistSig == None:
        coefflistSig = root.RooArgList(coeff)
      else:
        coefflistSig.add(coeff)
       
    pdflistSig = None
    for pdf in signalPdfList:
      if pdflistSig == None:
        setPDFfromFR(self.fitresult,
                     pdf,
                     self.datasets[dname])
        pdflistSig = root.RooArgList(pdf)
      else:
        setPDFfromFR(self.fitresult,
                     pdf,
                     self.datasets[dname])
        pdflistSig.add(pdf)

    signalPdfSum = root.RooAddPdf("signalPdfSum","signalPdfSum",
                                  pdflistSig,coefflistSig)



    coefflistBak = None
    for coeff in bakPdfCoeffList:
      if coefflistBak == None:
        coefflistBak = root.RooArgList(coeff)
      else:
       coefflistBak.add(coeff)
       
    pdflistBak = None
    for pdf in bakPdfList:
      if pdflistBak == None:
        setPDFfromFR(self.fitresult,
                     pdf,
                     self.datasets[dname])
        pdflistBak = root.RooArgList(pdf)
      else:
        setPDFfromFR(self.fitresult,
                     pdf,
                     self.datasets[dname])
        pdflistBak.add(pdf)


    bakPdfSum = root.RooAddPdf("bakPdfSum","bakPdfSum",
                               pdflistBak,coefflistBak)

    
    mMin = 110.
    mMax = 160.
    dimuonMass.setRange("whole",mMin,mMax)
    # 1 GeV bin
    dimuonMass.setBins( int(mMax-mMin) )

    #Plot Time
    rmp = ModelPlotter(dimuonMass,
                       bakPdfSum,
                       datasetSoSB,
                       self.fitresult,
                       title,self.energyStr,self.lumi,
                       None,None,
                       self.signalInject, # sum pdf already norm to num events (you can check it if you put 1)
                       signalPdfSum,
                       "Signal m_{H}=125 GeV #times %d" % ( self.signalInject )
                       )

    savename = outDir + name
    #print savename, savename, savename
    rmp.draw(savename+"_"+self.energyStr+self.massStr+"_SoSB")
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

  plotRange = [110.,160]
  normRange = [110.,120.,130.,160]

  rebin=1

  shapePlotterList = []
  for fn in glob.glob(dataDir+"*.root"):
    if re.search("P[\d.]+TeV",fn):
        continue

    if ("125" not in fn):
        continue
    if ("CombSplitAll" not in fn):
      continue
    
    print fn
    s = ShapePlotter(fn,outDir,fitDir,titleMap,signalInject=args.signalInject,binWidthOverride=args.binWidthOverride)
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

  # merge
  print "Drawing the full combination"
  sp_merged.drawSoSB("mMuMu_CombSplitAll7P8",
                     "S/(S+B) Weighted",
                     [], outDir)

