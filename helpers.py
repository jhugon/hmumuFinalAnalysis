
import ROOT as root
from ROOT import gStyle as gStyle
#from ROOT import RooRealVar, RooGaussian, RooArgList, RooDataHist
import re
import csv
import glob
from math import exp
from math import sqrt
from math import log10
import math
import numpy
import scipy
import scipy.stats
import array
import os
import sys
import time
import datetime
#import matplotlib.pyplot as mpl

#PRELIMINARYSTRING="CMS Internal"
PRELIMINARYSTRING="CMS Preliminary"
#PRELIMINARYSTRING="CMS"

def getDataStage2Directory():
  hostname = os.uname()[1]
  if "melrose" in hostname:
    result = "/raid/raid8/jhugon/higgsSamples/stage2/"
  elif "uftrig" in hostname:
    result = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
  elif "cern" in hostname: # Works on lxplus and lxbatch
    result = "/afs/cern.ch/work/j/jhugon/public/hmumuNtuplesLevel2/unzipped/"
  elif "cyril" == hostname:
    result = "/opt/hep/data/hmumu/stage2/"
  else: # Works on ihepa computers and UF HPC
    result = "/cms/data/store/user/jhugon/hmumu/stage2/"
  #print "Using Input Directory: ",result
  return result

def sortCatNames(l):
  orderDef = [
    "CombSplitAll",
    "Jets01SplitCatAll",
    "Jet2SplitCutsGFSplit",
    "Jets01PassCatAll" ,
    "Jets01FailCatAll" ,

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
    "Jets01FailPtG10EE",

    "Jet2CutsVBFPass",
    "Jet2CutsGFPass",
    "Jet2CutsFailVBFGF",
  ]
  return sorted(l,key=lambda x: orderDef.index(x))

def getOrdinalStr(inInt):
  result = str(inInt)
  if result[-1] == "1":
    result += "st"
  elif result[-1] == "2":
    result += "nd"
  elif result[-1] == "3":
    result += "rd"
  else:
    result += "th"
  return result


def drange(start, stop, step):
  r = start
  while r < stop:
    yield r
    r += step

def revdrange(start, stop, step):
  r = start
  while r > stop:
    yield r
    r -= step

# calculate FWHM
def calcFWHM(pdf,obs,min,max,step):

  var = pdf.getObservables(root.RooArgSet(obs)).first();

  ymaxVal = float(0)
  xmaxVal = float(0)

  # find the maximum value
  for x in drange(min,max,step):
   
    var.setVal(x)
    pdfVal = pdf.getVal(root.RooArgSet(var)) 

    if (pdfVal > ymaxVal):
       xmaxVal = x
       ymaxVal = pdfVal
       
    #print "x=%s, pdfVal=%s" % (x,pdfVal)

  #print "xMax=%s, ymaxVal=%s\n\n\n\n" % (xmaxVal,ymaxVal)


  # find lower boundary with y=max/2
  xLow = float(0)
  for x in drange(min,max,step):
   
    var.setVal(x)
    pdfVal = pdf.getVal(root.RooArgSet(var)) 

    #print "x=%s, pdfVal=%s, ymaxVal/2.=%s" % (x,pdfVal, ymaxVal/2.)
    if (pdfVal > ymaxVal/2. and xLow==0):
       xLow = x

  #print "xLow=%s" % xLow
  

  # find higher boundary with y=max/2
  xHigh = float(0)
  for x in revdrange(max,min,step):
   
    var.setVal(x)
    pdfVal = pdf.getVal(root.RooArgSet(var)) 

    if (pdfVal > ymaxVal/2. and xHigh==0):
       xHigh = x

  #print "xHigh=%s" % xHigh
  
  return (xHigh-xLow)


  
#fr stands for FitResults
def setAddPDFfromFR(fr,PDF,data):

  # print it for debugging
  #fr.Print()

  fitpars = {}
  fiterrs = {}

  for i in range(0,fr.floatParsFinal().getSize()):
    
    parName = fr.floatParsFinal().at(i).GetName()
    parValue= fr.floatParsFinal().at(i).getVal()
    parErr  = fr.floatParsFinal().at(i).getError()
    
    # print it for debugging
    #print i, parName, parValue, parErr 

    fitpars[parName] = float(parValue)
    fiterrs[parName] = float(parErr)

    # print it for debugging
    #print  fitpars

    pdfList  =  PDF.pdfList()
    coefList =  PDF.coefList()
    
    # print it for debugging
    #print "pdfSize  = ", pdfList.getSize() 
    #print "coefSize = ", coefList.getSize()
    
    for c in range(0,coefList.getSize()):
        coefName = coefList[0].GetName()
        print coefName
        if (coefName in fitpars.keys()):
            coefList[0].setVal  ( fitpars[coefName] )
            coefList[0].setError( fiterrs[coefName] )
        else:
            print "Potential Problem: No Coefficient Matching Found"
    
    
    for p in range(0,pdfList.getSize()):
        pdfName = pdfList[p].GetName()
        print pdfName
    
        # get the parameters
        parameters = pdfList[p].getParameters(data)
    
        # get the iterator
        iterator =  parameters.createIterator()
    
        cond = True
    
        while (cond):
            par = iterator.Next()
    
            if (par == None):
                cond = False
    
            else:
                parName = par.GetName()
                print 'the parameter is ', parName
                if (parName in fitpars.keys()):
                    par.setVal  ( fitpars[parName] )
                    par.setError( fiterrs[parName] )
                else:
                    print "Potential Problem: No Coefficient Matching Found"


def setPDFfromFR(fr,PDF,data):

  # print it for debugging
  #fr.Print()

  fitpars = {}
  fiterrs = {}

  for i in range(0,fr.floatParsFinal().getSize()):
    
    parName = fr.floatParsFinal().at(i).GetName()
    parValue= fr.floatParsFinal().at(i).getVal()
    parErr  = fr.floatParsFinal().at(i).getError()
    
    # print it for debugging
    #print i, parName, parValue, parErr 

    fitpars[parName] = float(parValue)
    fiterrs[parName] = float(parErr)

    # print it for debugging
    #print  fitpars

    # get the parameters
    parameters = PDF.getParameters(data)
    
    # get the iterator
    iterator =  parameters.createIterator()
    
    cond = True
    
    while (cond):
        par = iterator.Next()
    
        if (par == None):
            cond = False
    
        else:
            parName = par.GetName()
            #print 'the parameter is ', parName
            if (parName in fitpars.keys()):
                par.setVal  ( fitpars[parName] )
                par.setError( fiterrs[parName] )
            #else:
                #print "Potential Problem: No Coefficient Matching Found"

def rooObjPrintVars(obj):
  obj.Print()
  for i in rooArgSet2List(obj.getVariables()):
    i.Print()

def rooArgSet2List(x):
  itr = x.createIterator()
  result = []
  while True:
    ele = itr.Next()
    if ele:
      result.append(ele)
    else:
      break
  return result

def rooArgSet2Dict(x):
  itr = x.createIterator()
  result = {}
  while True:
    ele = itr.Next()
    if ele:
      result[ele.GetName()] = ele
    else:
      break
  return result

def rooPdfNFreeParams(pdf,data):
  paramList = rooArgSet2List(pdf.getParameters(data))
  n = 0
  for i in paramList:
    if not i.isConstant():
        n += 1
  return n

def doubleGauss(x,par):
  meanG1  = par[0]
  widthG1 = par[1]
  meanG2  = par[2]
  widthG2 = par[3]
  mixGG   = par[4]
  scale   = par[5]
  
  #if (par[1] != 0.0):
  
  arg1 = (x[0]-meanG1)/widthG1
  arg2 = (x[0]-meanG2)/widthG2
  
  gauss1 = exp(-0.5*arg1*arg1)
  gauss2 = exp(-0.5*arg2*arg2)
  dgauss = (1-mixGG)*gauss1 + mixGG*gauss2 
  
  return scale*dgauss
  #return meanG1 + widthG1*x[0]
  
def fit2DResHist(hist,color):
  histName = hist.GetName()
  hist.FitSlicesY()
  outMeanHist = root.gDirectory.Get(histName+"_1")
  outSigmaHist = root.gDirectory.Get(histName+"_2")

  outMeanHist.SetName(histName+"FitMean")
  outSigmaHist.SetName(histName+"FitSigma")
  outMeanHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
  outSigmaHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
  outMeanHist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
  outSigmaHist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())

  outSigOverMeanHist = outSigmaHist.Clone(histName+"FitSigmaOverMean")
  outSigOverMeanHist.Divide(outMeanHist)

  outMeanHist.SetMarkerColor(color)
  outMeanHist.SetLineColor(color)
  outSigmaHist.SetMarkerColor(color)
  outSigmaHist.SetLineColor(color)
  outSigOverMeanHist.SetMarkerColor(color)
  outSigOverMeanHist.SetLineColor(color)

  outMeanHist.SetTitle("")
  outSigmaHist.SetTitle("")
  outSigOverMeanHist.SetTitle("")

  return outMeanHist, outSigmaHist, outSigOverMeanHist

class Fit1DResHist:
  def __init__(self,hist1,hist2,canvas,fitList=[0.9,1.1],xlimits=[0.5,1.5],captions=[]):
    histName1 = hist1.GetName()

    self.captions = captions

    self.canvas = canvas
    canvas.cd()

    self.hist1 = hist1
    self.hist2 = hist2

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextSize(0.035)
    tlatex.SetTextAlign(12)
    self.tlatex = tlatex
  
    gaus1 = root.TF1("gaus1"+histName1,"gaus",fitList[0],fitList[1])
    gaus2 = root.TF1("gaus2"+histName1,"gaus",fitList[0],fitList[1])
    gaus1.SetLineColor(root.kBlack)
    gaus2.SetLineColor(root.kRed+1)
    self.gaus1 = gaus1
    self.gaus2 = gaus2

    hist1.Fit(gaus1,"WLMEIQ","",fitList[0],fitList[1])
    hist2.Fit(gaus2,"WLMEIQ","",fitList[0],fitList[1])

    hist1.GetXaxis().SetRangeUser(*xlimits)
    hist1.GetXaxis().SetTitle("p_{T}^{X}/p_{T}^{True}")
    
    hist1.SetLineColor(1)
    hist1.SetMarkerColor(1)
    hist2.SetLineColor(root.kRed+1)
    hist2.SetMarkerColor(root.kRed+1)
    hist1.Draw()
    hist2.Draw("same")
    
    tlatex.SetTextAlign(12)
    tlatex.SetTextSize(0.022)
    #tlatex.DrawLatex(0.65,0.55,"Mean_{RECO} = %.3f" % 
    #        (gaus1.GetParameter(1)))
    #tlatex.DrawLatex(0.65,0.5,"Mean_{Regr.} = %.3f" % 
    #        (gaus2.GetParameter(1)))
    #tlatex.DrawLatex(0.65,0.45,"#sigma_{RECO} = %.3f" % 
    #        (gaus1.GetParameter(2)))
    #tlatex.DrawLatex(0.65,0.4,"#sigma_{Regr.} = %.3f" % 
    #        (gaus2.GetParameter(2)))
    tlatex.DrawLatex(0.65,0.55,"Mean_{RECO} = %.3f #pm %.3f" % 
            (gaus1.GetParameter(1),gaus1.GetParError(1)))
    tlatex.DrawLatex(0.65,0.5,"Mean_{Regr.} = %.3f #pm %.3f" % 
            (gaus2.GetParameter(1),gaus2.GetParError(1)))
    tlatex.DrawLatex(0.65,0.45,"#sigma_{RECO} = %.3f #pm %.3f" % 
            (gaus1.GetParameter(2),gaus1.GetParError(2)))
    tlatex.DrawLatex(0.65,0.4,"#sigma_{Regr.} = %.3f #pm %.3f" % 
            (gaus2.GetParameter(2),gaus2.GetParError(2)))
    
    tlatex.SetTextSize(0.035)
    i = 0.85
    mini = 0.3
    for el in captions:
        tlatex.DrawLatex(0.2,i,el)
        i = i - 0.05
        if(i<mini):
            break

def fit2DResHistManual(hist,color,fitDictList=[],xlimits=[0.0,2.0]):
  histName = hist.GetName()
  histSaveName = "slices_"+histName+".pdf"

  outMeanHist = hist.ProjectionX()
  outMeanHist.Clear()
  outMeanHist.Sumw2()
  outMeanHist.SetName(histName+"FitMean")
  outSigmaHist = outMeanHist.Clone(histName+"FitSigma")
  outRMSHist = outMeanHist.Clone(histName+"FitRMS")

  tmpCanvas = root.TCanvas("tmpCanvas"+histName)
  tmpCanvas.cd()
  tmpCanvas.SetGridx(1)
  tlatex = root.TLatex()
  tlatex.SetNDC()
  tlatex.SetTextSize(0.035)
  tlatex.SetTextAlign(12)

  nBins = hist.GetXaxis().GetNbins()
  for i in range(0,nBins+2):
    sliceHist = getXBinHist(hist,i)

    mean = sliceHist.GetMean()
    meanErr = sliceHist.GetMeanError()
    rms = sliceHist.GetRMS()
    rmsErr = sliceHist.GetRMSError()

    maxBin = sliceHist.GetMaximumBin()
    xMax = sliceHist.GetXaxis().GetBinCenter(maxBin)
    gaus = root.TF1("gaus","gaus");
    if fitDictList==[]:
      sliceHist.Fit(gaus,"WLMEIQ","",xMax-0.2,xMax+0.2)
    else:
      sliceHist.Fit(gaus,"WLMEIQ","",fitDictList[i][0],fitDictList[i][1])
    sliceHist.SetTitle("pt: "+str(hist.GetXaxis().GetBinCenter(i))+"   iBin: "+str(i))
    sliceHist.GetXaxis().SetRangeUser(*xlimits)
    sliceHist.Draw()

    # 1 is mean, 2 is sigma
    fitMean = gaus.GetParameter(1)
    fitMeanErr = gaus.GetParError(1)
    fitSig = gaus.GetParameter(2)
    fitSigErr = gaus.GetParError(2)
    fitChi2 = gaus.GetChisquare()
    fitNDF = gaus.GetNDF()
    if fitMean ==0.0:
        fitMean = 1.0e-10
    if fitSig ==0.0:
        fitSig = 1.0e-10
    if rms ==0.0:
        rms = 1.0e-10
    outMeanHist.SetBinContent(i,fitMean)
    outMeanHist.SetBinError(i,fitMeanErr)
    outSigmaHist.SetBinContent(i,fitSig/fitMean)
    outSigmaHist.SetBinError(i,fitSig/fitMean*sqrt((fitSigErr/fitSig)**2+(fitMeanErr/fitMean)**2))
    outRMSHist.SetBinContent(i,rms/fitMean)
    outRMSHist.SetBinError(i,rms/fitMean*sqrt((rmsErr/rms)**2+(fitMeanErr/fitMean)**2))

    """
    outMeanHist.SetBinContent(i,mean)
    outMeanHist.SetBinError(i,meanErr)
    outSigmaHist.SetBinContent(i,rms)
    outSigmaHist.SetBinError(i,rmsErr)
    """

    """
    ptDiff = RooRealVar("ptDiff","pt - pt",hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    mean1 = RooRealVar("mean1","<pt-pt>",0.0)
    sigma1 = RooRealVar("sigma1","#sigma(pt-pt)",5.0)
    mean1.setConstant(False)
    sigma1.setConstant(False)

    frame1 = ptDiff.frame()

    gaus1 = RooGaussian("guas1","Guassian",ptDiff,mean1,sigma1)
    data = RooDataHist("data","Data",RooArgList(ptDiff),sliceHist)
    gaus1.fitTo(data)

    data.plotOn(frame1)
    gaus1.plotOn(frame1)

    outMeanHist.SetBinContent(i,mean1.GetVal())
    outMeanHist.SetBinError(i,mean1.GetError())
    outSigmaHist.SetBinContent(i,sigma1.GetVal())
    outSigmaHist.SetBinError(i,sigma1.GetError())

    frame1.SetTitle("pt: "+str(hist.GetXaxis().GetBinCenter(i))+"   iBin: "+str(i))
    frame1.GetXaxis().SetRangeUser(-100,100)
    frame1.Draw()
    """

    tlatex.DrawLatex(0.6,0.75,"Mean: %.3f #pm %.3f" % (fitMean,fitMeanErr))
    tlatex.DrawLatex(0.6,0.8,"#sigma: %.3f #pm %.3f" % (fitSig,fitSigErr))
    tlatex.DrawLatex(0.6,0.7,"#chi/NDF: %.3f / %.3f" % (fitChi2,fitNDF))
    if(i == 0):
	  tmpCanvas.SaveAs(histSaveName+"(")
	  tmpCanvas.Clear()
    elif(i == nBins+1):
	  tmpCanvas.SaveAs(histSaveName+")")
	  tmpCanvas.Clear()
    else:
	  tmpCanvas.SaveAs(histSaveName)
	  tmpCanvas.Clear()

  outMeanHist.SetMarkerColor(color)
  outMeanHist.SetLineColor(color)
  outSigmaHist.SetMarkerColor(color)
  outSigmaHist.SetLineColor(color)
  outRMSHist.SetMarkerColor(color)
  outRMSHist.SetLineColor(color)
  outMeanHist.SetTitle("")
  outSigmaHist.SetTitle("")
  outRMSHist.SetTitle("")

  #outMeanHist.Print()
  #outSigmaHist.Print()
  #outRMSHist.Print()

  return outMeanHist, outSigmaHist, outRMSHist

def getXBinHist(inHist, xBin):
  outHist = inHist.ProjectionY()
  outHist.Clear()
  outHist.SetName(inHist.GetName()+"XSliceBin"+str(xBin))
  outHist.Sumw2()
  nBins = outHist.GetXaxis().GetNbins()
  for i in range(0,nBins+2):
    outHist.SetBinContent(i,inHist.GetBinContent(xBin,i))
    outHist.SetBinError(i,inHist.GetBinError(xBin,i))
  return outHist

def getYBinHist(inHist, yBin):
  outHist = inHist.ProjectionX()
  outHist.Clear()
  outHist.SetName(inHist.GetName()+"YSliceBin"+str(yBin))
  outHist.Sumw2()
  nBins = outHist.GetYaxis().GetNbins()
  for i in range(0,nBins+2):
    outHist.SetBinContent(i,inHist.GetBinContent(i,yBin))
    outHist.SetBinError(i,inHist.GetBinError(i,yBin))
  return outHist

def fitSlicesTopRooFit(hist):
  #mTopBest = 172.9
  mass = root.RooRealVar("mass","m_{jjj} of Top Candidate [GeV]",20,250)
  #mtParam = root.RooRealVar("mtParam","m_{t} Parameter",mTopBest);
  #wtParam = root.RooRealVar("wtParam","#Gamma_{t} Parameter",2.0);

  smearMeanParam = root.RooRealVar("smearMeanParam","m_{t} Smearing Mean",0.0);
  smearSigParam = root.RooRealVar("smearSigParam","m_{t} Smearing Sigma",0.5);

  smearMeanParam.setConstant(False);
  smearSigParam.setConstant(False);

  #bwPdf = root.RooBreitWigner("bwPdf","m_{t} Breit-Wigner Pdf",
  #						mass,mtParam,wtParam);
  smearPdf = root.RooGaussian("smearPdf","m_{t} Smearing PDF",
						mass,smearMeanParam,smearSigParam);
  #convPdf = root.RooFFTConvPdf("convPdf","m_{t} Convolution of BW with Gaussian",
  #						mass,bwPdf,smearPdf);

  outMeanGraph = root.TGraphErrors()
  outSigGraph = root.TGraphErrors()
  if(hist.GetDimension()==1):
      
      mass = root.RooRealVar("mass","m_{jjj} of Top Candidate [GeV]",hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax())
      smearMeanParam = root.RooRealVar("smearMeanParam","m_{t} Smearing Mean",0.0);
      smearSigParam = root.RooRealVar("smearSigParam","m_{t} Smearing Sigma",0.5);
      smearMeanParam.setConstant(False);
      smearSigParam.setConstant(False);
      smearPdf = root.RooGaussian("smearPdf","m_{t} Smearing PDF",
						mass,smearMeanParam,smearSigParam);

      data = root.RooDataHist("data",
  			"t Candidate m_{jjj} from Measured b-jet and Fit W-mass",
  			root.RooArgList(mass),hist)
  
      smearPdf.fitTo(data)
                                                                            
      smearMean = smearMeanParam.getVal()
      smearMeanErr = smearMeanParam.getError()
      smearSig = smearSigParam.getVal()
      smearSigErr = smearSigParam.getError()
      return smearMean, smearMeanErr, smearSig, smearSigErr

  elif(hist.GetDimension()==2):
    iPoint=0;
    nBinsX = hist.GetXaxis().GetNbins()
    for iBinX in range(1,nBinsX+1):
      dataHist = getXBinHist(hist,iBinX);
      dataHist.Rebin(2)
      data = root.RooDataHist("data",
  			"t Candidate m_{jjj} from Measured b-jet and Fit W-mass",
  			root.RooArgList(mass),dataHist)
  
      print("justin bin: "+str(iBinX))
      dataHist.Print("v")
      data.Print("v")
      fitResult = convPdf.fitTo(data,root.RooFit.Save())
                                                                            
      smearMean = smearMeanParam.getVal()
      smearMeanErr = smearMeanParam.getError()
      smearSig = smearSigParam.getVal()
      smearSigErr = smearSigParam.getError()
  
      xPos = hist.GetXaxis().GetBinCenter(iBinX)
      xErr = hist.GetXaxis().GetBinWidth(iBinX)/2.0
  
      #should it be divided by mTopBest or xPos?
      smearSig = smearSig/xPos
      smearSigErr = smearSigErr/xPos
      smearMean = smearMean/xPos
      smearMeanErr = smearMeanErr/xPos
      """
      outMeanGraph.SetPoint(iPoint,xPos,dataHist.GetMean())
      outMeanGraph.SetPointError(iPoint,xErr,dataHist.GetMeanError())
      smearSig = dataHist.GetRMS()/mTopBest
      smearSigErr = dataHist.GetRMSError()/mTopBest
      """
  
      outMeanGraph.SetPoint(iPoint,xPos,smearMean)
      outMeanGraph.SetPointError(iPoint,xErr,smearMeanErr)
      outSigGraph.SetPoint(iPoint,xPos,smearSig)
      outSigGraph.SetPointError(iPoint,xErr,smearSigErr)
  
      iPoint = iPoint+1

  return outMeanGraph, outSigGraph

def divideYValByXVal(hist):
    nBinsX = hist.GetXaxis().GetNbins()
    for iBinX in range(1,nBinsX+1):
	binVal = hist.GetBinContent(iBinX)
	binErrVal = hist.GetBinError(iBinX)
	xVal = hist.GetXaxis().GetBinCenter(iBinX)
	hist.SetBinContent(iBinX,binVal/xVal)
	hist.SetBinError(iBinX,binErrVal/xVal)

def drawMVAHist(tfile, histToGetRE, meanOrsigmaString):
  """
   Draws a hist from the TMVA out file.  Must already have an active TCanvas
  """
  leg = root.TLegend(0.70,0.6,0.9,0.9)
  leg.SetFillColor(root.kWhite)
  leg.SetLineColor(root.kWhite)
  colors = [root.kBlack,root.kBlue+1,root.kRed+1,root.kBlue+1,root.kOrange+1]
  iColor = 0
  for dirName in tfile.GetListOfKeys():
   if(re.match(r"Method_.*",dirName.GetName())):
     for subDirName in dirName.ReadObj().GetListOfKeys():
       for histName in subDirName.ReadObj().GetListOfKeys():
	if(re.match(histToGetRE,histName.GetName())):
	  hist = histName.ReadObj()
	  #hist.GetXaxis().Rebin(2)
	  hM, hS, hSOM = fit2DResHist(hist)
	  hM.SetLineColor(colors[iColor % len(colors)])
	  hM.SetMarkerColor(colors[iColor % len(colors)])
	  hS.SetLineColor(colors[iColor % len(colors)])
	  hS.SetMarkerColor(colors[iColor % len(colors)])
	  hSOM.SetLineColor(colors[iColor % len(colors)])
	  hSOM.SetMarkerColor(colors[iColor % len(colors)])
	  tmpLegLabel = re.sub(r"MVA_","",histName.GetName())
	  tmpLegLabel = re.sub(r"test.*","",tmpLegLabel)
	  leg.AddEntry(hM,tmpLegLabel,"lep")
	  divideYValByXVal(hS)
	  
	  if(iColor==0):
	    if(meanOrsigmaString == "mean"):
		hM.SetTitle("")
		hM.GetYaxis().SetRangeUser(-10,10)
		hM.Draw()
	    else:
		hS.SetTitle("")
		hS.GetYaxis().SetRangeUser(0.0,0.6)
		hS.Draw()
	  else:
	    if(meanOrsigmaString == "mean"):
		hM.Draw("same")
	    else:
		hS.Draw("same")
	  iColor += 1
  leg.Draw()
  leg.Print()
  print("made it to leg.Draw")

def makeResPlotFromMVATree(tfile,xname,yname,truename,cuts,doDivide):
  histName2D = "mvaResTmpHist"
  tree = tfile.Get("TestTree")
  #tree = tfile.Get("TrainTree")
  drawStr = ""
  if doDivide:
    drawStr = "("+yname+"-"+truename+")/"+truename+":"+xname
  else:
    drawStr = yname+"-"+truename+":"+xname
  drawStr += ">>"
  drawStr += histName2D
  print(drawStr)
  tree.Draw(drawStr,cuts)
  tmpHist = root.gDirectory.Get(histName2D)
  tmpM, tmpS, tmpSOM = fit2DResHist(tmpHist,root.kBlack)
  if(doDivide):
    tmpM.GetYaxis().SetRangeUser(-1,1)
    tmpS.GetYaxis().SetRangeUser(0,0.5)
  else:
    tmpM.GetYaxis().SetRangeUser(-50,50)
    tmpS.GetYaxis().SetRangeUser(0,100)
  tmpM.SetTitle("")
  tmpS.SetTitle("")
  return tmpM, tmpS

def setNormalColorTable():
  rArray = array.array('d',[0.0,1.0,1.0])
  gArray = array.array('d',[1.0,1.0,0.0])
  bArray = array.array('d',[0.0,0.0,0.0])
  stopArray = array.array('d',[0.,0.5,1.])
  nTabColors = 500
  root.TColor.CreateGradientColorTable(len(stopArray),
            stopArray,rArray,gArray,bArray,nTabColors
         )
def setInvertColorTable():
  rArray = array.array('d',[1.0,1.0,0.0])
  gArray = array.array('d',[0.0,1.0,1.0])
  bArray = array.array('d',[0.0,0.0,0.0])
  stopArray = array.array('d',[0.,0.5,1.])
  nTabColors = 500
  root.TColor.CreateGradientColorTable(len(stopArray),
            stopArray,rArray,gArray,bArray,nTabColors
         )

def setStyle():
  gStyle.SetCanvasColor(0)
  gStyle.SetCanvasBorderSize(10)
  gStyle.SetCanvasBorderMode(0)
  gStyle.SetCanvasDefH(700)
  gStyle.SetCanvasDefW(700)

  gStyle.SetPadColor       (0)
  gStyle.SetPadBorderSize  (10)
  gStyle.SetPadBorderMode  (0)
  gStyle.SetPadBottomMargin(0.13)
  gStyle.SetPadTopMargin   (0.08)
  gStyle.SetPadLeftMargin  (0.15)
  gStyle.SetPadRightMargin (0.05)
  gStyle.SetPadGridX       (0)
  gStyle.SetPadGridY       (0)
  gStyle.SetPadTickX       (1)
  gStyle.SetPadTickY       (1)

  gStyle.SetFrameFillStyle ( 0)
  gStyle.SetFrameFillColor ( 0)
  gStyle.SetFrameLineColor ( 1)
  gStyle.SetFrameLineStyle ( 0)
  gStyle.SetFrameLineWidth ( 1)
  gStyle.SetFrameBorderSize(10)
  gStyle.SetFrameBorderMode( 0)

  gStyle.SetNdivisions(505)

  gStyle.SetLineWidth(2)
  gStyle.SetHistLineWidth(2)
  gStyle.SetFrameLineWidth(2)
  gStyle.SetLegendFillColor(root.kWhite)
  gStyle.SetLegendFont(42)
  gStyle.SetMarkerSize(1.2)
  gStyle.SetMarkerStyle(20)
  gStyle.SetHistLineColor(1)
 
  gStyle.SetLabelSize(0.040,"X")
  gStyle.SetLabelSize(0.040,"Y")

  gStyle.SetLabelOffset(0.010,"X")
  gStyle.SetLabelOffset(0.010,"Y")
 
  gStyle.SetLabelFont(42,"X")
  gStyle.SetLabelFont(42,"Y")
 
  gStyle.SetTitleBorderSize(0)
  gStyle.SetTitleFont(42)
  gStyle.SetTitleFont(42,"X")
  gStyle.SetTitleFont(42,"Y")

  gStyle.SetTitleSize(0.045,"X")
  gStyle.SetTitleSize(0.045,"Y")
 
  gStyle.SetTitleOffset(1.4,"X")
  gStyle.SetTitleOffset(1.4,"Y")
 
  gStyle.SetTextSize(0.055)
  gStyle.SetTextFont(42)
 
  gStyle.SetOptStat(0)
  setNormalColorTable()
  
setStyle()

def setHistTitles(hist,xlabel,ylabel):
    hist.GetXaxis().SetTitle(xlabel)
    hist.GetYaxis().SetTitle(ylabel)

def makeWeightHist(f1,canvas,leg):
  firstHist = True
  canvas.cd()
  canvas.SetLogy()
  colorsList = [1,2,3,4,5,6,7,8]
  nColors = len(colorsList)
  iDir = 0
  leg.Clear()
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  tmpList = []
  for dirName in f1.GetListOfKeys():
    tmpList.append(dirName)
  tmpList.reverse()
  for dirName in tmpList:
    print(dirName.GetName())
    if(re.search(r"data",dirName.GetName())):
	continue
    directory = dirName.ReadObj()
    for histKey in directory.GetListOfKeys():
      if(histKey.GetName()=="hWeight"):
        hist = histKey.ReadObj()
	hist.UseCurrentStyle()
	hist.SetLineColor(colorsList[iDir % nColors])
	hist.SetMarkerColor(colorsList[iDir % nColors])
	allIntegral = hist.Integral(0,hist.GetNbinsX()+1)
	integral = hist.Integral()
	if integral > 0.0:
	  print("Fraction Outside of bounds: %f" % (allIntegral/integral-1.0))
	  #hist.Scale(1.0/allIntegral)
	  hist.Scale(1.0/integral)
	else:
	  leg.AddEntry(hist,dirName.GetName(),"lep")
	if(firstHist):
	  firstHist=False
	  hist.GetYaxis().SetTitle("Fraction of Events")
	  hist.GetXaxis().SetTitle("Event Weight")
	  #hist.GetXaxis().SetRangeUser(0.0,1.0)
	  hist.Draw()
	else:
	  hist.Draw("same")
    iDir += 1
  leg.Draw("same")

class DataMCStack:
  def __init__(self, mcHistList, dataHist, canvas, xtitle, ytitle="", drawStack=True,nDivX=7,xlimits=[],showOverflow=False,lumi=5.0,logy=False,signalsNoStack=[],showCompatabilityTests=True,integralPlot=False,energyStr="8TeV",ylimits=[],ylimitsRatio=[],pullType="",doMCErrors=False,showPullStats=False,yMaxVals=[],yMaxXRanges=[],mcVariations=None,scaleMC2Data=False):
    nBinsX = dataHist.GetNbinsX()
    self.xlimits = xlimits
    self.ylimits = ylimits
    self.logy = logy
    self.nBinsX = nBinsX
    self.dataHist = dataHist
    self.canvas = canvas
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)
    self.mcVarHist = None
    setYLimitsAuto = getattr(self,"setYLimitsAuto")
    if ytitle=="":
      ytitle="Events/%s" % (getBinWidthStr(dataHist))
    for mcHist in mcHistList:
      #print("nBinsX data: %i, mc: %i" % (nBinsX,mcHist.GetNbinsX()))
      assert(nBinsX == mcHist.GetNbinsX())
    for sigHist in signalsNoStack:
      assert(nBinsX == sigHist.GetNbinsX())

    if integralPlot:
      dataHist = getIntegralHist(dataHist,True)
      self.dataHist = dataHist
      newMcHistList = []
      for i in mcHistList:
        newMcHistList.append(getIntegralHist(i))
      mcHistList = newMcHistList
      newSigHistList = []
      for i in signalsNoStack:
        newSigHistList.append(getIntegralHist(i))
      signalsNoStack = newSigHistList
      ytitle = "Integral of "+ytitle+" #geq X"
    self.signalsNoStack = signalsNoStack
    self.mcHistList = mcHistList
    self.dataHist = dataHist

    self.nDataEvents = dataHist.Integral(0,dataHist.GetNbinsX()+1)
    self.mc2DataSF = 1.
    if scaleMC2Data:
      tmpMCSum = 0.
      for mcHist in mcHistList:
        tmpMCSum += mcHist.Integral(0,mcHist.GetNbinsX()+1)
      self.mc2DataSF = float(self.nDataEvents)/tmpMCSum
      print("DataMC SF: %.2f" % self.mc2DataSF)

    # Make MC Stack/sumHist
    self.stack = root.THStack()
    self.mcSumHist = dataHist.Clone("mcSumHist"+dataHist.GetName())
    self.mcSumHist.Reset()
    for mcHist in mcHistList:
      mcHist.SetMaximum(1e12)
      mcHist.SetMinimum(1e-12)
      mcHist.SetLineColor(mcHist.GetFillColor())
      if showOverflow:
        showHistOverflow(mcHist)
      mcHist.Scale(self.mc2DataSF)
      self.mcSumHist.Add(mcHist)
      self.stack.Add(mcHist)

    if showOverflow:
        showHistOverflow(dataHist)

    self.doMCVariations(mcVariations)

    self.mcSumHist.SetFillColor(root.kGray+3)
    self.mcSumHist.SetFillStyle(3254)
    self.mcSumHist.SetMarkerSize(0)
    if doMCErrors and drawStack:
        self.mcSumHist.SetLineStyle(0)

    self.nMCEvents = self.mcSumHist.Integral(0,self.mcSumHist.GetNbinsX()+1)

    # Get chi^2 Prob Data/MC
    self.normchi2 = dataHist.Chi2Test(self.mcSumHist,"UW CHI2/NDF")
    self.chi2Prob = dataHist.Chi2Test(self.mcSumHist,"UW")
    self.KSProb = dataHist.KolmogorovTest(self.mcSumHist)
    if self.mcVarHist != None:
      self.normchi2 = dataHist.Chi2Test(self.mcVarHist,"UW CHI2/NDF")
      self.chi2Prob = dataHist.Chi2Test(self.mcVarHist,"UW")
      self.KSProb = dataHist.KolmogorovTest(self.mcVarHist)
    if self.chi2Prob < 1e-20:
        self.chi2Prob = 0.0
    if self.KSProb < 1e-20:
        self.KSProb = 0.0

    # Make Pull Hist
    self.pullList = []
    self.pullHist = dataHist.Clone("pullHist"+dataHist.GetName())
    self.pullHist.Reset()
    self.oneGraph = root.TGraph()
    self.oneGraph.SetLineWidth(2)
    self.oneGraph.SetLineStyle(2)
    iGraph = 0
    for i in range(0,self.pullHist.GetNbinsX()+2):
      nData = dataHist.GetBinContent(i)
      nMC = self.mcSumHist.GetBinContent(i)
      error = dataHist.GetBinError(i)
      errorMC = self.mcSumHist.GetBinError(i)
      if self.mcVarHist != None:
        errorMC = self.mcVarHist.GetBinError(i)
      pull = 0.0
      ratio = 0.0
      ratioErr = 0.0
      self.oneGraph.SetPoint(iGraph,dataHist.GetXaxis().GetBinCenter(i),1.0)
      iGraph += 1
      if error != 0.0:
        if pullType=="adrian1":
          pull = (nData -nMC)/nData
        else:
          pull = (nData -nMC)/error
      if pullType=="pullMC":
        if errorMC != 0.0:
          pull = (nData -nMC)/errorMC
        else:
          pull = 0.0
      if nMC != 0.0:
        ratio = nData/nMC
        ratioErr = error/nMC
      if pullType=="ratio":
        self.pullHist.SetBinContent(i,ratio)
        self.pullHist.SetBinError(i,ratioErr)
        #print("nData: {0:.2f} +/- {1:.2f}, nMC: {2:.2f}, ratio: {3:.2f} +/- {4:.2f}".format(nData,error,nMC,ratio,ratioErr))
      else:
        self.pullHist.SetBinContent(i,pull)
        #print("nData: %f, nMC: %f, error: %f, pull: %f" % (nData,nMC,error,pull))
      #pullDistribution
      if pullType == "pullMC":
        if errorMC != 0.0:
          self.pullList.append((nData -nMC)/errorMC)
      else:
        if error != 0.0:
          self.pullList.append((nData -nMC)/error)
    #print getattr(self,"getPullDistributionParams")(self.pullList)

    #Find Maximum y-value
    if xlimits != []:
      self.mcSumHist.GetXaxis().SetRangeUser(*xlimits)
      self.dataHist.GetXaxis().SetRangeUser(*xlimits)
    mcMax = self.mcSumHist.GetMaximum()
    if self.mcVarHist != None:
      mcMax = self.mcSumHist.GetMaximum()
    dataMaxBin = self.dataHist.GetMaximumBin()
    dataMax = dataHist.GetBinContent(dataMaxBin)+dataHist.GetBinError(dataMaxBin)
    ymax = 0.0
    if mcMax > dataMax:
       ymax = mcMax
    else:
       ymax = dataMax
    self.ymax = ymax
  
    #Setup Canvas
    canvas.cd()
    self.pad1Top = 0.98
    self.pad1Bot = 0.30
    self.pad1Right = 0.98
    self.pad1Left = 0.02
    pad1 = root.TPad("pad1"+dataHist.GetName(),"",0.02,0.30,0.98,0.98,0)
    pad2 = root.TPad("pad2"+dataHist.GetName(),"",0.02,0.01,0.98,0.29,0)
    self.pad1 = pad1
    self.pad2 = pad2
  
    pad1.SetBottomMargin(0.005);
    pad2.SetTopMargin   (0.005);
    pad2.SetBottomMargin(0.33);
    """
    pad1.SetBottomMargin(0.01);
    pad2.SetTopMargin   (0.3);
    pad2.SetBottomMargin(0.33);
    """
    canvas.SetLogy(0)
    pad2.SetLogy(0)
    if logy:
        pad1.SetLogy(1)
    else:
        pad1.SetLogy(0)
  
    pad1.Draw() # Projections pad
    pad2.Draw() # Residuals   pad

    pad1Width = pad1.XtoPixel(pad1.GetX2())
    pad1Height = pad1.YtoPixel(pad1.GetY1())
    pad2Height = pad2.YtoPixel(pad2.GetY1())
    #pad1ToPad2FontScalingFactor = float(pad1Width)/pad2Height
    pad1ToPad2FontScalingFactor = float(pad1Height)/pad2Height
  
    # Main Pad
    pad1.cd();
    xAxis = None
    yAxis = None
    histForAxis = None
    if len(self.ylimits)==2:
      ylimits[0] += 1e-3
      histForAxis = root.TH2F(dataHist.GetName()+"ForAxis","",1,xlimits[0],xlimits[1],1,self.ylimits[0],self.ylimits[1])
    elif self.logy:
      histForAxis = root.TH2F(dataHist.GetName()+"ForAxis","",1,xlimits[0],xlimits[1],1,0.1,ymax*2.0)
    else:
      histForAxis = root.TH2F(dataHist.GetName()+"ForAxis","",1,xlimits[0],xlimits[1],1,1e-3,ymax*1.05)
    self.histForAxis = histForAxis
    self.histForAxis.Draw()
    self.mcSumHist.Draw("e1same")
    #self.canvas.SaveAs("debug.png")
    if len(self.ylimits)!=2:
      setYLimitsAuto(yMaxXRanges,yMaxVals,self.ymax)
    self.histForAxis.Draw()
    self.histForAxis.GetXaxis().SetTitle("")
    self.histForAxis.GetXaxis().SetLabelSize(0)
    self.histForAxis.GetYaxis().SetTitle(ytitle)
    self.histForAxis.GetYaxis().SetLabelSize(0.050)
    self.histForAxis.GetYaxis().SetTitleSize(0.055)
    self.histForAxis.GetXaxis().SetNdivisions(nDivX)
    self.histForAxis.GetXaxis().SetTitleColor(0)
    self.histForAxis.GetXaxis().SetLabelColor(0)
    if drawStack:
      self.stack.Draw("hist same")
      if doMCErrors:
        if self.mcVarHist != None:
          self.mcVarHist.Draw("e2same")
        self.mcSumHist.Draw("e2same")
      pad1.Update()
    else:
      self.mcSumHist.SetFillColor(856)
      self.mcSumHist.SetLineColor(856)
      self.mcSumHist.SetMarkerColor(856)
      self.mcSumHist.SetFillStyle(1001)
      self.mcSumHist.Draw("histo b")
    for sigHist in signalsNoStack:
      sigHist.Draw("histo same")
    dataHist.Draw("pe same")

    pad1.RedrawAxis() # Updates Axis Lines
  
    # Pulls Pad
    pad2.cd()
    self.pullHist.SetTitle("")
    if xlimits != []:
      self.pullHist.GetXaxis().SetRangeUser(*xlimits)
    self.pullHist.GetXaxis().SetTitle(xtitle)
    self.pullHist.GetXaxis().CenterTitle(1)
    self.pullHist.GetXaxis().SetNdivisions(nDivX)
    self.pullHist.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
    self.pullHist.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist.SetLineColor(root.kBlue)
    self.pullHist.SetLineStyle(1)
    self.pullHist.SetLineWidth(2)
    if pullType=="adrian1":
      self.pullHist.GetYaxis().SetTitle("#frac{Data-MC}{Data}")
    elif pullType=="pullMC":
      self.pullHist.GetYaxis().SetTitle("#frac{Data-MC}{\sigma_{MC}}")
    else:
      self.pullHist.GetYaxis().SetTitle("#frac{Data-MC}{\sigma_{Data}}")
    self.pullHist.GetYaxis().SetTitleSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().SetLabelSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().CenterTitle(1)
    self.pullHist.GetXaxis().SetTitleOffset(0.75*self.pullHist.GetXaxis().GetTitleOffset())
    self.pullHist.GetYaxis().SetTitleOffset(0.70)
    self.pullHist.SetFillColor(856)
    self.pullHist.SetFillStyle(1001)
    if len(ylimitsRatio) == 2:
      ylimitsRatio[0] += 1e-3
      ylimitsRatio[1] -= 1e-3
      self.pullHist.GetYaxis().SetRangeUser(*ylimitsRatio)

    if pullType=="ratio":
      #pad2.SetGridy(1)
      self.pullHist.GetYaxis().SetTitle("#frac{Data}{MC}")
      self.pullHist.Draw("")
      self.oneGraph.Draw()
      self.pullHist.Draw("same")
    else:
      self.pullHist.Draw("histo")

    if showCompatabilityTests:
      self.problatex = root.TLatex()
      self.problatex.SetNDC()
      self.problatex.SetTextFont(root.gStyle.GetLabelFont())
      self.problatex.SetTextSize(self.pullHist.GetYaxis().GetLabelSize())
      self.problatex.SetTextAlign(12)
      yToDraw = 0.41 #bottom
      yToDraw = 0.92 #top
      #self.problatex.DrawLatex(0.18,yToDraw,"KS Prob: {0:.3g}".format(self.KSProb))
      self.problatex.DrawLatex(0.18,yToDraw,"#chi^{2}/NDF: %.3g" % (self.normchi2))
      self.problatex.DrawLatex(0.18,yToDraw-0.08,"#chi^{2}  Prob: %.3g" % (self.chi2Prob))

    pad2.Update()
    pad2.GetFrame().DrawClone()
    pad2.RedrawAxis() # Updates Axis Lines
  
    canvas.cd()
    self.tlatex.DrawLatex(0.33,0.96,PRELIMINARYSTRING)
    self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=%s, L=%.1f fb^{-1}" % (energyStr,lumi))

  def getPullDistributionParams(self,pullList):
    pull = root.RooRealVar("pull","pull",-20,20)
    mean = root.RooRealVar("mean","pull Mean",0.0,-20,20)
    sigma = root.RooRealVar("sigma","pull sigma",1.0,0.01,20)
    self.pullGaus = root.RooGaussian("pullGaus","pullGaus",pull,mean,sigma)
    self.pullDS = root.RooDataSet("pullDS","pullDS",root.RooArgSet(pull))
    for i in pullList:
      pull.setVal(i)
      self.pullDS.add(root.RooArgSet(pull))
    self.pullFR = self.pullGaus.fitTo(self.pullDS,PRINTLEVEL)
    self.pullMean = mean
    self.pullSigma = sigma
    meanStr = "<Pull> = %.2f #pm %.2f" % (mean.getVal(), mean.getError())
    sigmaStr = "#sigma(Pull) = %.2f #pm %.2f" % (sigma.getVal(), sigma.getError())

    frame = pull.frame(root.RooFit.Bins(20))
    self.pullDS.plotOn(frame)
    self.pullGaus.plotOn(frame)
    frame.Draw()
    self.canvas.SaveAs("pullDist"+self.dataHist.GetName()+".png")
    return meanStr, sigmaStr

  def getXNDC(self,x):
    minX = self.pad1.GetX1()
    maxX = self.pad1.GetX2()
    result=(x-minX)/(maxX-minX)
    return result
  def getYNDC(self,y):
    minY = self.pad1.GetY1()
    maxY = self.pad1.GetY2()
    result=(y-minY)/(maxY-minY)
    return result
  def getXUser(self,x):
    minX = self.pad1.GetX1()
    maxX = self.pad1.GetX2()
    result=x*(maxX-minX)+minX
    return result
  def getYUser(self,y):
    minY = self.pad1.GetY1()
    maxY = self.pad1.GetY2()
    result=y*(maxY-minY)+minY
    #print "running getYUser with: %.2f" % y
    #print "  minY: %.2f" % minY
    #print "  maxY: %.2f" % maxY
    #print "  result: %.2f" % result
    return result
  def setYLimitsAuto(self,rangesNDC,yNDCLimits,yMaxCurrent):
    #self.canvas.SaveAs("before_"+str(int(time.time()*100))+".png")
    #print("Running setYLimitsAuto...")
    self.pad1.Update()
    self.canvas.Update()
    getXUser = getattr(self,"getXUser")
    getYUser = getattr(self,"getYUser")
    setYLimitsAuto = getattr(self,"setYLimitsAuto")
    self.pad1.cd()
    ranges = [[getXUser(i[0]),getXUser(i[1])] for i in rangesNDC]
    yLimitsScaleFactor = 1.0
    if self.logy:
      yLimitsScaleFactor = 0.75
    yLimits = [getYUser(i)*yLimitsScaleFactor for i in yNDCLimits]
    maxPoints = []
    xAxis = self.mcSumHist.GetXaxis()
    #print("yMaxCurrent: %.2f " % (yMaxCurrent))
    for r,yLim in zip(ranges,yLimits):
      maxY = 0.0
      for i in range(1,xAxis.GetNbins()+1):
        if xAxis.GetBinUpEdge(i) >= r[0] and xAxis.GetBinLowEdge(i) <= r[1]:
          y = self.mcSumHist.GetBinContent(i)
          yErrTmp = self.mcSumHist.GetBinError(i)
          yErr2Tmp = 0.
          if self.mcVarHist != None:
            yErr2Tmp = self.mcVarHist.GetBinError(i)
          y += max(yErrTmp,yErr2Tmp)
          maxY = max(y,maxY)
      maxPoints += [maxY]
    rescale = 0.0
    if self.logy:
      newMaxPoints = []
      for x in maxPoints:
        if x>0.:
          newMaxPoints += [log10(x)]
        else:
          newMaxPoints += [0.]
      maxPoints = newMaxPoints
    for yLim,maxY in zip(yLimits,maxPoints):
      #print("yLim: %.2f maxY: %.2f" % (yLim, maxY))
      if maxY > yLim:
        rescaleTmp = (maxY/yLim)
        if rescaleTmp > rescale:
          rescale = rescaleTmp
    if rescale == 0.0:
        self.ymax = yMaxCurrent*1.1
        return
    if self.logy:
      rescale = 10**rescale*5.
    #print(rescale)
    newYMax = yMaxCurrent*rescale*1.5
    newYMin = 1e-3
    if self.logy:
      newYMin = 0.1
    self.histForAxis = root.TH2F(self.histForAxis.GetName()+"ForAxis","",1,self.xlimits[0],self.xlimits[1],1,newYMin,newYMax)
    self.histForAxis.Draw("")
    self.mcSumHist.Draw("e1 same")
    #self.canvas.SaveAs("after_"+str(int(time.time()*100))+".png")
    setYLimitsAuto(rangesNDC,yNDCLimits,newYMax)

  def doMCVariations(self,mcVariations):
    self.mcVarHist = None
    if mcVariations==None:
      return
    for key in mcVariations:
      for hist in mcVariations[key]:
        hist.Scale(self.mc2DataSF)
    errorTypes = set()
    for key in mcVariations:
      key = re.sub("Up$","",key)
      key = re.sub("Down$","",key)
      if not key in errorTypes:
        errorTypes.add(key)
    mcSumVariations = {}
    for key in mcVariations:
      if len(mcVariations[key])==0:
        continue
      sumHist = mcVariations[key][0].Clone()
      sumHist.Reset()
      for h in mcVariations[key]:
        sumHist.Add(h)
      mcSumVariations[key] = sumHist
    self.mcVarHist = self.mcSumHist.Clone(self.mcSumHist.GetName()+"_mcVariations")
    for iBin in range(1,self.mcVarHist.GetNbinsX()+1):
      nom = self.mcVarHist.GetBinContent(iBin)
      err2 = self.mcVarHist.GetBinError(iBin)**2
      for eBase in errorTypes:
        errUp = mcSumVariations[eBase+"Up"].GetBinContent(iBin)
        errDown = mcSumVariations[eBase+"Down"].GetBinContent(iBin)
        errUp = abs(nom-errUp)
        errDown = abs(nom-errDown)
        if errUp > errDown:
            err2 += errUp**2
        else:
            err2 += errDown**2
      err = sqrt(err2)
      self.mcVarHist.SetBinError(iBin,err)
    self.mcVarHist.SetFillColor(root.kRed)
    self.mcVarHist.SetFillStyle(3245)
    self.mcVarHist.SetMarkerSize(0)
    self.mcVarHist.SetLineStyle(0)


class CompareTwoHists:
  def __init__(self, hist1,hist2, canvas, xtitle, ytitle="Events",nDivX=7,nDivPullY=5,xlimits=[],ylimits=[],pullHistRangeY=[0.0,2.0],energyStr="8TeV",lumi=19.4):
    nBinsX = hist1.GetNbinsX()
    assert(nBinsX == hist2.GetNbinsX())
    self.nBinsX = nBinsX
    self.hist1 = hist1
    self.hist2 = hist2
    self.canvas = canvas
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)

    if xlimits != []:
      self.hist1.GetXaxis().SetRangeUser(*xlimits)
      self.hist2.GetXaxis().SetRangeUser(*xlimits)
  
    # Make Pull Hist
    self.pullHist = hist1.Clone("pullHist"+hist1.GetName())
    self.pullHist.Reset()
    self.pullHist.SetLineColor(hist2.GetLineColor())
    self.pullHist.SetMarkerColor(hist2.GetMarkerColor())
    self.pullErrorBand = root.TGraphAsymmErrors()
    self.pullErrorBand.SetFillColor(856)
    self.pullErrorBand.SetFillStyle(1001)
    self.pullErrorBand.SetLineStyle(2)
    self.pullErrorBand.SetLineColor(root.kBlack)
    for i in range(0,nBinsX+2):
      nhist1 = hist1.GetBinContent(i)
      nhist2 = hist2.GetBinContent(i)
      nhist1Err = hist1.GetBinError(i)
      nhist2Err = hist2.GetBinError(i)
      ratio = 0.0
      ratioErr = 0.0
      if nhist1 != 0.0 and nhist2 != 0.0:
        ratio = nhist2/nhist1
        ratioErr = nhist2/nhist1 * sqrt((nhist1Err/nhist1)**2+(nhist2Err/nhist2)**2)
      self.pullHist.SetBinContent(i,ratio)
      self.pullHist.SetBinError(i,ratioErr)
      tmpAxis = hist1.GetXaxis()
      self.pullErrorBand.SetPoint(i,tmpAxis.GetBinCenter(i),1.0)
      if nhist1 != 0.0:
        self.pullErrorBand.SetPointError(i,tmpAxis.GetBinLowEdge(i),tmpAxis.GetBinUpEdge(i),
                                                            nhist1Err/nhist1,nhist1Err/nhist1)
      else:
        self.pullErrorBand.SetPointError(i,tmpAxis.GetBinLowEdge(i),tmpAxis.GetBinUpEdge(i),0.,0.)
      #print("nData: %f, nMC: %f, error: %f, pull: %f" % (nData,nMC,error,pull))

    firstVizBin = self.pullHist.GetXaxis().GetFirst()
    lastVizBin = self.pullHist.GetXaxis().GetLast()
    for i in range(0,firstVizBin):
      self.pullErrorBand.SetPointEYhigh(i,self.pullErrorBand.GetErrorYhigh(firstVizBin))
      self.pullErrorBand.SetPointEYlow(i,self.pullErrorBand.GetErrorYlow(firstVizBin))
    for i in range(lastVizBin+1,nBinsX+2):
      self.pullErrorBand.SetPointEYhigh(i,self.pullErrorBand.GetErrorYhigh(lastVizBin))
      self.pullErrorBand.SetPointEYlow(i,self.pullErrorBand.GetErrorYlow(lastVizBin))

    #Find Maximum y-value
    max1 = hist1.GetMaximum()
    max2 = hist2.GetMaximum()
    max1Bin = hist1.GetMaximumBin()
    max2Bin = hist2.GetMaximumBin()
    max1 = hist1.GetBinContent(max1Bin)+hist1.GetBinError(max1Bin)
    max2 = hist2.GetBinContent(max2Bin)+hist2.GetBinError(max2Bin)
    ymax = 0.0
    if max1 > max2:
       ymax = max1
    else:
       ymax = max2
    self.ymax = ymax

    #Setup Canvas
    canvas.cd()
    canvas.Clear()
    pad1 = root.TPad("pad1"+hist1.GetName(),"",0.02,0.30,0.98,0.98,0)
    pad2 = root.TPad("pad2"+hist1.GetName(),"",0.02,0.01,0.98,0.29,0)
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
  
    # Main Pad
    pad1.cd();
    self.hist2.GetXaxis().SetTitle("")
    self.hist2.GetXaxis().SetLabelSize(0)
    self.hist2.GetYaxis().SetTitle(ytitle)
    self.hist2.GetYaxis().SetLabelSize(0.050)
    self.hist2.GetYaxis().SetTitleSize(0.055)
    self.hist2.GetXaxis().SetNdivisions(nDivX)
    if ylimits==[]:
      self.hist2.GetYaxis().SetRangeUser(0.0,ymax*1.04)
    else:
      self.hist2.GetYaxis().SetRangeUser(*ylimits)
    self.hist2.Draw("")
    self.hist1.Draw("same")
  
    # Pulls Pad
    pad2.cd()
    self.pullHist.SetTitle("")
    if xlimits != []:
      self.pullHist.GetXaxis().SetRangeUser(*xlimits)
    self.pullHist.GetXaxis().SetTitle(xtitle)
    self.pullHist.GetXaxis().CenterTitle(1)
    self.pullHist.GetXaxis().SetNdivisions(nDivX)
    self.pullHist.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
    self.pullHist.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().SetTitle("Ratio")
    self.pullHist.GetYaxis().SetTitleSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().SetLabelSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist.GetYaxis().CenterTitle(1)
    self.pullHist.GetYaxis().SetTitleOffset(0.5)
    self.pullHist.GetYaxis().SetNdivisions(nDivPullY)
    self.pullHist.GetYaxis().SetRangeUser(*pullHistRangeY)
    self.pullHist.GetXaxis().SetTitleOffset(0.75*self.pullHist.GetXaxis().GetTitleOffset())
    self.pullHist.GetYaxis().SetTitleOffset(0.70)
    self.pullHist.SetFillColor(856)
    self.pullHist.SetFillStyle(1001)
    self.pullHist.Draw("")
    self.pullErrorBand.Draw("3")
    self.pullErrorBand.Draw("LX")
    #self.pullErrorBandLine = self.pullErrorBand.Clone(self.pullErrorBand.GetName()+"Line")
    #self.pullErrorBandLine.SetFillStyle(0)
    #self.pullErrorBandLine.Draw("same HIST L")
    self.pullHist.Draw("same")

    pad1.RedrawAxis() # Updates Axis Lines
    pad2.RedrawAxis() # Updates Axis Lines
  
    #canvas.cd()
    pad1.cd()
    self.tlatex.DrawLatex(0.33,0.96,PRELIMINARYSTRING)
    self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=%s, L=%.1f fb^{-1}" % (energyStr,lumi))

class CompareTwoHistsAndData:
  def __init__(self, hist1,hist2, data, canvas, xtitle, ytitle="Events",nDivX=7,nDivPullY=5,xlimits=[],ylimits=[],pullHistRangeY=[0.0,2.0],isPreliminary=True,is7TeV=False,lumi=5.0,logy=False,integralPlot=False,energyStr="8TeV"):
    nBinsX = hist1.GetNbinsX()
    assert(nBinsX == hist2.GetNbinsX())
    assert(nBinsX == data.GetNbinsX())
    self.nBinsX = nBinsX
    self.hist1 = hist1
    self.hist2 = hist2
    self.data = data
    self.canvas = canvas
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)

    if xlimits != []:
      self.hist1.GetXaxis().SetRangeUser(*xlimits)
      self.hist2.GetXaxis().SetRangeUser(*xlimits)
  
    # Make Pull Hist
    self.pullHist1 = hist1.Clone("pullHist1"+data.GetName())
    self.pullHist1.Reset()
    self.pullHist2 = hist2.Clone("pullHist2"+data.GetName())
    self.pullHist2.Reset()
    for i in range(1,self.pullHist1.GetNbinsX()):
      nData = data.GetBinContent(i)
      nMC1 = hist1.GetBinContent(i)
      nMC2 = hist2.GetBinContent(i)
      error = data.GetBinError(i)
      pull1 = 0.0
      pull2 = 0.0
      if error != 0.0:
        pull1 = (nMC1 - nData)/error
        pull2 = (nMC2 - nData)/error
      self.pullHist1.SetBinContent(i,pull1)
      self.pullHist2.SetBinContent(i,pull2)

    #Find Maximum y-value
    max1 = hist1.GetMaximum()
    max2 = hist2.GetMaximum()
    max1Bin = hist1.GetMaximumBin()
    max2Bin = hist2.GetMaximumBin()
    max1 = hist1.GetBinContent(max1Bin)+hist1.GetBinError(max1Bin)
    max2 = hist2.GetBinContent(max2Bin)+hist2.GetBinError(max2Bin)
    ymax = 0.0
    if max1 > max2:
       ymax = max1
    else:
       ymax = max2
    self.ymax = ymax

    #Find min/max pulls
    max1 = self.pullHist1.GetMaximum()
    max2 = self.pullHist2.GetMaximum()
    pullmax = 0.0
    if max1 > max2:
       pullmax = max1
    else:
       pullmax = max2
    self.pullmax = pullmax
    min1 = self.pullHist1.GetMinimum()
    min2 = self.pullHist2.GetMinimum()
    pullmin = 0.0
    if min1 < min2:
       pullmin = min1
    else:
       pullmin = min2
    self.pullmin = pullmin

    #Setup Canvas
    canvas.cd()
    canvas.Clear()
    pad1 = root.TPad("pad1"+hist1.GetName(),"",0.02,0.30,0.98,0.98,0)
    pad2 = root.TPad("pad2"+hist1.GetName(),"",0.02,0.01,0.98,0.29,0)
  
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
  
    # Main Pad
    pad1.cd();
    self.hist2.SetTitle("")
    self.hist2.GetXaxis().SetTitle("")
    self.hist2.GetXaxis().SetLabelSize(0)
    self.hist2.GetYaxis().SetTitle(ytitle)
    self.hist2.GetYaxis().SetLabelSize(0.050)
    self.hist2.GetYaxis().SetTitleSize(0.055)
    self.hist2.GetXaxis().SetNdivisions(nDivX)
    if ylimits==[]:
      self.hist2.GetYaxis().SetRangeUser(0.0,ymax*1.04)
    else:
      self.hist2.GetYaxis().SetRangeUser(*ylimits)
    self.hist2.SetFillStyle(0)
    self.hist1.SetFillStyle(0)
    self.hist2.Draw("hist")
    self.hist1.Draw("hist same")
    self.data.Draw("pe same")
  
    # Pulls Pad
    pad2.cd()
    self.pullHist1.SetTitle("")
    if xlimits != []:
      self.pullHist1.GetXaxis().SetRangeUser(*xlimits)
    self.pullHist0 = self.pullHist1.Clone("pullHist0")
    self.pullHist0.Reset()
    self.pullHist0.SetLineColor(1)
    self.pullHist0.SetLineStyle(2)
    self.pullHist0.SetFillStyle(0)
    self.pullHist1.GetXaxis().SetTitle(xtitle)
    self.pullHist1.GetXaxis().CenterTitle(1)
    self.pullHist1.GetXaxis().SetNdivisions(nDivX)
    self.pullHist1.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetYaxis().SetTitle("#frac{MC-Data}{Error}")
    self.pullHist1.GetYaxis().SetTitleSize(0.050*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetYaxis().SetLabelSize(0.040*pad1ToPad2FontScalingFactor)
    self.pullHist1.GetYaxis().CenterTitle(1)
    self.pullHist1.GetYaxis().SetTitleOffset(0.5)
    self.pullHist1.GetYaxis().SetNdivisions(nDivPullY)
    self.pullHist1.GetYaxis().SetRangeUser(pullmin*0.90,pullmax*1.1)
    self.pullHist1.GetXaxis().SetTitleOffset(0.75*self.pullHist1.GetXaxis().GetTitleOffset())
    #self.pullHist1.GetYaxis().SetTitleOffset(0.70)
    self.pullHist1.SetFillStyle(0)
    self.pullHist2.SetFillStyle(0)
    self.pullHist1.Draw("hist")
    self.pullHist0.Draw("hist same")
    self.pullHist1.Draw("hist same")
    self.pullHist2.Draw("hist same")
    pad2.Update()
    pad2.GetFrame().DrawClone()
  
    canvas.cd()
    self.tlatex.DrawLatex(0.33,0.96,PRELIMINARYSTRING)
    self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=%s, L=%.1f fb^{-1}" % (energyStr,lumi))

def jacknife(a,func,d=1):
  """
    Returns jacknife estimators of f, error on f, and bias of f 
    a must be a 1D list or array
    f must be a function that takes a list or array
    d is the number of samples to delete (as in delete-d jacknife)
  """
  a = numpy.array(a)
  estimate = func(a)
  n = len(a)
  if d < 1 or d > n:
    print("Error: jacknife: d="+str(d)+" is is out of bounds for n: "+str(n))
    sys.exit(1)
  elif d == 1:
    jnEstimates = numpy.zeros(n)
    bools = numpy.ones(n,dtype=numpy.bool)
    for i in range(n):
      bools[i] = False
      jnEstimates[i] = func(a[bools])
      bools[i] = True
    jnEstimate = numpy.sum(jnEstimates)/n
    jnBias = estimate-jnEstimate
    jnError = jnEstimates-jnEstimate
    jnError = numpy.power(jnError,2)
    jnError = numpy.sum(jnError)
    jnError *= (n-1)/float(n)
    jnError = numpy.sqrt(jnError)
    return jnEstimate, jnError, jnBias
  else:
    return None

def makeBootstrapHist(hist,outHist,entries=None):
 outHist.Reset()
 samples = entries
 if samples == None:
   integral = hist.Integral()
 for i in range(samples):
   outHist.Fill(hist.GetRandom())

def getMedianAndQuantileInterval(hist,amountForQuantile,doErrors=False):
    """ 
        Takes and input histogram, and an amount for the quantile 
        returns a list with: 
          [amountForQuantile, median, and 1.0-amountForQuantile] quantiles
    """
    quantilesToGet = array.array('d',[amountForQuantile,0.5,1.0-amountForQuantile])
    if doErrors:
      lowErrors = numpy.zeros(1000)
      medianErrors = numpy.zeros(1000)
      highErrors = numpy.zeros(1000)
      errors = []
      for i in range(1000):
        tmpHist = hist.Clone("tmpForQuantiles")
        makeBootstrapHist(hist,tmpHist,int(hist.GetEntries()))
        quantiles = array.array('d',[0.0,0.0,0.0])
        nQuantiles = tmpHist.GetQuantiles(3,quantiles,quantilesToGet)
        if(nQuantiles != 3):
          raise Exception("ROOT Hist Quantile Estimation didn't work!!")
        lowErrors[i] = quantiles[0]
        medianErrors[i] = quantiles[1]
        highErrors[i] = quantiles[2]
      quantiles = array.array('d',[0.0,0.0,0.0])
      nQuantiles = hist.GetQuantiles(3,quantiles,quantilesToGet)
      if(nQuantiles != 3):
          raise Exception("ROOT Hist Quantile Estimation didn't work!!")
      medianErros = []
      highErrors = []
      errors.append(numpy.std(lowErrors))
      errors.append(numpy.std(medianErrors))
      errors.append(numpy.std(highErrors))
      return list(quantiles), errors
    else:
      quantiles = array.array('d',[0.0,0.0,0.0])
      nQuantiles = hist.GetQuantiles(3,quantiles,quantilesToGet)
      if(nQuantiles != 3):
          raise Exception("ROOT Hist Quantile Estimation didn't work!!")
      return list(quantiles)

def sqrtThisHistogram(hist):
    """
        Sqrt's bin contents 
        of the input hist bin contents, properly treating the errors.
    """
    nBins = hist.GetNbinsX()

    for i in range(nBins+2):
      y = hist.GetBinContent(i)
      yErr = hist.GetBinError(i)
      if y < 0.0:
        print("Warning sqrtThisHIstogram: hist named %s bin %i has negative y value %f" % (hist.GetName(),i,y))
        hist.SetBinContent(i,0.0)
        hist.SetBinError(i,0.0)
        continue
      if y == 0.0:
        print("Warning sqrtThisHIstogram: hist named %s bin %i has zero y value" % (hist.GetName(),i))
        hist.SetBinContent(i,0.0)
        hist.SetBinError(i,0.0)
        continue
      hist.SetBinContent(i,sqrt(y))
      hist.SetBinError(i,yErr/(sqrt(2*y)))

def getSqrtCopyOfHistogram(hist):
    """
        Reterns a histogram of where the bin contents are the sqrt
        of the input hist bin contents, properly treating the errors.
    """
    outHist = hist.Clone(hist.GetName()+"SqrtHist")
    sqrtThisHistogram(outHist)
    return outHist

def drawSilly(isPreliminary=True,is7TeV=False):
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(0.05)
    tlatex.SetTextAlign(22)
    if isPreliminary:
      tlatex.DrawLatex(0.33,0.96,"CMS Preliminary")
    if is7TeV:
      tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=8 TeV, L=4.7 fb^{-1}")

def normalizeHist(hist):
  integral = hist.Integral(0,hist.GetNbinsX()+1)
  if integral != 0.0:
    hist.Scale(1.0/integral)

def showHistOverflow(hist):
  nBins = hist.GetNbinsX()

  overflow = hist.GetBinContent(nBins+1)
  overflowErr = hist.GetBinError(nBins+1)
  lastBin = hist.GetBinContent(nBins)
  lastBinErr = hist.GetBinError(nBins)

  hist.SetBinContent(nBins,lastBin+overflow)
  hist.SetBinError(nBins,sqrt(lastBinErr**2+overflowErr**2))

  underflow = hist.GetBinContent(0)
  underflowErr = hist.GetBinError(0)
  firstBin = hist.GetBinContent(1)
  firstBinErr = hist.GetBinError(1)

  hist.SetBinContent(1,firstBin+underflow)
  hist.SetBinError(1,sqrt(firstBinErr**2+underflowErr**2))

class PlotOfSlices:
  def __init__(self, hist2D, xtitle, ytitle, canvas, xlimits=[], ylimits=[],sliceLabelPrefix="",isPreliminary=True,is7TeV=False):
    canvas.cd(0)
    canvas.Clear()
    nBinsX = hist2D.GetNbinsX()
    nBinsY = hist2D.GetNbinsY()
    self.nBinsX = nBinsX
    self.nBinsY = nBinsY
    self.hist2D = hist2D
    self.canvas = canvas
    self.sliceLabelPrefix = sliceLabelPrefix
    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(0.05)
    self.tlatex.SetTextAlign(22)
    self.histList = []

    colorsListTmp = [root.kRed+1,root.kBlue+1,root.kGreen+1,root.kCyan,root.kMagenta+1]
    self.colorsList=[]
    for i in [0,-11,+2]:
        for j in range(len(colorsListTmp)):
            self.colorsList.append(colorsListTmp[j]+i)

    if xlimits != []:
      self.hist2D.GetXaxis().SetRangeUser(*xlimits)
    if ylimits != []:
      self.hist2D.GetYaxis().SetRangeUser(*ylimits)

    ymax = 0.0
    for i in range(nBinsX+2):
        tmpHist = root.TH1F(hist2D.GetName()+"_slice"+str(i),"",
                            nBinsY,hist2D.GetYaxis().GetXbins().GetArray())
        for j in range(nBinsY+2):
            tmpHist.SetBinContent(j,hist2D.GetBinContent(i,j))
        tmpMax = tmpHist.GetMaximum()
        if tmpMax > ymax:
            ymax = tmpMax
        tmpHist.SetLineColor(self.colorsList[i])
        self.histList.append(tmpHist)
    
    firstHist = self.histList[0]
    firstHist.SetTitle("")
    firstHist.GetXaxis().SetTitle(xtitle)
    firstHist.GetYaxis().SetTitle(ytitle)
    if xlimits != []:
        firstHist.GetXaxis().SetRangeUser(*xlimits)
    if ylimits==[]:
        firstHist.GetYaxis().SetRangeUser(0.0,ymax*1.05)
    else:
        firstHist.GetYaxis().SetRangeUser(*ylimits)

    firstHist.Draw("")
    for hist in self.histList[1:]:
        hist.Draw("same")

    if isPreliminary:
      self.tlatex.DrawLatex(0.33,0.96,"CMS Preliminary")
    if is7TeV:
      self.tlatex.DrawLatex(0.75,0.96,"#sqrt{s}=8 TeV, L=4.7 fb^{-1}")

    ## Lgend Part

    leg = root.TLegend(0.6,0.3,0.9,0.9)
    leg.SetLineColor(root.kWhite)
    leg.SetFillColor(root.kWhite)
    self.leg = leg
    xAxis = self.hist2D.GetXaxis()
    xBin = 0
    for hist in self.histList:
      tmpLabel = ""
      if xBin == 0:
        tmpLabel = "%s [0.0,%.1f]" % (sliceLabelPrefix,xAxis.GetBinUpEdge(xBin))
      elif xBin == nBinsX+1:
        tmpLabel = "%s [%.1f,#infty]" % (sliceLabelPrefix,xAxis.GetBinLowEdge(xBin))
      else:
        tmpLabel = "%s [%.1f,%.1f]" % (sliceLabelPrefix,xAxis.GetBinLowEdge(xBin),xAxis.GetBinUpEdge(xBin))
      leg.AddEntry(hist,tmpLabel,"l")
      xBin += 1
    leg.Draw("same")

def getIntegralHist(hist,setErrors=True):
  result = hist.Clone(hist.GetName()+"_Integral")
  nBins = result.GetNbinsX()
  for i in range(nBins+1):
    sumw = 0.0
    sumw2 = 0.0
    for j in range(i,nBins+2):
      sumw += result.GetBinContent(j)
      sumw2 += (result.GetBinError(j))**2
    result.SetBinContent(i,sumw)
    if setErrors:
        result.SetBinError(i,sumw2**0.5)
  return result

def hist2to1(hist):
  assert(hist.InheritsFrom("TH1"))
  result = None
  nBinsX = hist.GetNbinsX()
  nBinsY = hist.GetNbinsY()
  totalBins = (nBinsX+2)*(nBinsY+2) - 2 #include underflow/overflow
  if hist.InheritsFrom("TH2F"):
    result = root.TH1F(hist.GetName()+"_1d","",totalBins,0,totalBins)
  elif hist.InheritsFrom("TH2D"):
    result = root.TH1D(hist.GetName()+"_1d","",totalBins,0,totalBins)
  else:
    print("Error: hist2to1: Input hist must be TH2F or TH2D, exiting.")
    sys.exit(1)
  k = 0
  for i in range(nBinsX+2):
    for j in range(nBinsY+2):
      tmp = hist.GetBinContent(i,j)
      tmpErr = hist.GetBinError(i,j)
      result.SetBinContent(k,tmp)
      result.SetBinError(k,tmpErr)
      k += 1
  return result

def hist2to1CollapseY(hist,xcuts=[]):
  assert(hist.InheritsFrom("TH1"))
  result = None
  nBinsX = hist.GetNbinsX()
  nBinsY = hist.GetNbinsY()
  ymin = hist.GetYaxis().GetXmin()
  ymax = hist.GetYaxis().GetXmax()
  totalBins = (nBinsX+2)*(nBinsY+2) - 2 #include underflow/overflow
  if hist.InheritsFrom("TH2F"):
    result = root.TH1F(hist.GetName()+"_1d","",nBinsY,ymin,ymax)
  elif hist.InheritsFrom("TH2D"):
    result = root.TH1D(hist.GetName()+"_1d","",nBinsY,ymin,ymax)
  else:
    print("Error: hist2to1CollapseY: Input hist must be TH2F or TH2D, exiting.")
    sys.exit(1)
  minBinX = 0
  maxBinX = nBinsX+2
  if len(xcuts)==2:
    minBinX = hist.GetXaxis().FindBin(xcuts[0])
    maxBinX = hist.GetXaxis().FindBin(xcuts[1])
    if hist.GetXaxis().GetBinCenter(maxBinX)> xcuts[1]:
        maxBinX -= 1
  for j in range(nBinsY+2):
    tmpSum = 0.0
    tmpSumErr2 = 0.0
    for i in range(minBinX,maxBinX):
      tmp = hist.GetBinContent(i,j)
      tmpErr = hist.GetBinError(i,j)
      tmpSum += tmp
      tmpSumErr2 += tmpErr*tmpErr
    result.SetBinContent(j,tmpSum)
    result.SetBinError(j,sqrt(tmpSumErr2))
  return result

def shrinkTH1(hist,xlow,xhigh,deleteOld=False):
  assert(hist.InheritsFrom("TH1"))
  taxis=hist.GetXaxis()
  oldXlow=taxis.GetXmin()
  oldXhigh=taxis.GetXmax()
  assert(xlow >= oldXlow)
  assert(xhigh <= oldXhigh)
  lowBin = taxis.FindBin(xlow)
  highBin = taxis.FindBin(xhigh)
  if taxis.GetBinLowEdge(highBin)==float(xhigh):
    highBin -= 1
  xlow = taxis.GetBinLowEdge(lowBin)
  xhigh = taxis.GetBinUpEdge(highBin)
  oldN = hist.GetNbinsX()
  newN = int((xhigh-xlow)/(oldXhigh-oldXlow)*oldN)
  name = hist.GetName()
  title = hist.GetTitle()
  hist.SetName(name+"_Old")
  newHist = root.TH1F(name,title,newN,xlow,xhigh)
  newHist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
  newHist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
  for iOld,iNew in zip(range(lowBin,highBin+1),range(1,newN+1)):
    newHist.SetBinContent(iNew,hist.GetBinContent(iOld))
    newHist.SetBinError(iNew,hist.GetBinError(iOld))
  if deleteOld:
    hist.Delete()
  return newHist

def linearChi2(xList,yList,order,yErrList=None,funcName="poly"):
  assert(len(xList)==len(yList))
  xList = numpy.array(xList)
  yList = numpy.array(yList)

  poly = None
  stats = None
  weights = None
  if yErrList == None:
    weights = 1.0/numpy.sqrt(yList)
  else:
    assert(len(yErrList)==len(yList))
    yErrList = numpy.array(yErrList)
    weights = 1.0/yErrList
  if funcName == "poly":
    poly,stats = numpy.polynomial.polyfit(xList,yList,order,full=True,w=weights)
  else:
    poly,stats = numpy.polynomial.chebfit(xList,yList,order,full=True,w=weights)

#  # For Debug
#  fig = mpl.figure()
#  ax = fig.add_subplot(111)
#  ax.errorbar(xList,yList,yerr=1.0/weights,linestyle="None",color="k")
#
#  polyInst = None
#
#  if funcName == "poly":
#    polyInst = numpy.polynomial.Polynomial(poly)
#  else:
#    polyInst = numpy.polynomial.Chebyshev(poly)
#
#  yy = polyInst(xList)
#  ax.plot(xList,yy,"-b")
#  
#  fig.savefig("debugLinearChi2.png")

  return poly

def linearChi2TH1(hist,order,funcName="poly"):
  nBinsX = hist.GetNbinsX()
  x = numpy.zeros(nBinsX)
  y = numpy.zeros(nBinsX)
  err = numpy.zeros(nBinsX)
  for i in range(1,nBinsX+1):
    x[i-1] = hist.GetXaxis().GetBinCenter(i)
    y[i-1] = hist.GetBinContent(i)
    err[i-1] = hist.GetBinError(i)
  print x
  print y
  result = linearChi2(x,y,order,yErrList=err,funcName=funcName)
  print result
  return result

def toyHistogram(hist):
  nBins = hist.GetNbinsX()
  random = root.TRandom3()
  for i in range(nBins+2):
    mean = hist.GetBinContent(i)
    n = random.Poisson(mean)
    err = sqrt(n)
    hist.SetBinContent(i,n)
    hist.SetBinError(i,err)

def getXbinsHighLow(hist,low,high):
  axis = hist.GetXaxis()
  xbinLow = axis.FindBin(low)
  xbinHigh = axis.FindBin(high)
  #print("xbinhigh: {0}, {1}, {2}".format(xbinHigh,axis.GetBinLowEdge(xbinHigh),float(high)))
  if axis.GetBinLowEdge(xbinHigh)==float(high):
    xbinHigh -= 1
  return xbinLow, xbinHigh

def getIntegralAll(hist,boundaries=[]):
  xbinLow = None
  xbinHigh = None
  if len(boundaries)==0:
    xbinLow = 0
    xbinHigh = hist.GetXaxis().GetNbins()
  elif len(boundaries)==2:
    xbinLow, xbinHigh = getXbinsHighLow(hist,boundaries[0],boundaries[1])
  else:
    return -1
  if hist.InheritsFrom("TH2"):
    nBinsY = hist.GetYaxis().GetNbins()
    return hist.Integral(xbinLow,xbinHigh,0,nBinsY+1)
  elif hist.InheritsFrom("TH1"):
    return hist.Integral(xbinLow,xbinHigh)
  else:
    return -1

def getIntegralLowHigh(hist,lowBoundaries,highBoundaries):
  lowInt = getIntegralAll(hist,lowBoundaries)
  highInt = getIntegralAll(hist,highBoundaries)
  return lowInt+highInt

def sqrtTH1(hist):
  nBins = hist.GetNbinsX()
  for i in range(nBins+2):
    n = hist.GetBinContent(i)
    nErr = hist.GetBinError(i)
    if n < 0.0:
      n = 0.0
    hist.SetBinContent(i,sqrt(n))
    hist.SetBinError(i,sqrt(nErr))

class CrossSecsErrs:
  def __init__(self,csvDict):
    self.data = csvDict
    self.nominal = {}
    self.err = {}
    self.errUp = {}
    self.errDown = {}
    self.lnN = {}
    for key in self.data:
      self.nominal[key] = self.data[key][0]
      self.errUp[key] = self.data[key][1]/100.
      self.errDown[key] = self.data[key][2]/100.
      self.err[key] = max(abs(self.data[key][1]),abs(self.data[key][2]))/100.
      self.lnN[key] = self.err[key] + 1.0
  def __getitem__(self,key):
    return self.extrap(self.nominal,key)
  def getLnN(self,key):
    return self.extrap(self.lnN,key)
  def has_key(self,key):
    return self.data.has_key(key)
  def keys(self):
    return self.data.keys()
  def extrap(self,dict,mass):
    if dict.has_key(mass):
      return dict[mass]
    massStr = mass
    mass = float(mass)
    dPos = []
    dNeg = []
    massKeys = dict.keys()
    for iMassStr in massKeys:
      iMass = float(iMassStr)
      if iMass >= mass:
        dPos.append(iMass-mass)
        dNeg.append(1e8)
      else:
        dPos.append(1e8)
        dNeg.append(mass-iMass)
    i1 = dPos.index(min(dPos))
    i2 = dNeg.index(min(dNeg))
    m1 = float(massKeys[i1])
    m2 = float(massKeys[i2])
    val1 = float(dict[massKeys[i1]])
    val2 = float(dict[massKeys[i2]])
    slope = (val2-val1)/(m2-m1)
    return val1 + slope*(mass-m1)

  

def readCSVXS(filename):
  f = open(filename)
  rd = csv.reader(f)
  result = {}
  for row in rd:
    if len(row) == 0:
        continue
    if len(row[0]) == 0:
        continue
    if re.search(r"[^\d.\s]",row[0]):
        continue
    mass = float(row[0])
    prec = "0"
    if mass % 1 > 0:
        prec = '1'
    result[("%."+prec+"f") % (mass)] = [float(i) for i in row[1:]]
  f.close()
  return CrossSecsErrs(result)

def getRooBinningFromTH1(hist):
  nbins = hist.GetNbinsX()
  xmax = hist.GetXaxis().GetXmax()
  xmin = hist.GetXaxis().GetXmin()
  print nbins,xmax,xmin
  return root.RooFit.RooBinning(nbins,xmin,xmax)

def getBinningFromTH1(hist,newName):
  nbins = hist.GetNbinsX()
  xmax = hist.GetXaxis().GetXmax()
  xmin = hist.GetXaxis().GetXmin()
  print nbins,xmax,xmin
  return newName,newName,nbins,xmin,xmax

def getRooVarRange(variable,name):
  assert(variable.InheritsFrom("RooRealVar"))
  assert(type(name)==str)
  binning = variable.getBinning(name)
  return binning.lowBound(), binning.highBound()

def saveAs(canvas,name):
  canvas.SaveAs(name+".png")
  canvas.SaveAs(name+".pdf")
  canvas.SaveAs(name+".eps")
  canvas.SaveAs(name+".root")
  canvas.SaveAs(name+".C")

def setLegPos(leg,legPos):
  leg.SetX1NDC(legPos[0])
  leg.SetX2NDC(legPos[2])
  leg.SetY1NDC(legPos[1])
  leg.SetY2NDC(legPos[3])

def getBinWidthStr(hist):
    binWidth = (hist.GetXaxis().GetXmax()-hist.GetXaxis().GetXmin())/hist.GetXaxis().GetNbins()
    binWidthPrec = "0"
    if binWidth % 1 > 0.0:
      binWidthPrec = "1"
      if binWidth*10 % 1 > 0.0:
        binWidthPrec = "2"
    return ("%."+binWidthPrec+"f") % (binWidth)

def getEfficiencyInterval(passed,total):
  eff = root.TEfficiency()
  nom = 0.
  quant=0.682689492137
  if total>0:
    nom = float(passed)/total
  low = eff.ClopperPearson(int(total),int(passed),quant,False)
  high = eff.ClopperPearson(int(total),int(passed),quant,True)
  return [low,nom,high]

class EfficiencyReader:
  def __init__(self,fileDir="EffMassScan/"):
    self.data = {}
    self.fileDir = fileDir
    for ifn in glob.glob(fileDir+"*.txt"):
      fmatch = re.match(r".*/Eff(.+)_(.+)Higgs([0-9.]+).txt",ifn)
      if not fmatch:
        print("Warning: EfficiencyReader: filename: %s isn't recognized" % (ifn))
        continue
      energy = fmatch.group(1)
      prodMode = fmatch.group(2)
      mass = fmatch.group(3)
      if not self.data.has_key(energy):
        self.data[energy] = {}
      if not self.data[energy].has_key(prodMode):
        self.data[energy][prodMode] = {}
      f = open(ifn)
      for iline in f:
        lmatch = re.match(r"([\w]+)[\s]+([\d.]+)[\s]+([\d.]+)",iline)
        if not fmatch:
          print("Warning: EfficiencyReader: text line isn't recognized: \n%s\n in file: %s" % (iline,ifn))
          continue
        category = lmatch.group(1)
        efficiency = float(lmatch.group(2))
        efficiencyError = float(lmatch.group(3))
        if not self.data[energy][prodMode].has_key(category):
          self.data[energy][prodMode][category] = {}
        if not self.data[energy][prodMode][category].has_key(mass):
          self.data[energy][prodMode][category][mass] = {}
        self.data[energy][prodMode][category][mass]['eff'] = efficiency
        self.data[energy][prodMode][category][mass]['effErr'] = efficiencyError
        
  def __getitem__(self,key):
    return self.data[key]
  def __call__(self,energy,prodMode,category,mass):
    eff = self.data[energy][prodMode][category][mass]['eff']
    err = self.data[energy][prodMode][category][mass]['effErr']
    return eff, err
  def getEfficiency(self,energy,prodMode,category,mass):
    return self.data[energy][prodMode][category][mass]['eff']
  def getEfficiencyError(self,energy,prodMode,category,mass):
    return self.data[energy][prodMode][category][mass]['effErr']
  def __str__(self):
    result = ""
    for energy in self.data:
     for mode in self.data[energy]:
      result += "%s %s:\n" % (mode,energy)
      sortedCats = sorted(self.data[energy][mode].keys())
      for cat in sortedCats:
        result += "  %s:\n" % (cat)
        sortedMasses = sorted(self.data[energy][mode][cat].keys())
        for mass in sortedMasses:
          eff = self.data[energy][mode][cat][mass]['eff']
          err = self.data[energy][mode][cat][mass]['effErr']
          result += "    %5s:  %8.2f%%  +/-  %8.2f%%\n" % (mass,eff*100.,err*100.)
    return result
  def __repr__(self):
    return str(self)
  def plot(self,folderBase):
    canvas = root.TCanvas("canvas")
    for energy in self.data:
     for mode in self.data[energy]:
      sortedCats = sorted(self.data[energy][mode].keys())
      for cat in sortedCats:
        sortedMasses = sorted(self.data[energy][mode][cat].keys())
        fileNameOut = folderBase+"/png/eff_"+mode+energy+"_"+cat
        graph = root.TGraphErrors()
        iPoint = 0
        eff = None
        err = None
        ymax = 0.0
        ymin = 1.0
        for mass in sortedMasses:
          x = float(mass)
          eff = self.data[energy][mode][cat][mass]['eff']
          err = self.data[energy][mode][cat][mass]['effErr']
          graph.SetPoint(iPoint,x,eff)
          graph.SetPointError(iPoint,0.,err)
          iPoint += 1
          if ymax < eff+err:
            ymax = eff+err
          if ymin > eff+err:
            ymin = eff+err
        setHistTitles(graph,"m_{H} [GeV/c^{2}]","Efficiency #times Acceptance")
        graph.SetTitle("%s %s %s:\n" % (mode,energy,cat))
        fitListPol=[0.1,0.,0.]
        #polim = root.TF1("polim","pol2",fitListPol[0],fitListPol[1], fitListPol[2])
        polim = root.TF1("polim","pol2")
        polim.SetLineColor(root.kBlack)
        self.polim = polim
        #graph.Fit(polim,"WLMEIQ","",fitListPol[0],fitListPol[1], fitListPol[2])
        graph.Fit(polim,"LE")

        graph.Draw("ape")
        ywidth = 1.0 * (ymax-ymin)
        ymax = ywidth+ymax
        ymin = ymin-ywidth
        if ymax > 1.0:
            ymax = 1.0
        ymin = 0.0
        graph.GetYaxis().SetRangeUser(ymin,ymax)
        graph.Draw("ape")
        canvas.SaveAs(fileNameOut+".png")
        canvas.SaveAs(fileNameOut+".pdf")
        canvas.SaveAs(fileNameOut+".root")

        # save the output of the fit
        outputfile = folderBase + '/effExtrapolation_' + mode + '_' + energy + '_' +  cat + '.txt' 
        outfile = open(outputfile,'w')

        outfile.write('#par0 err_par0 par1 err_par1 par2 err_par2 chisquare ndf\n')  
        outfile.write('{0:.10g} {1:.10g} {2:.10g} {3:.10g} {4:.10g} {5:.10g} {6:.4g} {7}\n'.format(polim.GetParameter(0),
                                                                                                   polim.GetParError(0),
                                                                                                   polim.GetParameter(1),
                                                                                                   polim.GetParError(1),
                                                                                                   polim.GetParameter(2),
                                                                                                   polim.GetParError(2),
                                                                                                   polim.GetChisquare(),
                                                                                                   polim.GetNDF()
                                                                                                   )
                      ) 
        
        outfile.close()
        
        

def fitDGFindQuantiles(hist,level):
    if not hist.InheritsFrom("TH1"):
      print("Error: fitDGFindQuantiles: input not TH1*. Exiting.")
      sys.exit(1)
    channelName = "silly"
    name = "silly"
    minMass = 110
    maxMass = 170
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    rooDataset = root.RooDataHist("DGDataSet","DGDataSet",root.RooArgList(mMuMu),hist)
    
    debug = ""
    debug += "### makePDFSigDG: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,mMuMu.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    rooParamList = []

    meanG1 = root.RooRealVar(channelName+"_"+name+"_MeanG1",
                             channelName+"_"+name+"_MeanG1", 
                             123.,100.,150.)
    meanG2 = root.RooRealVar(channelName+"_"+name+"_MeanG2",
                             channelName+"_"+name+"_MeanG2", 
                             125.,100.,150.)
    
    widthG1 = root.RooRealVar(channelName+"_"+name+"_WidthG1",
                             channelName+"_"+name+"_WidthG1", 
                             5.,2.,10.)
    widthG2 = root.RooRealVar(channelName+"_"+name+"_WidthG2",
                              channelName+"_"+name+"_WidthG2", 
                              1.,0.5,5.)
    
    mixGG = root.RooRealVar(channelName+"_"+name+"_mixGG",
                            channelName+"_"+name+"_mixGG", 
                            0.5,0.,1.)
    gaus1 = root.RooGaussian(channelName+"_"+name+"_gaus1",
                             channelName+"_"+name+"_gaus1",
                             mMuMu,meanG1,widthG1)
    gaus2 = root.RooGaussian(channelName+"_"+name+"_gaus2",
                             channelName+"_"+name+"_gaus2",
                             mMuMu,meanG2,widthG2)
    pdfMmumu = root.RooAddPdf(name,
                              name,
                              gaus1,gaus2,mixGG)
    #workspaceImportFn(pdfMmumu)
    rooParamList += [meanG1,meanG2,widthG1,widthG2,mixGG]
    
    PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
    fr.SetName(name+"_fitResult")

    pdftf = pdfMmumu.asTF(root.RooArgList(mMuMu))


    # Debugging
    #canvas = root.TCanvas("sillyCanvas")
    #pdftf.Draw()
    #now = datetime.datetime.now().isoformat()
    #canvas.SaveAs("debug_"+now+".png")

    quants =  getMedianAndQuantileInterval(pdftf,level)
    return quants

class RooCompareModels:
  def __init__(self,xVar,data,pdfList,frList,pdfTitleList,title,energyStr,lumi,canvas=None):
    self.xVar = xVar
    self.pdfList = pdfList
    self.data = data
    self.frList = frList
    self.pdfTitleList = pdfTitleList
    self.title = title
    self.energyStr = energyStr
    self.lumi = lumi

    self.lumiStr = "L = {0:.1f} fb^{{-1}}".format(lumi)

    self.xtitle = xVar.GetTitle()

    binning = xVar.getBinning()
    self.binning = binning
    nBins = self.binning.numBins()
    binWidth = (self.binning.highBound()-self.binning.lowBound())/nBins
    self.binWidth = binWidth
    xVar.setRange("RMPRange",self.binning.lowBound(),self.binning.highBound())
    rangeArg = root.RooFit.Range("RMPRange")
    self.rangeArg = rangeArg

    nowStr = str(int(time.time()*1e6))
    self.nowStr = nowStr
    if canvas == None:
      canvas = root.TCanvas("canvas"+nowStr)
    self.canvas = canvas

    self.colors = [root.kBlue,root.kRed,root.kGreen,root.kCyan,root.kMagenta,root.kOrange-3,root.kViolet-6]
    
  def draw(self,saveName):
    canvas = self.canvas
    data = self.data
    xVar = self.xVar
    xtitle = self.xtitle

    rangeArg = self.rangeArg
    binningArg = root.RooFit.Binning("")
    lineWidthArg = root.RooFit.LineWidth(2)
    lineDrawOptArg = root.RooFit.DrawOption("L")
    graphDrawOptArg = root.RooFit.DrawOption("PEZ")


    frame       = xVar.frame(root.RooFit.Title(""))
    data.plotOn(frame,      graphDrawOptArg,binningArg)

    leg = root.TLegend(0.55,0.60,0.9,0.9)
    leg.SetFillColor(0)
    leg.SetLineColor(0)

    fakeGraphs = []
    for iPdf, pdf in enumerate(self.pdfList):
      fr = self.frList[iPdf]
      pdfTitle = self.pdfTitleList[iPdf]
      color = self.colors[iPdf]

      #Set the PDF pars value from the FitResults
      setPDFfromFR(fr,pdf,data)

      lineColorArg = root.RooFit.LineColor(color)
      pdf.plotOn(frame,lineDrawOptArg,lineColorArg,lineWidthArg,rangeArg)

      fakeG = root.TGraph()
      fakeG.SetLineColor(color)
      fakeG.SetLineWidth(2)
      leg.AddEntry(fakeG,pdfTitle,"l")
      fakeGraphs.append(fakeG)

    frame.SetTitle("")
    frame.Draw()

    unitMatch =  re.search(r"GeV([\s]*/[\s]*c\^\{2\}|[\s]*/[\s]*c)?",xtitle)
    units = ""
    if unitMatch:
      units = " "+unitMatch.group(0)
    frame.SetXTitle(xtitle)
    frame.SetYTitle("Events/"+str(self.binWidth)+units)
    leg.Draw()
    energyLumiStr = "#sqrt{{s}} = {0}, L = {1:.1f} fb^{{-1}}".format(self.energyStr.replace("TeV"," TeV"),self.lumi)
    drawStandardCaptions(canvas,self.title,energyLumiStr,preliminaryString="CMS Internal")

    saveAs(canvas,saveName)

  def drawCurveHists(self,saveName):
    canvas = self.canvas
    hists = []
    leg = root.TLegend(0.55,0.60,0.9,0.9)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    xtitle = self.xtitle
    unitMatch =  re.search(r"GeV([\s]*/[\s]*c\^\{2\}|[\s]*/[\s]*c)?",xtitle)
    units = ""
    if unitMatch:
      units = " "+unitMatch.group(0)
    for i,model in enumerate(self.pdfList):
      hist = self.getHistFromModel(model)
      setHistTitles(hist,xtitle,"Events/"+str(self.binWidth)+units)
      color = self.colors[i]
      pdfTitle = self.pdfTitleList[i]
      hist.SetLineColor(color)
      leg.AddEntry(hist,pdfTitle,"l")
      if i == 0:
        hist.Draw()
      else:
        hist.Draw("SAME")
      hists.append(hist)
    leg.Draw()
    energyLumiStr = "#sqrt{{s}} = {0}, L = {1:.1f} fb^{{-1}}".format(self.energyStr.replace("TeV"," TeV"),self.lumi)
    drawStandardCaptions(canvas,self.title,energyLumiStr,preliminaryString="CMS Internal")
    canvas.RedrawAxis()
    saveAs(canvas,saveName)

  def drawPullHists(self,saveName):
    canvas = self.canvas
    hists = []
    leg = root.TLegend(0.55,0.60,0.9,0.9)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    xtitle = self.xtitle
    unitMatch =  re.search(r"GeV([\s]*/[\s]*c\^\{2\}|[\s]*/[\s]*c)?",xtitle)
    units = ""
    if unitMatch:
      units = " "+unitMatch.group(0)
    axisHist = root.TH2F("axisHist","",1,110,160,1,-3,8)
    setHistTitles(axisHist,xtitle,"#frac{Data-Fit}{#sqrt{Fit}}")
    zeroGraph = root.TGraph()
    zeroGraph.SetLineStyle(3)
    zeroGraph.SetPoint(0,110,0)
    zeroGraph.SetPoint(1,160,0)
    axisHist.Draw()
    zeroGraph.Draw("L")
    for i,model in enumerate(self.pdfList):
      hist = self.getPullHistFromModel(model)
      setHistTitles(hist,xtitle,"#frac{Data-Fit}{#sqrt{Fit}}")
      color = self.colors[i]
      pdfTitle = self.pdfTitleList[i]
      hist.SetLineColor(color)
      leg.AddEntry(hist,pdfTitle,"l")
      hists.append(hist)
    for hist in reversed(hists):
      hist.Draw("same")
    leg.Draw()
    energyLumiStr = "#sqrt{{s}} = {0}, L = {1:.1f} fb^{{-1}}".format(self.energyStr.replace("TeV"," TeV"),self.lumi)
    drawStandardCaptions(canvas,self.title,energyLumiStr,preliminaryString="CMS Internal")
    canvas.RedrawAxis()
    saveAs(canvas,saveName)
        
  def drawDiff(self,iModel,saveName):
    canvas = self.canvas
    hists = []
    leg = root.TLegend(0.55,0.60,0.9,0.9)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    xtitle = self.xtitle
    unitMatch =  re.search(r"GeV([\s]*/[\s]*c\^\{2\}|[\s]*/[\s]*c)?",xtitle)
    units = ""
    if unitMatch:
      units = " "+unitMatch.group(0)
    axisHist = root.TH2F("axisHist","",1,110,160,1,-75,100)
    #axisHist = root.TH2F("axisHist","",1,110,160,1,-15,20)
    #axisHist = root.TH2F("axisHist","",1,110,160,1,-2,6)
    setHistTitles(axisHist,xtitle,"PDF - "+self.pdfList[iModel].GetName()+" [Events/"+str(self.binWidth)+units+"]")
    zeroGraph = root.TGraph()
    zeroGraph.SetLineStyle(3)
    zeroGraph.SetPoint(0,110,0)
    zeroGraph.SetPoint(1,160,0)
    for i,model in enumerate(self.pdfList):
      hist = self.getHistFromModel(model)
      setHistTitles(hist,xtitle,"PDF - "+self.pdfList[iModel].GetName()+" [Events/"+str(self.binWidth)+units+"]")
      color = self.colors[i]
      pdfTitle = self.pdfTitleList[i]
      hist.SetLineColor(color)
      hists.append(hist)
      if i != iModel:
        leg.AddEntry(hist,pdfTitle,"l")
    iModelHist = hists[iModel]
    drawn = False
    axisHist.Draw()
    zeroGraph.Draw("L")
    drawn = True
    for i,hist in enumerate(hists):
      if i == iModel:
        continue
      hist.Add(iModelHist,-1)
      if drawn:
        hist.Draw("SAME")
      else:
        hist.Draw()
        drawn = True
    leg.Draw()
    energyLumiStr = "#sqrt{{s}} = {0}, L = {1:.1f} fb^{{-1}}".format(self.energyStr.replace("TeV"," TeV"),self.lumi)
    drawStandardCaptions(canvas,self.title,energyLumiStr,preliminaryString="CMS Internal")
    canvas.RedrawAxis()
    saveAs(canvas,saveName)
        
  def getHistFromModel(self,model):
    nowStr = str(int(time.time()*1e6))
    frame = self.xVar.frame(root.RooFit.Name("frameToGetHist"+nowStr),self.rangeArg)
    histPlotName = self.data.GetName()+nowStr
    pdfPlotName = model.GetName()+nowStr
    self.data.plotOn(frame,root.RooFit.Name(histPlotName))
    model.plotOn(frame,root.RooFit.Name(pdfPlotName),self.rangeArg)

    hist = frame.findObject(histPlotName)
    curve = frame.findObject(pdfPlotName)
    assert(hist)
    assert(curve)

    nBins = self.binning.numBins()
    lowB  = self.binning.lowBound() 
    highB = self.binning.highBound() 

    myCurveHist =  root.TH1F("myCurveHist_"+model.GetName()+"_"+nowStr,"",
                          int(nBins), lowB, highB
                          )
      
    x = root.Double(0.)
    y = root.Double(0.)

    curve.GetPoint(0,x,y) # Get Curve Start X
    xCurveMin = float(x)
    curve.GetPoint(curve.GetN()-1,x,y) # Get Curve End X
    xCurveMax = float(x)
    iBin = 1
    for i in range(1,self.binning.numBins()+1):
      hist.GetPoint(i-1,x,y)
      if (float(x) < lowB or float(x) > highB):
        continue
      if x > xCurveMin and x < xCurveMax:
        curvePoint = curve.interpolate(x)
        myCurveHist.SetBinContent(iBin,curvePoint)
      iBin += 1

    return myCurveHist

  def getPullHistFromModel(self,model):
    nowStr = str(int(time.time()*1e6))
    frame = self.xVar.frame(root.RooFit.Name("frameToGetHist"+nowStr),self.rangeArg)
    histPlotName = self.data.GetName()+nowStr
    pdfPlotName = model.GetName()+nowStr
    self.data.plotOn(frame,root.RooFit.Name(histPlotName))
    model.plotOn(frame,root.RooFit.Name(pdfPlotName),self.rangeArg)

    hist = frame.findObject(histPlotName)
    curve = frame.findObject(pdfPlotName)
    assert(hist)
    assert(curve)

    nBins = self.binning.numBins()
    lowB  = self.binning.lowBound() 
    highB = self.binning.highBound() 

    myPullHist =  root.TH1F("myPullHist_"+model.GetName()+"_"+nowStr,"",
                          int(nBins), lowB, highB
                          )
      
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
      if x > xCurveMin and x < xCurveMax:
        curvePoint = curve.interpolate(x)
        if curvePoint == 0.:
          pull = 0.
        else:
          pull -= curvePoint
          pull /= sqrt(curvePoint)
      else:
        pull = 0.
      myPullHist.SetBinContent(iBin,pull)
      iBin += 1

    return myPullHist

  def printGOFTable(self,log=sys.stdout):
    log.write("% GOF Comparison: "+self.title+" "+self.energyStr+"\n")
    log.write(r"\begin{tabular}{|l|r|r|} \hline"+"\n")
    log.write(r"\multicolumn{3}{|c|}{\textbf{"+self.title+" "+self.energyStr.replace("TeV"," TeV")+r"}} \\ \hline"+"\n")
    log.write(r"PDF & $\chi^2/$NDF & $\chi^2$ Probability \\ \hline \hline"+"\n")
    dataHist = self.data.binnedClone()
    for i,model in enumerate(self.pdfList):
      chi2,ndf,chi2Prob = rooCalcChi2(model,dataHist)
      pdfTitle = self.pdfTitleList[i]
      if "#" in pdfTitle or "^" in pdfTitle or "_" in pdfTitle:
        pdfTitle = r"$\mathrm{"+pdfTitle+"}$"
        pdfTitle = pdfTitle.replace("#",'\\')
      log.write(r"{0} & {1:.1f}/{2} & {3:.3g} \\ \hline".format(pdfTitle,chi2,ndf,chi2Prob)+"\n")
    log.write(r"\end{tabular}"+"\n\n")

  def printDiffTable(self,iModel,log=sys.stdout):
    hists = []
    for i,model in enumerate(self.pdfList):
      hist = self.getHistFromModel(model)
      hists.append(hist)
    iModelHist = hists[iModel]
    nBins = iModelHist.GetNbinsX()
    ndf = nBins
    iModelTitle = self.pdfTitleList[iModel]
    if "#" in iModelTitle or "^" in iModelTitle or "_" in iModelTitle:
      iModelTitle = r"$\bm{\mathrm{"+iModelTitle+"}}$"
      iModelTitle = iModelTitle.replace("#",'\\')
    log.write("% Difference Between PDFs: "+self.title+" "+self.energyStr+"\n")
    log.write(r"\begin{tabular}{|l|r|r|} \hline"+"\n")
    log.write(r"\multicolumn{3}{|c|}{\textbf{PDF Fit Difference with "+iModelTitle+r"}} \\"+"\n")
    log.write(r"\multicolumn{3}{|c|}{\textbf{"+self.title+" "+self.energyStr.replace("TeV"," TeV")+r"}} \\ \hline"+"\n")
    #log.write(r"PDF & $\chi^2/$NDF & $\chi^2$ Probability \\ \hline \hline"+"\n")
    log.write(r"PDF & $\chi^2$ & $N_{bins}$ \\ \hline \hline"+"\n")
    for i,hist in enumerate(hists):
      if i == iModel:
        continue
      chi2 = 0.
      for iBin in range(1,nBins+1):
        nEvt = hist.GetBinContent(iBin)
        nEvtIModel = iModelHist.GetBinContent(iBin)
        chi2 += (nEvt-nEvtIModel)**2/nEvtIModel
      pdfTitle = self.pdfTitleList[i]
      if "#" in pdfTitle or "^" in pdfTitle or "_" in pdfTitle:
        pdfTitle = r"$\mathrm{"+pdfTitle+"}$"
        pdfTitle = pdfTitle.replace("#",'\\')
      #chi2Prob = scipy.stats.chi2.sf(chi2,ndf)
      #log.write(r"{0} & {1:.1f}/{2} & {3:.3g} \\ \hline".format(pdfTitle,chi2,ndf,chi2Prob)+"\n")
      log.write(r"{0} & {1:.2f} & {2} \\ \hline".format(pdfTitle,chi2,ndf)+"\n")
    log.write(r"\end{tabular}"+"\n\n")

class RooModelPlotter:
  def __init__(self,xVar,pdf,data,fr,title,energyStr,lumi,
                backgroundPDFName=None,signalPDFName=None,pdfDotLineName=None,
                nSignal=0,signalPdf=None,
                RangeName="",
                canvas=None,
                caption1="",caption2="",caption3="",caption4="",
                legEntryData="Data",legEntryModel="Background Model",legEntrySignal="Signal",
                pullsYLabel="#frac{Data-Fit}{#sqrt{Fit}}",
                extraPDFs=[],
                extraLegEntries=[],
                extraPDFDotLineNames=[],
                doLinearErrs=True
              ):
    self.xVar = xVar
    self.pdf = pdf
    self.pdfName = pdf.GetName()
    self.data = data
    self.fr = fr
    self.title = title
    self.energyStr = energyStr
    self.lumi = lumi
    self.backgroundPDFName = backgroundPDFName
    self.signalPDFName = signalPDFName
    self.pdfDotLineName = pdfDotLineName
    self.nSignal = nSignal
    self.pullsYLabel = pullsYLabel
    nowStr = str(int(time.time()*1e6))
    self.nowStr = nowStr
    self.caption1 = caption1
    self.caption2 = caption2
    self.caption3 = caption3
    self.caption4 = caption4
    self.myCurveHist = None
    self.myDataHist = None
    assert(len(extraPDFs) == len(extraLegEntries))
    self.extraPDFs = extraPDFs
    #self.extraPDFColors = [root.kGreen+3,root.kMagenta+3,root.kOrange-3]
    self.extraPDFColors = [root.kSpring-7,root.kMagenta+3,root.kOrange-3]
    assert(len(self.extraPDFs)<=len(self.extraPDFColors))
    self.extraLegEntries = extraLegEntries
    self.extraPDFDotLineNames = extraPDFDotLineNames
    self.doLinearErrs = doLinearErrs

    self.legEntryData = legEntryData
    self.legEntryModel = legEntryModel

    self.lumiStr = "L = {0:.1f} fb^{{-1}}".format(lumi)

    xtitle = xVar.GetTitle()

    binning = xVar.getBinning()
    self.binning = binning
    nBins = self.binning.numBins()
    binWidth = (self.binning.highBound()-self.binning.lowBound())/nBins

    if canvas == None:
      canvas = root.TCanvas("canvas"+nowStr)
    self.canvas = canvas

    #Set the PDF pars value from the FitResults
    setPDFfromFR(self.fr,self.pdf,self.data)
    #setAddPDFfromFR(self.fr,self.pdf,self.data)

    xtitle = xVar.GetTitle()

    binning = xVar.getBinning()
    self.binning = binning
    nBins = self.binning.numBins()
    binWidth = (self.binning.highBound()-self.binning.lowBound())/nBins

    if canvas == None:
      canvas = root.TCanvas("canvas"+nowStr)
    self.canvas = canvas
    self.canvas2 = root.TCanvas("canvas2"+nowStr,"",
                        2*root.gStyle.GetCanvasDefW(),
                        root.gStyle.GetCanvasDefH()
                        )
    self.canvas2.SetMargin(0,0,0,0)
    self.canvas2.Divide(2,1,0,0,0)

    self.tlatex = root.TLatex()
    self.tlatex.SetNDC()
    self.tlatex.SetTextFont(root.gStyle.GetLabelFont())
    self.tlatex.SetTextSize(root.gStyle.GetLabelSize())
    self.tlatex.SetTextAlign(22)

    errVisArg = root.RooFit.VisualizeError(fr,1,self.doLinearErrs)
    errColorArg = root.RooFit.FillColor(root.kCyan)
    lineColorArg = root.RooFit.LineColor(root.kBlue)
    lineWidthArg = root.RooFit.LineWidth(2)
    lineDrawOptArg = root.RooFit.DrawOption("L")
    graphDrawOptArg = root.RooFit.DrawOption("PEZ")
    sigLineColorArg = root.RooFit.LineColor(root.kRed)
    dotLineStyleArg = root.RooFit.LineStyle(2)

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
    # simply cloning the frame did not work for me... :-(
    # working around it
    frameBkgSub = xVar.frame(root.RooFit.Title("BkgSub"))

    if (RangeName):
      self.rangename = RangeName
      frame = xVar.frame(root.RooFit.Title(""),root.RooFit.Range(self.rangename))
      frameBkgSub = xVar.frame(root.RooFit.Title("BkgSub"),root.RooFit.Range(self.rangename))
      #frame.SetMaximum(6000)

    self.frame       = frame
    self.frameBkgSub = frameBkgSub
    
    data.plotOn(frame,      graphDrawOptArg,binningArg)
      
    #print "backgroundPDFName = %s" % backgroundPDFName 
    if backgroundPDFName != None:
      bakCompArg = root.RooFit.Components(backgroundPDFName)
      pdf.plotOn(frame,errVisArg,errColorArg,bakCompArg,rangeArg,tmpBakPDFErrorNameArg)
      for iColor,extraPdf in enumerate(self.extraPDFs):
        extraPdf.plotOn(frame,lineDrawOptArg,root.RooFit.LineColor(self.extraPDFColors[iColor]),lineWidthArg,rangeArg)
      pdf.plotOn(frame,lineDrawOptArg,lineColorArg,lineWidthArg,bakCompArg,rangeArg,tmpBakPDFNameArg)
    else:
      pdf.plotOn(frame,errVisArg,errColorArg,rangeArg,tmpBakPDFErrorNameArg)
      for iColor,extraPdf in enumerate(self.extraPDFs):
        extraPdf.plotOn(frame,lineDrawOptArg,root.RooFit.LineColor(self.extraPDFColors[iColor]),lineWidthArg,rangeArg)
        if len(self.extraPDFDotLineNames) == len(self.extraPDFs):
          dotLineCompArg = root.RooFit.Components(self.extraPDFDotLineNames[iColor])
          extraPdf.plotOn(frame,lineDrawOptArg,root.RooFit.LineColor(self.extraPDFColors[iColor]),lineWidthArg,dotLineStyleArg,dotLineCompArg,rangeArg)
      if self.pdfDotLineName != None:
        dotLineCompArg = root.RooFit.Components(self.pdfDotLineName)
        pdf.plotOn(frame,lineDrawOptArg,lineColorArg,lineWidthArg,dotLineStyleArg,dotLineCompArg,rangeArg)
      pdf.plotOn(frame,lineDrawOptArg,lineColorArg,lineWidthArg,rangeArg,tmpBakPDFNameArg)

    data.plotOn(frame,graphDrawOptArg,binningArg,tmpDataHistNameArg)

    # filling it with white points for the sig PDF normalization
    data.plotOn(frameBkgSub,graphDrawOptArg,
                root.RooFit.MarkerColor(root.kWhite),
                root.RooFit.LineColor(root.kWhite),
                binningArg,tmpDataHistNameArg)

    frame.SetTitle("")
    frame.GetXaxis().SetLabelSize(0)
    frame.GetYaxis().SetLabelSize(0.050)
    frame.GetYaxis().SetTitleSize(0.055*1.2)
    frame.GetYaxis().SetTitleOffset(
        0.85*frame.GetYaxis().GetTitleOffset()
        )

    frameBkgSub.SetTitle("")
    frameBkgSub.GetYaxis().SetLabelSize(0.04)
    frameBkgSub.GetYaxis().SetTitleSize(0.055*0.9)
    frameBkgSub.GetYaxis().SetTitleOffset(
        1.1*frameBkgSub.GetYaxis().GetTitleOffset()
        )

    
    unitMatch =  re.search(r"GeV([\s]*/[\s]*c\^\{2\}|[\s]*/[\s]*c)?",xtitle)
    units = ""
    if unitMatch:
      units = " "+unitMatch.group(0)
    frame.SetYTitle("Events/"+str(binWidth)+units)
    frameBkgSub.SetYTitle("Events/"+str(binWidth)+units)
    #if (RangeName):
      #self.frame.SetMaximum(6000)
      #self.frame.SetMinimum(1000)
      #self.frame.GetYaxis().SetRangeUser(1500,8000)
    #else:
    #  self.frame.SetMinimum(0)

    dataHist = None
    if type(self.data) == root.RooDataHist:
      dataHist = self.data
    else:
      dataHist = self.data.binnedClone()
    self.chi2 = rooCalcChi2(self.pdf,dataHist)

    # Pulls Frame
    pullsHist = self.makePullPlotHist(frame,tmpDataHistName,tmpBakPDFName)
    self.pullsHist = pullsHist
    pullsHist.SetLineColor(root.kBlue)
    pullsHist.SetLineWidth(2)
    pullsHist.SetFillColor(856)
    pullsHist.SetFillStyle(1001)
    setHistTitles(pullsHist,xtitle,self.pullsYLabel)

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

    # Bkg Sub Hist
    bkgSubHist = self.makeBkgSubHist(frame,tmpDataHistName,tmpBakPDFName)
    bkgSubHist.SetLineColor(root.kGray+2)
    bkgSubHist.SetLineWidth(2)
    bkgSubHist.SetFillColor(root.kGray)
    bkgSubHist.SetFillStyle(1001)
    setHistTitles(bkgSubHist,xtitle,"Data-Bkg Fit Model")
    self.bkgSubHist = bkgSubHist

    bkgSubDataHist = root.RooDataHist("data_men_fit","",
                                      root.RooArgList(xVar),bkgSubHist)    

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
        sigPdfToDraw.plotOn(frameBkgSub,lineDrawOptArg,sigLineColorArg,lineWidthArg,componentToDraw)
      
        
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
    self.phonyExtraPDFLegHists = []
    for iColor in range(len(self.extraPDFs)):
        tmpHist = root.TH1F("phonyExtraPDFHist"+str(iColor)+nowStr,"",1,0,1)
        tmpHist.SetLineColor(self.extraPDFColors[iColor])
        tmpHist.SetLineWidth(2)
        self.phonyExtraPDFLegHists.append(tmpHist)
    
    legPos = [0.65,0.65,1.0-gStyle.GetPadRightMargin()-0.01,1.0-gStyle.GetPadTopMargin()-0.01]
    #legPos = [0.73,0.65,1.0-gStyle.GetPadRightMargin()-0.01,1.0-gStyle.GetPadTopMargin()-0.01]
    self.legPos = legPos
    self.leg = root.TLegend(*legPos)
    self.leg.SetFillColor(0)
    self.leg.SetLineColor(0)
    self.leg.AddEntry(self.phonyDatLegHist,self.legEntryData,"lp")
    self.leg.AddEntry(self.phonyFitLegHist,self.legEntryModel,"lf")
    for hist, legEntry in zip(self.phonyExtraPDFLegHists,self.extraLegEntries):
      self.leg.AddEntry(hist,legEntry,"l")

    legPosBkgSub = [0.64,0.72,0.92,0.88]
    self.legBkgSub = root.TLegend(*legPosBkgSub)
    self.legBkgSub.SetFillColor(0)
    self.legBkgSub.SetLineColor(0)
    self.legBkgSub.AddEntry(self.bkgSubHist,"Data - Background Model","lf")

    if signalPdf != None:
      if signalLegEntry != None:
        self.leg.AddEntry(self.phonySigLegHist,signalLegEntry,"l")
        self.legBkgSub.AddEntry(self.phonySigLegHist,signalLegEntry,"l")
      else:
        self.leg.AddEntry(self.phonySigLegHist,"Signal","l")
        self.legBkgSub.AddEntry(self.phonySigLegHist,"Signal","l")


  def draw(self,filenameNoExt,canvas=None,motherPad=None):
    if canvas==None:
        canvas=self.canvas
    if motherPad==None:
        motherPad=self.canvas
    nowStr = self.nowStr
    motherPad.SetLogy(0)
    motherPad.cd()
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
    motherPadToPad1FontScalingFactor = float(motherPad.YtoPixel(motherPad.GetY1()))/pad1.YtoPixel(pad1.GetY1())
    motherPadToPad2FontScalingFactor = float(motherPad.YtoPixel(motherPad.GetY1()))/pad2.YtoPixel(pad2.GetY1())
  
    # Main Pad
    pad1.cd();
    self.frame.Draw()
    self.addPDFNormError(pad1)
    self.leg.Draw()

    # Pulls Pad
    pad2.cd();
    self.pullsHist.Draw("")
  
    # Text
    self.pad1.cd()
    #self.tlatex.SetTextSize(root.gStyle.GetLabelSize())
    self.tlatex.SetTextSize(0.04*motherPadToPad1FontScalingFactor)
    self.tlatex.SetTextAlign(12)
    self.tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    self.tlatex.SetTextAlign(32)
    #self.tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.title)
    caption1 = self.caption1
    if caption1 != "":
        caption1 += ": "
    self.tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption1+self.title)

    self.tlatex.SetTextAlign(32)
    if (self.lumi != 0):
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.820,self.lumiStr)
    if (self.lumi == 0):
      self.tlatex.SetTextSize(0.04)
      self.tlatex.DrawLatex(self.legPos[0]-0.03,0.850,"#sqrt{s}=7 TeV L =  5.0 fb^{-1} ")
      self.tlatex.DrawLatex(self.legPos[0]-0.03,0.770,"#sqrt{s}=8 TeV L = 19.8 fb^{-1}")
            
    energyStr = self.energyStr
    if re.search(r"[\d]TeV",energyStr):
      energyStr = energyStr.replace("TeV"," TeV")
    if (self.energyStr != ""):
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.875,"#sqrt{s} = "+self.energyStr)

    if self.caption2 != "":
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.76,self.caption2)

    if self.caption3 != "":
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.70,self.caption3)

    if self.caption4 != "":
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.64,self.caption4)

    pad2.cd()
    self.tlatex.SetTextSize(self.pullsHist.GetYaxis().GetLabelSize())
    self.tlatex.SetTextAlign(12)
    #self.tlatex.DrawLatex(0.18,0.41,"#chi^{{2}}/NDF: {0:.3g}".format(self.chi2))
    self.tlatex.DrawLatex(0.18,0.41,"#chi^{{2}}/NDF = {0:.1f}/{1} = {2:.3g}; Probability: {3:.3g}".format(self.chi2[0],self.chi2[1],self.chi2[0]/self.chi2[1],self.chi2[2]))

    if (filenameNoExt != ""):
      saveAs(canvas,filenameNoExt)

    #### Save a copy of the pull histogram
    ##pullHistFile = root.TFile(filenameNoExt+"_pullHist.root","RECREATE")
    ##self.pullsHist.Write("pullHist")
    ##if self.myDataHist != None:
    ##  self.myDataHist.Write("dataHist")
    ##if self.myCurveHist != None:
    ##  self.myCurveHist.Write("fitHist")
    ##pullHistFile.Close()

  def drawWithParams(self,filenameNoExt,paramsToPlot=None):
    rightPad = self.canvas2.cd(2)
    paramPave = root.TPaveText(0,0,1,1,"NDC")
    paramPave.SetFillColor(0)
    paramPave.SetLineColor(0)
    paramPave.AddText("Fitted Parameters")
    fpf = self.fr.floatParsFinal()
    for i in range(fpf.getSize()):
      param = fpf.at(i)
      if paramsToPlot != None:
        doContinue = True
        for j in paramsToPlot:
          match = re.search(j,param.GetName())
          if match:
            doContinue = False
            break
        if doContinue:
          continue
      paramName = param.GetTitle()
      paramPave.AddText("{0:8} {1:8.3g} +/- {2:5.3g}".format(paramName,param.getVal(),param.getError()))
    paramPave.Draw()
    leftPad = self.canvas2.cd(1)
    self.draw(filenameNoExt,canvas=self.canvas2,motherPad=leftPad)
    self.paramPave = paramPave

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
    #tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.title)
    caption1 = self.caption1
    if caption1 != "":
        caption1 += ": "
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption1+self.title)
    tlatex.DrawLatex(0.98-gStyle.GetPadRightMargin(),0.875,"#sqrt{s} = "+self.energyStr)
    tlatex.DrawLatex(0.98-gStyle.GetPadRightMargin(),0.825,self.lumiStr)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.875,"#chi^{{2}}/NDF = {0:.2g}".format(float(chi2)/ndf))
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.825,"#mu = {0:.2f} #pm {1:.2f}".format(mean,meanErr))
    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.775,"#sigma = {0:.2f} #pm {1:.2f}".format(sigma,sigmaErr))
    saveAs(self.canvas,filenameNoExt)


  def makeBkgSubHist(self,frame,histPlotName,pdfPlotName):

    """
    Makes bkg subtracted hist which is (data-fit)
    where fit is the average value of the PDF within the data histogram bin.
    """
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


    # Histograms to compute the bkg subtracted histogram
    bkgSubHistData = root.TH1F("bkgSubHistData_"+histPlotName+"_"+pdfPlotName,"",
                               int(nBins), lowB, highB
                               )

    bkgSubHistFit = root.TH1F("bkgSubHistFit_"+histPlotName+"_"+pdfPlotName,"",
                              int(nBins), lowB, highB
                              )


    x = root.Double(0.)
    y = root.Double(0.)

    curve.GetPoint(0,x,y) # Get Curve Start X
    xCurveMin = float(x)
    curve.GetPoint(curve.GetN()-1,x,y) # Get Curve End X
    xCurveMax = float(x)

    iBin = 1
    for i in range(1,self.binning.numBins()+1):
      hist.GetPoint(i-1,x,y)
      diff = float(y)
      if (float(x) < lowB or float(x) > highB):
        continue
      #print("hist bin: %10i, x: %10.2f, y: %10.2f" % (iBin,float(x),float(y)))
      if x > xCurveMin and x < xCurveMax:
        curvePoint = curve.interpolate(x)
        diff -= curvePoint
        #print(" curve interpolation: %10.2f" % (curvePoint))
      else:
        diff = 0.
        #print(" Warning: x outside of curve range: [ %10.2f %10.2f ]" % (xCurveMin,xCurveMax))
      #print(" diff: %10.2f" % (diff))

      bkgSubHistData.SetBinContent(iBin,float(y))

      bkgSubHistFit.SetBinContent(iBin,curvePoint)
      bkgSubHistFit.SetBinError(iBin,0) # assuming the error on the fit is negligible

      iBin += 1

    # Construct the bkg subtracted histogram
    bkgSubHist = bkgSubHistData.Clone("bkgSubHist")
    bkgSubHist.Sumw2()
    bkgSubHist.Add(bkgSubHistFit,-1)
    
    return bkgSubHist
  

  def drawBkgSub(self,filenameNoExt):

    self.canvas.SetLogy(0)
    self.canvas.cd()
  
    self.frameBkgSub.SetMinimum(-100)
    self.frameBkgSub.SetMaximum(+200)
    self.frameBkgSub.Draw()
    self.bkgSubHist.Draw("histo same")
    self.bkgSubHist.Draw("pe same")
    self.frameBkgSub.Draw("same")
    self.legBkgSub.Draw()

    # Text
    self.tlatex.SetTextSize(0.04)
    self.tlatex.SetTextAlign(12)
    self.tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    self.tlatex.SetTextAlign(32)
    self.tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.title)

    self.tlatex.SetTextAlign(32)
    if (self.lumi != 0):
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.820,self.lumiStr)
    if (self.lumi == 0):
      self.tlatex.SetTextSize(0.04)
      self.tlatex.DrawLatex(self.legPos[0]-0.03,0.850,"#sqrt{s}=7 TeV L =  5.0 fb^{-1} ")
      self.tlatex.DrawLatex(self.legPos[0]-0.03,0.770,"#sqrt{s}=8 TeV L = 19.8 fb^{-1}")
            
    energyStr = self.energyStr
    if re.search(r"[\d]TeV",energyStr):
      energyStr = energyStr.replace("TeV"," TeV")
    if (self.energyStr != ""):
      self.tlatex.DrawLatex(self.legPos[0]-0.01,0.875,"#sqrt{s} = "+self.energyStr)

    saveAs(self.canvas,filenameNoExt)


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
    myCurveHist =  root.TH1F("myCurveHist_"+histPlotName+"_"+pdfPlotName,"",
                          int(nBins), lowB, highB
                          )
    myDataHist =  root.TH1F("myDataHist_"+histPlotName+"_"+pdfPlotName,"",
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
      myDataHist.SetBinContent(iBin,y)
      #print("hist bin: %10i, x: %10.2f, y: %10.2f" % (iBin,float(x),float(y)))
      if x > xCurveMin and x < xCurveMax:
        curvePoint = curve.interpolate(x)
        myCurveHist.SetBinContent(iBin,curvePoint)
        if curvePoint == 0.:
          pull = 0.
        else:
          pull -= curvePoint
          pull /= sqrt(curvePoint)
        #print(" curve interpolation: %10.2f" % (curvePoint))
      else:
        pull = 0.
        #print(" Warning: x outside of curve range: [ %10.2f %10.2f ]" % (xCurveMin,xCurveMax))
      #print(" pull: %10.2f" % (pull))
      pullsHist.SetBinContent(iBin,pull)
      iBin += 1

    self.myCurveHist = myCurveHist
    self.myDataHist = myDataHist
      
    return pullsHist

  def addPDFNormError(self,pad):
    bakCurve = None
    bakErrCurve = None
    for i in pad.GetListOfPrimitives():
      name = i.GetName()
      if re.match(self.pdfName+r"[\d]+",name):
        bakCurve = pad.FindObject(name)
      elif re.match(self.pdfName+r"Error[\d]+",name):
        bakErrCurve = pad.FindObject(name)
    assert(bakCurve)
    assert(bakErrCurve)
    assert(bakCurve.GetN()*2 == bakErrCurve.GetN())
    nPoints = bakCurve.GetN()
    x = root.Double(0.)
    y = root.Double(0.)
    xUp = root.Double(0.)
    yUp = root.Double(0.)
    xDown = root.Double(0.)
    yDown = root.Double(0.)
    for i in range(nPoints):
      bakCurve.GetPoint(i,x,y)
      bakErrCurve.GetPoint(i,xDown,yDown)
      iUp = 2*nPoints-i-1
      bakErrCurve.GetPoint(iUp,xUp,yUp)
      normErr = self.normErr*y
      errUp = float(yUp)-float(y)
      errDown = float(yDown)-float(y)
      errUp2 = errUp**2+normErr**2
      errDown2 = errDown**2+normErr**2
      if errUp < 0:
        errUp = - sqrt(errUp2)
      else:
        errUp = sqrt(errUp2)
      if errDown < 0:
        errDown = - sqrt(errDown2)
      else:
        errDown = sqrt(errDown2)
      bakErrCurve.SetPoint(i,xDown,y+errDown)
      bakErrCurve.SetPoint(iUp,xUp,y+errUp)

def treeCut(category,cutString,eventWeights=True,muonRequirements=True,KDString="KD"):
  global MENormDict
  result = cutString
  if 'KD' in result:
    result = result.replace("KD", KDString)
  if len(result)==0:
    result = "1"
  if len(category) > 0:
    mask = 0
    #if "VBFPresel" in category:
    #    result += " && ((1 & eventType) > 0) && ptMiss < 40."
    #if "VBFBDT" in category:
    #    result += " && ((1 & eventType) > 0) && ptMiss < 40. && bdtVBF>0.0"
    #if "IncPresel" in category:
    #    result += " && ((4 & eventType) > 0)"
    #if "IncBDT" in category:
    #    result += " && ((8 & eventType) > 0)"
    #if "Jet0" in category:
    #    result += " && nJets == 0"
    #if "Jet1" in category:
    #    result += " && nJets == 1"
    #if "Jet2" in category:
    #    result += " && nJets >= 2"
    if "NotBB" in category:
        result += " && ((1024 & eventType) > 0)"
    elif "BB" in category:
        result += " && ((16 & eventType) > 0)"
    if "BO" in category:
        result += " && ((32 & eventType) > 0)"
    if "BE" in category:
        result += " && ((64 & eventType) > 0)"
    if "OO" in category:
        result += " && ((128 & eventType) > 0)"
    if "OE" in category:
        result += " && ((256 & eventType) > 0)"
    if "EE" in category:
        result += " && ((512 & eventType) > 0)"
    #if "FF" in category:
    #    result += " && ((512 & eventType) > 0 || (256 & eventType) > 0)"
    #if "CC" in category:
    #    result += " && ((128 & eventType) > 0 || (64 & eventType) > 0)"
    #if "FC" in category:
    #    result += " && ((128 & eventType) > 0 || (64 & eventType) > 0 || (512 & eventType) > 0 || (256 & eventType) > 0)"
    #if "PtG10" in category:
    #    result += " && (dimuonPt > 10.)"
    #if "PtG20" in category:
    #    result += " && (dimuonPt > 20.)"
    #if "PtG50" in category:
    #    result += " && (dimuonPt > 50.)"
    #if "PtL10" in category:
    #    result += " && (dimuonPt <= 10.)"
    #if "PtL20" in category:
    #    result += " && (dimuonPt <= 20.)"
    #if "PtL50" in category:
    #    result += " && (dimuonPt <= 50.)"
    #if "VBFCutBased" in category:
    #    result += " && ((1 & eventType) > 0) && deltaEtaJets>3.5 && dijetMass>500. && ptMiss<40."
  if muonRequirements:
    #result += " && muonLead_pt>25. && muonSub_pt>25."
    #result += " && muonLead_passPFRelIso && muonSub_passPFRelIso && (muonLead_isHltMatched || muonSub_isHltMatched)"
    result += " && muonLead_passPFRelIso && muonSub_passPFRelIso && ("#(muonLead_isHltMatched || muonSub_isHltMatched)"
    result += " (muonLead_pt>25 && muonLead_isHltMatched) || "
    result += " (muonSub_pt >25 &&  muonSub_isHltMatched) "
    result += ")"

  if eventWeights:
    result = "("+result+")*puWeight"
  return result

def drawStandardCaptions(canvas,caption1,caption2="",caption3="",caption4="",caption5="",preliminaryString=PRELIMINARYSTRING):
  tlatex = root.TLatex()
  tlatex.SetNDC()

  tlatex.SetTextFont(root.gStyle.GetLabelFont())
  tlatex.SetTextSize(0.04)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,preliminaryString)

  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption1)
  tlatex.SetTextAlign(12)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.88,caption2)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.82,caption3)
  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.76,caption4)
  tlatex.SetTextAlign(32)
  tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.88,caption5)
  return tlatex

def copyTreeBranchToNewNameTree(tree,oldBranchName,newBranchName):
  """
  Returns a new tree with the contents of oldBranchName in the old tree, but with
  the branch name newBranchName

  Assumes both branches just contain floats!!
  """
  result = root.TTree(tree.GetName()+"_newNames","")

  newVal = array.array( 'f', [ 0. ] ) # one element array so we get a pointer to the value
  newBranch = result.Branch( newBranchName, newVal, newBranchName+'/F' )

  nEntries = tree.GetEntries()

  for i in range(nEntries):
    tree.GetEntry(i)
    oldVal = getattr(tree,oldBranchName)
    newVal[0] = oldVal
    #print i,newVal[0]
    result.Fill()
  newBranch.SetAddress(0)
  return result
  
def rooDebugFR(fr,printCorrelation=True):
  result = ""
  result += "{0:20}: {1:10.3g}\n".format("Status",fr.status())
  result += "{0:20}: {1:10.3g}\n".format("EDM",fr.edm())
  result += "{0:20}: {1:10.3g}\n".format("NLL",fr.minNll())
  result += "{0:20}: {1:10.3g}\n".format("Cov. Matrix Qual",fr.edm())
  result += "{0:20}: {1:10}\n".format("N Invalid NLL",fr.numInvalidNLL())
  result += "Status History:\n"
  for i in range(fr.numStatusHistory()):
    result += "  {1} ({0})\n".format(fr.statusCodeHistory(i),fr.statusLabelHistory(i))

  rowString = ""
  rowString += "{0:30}".format("Variable Name")
  rowString += "{0:10} +/- {1:<10}  ".format("Fit Value","HESSE Err")
  rowString += " {0:<10} {1:<10}  ".format("High Err","Low Err")

  rowString += "{0:20}".format("[{0}]".format("Variable Limits"))
  result += rowString + "\n"

  fpf_s = fr.floatParsFinal()
  for i in range(fpf_s.getSize()):
    nuis_s = fpf_s.at(i)
    name   = nuis_s.GetName();
    # .getVal() .getError() .getMin() .getMax()

    rowString = ""
    rowString += "{0:30}".format(re.sub(r".*TeV_","",name))
    rowString += "{0:10.3g} +/- {1:<10.3g}  ".format(nuis_s.getVal(),nuis_s.getError())
    rowString += "+{0:<10.3g} {1:<10.3g}  ".format(nuis_s.getErrorHi(),nuis_s.getErrorLo())

    rowString += "{0:20}".format("[{0:.1g},{1:.1g}]".format(nuis_s.getMin(),nuis_s.getMax()))
    result += rowString + "\n"
  constPars = fr.constPars()
  for i in range(constPars.getSize()):
    nuis_s = constPars.at(i)
    name   = nuis_s.GetName();
    # .getVal() .getError() .getMin() .getMax()

    rowString = ""
    rowString += "{0:30}".format(re.sub(r".*TeV_","",name))
    rowString += "{0:10.3g}  FIXED  ".format(nuis_s.getVal())
    result += rowString + "\n"
    #print name, nuis_s.getVal()
  if printCorrelation:
    result += "Correlations:\n"

    rowString = " "*12
    for j in range(fpf_s.getSize()): # Columns
      name = fpf_s.at(j).GetName()
      if '_' in name:
        name = name[-[i for i in reversed(name)].index("_"):]
      rowString += "{0:10}  ".format(name)
    result += rowString + "\n"
    for i in range(fpf_s.getSize()): # Rows
      name = fpf_s.at(i).GetName()
      if '_' in name:
        name = name[-[k for k in reversed(name)].index("_"):]
      rowString = "{0:10}  ".format(name)
      for j in range(fpf_s.getSize()): # Columns
        corr = fr.correlation(fpf_s.at(i),fpf_s.at(j))
        rowString += "{0:10.3g}  ".format(corr)
      result += rowString + "\n"
  return result 

def rooCalcChi2(pdf,dataHist):
  """
  Returns (chi^2,ndf,chi^2 Probablility)
  """
  assert(type(dataHist) == root.RooDataHist)
  observables = rooArgSet2List(pdf.getObservables(dataHist))
  assert(len(observables)==1)
  binning = observables[0].getBinning()
  nBins = binning.numBins()
  pdfParams = rooArgSet2List(pdf.getParameters(dataHist))
  nparams = len(pdfParams)
  nparamsFree = 0
  for p in pdfParams:
    if not p.isConstant():
      nparamsFree += 1
  ndf = nBins - nparamsFree
  #for errsForChi2, errsForChi2Label in zip([root.RooAbsData.Expected,root.RooAbsData.SumW2,root.RooAbsData.Poisson],["Errors from PDF","Errors Data Weights^2","Errors Poisson"]):
  #for errsForChi2, errsForChi2Label in zip([getattr(root.RooAbsData,"None"),root.RooAbsData.Auto,root.RooAbsData.SumW2,root.RooAbsData.Poisson],["Errors None","Errors Auto","Errors Data Weights^2","Errors Poisson"]):
  chi2Var = pdf.createChi2(dataHist,root.RooFit.DataError(root.RooAbsData.Poisson))
  chi2 = chi2Var.getVal()
  chi2Prob = scipy.stats.chi2.sf(chi2,ndf)
  return (chi2,ndf,chi2Prob)

def rooDebugChi2(pdf,data):
  result = ""
  observables = rooArgSet2List(pdf.getObservables(data))
  assert(len(observables)==1)
  binning = observables[0].getBinning()
  nBins = binning.numBins()
  pdfParams = rooArgSet2List(pdf.getParameters(data))
  nparams = len(pdfParams)
  nparamsFree = 0
  for p in pdfParams:
    if not p.isConstant():
      nparamsFree += 1
  ndf = nBins - nparamsFree
  #for errsForChi2, errsForChi2Label in zip([root.RooAbsData.Expected,root.RooAbsData.SumW2,root.RooAbsData.Poisson],["Errors from PDF","Errors Data Weights^2","Errors Poisson"]):
  #for errsForChi2, errsForChi2Label in zip([getattr(root.RooAbsData,"None"),root.RooAbsData.Auto,root.RooAbsData.SumW2,root.RooAbsData.Poisson],["Errors None","Errors Auto","Errors Data Weights^2","Errors Poisson"]):
  for errsForChi2, errsForChi2Label in zip([root.RooAbsData.Poisson],["Errors Poisson"]):
    chi2Var = pdf.createChi2(data,root.RooFit.DataError(errsForChi2))
    chi2 = chi2Var.getVal()
    chi2Prob = scipy.stats.chi2.sf(chi2,ndf)
    normChi2 = chi2/ndf
    result += "{0}:\n".format("Chi2 Info Using "+errsForChi2Label+" "+str(errsForChi2))
    result += "{0:20}: {1:<10.3g}\n".format("chi2",chi2)
    result += "{0:20}: {1:<10.3g}\n".format("chi2/ndf",normChi2)
    result += "{0:20}: {1:<10.3g}\n".format("chi2 Prob",chi2Prob)
    result += "{0:20}: {1:<10.3g}\n".format("ndf",ndf)
    result += "{0:20}: {1:<10.3g}\n".format("nBins",nBins)
    result += "{0:20}: {1:<10.3g}\n".format("nParams Free",nparamsFree)
    result += "{0:20}: {1:<10.3g}\n".format("nParams",nparams)
  return result

class RooPredictionIntervalPlotter:
  def __init__(self,xVar,pdf,data,fr,nToys=1000):
    self.xVar = xVar
    self.xBinning = self.xVar.getBinning()
    self.nBins = self.xBinning.numBins()
    self.xMin = self.xBinning.lowBound()
    self.xMax = self.xBinning.highBound()
    self.binWidth = self.xBinning.averageBinWidth()
    self.pdf = pdf
    self.data = data
    self.fr = fr
    self.nToys = nToys
    self.nData = data.sumEntries()
    self.pdfE = pdf
    self.toyHists = []
    if not self.pdfE.InheritsFrom("RooExtendPdf"):
        self.nDataVar = root.RooRealVar(pdf.GetName()+"_extendPdf_nDataVar","",self.nData)
        self.pdfE = root.RooExtendPdf(pdf.GetName()+"_extendPdf","",self.pdf,self.nDataVar)
    self.studyPdf, self.paramPdf = self.makeStudyPdf(self.pdfE,fr)
    self.constraintParamSet = root.RooArgSet(*rooArgSet2List(fr.floatParsFinal()))
    self.mcStudy = root.RooMCStudy(self.studyPdf,root.RooArgSet(self.xVar),
                                    root.RooFit.Binned(True),
                                    root.RooFit.Silence(),
                                    root.RooFit.Extended(),
                                    root.RooFit.FitOptions(
                                            root.RooFit.Save(True),
                                            root.RooFit.PrintEvalErrors(0)
                                    ),
                                    root.RooFit.Constrain(self.constraintParamSet),
    )
    self.mcStudy.generate(self.nToys,0,True)
    #self.mcStudy.generateAndFit(self.nToys,0,True)  ## Helps to be able to debug constriant parameters
    #fitParamDataSet = self.mcStudy.fitParDataSet()
    for iToy in range(nToys):
      toyDataSet = self.mcStudy.genData(iToy)
      toyHist = toyDataSet.createHistogram(self.pdf.GetName()+"toyPredTolHist{0}".format(iToy),self.xVar,
                        root.RooFit.Binning(self.nBins,self.xMin,self.xMax)
            )
      self.toyHists.append(toyHist)

  def drawPrediction(self,sigmas=1,color=root.kRed-9,drawOpt="2"):
    predGraph = root.TGraphAsymmErrors()
    predGraph.SetMarkerColor(color)
    predGraph.SetLineColor(color)
    predGraph.SetFillColor(color)
    for iBin in range(self.nBins):
      x = self.xBinning.binCenter(iBin)
      xLow = x-self.xBinning.binLow(iBin)
      xHigh = self.xBinning.binHigh(iBin)-x
      yToyList = numpy.zeros(self.nToys)
      for iToy in range(self.nToys):
        yToy = self.toyHists[iToy].GetBinContent(iBin+1)
        yToyList[iToy] = yToy
      quantiles = numpy.percentile(yToyList,
                    [
                        100.*scipy.stats.norm.cdf(-sigmas),
                        50.,
                        100.*scipy.stats.norm.cdf(sigmas)
                    ]
      )
      y = quantiles[1]
      yLow = y-quantiles[0]
      yHigh = quantiles[2]-y
      predGraph.SetPoint(iBin,x,y)
      predGraph.SetPointError(iBin,xLow,xHigh,yLow,yHigh)
    self.predGraph = predGraph
    self.predGraph.Draw(drawOpt)
    # For debugging
    #for hist in self.toyHists:
    #    hist.SetLineColor(root.kGreen+1)
    #    hist.Draw("same hist")

  def drawStatOnly(self,sigmas=1,color=root.kMagenta-9,drawOpt="2"):
    statErrGraph = root.TGraphAsymmErrors()
    statErrGraph = root.TGraphAsymmErrors()
    statErrGraph.SetMarkerColor(color)
    statErrGraph.SetLineColor(color)
    statErrGraph.SetFillColor(color)
    for iBin in range(self.nBins):
      x = self.xBinning.binCenter(iBin)
      xLow = x-self.xBinning.binLow(iBin)
      xHigh = self.xBinning.binHigh(iBin)-x
      yToyList = numpy.zeros(self.nToys)
      for iToy in range(self.nToys):
        yToy = self.toyHists[iToy].GetBinContent(iBin+1)
        yToyList[iToy] = yToy
      y = numpy.median(yToyList)
      statErrGraph.SetPoint(iBin,x,y)
      statErrGraph.SetPointError(iBin,xLow,xHigh,sqrt(y),sqrt(y))
    self.statErrGraph = statErrGraph
    self.statErrGraph.Draw(drawOpt)
    # For debugging
    #for hist in self.toyHists:
    #    hist.SetLineColor(root.kGreen+1)
    #    hist.Draw("same hist")


  def makeStudyPdf(self,pdf,fr):
    fitFloatParams = rooArgSet2List(fr.floatParsFinal())
    paramPdf = fr.createHessePdf(root.RooArgSet(*fitFloatParams))
    #paramPdf = root.RooGaussian("sillyJustin","",fitFloatParams[0],root.RooFit.RooConst(fitFloatParams[0].getVal()),root.RooFit.RooConst(fitFloatParams[0].getVal()*0.5))
    result = root.RooProdPdf(pdf.GetName()+"withConstrainedHessePdf","",pdf,paramPdf)
    return result, paramPdf

def rooLinearErrorPropagation(pdf,fr,xVar,parameters,nEvents,nPoints=100):
  errGraph = root.TGraphAsymmErrors()
  binning = xVar.getBinning()
  xMin = binning.lowBound()
  xMax = binning.highBound()
  binWidth = binning.averageBinWidth()
  pointWidth = (xMax-xMin)/(nPoints-1)
  observables = root.RooArgSet(xVar)
  parList = rooArgSet2List(parameters)
  parListFR = rooArgSet2List(fr.floatParsFinal())
  assert(len(parList) == len(parListFR))
  for par,parFR in zip(parList,parListFR):
    par.setVal(parFR.getVal())
  relStatUnc = nEvents**(-0.5)
  for iPoint in range(nPoints):
    xNow = iPoint*pointWidth
  
    rangeName = "linearErrPropRange_{0}".format(iPoint)
    xVar.setRange(rangeName,xNow-binWidth/2.,xNow+binWidth/2.)
    pdfInt =  pdf.createIntegral(observables,observables,rangeName)
    pdfVal = pdfInt.getVal()
    errGraph.SetPoint(iPoint,xNow,nEvents * pdfVal)
  
    pdfErrParList = []
    for par,parFR in zip(parList,parListFR):
      originalVal = parFR.getVal()
      paramUnc = parFR.getError()
  
      par.setVal(originalVal+paramUnc)
      pdfErrVal = pdfInt.getVal()
      pointErr = abs(pdfErrVal-pdfVal)*nEvents
      par.setVal(originalVal-paramUnc)
      pointErr = abs(pdfErrVal-pdfVal)*nEvents
      pointErr = max(abs(pdfErrVal-pdfVal)*nEvents,pointErr)
  
      pdfErrParList.append(pointErr)
      par.setVal(originalVal)
  
    totalVariance = pdfVal*relStatUnc # Include Stat unc of normalization
    for i in range(len(parList)):
      iParam = parList[i]
      iPdfErr = pdfErrParList[i]
      for j in range(i,len(parList)):
        jParam = parList[j]
        jPdfErr = pdfErrParList[j]
        totalVariance += iPdfErr*jPdfErr*fr.correlation(iParam,jParam)
      
    totalUncertainty = sqrt(totalVariance)
    errGraph.SetPointError(iPoint,binWidth/2.,binWidth/2.,totalUncertainty,totalUncertainty)
  return errGraph

if __name__ == "__main__":

  root.gROOT.SetBatch(True)
  print("Running helpers.py")
  

  zs = []
  f = numpy.mean
  for iToy in range(100):
    sample = numpy.random.normal(0.,1.,(100))
    est = f(sample)
    jnEst, jnErr, jnBias = jacknife(sample,f)
    zs.append((est)/jnErr)

  print("Mean and std Z-score for 100 experiments: ",numpy.mean(zs),numpy.std(zs))
  #print "median: ",numpy.median(sample)
  #print jacknife(sample,numpy.median)

