#!/usr/bin/env python

from helpers import *
from xsec import *
import math
import os.path
import glob

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

def convertErrors(hist):
  nBinsX = hist.GetNbinsX()
  for i in range(0,nBinsX+2):
    x = hist.GetBinContent(i)
    hist.SetBinError(i,math.sqrt(x))

def findMaxInRange(dataHist,xName,minX,maxX):
  result = 0.0
  for i in range(dataHist.numEntries()):
    argSet = dataHist.get(i)
    xVar = argSet.find(xName)
    x = xVar.getVal()
    if x>maxX or x<minX:
        continue
    y = dataHist.weight()
    if y>result:
      result = y
  return result

class ShapePlotter:
  def __init__(self,filename,outDir,titleMap,rebin=1,doSignalScaling=True,xlimits=[]):
    self.padList = []
    self.pullsList = []
    self.legList = []
    self.titleMap = titleMap
    self.filename = filename
    self.textFileName = os.path.splitext(filename)[0]+".txt"
    self.textData = getattr(self,"readCard")(self.textFileName)
    self.colors = [root.kRed-9, root.kGreen-9, root.kBlue-9, root.kMagenta-9, root.kCyan-9]
    self.fillStyles = [3004,3005,3003,3006,3007]

    self.canvas = root.TCanvas()

    self.lumi = -1
    self.lumiStr = ""
    self.energyStr = ""
    tmpMatch = re.search(r"([\w]*)_(.+)_([.0-9]+)\.root",filename)
    if tmpMatch:
      self.lumi = int(float(tmpMatch.group(3)))
      self.lumiStr = "L = {0} fb^{{-1}}".format(self.lumi)
      self.energyStr = tmpMatch.group(2)

    self.data = {}
    self.f = root.TFile(filename)
    for channelKey in self.f.GetListOfKeys():
      if channelKey.GetClassName() != "RooWorkspace":
        continue
      channelName = channelKey.GetName()
      channelWS = channelKey.ReadObj()
      mMuMu = channelWS.var("mMuMu")
      mMuMu.SetTitle("m_{#mu#mu} [GeV]")
      bakPDF = channelWS.pdf("bak")
      data_obs = channelWS.data("data_obs")
      maxVal = None

      #Data Time
      if len(xlimits) ==2:
         frame = mMuMu.frame(root.RooFit.Range(*xlimits))
         maxVal = findMaxInRange(data_obs,mMuMu.GetName(),xlimits[0],xlimits[1])
         frame.SetMaximum(maxVal)
      else:
         frame = mMuMu.frame()
      frame.SetTitle("")
      data_obs.plotOn(frame)
      bakPDF.plotOn(frame,root.RooFit.Range(*xlimits))
      #bakPDF.plotOn(frame,root.RooFit.LineStyle(2),root.RooFit.Range(*xlimits))
      chi2 = bakPDF.createChi2(data_obs)
      getattr(self,"draw")(frame,xlimits,channelName,maxVal)
      saveAs(self.canvas,outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName)

      #Templates Time
      for processName in self.textData[channelName]:
        template = channelWS.data(processName+"_Template")
        pdf = channelWS.pdf(processName)
        tmpXlims = xlimits
        if "125" in processName:
            tmpXlims = [110.,140.]
        if len(tmpXlims) ==2:
           frame = mMuMu.frame(root.RooFit.Range(*tmpXlims))
           maxVal = findMaxInRange(template,mMuMu.GetName(),tmpXlims[0],tmpXlims[1])
           frame.SetMaximum(maxVal)
        else:
           frame = mMuMu.frame()
        frame.SetTitle("")
        template.plotOn(frame)
        pdf.plotOn(frame,root.RooFit.Range(*tmpXlims))
        #if processName == "bak":
        #  bakPDF.plotOn(frame,root.RooFit.LineStyle(2),root.RooFit.Range(*xlimits))
        chi2 = pdf.createChi2(template)
        getattr(self,"draw")(frame,tmpXlims,channelName,maxVal)
        saveAs(self.canvas,outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName+"_"+processName)

  def readCard(self,fn):
    f = open(fn)
    foundBin = False
    binList = []
    processList = []
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
    binList = [x.replace(r" ","") for x in binList]
    processList = [x.replace(r" ","") for x in processList]
    result = {}
    for i in binList:
      if not result.has_key(i):
        result[i] = []
    for b,p in zip(binList,processList):
      result[b].append(p)
    return result

  def getRatioGraph(self,frame):
    def findI(graph,x):
      xList = []
      X = root.Double()
      Y = root.Double()
      for i in range(graph.GetN()):
        graph.GetPoint(i,X,Y)
        xList.append(float(X))
      return min(range(len(xList)), key=lambda i: xList[i])
    curve = None
    hist = None
    for i in range(2):
        obj = frame.getObject(i)
        if obj.InheritsFrom("RooCurve"):
            curve = obj
        elif obj.InheritsFrom("RooHist"):
            hist = obj
    assert(hist != None)
    assert(curve != None)
    #frame.Print()
    #curve.Print()
    #hist.Print()
    ratioPlot = root.TGraphAsymmErrors()
    histX = root.Double()
    histY = root.Double()
    curveX = root.Double()
    curveY = root.Double()
    for i in range(hist.GetN()):
      hist.GetPoint(i,histX,histY)
      #curveI = findI(curve,histX)
      curveI = curve.findPoint(histX,1.0)
      if curveI < 0:
        continue
      curve.GetPoint(curveI,curveX,curveY)
      ratio = 0.0
      ratioErrUp = 0.0
      ratioErrDown = 0.0
      if curveY != 0.0 and histY != 0.0:
        histYErrUp = hist.GetErrorYhigh(i)/histY
        histYErrDown = hist.GetErrorYlow(i)/histY
        ratio = histY/curveY
        ratioErrUp = (histY+histYErrUp)/curveY - ratio
        ratioErrDown = ratio - (histY-histYErrDown)/curveY
        #if "bak" in curve.GetName():
          #print("{:<5.1f}: hx: {:<10.1f} cx: {:<10.1f} ratio: {:<10.3f}".format(histX,histX,curveX,ratio))
          #print("{:<5.1f}: h: {:<10.1f} c: {:<10.1f} ratio: {:<10.3f}".format(histX,histY,curveY,ratio))
      else:
        continue
      ratioPlot.SetPoint(i,histX,ratio)
      ratioPlot.SetPointError(i,0.,0.,ratioErrUp,ratioErrDown)
    return ratioPlot

  def draw(self,frame,xlimits,channelName,maxVal=None):
      dataLabel = "MC Data"

      #Setup Canvas
      self.canvas.cd()
      self.canvas.Clear()
      self.canvas.SetLogy(0)
      pad1 = root.TPad("pad1"+frame.GetName(),"",0.02,0.30,0.98,0.98,0)
      pad2 = root.TPad("pad2"+frame.GetName(),"",0.02,0.01,0.98,0.29,0)
    
      pad1.SetBottomMargin(0.005);
      pad2.SetTopMargin   (0.005);
      pad2.SetBottomMargin(0.33);
      pad2.SetLogy(0)
      pad1.SetLogy(0)
    
      pad1.Draw() # Projections pad
      pad2.Draw() # Residuals   pad
  
      pad1Width = pad1.XtoPixel(pad1.GetX2())
      pad1Height = pad1.YtoPixel(pad1.GetY1())
      pad2Height = pad2.YtoPixel(pad2.GetY1())
      pad1ToPad2FontScalingFactor = float(pad1Height)/pad2Height
      canvasToPad1FontScalingFactor = float(self.canvas.YtoPixel(self.canvas.GetY1()))/pad1.YtoPixel(pad1.GetY1())
      canvasToPad2FontScalingFactor = float(self.canvas.YtoPixel(self.canvas.GetY1()))/pad2.YtoPixel(pad2.GetY1())
    
      # Main Pad
      pad1.cd();
      frame.Draw()
      frame.GetYaxis().SetTitle("Events/Bin")
      frame.GetYaxis().SetTitleSize(gStyle.GetTitleSize("Y")*canvasToPad1FontScalingFactor)
      frame.GetYaxis().SetLabelSize(gStyle.GetLabelSize("Y")*canvasToPad1FontScalingFactor)
      frame.GetYaxis().SetTitleOffset(0.9*gStyle.GetTitleOffset("Y"))
      frame.GetXaxis().SetLabelOffset(0.70)
      frame.Draw("")
      if maxVal != None:
        frame.GetYaxis().SetRangeUser(0.0,maxVal*1.05)
      pad1.Update()
      pad1.RedrawAxis() # Updates Axis Lines

      # Pulls Pad
      pad2.cd()

      pulls = getattr(self,"getRatioGraph")(frame)

      pulls.SetLineStyle(1)
      pulls.SetLineColor(1)
      pulls.SetMarkerColor(1)
      pulls.SetTitle("")
      pulls.GetXaxis().SetRangeUser(*xlimits)
      pulls.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
      pulls.GetXaxis().CenterTitle(1)
      pulls.GetYaxis().SetNdivisions(5)
      pulls.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
      pulls.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
      pulls.GetYaxis().SetTitle("#frac{"+dataLabel+"}{Fit}")
      pulls.GetYaxis().SetTitleSize(0.045*pad1ToPad2FontScalingFactor)
      pulls.GetYaxis().SetLabelSize(0.045*pad1ToPad2FontScalingFactor)
      pulls.GetYaxis().CenterTitle(1)
      pulls.GetXaxis().SetTitleOffset(0.75*pulls.GetXaxis().GetTitleOffset())
      pulls.GetYaxis().SetTitleOffset(0.70)

      pulls.GetXaxis().SetLabelSize(gStyle.GetLabelSize("X")*canvasToPad2FontScalingFactor)

      pulls.Draw("ape")

      ## Pretty Stuff
  
      normchi2 = frame.chiSquare()
      problatex = root.TLatex()
      problatex.SetNDC()
      problatex.SetTextFont(frame.GetXaxis().GetLabelFont())
      problatex.SetTextSize(pulls.GetYaxis().GetLabelSize())
      problatex.SetTextAlign(12)
      problatex.DrawLatex(0.18,0.41,"#chi^{{2}}/NDF: {0:.3g}".format(normchi2))
  
      pad2.Update()
      pad2.GetFrame().DrawClone()
      pad2.RedrawAxis() # Updates Axis Lines
    
      pad1.cd()

      curve = None
      hist = None
      for i in range(2):
        obj = frame.getObject(i)
        if obj.InheritsFrom("RooCurve"):
            curve = obj
        elif obj.InheritsFrom("RooHist"):
            hist = obj
      assert(hist != None)
      assert(curve != None)

      legPos = [0.65,0.65,1.0-gStyle.GetPadRightMargin()-0.01,1.0-gStyle.GetPadTopMargin()-0.01]
      leg = root.TLegend(*legPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      leg.AddEntry(hist,dataLabel,"pe")
      leg.AddEntry(curve,"Model","l")
      #leg.AddEntry(combinedSigErrorGraph,"SM Higgs #times {0:.1f}".format(self.limit),"lf")
      leg.Draw()

      tlatex = root.TLatex()
      tlatex.SetNDC()
      tlatex.SetTextFont(root.gStyle.GetLabelFont())
      #tlatex.SetTextSize(0.05)
      tlatex.SetTextSize(0.04*canvasToPad1FontScalingFactor)
      tlatex.SetTextAlign(12)
      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.titleMap[channelName])
      tlatex.DrawLatex(legPos[0]-0.01,0.875,self.lumiStr)
      tlatex.DrawLatex(legPos[0]-0.01,0.82,"#sqrt{s}="+self.energyStr)

      self.padList.extend([pad1,pad2])
      self.pullsList.append(pulls)
      self.legList.append(leg)

titleMap = {
  "AllCat":"All Categories Comb.",
  "IncCat":"Inclusive Categories Comb.",
  "VBFCat":"VBF Categories Comb.",

  "IncPresel":"Inclusive Preselection",
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
  "IncBDTCut":"Inclusive BDT Cut",
  "VBFBDTCut":"VBF BDT Cut",

  "BDTCutCat":"BDT Cut Cat. Combination",
  "IncBDTCutCat":"Inclusive BDT Cut",
  "VBFBDTCutCat":"VBF BDT Cut",

  "IncPreselCat":"Inclusive Cat. Preselection",
  "VBFPreselCat":"VBF Cat. Preselection",

  "IncBDTCutBB":"Inclusive BDT Cut BB",
  "IncBDTCutBO":"Inclusive BDT Cut BO",
  "IncBDTCutBE":"Inclusive BDT Cut BE",
  "IncBDTCutOO":"Inclusive BDT Cut OO",
  "IncBDTCutOE":"Inclusive BDT Cut OE",
  "IncBDTCutEE":"Inclusive BDT Cut EE",
  "IncBDTCutNotBB":"Inclusive BDT Cut !BB",
  "VBFBDTCutBB":"VBF BDT Cut BB",
  "VBFBDTCutNotBB":"VBF BDT Cut !BB",
  "IncPreselBB":"Inclusive Preselection BB",
  "IncPreselBO":"Inclusive Preselection BO",
  "IncPreselBE":"Inclusive Preselection BE",
  "IncPreselOO":"Inclusive Preselection OO",
  "IncPreselOE":"Inclusive Preselection OE",
  "IncPreselEE":"Inclusive Preselection EE",
  "IncPreselNotBB":"Inclusive Preselection !BB",
  "VBFPreselBB":"VBF Preselection BB",
  "VBFPreselNotBB":"VBF Preselection !BB"
}
        
if __name__ == "__main__":
  dataDir = "statsCards/"
  outDir = "shapes/"

  plotRange= [105,160]
  #plotRange= []

  rebin=1

  for fn in glob.glob(dataDir+"*20.root")+glob.glob(dataDir+"*5.05.root"):
    s = ShapePlotter(fn,outDir,titleMap,rebin,xlimits=plotRange)
