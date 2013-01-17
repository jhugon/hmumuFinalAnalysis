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

def getGraphIntegral(graph):
  X = root.Double()
  Y = root.Double()
  result = 0.0
  for i in range(graph.GetN()):
    graph.GetPoint(i,X,Y)
    result += float(Y)
  return result

class ShapePlotter:
  def __init__(self,filename,outDir,titleMap,rebin=1,doSignalScaling=True,xlimits=[]):
    self.histList = []
    self.padList = []
    self.pullsList = []
    self.legList = []
    self.modelLines = []
    self.titleMap = titleMap
    self.filename = filename
    self.xlimits = xlimits
    self.textFileName = os.path.splitext(filename)[0]+".txt"
    self.processNameMap, self.params = getattr(self,"readCard")(self.textFileName)
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
      bakPDF = channelWS.pdf("bak")
      data_obs = channelWS.data("data_obs")
      rooDataTitle = data_obs.GetTitle()
      if rebin != 1:
        tmpHist = data_obs.createHistogram("mMuMu")
        tmpHist.Rebin(rebin)
        data_obs =  root.RooDataHist(data_obs.GetName(),data_obs.GetTitle(),
                                                            root.RooArgList(mMuMu),tmpHist)

      #Data Time
      dataGraph, bakPDFGraph, pullsGraph,chi2 = getattr(self,"makeTGraphs")(bakPDF,data_obs,mMuMu)
      getattr(self,"draw")(channelName,dataGraph,bakPDFGraph,pullsGraph,chi2,rooDataTitle)
      saveAs(self.canvas,outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName)

      #Templates Time
      for processName in self.processNameMap[channelName]:
        template = channelWS.data(processName+"_Template")
        rooDataTitle = template.GetTitle()
        if rebin != 1:
            tmpHist = template.createHistogram("mMuMu")
            tmpHist.Rebin(rebin)
            template =  root.RooDataHist(template.GetName(),template.GetTitle(),
                                                            root.RooArgList(mMuMu),tmpHist)
        pdf = channelWS.pdf(processName)
        dataGraph, pdfGraph, pullsGraph,chi2 = getattr(self,"makeTGraphs")(pdf,template,mMuMu)
        getattr(self,"draw")(channelName,dataGraph,pdfGraph,pullsGraph,chi2,rooDataTitle)
        saveAs(self.canvas,outDir+os.path.splitext(os.path.split(self.filename)[1])[0]+'_'+channelName+"_"+processName)

  def readCard(self,fn):
    f = open(fn)
    foundBin = False
    binList = []
    processList = []
    rateList = []
    paramMap = {}
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
    binList = [x.replace(r" ","") for x in binList]
    processList = [x.replace(r" ","") for x in processList]
    result = {}
    for i in binList:
      if not result.has_key(i):
        result[i] = {}
    for b,p,r in zip(binList,processList,rateList):
      result[b][p] = r
    return result, paramMap

  def getRatioGraph(self,hist,curve):
    def findI(graph,x):
      xList = []
      X = root.Double()
      Y = root.Double()
      for i in range(graph.GetN()):
        graph.GetPoint(i,X,Y)
        xList.append(abs(float(X)-float(x)))
      bestI = min(range(len(xList)), key=lambda i: xList[i])
      if xList[bestI]>1.0:
        return -1
      return bestI
    ratioPlot = root.TGraphAsymmErrors()
    histX = root.Double()
    histY = root.Double()
    curveX = root.Double()
    curveY = root.Double()
    for i in range(hist.GetN()):
      hist.GetPoint(i,histX,histY)
      curveI = findI(curve,histX)
      if curveI < 0:
        continue
      curve.GetPoint(curveI,curveX,curveY)
      ratio = 0.0
      ratioErrUp = 0.0
      ratioErrDown = 0.0
      if curveY != 0.0 and histY != 0.0:
        histYErrUp = hist.GetErrorYhigh(i)
        histYErrDown = hist.GetErrorYlow(i)
        ratio = histY/curveY
        ratioErrUp = histYErrUp/curveY
        ratioErrDown = histYErrDown/curveY
        #if "bak" in curve.GetName():
          #print("hx: {:<10.1f} cx: {:<10.1f} hy: {:<10.1f} cy: {:<10.1f} ratio: {:<10.3f} hyErrUp: {:<10.1f}".format(histX,curveX,histY,curveY,ratio,hist.GetErrorYhigh(i)))
          #print("{:<5.1f}: hx: {:<10.1f} cx: {:<10.1f} ratio: {:<10.3f}".format(histX,histX,curveX,ratio))
          #print("{:<5.1f}: h: {:<10.1f} c: {:<10.1f} ratio: {:<10.3f}".format(histX,histY,curveY,ratio))
      else:
        continue
      ratioPlot.SetPoint(i,histX,ratio)
      ratioPlot.SetPointError(i,0.,0.,ratioErrUp,ratioErrDown)
    return ratioPlot

  def makeTGraphs(self,pdf,data,observable):
    avgWeight = 1.0
    isRealData = (data.GetTitle() == "Real Observed Data")
    if (not isRealData) and (data.GetName() == "data_obs" or "bak" in data.GetName()):
      for i in backgroundList:
        if "DY" in i:
          datasetName = i + "_" + self.energyStr
          avgWeight = self.lumi*1000.0*xsec[datasetName]/nEventsMap[datasetName]
    frame = None
    if len(self.xlimits) == 2:
      frame = observable.frame(root.RooFit.Range(*self.xlimits))
    else:
      frame = observable.frame()
    pdfParams = pdf.getParameters(data)
    itr = pdfParams.createIterator()
    curveNomName = pdf.GetName()+"_CurveNom"
    curveDataName = data.GetName()+"_Curve"
    rng = root.RooFit.Range(*self.xlimits)
    if "Hmumu125" in data.GetName():
      rng = root.RooFit.Range(110.,140.)
    data.plotOn(frame,root.RooFit.Name(curveDataName))
    pdf.plotOn(frame,root.RooFit.Name(curveNomName),rng)
    pulls = frame.pullHist(curveDataName,curveNomName)
    chi2 = frame.chiSquare()
    curveNames = []
    for iParam in range(pdfParams.getSize()):
      param = itr.Next()
      paramName = param.GetName()
      if self.params.has_key(paramName):
        nominal,err = self.params[paramName]
        nominal = float(nominal)
        err = float(err)
        tmpName = pdf.GetName()+"_Curve{0}p".format(iParam)
        param.setVal(nominal+err)
        pdf.plotOn(frame,root.RooFit.Name(tmpName),rng)
        curveNames.append(tmpName)
        tmpName = pdf.GetName()+"_Curve{0}m".format(iParam)
        param.setVal(nominal-err)
        pdf.plotOn(frame,root.RooFit.Name(tmpName),rng)
        curveNames.append(tmpName)
    varUp = []
    varDown = []
    for curveName in curveNames:
      curve = frame.findObject(curveName)
      if curveName[len(curveName)-1]=='p':
        varUp.append(curve)
      else:
        varDown.append(curve)
    curveNom = frame.findObject(curveNomName)
    curveData = frame.findObject(curveDataName)
    modelGraph = root.TGraphAsymmErrors()
    dataGraph = root.TGraphAsymmErrors()
    pullsGraph = root.TGraphAsymmErrors()
    iPoint = 0
    xNom = root.Double()
    yNom = root.Double()
    xErr = root.Double()
    yErr = root.Double()
    for i in range(curveNom.GetN()):
      curveNom.GetPoint(i,xNom,yNom)
      errUp = 0.0
      errDown = 0.0
      badPoint = False
      for errCurve in varUp:
        iErr = errCurve.findPoint(xNom,1.0)
        if iErr < 0:
          badPoint = True
        errCurve.GetPoint(iErr,xErr,yErr)
        errUp += (yErr-float(yNom))**2
      for errCurve in varDown:
        iErr = errCurve.findPoint(xNom,1.0)
        if iErr < 0:
          badPoint = True
        errCurve.GetPoint(iErr,xErr,yErr)
        errDown += (yErr-float(yNom))**2
      if badPoint:
        continue
      errUp = sqrt(errUp)
      errDown = sqrt(errDown)
      modelGraph.SetPoint(iPoint,xNom,yNom)
      modelGraph.SetPointError(iPoint,0.,0.,errUp,errDown)
      iPoint+=1
    iPoint = 0
    for i in range(curveData.GetN()):
      curveData.GetPoint(i,xNom,yNom)
      if yNom <=0.0:
        continue
      dataGraph.SetPoint(iPoint,xNom,yNom)
      errUp = sqrt(yNom)*sqrt(avgWeight)
      errDown = sqrt(yNom)*sqrt(avgWeight)
      dataGraph.SetPointError(iPoint,0.,0.,errUp,errDown)
      iPoint+=1
    iPoint = 0
    for i in range(pulls.GetN()):
      pulls.GetPoint(i,xNom,yNom)
      pullsGraph.SetPoint(iPoint,xNom,yNom)
      errUp = pulls.GetErrorXhigh(i)
      errDown = pulls.GetErrorXlow(i)
      pullsGraph.SetPointError(iPoint,0.,0.,errUp,errDown)
      iPoint+=1
    dataGraph.SetName(curveData.GetName()+"_TGraph")
    modelGraph.SetName(curveNom.GetName()+"_TGraph")
    pullsGraph.SetName(curveData.GetName()+"_pullsTGraph")
    dataGraph.SetMarkerColor(1)
    modelGraph.SetLineColor(root.kBlue+1)
    modelGraph.SetFillColor(root.kCyan)
    modelGraph.SetFillStyle(1)
    return dataGraph, modelGraph, pullsGraph, chi2

  def draw(self,channelName,data,model,pulls,chi2,rooDataTitle):
      drawXLimits = self.xlimits
      if "Hmumu125" in data.GetName():
        drawXLimits = [110.,140.]
      def getMaxYVal(self,graph):
        l = []
        x = root.Double()
        y = root.Double()
        for i in range(graph.GetN()):
          graph.GetPoint(i,x,y)
          if x < drawXLimits[0]:
                continue
          if x > drawXLimits[1]:
                continue
          l.append(float(y))
        if len(l) == 0:
            return -1.0
        return max(l)
      def copyGraphNoErrs(graph,outGraph):
        x = root.Double()
        y = root.Double()
        for i in range(graph.GetN()):
          graph.GetPoint(i,x,y)
          outGraph.SetPoint(i,x,y)
        outGraph.SetLineColor(graph.GetLineColor())
        outGraph.SetLineStyle(graph.GetLineStyle())
        outGraph.SetLineWidth(graph.GetLineWidth())
        outGraph.SetMarkerColor(graph.GetMarkerColor())
        outGraph.SetMarkerStyle(graph.GetMarkerStyle())
        outGraph.SetMarkerSize(graph.GetMarkerSize())
        self.modelLines.append(outGraph)
      def makeRatioModelPlots(graph,outRatio,outOne):
        for i in range(graph.GetN()):
          x = root.Double()
          y = root.Double()
          graph.GetPoint(i,x,y)
          outOne.SetPoint(i,x,1.0)
          outRatio.SetPoint(i,x,1.0)
          if y == 0.0:
            continue
          errUp = graph.GetErrorYhigh(i)
          errDown = graph.GetErrorYlow(i)
          outRatio.SetPointError(i,0.,0.,errDown/float(y),errUp/float(y))
        outOne.SetLineColor(graph.GetLineColor())
        outOne.SetLineStyle(graph.GetLineStyle())
        outOne.SetLineWidth(graph.GetLineWidth())
        outRatio.SetFillColor(graph.GetFillColor())
        outRatio.SetFillStyle(graph.GetFillStyle())
        self.modelLines.append(outOne)
        self.modelLines.append(outRatio)
      def getHistArgsFromGraph(graph,drawXLimits):
        if graph.GetN() < 2:
          return graph.GetName(),graph.GetTitle(),100,drawXLimits[0],drawXLimits[1]
        diffs = []
        x = root.Double()
        y = root.Double()
        for i in range(graph.GetN()-1):
          graph.GetPoint(i,x,y)
          x1 = float(x)
          graph.GetPoint(i+1,x,y)
          x2 = float(x)
          diffs.append(x2-x1)
        #return max(set(diffs), key=diffs.count)
        nBins = int((drawXLimits[1]-drawXLimits[0])/min(diffs))
        return graph.GetName(),graph.GetTitle(),nBins,drawXLimits[0],drawXLimits[1]
      def getHistFromGraph(graph,hist):
        x = root.Double()
        y = root.Double()
        for i in range(graph.GetN()-1):
          graph.GetPoint(i,x,y)
          yErr = graph.GetErrorYhigh(i)
          iBin = hist.GetXaxis().FindBin(x)
          hist.SetBinContent(iBin,y)
          hist.SetBinError(iBin,yErr)
        hist.SetMarkerColor(1)
        hist.SetLineColor(1)

      th1Args = getHistArgsFromGraph(data,drawXLimits)
      dataHist = root.TH1F(*th1Args)
      self.histList.append(dataHist)
      pullHist = root.TH1F(*th1Args)
      self.histList.append(pullHist)

      getHistFromGraph(data,dataHist)

      maxVal = getMaxYVal(self,data)
      maxVal = max(getMaxYVal(self,model),maxVal)
      dataLabel = "FullSim Data"
      if rooDataTitle == "Toy Data":
        dataLabel = "Toy Data"
      elif rooDataTitle == "Real Observed Data":
        dataLabel = "Data"

      #Setup Canvas
      self.canvas.cd()
      self.canvas.Clear()
      self.canvas.SetLogy(0)
      pad1 = root.TPad("pad1"+data.GetName(),"",0.02,0.30,0.98,0.98,0)
      pad2 = root.TPad("pad2"+data.GetName(),"",0.02,0.01,0.98,0.29,0)
    
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
      dataHist.Draw("")
      dataHist.GetYaxis().SetTitle("Events/Bin")
      dataHist.GetYaxis().SetTitleSize(gStyle.GetTitleSize("Y")*canvasToPad1FontScalingFactor)
      dataHist.GetYaxis().SetLabelSize(gStyle.GetLabelSize("Y")*canvasToPad1FontScalingFactor)
      dataHist.GetYaxis().SetTitleOffset(0.9*gStyle.GetTitleOffset("Y"))
      dataHist.GetXaxis().SetLabelOffset(0.70)
      dataHist.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
      #dataHist.GetXaxis().SetRangeUser(*drawXLimits)
      if maxVal != None:
        dataHist.GetYaxis().SetRangeUser(0.0,maxVal*1.05)
      modelLine = root.TGraph()
      copyGraphNoErrs(model,modelLine)
      model.Draw("3")
      modelLine.Draw("l")
      dataHist.Draw("same")
      pad1.Update()
      pad1.RedrawAxis() # Updates Axis Lines

      # Pulls Pad
      pad2.cd()

      ratioGraph = self.getRatioGraph(data,model)
      self.pullsList.append(ratioGraph)
      doRatio = True
      if doRatio:
        pulls =  ratioGraph
      getHistFromGraph(pulls,pullHist)
      pullHist.SetLineStyle(1)
      pullHist.SetLineColor(1)
      pullHist.SetMarkerColor(1)
      pullHist.SetTitle("")
      pullHist.GetXaxis().SetRangeUser(*drawXLimits)
      pullHist.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
      pullHist.GetXaxis().CenterTitle(1)
      pullHist.GetYaxis().SetNdivisions(5)
      pullHist.GetXaxis().SetTitleSize(0.055*pad1ToPad2FontScalingFactor)
      pullHist.GetXaxis().SetLabelSize(0.050*pad1ToPad2FontScalingFactor)
      pullHist.GetYaxis().SetTitle("#frac{"+dataLabel+"-Fit}{#sigma_{"+dataLabel+"}}")
      if doRatio:
        pullHist.GetYaxis().SetTitle("#frac{"+dataLabel+"}{Fit}")
        pullHist.GetYaxis().SetRangeUser(0,2)
      pullHist.GetYaxis().SetTitleSize(0.045*pad1ToPad2FontScalingFactor)
      pullHist.GetYaxis().SetLabelSize(0.045*pad1ToPad2FontScalingFactor)
      pullHist.GetYaxis().CenterTitle(1)
      pullHist.GetXaxis().SetTitleOffset(0.75*pullHist.GetXaxis().GetTitleOffset())
      pullHist.GetYaxis().SetTitleOffset(0.70)

      pullHist.GetXaxis().SetLabelSize(gStyle.GetLabelSize("X")*canvasToPad2FontScalingFactor)

      pullHist.Draw("")
      if doRatio:
        modelRatio = root.TGraphAsymmErrors()
        modelOne = root.TGraphAsymmErrors()
        makeRatioModelPlots(model,modelRatio,modelOne)
        modelRatio.Draw("3")
        modelOne.Draw("l")
        pullHist.Draw("same")


      ## Pretty Stuff
  
      normchi2 = chi2
      problatex = root.TLatex()
      problatex.SetNDC()
      problatex.SetTextFont(dataHist.GetXaxis().GetLabelFont())
      problatex.SetTextSize(pullHist.GetYaxis().GetLabelSize())
      problatex.SetTextAlign(12)
      problatex.DrawLatex(0.18,0.41,"#chi^{{2}}/NDF: {0:.3g}".format(normchi2))
  
      pad2.Update()
      pad2.RedrawAxis() # Updates Axis Lines
    
      pad1.cd()

      legPos = [0.65,0.65,1.0-gStyle.GetPadRightMargin()-0.01,1.0-gStyle.GetPadTopMargin()-0.01]
      leg = root.TLegend(*legPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      leg.AddEntry(dataHist,dataLabel,"pe")
      leg.AddEntry(model,"Model","lf")
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

  plotRange= [110,160]
  #plotRange= []

  rebin=1

  shapePlotterList = []
  #for fn in glob.glob(dataDir+"*20.root")+glob.glob(dataDir+"*5.05.root"):
  for fn in glob.glob(dataDir+"*.root"):
    s = ShapePlotter(fn,outDir,titleMap,rebin,xlimits=plotRange)
    shapePlotterList.append(s)
