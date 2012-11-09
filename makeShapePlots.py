#!/usr/bin/env python

from helpers import *
import math
import os.path
import glob

root.gErrorIgnoreLevel = root.kWarning
root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

def convertErrors(hist):
  nBinsX = hist.GetNbinsX()
  for i in range(0,nBinsX+2):
    x = hist.GetBinContent(i)
    hist.SetBinError(i,math.sqrt(x))

class ShapePlotter:
  def __init__(self,filename,titleMap,rebin=1):
    self.titleMap = titleMap
    self.filename = filename
    self.data = {}
    self.f = root.TFile(filename)
    for channelKey in self.f.GetListOfKeys():
      if channelKey.GetClassName() != "TDirectoryFile":
        continue
      channelFolder = channelKey.ReadObj()
      debugFolderKey = channelFolder.GetKey("debug")
      try:
       if debugFolderKey.GetClassName() != "TDirectoryFile":
          print("Warning: loadData: {0}/{1} in {2} is not a directory".format(channelFolder.GetName(),debugFolder.GetName(),filename))
          continue
      except Exception as e:
          print("Warning: loadData: {0}/{1} in {2} Does not exist".format(channelFolder.GetName(),"debug",filename))
          continue
      debugFolder = debugFolderKey.ReadObj()
      self.data[channelKey.GetName()] = {}
      for key in debugFolder.GetListOfKeys():
        if re.match(r"TH1.*",key.GetClassName()) and (re.match(r"bak.*",key.GetName()) or key.GetName() == "data_obs" or key.GetName() == "sig"):
          hist = key.ReadObj()
          hist = hist.Clone(channelKey.GetName()+"_"+key.GetName())
          hist.Rebin(rebin)
          self.data[channelKey.GetName()][key.GetName()] = hist

    self.lumi = -1
    self.lumiStr = ""
    self.energyStr = ""
    tmpMatch = re.search(r"([\w]*)_(.+)_([0-9]+)\.root",filename)
    if tmpMatch:
      self.lumi = int(tmpMatch.group(3))
      self.lumiStr = "L = {0} fb^{{-1}}".format(self.lumi)
      self.energyStr = tmpMatch.group(2)

    self.colors = [root.kRed-9, root.kGreen-9, root.kBlue-9, root.kMagenta-9, root.kCyan-9]
    self.fillStyles = [3004,3005,3003,3006,3007]

    self.limit = 1.0
    #if len(self.data)==1:
    if True:
        scaleSignal = getattr(self,"scaleSignal")
        scaleSignal()
    self.graphsListList = []
    self.padList = []

  def scaleSignal(self):
    limitfname = ""
    try:
      cardTime = os.path.getmtime(self.filename)
      limitfname = os.path.splitext(self.filename)[0] + ".txt.out"
      limitfname = limitfname.replace("statsCards","statsInput")
      limitTime = os.path.getmtime(limitfname)
    except Exception as e:
      print("Warning: Couldn't find limit file to scale signal by: {0}".format(limitfname))
      #print("Error: {0}, {1}".format(e.errno,e.strerror))
      return
    if limitTime > cardTime:
      try:
        f = open(limitfname)
        for line in f:
          #match = re.search(r"Observed[\s]Limit:[^.\d]*< ([.\deE]+)",line)
          match = re.search(r"Expected.*50\.0.*< ([.\deE]+)",line)
          if match:
            self.limit = float(match.group(1))
            break
      except Exception as e:
        print("Warning: Couldn't find limit file to scale signal by: {0}".format(limitfname))
        #print("Error: {0}, {1}".format(e.errno,e.strerror))
        return
      if self.limit != 1.0:
        for channelName in self.data:
          channel = self.data[channelName]
          if channel.has_key("sig"):
            channel["sig"].Scale(self.limit)

  def hist2Graph(self,hist,outGraph):
    assert(hist.InheritsFrom("TH1"))
    assert(outGraph.InheritsFrom("TGraphAsymmErrors") or outGraph.InheritsFrom("TGraphErrors") 
            or outGraph.InheritsFrom("TGraph"))
    iGraph = 0
    for i in range(1,hist.GetNbinsX()+1):
      x = hist.GetXaxis().GetBinCenter(i)
      y = hist.GetBinContent(i)
      if y<= 0:
        continue
      outGraph.SetPoint(iGraph,x,y)
      iGraph += 1

  def makeGraph(self,nominal,up,down,outGraph):
    assert(nominal.InheritsFrom("TH1"))
    assert(up.InheritsFrom("TH1"))
    assert(down.InheritsFrom("TH1"))
    assert(outGraph.InheritsFrom("TGraphAsymmErrors"))
    assert(nominal.GetNbinsX()==up.GetNbinsX())
    assert(nominal.GetNbinsX()==down.GetNbinsX())
    
    iGraph = 0
    for i in range(1,nominal.GetNbinsX()+1):
      x = nominal.GetXaxis().GetBinCenter(i)
      yu = up.GetBinContent(i)
      yd = down.GetBinContent(i)
      if yd > yu:
        tmp = yu
        yu = yd
        yd = tmp
      #y = (yu+yd)/2.0
      y = nominal.GetBinContent(i)
      yu = yu-y
      yd = y-yd
      if y<=0:
        continue
      #print("x,y: {0:.3g},{1:.3g} + {2:.3g} - {3:.3g}".format(x,y, yu, yd))
      outGraph.SetPoint(iGraph,x,y)
      outGraph.SetPointEYhigh(iGraph,yu)
      outGraph.SetPointEYlow(iGraph,yd)
      iGraph += 1

  def combineErrors(self,listOfGraphs,outGraph):
    N = listOfGraphs[0].GetN()
    for g in listOfGraphs:
      assert(g.GetN() == N)
    for i in range(N):
      err2Up = 0.0
      err2Down = 0.0
      y = root.Double()
      x = root.Double()
      for g in listOfGraphs:
        err2Up += g.GetErrorYhigh(i)**2
        err2Down += g.GetErrorYlow(i)**2
        g.GetPoint(i,x,y)
      if y <= 0.0:
        continue
      outGraph.SetPoint(i,x,y)
      outGraph.SetPointError(i,0.0,0.0,math.sqrt(err2Down),math.sqrt(err2Up))

  def getPlots(self,channelName):
      channel = self.data[channelName]
      nominal = channel["bak"]
      obs = channel["data_obs"]
      convertErrors(obs)
      sig = channel["sig"]
      sig.Add(nominal)
      sigGraph = root.TGraph()
      self.hist2Graph(sig,sigGraph)
      sigGraph.SetLineColor(root.kRed)
      paramSet = set()
      for histName in channel:
        if histName != "bak" and histName != "data_obs" and histName != "sig":
          tmp = histName.replace("bak_","")
          tmp = tmp.replace("Up","")
          tmp = tmp.replace("Down","")
          if not (tmp in paramSet):
            paramSet.add(tmp)
      obs.SetFillStyle(0)
      obs.SetLineStyle(1)
      obs.SetLineStyle(1)
      obs.SetLineColor(1)
      obs.SetTitle("")
      obs.GetXaxis().SetTitle("")
      if len(plotRange) ==2:
        obs.GetXaxis().SetRangeUser(*plotRange)
      obs.GetYaxis().SetTitle("Events/Bin")
      graphs = []
      for param,col,sty in zip(paramSet,self.colors,self.fillStyles):
        tmp = root.TGraphAsymmErrors()
        up = channel["bak_"+param+"Up"]
        down = channel["bak_"+param+"Down"]
        self.makeGraph(nominal,up,down,tmp)
        #tmp.SetFillStyle(sty)
        #tmp.SetFillColor(col)
        #tmp.SetLineColor(col)
        #tmp.SetMarkerColor(col)
        #tmp.Draw("4")
        graphs.append(tmp)
      self.graphsListList.append(graphs)
      combinedErrorGraph = root.TGraphAsymmErrors()
      combinedErrorGraph.SetFillColor(root.kCyan)
      combinedErrorGraph.SetLineColor(root.kBlue+1)
      self.combineErrors(graphs,combinedErrorGraph)
      combinedErrorGraph.Draw("3")
      sigGraph.Draw("C")
      combinedErrorGraph.Draw("CX")

      normchi2 = nominal.Chi2Test(obs,"UW CHI2/NDF")
      print("Warning: chi2 is computed for weighted data!!!! Not suitable for real data!!!")

      # Make Pull Hist
      pullYmax = 0.0
      pullYmin = 10.0
      pullHist = obs.Clone("pullHist_"+obs.GetName())
      pullHist.Reset()
      iGraph = 0
      pullGraph = root.TGraphAsymmErrors()
      for i in range(1,pullHist.GetNbinsX()):
        x = obs.GetXaxis().GetBinCenter(i)
        nObs = obs.GetBinContent(i)
        nNom = nominal.GetBinContent(i)
        nObsErr = obs.GetBinError(i)
        nNomErr = nominal.GetBinError(i)
        ratio = 0.0
        ratioErrUp = 0.0
        ratioErrDown = 0.0
        if nObs != 0.0 and nNom != 0.0:
          ratio = nObs/nNom
          ratioErrUp = (nObs+nObsErr)/nNom - ratio
          ratioErrDown = ratio - (nObs-nObsErr)/nNom
        pullGraph.SetPoint(iGraph,x,ratio)
        pullGraph.SetPointError(iGraph,0.,0.,ratioErrUp,ratioErrDown)
        pullHist.SetBinContent(i,ratio)
        if pullYmax < ratio + ratioErrUp:
            pullYmax = ratio + ratioErrUp
        if pullYmin > ratio - ratioErrDown:
            pullYmin = ratio - ratioErrDown
        iGraph += 1

      pullErrs = root.TGraphAsymmErrors()
      for i in range(0,combinedErrorGraph.GetN()):
        x = root.Double()
        y = root.Double()
        combinedErrorGraph.GetPoint(i,x,y)
        yErrUp = combinedErrorGraph.GetErrorYhigh(i)/y
        yErrDown = combinedErrorGraph.GetErrorYlow(i)/y
        pullErrs.SetPoint(i,x,1.0)
        pullErrs.SetPointError(i,0.0,0.0,yErrDown,yErrUp)
        if pullYmax < 1. + yErrUp:
            pullYmax = 1. + yErrUp
        if pullYmin > 1. - yErrDown:
            pullYmin = 1. - ratioErrDown
      pullHist.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
      pullHist.GetYaxis().SetRangeUser(pullYmin,pullYmax)
      pullErrs.SetLineColor(root.kBlue)
      pullErrs.SetFillColor(root.kCyan)

      return obs, combinedErrorGraph, sigGraph, pullHist, pullGraph, pullErrs, normchi2

  def makePlot(self,outDir,plotRange=[]):
    dataLabel = "MC Data"
    canvas = root.TCanvas("canvas")
    for channelName in self.data:
      canvas.Clear()
      obs, combinedErrorGraph, sigGraph, pulls, pullGraph, pullErrGraph, normchi2 = self.getPlots(channelName)

      #Setup Canvas
      canvas.cd()
      pad1 = root.TPad("pad1"+obs.GetName(),"",0.02,0.30,0.98,0.98,0)
      pad2 = root.TPad("pad2"+obs.GetName(),"",0.02,0.01,0.98,0.29,0)
    
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
      if False:
          pad1.SetLogy(1)
      else:
          pad1.SetLogy(0)
    
      pad1.Draw() # Projections pad
      pad2.Draw() # Residuals   pad
  
      pad1Width = pad1.XtoPixel(pad1.GetX2())
      pad1Height = pad1.YtoPixel(pad1.GetY1())
      pad2Height = pad2.YtoPixel(pad2.GetY1())
      pad1ToPad2FontScalingFactor = float(pad1Height)/pad2Height
      canvasToPad1FontScalingFactor = float(canvas.YtoPixel(canvas.GetY1()))/pad1.YtoPixel(pad1.GetY1())
      canvasToPad2FontScalingFactor = float(canvas.YtoPixel(canvas.GetY1()))/pad2.YtoPixel(pad2.GetY1())
    
      # Main Pad
      pad1.cd();
      #obs.GetYaxis().SetTitleSize(0.045)
      #obs.GetYaxis().SetLabelSize(0.045)
      obs.GetYaxis().SetTitleSize(gStyle.GetTitleSize("Y")*canvasToPad1FontScalingFactor)
      obs.GetYaxis().SetLabelSize(gStyle.GetLabelSize("Y")*canvasToPad1FontScalingFactor)
      obs.GetYaxis().SetTitleOffset(0.9*gStyle.GetTitleOffset("Y"))
      obs.GetXaxis().SetLabelOffset(0.70)
      obs.Draw("")
      combinedErrorGraph.Draw("3")
      sigGraph.Draw("C")
      combinedErrorGraph.Draw("CX")
      obs.Draw("same")
      pad1.Update()
  
      pad1.RedrawAxis() # Updates Axis Lines

      # Pulls Pad
      pad2.cd()
      pulls.SetTitle("")
      if plotRange != []:
        pulls.GetXaxis().SetRangeUser(*plotRange)
      #pulls.GetXaxis().SetTitle(xtitle)
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

      pulls.Draw("p")
      #pulls.GetYaxis().SetRangeUser(0,2)
      pulls.Draw("p")
      pullGraph.Draw("p")
      pullErrGraph.Draw("3")
      pullErrGraph.Draw("CX")
      pullGraph.Draw("p")
  
      problatex = root.TLatex()
      problatex.SetNDC()
      problatex.SetTextFont(obs.GetXaxis().GetLabelFont())
      problatex.SetTextSize(pulls.GetYaxis().GetLabelSize())
      problatex.SetTextAlign(12)
      problatex.DrawLatex(0.18,0.41,"#chi^{{2}}/NDF: {0:.3g}".format(normchi2))
  
      pad2.Update()
      pad2.GetFrame().DrawClone()
      pad2.RedrawAxis() # Updates Axis Lines
    

      """
      obs.Draw("")
      combinedErrorGraph.Draw("3")
      sigGraph.Draw("C")
      combinedErrorGraph.Draw("CX")
      obs.Draw("same")

      legPos = [0.65,0.65,1.0-gStyle.GetPadRightMargin()-0.01,1.0-gStyle.GetPadTopMargin()-0.01]
      leg = root.TLegend(*legPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      leg.AddEntry(obs,"MC Data","pe")
      leg.AddEntry(combinedErrorGraph,"Background Model","lf")
      leg.AddEntry(sigGraph,"SM Higgs #times {0:.1f}".format(self.limit),"l")
      leg.Draw()

      tlatex = root.TLatex()
      tlatex.SetNDC()
      tlatex.SetTextFont(root.gStyle.GetLabelFont())
      #tlatex.SetTextSize(0.05)
      tlatex.SetTextSize(0.04)
      tlatex.SetTextAlign(12)
      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.titleMap[channelName])
      tlatex.DrawLatex(legPos[0]-0.01,0.875,self.lumiStr)
      tlatex.DrawLatex(legPos[0]-0.01,0.83,"#sqrt{s}=8 TeV")
      tlatex.SetTextSize(0.035)
      tlatex.SetTextAlign(11)
      tlatex.DrawLatex(gStyle.GetPadLeftMargin()+0.025,gStyle.GetPadBottomMargin()+0.025,
        "#chi^{{2}}/NDF = {0:.3g}".format(normchi2)
      )
        
      canvas.RedrawAxis()
      canvas.SaveAs(outDir+"/"+channelName+".png")
      """

      pad1.cd()

      legPos = [0.65,0.65,1.0-gStyle.GetPadRightMargin()-0.01,1.0-gStyle.GetPadTopMargin()-0.01]
      leg = root.TLegend(*legPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      leg.AddEntry(obs,dataLabel,"pe")
      leg.AddEntry(combinedErrorGraph,"Background Model","lf")
      leg.AddEntry(sigGraph,"SM Higgs #times {0:.1f}".format(self.limit),"l")
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

      canvas.SaveAs(outDir+"/"+channelName+"_"+self.energyStr+".png")

      self.padList.extend([canvas,pad1,pad2])

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

  "BDTSig80":"BDT Cut Combination",
  "IncBDTSig80":"Inclusive BDT Cut",
  "VBFBDTSig80":"VBF BDT Cut"
}
        
if __name__ == "__main__":
  dataDir = "statsCards/"
  outDir = "shapes/"

  plotRange= [115,150]
  #plotRange= []

  rebin=1

  for fn in glob.glob(dataDir+"*20.root"):
    #print fn
    if fn.count("Cat")>0:
      continue
    if fn.count("/BDTSig80")>0:
      continue
    s = ShapePlotter(fn,titleMap,rebin)
    s.makePlot(outDir,plotRange)
