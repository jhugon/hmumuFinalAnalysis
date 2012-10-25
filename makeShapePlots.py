#!/usr/bin/python

from helpers import *

root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

dataDir = "statsCards/"
outDir = ""

class ShapePlotter:
  def __init__(self,filename,titleMap):
    self.titleMap = titleMap
    self.filename = filename
    self.data = {}
    self.f = root.TFile(filename)
    for channelKey in self.f.GetListOfKeys():
      if channelKey.GetClassName() != "TDirectoryFile":
        continue
      channelFolder = channelKey.ReadObj()
      debugFolderKey = channelFolder.GetKey("debug")
      if debugFolderKey.GetClassName() != "TDirectoryFile":
        print("Warning: loadData: {}/{} in {} is not a directory".format(channelFolder.GetName(),debugFolder.GetName(),filename))
        continue
      debugFolder = debugFolderKey.ReadObj()
      self.data[channelKey.GetName()] = {}
      for key in debugFolder.GetListOfKeys():
        if re.match(r"TH1.*",key.GetClassName()) and (re.match(r"bak.*",key.GetName()) or key.GetName() == "data_obs" or key.GetName() == "sig"):
          hist = key.ReadObj()
          hist = hist.Clone(channelKey.GetName()+"_"+key.GetName())
          self.data[channelKey.GetName()][key.GetName()] = hist

    self.lumi = -1
    self.lumiStr = ""
    tmpMatch = re.search(r"([\w]*)_([0-9]+)\.root",filename)
    if tmpMatch:
      self.lumi = int(tmpMatch.group(2))
      self.lumiStr = "L = {} fb^{{-1}}".format(self.lumi)


    self.colors = [root.kRed-9, root.kGreen-9, root.kBlue-9, root.kMagenta-9, root.kCyan-9]
    self.fillStyles = [3004,3005,3003,3006,3007]

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
      y = (yu+yd)/2.0
      yu = yu-y
      yd = y-yd
      #print("x,y: {:.3g},{:.3g} + {:.3g} - {:.3g}".format(x,y, yu, yd))
      outGraph.SetPoint(iGraph,x,y)
      outGraph.SetPointEYhigh(iGraph,yu)
      outGraph.SetPointEYlow(iGraph,yd)
      iGraph += 1

  def makePlot(self):
    canvas = root.TCanvas("canvas")
    for channelName in self.data:
      canvas.Clear()
      channel = self.data[channelName]
      nominal = channel["bak"]
      obs = channel["data_obs"]
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
      obs.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
      obs.GetYaxis().SetTitle("Events/Bin")
      nominal.SetFillStyle(0)
      nominal.SetLineStyle(1)
      nominal.SetLineColor(root.kRed+1)
      obs.Draw("")
      nominal.Draw("hist L same")
      graphs = []
      for param,col,sty in zip(paramSet,self.colors,self.fillStyles):
        tmp = root.TGraphAsymmErrors()
        up = channel["bak_"+param+"Up"]
        down = channel["bak_"+param+"Down"]
        self.makeGraph(nominal,up,down,tmp)
        tmp.SetFillStyle(sty)
        tmp.SetFillColor(col)
        tmp.SetLineColor(col)
        tmp.SetMarkerColor(col)
        tmp.Draw("4")
        graphs.append(tmp)
      nominal.Draw("hist L same")
      obs.Draw("same")

      tlatex = root.TLatex()
      tlatex.SetNDC()
      tlatex.SetTextFont(root.gStyle.GetLabelFont())
      #tlatex.SetTextSize(0.05)
      tlatex.SetTextSize(0.04)
      tlatex.SetTextAlign(12)
      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
      tlatex.SetTextAlign(32)
      tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,self.titleMap[channelName]+", "+self.lumiStr)
        
      canvas.SaveAs(channelName+".png")

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
  "VBFVeryTight":"VBFVT"
}
        
s = ShapePlotter(dataDir+"IncPresel_20.root",titleMap)
#print s.data
s.makePlot()
