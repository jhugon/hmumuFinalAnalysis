#!/usr/bin/python

from helpers import *

root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

dataDir = "statsCards/"
outDir = ""

class ShapePlotter:
  def __init__(self,filename):
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
        if re.match(r"TH1.*",key.GetClassName()) and (re.match(r"bak.*",key.GetName()) or key.GetName() == "data_obs"):
          hist = key.ReadObj()
          hist = hist.Clone(channelKey.GetName()+"_"+key.GetName())
          hist.Print()
          self.data[channelKey.GetName()][key.GetName()] = hist

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
      print("x,y: {:.3g},{:.3g} + {:.3g} - {:.3g}".format(x,y, yu, yd))
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
        if histName != "bak" and histName != "data_obs":
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
        
      canvas.SaveAs(channelName+".png")
    
        
s = ShapePlotter(dataDir+"IncPresel_20.root")
print s.data
s.makePlot()
