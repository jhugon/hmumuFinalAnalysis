#!/usr/bin/env python

from xsec import *
from helpers import *

root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

dataDir = "input/"
outDir = "output/"

histPrefix = ""

histlist = [
    "BDTHistMuonOnlyVMass",
    "BDTHistVBFVMass",
    "likelihoodHistMuonOnlyVMass",
    "likelihoodHistVBFVMass"#,
    #"yVmDiMu",
    #"ptVmDiMu",
    #"phiVmDiMu"
]

drawOpt = "colz"

xlabels = {
"yVptDiMu":"Dimuon p_{T} [GeV]",
"BDTHistMuonOnlyVMass":"Dimuon Mass [GeV]",
"BDTHistVBFVMass":"Dimuon Mass [GeV]",
"likelihoodHistMuonOnlyVMass":"Dimuon Mass [GeV]",
"likelihoodHistVBFVMass":"Dimuon Mass [GeV]",
    "yVmDiMu":"Dimuon Mass [GeV]",
    "ptVmDiMu":"Dimuon Mass [GeV]",
    "phiVmDiMu":"Dimuon Mass [GeV]"
}

ylabels = {
"yVptDiMu":"Dimuon |y|",
"BDTHistMuonOnlyVMass":"BDT Output",
"BDTHistVBFVMass":"BDT Output",
"likelihoodHistMuonOnlyVMass":"Likelihood",
"likelihoodHistVBFVMass":"Likelihood",
    "yVmDiMu":"Dimuon |y|",
    "ptVmDiMu":"Dimuon p_{T} [GeV]",
    "phiVmDiMu":"Dimuon #phi"
}

xranges = {
"yVptDiMu":[0,200],
"BDTHistMuonOnlyVMass":[70,150],
"BDTHistVBFVMass":[70,150],
"likelihoodHistMuonOnlyVMass":[70,150],
"likelihoodHistVBFVMass":[70,150],
    "yVmDiMu":[70,150],
    "ptVmDiMu":[70,150],
    "phiVmDiMu":[70,150]
}

yranges = {
"yVptDiMu":[0,2.1],
"BDTHistMuonOnlyVMass":[-1,0],
"BDTHistVBFVMass":[-0.5,0.5],
"likelihoodHistMuonOnlyVMass":[-0.5,0.5],
"likelihoodHistVBFVMass":[-0.5,0.5],
    "yVmDiMu":[0,2.2],
    "ptVmDiMu":[0,250],
    "phiVmDiMu":[0,3.2]
}

rebins = {
"yVptDiMu":[1,2],
"BDTHistMuonOnlyVMass":[10,20],
"BDTHistVBFVMass":[10,20],
"likelihoodHistMuonOnlyVMass":[10,20],
"likelihoodHistVBFVMass":[10,20],
    "yVmDiMu":[10,10],
    "ptVmDiMu":[10,10],
    "phiVmDiMu":[10,10]
}

fSigList = []
for fname in signalList:
  fSigList.append((root.TFile(dataDir+fname+".root"),xsec[fname],nEventsMap[fname]))

fBakList = []
for fname in backgroundList:
  fBakList.append((root.TFile(dataDir+fname+".root"),xsec[fname],nEventsMap[fname]))

fDataList = []
for fname in dataList:
  fDataList.append(root.TFile(dataDir+fname+".root"))

for hist in histlist:
  hSig = None
  hBak = None
  hData = None
  first = True
  for tup in fSigList:
    f = tup[0]
    xsec = tup[1]
    nEvents = tup[2]
    if first:
      hSig = f.Get(histPrefix+hist)
      hSig.Scale(LUMI*xsec/nEvents)
      first = False
    else:
      tmp = f.Get(histPrefix+hist)
      tmp.Scale(LUMI*xsec/nEvents)
      hSig.Add(tmp)
  first = True
  for tup in fBakList:
    f = tup[0]
    xsec = tup[1]
    nEvents = tup[2]
    if first:
      hBak = f.Get(histPrefix+hist)
      hBak.Scale(LUMI*xsec/nEvents)
      first = False
    else:
      tmp = f.Get(histPrefix+hist)
      tmp.Scale(LUMI*xsec/nEvents)
      hBak.Add(tmp)
  first = True
  for f in fDataList:
    if first:
      hData = f.Get(hist)
      first = False
    else:
      tmp = f.Get(hist)
      hData.Add(tmp)
  
  nSlices = 40
  rebinX = hBak.GetNbinsX()/nSlices
  rebinY = 100
  rebinY = 20
  hBak.RebinX(rebinX)
  hSig.RebinX(rebinX)
  hData.RebinX(rebinX)
  hBak.RebinY(rebinY)
  hSig.RebinY(rebinY)
  hData.RebinY(rebinY)
  print("{}: nBins: {}, {}".format(hist,hBak.GetNbinsX(),hBak.GetNbinsY()))

  slicesSig = []
  slicesBak = []
  slicesData = []
  colors = [root.kRed,root.kBlue,root.kGreen,root.kPink,root.kCyan,root.kMagenta]
  colors.extend([x+2 for x in colors]+[x-7 for x in colors])
  colors.extend([x+2 for x in colors]+[x-7 for x in colors])
  maxValue = 0.0
  legendEntries = []
  leg = root.TLegend(0.7,0.65,0.89,0.89)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  for i in range(1,nSlices):
    if hBak.GetXaxis().GetBinLowEdge(i)<110.0:
      continue
    elif hBak.GetXaxis().GetBinUpEdge(i)>150.0:
      continue
    legendEntry = "m_{{#mu#mu}} #in [{:.0f},{:.0f}]".format(hSig.GetXaxis().GetBinLowEdge(i),hSig.GetXaxis().GetBinUpEdge(i))

    tmp = getXBinHist( hSig, i)
    tmp.SetMarkerColor(colors[i])
    tmp.SetMarkerStyle(33)
    tmp.SetLineWidth(0)
    tmp.SetMarkerSize(2.5)
    if tmp.GetMaximum()>maxValue:
        maxValue = tmp.GetMaximum()
    slicesSig.append(tmp)
    leg.AddEntry(tmp,legendEntry,"p")
     
    tmp = getXBinHist( hBak, i)
    tmp.SetMarkerColor(colors[i])
    tmp.SetMarkerStyle(22)
    tmp.SetLineWidth(0)
    tmp.SetMarkerSize(2.5)
    if tmp.GetMaximum()>maxValue:
        maxValue = tmp.GetMaximum()
    slicesBak.append(tmp)

    tmp = getXBinHist( hData, i)
    tmp.SetMarkerColor(colors[i])
    tmp.SetMarkerStyle(20)
    tmp.SetMarkerSize(2.5)
    tmp.SetLineWidth(0)
    if tmp.GetMaximum()>maxValue:
        maxValue = tmp.GetMaximum()
    slicesData.append(tmp)

  canvas = root.TCanvas()
  canvas.SetLogy(1)
  #slicesBak[0].GetYaxis().SetRangeUser(1e-6,maxValue)
  slicesBak[0].GetXaxis().SetTitle(ylabels[hist])
  slicesBak[0].GetYaxis().SetTitle("Counts")
  slicesBak[0].Draw("hist M")
  for i in slicesBak:
    i.Draw("hist same M")
  for i in slicesSig:
    i.Draw("hist same M")
  #for i in slicesData:
  #  i.Draw("same")

  leg.Draw()

  saveAs(canvas,outDir+"slices_"+hist)

  ####################################
  # background shapes

  canvas.SetLogy(0)

  maxValue = 0.0
  for i in slicesBak:
    normalizeHist(i)
    i.SetLineStyle(1)
    i.SetLineWidth(2)
    i.SetMarkerStyle(0)
    i.SetLineColor(i.GetMarkerColor())
    if i.GetMaximum()>maxValue:
        maxValue = i.GetMaximum()

  slicesBak[0].GetYaxis().SetRangeUser(0,maxValue*1.1)
  slicesBak[0].GetYaxis().SetTitle("Normalized Counts")
  slicesBak[0].Draw("hist L")
  for i in slicesBak:
    i.Draw("hist same L")
  leg.Draw()
  saveAs(canvas,outDir+"slices_bak_"+hist)


  ################################################
  ################################################
  ################################################
  ################################################

  slicesSig = []
  slicesBak = []
  slicesData = []
  maxValue = 0.0
  legendEntries = []
  leg = root.TLegend(0.7,0.65,0.89,0.89)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  print hSig.GetNbinsY()
  for i in range(10,hSig.GetNbinsY()-5):
    legendEntry = "Y #in [{:.2f},{:.2f}]".format(hSig.GetYaxis().GetBinLowEdge(i),hSig.GetYaxis().GetBinUpEdge(i))

    tmp = getYBinHist( hSig, i)
    tmp.SetMarkerColor(colors[i])
    tmp.SetMarkerStyle(33)
    tmp.SetLineWidth(0)
    tmp.SetMarkerSize(2.5)
    if tmp.GetMaximum()>maxValue:
        maxValue = tmp.GetMaximum()
    slicesSig.append(tmp)
    leg.AddEntry(tmp,legendEntry,"p")
     
    tmp = getYBinHist( hBak, i)
    tmp.SetMarkerColor(colors[i])
    tmp.SetMarkerStyle(22)
    tmp.SetLineWidth(0)
    tmp.SetMarkerSize(2.5)
    if tmp.GetMaximum()>maxValue:
        maxValue = tmp.GetMaximum()
    slicesBak.append(tmp)

    tmp = getYBinHist( hData, i)
    tmp.SetMarkerColor(colors[i])
    tmp.SetMarkerStyle(20)
    tmp.SetMarkerSize(2.5)
    tmp.SetLineWidth(0)
    if tmp.GetMaximum()>maxValue:
        maxValue = tmp.GetMaximum()
    slicesData.append(tmp)

  canvas = root.TCanvas()
  canvas.SetLogy(1)
  slicesBak[0].GetYaxis().SetRangeUser(1e-6,maxValue)
  slicesBak[0].GetXaxis().SetTitle(xlabels[hist])
  slicesBak[0].GetYaxis().SetTitle("Counts")
  slicesBak[0].Draw("")
  for i in slicesBak:
    i.Draw("same")
  for i in slicesSig:
    i.Draw("same")
  #for i in slicesData:
  #  i.Draw("same")

  leg.Draw()

  saveAs(canvas,outDir+"slices_y_"+hist)


#class PlotOfSlices:
#  def __init__(self, hist2D, xtitle, ytitle, canvas, xlimits=[], ylimits=[],sliceLabelPrefix="",isPreliminary=True,is7TeV=False):

#  PlotOfSlices(hSig,"mass","ytitle",canvas)
#  saveAs(canvas,outDir+"test_slices_y_"+hist)

