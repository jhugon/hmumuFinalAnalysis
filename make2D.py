#!/usr/bin/env python

from xsec import *

root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

dataDir = "input/"
outDir = "output/"

histPrefix = ""

histlist = [
    "yVptDiMu",
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

def CheckForInfSoB(sig,bak):
  nBinsX = sig.GetNbinsX()
  nBinsY = sig.GetNbinsY()
  assert(nBinsX == bak.GetNbinsX())
  assert(nBinsY == bak.GetNbinsY())

  for i in range(nBinsX+2):
    for j in range(nBinsY+2):
      s = sig.GetBinContent(i,j)
      b = bak.GetBinContent(i,j)
      if s>0.0 and b==0:
        print("SoB is inf for hist: {}, s: {}, b: {} x: {} y: {}".format(
         sig.GetName(),
         s,
         b,
         sig.GetXaxis().GetBinCenter(i),
         sig.GetYaxis().GetBinCenter(j)
        ))
  

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
  
  hSig.SetTitle("Signal")
  hSig.GetXaxis().SetTitle(xlabels[hist])
  hSig.GetYaxis().SetTitle(ylabels[hist])
  hSig.Rebin2D(*(rebins[hist]))
  hSig.GetXaxis().SetRangeUser(*(xranges[hist]))
  hSig.GetYaxis().SetRangeUser(*(yranges[hist]))

  hBak.SetTitle("Background")
  hBak.GetXaxis().SetTitle(xlabels[hist])
  hBak.GetYaxis().SetTitle(ylabels[hist])
  hBak.Rebin2D(*(rebins[hist]))
  hBak.GetXaxis().SetRangeUser(*(xranges[hist]))
  hBak.GetYaxis().SetRangeUser(*(yranges[hist]))
  
  hData.SetTitle("2012 Data {} fb^{{-1}}".format(LUMI))
  hData.GetXaxis().SetTitle(xlabels[hist])
  hData.GetYaxis().SetTitle(ylabels[hist])
  hData.Rebin2D(*(rebins[hist]))
  hData.GetXaxis().SetRangeUser(*(xranges[hist]))
  hData.GetYaxis().SetRangeUser(*(yranges[hist]))

  #CheckForInfSoB(hSig,hBak)

  hSoB = hSig.Clone(hist+"_SoB")
  hSoB.Divide(hBak)
  
  canvas = root.TCanvas()
  
  hSig.Draw(drawOpt)
  canvas.SaveAs(outDir+hist+"Sig.png")
  canvas.SaveAs(outDir+hist+"Sig.pdf")
  
  hBak.Draw(drawOpt)
  canvas.SaveAs(outDir+hist+"Bak.png")
  canvas.SaveAs(outDir+hist+"Bak.pdf")

  hData.Draw(drawOpt)
  canvas.SaveAs(outDir+hist+"Data.png")
  canvas.SaveAs(outDir+hist+"Data.pdf")

  hSoB.Draw(drawOpt)
  canvas.SaveAs(outDir+hist+"SoB.png")
  canvas.SaveAs(outDir+hist+"SoB.pdf")

  leg = root.TLegend(0.7,0.7,0.88,0.88)
  leg.SetFillColor(0)
  leg.SetLineColor(0)
  leg.AddEntry(hSig,"Signal","l")
  leg.AddEntry(hBak,"Background","l")

  hBak.SetTitle("")
  hBak.SetLineColor(root.kRed+1)
  hBak.SetFillColor(root.kRed+1)
  hSig.SetLineColor(root.kGreen+1)
  hSig.SetFillColor(root.kGreen+1)
  hBak.Draw("cont3")
  hSig.Draw("cont3 same")
  leg.Draw()
  canvas.SaveAs(outDir+hist+"Cont.png")
  canvas.SaveAs(outDir+hist+"Cont.pdf")

  leg.Clear()
  leg.AddEntry(hSig,"Signal","f")
  leg.AddEntry(hBak,"Background","f")
  
  hSig.Scale(hBak.Integral()/hSig.Integral())
  hBak.Draw("box")
  hSig.Draw("box same")
  leg.Draw()
  canvas.SaveAs(outDir+hist+"Box.png")
  canvas.SaveAs(outDir+hist+"Box.pdf")
  
  hBak.Draw("box")
  hSig.Draw("box same")
  leg.Draw()
  canvas.SaveAs(outDir+hist+"Box.png")
  canvas.SaveAs(outDir+hist+"Box.pdf")
  
  hBak.Draw("lego")
  hSig.Draw("lego same")
  leg.Draw()
  canvas.SaveAs(outDir+hist+"lego.png")
  canvas.SaveAs(outDir+hist+"lego.pdf")
  
