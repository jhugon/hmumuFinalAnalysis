#!/usr/bin/python

import ROOT as root

root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

dataDir = "input/window4GeV/"
outDir = "output/"

sigFname = "ggHmumu"
bakFname = "DYJetsToLL"

fSig = root.TFile(dataDir+sigFname+".root")
fBak = root.TFile(dataDir+bakFname+".root")

hSig = fSig.Get("yVptDiMu")
hBak = fBak.Get("yVptDiMu")

xlims = [0,200]
xlims = [0,100]
ylims = [0,2.5]
rebin = [1,2]
drawOpt = "colz"

hSig.SetTitle("")
hSig.GetXaxis().SetTitle("Dimuon p_{T} [GeV]")
hSig.GetYaxis().SetTitle("Dimuon |y|")
hSig.Rebin2D(*rebin)
hSig.GetXaxis().SetRangeUser(*xlims)
hSig.GetYaxis().SetRangeUser(*ylims)

hBak.SetTitle("")
hBak.GetXaxis().SetTitle("Dimuon p_{T} [GeV]")
hBak.GetYaxis().SetTitle("Dimuon |y|")
hBak.Rebin2D(*rebin)
hBak.GetXaxis().SetRangeUser(*xlims)
hBak.GetYaxis().SetRangeUser(*ylims)

canvas = root.TCanvas()

hSig.Draw(drawOpt)
canvas.SaveAs(outDir+"ptVySig.png")
canvas.SaveAs(outDir+"ptVySig.pdf")

hBak.Draw(drawOpt)
canvas.SaveAs(outDir+"ptVyBak.png")
canvas.SaveAs(outDir+"ptVyBak.pdf")
