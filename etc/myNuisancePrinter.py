#!/usr/bin/env python

import sys
import ROOT
root = ROOT
gStyle = root.gStyle

ROOT.gROOT.SetBatch(True)

if sys.argv < 2:
  print("Error: Requires filename argument")
  sys.exit(1)


import os,sys,inspect

def ensureImportInCurrentPath():
  currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
  parentdir = os.path.dirname(currentdir)
  sys.path.insert(0,currentdir) 
  sys.path.insert(0,parentdir) 

ensureImportInCurrentPath()
import helpers


file = ROOT.TFile(sys.argv[1])
if file == None: raise RuntimeError, "Cannot open file %s" % "mlfit.root"
fit_s  = file.Get("fit_s")
fit_b  = file.Get("fit_b")
prefit = file.Get("nuisances_prefit")
if fit_s == None or fit_s.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the signal fit 'fit_s'"     % args[0]
if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % args[0]
if prefit == None or prefit.ClassName() != "RooArgSet":    raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'"  % args[0]

graphAllGaus = root.TGraphErrors()
graphCMSCommon = root.TGraphErrors()
allGausNames = []
allCMSCommonNames = []
iAll = 0
iCommon = 0

table = {}
fpf_b = fit_b.floatParsFinal()
fpf_s = fit_s.floatParsFinal()
for i in range(fpf_s.getSize()):
    nuis_s = fpf_s.at(i)
    name   = nuis_s.GetName();
    nuis_b = fpf_b.find(name)
    nuis_p = prefit.find(name)
    row = []
    flag = False;
    # .getVal() .getError() .getMin() .getMax()

    rowString = ""
    rowString += "{0:30}".format(name)
    rowString += "{0:8.3f} +/- {1:5.3f}  ".format(nuis_s.getVal(),nuis_s.getError())

    rowString += "{0:15}".format("[{0:.1f},{1:.1f}]".format(nuis_s.getMin(),nuis_s.getMax()))

    if nuis_p:
      rowString += " Gaus Constraint: {0:8.3f} +/- {1:5.3f}".format(nuis_p.getVal(),nuis_p.getError())

    print(rowString)

    if nuis_p:
      if nuis_p.getError() == 0.:
        #nuis_p.Print()
        #nuis_s.Print()
        pass
      elif not "Mean" in name and not "Width" in name:
        pullVal = (nuis_s.getVal()-nuis_p.getVal())/nuis_p.getError()
        pullErr = nuis_s.getError()/nuis_p.getError()
        graphAllGaus.SetPoint(iAll,pullVal,iAll)
        graphAllGaus.SetPointError(iAll,pullErr,0.)
        allGausNames.append(name)
        iAll += 1
        if "CMS_res" in name or "CMS_scale" in name or "QCDscale" in name or "lumi" in name or name == "pdf_gg" or name == "pdf_qqbar":
          graphCMSCommon.SetPoint(iCommon,pullVal,iCommon)
          graphCMSCommon.SetPointError(iCommon,pullErr,0.)
          allCMSCommonNames.append(name)
          iCommon += 1

## Plot all nuisances
#root.gStyle.SetPadLeftMargin(2.5*root.gStyle.GetPadLeftMargin())
root.gStyle.SetPadLeftMargin(0.45)
canvas = root.TCanvas()
nPoints = graphCMSCommon.GetN()
axisHist = root.TH2F("axisHistCommon","",1,-1.5,1.5,nPoints,-0.5,nPoints-0.5)
axisHist.GetXaxis().SetTitle("Nuisance Parameter Pull")
yaxis = axisHist.GetYaxis()
for iBin, label in enumerate(allCMSCommonNames):
  yaxis.SetBinLabel(iBin+1,label)
yaxis.SetLabelSize(yaxis.GetLabelSize()*1.5)
axisHist.Draw()
zeroGraph = root.TGraph()
zeroGraph.SetPoint(0,0.,-0.5)
zeroGraph.SetPoint(1,0.,nPoints-0.5)
zeroGraph.SetLineStyle(2)
zeroGraph.Draw("L")
#bandGraph = root.TGraphErrors()
#bandGraph.SetPoint(0,0.,0.001)
#bandGraph.SetPointError(0,1.,0.)
#bandGraph.SetPoint(1,0.,nPoints-0.001)
#bandGraph.SetPointError(1,1.,0.)
#bandGraph.SetFillColor(root.kGreen)
#bandGraph.SetFillStyle(1)
#bandGraph.Draw("3")
bandBox = root.TBox(-1.,-0.5+0.0001,1.0,nPoints-0.5001)
bandBox.SetLineColor(root.kGreen)
bandBox.SetFillColor(root.kGreen)
bandBox.Draw()
graphCMSCommon.Draw("PE")
tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.04)
tlatex.SetTextAlign(12)
tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
tlatex.SetTextAlign(32)
tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,"H#rightarrow#mu#mu")
canvas.RedrawAxis()
canvas.SaveAs("CMSCommonNuisances.png")
canvas.SaveAs("CMSCommonNuisances.pdf")

## Plot all nuisances
#root.gStyle.SetPadLeftMargin(2.5*root.gStyle.GetPadLeftMargin())
root.gStyle.SetPadLeftMargin(0.45)
nPoints = graphAllGaus.GetN()
axisHist = root.TH2F("axisHistAll","",1,-1.5,1.5,nPoints,-0.5,nPoints-0.5)
axisHist.GetXaxis().SetTitle("Nuisance Parameter Pull")
yaxis = axisHist.GetYaxis()
for iBin, label in enumerate(allGausNames):
  yaxis.SetBinLabel(iBin+1,label)
yaxis.SetLabelSize(yaxis.GetLabelSize()*0.9)
axisHist.Draw()
zeroGraph = root.TGraph()
zeroGraph.SetPoint(0,0.,-0.5)
zeroGraph.SetPoint(1,0.,nPoints-0.5)
zeroGraph.SetLineStyle(2)
zeroGraph.Draw("L")
bandBox = root.TBox(-1.,-0.5+0.0001,1.0,nPoints-0.5001)
bandBox.SetLineColor(root.kGreen)
bandBox.SetFillColor(root.kGreen)
bandBox.Draw()
graphAllGaus.Draw("PE")
tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(root.gStyle.GetLabelFont())
tlatex.SetTextSize(0.04)
tlatex.SetTextAlign(12)
tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
tlatex.SetTextAlign(32)
tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,"H#rightarrow#mu#mu")
canvas.RedrawAxis()
canvas.SaveAs("AllNuisances.png")
canvas.SaveAs("AllNuisances.pdf")
