#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os

infilename = "matches.root"
outdir = "output/"
CUTS="1"

nBinsDiff=20
minDiff=-5.
maxDiff=5.

RUNPERIOD="8TeV"
LUMI=lumiDict[RUNPERIOD]
LEGDRAWSTRING="lep"

urLegendPos = [0.68,0.65,0.88,0.88]
lcLegendPos = [0.45,0.25,0.65,0.48]
lrLegendPos = [0.68,0.25,0.88,0.48]
ulLegendPos = [0.25,0.65,0.45,0.88]
stdLegendPos = urLegendPos

def getHist(tree,var,cuts,name,nbins=50,minX=110,maxX=160):
   drawStr = var+" >> "+name+"("+str(nbins)+","+str(minX)+","+str(maxX)+")"
   tree.Draw(drawStr,cuts)
   tmp = root.gDirectory.Get(name)
   if type(tmp) != root.TH1F:
     print("Warning: loading histogram: %s from var %s: Object type is not TH1F!!" % (name,var))
     print("  Draw string: '%s'" % drawStr)
     print("  cut string: '%s'" % cuts)
   tmp.UseCurrentStyle()
   tmp.Sumw2()
   tmp.SetTitle("")
   return tmp
def getHist2D(tree,varx,vary,cuts,name,nbinsY,minY,maxY,nbinsX=50,minX=110,maxX=160):
   drawStr = varx+":"+vary+" >> "+name+"("+str(nbinsX)+","+str(minX)+","+str(maxX)+","+str(nbinsY)+","+str(minY)+","+str(maxY)+")"
   tree.Draw(drawStr,cuts)
   tmp = root.gDirectory.Get(name)
   if type(tmp) != root.TH2F:
     print("Warning: loading histogram: %s from varx %s, vary %s: Object type is not TH2F!!" % (name,varx,vary))
     print("  Draw string: '%s'" % drawStr)
     print("  cut string: '%s'" % cuts)
   tmp.UseCurrentStyle()
   tmp.Sumw2()
   tmp.SetTitle("")
   return tmp
def setHistColor(hist,col):
  hist.SetLineColor(col)
  hist.SetMarkerColor(col)
  hist.SetFillColor(col)
def newLeg():
  leg = root.TLegend(*stdLegendPos)
  leg.SetLineColor(0)
  leg.SetFillColor(0)
  return leg

#######################################
root.gROOT.SetBatch(True)
setStyle()
#######################################

canvas = root.TCanvas("canvas")
tlatex = root.TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(gStyle.GetLabelFont("X"))
tlatex.SetTextSize(gStyle.GetLabelSize("X"))
tlatex.SetTextAlign(13)
latexSize = tlatex.GetTextSize()

f = root.TFile(infilename)
matchTree = f.Get("match")
only1Tree = f.Get("only1")
only2Tree = f.Get("only2")
all1Tree = f.Get("all1")
all2Tree = f.Get("all2")

tf1 = root.TF1("tf1","gaus",-5,5)

for case in ["","BB"]:
  title = "Inclusive"
  caseCutStr=""
  caseCutStr1=""
  caseCutStr2=""
  if case =="":
    pass
  else:
    title += " "+case
  if case == "BB":
    caseCutStr1=" && ((16 & eventType1) > 0)"
    caseCutStr2=" && ((16 & eventType2) > 0)"
    caseCutStr=" && ((16 & eventType1) > 0) && ((16 & eventType2) > 0)"

  print CUTS+caseCutStr

  name = "compareMass"+case
  all1Hist = getHist(all1Tree,"dimuonMass1",CUTS+caseCutStr1,name+"1")
  all2Hist = getHist(all2Tree,"dimuonMass_noMuscle2",CUTS+caseCutStr2,name+"2")
  setHistTitles(all1Hist,"m_{#mu#mu} GeV/c^{2}","Events / 1 GeV/c^{2}")
  setHistTitles(all2Hist,"m_{#mu#mu} GeV/c^{2}","Events / 1 GeV/c^{2}")
  setHistColor(all1Hist,root.kBlue)
  setHistColor(all2Hist,root.kRed)
  leg = newLeg()
  leg.AddEntry(all1Hist,"Prompt","lep")
  leg.AddEntry(all2Hist,"22 Jan ReReco","lep")
  all1Hist.Draw()
  all2Hist.Draw("same")
  leg.Draw()
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  saveAs(canvas,outdir+name)

  #canvas.SetLogy(1)
  name = "compareMassMore"+case
  match1Hist = getHist(matchTree,"dimuonMass1",CUTS+caseCutStr1,"matchh"+"1")
  match2Hist = getHist(matchTree,"dimuonMass_noMuscle2",CUTS+caseCutStr2,"matchh"+"2")
  only1Hist = getHist(only1Tree,"dimuonMass1",CUTS+caseCutStr1,"onlyh"+"1")
  only2Hist = getHist(only2Tree,"dimuonMass_noMuscle2",CUTS+caseCutStr2,"onlyh"+"2")
  setHistColor(match1Hist,root.kBlue+1)
  setHistColor(match2Hist,root.kRed+1)
  setHistColor(only1Hist,root.kBlue+6)
  setHistColor(only2Hist,root.kRed+6)
  leg.AddEntry(match1Hist,"Prompt Match","lep")
  leg.AddEntry(match2Hist,"ReReco Match","lep")
  leg.AddEntry(only1Hist,"Prompt Only","lep")
  leg.AddEntry(only2Hist,"ReReco Only","lep")
  all1Hist.GetXaxis().SetRangeUser(120,130)
  all1Hist.GetYaxis().SetRangeUser(0.,all1Hist.GetMaximum()*1.05)
  all1Hist.Draw()
  all2Hist.Draw("same")
  match1Hist.Draw("same")
  match2Hist.Draw("same")
  only1Hist.Draw("same")
  only2Hist.Draw("same")
  leg.Draw()
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  saveAs(canvas,outdir+name)
  canvas.SetLogy(0)
  
  name = "massDiffAll"+case
  massDiffAllHist = getHist(matchTree,"(dimuonMass_noMuscle2-dimuonMass1)/dimuonMass_noMuscle2*100.",CUTS+caseCutStr+" && dimuonMass1>110. && dimuonMass1<160. && dimuonMass_noMuscle2>110. && dimuonMass_noMuscle2<160.",name,40,-5,5)
  setHistTitles(massDiffAllHist,"(m_{#mu#mu}(ReReco)-m_{#mu#mu}(Prompt MuScle))/m_{#mu#mu}(ReReco) %","Events / 0.25%")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"110 GeV < m_{#mu#mu} < 160 GeV")
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)
  
  name = "massDiff120130"+case
  massDiffAllHist = getHist(matchTree,"dimuonMass_noMuscle2-dimuonMass1",CUTS+caseCutStr+" && dimuonMass1>120. && dimuonMass1<130. && dimuonMass_noMuscle2>120. && dimuonMass_noMuscle2<130.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(ReReco)-m_{#mu#mu}(Prompt MuScle) GeV/c^{2}","Events / 0.5 GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"120 GeV < m_{#mu#mu} < 130 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)
  
  name = "massDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"dimuonMass_noMuscle2-dimuonMass1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(ReReco)-m_{#mu#mu}(Prompt MuScle) GeV/c^{2}","Events / 0.5 GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  print "Prompt Match 125 126: %.1f" % massDiffAllHist.Integral()
  saveAs(canvas,outdir+name)

  name = "massDiffReReco125126"+case
  massDiffAllHist = getHist(matchTree,"dimuonMass_noMuscle2-dimuonMass1",CUTS+caseCutStr+" && dimuonMass_noMuscle2>=125. && dimuonMass_noMuscle2<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(ReReco)-m_{#mu#mu}(Prompt MuScle) GeV/c^{2}","Events / 0.5 GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{ReReco} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  print "Rereco Match 125 126: %.1f" % massDiffAllHist.Integral()
  saveAs(canvas,outdir+name)

  name = "massDiffPromptNoMuscle125126"+case
  massDiffAllHist = getHist(matchTree,"dimuonMass_noMuscle2-dimuonMass_noMuscle1",CUTS+caseCutStr+" && dimuonMass_noMuscle1>=125. && dimuonMass_noMuscle1<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(ReReco)-m_{#mu#mu}(Prompt Not MuScle) GeV/c^{2}","Events / 0.5 GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  print "Prompt Match 125 126: %.1f" % massDiffAllHist.Integral()
  saveAs(canvas,outdir+name)

  name = "massDiffReRecoNoMuscle125126"+case
  massDiffAllHist = getHist(matchTree,"dimuonMass_noMuscle2-dimuonMass_noMuscle1",CUTS+caseCutStr+" && dimuonMass_noMuscle2>=125. && dimuonMass_noMuscle2<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(ReReco)-m_{#mu#mu}(Prompt Not MuScle) GeV/c^{2}","Events / 0.5 GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{ReReco} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  print "Rereco Match 125 126: %.1f" % massDiffAllHist.Integral()
  saveAs(canvas,outdir+name)

  name = "massDiffPromptOnly125126"+case
  massDiffAllHist = getHist(matchTree,"dimuonMass1-dimuonMass_noMuscle1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(Prompt MuScle)-m_{#mu#mu}(Prompt) GeV/c^{2}","Events / 0.5 GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "massDiffPromptOnly125126"+case
  massDiffAllHist = getHist(matchTree,"dimuonMass1-dimuonMass_noMuscle1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(Prompt MuScle)-m_{#mu#mu}(Prompt) GeV/c^{2}","Events / 0.5 GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)


#################################################
# pt, Eta, and Phi for Prompt 125-126

  name = "ptLeadDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"muonLead_pt_noMuscle2-muonLead_pt1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"Leading Muon p_{T}^{ReReco}-p_{T}^{Prompt MuScle} GeV/c","Events / 0.5 GeV/c")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "ptSubDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"muonSub_pt_noMuscle2-muonSub_pt1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,nBinsDiff,minDiff,maxDiff)
  setHistTitles(massDiffAllHist,"Sub-leading Muon p_{T}^{ReReco}-p_{T}^{Prompt MuScle} GeV/c","Events / 0.5 GeV/c")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "ptLeadRelDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"(muonLead_pt_noMuscle2-muonLead_pt1)/muonLead_pt_noMuscle2",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,20,-0.1,0.1)
  setHistTitles(massDiffAllHist,"Leading Muon (p_{T}^{ReReco}-p_{T}^{Prompt MuScle})/p_{T}^{ReReco}","Events / 0.01")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.3f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "ptSubRelDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"(muonSub_pt_noMuscle2-muonSub_pt1)/muonSub_pt_noMuscle2",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,20,-0.1,0.1)
  setHistTitles(massDiffAllHist,"Sub-leading Muon (p_{T}^{ReReco}-p_{T}^{Prompt MuScle})/p_{T}^{ReReco}","Events / 0.01")
  setHistColor(massDiffAllHist,root.kBlack)
  #massDiffAllHist.GetXaxis().SetTitleSize(massDiffAllHist.GetXaxis().GetTitleSize()*0.7)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.5)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.3f" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "etaLeadDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"muonLead_eta2-muonLead_eta1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,20,-0.001,0.001)
  setHistTitles(massDiffAllHist,"Leading Muon #eta^{ReReco}-#eta^{Prompt}","Events / 0.001")
  setHistColor(massDiffAllHist,root.kBlack)
  #massDiffAllHist.GetXaxis().SetTitleSize(massDiffAllHist.GetXaxis().GetTitleSize()*0.7)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.25)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2e" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "etaSubDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"muonSub_eta2-muonSub_eta1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,20,-0.001,0.001)
  setHistTitles(massDiffAllHist,"Sub-leading Muon #eta^{ReReco}-#eta^{Prompt}","Events / 0.001")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.25)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2e" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "phiLeadDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"muonLead_phi2-muonLead_phi1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,20,-0.001,0.001)
  setHistTitles(massDiffAllHist,"Leading Muon #phi^{ReReco}-#phi^{Prompt}","Events / 10^{-4}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.25)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2e" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

  name = "phiSubDiffPrompt125126"+case
  massDiffAllHist = getHist(matchTree,"muonSub_phi2-muonSub_phi1",CUTS+caseCutStr+" && dimuonMass1>=125. && dimuonMass1<126.",name,20,-0.001,0.001)
  setHistTitles(massDiffAllHist,"Sub-leading Muon #phi^{ReReco}-#phi^{Prompt}","Events / 10^{-4}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(11)
  tlatex.DrawLatex(root.gStyle.GetPadLeftMargin(),1.02-root.gStyle.GetPadTopMargin(),"125 GeV < m_{#mu#mu}^{Prompt} < 126 GeV")
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  tlatex.SetTextAlign(33)
  tlatex.SetTextSize(latexSize*1.25)
  tlatex.DrawLatex(0.97-root.gStyle.GetPadRightMargin(),0.96-root.gStyle.GetPadTopMargin(),"RMS: %.2e" % massDiffAllHist.GetRMS())
  tlatex.SetTextSize(latexSize)
  saveAs(canvas,outdir+name)

#######################################################################################

  name = "massDiffvMPrompt"+case
  massDiffAllHist = getHist2D(matchTree,"dimuonMass_noMuscle2-dimuonMass1","dimuonMass1",CUTS+caseCutStr,name,nBinsDiff,minDiff,maxDiff,10,110,160)
  massDiffAllHist = massDiffAllHist.ProfileX("_pfx",1,-1,'s')
  setHistTitles(massDiffAllHist,"m_{#mu#mu}(Prompt MuScle) GeV/c^{2}","m_{#mu#mu}(ReReco)-m_{#mu#mu}(Prompt MuScle) GeV/c^{2}")
  setHistColor(massDiffAllHist,root.kBlack)
  massDiffAllHist.Draw()
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  saveAs(canvas,outdir+name)

  name = "ptDiffvPhiPrompt125126"+case
  massDiffAllHist = getHist2D(matchTree,"muonLead_pt_noMuscle2-muonLead_pt1","muonLead_phi2",CUTS+caseCutStr,name,100,-5,5,10,-3.2,3.2)
  massDiffAllHist = massDiffAllHist.ProfileX("_pfx",1,-1,'s')
  newHist = root.TH1F("newHist"+name+case,"",10,-3.2,3.2)
  setHistTitles(newHist,"Leading Muon #phi","RMS of Leading Muon (p_{T}^{ReReco}-p_{T}^{Prompt MuScle})/p_{T}^{ReReco}")
  setHistColor(newHist,root.kBlack)
  newHist.SetFillStyle(0)
  for i in range(massDiffAllHist.GetNbinsX()):
    newHist.SetBinContent(i,massDiffAllHist.GetBinError(i))
  newHist.Draw()
  tlatex.SetTextAlign(31)
  tlatex.DrawLatex(1.0-root.gStyle.GetPadRightMargin(),1.02-root.gStyle.GetPadTopMargin(),title)
  saveAs(canvas,outdir+name)

