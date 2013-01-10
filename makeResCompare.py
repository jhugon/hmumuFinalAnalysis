#!/usr/bin/env python

from helpers import *
import math
import os.path
import glob
import numpy
import matplotlib.pyplot as mpl

root.gROOT.SetBatch(True)
root.gStyle.SetOptStat(0)

#from ROOT import gSystem
#gSystem.Load('libRooFit')

from makeCards import makePDFSig

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

class ResCompare:
  def __init__(self,fileNames,titles,categories):
    assert(len(fileNames) == len(titles))
    self.files = []
    self.histsCatDict = {}
    self.titles = titles
    self.nFiles = len(fileNames)
    self.categories = categories

    self.colors = [1,root.kRed+1,root.kBlue+1,root.kGreen+1,root.kCyan,root.kPink]

    for i in categories:
      self.histsCatDict[i] = []

    for f in fileNames:
      tmpFile = root.TFile(f)
      self.files.append(tmpFile)
      for i in categories:
        strToGet = i + '/mDiMu'
        strToGet = os.path.normpath(strToGet)
        if strToGet[0] == '/':
            strToGet = strToGet[1:]
        tmpHist = tmpFile.Get(strToGet)
        self.histsCatDict[i].append(tmpHist)
    
    minMass = 110
    maxMass = 140

    self.workspace = root.RooWorkspace("w")
    wImport = getattr(self.workspace,"import")
    mMuMu = root.RooRealVar("mMuMu","m_{#mu#mu} [GeV]",minMass,maxMass)
    wImport(mMuMu)
    self.mMuMu = mMuMu
    for i in categories:
     for j in range(self.nFiles):
       makePDFSig(i+str(j),self.histsCatDict[i][j],mMuMu,minMass,maxMass,wImport,i)

    #self.workspace.Print()

    self.canvas = root.TCanvas("canvas")
    
  def plotData(self,saveNameBase):
    datasets = []
    for i in self.categories:
      self.canvas.Clear()
      frame = self.mMuMu.frame()
      frame.SetTitle("")
      frame.SetYTitle("Signal MC Events")
      urLegendPos = [0.70,0.67,0.9,0.9]
      leg = root.TLegend(*urLegendPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      for j in range(self.nFiles):
        tmpDataset = self.workspace.data(i+str(j)+"_Template")
        rooLCol = root.RooFit.LineColor(self.colors[j])
        rooMCol = root.RooFit.MarkerColor(self.colors[j])
        rooNameStr = "Curve_"+i+str(j)+"_Template"
        rooName = root.RooFit.Name(rooNameStr)
        tmpDataset.plotOn(frame,rooLCol,rooMCol,rooName)
        tmpDatasetH = frame.getHist(rooNameStr)
        leg.AddEntry(tmpDatasetH,self.titles[j],"p")
        datasets.append(tmpDataset)
      frame.Draw()
      leg.Draw()
      tlatex = root.TLatex()
      tlatex.SetNDC()
      tlatex.SetTextFont(root.gStyle.GetLabelFont())
      tlatex.SetTextSize(0.05)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.33,0.96,"CMS Internal")
      saveAs(self.canvas,saveNameBase+i)

  def plotPDF(self,saveNameBase):
    curves = []
    for i in self.categories:
      self.canvas.Clear()
      frame = self.mMuMu.frame()
      frame.SetTitle("")
      frame.SetYTitle("Signal MC Events")
      urLegendPos = [0.70,0.67,0.9,0.9]
      leg = root.TLegend(*urLegendPos)
      leg.SetFillColor(0)
      leg.SetLineColor(0)
      for j in range(self.nFiles):
        tmpCurve = self.workspace.pdf(i+str(j))
        rooLCol = root.RooFit.LineColor(self.colors[j])
        rooNameStr = "Curve_"+i+str(j)
        rooName = root.RooFit.Name(rooNameStr)
        tmpCurve.plotOn(frame,rooLCol,rooName)
        tmpCurveH = frame.getCurve(rooNameStr)
        leg.AddEntry(tmpCurveH,self.titles[j],"l")
        curves.append(tmpCurve)
      frame.Draw()
      leg.Draw()
      tlatex = root.TLatex()
      tlatex.SetNDC()
      tlatex.SetTextFont(root.gStyle.GetLabelFont())
      tlatex.SetTextSize(0.05)
      tlatex.SetTextAlign(22)
      tlatex.DrawLatex(0.33,0.96,"CMS Internal")

      saveAs(self.canvas,saveNameBase+i)

  def plotDataPDF(self,saveNameBase):
    datasets = []
    for i in self.categories:
      for j in range(self.nFiles):
        self.canvas.Clear()
        frame = self.mMuMu.frame()
        frame.SetTitle("")
        frame.SetYTitle("Signal MC Events")
        tmpDataset = self.workspace.data(i+str(j)+"_Template")
        rooNameStr = "Curve_"+i+str(j)+"_Template"
        rooName = root.RooFit.Name(rooNameStr)
        tmpDataset.plotOn(frame,rooName)
        datasets.append(tmpDataset)

        tmpCurve = self.workspace.pdf(i+str(j))
        rooNameStr = "Curve_"+i+str(j)
        rooName = root.RooFit.Name(rooNameStr)
        tmpCurve.plotOn(frame,rooName)

        frame.Draw()

        chi2ondf =  frame.chiSquare()
        tlatex = root.TLatex()
        tlatex.SetNDC()
        tlatex.SetTextFont(root.gStyle.GetLabelFont())
        tlatex.SetTextSize(0.05)
        tlatex.SetTextAlign(22)
        tlatex.DrawLatex(0.75,0.85,"#chi^2/NDF = {0:.2f}".format(chi2ondf))

        tlatex.DrawLatex(0.33,0.96,"CMS Internal")
        
        saveAs(self.canvas,saveNameBase+i+self.titles[j])

  def printVars(self,saveNameBase,tex=False,varNames=["Mean","Width","Alpha","n"]):
    doNum = True
    for i in self.categories:
      output = " "*10
      if tex:
        output += r" & "
      for v in varNames:
        output += "{0:<10}".format(v)
        if tex:
          output += " & "
      if doNum:
        output += "{0:<14}".format("Events Pass")
        if tex:
          output += " & "
        output += "{0:<14}".format(r"Events $\in [-1\sigma,+1\sigma]$")
        if tex:
          output += " & "
      if tex:
        output = output[:len(output)-3]
        output += r" \\ \hline \hline"
      output += "\n"
      for j in range(self.nFiles):
        output += "{0:<10}".format(self.titles[j])
        if tex:
          output += r" & "
        for v in varNames:
          tmpVar = self.workspace.var(i+"_"+i+str(j)+"_"+v)
          val = tmpVar.getVal()
          output += "{0:<10.2f}".format(val)
          if tex:
            output += " & "
        if doNum:
          tmpDataset = self.workspace.data(i+str(j)+"_Template")
          xName = self.mMuMu.GetName()
          meanMass = self.workspace.var(i+"_"+i+str(j)+"_"+"Mean").getVal()
          widthMass = self.workspace.var(i+"_"+i+str(j)+"_"+"Width").getVal()
          nEvts = tmpDataset.sumEntries("1")
          output += "{0:<14.2f}".format(nEvts)
          nEvts = tmpDataset.sumEntries("{0}<{1} && {0}>{2}".format(xName,meanMass+widthMass,meanMass-widthMass))
          if tex:
            output += " & "
          output += "{0:<14.2f}".format(nEvts)
          if tex:
            output += " & "
        if tex:
          output = output[:len(output)-3]
          output += r" \\ \hline"
        output += "\n"
      if tex:
        columnFormat = "|l|" + "c|"*len(varNames)
        if doNum:
          columnFormat += "c|c|"
        output = r"\begin{tabular}{"+columnFormat+"} \hline \n" + output + r"\end{tabular}"+"\n"
      if saveNameBase != None:
        ext = ".txt"
        if tex:
          ext = ".tex"
        f = open(saveNameBase+i+ext,"w")
        f.write(output)
        f.close()
      else:
        print(i)
        print("")
        print(output)
    if tex and saveNameBase != None:
      f = open(saveNameBase+"Test.tex","w")
      testStr = r"""
\documentclass[12pt,a4paper]{article}
\usepackage{lscape}
\begin{document}
%\tiny
\small
\begin{center}

"""
      for i in self.categories:
        testStr += r"\vspace{2ex}{\large "+i +"} \\vspace{1ex}\n \\\\ \n"
        testStr += r"\input{"+saveNameBase+i+"}\n \\\\ \n"
      testStr += r"""

\end{center}
\end{document}
"""
      f.write(testStr)
      f.close()

if __name__ == "__main__":

  #categories = ["IncPresel","VBFPresel","IncBDTCut","VBFBDTCut"]
  categories = ["IncPresel","VBFPresel"]

  infiles = []
  titles = []
  infiles.append("input/ggHmumu125_8TeV.root")
  infiles.append("input/smearing/ggHmumu125_8TeV.root")
  infiles.append("input/rochester/ggHmumu125_8TeV.root")
  infiles.append("input/muscle/ggHmumu125_8TeV.root")
  titles.append("Default")
  titles.append("Smearing")
  titles.append("Rochester")
  titles.append("MuScle")

  #infiles.append("input/vbfHmumu125_8TeV.root")
  #infiles.append("input/smearing/vbfHmumu125_8TeV.root")
  #infiles.append("input/rochester/vbfHmumu125_8TeV.root")
  #infiles.append("input/muscle/vbfHmumu125_8TeV.root")

  rs = ResCompare(infiles,titles,categories)
  rs.plotData("resData")
  rs.plotPDF("resShape")
  rs.printVars("resTable",True)
  rs.plotDataPDF("resData")
