#! /usr/bin/env python

from ROOT import gSystem

import datetime
import sys
import os
import re
import math
import cPickle
import ROOT as root
root.gSystem.Load('libRooFit')
root.gROOT.SetBatch(True)
import scipy.stats

from helpers import *
import makeCards

from numpy import mean, median, corrcoef, percentile
from numpy import std as stddev

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

canvas = root.TCanvas()

class BiasStudy:
  #def __init__(self,category,dataFiles,energyStr):
  def __init__(self,category,dataFileNames,energyStr,nToys=10,pklOutFnBase="output/biasData",inputPkl=None):
    self.dataFileNames = dataFileNames
    self.sigMasses = range(115,156,5)
    self.sigMasses = [125,150]
    ## Try to load data from pkl file
    if inputPkl != None:
      try:
        inputPklF = open(inputPkl)
        self.data = cPickle.load(inputPklF)
        inputPklF.close()
        self.refPdfNameList = self.data['meta']['refPdfNameList']
        self.pdfAltNamesDict = self.data['meta']['pdfAltNamesDict']
        self.nToys = self.data['meta']['nToys']
        self.nData = self.data['meta']['nData']
        self.catName = self.data['meta']['catName']
        self.energyStr = self.data['meta']['energyStr']
        energyStr = self.energyStr
        self.sigMasses = self.data['meta']['sigMasses']
      except Exception, err:
        print("Error loading data from pkl file: "+str(inputPkl))
        print(err)
        self.data = None
    else:
      self.catName = category[0]
      self.catCuts = category[1]
      self.energyStr = energyStr
      self.nToys = nToys
      self.pklOutFn = pklOutFnBase+"_"+self.catName+"_"+energyStr+".pkl"
      self.data = None
    catName = self.catName
    canvas = root.TCanvas()
    self.canvas = canvas
    ## Run
    if self.data==None:
      self.refPdfNameList = [
      #    "ExpLog",
      #    "MOverSq",
          "Old",
          "ExpMOverSq",
      ]
      self.pdfAltNamesDict = {
          "ExpLog":["ExpMOverSq"],
          "MOverSq":["ExpMOverSq"],
          "Old":["ExpMOverSq"],
          "ExpMOverSq":[          
          #                  "ExpLog",
          #                  "MOverSq",
                            "Old",
                        ],
      }

      self.data = {'meta':{}}
      data = self.data
      data['meta']['sigMasses'] = self.sigMasses
      data['meta']['refPdfNameList'] = self.refPdfNameList
      data['meta']['pdfAltNamesDict'] = self.pdfAltNamesDict
      data['meta']['nToys'] = self.nToys
      data['meta']['catName'] = self.catName
      data['meta']['energyStr'] = self.energyStr
      self.dataTree = root.TChain()
      self.iPklAutoSave = 1
      for i in self.dataFileNames:
        self.dataTree.Add(i+"/outtree"+self.catName)
      self.dataTree.SetCacheSize(10000000);
      self.dataTree.AddBranchToCache("*");
      for refPdfName in self.refPdfNameList:
        self.runStudy(refPdfName)
      pklFile = open(self.pklOutFn,'w')
      cPickle.dump(data,pklFile)
      pklFile.close()

    outStr = "#"*80+"\n"
    outStr = "#"*80+"\n"
    outStr += "\n"+self.catName +"  "+self.energyStr + "\n\n"

    outStr += "nToys: {0}\n".format(self.nToys)
    outStr += "nData Events: {0}\n".format(self.nData)
    outStr += "\n"

    data = self.data

    for refPdfName in self.refPdfNameList:
      outStr +=  "Reference PDF: "+str(refPdfName)+'\n'
      for hmass in self.sigMasses:
        outStr +=  "mass: "+str(hmass)+'\n'
        dataH = data[refPdfName][hmass]
        shapiroStat, shapiroP = scipy.stats.shapiro(dataH['zTrue'])
        sumChi2 = sum(dataH['chi2True'])
        sumNDF = sum(dataH['ndfTrue'])
        outStr +=  "  True Z Scores:   {0:.2f} +/- {1:.2f}  Median: {2:.2f}    S-W Normal p-Val: {3:.3g}\n".format(mean(dataH['zTrue']),stddev(dataH['zTrue']),median(dataH['zTrue']),shapiroP)
        outStr +=  "  True Fit Prob:   {0:.3g},                              chi2: {1:.2f}  NDF: {2}\n".format(scipy.stats.chi2.sf(sumChi2,sumNDF),sumChi2,sumNDF)
        outStr +=  "  All Pulls:       {0:.2f} +/- {1:.2f}  Median: {2:.2f}\n".format(mean(dataH['pullAll']),stddev(dataH['pullAll']),median(dataH['pullAll']))
        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          dataHA = dataH[pdfAltName]
          shapiroStat, shapiroP = scipy.stats.shapiro(dataHA['z'])
          sumChi2 = sum(dataHA['chi2'])
          sumNDF = sum(dataHA['ndf'])
          outStr +=  "  "+pdfAltName+":\n"
          outStr +=  "    Z Scores:      {0:.2f} +/- {1:.2f}  Median: {2:.2f}    S-W Normal p-Val: {3:.3g}\n".format(mean(dataHA['z']),stddev(dataHA['z']),median(dataHA['z']),shapiroP)
          outStr +=  "    Fit Prob:      {0:.3g},                              chi2: {1:.2f}  NDF: {2}\n".format(scipy.stats.chi2.sf(sumChi2,sumNDF),sumChi2,sumNDF)
          outStr +=  "    Pulls:         {0:.2f} +/- {1:.2f}  Median: {2:.2f}\n".format(mean(dataHA['pull']),stddev(dataHA['pull']),median(dataHA['pull']))
    print outStr
    self.outStr = outStr

  def runStudy(self,truePdfName):
      data = self.data
      catName = self.catName
      energyStr = self.energyStr
      truePdfFunc = getattr(makeCards,"makePDFBak"+truePdfName)
      pdfAltNameList = self.pdfAltNamesDict[truePdfName]
      pdfAltFuncList = [ getattr(makeCards,"makePDFBak"+i) for i in pdfAltNameList]

      dimuonMass = root.RooRealVar("dimuonMass","m [GeV/c^{2}]",110.,170.)
      dimuonMass.setRange("low",110,120) # Silly ranges for old fit functionality
      dimuonMass.setRange("high",130,170)
      dimuonMass.setRange("signal",120,130)
      dimuonMass.setRange("signalfit",110,140)
      wTrue = root.RooWorkspace("wTrue")
      wTrueToy = root.RooWorkspace("wTrueToy")
      wTrueImport = getattr(wTrue,"import")
      wTrueToyImport = getattr(wTrueToy,"import")

      # Hack to Make makePDFBakOld work
      minMassZ = 88.
      maxMassZ = 94.
      dimuonMassZ = root.RooRealVar("dimuonMass","dimuonMass",minMassZ,maxMassZ)
      dimuonMassZ = dimuonMassZ

      ### Load data
      
      realData = root.RooDataSet("realData"+catName+energyStr,
                                      "realData"+catName+energyStr,
                                          self.dataTree,root.RooArgSet(dimuonMass)
                                        )
      nData = realData.sumEntries()
      self.nData = nData
      data['meta']['nData'] = self.nData
      realDataZ = root.RooDataSet("realDataZ"+catName+energyStr,
                                      "realDataZ"+catName+energyStr,
                                          self.dataTree,root.RooArgSet(dimuonMassZ)
                                        )

      ### Make Bak Pdfs

      truePdfFunc(truePdfName+catName+energyStr,realData,dimuonMass,110,170,wTrueImport,dimuonMassZ,realDataZ)
      truePdf = wTrue.pdf("bak")
      truePdf.SetName(truePdfName)
      truePdf.SetTitle("True PDF ")
      truePdf.fitTo(realData,
                             PRINTLEVEL
                           )

      trueToyPdfName = "trueToy"+catName+energyStr
      truePdfFunc(trueToyPdfName,realData,dimuonMass,110,170,wTrueToyImport,dimuonMassZ,realDataZ)
      trueToyPdf = wTrueToy.pdf("bak")
      trueToyPdf.SetName(trueToyPdfName)

      # Make sure Voigt params are set to True vals and constant
      if truePdfName == "Old":
        for xTrue in rooArgSet2List(truePdf.getParameters(realData)):
          if not ("voit" in xTrue.GetName()):
            continue
          for xToy in rooArgSet2List(trueToyPdf.getParameters(realData)):
            trueMatch = re.match(r".*(_voit.*)",xTrue.GetName()) 
            toyMatch = re.match(r".*(_voit.*)",xToy.GetName()) 
            assert(trueMatch)
            if not toyMatch:
                continue
            trueBaseName = trueMatch.group(1)
            toyBaseName = toyMatch.group(1)
            if not ( trueBaseName == toyBaseName ):
              continue
            xToy.setVal(xTrue.getVal())
            xTrue.setConstant(True)
            xToy.setConstant(True)

      pdfAltList = []
      pdfAltwList = []
      for pdfAltName,pdfAltFunc in zip(pdfAltNameList,pdfAltFuncList):
        pdfName = "alt"+catName+energyStr+"_"+pdfAltName
        wAlt = root.RooWorkspace("wAlt"+catName+energyStr+"_"+pdfAltName)
        pdfAltFunc(pdfName,realData,dimuonMass,110,170,getattr(wAlt,"import"),dimuonMassZ,realDataZ)
        altPdf = wAlt.pdf("bak")
        altPdf.SetName(pdfName)
        # Make sure Voigt params are constant
        if pdfAltName == "Old":
          for x in rooArgSet2List(altPdf.getParameters(realData)):
            if "voit" in x.GetName():
              x.setConstant(True)
        pdfAltList.append(altPdf)
        pdfAltwList.append(wAlt)

      nBakVar = root.RooRealVar("nBak","N_{B}",nData/2.,nData*2)

      ### Now load Signal PDFs
      sigMasses = self.sigMasses
      sigPdfs = []
      wSigs = []
      for hmass in sigMasses:
        wSig = root.RooWorkspace("signal"+catName+energyStr+str(hmass))
        makeCards.makePDFSigNew(catName+energyStr,"sig_ggH",dimuonMass,float(hmass),
                                getattr(wSig,"import")
                               )
        sigPdf = wSig.pdf("ggH")
        sigPdf.SetName("sigPDF_"+str(hmass)+"_"+catName+energyStr)
        sigPdfs.append(sigPdf)
        wSigs.append(wSig)

      nSigVar = root.RooRealVar("nSig","N_{S}",-nData/4.,nData/4)

      ### Make results data structure and begin log
      data[truePdfName] = {}
      for hmass in sigMasses:
        data[truePdfName][hmass] = {}
        data[truePdfName][hmass]['zTrue'] = []
        data[truePdfName][hmass]['nTrue'] = []
        data[truePdfName][hmass]['chi2BOnly'] = []
        data[truePdfName][hmass]['ndfBOnly'] = []
        data[truePdfName][hmass]['chi2True'] = []
        data[truePdfName][hmass]['ndfTrue'] = []
        data[truePdfName][hmass]['errTrue'] = []
        data[truePdfName][hmass]['pullAll'] = []
        for pdfAltName in pdfAltNameList:
          data[truePdfName][hmass][pdfAltName] = {'z':[],'pull':[],'n':[],'err':[],'chi2':[],'ndf':[]}

      ### Toy Loop

      for iToy in range(self.nToys):
        toyData = truePdf.generate(root.RooArgSet(dimuonMass),int(nData))
        toyData.SetName("toyData"+catName+energyStr+str(iToy))
        toyDataHist = toyData.binnedClone("toyDataHist"+catName+energyStr+str(iToy))
        plotThisToy = (iToy % 10 == 0)
        saveThisData = (iToy % 100 == 99)
        for hmass,sigPdf in zip(sigMasses,sigPdfs):
          frame = None 
          if plotThisToy:
            frame = dimuonMass.frame()
            toyData.plotOn(frame)
          # Check Chi^2 for ref background only fit
          trueToyPdf.fitTo(toyData,
                             PRINTLEVEL
                           )
          chi2TrueToyVar = trueToyPdf.createChi2(toyDataHist)
          ndfTrue = dimuonMass.getBins() - 1  # b/c roofit normalizes
          ndfTrue -= trueToyPdf.getParameters(toyDataHist).getSize()
          data[truePdfName][hmass]['chi2BOnly'].append(chi2TrueToyVar.getVal())
          data[truePdfName][hmass]['ndfBOnly'].append(ndfTrue)
          # Now Create S+B PDF and fit for real
          trueToySBPdf = root.RooAddPdf("SBTrue"+catName+energyStr+str(hmass)+"_"+str(iToy),"",
                              root.RooArgList(trueToyPdf,sigPdf),
                              root.RooArgList(nBakVar,nSigVar)
                          )
          trueToySBPdf.fitTo(toyData,
                             PRINTLEVEL
                           )
          chi2TrueToyVar = trueToySBPdf.createChi2(toyDataHist)
          ndfTrue = dimuonMass.getBins() - 1  # b/c roofit normalizes
          ndfTrue -= trueToySBPdf.getParameters(toyDataHist).getSize()
          if plotThisToy:
            trueToySBPdf.plotOn(frame,root.RooFit.LineColor(6))
            #trueToySBPdf.plotOn(frame,root.RooFit.Components(root.RooArgSet(trueToyPdf)),root.RooFit.LineColor(1),root.RooFit.LineStyle(2))
          nTrueToy = nSigVar.getVal()
          errTrueToy = nSigVar.getError()
          data[truePdfName][hmass]['nTrue'].append(nTrueToy)
          data[truePdfName][hmass]['errTrue'].append(errTrueToy)
          data[truePdfName][hmass]['chi2True'].append(chi2TrueToyVar.getVal())
          data[truePdfName][hmass]['ndfTrue'].append(ndfTrue)
          if errTrueToy != 0.:
            data[truePdfName][hmass]['zTrue'].append(nTrueToy/errTrueToy)
          for pdfAlt,pdfAltName,color in zip(pdfAltList,pdfAltNameList,range(2,len(pdfAltList)+2)):
              altSBPdf = root.RooAddPdf("SB"+pdfAltName+catName+energyStr+str(hmass)+"_"+str(iToy),"",
                              root.RooArgList(pdfAlt,sigPdf),
                              root.RooArgList(nBakVar,nSigVar)
                          )
              altSBPdf.fitTo(toyData,
                              PRINTLEVEL
                              )
              altChi2Var = altSBPdf.createChi2(toyDataHist)
              ndfAlt = dimuonMass.getBins() - 1  # b/c roofit normalizes
              ndfAlt -= altSBPdf.getParameters(toyDataHist).getSize()
              if plotThisToy:
                altSBPdf.plotOn(frame,root.RooFit.LineColor(color))
                #altSBPdf.plotOn(frame,root.RooFit.Components(root.RooArgSet(pdfAlt)),root.RooFit.LineColor(color),root.RooFit.LineStyle(2))
              nAlt = nSigVar.getVal()
              errAlt = nSigVar.getError()
              if errAlt == 0.:
                  pull = -1e9
              else:
                  pull = (nAlt-nTrueToy)/errAlt
              data[truePdfName][hmass][pdfAltName]['n'].append(nAlt)
              data[truePdfName][hmass][pdfAltName]['err'].append(errAlt)
              data[truePdfName][hmass][pdfAltName]['chi2'].append(altChi2Var.getVal())
              data[truePdfName][hmass][pdfAltName]['ndf'].append(ndfAlt)
              data[truePdfName][hmass][pdfAltName]['z'].append(nAlt/errAlt)
              data[truePdfName][hmass][pdfAltName]['pull'].append(pull)
              data[truePdfName][hmass]['pullAll'].append(pull)
          if plotThisToy:
            frame.Draw()
            canvas.SaveAs("output/debug_"+catName+"_"+energyStr+"_"+str(hmass)+"_"+str(iToy)+".png")
          if saveThisData:
            pklTmpFile = open(self.pklOutFn+"."+str(self.iPklAutoSave),'w')
            self.iPklAutoSave += 1
            cPickle.dump(data,pklTmpFile)
            pklTmpFile.close()
    
          #if truePdfName == "Old":
          #  print "**************************************************************"
          #  print "True PDF Parameters:"
          #  for x in rooArgSet2List(truePdf.getParameters(realData)):
          #      x.Print()
          #  print "True Toy PDF Parameters:"
          #  for x in rooArgSet2List(trueToyPdf.getParameters(realData)):
          #      x.Print()
          #  print "**************************************************************"
        del toyData
        del toyDataHist

  def plot(self,outputPrefix):
    self.pdfTitleMap = {
        "ExpLog":"Exp(p_{1}m^{2}+p_{2}m+p_{3}ln(m))",
        "MOverSq":"#frac{m}{(m-p_{1})^{2}}",
        "Old":"Voigtian+Exp",
        "ExpMOverSq":"#frac{Exp(p_{1}m)}{(m-p_{2})^{2}}",
    }
    titleMap = {
      "Jets01PassPtG10BB": "0,1-Jet Tight BB",
      "Jets01PassPtG10BO": "0,1-Jet Tight BO",
      "Jets01PassPtG10BE": "0,1-Jet Tight BE",
      "Jets01PassPtG10OO": "0,1-Jet Tight OO",
      "Jets01PassPtG10OE": "0,1-Jet Tight OE",
      "Jets01PassPtG10EE": "0,1-Jet Tight EE",
      "Jets01PassCatAll" : "0,1-Jet Tight Combination",
                            
      "Jets01FailPtG10BB": "0,1-Jet Loose BB",
      "Jets01FailPtG10BO": "0,1-Jet Loose BO",
      "Jets01FailPtG10BE": "0,1-Jet Loose BE",
      "Jets01FailPtG10OO": "0,1-Jet Loose OO",
      "Jets01FailPtG10OE": "0,1-Jet Loose OE",
      "Jets01FailPtG10EE": "0,1-Jet Loose EE",
      "Jets01FailCatAll" : "0,1-Jet Loose Combination",
    
      "Jet2CutsVBFPass":"2-Jet VBF Tight",
      "Jet2CutsGFPass":"2-Jet GF Tight",
      "Jet2CutsFailVBFGF":"2-Jet Loose",
    }

    canvas = self.canvas
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(0.04)
    caption = titleMap[self.catName]
    caption2 = ""
    caption3 = ""
    caption4 = ""
    iHist = 0

    for refPdfName in self.refPdfNameList:
      ##### Pull plots 1D
      for hmass in self.sigMasses:
        if len(self.pdfAltNamesDict[refPdfName])>1:
          hist = root.TH1F("hist"+str(iHist),"",30,-3,3)
          setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","N_{Toys}")
          iHist += 1
          for pull in self.data[refPdfName][hmass]['pullAll']:
              hist.Fill(pull)
          medianPull = median(self.data[refPdfName][hmass]['pullAll'])
          hist.Draw()
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"All Alternate PDFs")
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.60,"Median: {0:.2f}".format(medianPull))
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"{0:.2f}".format(medianPull))
          line = self.setYMaxAndDrawVertLines(hist,medianPull)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_AllPulls_Ref"+refPdfName)
          canvas.Clear()

        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          hist = root.TH1F("hist"+str(iHist),"",30,-3,3)
          setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","N_{Toys}")
          iHist += 1
          for pull in self.data[refPdfName][hmass][pdfAltName]['pull']:
              hist.Fill(pull)
          medianPull = median(self.data[refPdfName][hmass][pdfAltName]['pull'])
          hist.Draw()
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.60,"Median: {0:.2f}".format(medianPull))
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          line = self.setYMaxAndDrawVertLines(hist,medianPull)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Pulls_Ref"+refPdfName+"_Alt"+pdfAltName)
          canvas.Clear()

      ##### Median pull plots v. mass
      for pdfAltName in self.pdfAltNamesDict[refPdfName]:
        minx = 115
        maxx = 155
        axisHist = root.TH2F("axishist"+str(iHist),"",1,minx,maxx,1,-1,1)
        setHistTitles(axisHist,"M_{H} [GeV/c^{2}]","Median[(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)]")
        graph = root.TGraph()
        graphBand = root.TGraphErrors()
        graphBand.SetPoint(0,minx,0.)
        graphBand.SetPointError(0,0.,0.14)
        graphBand.SetPoint(1,maxx,0.)
        graphBand.SetPointError(1,0.,0.14)
        graphBand.SetFillStyle(1)
        graphBand.SetFillColor(root.kGreen-9)
        iHist += 1
        iGraph = 0
        for hmass in self.sigMasses:
          medianPull = median(self.data[refPdfName][hmass][pdfAltName]['pull'])
          graph.SetPoint(iGraph,hmass,medianPull)
          iGraph += 1
        axisHist.Draw()
        graphBand.Draw("3")
        graph.Draw("LP")
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
        tlatex.SetTextAlign(32)
        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
        canvas.RedrawAxis()
        saveAs(canvas,outputPrefix+self.catName+"_Pulls_Ref"+refPdfName+"_Alt"+pdfAltName)
        canvas.Clear()

      ##### Chi2 Prob Plots
      for hmass in self.sigMasses:
        hist = root.TH1F("hist"+str(iHist),"",20,0,1)
        setHistTitles(hist,"#chi^{2} p-Value of Fit","N_{Toys}")
        iHist += 1
        for chi2,ndf in zip(self.data[refPdfName][hmass]['chi2True'],self.data[refPdfName][hmass]['ndfTrue']):
            hist.Fill(scipy.stats.chi2.sf(chi2,ndf))
        hist.Draw()
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"m_{H} = "+str(hmass)+" GeV/c^{2}")
        tlatex.SetTextAlign(32)
        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
        self.setYMaxAndDrawVertLines(hist,None)
        canvas.RedrawAxis()
        saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2_Ref"+refPdfName)
        canvas.Clear()

        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          hist = root.TH1F("hist"+str(iHist),"",20,0,1)
          setHistTitles(hist,"#chi^{2} p-Value of Fit","N_{Toys}")
          iHist += 1
          for chi2,ndf in zip(self.data[refPdfName][hmass][pdfAltName]['chi2'],self.data[refPdfName][hmass][pdfAltName]['ndf']):
            hist.Fill(scipy.stats.chi2.sf(chi2,ndf))
          hist.Draw()
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          self.setYMaxAndDrawVertLines(hist,None)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2_Ref"+refPdfName+"_Alt"+pdfAltName)
          canvas.Clear()

      ##### Chi2vPull
      for hmass in self.sigMasses:
        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          hist = root.TH2F("hist"+str(iHist),"",12,-2,2,5,0,1)
          setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","#chi^{2} p-Value of Fit")
          iHist += 1
          chi2pVals = scipy.stats.chi2.sf(self.data[refPdfName][hmass][pdfAltName]['chi2'],self.data[refPdfName][hmass][pdfAltName]['ndf'])
          for pull,chi2pVal in zip(self.data[refPdfName][hmass][pdfAltName]['pull'],chi2pVals):
            hist.Fill(pull,chi2pVal)
          hist.Draw('col')
          xLine = median(self.data[refPdfName][hmass][pdfAltName]['pull'])
          line = root.TLine()
          line.SetLineColor(root.kBlue)
          line.SetLineWidth(2)
          line.SetLineStyle(2)
          line.DrawLine(xLine,0,xLine,1)
          line.SetLineStyle(1)
          for iY in range(1,hist.GetNbinsY()+1):
            binPullList = []
            binLowVal = hist.GetYaxis().GetBinLowEdge(iY)
            binHighVal = hist.GetYaxis().GetBinUpEdge(iY)
            for iEntry in range(len(self.data[refPdfName][hmass][pdfAltName]['pull'])):
              if chi2pVals[iEntry] >= binLowVal and chi2pVals[iEntry] < binHighVal:
                binPullList.append(self.data[refPdfName][hmass][pdfAltName]['pull'][iEntry])
            binPullMed = median(binPullList)
            line.DrawLine(binPullMed,binLowVal,binPullMed,binHighVal)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2vPull_Ref"+refPdfName+"_Alt"+pdfAltName)
          canvas.Clear()

      ##### PVdeltaNsigAlt
      for hmass in self.sigMasses:
        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          minY = percentile(self.data[refPdfName][hmass][pdfAltName]['err'],5.)
          maxY = percentile(self.data[refPdfName][hmass][pdfAltName]['err'],95.)
          hist = root.TH2F("hist"+str(iHist),"",12,-2,2,30,minY*0.6,maxY*1.5)
          setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","#DeltaN_{sig}(Alt)")
          iHist += 1
          for pull,err in zip(self.data[refPdfName][hmass][pdfAltName]['pull'],self.data[refPdfName][hmass][pdfAltName]['err']):
            hist.Fill(pull,err)
          hist.Draw('col')
          xLine = median(self.data[refPdfName][hmass][pdfAltName]['pull'])
          line = root.TLine()
          line.SetLineColor(root.kBlue)
          line.SetLineWidth(2)
          line.SetLineStyle(2)
          line.DrawLine(xLine,hist.GetYaxis().GetBinLowEdge(1),xLine,hist.GetYaxis().GetBinUpEdge(hist.GetNbinsY()+1))
          line.SetLineStyle(1)
          for iY in range(1,hist.GetNbinsY()+1):
            binPullList = []
            binLowVal = hist.GetYaxis().GetBinLowEdge(iY)
            binHighVal = hist.GetYaxis().GetBinUpEdge(iY)
            for iEntry in range(len(self.data[refPdfName][hmass][pdfAltName]['pull'])):
              if self.data[refPdfName][hmass][pdfAltName]['err'][iEntry] >= binLowVal and self.data[refPdfName][hmass][pdfAltName]['err'][iEntry] < binHighVal:
                binPullList.append(self.data[refPdfName][hmass][pdfAltName]['pull'][iEntry])
            binPullMed = median(binPullList)
            line.DrawLine(binPullMed,binLowVal,binPullMed,binHighVal)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_errNvPull_Ref"+refPdfName+"_Alt"+pdfAltName)
          canvas.Clear()

      ##### Z Plots
      for hmass in self.sigMasses:
        hist = root.TH1F("hist"+str(iHist),"",50,-3,3)
        setHistTitles(hist,"N_{sig}(Ref)/\DeltaN_{sig}(Ref)","N_{Toys}")
        iHist += 1
        for nsigref in self.data[refPdfName][hmass]['zTrue']:
            hist.Fill(nsigref)
        hist.Draw()
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"m_{H} = "+str(hmass)+" GeV/c^{2}")
        tlatex.SetTextAlign(32)
        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
        tmpDat = self.data[refPdfName][hmass]['zTrue']
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.1f}".format(median(tmpDat)))
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.1f}".format(mean(tmpDat)))
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.1f}".format(stddev(tmpDat)))
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"N_{{out of hist}}: {0:.0f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
        self.setYMaxAndDrawVertLines(hist,None)
        canvas.RedrawAxis()
        saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Z_Ref"+refPdfName)
        canvas.Clear()

        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          hist = root.TH1F("hist"+str(iHist),"",50,-3,3)
          setHistTitles(hist,"N_{sig}(Alt)/#DeltaN_{sig}(Alt)","N_{Toys}")
          iHist += 1
          for nsigalt in self.data[refPdfName][hmass][pdfAltName]['z']:
            hist.Fill(nsigalt)
          hist.Draw()
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          tmpDat = self.data[refPdfName][hmass][pdfAltName]['z']
          tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.1f}".format(median(tmpDat)))
          tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.1f}".format(mean(tmpDat)))
          tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.1f}".format(stddev(tmpDat)))
          tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"N_{{out of hist}}: {0:.1f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
          self.setYMaxAndDrawVertLines(hist,None)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Z_Ref"+refPdfName+"_Alt"+pdfAltName)
          canvas.Clear()

    #  ##### NSigRef/Alt Plots
    #  for hmass in self.sigMasses:
    #    hist = root.TH1F("hist"+str(iHist),"",50,-200,200)
    #    setHistTitles(hist,"N_{sig}(Ref)","N_{Toys}")
    #    iHist += 1
    #    for nsigref in self.data[refPdfName][hmass]['nTrue']:
    #        hist.Fill(nsigref)
    #    hist.Draw()
    #    tlatex.SetTextAlign(12)
    #    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #    tlatex.SetTextAlign(12)
    #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #    tlatex.SetTextAlign(32)
    #    tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #    tmpDat = self.data[refPdfName][hmass]['nTrue']
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.1f}".format(median(tmpDat)))
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.1f}".format(mean(tmpDat)))
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.1f}".format(stddev(tmpDat)))
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"N_{{out of hist}}: {0:.0f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
    #    self.setYMaxAndDrawVertLines(hist,None)
    #    canvas.RedrawAxis()
    #    saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_NSigRef_Ref"+refPdfName)
    #    canvas.Clear()

    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      hist = root.TH1F("hist"+str(iHist),"",50,-200,200)
    #      setHistTitles(hist,"N_{sig}(Alt)","N_{Toys}")
    #      iHist += 1
    #      for nsigalt in self.data[refPdfName][hmass][pdfAltName]['n']:
    #        hist.Fill(nsigalt)
    #      hist.Draw()
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #      tlatex.SetTextAlign(32)
    #      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #      tmpDat = self.data[refPdfName][hmass][pdfAltName]['n']
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.1f}".format(median(tmpDat)))
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.1f}".format(mean(tmpDat)))
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.1f}".format(stddev(tmpDat)))
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"N_{{out of hist}}: {0:.1f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
    #      self.setYMaxAndDrawVertLines(hist,None)
    #      canvas.RedrawAxis()
    #      saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_NSigAlt_Ref"+refPdfName+"_Alt"+pdfAltName)
    #      canvas.Clear()

    #  ##### deltaNSigRef/Alt Plots
    #  for hmass in self.sigMasses:
    #    hist = root.TH1F("hist"+str(iHist),"",50,0,100)
    #    setHistTitles(hist,"#DeltaN_{sig}(Ref)","N_{Toys}")
    #    iHist += 1
    #    for errref in self.data[refPdfName][hmass]['errTrue']:
    #        hist.Fill(errref)
    #    hist.Draw()
    #    tlatex.SetTextAlign(12)
    #    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #    tlatex.SetTextAlign(12)
    #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #    tlatex.SetTextAlign(32)
    #    tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #    self.setYMaxAndDrawVertLines(hist,None)
    #    canvas.RedrawAxis()
    #    saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_deltaNSigRef_Ref"+refPdfName)
    #    canvas.Clear()

    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      hist = root.TH1F("hist"+str(iHist),"",50,0,100)
    #      setHistTitles(hist,"#DeltaN_{sig}(Alt)","N_{Toys}")
    #      iHist += 1
    #      for erralt in self.data[refPdfName][hmass][pdfAltName]['err']:
    #        hist.Fill(erralt)
    #      hist.Draw()
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #      tlatex.SetTextAlign(32)
    #      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #      self.setYMaxAndDrawVertLines(hist,None)
    #      canvas.RedrawAxis()
    #      saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_deltaNSigAlt_Ref"+refPdfName+"_Alt"+pdfAltName)
    #      canvas.Clear()

    #  ##### Chi2/NDF Plots
    #  for hmass in self.sigMasses:
    #    hist = root.TH1F("hist"+str(iHist),"",20,0,2)
    #    setHistTitles(hist,"#chi^{2}/NDF","N_{Toys}")
    #    iHist += 1
    #    for chi2,ndf in zip(self.data[refPdfName][hmass]['chi2True'],self.data[refPdfName][hmass]['ndfTrue']):
    #        hist.Fill(chi2/ndf)
    #    hist.Draw()
    #    tlatex.SetTextAlign(12)
    #    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #    tlatex.SetTextAlign(12)
    #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #    tlatex.SetTextAlign(32)
    #    tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #    NDF = self.data[refPdfName][hmass]['ndfTrue'][0]
    #    tmpDat = [i/NDF for i in self.data[refPdfName][hmass]['chi2True']]
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.1f}".format(median(tmpDat)))
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.1f}".format(mean(tmpDat)))
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.1f}".format(stddev(tmpDat)))
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"N_{{out of hist}}: {0:.0f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
    #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.45,"NDF: {0:.0f}".format(NDF))
    #    self.setYMaxAndDrawVertLines(hist,None)
    #    canvas.RedrawAxis()
    #    saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2NDF_Ref"+refPdfName)
    #    canvas.Clear()

    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      hist = root.TH1F("hist"+str(iHist),"",20,0,2)
    #      setHistTitles(hist,"#chi^{2}/NDF","N_{Toys}")
    #      iHist += 1
    #      for chi2,ndf in zip(self.data[refPdfName][hmass][pdfAltName]['chi2'],self.data[refPdfName][hmass][pdfAltName]['ndf']):
    #        hist.Fill(chi2/ndf)
    #      hist.Draw()
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #      tlatex.SetTextAlign(32)
    #      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #      NDF = self.data[refPdfName][hmass][pdfAltName]['ndf'][0]
    #      tmpDat = [i/NDF for i in self.data[refPdfName][hmass][pdfAltName]['chi2']]
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Median: {0:.1f}".format(median(tmpDat)))
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Mean: {0:.1f}".format(mean(tmpDat)))
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"#sigma: {0:.1f}".format(stddev(tmpDat)))
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"N_{{out of hist}}: {0:.1f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.45,"NDF: {0:.0f}".format(NDF))
    #      self.setYMaxAndDrawVertLines(hist,None)
    #      canvas.RedrawAxis()
    #      saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2NDF_Ref"+refPdfName+"_Alt"+pdfAltName)
    #      canvas.Clear()

    #  ##### Chi2/NDFvPull
    #  for hmass in self.sigMasses:
    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      hist = root.TH2F("hist"+str(iHist),"",12,-2,2,10,0,2)
    #      setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","#chi^{2}/NDF")
    #      iHist += 1
    #      NDF = self.data[refPdfName][hmass][pdfAltName]['ndf'][0]
    #      chi2pVals = [i/NDF for i in self.data[refPdfName][hmass][pdfAltName]['chi2']]
    #      for pull,chi2pVal in zip(self.data[refPdfName][hmass][pdfAltName]['pull'],chi2pVals):
    #        hist.Fill(pull,chi2pVal)
    #      hist.Draw('col')
    #      xLine = median(self.data[refPdfName][hmass][pdfAltName]['pull'])
    #      line = root.TLine()
    #      line.SetLineColor(root.kBlue)
    #      line.SetLineWidth(2)
    #      line.SetLineStyle(2)
    #      line.DrawLine(xLine,0,xLine,1)
    #      line.SetLineStyle(1)
    #      for iY in range(1,hist.GetNbinsY()+1):
    #        binPullList = []
    #        binLowVal = hist.GetYaxis().GetBinLowEdge(iY)
    #        binHighVal = hist.GetYaxis().GetBinUpEdge(iY)
    #        for iEntry in range(len(self.data[refPdfName][hmass][pdfAltName]['pull'])):
    #          if chi2pVals[iEntry] >= binLowVal and chi2pVals[iEntry] < binHighVal:
    #            binPullList.append(self.data[refPdfName][hmass][pdfAltName]['pull'][iEntry])
    #        binPullMed = median(binPullList)
    #        line.DrawLine(binPullMed,binLowVal,binPullMed,binHighVal)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #      tlatex.SetTextAlign(32)
    #      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #      tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"NDF: {0:.0f}".format(NDF))
    #      canvas.RedrawAxis()
    #      saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2NDFvPull_Ref"+refPdfName+"_Alt"+pdfAltName)
    #      canvas.Clear()

    #  ##### Chi2/NDFvChi2/NDF
    #  for hmass in self.sigMasses:
    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      hist = root.TH2F("hist"+str(iHist),"",10,0,2,10,0,2)
    #      setHistTitles(hist,"Alternate #chi^2/NDF","Reference #chi^{2}/NDF")
    #      iHist += 1
    #      NDF = self.data[refPdfName][hmass][pdfAltName]['ndf'][0]
    #      NDFRef = self.data[refPdfName][hmass]['ndfTrue'][0]
    #      for chi2Alt,chi2True in zip(self.data[refPdfName][hmass][pdfAltName]['chi2'],self.data[refPdfName][hmass]['chi2True']):
    #        hist.Fill(chi2Alt/NDF,chi2True/NDFRef)
    #      hist.Draw('col')
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #      tlatex.SetTextAlign(32)
    #      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #      canvas.RedrawAxis()
    #      saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2NDFvChi2NDF_Ref"+refPdfName+"_Alt"+pdfAltName)
    #      canvas.Clear()

    #  ##### Chi2ProbvChi2Prob
    #  for hmass in self.sigMasses:
    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      hist = root.TH2F("hist"+str(iHist),"",10,0,1,10,0,1)
    #      setHistTitles(hist,"Alternate #chi^2 p-Value","Reference #chi^{2} p-Value")
    #      iHist += 1
    #      NDF = self.data[refPdfName][hmass][pdfAltName]['ndf'][0]
    #      NDFRef = self.data[refPdfName][hmass]['ndfTrue'][0]
    #      chi2pVals = [i/NDF for i in self.data[refPdfName][hmass][pdfAltName]['chi2']]
    #      for chi2Alt,chi2True in zip(self.data[refPdfName][hmass][pdfAltName]['chi2'],self.data[refPdfName][hmass]['chi2True']):
    #        hist.Fill(
    #                    scipy.stats.chi2.sf(chi2Alt,NDF),
    #                    scipy.stats.chi2.sf(chi2True,NDFRef),
    #                 )
    #      hist.Draw('col')
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #      tlatex.SetTextAlign(32)
    #      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #      canvas.RedrawAxis()
    #      saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2vChi2_Ref"+refPdfName+"_Alt"+pdfAltName)
    #      canvas.Clear()

    #  ##### NAltvNRef
    #  for hmass in self.sigMasses:
    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      hist = root.TH2F("hist"+str(iHist),"",20,-200,200,20,-200,200)
    #      setHistTitles(hist,"N_{sig}(Ref)","N_{sig}(Alt)")
    #      iHist += 1
    #      for nAlt,nRef in zip(self.data[refPdfName][hmass][pdfAltName]['n'],self.data[refPdfName][hmass]['nTrue']):
    #        hist.Fill(nRef,nAlt)
    #      hist.Draw('col')
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #      tlatex.SetTextAlign(12)
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #      tlatex.SetTextAlign(32)
    #      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #      canvas.RedrawAxis()
    #      saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_NAltvNRef_Ref"+refPdfName+"_Alt"+pdfAltName)
    #      canvas.Clear()

      ##### Chi2 B-Only Prob Plots
      for hmass in self.sigMasses:
        hist = root.TH1F("hist"+str(iHist),"",20,0,1)
        setHistTitles(hist,"#chi^{2} p-Value of Background Only Fit","N_{Toys}")
        iHist += 1
        for chi2,ndf in zip(self.data[refPdfName][hmass]['chi2BOnly'],self.data[refPdfName][hmass]['ndfBOnly']):
            hist.Fill(scipy.stats.chi2.sf(chi2,ndf))
        hist.Draw()
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"m_{H} = "+str(hmass)+" GeV/c^{2}")
        tlatex.SetTextAlign(32)
        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
        self.setYMaxAndDrawVertLines(hist,None)
        canvas.RedrawAxis()
        saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2BOnly_Ref"+refPdfName)
        canvas.Clear()

        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          hist = root.TH1F("hist"+str(iHist),"",20,0,1)
          setHistTitles(hist,"#chi^{2} p-Value of Fit","N_{Toys}")
          iHist += 1
          for chi2,ndf in zip(self.data[refPdfName][hmass][pdfAltName]['chi2'],self.data[refPdfName][hmass][pdfAltName]['ndf']):
            hist.Fill(scipy.stats.chi2.sf(chi2,ndf))
          hist.Draw()
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+self.pdfTitleMap[refPdfName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+self.pdfTitleMap[pdfAltName])
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          self.setYMaxAndDrawVertLines(hist,None)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Chi2_Ref"+refPdfName+"_Alt"+pdfAltName)
          canvas.Clear()

    ####### QQ plot
    #for hmass in self.sigMasses:
    #  allZTrues = []
    #  allZAlts = []
    #  for refPdfName in self.refPdfNameList:
    #    allZTrues.extend(self.data[refPdfName][hmass]['zTrue'])
    #    for pdfAltName in self.pdfAltNamesDict[refPdfName]:
    #      allZAlts.extend(self.data[refPdfName][hmass][pdfAltName]['z'])
    #  zRefQQ = scipy.stats.probplot(allZTrues)
    #  zAltQQ = scipy.stats.probplot(allZAlts)
    #  refGraph = root.TGraph()
    #  altGraph = root.TGraph()
    #  oneGraph = root.TGraph()
    #  oneGraph.SetPoint(0,-2,-2)
    #  oneGraph.SetPoint(1,2,2)
    #  refGraph.SetLineColor(root.kRed)
    #  altGraph.SetLineColor(root.kBlue)
    #  refGraph.SetMarkerColor(root.kRed)
    #  altGraph.SetMarkerColor(root.kBlue)
    #  iGraph = 0
    #  for x,y in zip(*zRefQQ[0]):
    #    refGraph.SetPoint(iGraph,x,y)
    #    iGraph += 1
    #  iGraph = 0
    #  for x,y in zip(*zAltQQ[0]):
    #    altGraph.SetPoint(iGraph,x,y)
    #    iGraph += 1
    #  hist = root.TH2F("hist"+str(iHist),"",1,-3.,3.,1,-3.,3.)
    #  iHist += 1
    #  setHistTitles(hist,"Gaussian Expected Quantiles","Gaussian Observed Quantiles")
    #  hist.Draw()
    #  oneGraph.Draw("L")
    #  refGraph.Draw("LP")
    #  altGraph.Draw("LP")
    #  tlatex.SetTextAlign(12)
    #  tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
    #  tlatex.SetTextAlign(12)
    #  tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"m_{H} = "+str(hmass)+" GeV/c^{2}")
    #  tlatex.SetTextAlign(32)
    #  tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
    #  canvas.RedrawAxis()
    #  saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_QQ")
    #  canvas.Clear()


  def setYMaxAndDrawVertLines(self,hist,x):
    ymax = 0.
    for i in range(1,hist.GetXaxis().GetNbins()+1):
      ymax = max(ymax,hist.GetBinContent(i))
    ymax * 1.5
    hist.GetYaxis().SetRangeUser(0.,ymax*1.5)
    if x == None:
        return
    line = root.TLine()
    line.SetLineColor(root.kRed)
    line.SetLineWidth(3)
    line.DrawLine(x,0,x,ymax*1.05)
    return line

if __name__ == "__main__":
  outDir = "output/"

  categories = []

  jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
  jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

  #categories += [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10BO",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #categories += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  #categories += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  categories += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  #categories += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

  dataDir = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
  dataFns8TeV = [
    "SingleMuRun2012Av1-22Jan2013",
    "SingleMuRun2012Bv1-22Jan2013",
    "SingleMuRun2012Cv1-22Jan2013",
    "SingleMuRun2012Dv1-22Jan2013",
    ]

  dataFns7TeV = [
    "SingleMuRun2011Av1",
    "SingleMuRun2011Bv1"
    ]
  dataFns7TeV = [dataDir+i+".root" for i in dataFns7TeV]
  dataFns8TeV = [dataDir+i+".root" for i in dataFns8TeV]

  logFile = open(outDir+"biasStudy.log",'w')
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
  logFile.write("# {0}\n\n".format(now))
  inputPklFiles = glob.glob(outDir+"*.pkl")
  if len(inputPklFiles)>0:
    for inputPkl in inputPklFiles:
      bs = BiasStudy(None,None,None,None,inputPkl=inputPkl)
      logFile.write(bs.outStr)
      bs.plot(outDir+"bias_")
  else:
    for category in categories:
      bs = BiasStudy(category,dataFns8TeV,"8TeV",100)
      logFile.write(bs.outStr)
      bs.plot(outDir+"bias_")
    
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
  logFile.write("\n\n# {0}\n".format(now))
  logFile.close()
  
