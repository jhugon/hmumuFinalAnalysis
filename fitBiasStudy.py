#! /usr/bin/env python

from ROOT import gSystem

import datetime
import sys
import os
import os.path
import re
import math
import cPickle
import ROOT as root
root.gSystem.Load('libRooFit')
root.gROOT.SetBatch(True)
import scipy.stats

from multiprocessing import Pool
import itertools
from itertools import repeat as itrRepeat

from helpers import *
import makeCards
import fitOrderChooser
from singleUseScripts.biasPklToMu import getSMSigCounts

from numpy import mean, median, corrcoef, percentile
from numpy import std as stddev

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

NPROCS = 1

TITLEMAP = {
  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
  #"Jets01PassCatAll" : "0,1-Jet Tight Combination",
                        
  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
  "Jets01FailPtG10EE": "0,1-Jet Loose EE",
  #"Jets01FailCatAll" : "0,1-Jet Loose Combination",

  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",
}

def medianAbsoluteDeviation(data,axis=None):
  return numpy.median(numpy.absolute(data - numpy.median(data, axis=axis)), axis=axis)
 
def madStddev(a,rescaleFactor=1.0/scipy.stats.norm.isf(1/4.)):
  """
  the default rescaleFactor is appropriate for the normal distribution
  """
  return medianAbsoluteDeviation(a)*rescaleFactor

def getOrdinalStr(inInt):
  result = str(inInt)
  if result[-1] == "1":
    result += "st"
  elif result[-1] == "2":
    result += "nd"
  elif result[-1] == "3":
    result += "rd"
  else:
    result += "th"
  return result

class PdfTitleMap(object):
  def __init__(self,data):
    self.data = data
  def __getitem__(self,key):
    orderMatch = re.match(r"([\d]+)(.+)",key)
    if orderMatch:
      keyNoOrder = orderMatch.group(2)
      order = orderMatch.group(1)
      order = getOrdinalStr(order)
      if keyNoOrder in self.data:
        return order+"-Order "+self.data[keyNoOrder]
      else:
        raise KeyError(keyNoOrder)
    elif key in self.data:
      return self.data[key]
    else:
      raise KeyError(key)
  def __setitem__(self,key,value):
    self.data[key] = value

PDFTITLEMAP = PdfTitleMap({
    "ExpLog":"Exp(p_{1}m^{2}+p_{2}m+p_{3}ln(m))",
    "MOverSq":"#frac{m}{(m-p_{1})^{2}}",
    "Old":"Voigtian+Exp",
    "ExpMOverSq":"#frac{Exp(p_{1}m)}{(m-p_{2})^{2}}",
    "ExpMOverSqP0":"#frac{Exp(-p_{1}^{2}m)}{(m-p_{2})}*(#frac{1}{m-p_{2}}+p_{3}^{2}m)",
    "ExpMOverSqP0New":"e^{-p_{1}^{2}m}/(m-p_{2})^{2}+p_{3}^{2}e^{-p_{1}^{2}m}",
    "Bernstein":"Bernstein",
    "Chebychev":"Chebychev",
    "Polynomial":"Polynomial",
    "SumExp":"Sum of Exponentials",
    "SumPow":"Sum of Power Functions",
    "Laurent":"Laurent",
    "ExpTimesBernstein":"Exp*Bernstein",
    "ExpTimesChebychev":"Exp*Chebychev",
    "ExpTimesPolynomial":"Exp*Polynomial",
    "MSSM":"Exp#times(Breit-Wigner+#frac{1}{m^{2}})",
    "VoigtPMm2":"Voigtian+#frac{1}{m^{2}}",
    "VoigtPExpMm2":"Voigtian+#frac{Exp}{m^{2}}",
})

def runStudy(iJob,iJobGroup,catName,energyStr,truePdfName,pdfAltNameList,dataFileNames,sigMasses,sigInject,toysPerJob):
      """
        Pure function so that we can do multiprocessing!!
      """
      plotEveryNToys = 10

      tmpJobStr = "_job"+str(iJob)
      if iJobGroup != None:
        tmpJobStr = "_jobGrp"+str(iJobGroup)+tmpJobStr

      randomGenerator = root.RooRandom.randomGenerator()
      iSeed = 10001+iJob
      if iJobGroup != None:
        iSeed += 1000*iJobGroup
      print "iSeed: {0}".format(iSeed)
      randomGenerator.SetSeed(iSeed)

      dataTree = root.TChain()
      for i in dataFileNames:
        dataTree.Add(i+"/outtree"+catName)
      dataTree.SetCacheSize(10000000);
      dataTree.AddBranchToCache("*");

      truePdfOrder = None
      truePdfNameNoOrder = truePdfName
      trueOrderMatch = re.match(r"([\d]+)(.+)",truePdfNameNoOrder)
      if trueOrderMatch:
        truePdfNameNoOrder = trueOrderMatch.group(2)
        truePdfOrder = int(trueOrderMatch.group(1))
      altPdfOrderList = []
      altPdfNameNoOrderList = []
      for i in pdfAltNameList:
        altOrderMatch = re.match(r"([\d]+)(.+)",i)
        if altOrderMatch:
          altPdfNameNoOrderList.append(altOrderMatch.group(2))
          altPdfOrderList.append(int(altOrderMatch.group(1)))
        else:
          altPdfNameNoOrderList.append(i)
          altPdfOrderList.append(None)

      truePdfFunc = None
      if truePdfNameNoOrder == "Bernstein" or truePdfNameNoOrder == "Chebychev" or truePdfNameNoOrder == "Polynomial" or truePdfNameNoOrder == "SumExp" or truePdfNameNoOrder == "SumPow" or truePdfNameNoOrder == "Laurent" or truePdfNameNoOrder == "ExpTimesBernstein" or truePdfNameNoOrder == "ExpTimesChebychev" or truePdfNameNoOrder == "ExpTimesPolynomial":
        truePdfFunc = getattr(fitOrderChooser,"makePDFBak"+truePdfNameNoOrder)
      else:
        truePdfFunc = getattr(makeCards,"makePDFBak"+truePdfNameNoOrder)
      pdfAltFuncList = []
      for i in altPdfNameNoOrderList:
        if i == "Bernstein" or i == "Chebychev" or i == "Polynomial" or i == "SumExp" or i == "SumPow" or i == "Laurent" or i == "ExpTimesBernstein" or i == "ExpTimesChebychev" or i == "ExpTimesPolynomial":
          pdfAltFuncList.append(getattr(fitOrderChooser,"makePDFBak"+i))
        else:
          pdfAltFuncList.append(getattr(makeCards,"makePDFBak"+i))

      #dimuonMass = root.RooRealVar("dimuonMass","m [GeV/c^{2}]",110.,170.)
      dimuonMass = root.RooRealVar("dimuonMass","m [GeV/c^{2}]",110.,160.)
      #dimuonMass.setBins(60)
      dimuonMass.setBins(50)
      dimuonMass.setRange("exprange",120,160)
      dimuonMass.setRange("whole",110,160)
      dimuonMass.setRange("low",110,120) # Silly ranges for old fit functionality
      #dimuonMass.setRange("high",130,170)
      dimuonMass.setRange("high",130,160)
      dimuonMass.setRange("signal",120,130)
      dimuonMass.setRange("signalfit",110,140)
      dimuonMass.setRange("annaRegion",123.5,127.5)
      dimuonMassArgSet = root.RooArgSet(dimuonMass)
      wTrue = root.RooWorkspace("wTrue")
      wTrueToy = root.RooWorkspace("wTrueToy")
      wTrueImport = getattr(wTrue,"import")
      wTrueToyImport = getattr(wTrueToy,"import")

      canvas = root.TCanvas("canvas"+catName+energyStr+truePdfName+str(iJob))
      tlatex = root.TLatex()
      tlatex.SetNDC()
      tlatex.SetTextFont(root.gStyle.GetLabelFont())
      tlatex.SetTextSize(0.04)

      # Hack to Make makePDFBakOld work
      minMassZ = 88.
      maxMassZ = 94.
      dimuonMassZ = root.RooRealVar("dimuonMass","dimuonMass",minMassZ,maxMassZ)

      ### Load data
      
      realData = root.RooDataSet("realData"+catName+energyStr,
                                      "realData"+catName+energyStr,
                                          dataTree,root.RooArgSet(dimuonMass)
                                        )
      nData = realData.sumEntries()
      realDataZ = root.RooDataSet("realDataZ"+catName+energyStr,
                                      "realDataZ"+catName+energyStr,
                                          dataTree,root.RooArgSet(dimuonMassZ)
                                        )

      ### Make Bak Pdfs

      trashParamList, trashBakNormTup, trashDebug, trueOrder = truePdfFunc(truePdfName+catName+energyStr,realData,dimuonMass,110,160,wTrueImport,dimuonMassZ,realDataZ,order=truePdfOrder)
      truePdf = wTrue.pdf("bak")
      truePdf.SetName(truePdfName)
      truePdf.SetTitle("True PDF ")
      truePdf.fitTo(realData,
                             PRINTLEVEL
                           )

      trueToyPdfName = "trueToy"+catName+energyStr
      trashParamList, trashBakNormTup, trashDebug, trueToyOrder = truePdfFunc(trueToyPdfName,realData,dimuonMass,110,160,wTrueToyImport,dimuonMassZ,realDataZ,order=truePdfOrder)
      trueToyPdf = wTrueToy.pdf("bak")
      trueToyPdf.SetName(trueToyPdfName)
      assert(trueOrder == trueToyOrder)
      nDataVar = root.RooFit.RooConst(nData)
      truePdfE = root.RooExtendPdf(truePdfName+"E","True PDF Extended",truePdf,nDataVar)

      # Debug plot for fit to data
      if toysPerJob > 2:
        frame = dimuonMass.frame()
        realDataHist = realData.binnedClone("realDataHist"+catName+energyStr+str(iJob))
        chi2RealDataVar = truePdf.createChi2(realDataHist)
        ndfRealData = dimuonMass.getBins() - 1  # b/c roofit normalizes
        ndfRealData -= rooPdfNFreeParams(truePdf,realDataHist)
        realData.plotOn(frame)
        truePdf.plotOn(frame,root.RooFit.Range('low,signal,high'),root.RooFit.NormRange('low,signal,high'))
        frame.Draw()
        frame.SetTitle("")
        frame.GetYaxis().SetTitle("Events / 1 GeV/c^{2}")
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Ref PDF: "+truePdfName)
        tlatex.SetTextAlign(32)
        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,catName+" "+energyStr)
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Ref. Fit to Obs. Data")
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.80,"Ref. GOF: {0:.2f}".format(scipy.stats.chi2.sf(chi2RealDataVar.getVal(),ndfRealData)))
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Ref. #chi^{{2}}/NDF: {0:.2f}".format(chi2RealDataVar.getVal()/ndfRealData))
        canvas.SaveAs("output/debug_RealData_"+truePdfName+"_"+catName+"_"+energyStr+tmpJobStr+".png")

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
      for pdfAltName,pdfAltNameNoOrder,pdfAltOrder,pdfAltFunc in zip(pdfAltNameList,altPdfNameNoOrderList,altPdfOrderList,pdfAltFuncList):
        pdfName = "alt"+catName+energyStr+"_"+pdfAltName
        wAlt = root.RooWorkspace("wAlt"+catName+energyStr+"_"+pdfAltName)
        pdfAltFunc(pdfName,realData,dimuonMass,110,160,getattr(wAlt,"import"),dimuonMassZ,realDataZ,order=pdfAltOrder)
        altPdf = wAlt.pdf("bak")
        altPdf.SetName(pdfName)
        # Make sure Voigt params are constant
        if pdfAltNameNoOrder == "Old":
          for x in rooArgSet2List(altPdf.getParameters(realData)):
            if "voit" in x.GetName():
              x.setConstant(True)
        pdfAltList.append(altPdf)
        pdfAltwList.append(wAlt)

      nBakVar = root.RooRealVar("nBak","N_{B}",nData/2.,nData*2)

      ### Now load Signal PDFs
      nSigVarBounds = nData/4
      nSigVar = root.RooRealVar("nSig","N_{S}",-nSigVarBounds,nSigVarBounds)
      sigPdfs = []
      sigPdfEs = []
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
        sigPdfE = root.RooExtendPdf(sigPdf.GetName()+"E",sigPdf.GetTitle()+" Extended",sigPdf,nSigVar)
        sigPdfEs.append(sigPdfE)

      ## Load the 1*SM N signal events
      nSigSMs = []
      for hmass in sigMasses:
        nSigSMs.append(getSMSigCounts(catName,hmass))

      ### Make results data structure and begin log
      data = {}
      data['meta'] = {'nData':nData,'sigInjectMu':sigInject}
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
        data[truePdfName][hmass]['nBakTrue'] = []
        data[truePdfName][hmass]['orderTrue'] = trueOrder
        for pdfAltName in pdfAltNameList:
          data[truePdfName][hmass][pdfAltName] = {'z':[],'pull':[],'n':[],'err':[],'chi2':[],'ndf':[]}
          data[truePdfName][hmass][pdfAltName]['nBak'] = []
          data[truePdfName][hmass][pdfAltName]['chi2BOnly'] = []
          data[truePdfName][hmass][pdfAltName]['ndfBOnly'] = []

      ### Toy Loop

      for iToy in range(toysPerJob):
        toyData = None
        toyDataHist = None
        if sigInject == 0.:
          toyData = truePdf.generate(root.RooArgSet(dimuonMass),int(nData))
          toyData.SetName("toyData"+catName+energyStr+str(iToy))
          toyDataHist = toyData.binnedClone("toyDataHist"+catName+energyStr+str(iToy))
        plotThisToy = (iToy % plotEveryNToys == 5)
        #plotThisToy = True
        for hmass,sigPdf,sigPdfE,nSigSM in zip(sigMasses,sigPdfs,sigPdfEs,nSigSMs):
          if sigInject != 0.:
            nSigVar.setVal(nSigSM*sigInject)
            truePdfPlusSigPdf = root.RooAddPdf("truePdfPlusSigPdf"+catName+energyStr+str(iToy),"",root.RooArgList(truePdfE,sigPdfE))
            toyData = truePdfPlusSigPdf.generate(root.RooArgSet(dimuonMass),int(nData))
            toyData.SetName("toyData"+catName+energyStr+str(iToy))
            toyDataHist = toyData.binnedClone("toyDataHist"+catName+energyStr+str(iToy))
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
          ndfTrue -= rooPdfNFreeParams(trueToyPdf,toyDataHist)
          chi2BOnlyVal = chi2TrueToyVar.getVal()
          ndfBOnlyVal = ndfTrue
          fracBakTrue = trueToyPdf.createIntegral(dimuonMassArgSet,root.RooFit.NormSet(dimuonMassArgSet),root.RooFit.Range("annaRegion"))
          nBakTrue = fracBakTrue.getVal()*nData
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
          ndfTrue -= rooPdfNFreeParams(trueToySBPdf,toyDataHist)
          nTrueToy = nSigVar.getVal() - nSigSM*sigInject
          errTrueToy = nSigVar.getError()
          if errTrueToy == 0.:
            continue
          if chi2TrueToyVar.getVal()==0.0:
            continue
          data[truePdfName][hmass]['nTrue'].append(nTrueToy/nSigSM)
          data[truePdfName][hmass]['errTrue'].append(errTrueToy/nSigSM)
          data[truePdfName][hmass]['chi2True'].append(chi2TrueToyVar.getVal())
          data[truePdfName][hmass]['ndfTrue'].append(ndfTrue)
          data[truePdfName][hmass]['zTrue'].append(nTrueToy/errTrueToy)
          data[truePdfName][hmass]['chi2BOnly'].append(chi2BOnlyVal)
          data[truePdfName][hmass]['ndfBOnly'].append(ndfBOnlyVal)
          data[truePdfName][hmass]['nBakTrue'].append(nBakTrue)
          if plotThisToy:
            trueToySBPdf.plotOn(frame,root.RooFit.LineColor(6))
          for pdfAlt,pdfAltName,color in zip(pdfAltList,pdfAltNameList,range(2,len(pdfAltList)+2)):
              ##### Get Background Only chi^2 and nBak
              pdfAlt.fitTo(toyData,
                                 PRINTLEVEL
                               )
              altChi2BOnlyVar = pdfAlt.createChi2(toyDataHist)
              ndfBOnly = dimuonMass.getBins() - 1  # b/c roofit normalizes
              ndfBOnly -= rooPdfNFreeParams(pdfAlt,toyDataHist)
              chi2BOnlyVal = altChi2BOnlyVar.getVal()
              fracBakAlt = pdfAlt.createIntegral(dimuonMassArgSet,root.RooFit.NormSet(dimuonMassArgSet),root.RooFit.Range("annaRegion"))
              nBakAlt = fracBakAlt.getVal()*nData
              ##### Do S+B Stuff
              altSBPdf = root.RooAddPdf("SB"+pdfAltName+catName+energyStr+str(hmass)+"_"+str(iToy),"",
                              root.RooArgList(pdfAlt,sigPdf),
                              root.RooArgList(nBakVar,nSigVar)
                          )
              altSBPdf.fitTo(toyData,
                              PRINTLEVEL
                              )
              altChi2Var = altSBPdf.createChi2(toyDataHist)
              ndfAlt = dimuonMass.getBins() - 1  # b/c roofit normalizes
              ndfAlt -= rooPdfNFreeParams(altSBPdf,toyDataHist)
              nAlt = nSigVar.getVal() - nSigSM*sigInject
              errAlt = nSigVar.getError()
              if errAlt == 0.:
                continue
              pull = (nAlt-nTrueToy)/errAlt
              data[truePdfName][hmass][pdfAltName]['n'].append(nAlt/nSigSM)
              data[truePdfName][hmass][pdfAltName]['err'].append(errAlt/nSigSM)
              data[truePdfName][hmass][pdfAltName]['chi2'].append(altChi2Var.getVal())
              data[truePdfName][hmass][pdfAltName]['ndf'].append(ndfAlt)
              data[truePdfName][hmass][pdfAltName]['z'].append(nAlt/errAlt)
              data[truePdfName][hmass][pdfAltName]['pull'].append(pull)
              data[truePdfName][hmass]['pullAll'].append(pull)
              data[truePdfName][hmass][pdfAltName]['chi2BOnly'].append(chi2BOnlyVal)
              data[truePdfName][hmass][pdfAltName]['ndfBOnly'].append(ndfBOnly)
              data[truePdfName][hmass][pdfAltName]['nBak'].append(nBakAlt)
              if plotThisToy:
                altSBPdf.plotOn(frame,root.RooFit.LineColor(color))
          if plotThisToy:
            frame.Draw()
            frame.SetTitle("")
            frame.GetYaxis().SetTitle("Events / 1 GeV/c^{2}")
            tlatex.SetTextAlign(12)
            tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
            tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Ref PDF: "+truePdfName)
            tlatex.SetTextAlign(32)
            tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,catName+" "+energyStr)
            tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Ref B-Only GOF: {0:.2f}".format(scipy.stats.chi2.sf(chi2BOnlyVal,ndfBOnlyVal)))
            tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.80,"Ref B-Only #chi^{{2}}/NDF: {0:.2f}".format(chi2BOnlyVal/ndfBOnlyVal))
            tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Ref S+B GOF: {0:.2f}".format(scipy.stats.chi2.sf(chi2TrueToyVar.getVal(),ndfTrue)))
            tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.70,"Ref S+B #chi^{{2}}/NDF: {0:.2f}".format(chi2TrueToyVar.getVal()/ndfTrue))
            canvas.SaveAs("output/debug_"+truePdfName+"_"+catName+"_"+energyStr+"_"+str(hmass)+tmpJobStr+"_Toy"+str(iToy)+".png")
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
      return data

def runStudyStar(argList):
  return runStudy(*argList)

################################################################################################
################################################################################################
################################################################################################

class BiasStudy:
  def __init__(self,category,dataFileNames,energyStr,sigMasses,refPdfNameList,pdfAltNamesDict,nToys=10,pklOutFnBase="output/biasData",inputPkl=None,processPool=None,iJobGroup=None,sigInject=0.):
    self.dataFileNames = dataFileNames
    self.sigMasses = sigMasses
    self.sigInject = sigInject
    self.iJobGroup = iJobGroup
    ## Try to load data from pkl file
    if inputPkl != None:
      if type(inputPkl) == str:
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
          self.sigInject = self.data['meta']['sigInjectMu']
        except Exception, err:
          print("Error loading data from pkl file: "+str(inputPkl))
          print(err)
          self.data = None
      elif type(inputPkl)==dict:
          self.data = inputPkl
          self.refPdfNameList = self.data['meta']['refPdfNameList']
          self.pdfAltNamesDict = self.data['meta']['pdfAltNamesDict']
          self.nToys = self.data['meta']['nToys']
          self.nData = self.data['meta']['nData']
          self.catName = self.data['meta']['catName']
          self.energyStr = self.data['meta']['energyStr']
          energyStr = self.energyStr
          self.sigMasses = self.data['meta']['sigMasses']
          self.sigInject = self.data['meta']['sigInjectMu']
      else:
          print("Error: unexpected type for input pickle filename or dict: "+type(inputPkl))
          print("Exiting.")
          sys.exit(1)
    else:
      self.catName = category[0]
      self.catCuts = category[1]
      self.energyStr = energyStr
      self.nToys = nToys
      tmpJobStr = ""
      if self.iJobGroup != None:
        tmpJobStr = "_jobGrp"+str(self.iJobGroup)
      self.pklOutFn = pklOutFnBase+"_"+self.catName+"_"+energyStr+tmpJobStr+".pkl"
      self.data = None
    catName = self.catName
    canvas = root.TCanvas()
    self.canvas = canvas
    ## Run
    if self.data==None:
      self.refPdfNameList = refPdfNameList
      self.pdfAltNamesDict = pdfAltNamesDict

      self.data = {'meta':{}}
      data = self.data
      data['meta']['sigMasses'] = self.sigMasses
      data['meta']['refPdfNameList'] = self.refPdfNameList
      data['meta']['pdfAltNamesDict'] = self.pdfAltNamesDict
      data['meta']['nToys'] = self.nToys
      data['meta']['catName'] = self.catName
      data['meta']['energyStr'] = self.energyStr
      data['meta']['sigInjectMu'] = self.sigInject
      self.iPklAutoSave = 1
      nProcesses = NPROCS
      nJobs = NPROCS
      for refPdfName,iRefPdfName in zip(self.refPdfNameList,range(len(self.refPdfNameList))):
        pdfAltNameList = self.pdfAltNamesDict[refPdfName]
        if processPool == None:
          mapResults = map(runStudyStar, itertools.izip(range(nJobs),itrRepeat(self.iJobGroup),itrRepeat(self.catName),itrRepeat(self.energyStr),itrRepeat(refPdfName),itrRepeat(pdfAltNameList),itrRepeat(self.dataFileNames),itrRepeat(self.sigMasses),itrRepeat(self.sigInject),itrRepeat(int(nToys/nJobs))))
        else:
          mapResults = processPool.map(runStudyStar, itertools.izip(range(nJobs),itrRepeat(self.iJobGroup),itrRepeat(self.catName),itrRepeat(self.energyStr),itrRepeat(refPdfName),itrRepeat(pdfAltNameList),itrRepeat(self.dataFileNames),itrRepeat(self.sigMasses),itrRepeat(self.sigInject),itrRepeat(int(nToys/nJobs))))
        for jobResults in mapResults:
          mergeDicts(data,jobResults)
        #if iRefPdfName != len(self.refPdfNameList)-1:
        #  pklFile = open(self.pklOutFn+"."+str(iRefPdfName),'w')
        #  cPickle.dump(data,pklFile)
        #  pklFile.close()
      pklFile = open(self.pklOutFn,'w')
      cPickle.dump(data,pklFile)
      pklFile.close()

    self.nData = self.data['meta']['nData']
    outStr = "#"*80+"\n"
    outStr = "#"*80+"\n"
    outStr += "\n"+self.catName +"  "+self.energyStr + "\n\n"

    outStr += "nToys: {0}\n".format(self.nToys)
    outStr += "nData Events: {0}\n".format(self.nData)
    outStr += "\n"

    data = self.data

    print outStr
    self.outStr = outStr

    self.pullSummaryDict = {}
    self.zSigmaSummaryDict = {}
    for refPdfName in self.refPdfNameList:
      self.pullSummaryDict[refPdfName] = {}
      self.zSigmaSummaryDict[refPdfName] = {}
      self.zSigmaSummaryDict[refPdfName]["zTrue"] = {}
      for hmass in self.sigMasses:
        self.zSigmaSummaryDict[refPdfName]['zTrue'][hmass] = stddev(data[refPdfName][hmass]['zTrue'])
        for pdfAltName in self.pdfAltNamesDict[refPdfName]:
          if not self.zSigmaSummaryDict[refPdfName].has_key(pdfAltName):
            self.zSigmaSummaryDict[refPdfName][pdfAltName] = {}
          if not self.pullSummaryDict[refPdfName].has_key(pdfAltName):
            self.pullSummaryDict[refPdfName][pdfAltName] = {}
          self.pullSummaryDict[refPdfName][pdfAltName][hmass] = median(data[refPdfName][hmass][pdfAltName]['pull'])
          self.pullSummaryDict[refPdfName]['orderRef'] = data[refPdfName][hmass]['orderTrue']
          self.zSigmaSummaryDict[refPdfName][pdfAltName][hmass] = stddev(data[refPdfName][hmass][pdfAltName]['z'])

  def plot(self,outputPrefix):

    canvas = self.canvas
    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    tlatex.SetTextSize(0.04)
    caption = TITLEMAP[self.catName]
    caption2 = ""
    caption3 = ""
    caption4 = ""
    iHist = 0

    for refPdfName in self.refPdfNameList:
      refPdfNameOrder = refPdfName
      if self.data[refPdfName][self.sigMasses[0]]['orderTrue'] != None:
        refPdfNameOrder = str(self.data[refPdfName][self.sigMasses[0]]['orderTrue'])+refPdfNameOrder
      ###### Pull plots 1D
      #for hmass in self.sigMasses:
      #  if len(self.pdfAltNamesDict[refPdfName])>1:
      #    hist = root.TH1F("hist"+str(iHist),"",30,-3,3)
      #    setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","N_{Toys}")
      #    iHist += 1
      #    for pull in self.data[refPdfName][hmass]['pullAll']:
      #        hist.Fill(pull)
      #    medianPull = median(self.data[refPdfName][hmass]['pullAll'])
      #    hist.Draw()
      #    tlatex.SetTextAlign(12)
      #    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
      #    tlatex.SetTextAlign(12)
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"All Alternate PDFs")
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.60,"Median: {0:.2f}".format(medianPull))
      #    tlatex.SetTextAlign(32)
      #    tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
      #    tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"{0:.2f}".format(medianPull))
      #    line = self.setYMaxAndDrawVertLines(hist,medianPull)
      #    canvas.RedrawAxis()
      #    saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_AllPulls_Ref"+refPdfName)
      #    canvas.Clear()

      #  for pdfAltName in self.pdfAltNamesDict[refPdfName]:
      #    hist = root.TH1F("hist"+str(iHist),"",30,-3,3)
      #    setHistTitles(hist,"(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)","N_{Toys}")
      #    iHist += 1
      #    for pull in self.data[refPdfName][hmass][pdfAltName]['pull']:
      #        hist.Fill(pull)
      #    medianPull = median(self.data[refPdfName][hmass][pdfAltName]['pull'])
      #    hist.Draw()
      #    tlatex.SetTextAlign(12)
      #    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
      #    tlatex.SetTextAlign(12)
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
      #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.60,"Median: {0:.2f}".format(medianPull))
      #    tlatex.SetTextAlign(32)
      #    tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
      #    line = self.setYMaxAndDrawVertLines(hist,medianPull)
      #    canvas.RedrawAxis()
      #    saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Pulls_Ref"+refPdfName+"_Alt"+pdfAltName)
      #    canvas.Clear()

#      ##### Median pull plots v. mass
#      for pdfAltName in self.pdfAltNamesDict[refPdfName]:
#        minx = 110
#        maxx = 160
#        axisHist = root.TH2F("axishist"+str(iHist),"",1,minx,maxx,1,-1,1)
#        setHistTitles(axisHist,"M_{H} [GeV/c^{2}]","Median[(N_{sig}(Alt)-N_{sig}(Ref))/#DeltaN_{sig}(Alt)]")
#        iHist += 1
#        graph = root.TGraph()
#        graphBand = root.TGraphErrors()
#        graphBand.SetPoint(0,minx,0.)
#        graphBand.SetPointError(0,0.,0.14)
#        graphBand.SetPoint(1,maxx,0.)
#        graphBand.SetPointError(1,0.,0.14)
#        graphBand.SetFillStyle(1001)
#        graphBand.SetFillColor(root.kGreen-9)
#        iHist += 1
#        iGraph = 0
#        for hmass in self.sigMasses:
#          medianPull = median(self.data[refPdfName][hmass][pdfAltName]['pull'])
#          graph.SetPoint(iGraph,hmass,medianPull)
#          iGraph += 1
#        axisHist.Draw()
#        graphBand.Draw("3")
#        graph.Draw("LP")
#        tlatex.SetTextAlign(12)
#        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
#        tlatex.SetTextAlign(12)
#        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
#        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
#        tlatex.SetTextAlign(32)
#        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
#        canvas.RedrawAxis()
#        saveAs(canvas,outputPrefix+self.catName+"_Pulls_Ref"+refPdfName+"_Alt"+pdfAltName)
#        canvas.Clear()
#
#      ##### Ref Z plots v. mass
#      minx = 110
#      maxx = 160
#      axisHist = root.TH2F("axishist"+str(iHist),"",1,minx,maxx,1,-5,5)
#      setHistTitles(axisHist,"M_{H} [GeV/c^{2}]","N_{sig}(Ref))/#DeltaN_{sig}(Ref)]")
#      iHist += 1
#      graph = root.TGraphErrors()
#      graphBand = root.TGraphErrors()
#      graphBand.SetPoint(0,minx,0.)
#      graphBand.SetPointError(0,0.,1.)
#      graphBand.SetPoint(1,maxx,0.)
#      graphBand.SetPointError(1,0.,1.)
#      graphBand.SetFillStyle(1001)
#      graphBand.SetFillColor(root.kGreen-9)
#      for iPoint,hmass in zip(range(len(self.sigMasses)),self.sigMasses):
#          zMeanTmp = mean(self.data[refPdfName][hmass]['zTrue'])
#          zSigmaTmp = stddev(self.data[refPdfName][hmass]['zTrue'])
#          graph.SetPoint(iPoint,hmass,zMeanTmp)
#          graph.SetPointError(iPoint,0.,zSigmaTmp)
#      axisHist.Draw()
#      graphBand.Draw("3")
#      graph.Draw("LP")
#      tlatex.SetTextAlign(12)
#      tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
#      tlatex.SetTextAlign(12)
#      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
#      tlatex.SetTextAlign(32)
#      tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
#      canvas.RedrawAxis()
#      saveAs(canvas,outputPrefix+self.catName+"_Z_Ref"+refPdfName)
#      canvas.Clear()
#
#      ##### Alt Z plots v. mass
#      for pdfAltName in self.pdfAltNamesDict[refPdfName]:
#        minx = 110
#        maxx = 160
#        axisHist = root.TH2F("axishist"+str(iHist),"",1,minx,maxx,1,-5,5)
#        setHistTitles(axisHist,"M_{H} [GeV/c^{2}]","N_{sig}(Alt))/#DeltaN_{sig}(Alt)]")
#        iHist += 1
#        graph = root.TGraphErrors()
#        graphBand = root.TGraphErrors()
#        graphBand.SetPoint(0,minx,0.)
#        graphBand.SetPointError(0,0.,1.)
#        graphBand.SetPoint(1,maxx,0.)
#        graphBand.SetPointError(1,0.,1.)
#        graphBand.SetFillStyle(1001)
#        graphBand.SetFillColor(root.kGreen-9)
#        for iPoint,hmass in zip(range(len(self.sigMasses)),self.sigMasses):
#            zMeanTmp = mean(self.data[refPdfName][hmass][pdfAltName]['z'])
#            zSigmaTmp = stddev(self.data[refPdfName][hmass][pdfAltName]['z'])
#            graph.SetPoint(iPoint,hmass,zMeanTmp)
#            graph.SetPointError(iPoint,0.,zSigmaTmp)
#        axisHist.Draw()
#        graphBand.Draw("3")
#        graph.Draw("LP")
#        tlatex.SetTextAlign(12)
#        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
#        tlatex.SetTextAlign(12)
#        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
#        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
#        tlatex.SetTextAlign(32)
#        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
#        canvas.RedrawAxis()
#        saveAs(canvas,outputPrefix+self.catName+"_Z_Ref"+refPdfName+"_Alt"+pdfAltName)
#        canvas.Clear()
#

      ##### Z Plots
      for hmass in self.sigMasses:
        hist = root.TH1F("hist"+str(iHist),"",50,-3,3)
        setHistTitles(hist,"N_{sig}(Ref)/\DeltaN_{sig}(Ref)","N_{Toys}")
        iHist += 1
        for nsigref in self.data[refPdfName][hmass]['zTrue']:
            hist.Fill(nsigref)
        fitFn = root.TF1("gaus"+str(hmass),"gaus",-3,3)
        fitResult = hist.Fit(fitFn,"LSMEQ") 
        hist.Draw()
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
        tlatex.SetTextAlign(12)
        tlatex.SetTextColor(root.kRed+1)
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
        tlatex.SetTextColor(1)
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"m_{H} = "+str(hmass)+" GeV/c^{2}")
        tlatex.SetTextAlign(32)
        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
        tmpDat = self.data[refPdfName][hmass]['zTrue']
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.68,"Mean: {0:.2f} #pm {1:.2f}".format(fitFn.GetParameter(1),fitFn.GetParError(1)))
        tlatex.SetTextColor(root.kRed+1)
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.58,"#sigma: {0:.2f} #pm {1:.2f}".format(fitFn.GetParameter(2),fitFn.GetParError(2)))
        tlatex.SetTextColor(1)
        #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"Median: {0:.2f}".format(median(tmpDat)))
        #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"MAD #sigma: {0:.2f}".format(madStddev(tmpDat)))
        #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.45,"Mean: {0:.2f}".format(mean(tmpDat)))
        #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.35,"#sigma: {0:.2f}".format(stddev(tmpDat)))
        #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.25,"N_{{out of hist}}: {0:.1f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
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
          fitResult = hist.Fit(fitFn,"LSMEQ") 
          hist.Draw()
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
          tlatex.SetTextAlign(12)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
          tlatex.SetTextColor(root.kRed+1)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
          tlatex.SetTextColor(1)
          tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.68,"m_{H} = "+str(hmass)+" GeV/c^{2}")
          tlatex.SetTextAlign(32)
          tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,caption)
          tmpDat = self.data[refPdfName][hmass][pdfAltName]['z']
          tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.68,"Mean: {0:.2f} #pm {1:.2f}".format(fitFn.GetParameter(1),fitFn.GetParError(1)))
          tlatex.SetTextColor(root.kRed+1)
          tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.58,"#sigma: {0:.2f} #pm {1:.2f}".format(fitFn.GetParameter(2),fitFn.GetParError(2)))
          tlatex.SetTextColor(1)
          #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.65,"Median: {0:.2f}".format(median(tmpDat)))
          #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.55,"MAD #sigma: {0:.2f}".format(madStddev(tmpDat)))
          #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.45,"Mean: {0:.2f}".format(mean(tmpDat)))
          #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.35,"#sigma: {0:.2f}".format(stddev(tmpDat)))
          #tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.25,"N_{{out of hist}}: {0:.1f}".format(hist.GetBinContent(0)+hist.GetBinContent(hist.GetNbinsX()+1)))
          self.setYMaxAndDrawVertLines(hist,None)
          canvas.RedrawAxis()
          saveAs(canvas,outputPrefix+self.catName+"_"+str(hmass)+"_Z_Ref"+refPdfName+"_Alt"+pdfAltName)
          canvas.Clear()

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
    #    tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
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
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Reference PDF: "+PDFTITLEMAP[refPdfNameOrder])
    #      tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.75,"Alternate PDF: "+PDFTITLEMAP[pdfAltName])
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

def mergeDicts(data,newData,multiJob=False):
  newDataKeys = newData.keys()
  dataKeys = data.keys()
  for key in newDataKeys:
    if key == 'meta':
      data[key]['nData'] = newData[key]['nData']
      if multiJob:
        data[key]['nToys'] += newData[key]['nToys']
    elif not (key in dataKeys):
      data[key] = newData[key]
    else:
      masses = sorted(newData[key].keys())
      massesOld = sorted(data[key].keys())
      assert(len(masses)==len(massesOld))
      for old,new in zip(masses,massesOld):
        if old != new:
          print "masses not equal:",old,new
          assert(False)
      for mass in masses:
        subkeys = sorted(newData[key][mass])
        subkeysOld = sorted(data[key][mass])
        assert(len(subkeys) == len(subkeysOld))
        for old,new in zip(subkeys,subkeysOld):
          if old != new:
            print "subkeys not equal:",old,new
            assert(False)
        for subkey in subkeys:
          if subkey == "orderTrue":
            assert(newData[key][mass][subkey] == data[key][mass][subkey])
          elif type(newData[key][mass][subkey])==list:
            assert(type(data[key][mass][subkey])==list)
            data[key][mass][subkey].extend(newData[key][mass][subkey])
          else:
            subsubkeys = sorted(newData[key][mass][subkey])
            subsubkeysOld = sorted(data[key][mass][subkey])
            assert(len(subsubkeys) == len(subsubkeysOld))
            for old,new in zip(subsubkeys,subsubkeysOld):
              if old != new:
                print "subsubkeys not equal:",old,new
                assert(False)
            for subsubkey in subsubkeys:
              assert(type(data[key][mass][subkey][subsubkey])==list)
              assert(type(newData[key][mass][subkey][subsubkey])==list)
              data[key][mass][subkey][subsubkey].extend(newData[key][mass][subkey][subsubkey])

def printBiasTable(dataCats,hmasses):
  catNames = sortCatNames(dataCats.keys())
  if len(catNames) == 0:
    return
  plainResult = ""
  latexResult = ""
  for catName in catNames:
    catTitle = TITLEMAP[catName]
    data = dataCats[catName]
    refNames = sorted(data.keys())
    allAltNames = set()
    for refName in refNames:
      allAltNames = allAltNames.union(data[refName].keys())
    allAltNames = sorted(list(allAltNames))
    allAltNames.pop(allAltNames.index("orderRef"))
    for altName in allAltNames:
      altTitle = PDFTITLEMAP[altName]
      if "#" in altTitle or "^" in altTitle:
        altTitle = '$'+altTitle+'$'
      altTitle = altTitle.replace("#","\\")
      plainResult += "############################################################\n"
      plainResult += catName+" Bias v. Mass for Alt: "+altName+"\n\n"
      plainResult += "\n{0:<10}".format("mH")+"{0:>15}\n".format("Reference")
      plainResult += "{0:<10}".format("")
      latexResult += "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
      latexResult += "%% "+catName+" Bias v. Mass for Alt: "+altName+"\n\n"
      latexResult += r"\begin{tabular}{|l|"+"r|"*len(refNames)+"} \\hline \n"
      latexResult += r"\multicolumn{"+str(len(refNames)+1)+r"}{|c|}{ \bf "+catTitle+' Bias for '+altTitle+r"} \\ \hline"+"\n"
      latexResult += r"\multicolumn{1}{|c|}{\multirow{2}{*}{$m_{H}$ [GeV/c$^{2}$]}} & \multicolumn{"+str(len(refNames))+r"}{c|}{Reference PDFs} \\ \cline{"+"2-{0}".format(len(refNames)+1)+"} "+"\n"
      for refName in refNames:
        plainResult += "{0:>15}".format(refName)
        nicePdfName = PDFTITLEMAP[refName]
        if "#" in nicePdfName:
          nicePdfName = '$'+nicePdfName+'$'
        latexResult += "& \multicolumn{1}{c|}{" +"{0:>15}".format(nicePdfName.replace("#","\\"))+"} "
      plainResult += "\n"
      latexResult += r"\\ \hline"+"\n"
      for hmass in hmasses:
        plainResult += "{0:<10.0f}".format(hmass)
        latexResult += "{0:10.0f} ".format(hmass)
        for refName in refNames:
           tmpBias = data[refName][altName][hmass]
           plainResult += "{0:>15.1%}".format(tmpBias)
           latexResult += ("& {0:15.1%} ".format(tmpBias)).replace("%","\%")
        plainResult += "\n"
        latexResult += r"\\ \hline"+"\n"
      plainResult += "\n\n"
      latexResult += "\\end{tabular}\n\n"
  print plainResult
  print latexResult

def printBiasSummary(dataCats):
  catNames = sorted(dataCats.keys())
  if len(catNames) == 0:
    return
  plainResult = ""
  latexResult = ""
  for catName in catNames:
    catTitle = TITLEMAP[catName]
    data = dataCats[catName]
    refNames = sorted(data.keys())
    allAltNames = set()
    for refName in refNames:
      allAltNames = allAltNames.union(data[refName].keys())
    allAltNames = sorted(list(allAltNames))
    allAltNames.pop(allAltNames.index("orderRef"))
    plainResult += "############################################################\n"
    plainResult += catName+" Maximum Bias\n\n"
    plainResult += "\n{0:<15}".format("Alternate")+"{0:>15}\n".format("Reference")
    plainResult += "{0:<15}".format("")
    latexResult += "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
    latexResult += "%% "+catName+" Maximum Bias\n\n"
    latexResult += r"\begin{tabular}{|l|"+"r|"*len(refNames)+"} \\hline \n"
    latexResult += r"\multicolumn{"+str(len(refNames)+1)+r"}{|c|}{ \bf "+catTitle+r" Maximum Bias} \\ \hline"+"\n"
    latexResult += r"\multicolumn{1}{|c|}{\multirow{3}{*}{Alternate PDFs}} & \multicolumn{"+str(len(refNames))+r"}{c|}{Reference PDFs} \\ \cline{"+"2-{0}".format(len(refNames)+1)+"} "+"\n"
    for refName in refNames:
      if data[refName]['orderRef'] != None:
        latexResult += "& \multicolumn{1}{c|}{"+getOrdinalStr(data[refName]['orderRef'])+"-Order} "
      else:
        latexResult += "&              ".format(getOrdinalStr(data[refName]['orderRef']))
    latexResult += r"\\ "+"\n"
    for refName in refNames:
      plainResult += "{0:>15}".format(refName)
      nicePdfName = PDFTITLEMAP[refName]
      if "#" in nicePdfName:
        nicePdfName = '$'+nicePdfName+'$'
      latexResult += "& \multicolumn{1}{c|}{" +"{0:>15}".format(nicePdfName.replace("#","\\"))+"} "
    plainResult += "\n"
    latexResult += r"\\ \hline"+"\n"
    for altName in allAltNames:
      plainResult += "{0:<15}".format(altName)
      nicePdfName = PDFTITLEMAP[altName]
      if "#" in nicePdfName:
        nicePdfName = '$'+nicePdfName+'$'
      latexResult += "{0:15} ".format(nicePdfName.replace("#","\\"))
      for refName in refNames:
          if not data[refName].has_key(altName):
            plainResult += "{0:>15}".format("-")
            latexResult += "& {0:15} ".format('-')
            continue
          maxBias = 0.
          absMaxBias = 0.
          for hmass in data[refName][altName]:
            tmpBias = data[refName][altName][hmass]
            tmpAbsBias = abs(tmpBias)
            if tmpAbsBias > absMaxBias:
              maxBias = tmpBias
              absMaxBias = tmpAbsBias
          plainResult += "{0:>15.1%}".format(maxBias)
          latexResult += ("& {0:15.1%} ".format(maxBias)).replace("%","\%")
      plainResult += "\n"
      latexResult += r"\\ \hline"+"\n"
    plainResult += "\n\n"
    latexResult += "\\end{tabular}\n\n"
  print plainResult
  print latexResult
          
      
def printDiagnosticSummary(dataCats,dataCatsZSig):
  catNames = sorted(dataCats.keys())
  if len(catNames) == 0:
    return
  assert(len(dataCats)==len(dataCatsZSig))
  plainResult = ""
  for catName in catNames:
    catTitle = TITLEMAP[catName]
    data = dataCats[catName]
    dataZSig = dataCatsZSig[catName]
    refNames = sorted(data.keys())
    allAltNames = set()
    for refName in refNames:
      allAltNames = allAltNames.union(data[refName].keys())
    allAltNames = sorted(list(allAltNames))
    allAltNames.pop(allAltNames.index("orderRef"))
    plainResult += "############################################################\n"
    plainResult += catName+" Maximum Bias\n\n"
    plainResult += "\n{0:<15}".format("Alternate")+"{0:>15}".format("Reference")
    for i in range(len(refNames)-1):
      plainResult += "{0:>15}".format("")
    plainResult += "{0:>15}".format("Z-Sigmas")
    plainResult += "\n"
    plainResult += "{0:>15}".format("")
    for refName in refNames:
      plainResult += "{0:>15}".format(refName,)
      nicePdfName = PDFTITLEMAP[refName]
    for refName in refNames:
       maxZSig = 0.
       for hmass in dataZSig[refName]['zTrue']:
         maxZSig = max(dataZSig[refName]['zTrue'][hmass],maxZSig)
       plainResult += "{0:>15.2f}".format(maxZSig)
    plainResult += "\n"
    for altName in allAltNames:
      plainResult += "{0:<15}".format(altName)
      nicePdfName = PDFTITLEMAP[altName]
      for refName in refNames:
          if not data[refName].has_key(altName):
            plainResult += "{0:>15}".format("-")
            continue
          maxBias = 0.
          absMaxBias = 0.
          for hmass in data[refName][altName]:
            tmpBias = data[refName][altName][hmass]
            tmpAbsBias = abs(tmpBias)
            if tmpAbsBias > absMaxBias:
              maxBias = tmpBias
              absMaxBias = tmpAbsBias
          plainResult += "{0:>15.1%}".format(maxBias)
      for refName in refNames:
          maxZSig = 0.
          for hmass in dataZSig[refName][altName]:
            maxZSig = max(dataZSig[refName][altName][hmass],maxZSig)
          plainResult += "{0:>15.2f}".format(maxZSig)
      plainResult += "\n"
    plainResult += "\n\n"
  print plainResult

def createSummaryMuDict(data):
  muDict = {}
  muDict = {'meta':data['meta']}
  pdfRefNames = data['meta']['refPdfNameList']
  altPdfName = 'MSSM'
  muDict['meta']['altPdf'] = altPdfName
  muDict['meta'].pop('pdfAltNamesDict')
  for refPdfName in pdfRefNames:
    muDict[refPdfName] = {}
    for hmass in data['meta']['sigMasses']:
      # muR is (N(alt)-N(ref))/N(SM) similar to what we used before
      # muAlt is N(alt)/N(SM), just in case you want it
      # muAltUnc is the uncertainty on N(alt) divided by N(SM)
      # muRef is N(ref)/N(SM), just in case you want it
      # muRefUnc is the uncertainty on N(ref) divided by N(SM)
      # You shouldn't need it
      muAltList = data[refPdfName][hmass][altPdfName]["n"]
      muAltUncList = data[refPdfName][hmass][altPdfName]["err"]
      muRefList = data[refPdfName][hmass]["nTrue"]
      muRefUncList = data[refPdfName][hmass]["errTrue"]
      muRList = [i-j for i,j in zip(muAltList,muRefList)]
      muR = median(muRList)
      muAlt = median(muAltList)
      muAltUnc = median(muAltUncList)
      muRef = median(muRefList)
      muRefUnc = median(muRefUncList)
      muDict[refPdfName][hmass] = {}
      muDict[refPdfName][hmass]["muR"]      = muR
      muDict[refPdfName][hmass]["muAlt"]    = muAlt
      muDict[refPdfName][hmass]["muAltUnc"] = muAltUnc
      muDict[refPdfName][hmass]["muRef"]    = muRef
      muDict[refPdfName][hmass]["muRefUnc"] = muRefUnc
  return muDict

def printBiasCombination(data):
  hmasses = None
  for catName in data:
    hmasses = data[catName]['meta']['sigMasses']
    continue
  print "Weighted Averages of Categories"
  combPullDict = {}
  refNames = sorted(data[catName].keys())
  refNames.remove('meta')
  for refName in refNames:
    print refName
    print "{0:<7} {1:>8} {2:>8} {3:>8}".format('hmass', 'Mu','deltaMu','Pull')
    for hmass in hmasses:
      normFactor = 0.
      sumOfWeightedMeas = 0.
      for catName in data:
        var = data[catName][refName][hmass]['muAltUnc']**2
        sumOfWeightedMeas += data[catName][refName][hmass]['muR'] / var
        normFactor += 1./var
      weightedMean = sumOfWeightedMeas/normFactor
      weightedStdDev = 1./sqrt(normFactor)
      print "{0:<7} {1:>8.2f} {2:>8.2f} {3:>8.2f}".format(hmass, weightedMean,weightedStdDev,weightedMean/weightedStdDev)
      if not combPullDict.has_key(hmass):
        combPullDict[hmass] = {}
      combPullDict[hmass][refName] = weightedMean/weightedStdDev

  print
  print

  rowStr = "{0:<8}"
  rowArgs = [r"m_H [\GeVcc{}]"]
  rowI = 1
  for refName in refNames:
    rowStr += " & {"+str(rowI)+":<14}"
    rowArgs.append(PDFTITLEMAP[refName])
    rowI += 1
  rowStr += r"\\ \hline \hline"
  print rowStr.format(*rowArgs)
  rowStr = "{0:<8}"
  rowArgs = ["hmass"]
  rowI = 1
  for refName in refNames:
    rowStr += " {"+str(rowI)+":<14}"
    rowArgs.append(refName)
    rowI += 1
  print rowStr.format(*rowArgs)
  for hmass in hmasses:
    rowStr = "{0:<8}"
    rowArgs = [hmass]
    rowI = 1
    for refName in refNames:
      rowStr += " & {"+str(rowI)+":<12.2f}"
      rowArgs.append(combPullDict[hmass][refName])
      rowI += 1
    print rowStr.format(*rowArgs)
    

if __name__ == "__main__":
  helpStr = "./fitBiasStudy.py [jobGroupNumber] [categoryName]\n  where jobGroupNumber is an int that will be added to the random number seed (*1000)\n    and the output pkl file name\n  if there is a jobGroupNumber, no plots or summary will be produced.\n  If categoryName is present, then only that category will be run,\n    otherwise a group of categories defined in the script will all be run."
  iJobGroup = None
  catToRun = None
  if len(sys.argv) > 3:
    print("Error: Too many arguments.  Run like:")
    print(helpStr)
    sys.exit(1)
  if len(sys.argv) >= 2:
    iJobGroup = int(sys.argv[1])
    print("iJobGroup: "+str(iJobGroup))
  if len(sys.argv) == 3:
    catToRun = sys.argv[2]
    print("Running only on category: "+catToRun)
    if not TITLEMAP.has_key(catToRun):
      print("Error: category name '"+catToRun+"' not recognized, exiting.")
      sys.exit(1)

  ############################################
  ### Define output directory

  outDir = "output/"

  ############################################
  ### Define number of toys to run over

  nToys = 1

  ############################################
  ### Define which reference functions to use

  refPdfNameList = [
          "Old",
          "ExpMOverSq",
          "Bernstein",
          "SumExp",
          "VoigtPExpMm2",
          "VoigtPMm2",
      #    "ExpMOverSqP0",
      #    "ExpLog",
      #    "MOverSq",
      #    "SumPow",
      #    "Laurent",
      #    "Chebychev",
      #    "Polynomial",
  ]
  ###############################################
  ### Define which alternate functions to test
  ### for each given reference function

  pdfAltNamesDict = {
      "ExpLog":["MSSM"],
      "MOverSq":["MSSM"],
      "Old":["MSSM"],
      "ExpMOverSq":["MSSM"],
      "ExpMOverSqP0":["MSSM"],
      "Bernstein":["MSSM"],          
      "SumExp":["MSSM"],          
      "SumPow":["MSSM"],          
      "Laurent":["MSSM"],          
      "Chebychev":["MSSM"],          
      "Polynomial":["MSSM"],          
      "VoigtPExpMm2":["MSSM"],
      "VoigtPMm2":["MSSM"],
  }

  ########################################
  ### Define which masses to run over

  #sigMasses = range(115,156,5)
  sigMasses = [115,120,125,130,135,140,145,150,155]

  sigInject = 0.

  ########################################

  jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
  jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

  categoriesAll = ["BB","BO","BE","OO","OE","EE"]
  categories = []

  ########################################
  ### Define which categories to run over

  #categories += [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10BO",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10BE",  "dimuonPt>10." +jet01PtCuts]]

  #categories += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #categories += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  categories += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  #categories += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  #categories += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

  ########################################
  ### Directory and file names

  dataDir = getDataStage2Directory()
  #dataDir = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
  #dataDir = "/afs/cern.ch/work/j/jhugon/public/hmumuNtuplesLevel2/unzipped/"
  #dataDir = "/cms/data/store/user/jhugon/hmumu/stage2/"
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

  #############################################

  if catToRun != None:
    categories = [[catToRun,""]]

  allSummaries = {}
  allZSigmaSummaries = {}
  muSummaries = {}
  tmpJobGroupStr = ""
  if iJobGroup != None:
    tmpJobGroupStr = "_jobGrp"+str(iJobGroup)
#  logFile = open(outDir+"biasStudy.log",'w')
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
#  logFile.write("# {0}\n\n".format(now))
  inputPklFiles = glob.glob(outDir+"*.pkl")
  tmpInputPklFiles = inputPklFiles
  if len(inputPklFiles)>0 and iJobGroup==None:
    foundJobGroupPkl = False
    foundNotJobGroupPkl = False
    for inputPkl in inputPklFiles:
      if "_jobGrp" in inputPkl:
        foundJobGroupPkl = True
      else:
        foundNotJobGroupPkl = True
    if foundNotJobGroupPkl and foundJobGroupPkl:
        print "Error: found .pkl files containing '_jobGrp' and not containing '_jobGrp'."
        print "  Can only process one or the other.  "
        print "  Please delete one or the other and try again. Exiting."
        sys.exit(1)
    if foundNotJobGroupPkl:
      for inputPkl in inputPklFiles:
        print "Running over input pkl file: "+inputPkl
        bs = BiasStudy(None,None,None,None,None,None,None,inputPkl=inputPkl)
#        logFile.write(bs.outStr)
        bs.plot(outDir+"bias_")
        allSummaries[bs.catName] = bs.pullSummaryDict
        allZSigmaSummaries[bs.catName] = bs.zSigmaSummaryDict
        muSummaries[bs.catName] = createSummaryMuDict(bs.data)
    else:
      # Identify basenames to combine job groups
      basenames = set()
      for inputPklFn in inputPklFiles:
        match = re.match(r"(.+jobGrp)[\d]+\.pkl",inputPklFn)
        assert(match)
        tmpBase = match.group(1)
        if not (tmpBase+"*.pkl") in basenames:
            basenames.add((tmpBase+"*.pkl"))
      for globStr in basenames:
        fns = glob.glob(globStr)
        resultData = None
        for tmpFn in fns:
          tmpF = open(tmpFn)
          tmpD = cPickle.load(tmpF)
          if resultData == None:
            resultData = tmpD 
          else:
            mergeDicts(resultData,tmpD,True)
          tmpF.close()
        bs = BiasStudy(None,None,None,None,None,None,None,inputPkl=resultData)
#        logFile.write(bs.outStr)
        bs.plot(outDir+"bias_")
        allSummaries[bs.catName] = bs.pullSummaryDict
        allZSigmaSummaries[bs.catName] = bs.zSigmaSummaryDict
        muSummaries[bs.catName] = createSummaryMuDict(bs.data)
  else:
    processPool = None
    if NPROCS > 1:
      processPool = Pool(processes=NPROCS)
    for category in categories:
      bs = BiasStudy(category,dataFns8TeV,"8TeV",sigMasses,refPdfNameList,pdfAltNamesDict,nToys,processPool=processPool,iJobGroup=iJobGroup,sigInject=sigInject)
#      logFile.write(bs.outStr)
      if iJobGroup == None:
        bs.plot(outDir+"bias_")
        allSummaries[bs.catName] = bs.pullSummaryDict
        allZSigmaSummaries[bs.catName] = bs.zSigmaSummaryDict
        muSummaries[bs.catName] = createSummaryMuDict(bs.data)
  printBiasTable(allSummaries,sigMasses)
  printBiasSummary(allSummaries)
  printDiagnosticSummary(allSummaries,allZSigmaSummaries)
    
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
#  logFile.write("\n\n# {0}\n".format(now))
#  logFile.close()
  printBiasCombination(muSummaries)
  
