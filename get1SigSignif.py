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

from helpers import *
import makeCards
import fitOrderChooser
from singleUseScripts.biasPklToMu import getSMSigCounts

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT


def runStudy(catName,energyStr,truePdfName,dataFileNames,sigMasses):
      """
        Pure function so that we can do multiprocessing!!
      """

      dataTree = root.TChain()
      for i in dataFileNames:
        dataTree.Add(i+"/outtree"+catName)
      dataTree.SetCacheSize(10000000);
      dataTree.AddBranchToCache("*");

      truePdfFunc = None
      if truePdfName == "Bernstein" or truePdfName == "Chebychev" or truePdfName == "Polynomial" or truePdfName == "SumExp" or truePdfName == "SumPow" or truePdfName == "Laurent" or truePdfName == "ExpTimesBernstein" or truePdfName == "ExpTimesChebychev" or truePdfName == "ExpTimesPolynomial":
        truePdfFunc = getattr(fitOrderChooser,"makePDFBak"+truePdfName)
      else:
        truePdfFunc = getattr(makeCards,"makePDFBak"+truePdfName)

      dimuonMass = root.RooRealVar("dimuonMass","m [GeV/c^{2}]",110.,160.)
      dimuonMass.setBins(50)
      dimuonMass.setRange("exprange",120,160)
      dimuonMass.setRange("whole",110,160)
      dimuonMass.setRange("low",110,120) # Silly ranges for old fit functionality
      dimuonMass.setRange("high",130,160)
      dimuonMass.setRange("signal",120,130)
      dimuonMass.setRange("signalfit",110,140)
      dimuonMass.setRange("annaRegion",123.5,127.5)
      dimuonMassArgSet = root.RooArgSet(dimuonMass)
      wTrue = root.RooWorkspace("wTrue")
      wTrueImport = getattr(wTrue,"import")

      canvas = root.TCanvas("canvas"+catName+energyStr+truePdfName)
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
      realDataHist = realData.binnedClone("realDataHist"+catName+energyStr)
      nData = realData.sumEntries()
      realDataZ = root.RooDataSet("realDataZ"+catName+energyStr,
                                      "realDataZ"+catName+energyStr,
                                          dataTree,root.RooArgSet(dimuonMassZ)
                                        )

      ### Make Bak Pdfs

      trashParamList, trashBakNormTup, trashDebug, trueOrder = truePdfFunc(truePdfName+catName+energyStr,realData,dimuonMass,110,160,wTrueImport,dimuonMassZ,realDataZ)
      truePdf = wTrue.pdf("bak")
      truePdf.SetName(truePdfName)
      truePdf.SetTitle("True PDF ")

      nDataVar = root.RooFit.RooConst(nData)
      nBakVar = root.RooRealVar("nBak","N_{B}",nData/2.,nData*2)
      truePdfE = root.RooExtendPdf(truePdfName+"E","True PDF Extended",truePdf,nBakVar)

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

      ### Now load Signal PDFs
      nSigVarBounds = nData/2.
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
        nSigSMs.append(getSMSigCounts(catName,hmass,energy=energyStr))

      result = {}

      ### Do S+B Fits
      for hmass,sigPdf,sigPdfE,nSigSM in zip(sigMasses,sigPdfs,sigPdfEs,nSigSMs):
        truePdfPlusSigPdf = root.RooAddPdf("truePdfPlusSigPdf"+catName+energyStr,"",root.RooArgList(truePdfE,sigPdfE))
        fr = truePdfPlusSigPdf.fitTo(realData,
                             PRINTLEVEL,
                             root.RooFit.Save(1)
                           )
        #frPars = fr.floatParsFinal()
        #for i in range(frPars.getSize()):
        #  frPars[i].Print()
        #nSigVar.Print()

        result[hmass] = nSigVar.getError()
        #result[hmass] = nSigVar.getError()/nSigSM

        # Debug plot for fit to data
        frame = dimuonMass.frame()
        chi2RealDataVar = truePdfPlusSigPdf.createChi2(realDataHist)
        ndfRealData = dimuonMass.getBins() - 1  # b/c roofit normalizes
        ndfRealData -= rooPdfNFreeParams(truePdfPlusSigPdf,realDataHist)
        realData.plotOn(frame)
        errVisArg = root.RooFit.VisualizeError(fr,1,True)
        errFillArg = root.RooFit.FillStyle(3001)
        truePdfPlusSigPdf.plotOn(frame,root.RooFit.Range('low,signal,high'),root.RooFit.NormRange('low,signal,high'),errVisArg,errFillArg,root.RooFit.FillColor(root.kGreen-7))
        truePdfPlusSigPdf.plotOn(frame,root.RooFit.Range('low,signal,high'),root.RooFit.NormRange('low,signal,high'),root.RooFit.Components(truePdf.GetName()),root.RooFit.LineStyle(2),root.RooFit.LineColor(root.kRed+1))
        truePdfPlusSigPdf.plotOn(frame,root.RooFit.Range('low,signal,high'),root.RooFit.NormRange('low,signal,high'))
        #truePdfPlusSigPdf.plotOn(frame,root.RooFit.Range('low,signal,high'),root.RooFit.NormRange('low,signal,high'),root.RooFit.Components(sigPdf.GetName()),root.RooFit.LineColor(root.kRed+1))
        
        frame.Draw()
        frame.SetTitle("")
        frame.GetYaxis().SetTitle("Events / 1 GeV/c^{2}")
        tlatex.SetTextAlign(12)
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
        tlatex.DrawLatex(0.02+gStyle.GetPadLeftMargin(),0.85,"Ref PDF: "+truePdfName)
        tlatex.SetTextAlign(32)
        tlatex.DrawLatex(0.99-gStyle.GetPadRightMargin(),0.96,catName+" "+energyStr)
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.85,"Ref. S+B Fit to Real Data")
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.80,"Ref. GOF: {0:.2f}".format(scipy.stats.chi2.sf(chi2RealDataVar.getVal(),ndfRealData)))
        tlatex.DrawLatex(0.97-gStyle.GetPadRightMargin(),0.75,"Ref. #chi^{{2}}/NDF: {0:.2f}".format(chi2RealDataVar.getVal()/ndfRealData))
        canvas.SaveAs("output/debug_oneSig_RealData_"+truePdfName+"_"+catName+"_"+energyStr+"_"+str(hmass)+".png")

      return result

if __name__ == "__main__":

  outDir = "output/"

  pdfName = "MSSM"

  ########################################
  ### Define which masses to run over

  #sigMasses = range(115,156,5)
  sigMasses = [115,120,125,130,135,140,145,150,155]

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

  categories += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  categories += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  categories += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  categories += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  categories += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

  ########################################
  ### Directory and file names

  dataDir = getDataStage2Directory()

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

  data = {}
  for energy,fns in zip(['7TeV','8TeV'],[dataFns7TeV,dataFns8TeV]):
    data[energy] = {}
    for category in categories:
      cat = category[0]
      data[energy][cat] = runStudy(cat,energy,pdfName,fns,sigMasses)
  outPklFile = open(outDir+"oneSig.pkl",'w')
  outTxtFile = open(outDir+"oneSig.txt",'w')
  cPickle.dump(data,outPklFile)
  outPklFile.close()
  for energy in ['7TeV','8TeV']:
    print energy
    outTxtFile.write(energy+'\n')
    for category in categories:
      print "  "+category[0]
      outTxtFile.write("  "+category[0]+'\n')
      for hmass in sigMasses:
        line =  "    {0}: {1:0.1f}".format(hmass,data[energy][category[0]][hmass])
        print line
        outTxtFile.write(line+'\n')
  outTxtFile.close()
