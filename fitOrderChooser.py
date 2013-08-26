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

from multiprocessing import Pool
import itertools
from itertools import repeat as itrRepeat

from helpers import *
from makeCards import Param
from xsec import lumiDict

from numpy import mean, median, corrcoef, percentile
from numpy import std as stddev

#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

titleMap = {

  "Jets01PassPtG10BB": "0,1-Jet Tight BB",
  "Jets01PassPtG10BO": "0,1-Jet Tight BO",
  "Jets01PassPtG10BE": "0,1-Jet Tight BE",
  "Jets01PassPtG10OO": "0,1-Jet Tight OO",
  "Jets01PassPtG10OE": "0,1-Jet Tight OE",
  "Jets01PassPtG10EE": "0,1-Jet Tight EE",
                        
  "Jets01FailPtG10BB": "0,1-Jet Loose BB",
  "Jets01FailPtG10BO": "0,1-Jet Loose BO",
  "Jets01FailPtG10BE": "0,1-Jet Loose BE",
  "Jets01FailPtG10OO": "0,1-Jet Loose OO",
  "Jets01FailPtG10OE": "0,1-Jet Loose OE",
  "Jets01FailPtG10EE": "0,1-Jet Loose EE",

  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",
}

def makePDFBernstein(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=3):
    debug = ""
    debug += "### makePDFBakExpMOverSq: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    rooParamList = []
    rooArgList = root.RooArgList()
    for i in range(order):
      tmpArg = root.RooRealVar(channelName+"_B"+str(i),"Bernstein Coefficient "+str(i), 0.0, 0., 1.)
      rooArgList.add(tmpArg)
      rooParamList.append(tmpArg)
  
    pdfMmumu = root.RooBernstein("bak","Bernstein Order: "+str(order),dimuonMass,rooArgList)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if False:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))

#    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    #for i in rooParamList:
    #  debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    #debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug

class OrderStudy:
  def __init__(self,catName,energyStr,dataFileNames,outPrefix):
      catName = catName[0]
      randomGenerator = root.RooRandom.randomGenerator()
      randomGenerator.SetSeed(10001)
      dataTree = root.TChain()
      for i in dataFileNames:
        dataTree.Add(i+"/outtree"+catName)
      dataTree.SetCacheSize(10000000);
      dataTree.AddBranchToCache("*");

      dimuonMass = root.RooRealVar("dimuonMass","M(#mu#mu) [GeV/c^{2}]",110.,170.)
      dimuonMass.setRange("low",110,120) # Silly ranges for old fit functionality
      dimuonMass.setRange("high",130,170)
      dimuonMass.setRange("signal",120,130)
      dimuonMass.setRange("signalfit",110,140)
      dimuonMass.setBins(60)

      # Hack to Make makePDFBakOld work
      minMassZ = 88.
      maxMassZ = 94.
      dimuonMassZ = root.RooRealVar("dimuonMass","dimuonMass",minMassZ,maxMassZ)
      dimuonMassZ = dimuonMassZ

      ### Load data
      realData = root.RooDataSet("realData"+catName+energyStr,
                                      "realData"+catName+energyStr,
                                          dataTree,root.RooArgSet(dimuonMass)
                                        )
      nData = realData.sumEntries()
      #realDataZ = root.RooDataSet("realDataZ"+catName+energyStr,
      #                                "realDataZ"+catName+energyStr,
      #                                    dataTree,root.RooArgSet(dimuonMassZ)
      #                                  )
      realDataZ=None

      realDataHist = realData.binnedClone()

      ### Make Bak Pdfs
      self.rmpList = []
      self.pdfList = []
      self.frList = []
      data = {}

      self.outStr = "################################################################"
      self.outStr += "\n\n"
      self.outStr += catName + energyStr
      self.outStr += "\n\n"
      for pdfBaseName in ["Bernstein"]:
        data[pdfBaseName] = {}
        for order in range(1,8):
          data[pdfBaseName][order] = {}
          w = root.RooWorkspace("w"+pdfBaseName+str(order))
          wImport = getattr(w,"import")
          pdfName = pdfBaseName+str(order)
          pdfFunc = globals()["makePDF"+pdfBaseName]
          pdfFunc(pdfName+catName+energyStr,realData,dimuonMass,110,170,wImport,dimuonMassZ,realDataZ,order=order)
          pdf = w.pdf("bak")
          fr = pdf.fitTo(realData,
                             root.RooFit.Save(True),
                             PRINTLEVEL
                           )

          rmp = RooModelPlotter(dimuonMass,pdf,realData,fr,
                            titleMap[catName],energyStr,lumiDict[energyStr],
                            caption2=pdfBaseName+" Order "+str(order)
                            )
          rmp.draw(outPrefix+"_"+catName+"_"+pdfBaseName+str(order))

          pdf.SetName(pdfName)

          chi2Var = pdf.createChi2(realDataHist)
          chi2 = chi2Var.getVal()
          ndf = dimuonMass.getBins() - 1  # b/c roofit normalizes
          ndf -= pdf.getParameters(realData).getSize()
          nll = fr.minNll()
          self.outStr+= "{0:15} chi2: {1:.2f} ndf: {2:.0f} nll: {3:.3g}\n".format(pdfName,chi2,ndf,nll)

          self.rmpList.append(rmp)
          self.pdfList.append(pdf)
          self.frList.append(fr)
          data[pdfBaseName][order]['chi2'] = chi2
          data[pdfBaseName][order]['ndf'] = ndf
          data[pdfBaseName][order]['nll'] = nll
          data[pdfBaseName][order]['ndfFunc'] = self.getNDF(pdfBaseName,order)


      for pdfBaseName in ["Bernstein"]:
        tableStr =  "\n\n"+pdfBaseName
        tableStr += "\n\n"+r"\begin{tabular}{|l|c|c|c|c|} \hline" + "\n"
        tableStr += r"Degree (d) & Goodness & NLL$_d$ & -2$\Delta$NLL(d+1,d) & $p_{\chi^2}$(d+1,d) \\"+" \n"
        tableStr += r" & of Fit  &  & &  \\ \hline \hline" + "\n"
        for order in sorted(data[pdfBaseName].keys()):
          tmpDat = data[pdfBaseName][order]
          tmpDatP1 = False
          if data[pdfBaseName].has_key(order+1):
            tmpDatP1 = data[pdfBaseName][order+1]
          tmpDat = data[pdfBaseName][order]
          gof = scipy.stats.chi2.sf(tmpDat['chi2'],tmpDat['ndf'])
          nll = tmpDat['nll']
          ndfFunc = tmpDat['ndfFunc']
          if tmpDatP1:
            nllP1 = tmpDatP1['nll']
            ndfFuncP1 = tmpDatP1['ndfFunc']
            deltaNdfFunc = ndfFuncP1-ndfFunc
            deltaNLL = -2.*(nllP1-nll)
            pDeltaNLL = scipy.stats.chi2.sf(deltaNLL,deltaNdfFunc)
            tableStr += "{0} & {1:.3g} & {2:.0f} & {3:.2f} & {4:.3g} ".format(order,gof,nll,deltaNLL,pDeltaNLL)+r" \\ \hline" + "\n"
          else:
            tableStr += "{0} & {1:.3g} & {2:.0f} & - & - ".format(order,gof,nll)+r" \\ \hline" + "\n"
        tableStr += r"\end{tabular}" + "\n\n"
        self.outStr += tableStr

      print self.outStr

  def getNDF(self,basename,order):
    if basename == "Bernstein":
        return order+1
    else:
        print "Error: getNDF: don't recognize function: "+basename
        sys.exit(1)
        
if __name__ == "__main__":
  canvas = root.TCanvas()
  outDir = "output/"

  categories = []

  jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
  jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

  categories += [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts]]
  categories += [["Jets01PassPtG10BO",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #categories += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  categories += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  categories += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  categories += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

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

  logFile = open(outDir+"orderStudy.log",'w')
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
  logFile.write("# {0}\n\n".format(now))
  inputPklFiles = glob.glob(outDir+"*.pkl")
  orderStudyList = []
  for category in categories:
    osy = OrderStudy(category,"8TeV",dataFns8TeV,outPrefix=outDir+"order_Shape")
    logFile.write(osy.outStr)
    orderStudyList.append(osy)
    
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
  logFile.write("\n\n# {0}\n".format(now))
  logFile.close()
  
