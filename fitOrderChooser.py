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

def getOrderToUseFromDict(dataOrds,name,mh):
  assert(dataOrds != None)
  assert(name != None)
  assert(mh != None)

  cat = None
  for key in dataOrds.keys():
    if key in name:
        cat = key
        break
  if cat == None:
    raise KeyError(name)

  data = dataOrds[cat]

  if data.has_key(mh):
    return data[mh]

  massList = data.keys()
  assert(len(massList)>1)

  massListP = []
  massListN = []
  massDiffListP = []
  massDiffListN = []

  for mass in massList:
    mdiff = mh - mass
    if mdiff < 0.:
      massDiffListN.append(-mdiff)
      massListN.append(mass)
    else:
      massDiffListP.append(mdiff)
      massListP.append(mass)
  minP = 0
  minN = 0
  if len(massListP)>0:
    minP = massListP[massDiffListP.index(min(massDiffListP))]
    minP = data[minP]
  if len(massListN)>0:
    minN = massListN[massDiffListN.index(min(massDiffListN))]
    minN = data[minN]
  result = max(minP,minN)

  return result

########################################################################################

def makePDFBakBernstein(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None,higgsMass=None):
    debug = ""
    debug += "### makePDFBakBernstein: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    # Bernstein Default Order Dict
    # For Window Width: 30.0 GeV
    defaultOrders = {
      'Jet2CutsFailVBFGF': {
        115:5,
        120:5,
        125:2,
        135:1,
        150:2,
        155:1,
      },
      'Jet2CutsGFPass': {
        115:3,
        120:3,
        125:2,
        135:1,
        150:1,
        155:1,
      },
      'Jet2CutsVBFPass': {
        115:2,
        120:2,
        125:1,
        135:1,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BB': {
        115:5,
        120:4,
        125:4,
        135:2,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BE': {
        115:6,
        120:4,
        125:3,
        135:2,
        150:2,
        155:2,
      },
      'Jets01PassPtG10BO': {
        115:6,
        120:4,
        125:4,
        135:2,
        150:2,
        155:2,
      },
    }

    if order == None:
      order = getOrderToUseFromDict(defaultOrders,name,higgsMass)

    rooParamList = []
    rooArgList = root.RooArgList()
    for i in range(order+1):
      tmpArg = root.RooRealVar(channelName+"_B"+str(i),"Bernstein Coefficient "+str(i), 0.0, 0., 1.)
      rooArgList.add(tmpArg)
      rooParamList.append(tmpArg)

    debug += "#    Bernstein Order: "+str(order)+"\n"
    debug += "#    pdfArgs: "+dimuonMass.GetName()+" "
    for i in rooParamList:
        debug += i.GetName()+" "
    debug += "\n"

    #print
    #print "Bernstein Order: ",order
    #for i in rooParamList:
    #    i.Print()
    #print
  
    pdfMmumu = root.RooBernstein("bak","Bernstein Order: "+str(order),dimuonMass,rooArgList)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
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

    return paramList, bakNormTup, debug, order

def makePDFBakChebychev(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None,higgsMass=None):
    debug = ""
    debug += "### makePDFBakChebychev: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    defaultOrders = None

    if order == None:
      order = getOrderToUseFromDict(defaultOrders,name,higgsMass)

    rooParamList = []
    rooArgList = root.RooArgList()
    for i in range(order):
      tmpArg = root.RooRealVar(channelName+"_B"+str(i),"Chebychev Coefficient "+str(i), 0.0, -1., 1.)
      rooArgList.add(tmpArg)
      rooParamList.append(tmpArg)

    debug += "#    Chebychev Order: "+str(order)+"\n"
    debug += "#    pdfArgs: "+dimuonMass.GetName()+" "
    for i in rooParamList:
        debug += i.GetName()+" "
    debug += "\n"

    print
    print "Chebychev Order: ",order
    for i in rooParamList:
        i.Print()
    print
  
    pdfMmumu = root.RooChebychev("bak","Chebychev Order: "+str(order),dimuonMass,rooArgList)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
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

    return paramList, bakNormTup, debug, order

def makePDFBakPolynomial(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None,higgsMass=None):
    debug = ""
    debug += "### makePDFBakPolynomial: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    defaultOrders = None
    if order == None:
      order = getOrderToUseFromDict(defaultOrders,name,higgsMass)

    rooParamList = []
    rooArgList = root.RooArgList()
    for i in range(order):
      tmpArg = root.RooRealVar(channelName+"_P"+str(i),"Polynomial Coefficient "+str(i), 0.0, -1., 1.)
      rooArgList.add(tmpArg)
      rooParamList.append(tmpArg)

    debug += "#    Polynomial Order: "+str(order)+"\n"
    debug += "#    pdfArgs: "+dimuonMass.GetName()+" "
    for i in rooParamList:
        debug += i.GetName()+" "
    debug += "\n"

    print
    print "Polynomial Order: ",order
    for i in rooParamList:
        i.Print()
    print
  
    pdfMmumu = root.RooPolynomial("bak","Polynomial Order: "+str(order),dimuonMass,rooArgList)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
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

    return paramList, bakNormTup, debug, order

def makePDFBakSumExp(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None,higgsMass=None):
    debug = ""
    debug += "### makePDFBakSumExp: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    # SumExp Default Order Dict
    # For Window Width: 30.0 GeV
    defaultOrders = {
      'Jet2CutsFailVBFGF': {
        115:2,
        120:2,
        125:2,
        135:1,
        150:1,
        155:1,
      },
      'Jet2CutsGFPass': {
        115:2,
        120:2,
        125:2,
        135:1,
        150:1,
        155:1,
      },
      'Jet2CutsVBFPass': {
        115:1,
        120:1,
        125:1,
        135:1,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BB': {
        115:2,
        120:2,
        125:2,
        135:2,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BE': {
        115:2,
        120:2,
        125:2,
        135:1,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BO': {
        115:3,
        120:2,
        125:2,
        135:2,
        150:1,
        155:1,
      },
    }

    if order == None:
      order = getOrderToUseFromDict(defaultOrders,name,higgsMass)

    rooParamList = []
    pyPdfList = []
    rooArgExpList = root.RooArgList()
    rooArgCoefList = root.RooArgList()
    rooExpPdfList = root.RooArgList()
    for i in range(order):
      tmpExpArg = root.RooRealVar(channelName+"_E"+str(i),"Exponential Parameter "+str(i), 0.0, -1., 0.)
      rooArgExpList.add(tmpExpArg)
      rooParamList.append(tmpExpArg)
      tmpExpPdf = root.RooExponential(channelName+"_ExpPdf"+str(i),"Exponential sub-Pdf "+str(i),dimuonMass,tmpExpArg)
      rooExpPdfList.add(tmpExpPdf)
      pyPdfList.append(tmpExpPdf)
      if i != 0:
        tmpCoefArg = root.RooRealVar(channelName+"_C"+str(i),"Exponential Coefficient "+str(i), 0.0, 0., 1.)
        rooArgCoefList.add(tmpCoefArg)
        rooParamList.append(tmpCoefArg)

    debug += "#    SumExp Order: "+str(order)+"\n"
    debug += "#    pdfArgs: "+dimuonMass.GetName()+" "
    for i in rooParamList:
        debug += i.GetName()+" "
    debug += "\n"

    print
    print "SumExp Order: ",order
    for i in rooParamList:
        i.Print()
    print
  
    pdfMmumu = root.RooAddPdf("bak","Sum of Exponentials Order: "+str(order),rooExpPdfList,rooArgCoefList)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
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

    return paramList, bakNormTup, debug, order

def makePDFBakSumPow(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None,higgsMass=None):
    debug = ""
    debug += "### makePDFBakSumPow: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    # SumPow Default Order Dict
    # For Window Width: 30.0 GeV
    defaultOrders = {
      'Jet2CutsFailVBFGF': {
        115:1,
        120:2,
        125:1,
        135:1,
        150:1,
        155:1,
      },
      'Jet2CutsGFPass': {
        115:1,
        120:1,
        125:1,
        135:1,
        150:1,
        155:1,
      },
      'Jet2CutsVBFPass': {
        115:1,
        120:1,
        125:1,
        135:1,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BB': {
        115:1,
        120:2,
        125:2,
        135:1,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BE': {
        115:1,
        120:2,
        125:2,
        135:1,
        150:1,
        155:1,
      },
      'Jets01PassPtG10BO': {
        115:1,
        120:2,
        125:2,
        135:2,
        150:1,
        155:1,
      },
    }

    if order == None:
      order = getOrderToUseFromDict(defaultOrders,name,higgsMass)

    rooParamList = []
    rooArgList = root.RooArgList(dimuonMass)
    coefIndexStrList = []
    iParam = 1
    pdfDefString = ""
    for i in range(order):
      tmpCoefArg = None
      tmpPowerArg = root.RooRealVar(channelName+"_P"+str(i),"Power Parameter "+str(i), 0.0, 0., 10.)
      rooArgList.add(tmpPowerArg)
      rooParamList.append(tmpPowerArg)
      iCoefParam = iParam
      iParam += 1
      if i != 0:
        tmpCoefArg = root.RooRealVar(channelName+"_C"+str(i),"Power Term Coefficient "+str(i), 0.0, 0., 1.)
        rooArgList.add(tmpCoefArg)
        rooParamList.append(tmpCoefArg)
        pdfDefString += "+@"+str(iParam)+"*"
        coefIndexStrList.append("@"+str(iParam))
        iParam += 1
      pdfDefString += "TMath::Power(@0,-@"+str(iCoefParam)+")"

    # now doing norm string to first term:
    normFirstTermStr = "(1"
    for iTerm in coefIndexStrList:
      normFirstTermStr+= "-"+iTerm
    normFirstTermStr += ")*"
    pdfDefString = normFirstTermStr+pdfDefString
  
    debug += "#    SumPow Order: "+str(order)+"\n"
    debug += "#    pdfDefString: "+pdfDefString+"\n"
    debug += "#    pdfArgs: "+dimuonMass.GetName()+" "
    for i in rooParamList:
        debug += i.GetName()+" "
    debug += "\n"

    print
    print "SumPow Order: ",order
    print "SumPow rooDefString: "+pdfDefString
    for i in rooParamList:
        i.Print()
    print
    rooArgList.Print()
    print

    pdfMmumu = root.RooGenericPdf("bak","Sum of Powers Order: "+str(order),pdfDefString,rooArgList)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
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

    return paramList, bakNormTup, debug, order

def makePDFBakLaurent(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None,higgsMass=None):
    debug = ""
    debug += "### makePDFBakLaurent: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    defaultOrders = None
    if order == None:
      order = getOrderToUseFromDict(defaultOrders,name,higgsMass)

    rooParamList = []
    rooTermList = root.RooArgList()
    pyTermList = []
    rooExtParList = root.RooArgList()
    for i in range(1,order+1):
      rooArgList = root.RooArgList(dimuonMass)
      tmpLPow = -4
      for j in range(1,i+1):
        tmpLPow += (-1)**j*(j-1)
      if i != 1:
        tmpCoefArg = root.RooRealVar(channelName+"_C"+str(i),"Laurent Coefficient for Power: "+str(tmpLPow)+" Temrm #: "+str(i), 0.0, 0., 1.)
        rooExtParList.add(tmpCoefArg)
        rooParamList.append(tmpCoefArg)
      pdfDefString = "TMath::Power(@0,"+str(tmpLPow)+")"
      pdfTerm = root.RooGenericPdf(name+"_LauTerm"+str(i),"Laurent Term: "+str(i),pdfDefString,rooArgList)
      pyTermList.append(pdfTerm)
      rooTermList.add(pdfTerm)

    debug += "#    Laurent Order: "+str(order)+"\n"
    debug += "#    pdfArgs: "+dimuonMass.GetName()+" "
    for i in rooParamList:
        debug += i.GetName()+" "
    debug += "\n"

    print
    print "Laurent Order: ",order
    for i in rooParamList:
        i.Print()
    print
    print "Laurent Term PDFs:"
    for i in pyTermList:
        i.Print()
    print

    pdfMmumu = root.RooAddPdf("bak","Laurent Order: "+str(order),rooTermList,rooExtParList)
    pdfMmumu.Print()

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
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

    return paramList, bakNormTup, debug, order

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

class OrderSummary(object):
  def __init__(self,order,gof,nll,nllp1,deltaNdfFunc):
    self.order = order
    self.gof = gof
    self.nll = nll
    self.nllp1 = nllp1
    self.m2DeltaNLL = -2.*(nllp1-nll)
    self.deltaNdfFunc = deltaNdfFunc
    self.pChi2 = scipy.stats.chi2.sf(self.m2DeltaNLL,self.deltaNdfFunc)

  def __str__(self):
    outStr = "{0:>4} {1:>10} {2:>10} {3:>10} {4:>10} ".format("d","GOF","NLL","-2DeltaNLL","pChi2")+"\n"
    outStr += "{0:4} {1:10.3g} {2:10.2f} {3:10.2f} {4:10.3g}".format(self.order,self.gof,-self.nll,self.m2DeltaNLL,self.pChi2)
    return outStr

class OrderStudy:
  def __init__(self,catName,energyStr,dataFileNames,outPrefix,pdfsToTry,ordersToTry,signalMassList,massWindowWidth):
      catName = catName[0]
      randomGenerator = root.RooRandom.randomGenerator()
      randomGenerator.SetSeed(10001)
      dataTree = root.TChain()
      for i in dataFileNames:
        dataTree.Add(i+"/outtree"+catName)
      dataTree.SetCacheSize(10000000);
      dataTree.AddBranchToCache("*");

      ### Make Bak Pdfs
      self.rmpList = []
      self.pdfList = []
      self.frList = []
      data = {}

      self.outStr = "################################################################"
      self.outStr += "\n\n"
      self.outStr += catName + energyStr
      self.outStr += "\n\n"
      self.latexStrDetail = ""
      self.outStrDetail = "##############################################################\n\n"
      self.outStrDetail += catName + energyStr + "\n\n"
      for pdfBaseName in pdfsToTry:
        data[pdfBaseName] = {}
        for sigMass in signalMassList:
          minMass = sigMass - massWindowWidth/2.
          maxMass = sigMass + massWindowWidth/2.

          dimuonMass = root.RooRealVar("dimuonMass","M(#mu#mu) [GeV/c^{2}]",minMass,maxMass)
          dimuonMass.setBins(int(massWindowWidth))
          ### Load data
          realData = root.RooDataSet("realData"+catName+energyStr,
                                          "realData"+catName+energyStr,
                                              dataTree,root.RooArgSet(dimuonMass)
                                            )
          nData = realData.sumEntries()
          dimuonMassZ = None
          realDataZ=None
          realDataHist = realData.binnedClone()
          data[pdfBaseName][sigMass] = {}
          for order in ordersToTry:
            data[pdfBaseName][sigMass][order] = {}
            w = root.RooWorkspace("w"+pdfBaseName+str(order))
            wImport = getattr(w,"import")
            pdfName = pdfBaseName+str(order)
            pdfFunc = globals()["makePDFBak"+pdfBaseName]
            tmpParamList,tmpNormTup,tmpDebug,tmpOrder = pdfFunc(pdfName+catName+energyStr,realData,dimuonMass,minMass,maxMass,wImport,dimuonMassZ,realDataZ,order=order)
            pdf = w.pdf("bak")
            fr = pdf.fitTo(realData, 
                               root.RooFit.Range(minMass,maxMass),
                               #root.RooFit.Hesse(True), 
                               #root.RooFit.Minos(True), # Doesn't Help, just makes it run longer
                               root.RooFit.Save(True),
                               PRINTLEVEL
                             )

            rmp = RooModelPlotter(dimuonMass,pdf,realData,fr,
                              titleMap[catName],energyStr,lumiDict[energyStr],
                              caption2=pdfBaseName+" Order "+str(order),
                              caption3="m_{{H}}={0:.0f} GeV/c^{{2}}".format(sigMass)
                              )
            rmp.draw(outPrefix+"_"+catName+"_"+pdfBaseName+"_"+str(sigMass)+"_"+str(order))
            floatParsFinal = fr.floatParsFinal()
            self.outStrDetail += "\n{0}\n".format(pdf.GetTitle())
            #self.outStrDetail += tmpDebug + "\n"
            for i in range(floatParsFinal.getSize()):
              parTitle = floatParsFinal.at(i).GetTitle()
              parVal = floatParsFinal.at(i).getVal()
              parErr = floatParsFinal.at(i).getError()
              self.outStrDetail += "  {0:30}: {1:10.3g} +/- {2:10.3g}\n".format(parTitle,parVal,parErr)

            pdf.SetName(pdfName)
          
            ndfFunc = rooPdfNFreeParams(pdf,realData)

            chi2Var = pdf.createChi2(realDataHist)
            chi2 = chi2Var.getVal()
            ndf = dimuonMass.getBins() - 1  # b/c roofit normalizes
            ndf -= ndfFunc
            nll = fr.minNll()
            #self.outStr+= "{0:15} chi2: {1:.2f} ndf: {2:.0f} nll: {3:.3g}\n".format(pdfName,chi2,ndf,nll)

            self.rmpList.append(rmp)
            self.pdfList.append(pdf)
            self.frList.append(fr)
            data[pdfBaseName][sigMass][order]['chi2'] = chi2
            data[pdfBaseName][sigMass][order]['ndf'] = ndf
            data[pdfBaseName][sigMass][order]['nll'] = nll
            data[pdfBaseName][sigMass][order]['ndfFunc'] = ndfFunc

      self.summary = {}

      self.latexStr = ""
      for pdfBaseName in pdfsToTry:
        self.summary[pdfBaseName] = {}
        for sigMass in signalMassList:
          tableStr =  catName+" "+energyStr+"\n"+pdfBaseName+" mH = {0:.1f}  massWindowWidth = {1:.2f} ".format(sigMass,massWindowWidth)
          tableStr += "\n\n"+r"\begin{tabular}{|l|c|c|c|c|} \hline" + "\n"
          tableStr += r"Order & Goodness & $-\ln\mathcal{L}$ & $-2\ln\lambda$ & $p_{\chi^2}$ \\"+" \n"
          #tableStr += r"Degree (d) & Goodness & NLL$_d$ & -2$\Delta$NLL(d+1,d) & $p_{\chi^2}$(d+1,d) \\"+" \n"
          tableStr += r" & of Fit  &  & &  \\ \hline \hline" + "\n"
          foundGoodOne = False
          tmpSummary = None
          for order in sorted(data[pdfBaseName][sigMass].keys()):
            tmpDat = data[pdfBaseName][sigMass][order]
            tmpDatP1 = False
            if data[pdfBaseName][sigMass].has_key(order+1):
              tmpDatP1 = data[pdfBaseName][sigMass][order+1]
            tmpDat = data[pdfBaseName][sigMass][order]
            gof = scipy.stats.chi2.sf(tmpDat['chi2'],tmpDat['ndf'])
            nll = tmpDat['nll']
            ndfFunc = tmpDat['ndfFunc']
            if tmpDatP1:
              nllP1 = tmpDatP1['nll']
              ndfFuncP1 = tmpDatP1['ndfFunc']
              deltaNdfFunc = ndfFuncP1-ndfFunc
              deltaNLL = -2.*(nllP1-nll)
              pDeltaNLL = scipy.stats.chi2.sf(deltaNLL,deltaNdfFunc)
              if not foundGoodOne and pDeltaNLL > 0.05:
                foundGoodOne = True
                tableStr += r"\bf {0} & \bf {1:.3g} & \bf {2:.2f} & \bf {3:.2f} & \bf {4:.3g} ".format(order,gof,-nll,deltaNLL,pDeltaNLL)+r" \\ \hline" + "\n"
                tmpSummary = OrderSummary(order,gof,nll,nllP1,deltaNdfFunc)
              else:
                tableStr += "{0} & {1:.3g} & {2:.2f} & {3:.2f} & {4:.3g} ".format(order,gof,-nll,deltaNLL,pDeltaNLL)+r" \\ \hline" + "\n"
            else:
              tableStr += "{0} & {1:.3g} & {2:.2f} & - & - ".format(order,gof,-nll)+r" \\ \hline" + "\n"
          tableStr += r"\end{tabular}" + "\n\n"
          self.latexStr += tableStr
          self.summary[pdfBaseName][sigMass] = tmpSummary
      for pdfBaseName in pdfsToTry:
        for sigMass in signalMassList:
          self.outStr +=  "\n\n"+catName+" "+energyStr+"\n"+pdfBaseName+" mH = {0:.1f}    massWindowWidth = {1:.2f}".format(sigMass,massWindowWidth)+"\n\n"
      
          self.outStr += "{0:>4} {1:>10} {2:>10} {3:>10} {4:>10} ".format("d","GOF","NLL","-2DeltaNLL","pChi2")+"\n"
          foundGoodOne = False
          for order in sorted(data[pdfBaseName][sigMass].keys()):
            tmpDat = data[pdfBaseName][sigMass][order]
            tmpDatP1 = False
            if data[pdfBaseName][sigMass].has_key(order+1):
              tmpDatP1 = data[pdfBaseName][sigMass][order+1]
            tmpDat = data[pdfBaseName][sigMass][order]
            gof = scipy.stats.chi2.sf(tmpDat['chi2'],tmpDat['ndf'])
            nll = tmpDat['nll']
            ndfFunc = tmpDat['ndfFunc']
            if tmpDatP1:
              nllP1 = tmpDatP1['nll']
              ndfFuncP1 = tmpDatP1['ndfFunc']
              deltaNdfFunc = ndfFuncP1-ndfFunc
              deltaNLL = -2.*(nllP1-nll)
              pDeltaNLL = scipy.stats.chi2.sf(deltaNLL,deltaNdfFunc)
              foundGoodStr = ""
              if not foundGoodOne and pDeltaNLL > 0.05:
                foundGoodStr = "*"
                foundGoodOne = True
              self.outStr += "{0:4} {1:10.3g} {2:10.2f} {3:10.2f} {4:10.3g} {5}".format(order,gof,-nll,deltaNLL,pDeltaNLL,foundGoodStr)+ "\n"
            else:
              self.outStr += "{0:4} {1:10.3g} {2:10.2f} {3:>10} {4:>10} ".format(order,gof,-nll,'-','-')+"\n"

      print self.outStr

  def getNDF(self,basename,order):
    if basename == "Bernstein":
        return order+1
    if basename == "Chebychev":
        return order+1
    if basename == "Polynomial":
        return order
    if basename == "SumExp":
        return 2*order
    if basename == "SumPow":
        return 2*order
    if basename == "Laurent":
        return order
    else:
        print "Error: getNDF: don't recognize function: "+basename
        sys.exit(1)

def summaryWriter(summary,windowSize):
  result = ""
  result += "#################################\n"
  result += "####### "+"{0:^17}".format("Window: {0:.1f}".format(windowSize))+" #######\n"
  result += "#################################\n"
  pdfNames = None
  categories = sorted(summary.keys())
  for category in categories:
    pdfNames =  sorted(summary[category].keys())
    break
  sigMasses = set()
  for pdfName in pdfNames:
    for category in categories:
      for sigMass in summary[category][pdfName].keys():
        if not sigMass in sigMasses:
          sigMasses.add(sigMass)
  sigMasses = sorted(list(sigMasses))
  for pdfName in pdfNames:
    result += "\n{0}\n".format(pdfName)
    result += "{0:20}".format("")
    for sigMass in sigMasses:
      result += "{0:>8.1f}".format(sigMass)
    result += "{0:>10}".format("WrstGOF")
    result += "\n"
    for category in categories:
      result += "{0:20}".format(category)
      gofList = []
      for sigMass in sigMasses:
        if orderSummaries[category][pdfName].has_key(sigMass)  and orderSummaries[category][pdfName][sigMass] != None:
          order =  orderSummaries[category][pdfName][sigMass].order
          gof =  orderSummaries[category][pdfName][sigMass].gof
          result +=  "{0:>8}".format(order)
          gofList.append(gof)
        else:
          result +=  "{0:>8}".format('-')
          gofList.append(99999999999.)
      result += "{0:>10.3f}".format(min(gofList))
      result += "\n"
  return result

def summaryLatex(summary,windowSize):
  result = ""
  pdfNames = None
  categories = sorted(summary.keys())
  for category in categories:
    pdfNames =  sorted(summary[category].keys())
    break
  sigMasses = set()
  for pdfName in pdfNames:
    for category in categories:
      for sigMass in summary[category][pdfName].keys():
        if not sigMass in sigMasses:
          sigMasses.add(sigMass)
  sigMasses = sorted(list(sigMasses))

  nSigMasses = len(sigMasses)
  nCols = nSigMasses+1
  result += "\n\n"+r"\begin{tabular}{|l|"+"c|"*nSigMasses+r"} \hline" + "\n"
  result += r"\multicolumn{"+str(nCols)+r"}{|c|}{ \bf Optimal Reference Function Order} \\ \hline"+"\n"
  result += r"% window width = {0} GeV/c^2".format(windowSize)+"\n"
  for pdfName in pdfNames:
    result += r"\multicolumn{"+str(nCols)+"}{|c|}{"+pdfName+r"} \\ \hline"+"\n"
    result += "{0:20} ".format("")
    for sigMass in sigMasses:
      result += "& {0:>8.1f} ".format(sigMass)
    result += r"\\ \hline"
    result += "\n"
    for category in categories:
      result += "{0:20} ".format(titleMap[category])
      for sigMass in sigMasses:
        if orderSummaries[category][pdfName].has_key(sigMass)  and orderSummaries[category][pdfName][sigMass] != None:
          order =  orderSummaries[category][pdfName][sigMass].order
          result +=  "& {0:>8} ".format(order)
        else:
          result +=  "& {0:>8} ".format('-')
      result += r"\\ \hline"
      result += "\n"
  result += "\end{tabular}"+"\n"
  return result

def summaryDictMaker(summary,windowSize):
  result = "\n########################################################\n"
  pdfNames = None
  categories = sorted(summary.keys())
  for category in categories:
    pdfNames =  sorted(summary[category].keys())
    break
  sigMasses = set()
  for pdfName in pdfNames:
    for category in categories:
      for sigMass in summary[category][pdfName].keys():
        if not sigMass in sigMasses:
          sigMasses.add(sigMass)
  sigMasses = sorted(list(sigMasses))
  # Base indentation
  ind = " "*4
  for pdfName in pdfNames:
    result += "\n"+ind+"# {0} Default Order Dict\n".format(pdfName)
    result += ind+"# For Window Width: {0} GeV\n".format(windowSize)
    result += ind+"defaultOrders = {\n"
    for category in categories:
      result += ind+"  '"+category+"': {\n"
      for sigMass in sigMasses:
        if orderSummaries[category][pdfName].has_key(sigMass)  and orderSummaries[category][pdfName][sigMass] != None:
          order =  orderSummaries[category][pdfName][sigMass].order
          result +=  ind+"    {0}:{1},\n".format(sigMass,order)
        else:
          result +=  ind+"    {0}:{1},\n".format(sigMass,'None')
      result += ind+"  },\n"
    result += ind+"}\n"
  return result
        
if __name__ == "__main__":
  canvas = root.TCanvas()
  outDir = "output/"

  #pdfsToTry = ["Bernstein","Chebychev","Polynomial","SumExp","SumPow","Laurent"]
  #pdfsToTry = ["SumExp","SumPow"]
  #ordersToTry= range(1,5)
  #pdfsToTry = ["Bernstein","Chebychev"]
  pdfsToTry = ["Bernstein"]
  #ordersToTry= range(1,9)
  ordersToTry= range(1,8)

  categories = []

  jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
  jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

  categoriesAll = ["BB","BO","BE","OO","OE","EE"]
  categories += [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts]]
  categories += [["Jets01PassPtG10BO",  "dimuonPt>10." +jet01PtCuts]]
  categories += [["Jets01PassPtG10BE",  "dimuonPt>10." +jet01PtCuts]]
  #categories += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #categories += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  categories += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  categories += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  categories += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]

  massWindow = 30.
  signalMasses = [115,120,125,135,150,155]

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
  logDetailFile = open(outDir+"orderStudyDetail.log",'w')
  texFile = open(outDir+"orderStudy.tex",'w')
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
  logFile.write("# {0}\n\n".format(now))
  inputPklFiles = glob.glob(outDir+"*.pkl")
  orderStudyList = []
  orderSummaries = {}
  for category in categories:
    osy = OrderStudy(category,"8TeV",dataFns8TeV,outDir+"order_Shape",pdfsToTry,ordersToTry,signalMasses,massWindow)
    logFile.write(osy.outStr)
    logFile.flush()
    logDetailFile.write(osy.outStrDetail)
    logDetailFile.flush()
    texFile.write(osy.latexStr)
    texFile.flush()
    orderStudyList.append(osy)
    orderSummaries[category[0]] = osy.summary

  logFile.write(summaryWriter(orderSummaries,massWindow))
  logFile.write(summaryDictMaker(orderSummaries,massWindow))
  texFile.write(summaryLatex(orderSummaries,massWindow))
    
  now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
  logFile.write("\n\n# {0}\n".format(now))
  logFile.close()
  logDetailFile.close()
  texFile.close()
