#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import sys

GLOBALCOUNTER=0

def absMax(x,y):
  absx = abs(x)
  absy = abs(y)
  if absx > absy:
    return x
  else:
    return y

class Dataset:
  def __init__(self,filename):
    self.filename = filename
    #print "filename: {0}".format(filename)
    if filename != "":
      self.rootFile = root.TFile(filename)
      self.tree = self.rootFile.Get("outtree")
      self.tree.SetCacheSize(10000000);
      self.tree.AddBranchToCache("*");
    self.datasetName = os.path.basename(filename)
    self.datasetName = self.datasetName.replace(".root","")

  def isZombie(self):
    return self.rootFile.IsZombie()

  def getYield(self,cuts,error=""):
    """
    cuts is the cut string
    error is the string for the error, including the 'Up' or 'Down'
    """
    global GLOBALCOUNTER
    #print("In datasetName: %s, loading cuts: %s error: %s" % (self.datasetName,cuts, error))
    tmpHistName = "effHist"+str(GLOBALCOUNTER)
    GLOBALCOUNTER += 1
    varToDraw = "dimuonMass"
    tmpCUTS = cuts

    if error != "":
      cutList = re.split('\W+', tmpCUTS)
      varList = re.split('\W+', varToDraw)
      cutSet = set()
      varSet = set()
      for i in cutList:
        if (not re.match(r'\d+',i)) and (not i==''):
           if not (i in cutSet):
             cutSet.add(i)
      for i in varList:
        if (not re.match(r'\d+',i)) and (not i==''):
           if not (i in varSet):
             varSet.add(i)
      cutList = list(cutSet)
      varList = list(varSet)
      newCutList = [i + "_"+error for i in cutList]
      newVarList = [i + "_"+error for i in varList]
      branchNameList = [i.GetName() for i in self.tree.GetListOfBranches()]
      for i in reversed(range(len(newCutList))):
        newVarTmp = newCutList[i]
        if not (newVarTmp in branchNameList):
           newCutList.pop(i)
           cutList.pop(i)
      for i in reversed(range(len(newVarList))):
        newVarTmp = newVarList[i]
        if not (newVarTmp in branchNameList):
           newVarList.pop(i)
           varList.pop(i)

      for i,iNew in zip(cutList,newCutList):
        tmpCUTS = tmpCUTS.replace(i,iNew)
      for i,iNew in zip(varList,newVarList):
        varToDraw = varToDraw.replace(i,iNew)

    drawStr = varToDraw+" >> "+tmpHistName
    #print drawStr
    cutStr = treeCut("",tmpCUTS)
    #print cutStr
    self.tree.Draw(drawStr,cutStr)
    tmp = root.gDirectory.Get(tmpHistName)
    if type(tmp) != root.TH1F:
      print("Warning: In datasetName: %s, loading histogram: %s: Object type is not TH1F!!" % (self.datasetName,cuts))
      return -1.
    counts = getIntegralAll(tmp)
    return counts

if __name__ == "__main__":
  dataDir = "input/V00-01-10/"
  outDir = "output/"
  root.gROOT.SetBatch(True)
  canvas = root.TCanvas('c')

  teff = root.TEfficiency()

  THRESHOLD = 0.001

  datasetNames = [
  "ggHmumu",
  "vbfHmumu"
  ]
  masses = [
    '125',
    '115',
    '135',
    '150',
  ]
  energies=[
  '7TeV',
  '8TeV'
  ]

  massCut = "dimuonMass > 110. && dimuonMass < 170."
  jet01Cuts = " && !(jetLead_pt>40. && jetSub_pt>30. && ptMiss<40.)"
  jet2Cuts = " && jetLead_pt>40. && jetSub_pt>30. && ptMiss<40."
  categoryDict = {
    "Jets01PassPtG10": massCut+jet01Cuts+" && dimuonPt>10.",
    "Jets01FailPtG10": massCut+jet01Cuts+" && !(dimuonPt>10.)",
    "Jet2CutsVBFPass": massCut+jet2Cuts+" && dijetMass > 650. && deltaEtaJets>3.5",
    "Jet2CutsGFPass": massCut+jet2Cuts+" && !(dijetMass > 650. && deltaEtaJets>3.5) && dijetMass>250. && dimuonPt>50.",
    "Jet2CutsFailVBFGF": massCut+jet2Cuts+" && !(dijetMass > 650. && deltaEtaJets>3.5) && !(dijetMass>250. && dimuonPt>50.)"
  }
  errors = ["JES","JER"]

  dataDict = {}
  for dsName in datasetNames:
    dataDict[dsName] = {}
    for energy in energies:
      dataDict[dsName][energy] = {}
#      print(dsName+" "+energy)
      for mass in masses:
        dataDict[dsName][energy][mass] = {}
#        print("  "+mass)
        filename = dataDir+dsName+mass+'_'+energy+".root"
        ds = Dataset(filename)
        for key in categoryDict:
          dataDict[dsName][energy][mass][key] = {}
          cuts = categoryDict[key]
          nominal = ds.getYield(cuts)
          for error in errors:
            maxDiff = 0.0
            negative = False
            for eSign in ["Up","Down"]:
              e = error+eSign
              varied = ds.getYield(cuts,e)
              diff = abs(varied-nominal)
              if diff > maxDiff:
                maxDiff = diff
              if eSign == "Up" and varied < nominal:
                negative = True
            relErr = maxDiff/nominal
            if negative:
              relErr *= -1.
#            print("    {0:20} Nom: {1:>7.2e} Stat Err: {2:>7.2%} {3:>10} Err: {4:<10.2%}".format(key,nominal,1./sqrt(nominal),error,relErr))
            dataDict[dsName][energy][mass][key][error] = relErr

  scriptString = ""
  
  for energy in energies:
    for dsName in datasetNames:
      print(dsName+" "+energy)
      for key in sorted(categoryDict.keys()):
        print("  "+key)
        for error in errors:
          maxErr = 0.
          for mass in masses:
            relErr = dataDict[dsName][energy][mass][key][error]
            maxErr = absMax(maxErr,relErr)
          print("    {0:20}  Err: {1:<10.2%}".format(error,maxErr))

  print "########################################\n"*3 + "\n"

  for error in errors:
    print("    self."+error+" = {")
    for dsName in datasetNames:
      dsMatch = re.match(r"([a-z]+)Hmumu.*",dsName)
      assert(dsMatch)
      print "      '"+dsMatch.group(1)+"' : {"
      for energy in energies:
        print "        '"+energy+"' : {"
        for key in sorted(categoryDict.keys()):
          maxErr = 0.
          for mass in masses:
            relErr = dataDict[dsName][energy][mass][key][error]
            maxErr = absMax(maxErr,relErr)
          if abs(maxErr) < THRESHOLD:
            print "          '"+key+"' : None,"
          else:
            if maxErr>0.:
              maxErr += 1.
            else:
              maxErr -= 1.
            print "          '"+key+"' : {0:.4f},".format(maxErr)
        print("          },")
      print("        },")
    # for wH and zH
    for dsName in datasetNames:
      dsPrint = "w"
      if "vbf" in dsName:
        dsPrint = "z"
      print "      '"+dsPrint+"' : {"
      for energy in energies:
        print "        '"+energy+"' : {"
        for key in sorted(categoryDict.keys()):
            print "          '"+key+"' : None,"
        print("          },")
      print("        },")
    print("      }")

  canvas = root.TCanvas("c1")
  #canvas.SetLogy(1)
  axes = root.TH2F("axesHist","",1,100.,160.,1,0.,10.)
  setHistTitles(axes,"m_{H} [GeV/c^{2}]","Error on Signal Yield [%]")

  graphs = []
  for dsName in datasetNames:
    for energy in energies:
      axes.Draw()
      for key,color in zip(sorted(categoryDict.keys()),[1,root.kRed,root.kBlue,root.kGreen]):
        for error,lineStyle in zip(errors,[1,2]):
          graph = root.TGraph()
          graph.SetLineColor(color)
          graph.SetLineStyle(lineStyle)
          graph.SetMarkerStyle(0)
          iPoint = 0
          for mass in sorted(masses,key=float):
            relErr = dataDict[dsName][energy][mass][key][error]
            graph.SetPoint(iPoint,float(mass),relErr*100.)
            iPoint += 1
          graph.Draw("l")
          graphs.append(graph)
  
      canvas.SaveAs("JetVarVMass_"+dsName+"_"+energy+".png")
  
