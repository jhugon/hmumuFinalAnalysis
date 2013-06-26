#!/usr/bin/env python

from xsec import *
from helpers import *
import ROOT as root
import os
import sys
from math import ceil

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

  def getYield(self,cuts,error="",extraWeight=None):
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
    if extraWeight:
      cutStr = "("+cutStr+")*"+extraWeight
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

  THRESHOLD = 0.01

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
  labelsDict = {
    "Jets01PassPtG10": "Non-VBF Presel. Tight",
    "Jets01FailPtG10": "Non-VBF Presel. Loose",
    "Jet2CutsVBFPass": "VBF Presel. VBF Tight",
    "Jet2CutsGFPass": "VBF Presel. GF Tight",
    "Jet2CutsFailVBFGF": "VBF Presel. Loose"
  }
  errors = ["JES","JER","PUID","MCStat"]

  effReader = EfficiencyReader()

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
            if error == "MCStat":
              dsNameForReader = dsName.replace("Hmumu","")
              nomEff, nomErr = effReader(energy,dsNameForReader,key,mass)
              relStatErr = nomErr/nomEff
              dataDict[dsName][energy][mass][key][error] = relStatErr
              continue
            maxDiff = 0.0
            negative = False
            for eSign in ["Up","Down"]:
              e = error+eSign
              varied = 0.
              if error == "PUID":
                if eSign=="Up":
                  varied = ds.getYield(cuts,extraWeight="puidUncWeight")
                else:
                  continue
              else:
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
    tableStr = ""
    extraCols = len(datasetNames)*len(energies)
    tableStr += r"\begin{tabular}{|l|"+"c|"*extraCols+r"} \hline" + "\n"
    tableStr += r"Category &"
    for ds in datasetNames:
      dsLabel = "GF"
      if "vbf" in ds:
        dsLabel = "VBF"
      for energy in energies:
        tableStr += r" %s %s &" % (dsLabel,energy.replace("TeV"," TeV"))
    tableStr = tableStr[:-1] + r"\\ \hline \hline" + "\n"
    for cat in sorted(categoryDict.keys()):
      label = labelsDict[cat]
      tableStr += label + " &"
      for ds in datasetNames:
        for energy in energies:
          yMax = 0.
          for mass in masses:
            y = dataDict[ds][energy][mass][cat][error]
            yMax = absMax(y,yMax)
          tableStr +=  " %.2f%% &" % (yMax*100.)
      tableStr = tableStr[:-1] + r" \\ \hline "+ "\n"
    tableStr += r"\end{tabular}" + "\n"

    tableStr = tableStr.replace(r"%",r"\%")
    print
    print error
    print
    print tableStr
    print

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
  tlatex = root.TLatex()
  tlatex.SetNDC()
  tlatex.SetTextFont(root.gStyle.GetLabelFont())
  tlatex.SetTextSize(0.04)
  tlatex.SetTextAlign(12)

  for dsName in datasetNames:
    for energy in energies:
      for error in errors:
        maxErr = 0.
        graphs = []
        #leg = root.TLegend(0.75,0.66,0.9,0.9)
        leg = root.TLegend(0.5,0.66,0.9,0.9)
        leg.SetFillColor(0)
        leg.SetLineColor(0)
        for key,color in zip(sorted(categoryDict.keys()),[1,root.kRed,root.kBlue,root.kGreen,root.kOrange+4]):
            label = labelsDict[key]
            graph = root.TGraphErrors()
            graph.SetLineColor(color)
            graph.SetMarkerColor(color)
            #graph.SetMarkerStyle(0)
            iPoint = 0
            for mass in sorted(masses,key=float):
              relErr = dataDict[dsName][energy][mass][key][error]
              graph.SetPoint(iPoint,float(mass),relErr*100.)
              maxErr = max(maxErr,relErr*100.)
              #relStatErr = dataDict[dsName][energy][mass][key]["MCStat"]
              #graph.SetPointError(iPoint,0.,relStatErr*relErr*100.)
              iPoint += 1
            graph.Draw("al")
            leg.AddEntry(graph,label,"ep")
            graphs.append(graph)
        maxErr = maxErr*1.8
        if maxErr > 30.:
          maxErr = ceil(maxErr/10) * 10
        elif maxErr > 16.:
          maxErr = ceil(maxErr/5) * 5
        elif maxErr > 6.:
          maxErr = ceil(maxErr/2) * 2
        elif maxErr > 1.6:
          maxErr = ceil(maxErr)
        elif maxErr > 0.6:
          maxErr = ceil(maxErr*5) / 5.
        elif maxErr > 0.2:
          maxErr = ceil(maxErr*10) / 10.
        elif maxErr > 0.03:
          maxErr = ceil(maxErr*100) / 100.
        axes = root.TH2F("axesHist"+dsName+energy+error,"",1,110.,160.,1,0.,maxErr)
        setHistTitles(axes,"m_{H} [GeV/c^{2}]","Error on Signal Yield [%]")
        axes.Draw()
        for graph in graphs:
          graph.Draw('elp')
        leg.Draw()
        tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,PRELIMINARYSTRING)
        tlatex.SetTextAlign(32)
        errorLabel = error
        if errorLabel=="MCStat":
          errorLabel = "MC Statistics"
        tlatex.DrawLatex(1.02-gStyle.GetPadRightMargin(),0.96,errorLabel)
        tlatex.SetTextAlign(12)
        dsLabel = "GF"
        if "vbf" in dsName:
          dsLabel = "VBF"
        captionStr = r"%s %s" % (dsLabel,energy.replace("TeV"," TeV"))
        tlatex.DrawLatex(0.04+gStyle.GetPadLeftMargin(),0.88,captionStr)
        saveAs(canvas,outDir+"Syst_"+error+"_"+dsName+"_"+energy)
  
