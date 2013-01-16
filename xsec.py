import re
import ROOT as root

lumiDict={}
lumiDict["8TeV"] = 12.034 #2012ABC
lumiDict["7TeV"] = 5.05 #2011AB

#LUMI=0.807 #2012A
#LUMI=4.421 #2012B
#LUMI=0.495 #2012Cv1
#LUMI=6.311 #2012Cv2


#LUMI=2.311 #2011A
#LUMI=2.739 #2011B

scaleHiggsBy = 1.0
#scaleHiggsBy = 50.0

xsec = {}
xsec["ggHmumu125_8TeV"] = 4.294e-3 * scaleHiggsBy
xsec["vbfHmumu125_8TeV"] = 3.347e-4 * scaleHiggsBy
xsec["zHmumu125_8TeV"] = 8.675e-5 * scaleHiggsBy
xsec["wHmumu125_8TeV"] = 1.533e-4 * scaleHiggsBy

xsec["DYJetsToLL_8TeV"] = 3503.71   ## madgraph
xsec["DY2JetsToLL_8TeV"] = 181.   ## madgraph
xsec["DY3JetsToLL_8TeV"] = 51.1   ## madgraph
xsec["DY4JetsToLL_8TeV"] = 23.04  ## madgraph
xsec["ttbar_8TeV"] = 225.197   ## madgraph

xsec["DYToMuMu_8TeV"] = 5745.25/3.0   ## powheg
xsec["DYToTauTau_8TeV"] = 5745.25/3.0   ## powheg
xsec["WW_8TeV"] =   54.838
xsec["WZ_8TeV"] =   33.21 
xsec["ZZ_8TeV"] =  17.654 
xsec["WJetsToLNu_8TeV"] = 36257.2
xsec["QCD_8TeV"] =  1.346e5

xsec["ggHmumu125_7TeV"] = 3.37e-3 * scaleHiggsBy
xsec["vbfHmumu125_7TeV"] = 2.65e-4 * scaleHiggsBy
xsec["zHmumu125_7TeV"] = 6.948e-5 * scaleHiggsBy
xsec["wHmumu125_7TeV"] = 1.26e-4 * scaleHiggsBy

xsec["DYJetsToLL_7TeV"] = 3048.   ## madgraph
xsec["ttbar_7TeV"] = 157.5   ## madgraph

xsec["DYToMuMu_7TeV"] = 1666.   ## powheg
xsec["DYToTauTau_7TeV"] = 1666.   ## powheg
xsec["WW_7TeV"] =  43.
xsec["WZ_7TeV"] =  18.2
xsec["ZZ_7TeV"] =  5.9
xsec["WJetsToLNu_7TeV"] = 27770.
xsec["QCD_7TeV"] =  84679.

#LHC Higgs XS WG: European Strat Group
xsec["ggHmumu125_14TeV"] = 50.35*2.2e-4 * scaleHiggsBy
xsec["vbfHmumu125_14TeV"] = 4.172*2.2e-4 * scaleHiggsBy
xsec["zHmumu125_14TeV"] = 0.883*2.2e-4 * scaleHiggsBy
xsec["wHmumu125_14TeV"] = 1.504*2.2e-4 * scaleHiggsBy

xsec["DYJetsToLL_14TeV"] = 6131.
xsec["ttbar_14TeV"] =   964.6 #LHC Higgs XS WG: European Strat Group

nEventsMap = {}
nEventsMap["ggHmumu125_8TeV"] = 9998
nEventsMap["vbfHmumu125_8TeV"] = 9990
nEventsMap["zHmumu125_8TeV"] = 10000
nEventsMap["wHmumu125_8TeV"] = 10000
nEventsMap["DYJetsToLL_8TeV"] = 29659503
nEventsMap["DY2JetsToLL_8TeV"] = 0.000001
nEventsMap["DY3JetsToLL_8TeV"] = 11015445
nEventsMap["DY4JetsToLL_8TeV"] = 0.000001
nEventsMap["ttbar_8TeV"] = 6923750
nEventsMap["DYToMuMu_8TeV"] = 48719386
nEventsMap["DYToTauTau_8TeV"] = 3295238
nEventsMap["WW_8TeV"] =  10000431
nEventsMap["WZ_8TeV"] = 9900283
nEventsMap["ZZ_8TeV"] = 9799908
nEventsMap["WJetsToLNu_8TeV"] =  55509905
nEventsMap["QCD_8TeV"] =  21384602

nEventsMap["ggHmumu125_7TeV"] = 10000
nEventsMap["vbfHmumu125_7TeV"] = 9994
nEventsMap["zHmumu125_7TeV"] = 10000
nEventsMap["wHmumu125_7TeV"] = 10000
nEventsMap["DYJetsToLL_7TeV"] = 36264432
nEventsMap["ttbar_7TeV"] = 25645835
nEventsMap["DYToMuMu_7TeV"] = 29243564
nEventsMap["DYToTauTau_7TeV"] = 13048745
nEventsMap["WW_7TeV"] =  4225916
nEventsMap["WZ_7TeV"] = 4265243
nEventsMap["ZZ_7TeV"] = 3991045
nEventsMap["WJetsToLNu_7TeV"] =  0.00000001
nEventsMap["QCD_7TeV"] =  0.00000001

nEventsMap["ggHmumu125_14TeV"] = nEventsMap["ggHmumu125_8TeV"]
nEventsMap["vbfHmumu125_14TeV"] = nEventsMap["vbfHmumu125_8TeV"]
nEventsMap["zHmumu125_14TeV"] = nEventsMap["zHmumu125_8TeV"]
nEventsMap["wHmumu125_14TeV"] = nEventsMap["wHmumu125_8TeV"]
nEventsMap["DYJetsToLL_14TeV"] = nEventsMap["DYJetsToLL_8TeV"]
nEventsMap["ttbar_14TeV"] = nEventsMap["ttbar_8TeV"]

backgroundList = [
#"DYToMuMu",
"DYJetsToLL",
#"WJetsToLNu",
"ttbar",
"WW",
"WZ",
"ZZ"
#"QCD"
]

signalList = [
"ggHmumu125",
"vbfHmumu125",
"wHmumu125",
"zHmumu125"
]

dataDict = {}

dataDict["8TeV"] = [
"SingleMuRun2012Av1",
"SingleMuRun2012Bv1",
"SingleMuRun2012Cv1",
"SingleMuRun2012Cv2"
]

dataDict["7TeV"] = [
"SingleMuRun2011Av1",
"SingleMuRun2011Bv1"
]

legendEntries = {}
legendEntries["DYJetsToLL"] = "DY+Jets"
legendEntries["DY2JetsToLL"] = "DY+Jets"
legendEntries["DY3JetsToLL"] = "DY+Jets"
legendEntries["DY4JetsToLL"] = "DY+Jets"
legendEntries["DYToMuMu"] = "DY->#mu#mu"
legendEntries["WJetsToLNu"] = "W\rightarrow\ell\nu+Jets"
legendEntries["ttbar"] = "t#bar{t}"
legendEntries["vbfHmumu125"] = "VBF H->#mu#mu"
legendEntries["ggHmumu125"] = "ggH->#mu#mu"
legendEntries["zHmumu125"] = "ZH, H->#mu#mu"
legendEntries["wHmumu125"] = "WH, H->#mu#mu"
legendEntries["ttbar"] = "t#bar{t}"
legendEntries["DYToTauTau"] = "DY->#tau#tau"
legendEntries["WW"] = "VV"
legendEntries["WZ"] = "VV"
legendEntries["ZZ"] = "VV"
legendEntries["QCD"] = "QCD"

legendEntries["7TeV"] = "CMS Data 2011"
legendEntries["8TeV"] = "CMS Data 2012"

colors = {}
colors["DYJetsToLL"] = root.kOrange
colors["DY2JetsToLL"] = root.kOrange
colors["DY3JetsToLL"] = root.kOrange
colors["DY4JetsToLL"] = root.kOrange
colors["DYToMuMu"] = root.kOrange
colors["WJetsToLNu"] = root.kCyan
colors["vbfHmumu125"] = root.kBlue
colors["ggHmumu125"] = root.kRed
colors["zHmumu125"] = root.kGreen+1
colors["wHmumu125"] = root.kOrange+7
colors["ttbar"] = root.kGreen-1
colors["DYToTauTau"] = root.kOrange+3 #brown
colors["WW"] = root.kPink+9
colors["WZ"] = root.kPink+9
colors["ZZ"] = root.kPink+9
colors["QCD"] = root.kSpring+8

efficiencyMap = {}
efficiencyMap["7TeV"] = 1.0
efficiencyMap["8TeV"] = 1.0
efficiencyMap["14TeV"] = 1.0

# Data/MC scale factors
mcPlotScaleFactorMap = {}
mcPlotScaleFactorMap["7TeV"] = 0.981 #IncPresel
mcPlotScaleFactorMap["8TeV"] = 1.0 #IncPresel
mcPlotScaleFactorMap["14TeV"] = 1.0

nuisanceMap = {}
nuisanceMap["lumi"] = {
  "vbfHmumu125_8TeV":0.044,
  "ggHmumu125_8TeV":0.044,
  "wHmumu125_8TeV":0.044,
  "zHmumu125_8TeV":0.044,

  "vbfHmumu125_7TeV":0.022,
  "ggHmumu125_7TeV":0.022,
  "wHmumu125_7TeV":0.022,
  "zHmumu125_7TeV":0.022,

  "vbfHmumu125_14TeV":0.044,
  "ggHmumu125_14TeV":0.044,
  "wHmumu125_14TeV":0.044,
  "zHmumu125_14TeV":0.044
}

nuisanceMap["br_Hmm"] = {
  "vbfHmumu125_8TeV":0.06,
  "ggHmumu125_8TeV":0.06,
  "wHmumu125_8TeV":0.06,
  "zHmumu125_8TeV":0.06,

  "vbfHmumu125_7TeV":0.06,
  "ggHmumu125_7TeV":0.06,
  "wHmumu125_7TeV":0.06,
  "zHmumu125_7TeV":0.06,

  "vbfHmumu125_14TeV":0.06,
  "ggHmumu125_14TeV":0.06,
  "wHmumu125_14TeV":0.06,
  "zHmumu125_14TeV":0.06,
}

nuisanceMap["xs_ggH"] = {
  "ggHmumu125_14TeV": 0.10,
  "ggHmumu125_8TeV": 0.147,
  "ggHmumu125_7TeV": 0.147
}
nuisanceMap["xs_vbfH"] = {
  "vbfHmumu125_14TeV": 0.019,
  "vbfHmumu125_8TeV": 0.03,
  "vbfHmumu125_7TeV": 0.024
}
nuisanceMap["xs_wH"] = {
  "wHmumu125_14TeV": 0.038,
  "wHmumu125_8TeV": 0.041,
  "wHmumu125_7TeV": 0.043
}
nuisanceMap["xs_zH"] = {
  "zHmumu125_14TeV": 0.046,
  "zHmumu125_8TeV": 0.051,
  "zHmumu125_7TeV": 0.051
}

nuisanceMap["muonRes"] = {
  "vbfHmumu125_8TeV":0.03,
  "ggHmumu125_8TeV":0.03,
  "wHmumu125_8TeV":0.03,
  "zHmumu125_8TeV":0.03,

  "vbfHmumu125_7TeV":0.03,
  "ggHmumu125_7TeV":0.03,
  "wHmumu125_7TeV":0.03,
  "zHmumu125_7TeV":0.03,

  "vbfHmumu125_14TeV":0.03,
  "ggHmumu125_14TeV":0.03,
  "wHmumu125_14TeV":0.03,
  "zHmumu125_14TeV":0.03,
}

def getLegendEntry(ds):
  return legendEntries[re.sub(r"_.*","",ds)]
def getColor(ds):
  return colors[re.sub(r"_.*","",ds)]
def appendPeriod(l,period):
  return [i+"_"+period for i in l]
def getPeriod(datasetName):
  match =  re.search(r"_(.*)",datasetName)
  if match:
    return match.group(1)
  else:
    return ""
##################################################

if __name__ == "__main__":
  if scaleHiggsBy != 1:  
    print("**** Higgs XSEC Scaled by Factor of: {} ****".format(scaleHiggsBy))
  print("Integrated Lumi for Datasets: [fb^-1]")
  for i in xsec: 
    print("{0:<15} {1:.3f}".format(i,nEventsMap[i]/xsec[i]/1000.))
  print("xsec/nEvent Scale Factors for Datasets: [fb/event]")
  for i in xsec: 
    print("{0:<15} {1:.3g}".format(i,xsec[i]/nEventsMap[i]*1000.))
