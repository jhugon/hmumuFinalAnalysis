import re
import ROOT as root
import helpers

lumiDict={}
#lumiDict["8TeV"] = 19.79 #2012ABCD 22Jan2013 HF-Lumi
lumiDict["8TeV"] = 19.712 #2012ABCD 22Jan2013 Pixel-Lumi
lumiDict["7TeV"] = 5.05 #2011AB

#LUMI=0.807 #2012A
#LUMI=4.421 #2012B
#LUMI=0.495 #2012Cv1
#LUMI=6.311 #2012Cv2


#LUMI=2.311 #2011A
#LUMI=2.739 #2011B

brDict = helpers.readCSVXS("etc/brHmm.txt")
ggHDict8 = helpers.readCSVXS("etc/8TeV-ggH.txt")
vbfHDict8 = helpers.readCSVXS("etc/8TeV-vbfH.txt")
wHDict8 = helpers.readCSVXS("etc/8TeV-WH.txt")
zHDict8 = helpers.readCSVXS("etc/8TeV-ZH.txt")

ggHDict7 = helpers.readCSVXS("etc/7TeV-ggH.txt")
vbfHDict7 = helpers.readCSVXS("etc/7TeV-vbfH.txt")
wHDict7 = helpers.readCSVXS("etc/7TeV-WH.txt")
zHDict7 = helpers.readCSVXS("etc/7TeV-ZH.txt")

class CrossSections:
  def __init__(self):
    self.data = {}
    self.br = helpers.readCSVXS("etc/brHmm.txt")
    self.vbf = {}
    self.gg = {}
    self.wh = {}
    self.zh = {}
    for e in ["8TeV","7TeV"]:
      self.vbf[e] = helpers.readCSVXS("etc/"+e+"-vbfH.txt")
      self.gg[e] = helpers.readCSVXS("etc/"+e+"-ggH.txt")
      self.wh[e] = helpers.readCSVXS("etc/"+e+"-WH.txt")
      self.zh[e] = helpers.readCSVXS("etc/"+e+"-ZH.txt")
  def keys(self):
    return self.data.keys()
  def __setitem__(self,key,value):
    self.data[key] = value
  def __getitem__(self,key):
    key = key.replace("ChangeEvents","")
    match = re.match(r"([a-z]+)Hmumu([\d.]+)_([\d]+TeV)",key)
    if match:
      prodMode = match.group(1)
      mass = match.group(2)
      energy = match.group(3)
      result = -1.
      if prodMode == "vbf":
        result = self.vbf[energy][mass]
      elif prodMode == "gg":
        result = self.gg[energy][mass]
      elif prodMode == "wh":
        result = self.wh[energy][mass]
      elif prodMode == "zh":
        result = self.zh[energy][mass]
      else:
        raise Exception("Higgs Production mode not recognized for: " + key)
      return result * self.br[mass]
    else:
      return self.data[key]

xsec = CrossSections()
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
xsec["DYJetsToLL_14TeV"] = 6131.
xsec["ttbar_14TeV"] =   964.6 #LHC Higgs XS WG: European Strat Group

nEventsMap = {}

#nEventsMap["ggHmumu125_8TeV"] = 9998
#nEventsMap["vbfHmumu125_8TeV"] = 9990

# For Different than changing sample
#nEventsMap["ggHmumu125_8TeV"] = 9999
#nEventsMap["vbfHmumu125_8TeV"] = 7992
#nEventsMap["ggHmumu125ChangeEvents_8TeV"] = 9999
#nEventsMap["vbfHmumu125ChangeEvents_8TeV"] = 7992

# For different large check samples
nEventsMap["ggHmumu125_8TeV"] = 49994
nEventsMap["vbfHmumu125_8TeV"] = 49894

nEventsMap["zhHmumu125_8TeV"] = 10000
nEventsMap["whHmumu125_8TeV"] = 10000
nEventsMap["DYJetsToLL_8TeV"] = 30086987
nEventsMap["DY2JetsToLL_8TeV"] = 0.000001
nEventsMap["DY3JetsToLL_8TeV"] = 0.000001
nEventsMap["DY4JetsToLL_8TeV"] = 0.000001
nEventsMap["ttbar_8TeV"] = 6921652 
nEventsMap["DYToMuMu_8TeV"] = 48179386
nEventsMap["DYToTauTau_8TeV"] = 3295238
nEventsMap["WW_8TeV"] =  5218045
nEventsMap["WZ_8TeV"] = 10000283
nEventsMap["ZZ_8TeV"] = 10320000
nEventsMap["WJetsToLNu_8TeV"] =  44371574
nEventsMap["QCD_8TeV"] =  20764602

#nEventsMap["ggHmumu125_7TeV"] = 10000
#nEventsMap["vbfHmumu125_7TeV"] = 9994

# For different large check samples
nEventsMap["ggHmumu125_7TeV"] = 49871
nEventsMap["vbfHmumu125_7TeV"] = 42861

nEventsMap["zhHmumu125_7TeV"] = 10000
nEventsMap["whHmumu125_7TeV"] = 10000
nEventsMap["DYJetsToLL_7TeV"] = 34334400
nEventsMap["ttbar_7TeV"] = 51360015
nEventsMap["DYToMuMu_7TeV"] = 23263564
nEventsMap["DYToTauTau_7TeV"] = 19937479
nEventsMap["WW_7TeV"] =  2145916
nEventsMap["WZ_7TeV"] = 4185243
nEventsMap["ZZ_7TeV"] = 4191045
nEventsMap["WJetsToLNu_7TeV"] =  0.00000001
nEventsMap["QCD_7TeV"] =  0.00000001

#nEventsMap["ggHmumu125_14TeV"] = nEventsMap["ggHmumu125_8TeV"]
#nEventsMap["vbfHmumu125_14TeV"] = nEventsMap["vbfHmumu125_8TeV"]
#nEventsMap["zHmumu125_14TeV"] = nEventsMap["zHmumu125_8TeV"]
#nEventsMap["wHmumu125_14TeV"] = nEventsMap["wHmumu125_8TeV"]
#nEventsMap["DYJetsToLL_14TeV"] = nEventsMap["DYJetsToLL_8TeV"]
#nEventsMap["ttbar_14TeV"] = nEventsMap["ttbar_8TeV"]

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
"whHmumu125",
"zhHmumu125"
]

dataDict = {}

#dataDict["8TeV"] = [
#    "SingleMuRun2012Av1",
#    "SingleMuRun2012Av1Recover",
#    "SingleMuRun2012Bv1",
#    "SingleMuRun2012Cv1",
#    "SingleMuRun2012Cv2",
#    "SingleMuRun2012D"
#]

# Test with the new ReReco
dataDict["8TeV"] = [
    "SingleMuRun2012Av1-22Jan2013",
    "SingleMuRun2012Bv1-22Jan2013",
    "SingleMuRun2012Cv1-22Jan2013",
    "SingleMuRun2012Dv1-22Jan2013"
]

dataDict["7TeV"] = [
"SingleMuRun2011Av1",
"SingleMuRun2011Bv1"
]

legendEntries = {}
legendEntries["DYJetsToLL"] = "Drell-Yan"
legendEntries["DY2JetsToLL"] = "Drell-Yan"
legendEntries["DY3JetsToLL"] = "Drell-Yan"
legendEntries["DY4JetsToLL"] = "Drell-Yan"
legendEntries["DYToMuMu"] = "DY#rightarrow#mu#mu"
legendEntries["WJetsToLNu"] = "W#rightarrow#ell#nu+Jets"
legendEntries["ttbar"] = "t#bar{t}"
legendEntries["vbfHmumu125"] = "SM VBF H#rightarrow#mu#mu"
legendEntries["ggHmumu125"] = "SM GF H#rightarrow#mu#mu"
legendEntries["vbfHmumu125ChangeEvents"] = "New VBF H#rightarrow#mu#mu"
legendEntries["ggHmumu125ChangeEvents"] = "New GF H#rightarrow#mu#mu"
legendEntries["zhHmumu125"] = "SM ZH H#rightarrow#mu#mu"
legendEntries["whHmumu125"] = "SM WH H#rightarrow#mu#mu"
legendEntries["ttbar"] = "t#bar{t}"
legendEntries["DYToTauTau"] = "D#rightarrow#tau#tau"
legendEntries["WW"] = "VV"
legendEntries["WZ"] = "VV"
legendEntries["ZZ"] = "VV"
legendEntries["QCD"] = "QCD"

legendEntries["7TeV"] = "Data"
legendEntries["8TeV"] = "Data"

colors = {}
colors["DYJetsToLL"] = root.kOrange
colors["DY2JetsToLL"] = root.kOrange
colors["DY3JetsToLL"] = root.kOrange
colors["DY4JetsToLL"] = root.kOrange
colors["DYToMuMu"] = root.kOrange
colors["WJetsToLNu"] = root.kCyan
colors["vbfHmumu125"] = root.kBlue
colors["ggHmumu125"] = root.kRed
colors["vbfHmumu125ChangeEvents"] = root.kBlue+3
colors["ggHmumu125ChangeEvents"] = root.kRed+2
colors["zhHmumu125"] = root.kGreen+1
colors["whHmumu125"] = root.kOrange+7
colors["ttbar"] = root.kGreen-1
colors["DYToTauTau"] = root.kOrange+3 #brown
#colors["WW"] = root.kPink+9
#colors["WZ"] = root.kPink+9
#colors["ZZ"] = root.kPink+9
colors["WW"] = root.kMagenta-9
colors["WZ"] = colors["WW"]
colors["ZZ"] = colors["WW"]
colors["QCD"] = root.kSpring+8

efficiencyMap = {}
efficiencyMap["7TeV"] = 1.0
efficiencyMap["8TeV"] = 1.0
efficiencyMap["14TeV"] = 1.0

# Data/MC scale factors
mcPlotScaleFactorMap = {}
mcPlotScaleFactorMap["7TeV"] = 1.0 #IncPresel
mcPlotScaleFactorMap["8TeV"] = 1.0 #IncPresel
mcPlotScaleFactorMap["14TeV"] = 1.0

# MEKD Normalization Factors
MENormDict = {}
MENormDict['7TeV'] = {}
MENormDict['8TeV'] = {}
MENormDict['7TeV']['sigME'] = 2104
MENormDict['7TeV']['bakME'] = 41.72
MENormDict['7TeV']['sigMEPdf'] = 0.0143
MENormDict['7TeV']['bakMEPdf'] = 0.01406
MENormDict['8TeV']['sigME'] = 2104
MENormDict['8TeV']['bakME'] = 41.62
MENormDict['8TeV']['sigMEPdf'] = 0.008808
MENormDict['8TeV']['bakMEPdf'] = 0.009364

# Background Parameterization Uncertainties for makeCards.py
# in terms of number of signal events for eac category
# For group 3 reference functions (including Bernstein)
# Using H->gamma gamma bias measure: N(alt)-N(ref)
# 7TeV 2-Jet VBF Tight corrected number
BakParameterizationUncDict = {'7TeV':{},'8TeV':{}}
BakParameterizationUncDict['7TeV']['Jets01PassPtG10BB'] = 22.50
BakParameterizationUncDict['7TeV']['Jets01PassPtG10BO'] = 42.42
BakParameterizationUncDict['7TeV']['Jets01PassPtG10BE'] = 16.62
BakParameterizationUncDict['7TeV']['Jets01PassPtG10OO'] = 11.48
BakParameterizationUncDict['7TeV']['Jets01PassPtG10OE'] = 26.54
BakParameterizationUncDict['7TeV']['Jets01PassPtG10EE'] = 11.37
BakParameterizationUncDict['7TeV']['Jets01FailPtG10BB'] = 17.15
BakParameterizationUncDict['7TeV']['Jets01FailPtG10BO'] = 18.94
BakParameterizationUncDict['7TeV']['Jets01FailPtG10BE'] = 19.05
BakParameterizationUncDict['7TeV']['Jets01FailPtG10OO'] = 19.12
BakParameterizationUncDict['7TeV']['Jets01FailPtG10OE'] = 16.13
BakParameterizationUncDict['7TeV']['Jets01FailPtG10EE'] = 5.64
BakParameterizationUncDict['7TeV']['Jet2CutsVBFPass'] = 0.52
BakParameterizationUncDict['7TeV']['Jet2CutsGFPass'] = 1.74
BakParameterizationUncDict['7TeV']['Jet2CutsFailVBFGF'] = 8.36
BakParameterizationUncDict['8TeV']['Jets01PassPtG10BB'] = 40.82
BakParameterizationUncDict['8TeV']['Jets01PassPtG10BO'] = 102.24
BakParameterizationUncDict['8TeV']['Jets01PassPtG10BE'] = 63.82
BakParameterizationUncDict['8TeV']['Jets01PassPtG10OO'] = 38.99
BakParameterizationUncDict['8TeV']['Jets01PassPtG10OE'] = 151.11
BakParameterizationUncDict['8TeV']['Jets01PassPtG10EE'] = 34.22
BakParameterizationUncDict['8TeV']['Jets01FailPtG10BB'] = 40.18
BakParameterizationUncDict['8TeV']['Jets01FailPtG10BO'] = 85.46
BakParameterizationUncDict['8TeV']['Jets01FailPtG10BE'] = 74.87
BakParameterizationUncDict['8TeV']['Jets01FailPtG10OO'] = 33.24
BakParameterizationUncDict['8TeV']['Jets01FailPtG10OE'] = 78.18
BakParameterizationUncDict['8TeV']['Jets01FailPtG10EE'] = 18.87
BakParameterizationUncDict['8TeV']['Jet2CutsVBFPass'] = 1.60
BakParameterizationUncDict['8TeV']['Jet2CutsGFPass'] = 11.82
BakParameterizationUncDict['8TeV']['Jet2CutsFailVBFGF'] = 25.27

class NuisanceMap:
  def __init__(self):
    self.data = {}
    self.br = helpers.readCSVXS("etc/brHmm.txt")
    self.vbf = {}
    self.gg = {}
    self.wh = {}
    self.zh = {}
    for e in ["8TeV","7TeV"]:
      self.vbf[e] = helpers.readCSVXS("etc/"+e+"-vbfH.txt")
      self.gg[e] = helpers.readCSVXS("etc/"+e+"-ggH.txt")
      self.wh[e] = helpers.readCSVXS("etc/"+e+"-WH.txt")
      self.zh[e] = helpers.readCSVXS("etc/"+e+"-ZH.txt")
    self.lumi = {
        #"14TeV" : 1.044,    # 2012 HF-Lumi
        #"8TeV" : 1.044,    # 2012 HF-Lumi
        "14TeV" : 1.026,    # 2012 Pixel-Lumi
        "8TeV" : 1.026,    # 2012 Pixel-Lumi
        "7TeV" : 1.022,
        }
    self.CMS_scale_j = {
      'gg' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0832,
          'Jet2CutsGFPass' : 1.0585,
          'Jet2CutsVBFPass' : 1.0759,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0799,
          'Jet2CutsGFPass' : 1.0477,
          'Jet2CutsVBFPass' : 1.0606,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'vbf' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0572,
          'Jet2CutsGFPass' : 1.0381,
          'Jet2CutsVBFPass' : 1.0358,
          'Jets01FailPtG10' : -1.0321,
          'Jets01PassPtG10' : -1.0308,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0482,
          'Jet2CutsGFPass' : 1.0249,
          'Jet2CutsVBFPass' : 1.0185,
          'Jets01FailPtG10' : -1.0278,
          'Jets01PassPtG10' : -1.0250,
          },
        },
      'wh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'zh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      }
    self.CMS_res_j = {
      'gg' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : -1.0166,
          'Jet2CutsGFPass' : -1.0126,
          'Jet2CutsVBFPass' : 1.0321,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : -1.0131,
          'Jet2CutsGFPass' : -1.0101,
          'Jet2CutsVBFPass' : -1.0273,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'vbf' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : -1.0111,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : -1.0167,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : -1.0162,
          'Jets01PassPtG10' : None,
          },
        },
      'wh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'zh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      }
    self.MCStat = {
      'gg' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0301,
          'Jet2CutsGFPass' : 1.0572,
          'Jet2CutsVBFPass' : 1.1628,
          'Jets01FailPtG10' : 1.0154,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0308,
          'Jet2CutsGFPass' : 1.0577,
          'Jet2CutsVBFPass' : 1.1068,
          'Jets01FailPtG10' : 1.0147,
          'Jets01PassPtG10' : None,
          },
        },
      'vbf' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0218,
          'Jet2CutsGFPass' : 1.0177,
          'Jet2CutsVBFPass' : 1.0208,
          'Jets01FailPtG10' : 1.0628,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0178,
          'Jet2CutsGFPass' : 1.0152,
          'Jet2CutsVBFPass' : 1.0135,
          'Jets01FailPtG10' : 1.0642,
          'Jets01PassPtG10' : None,
          },
        },
      'wh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'zh' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      }
    self.PDF = {
      'gg' : {
        '8TeV' : {
          'Jets01PassPtG10' : 1.0799,
          'Jets01FailPtG10' : 1.0844,
          'Jet2CutsVBFPass' : 1.0578,
          'Jet2CutsGFPass' : 1.0715,
          'Jet2CutsFailVBFGF' : 1.0809,
          },
        '7TeV' : {
          'Jets01PassPtG10' : 1.0830,
          'Jets01FailPtG10' : 1.0868,
          'Jet2CutsVBFPass' : 1.0593,
          'Jet2CutsGFPass' : 1.0743,
          'Jet2CutsFailVBFGF' : 1.0801,
          },
        },
      'vbf' : {
        '8TeV' : {
          'Jets01PassPtG10' : 1.0197,
          'Jets01FailPtG10' : 1.0186,
          'Jet2CutsVBFPass' : 1.0378,
          'Jet2CutsGFPass' : 1.0212,
          'Jet2CutsFailVBFGF' : 1.0196,
          },
        '7TeV' : {
          'Jets01PassPtG10' : 1.0220,
          'Jets01FailPtG10' : 1.0214,
          'Jet2CutsVBFPass' : 1.0448,
          'Jet2CutsGFPass' : 1.0299,
          'Jet2CutsFailVBFGF' : 1.0204,
          },
        },
      'wh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      'zh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      }
    self.QCDScale = {
      'gg' : {
        '8TeV' : {
          'Jets01PassPtG10' : -1.0121,
          'Jets01FailPtG10' : 1.0450,
          'Jet2CutsVBFPass' : -1.1585,
          'Jet2CutsGFPass' : 1.1515,
          'Jet2CutsFailVBFGF' : -1.0662,
          },
        '7TeV' : {
          'Jets01PassPtG10' : 1.0104,
          'Jets01FailPtG10' : 1.0417,
          'Jet2CutsVBFPass' : -1.2424,
          'Jet2CutsGFPass' : -1.1608,
          'Jet2CutsFailVBFGF' : 1.0626,
          },
        },
      'vbf' : {
        '8TeV' : {
          'Jets01PassPtG10' : -1.0139,
          'Jets01FailPtG10' : -1.1598,
          'Jet2CutsVBFPass' : 1.0313,
          'Jet2CutsGFPass' : -1.0372,
          'Jet2CutsFailVBFGF' : -1.0427,
          },
        '7TeV' : {
          'Jets01PassPtG10' : 1.0101,
          'Jets01FailPtG10' : 1.0571,
          'Jet2CutsVBFPass' : -1.0489,
          'Jet2CutsGFPass' : -1.0380,
          'Jet2CutsFailVBFGF' : 1.0653,
          },
        },
      'wh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      'zh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      }
    self.UEPS = {
      'gg' : {
        '8TeV' : {
          'Jets01PassPtG10' : -1.0618,
          'Jets01FailPtG10' : 1.1698,
          'Jet2CutsVBFPass' : 1.4228,
          'Jet2CutsGFPass' : 1.2317,
          'Jet2CutsFailVBFGF' : -1.0847,
          },
        '7TeV' : {
          'Jets01PassPtG10' : -1.0658,
          'Jets01FailPtG10' : 1.1715,
          'Jet2CutsVBFPass' : 1.6606,
          'Jet2CutsGFPass' : 1.2570,
          'Jet2CutsFailVBFGF' : -1.0669,
          },
        },
      'vbf' : {
        '8TeV' : {
          'Jets01PassPtG10' : -1.0309,
          'Jets01FailPtG10' : -1.0794,
          'Jet2CutsVBFPass' : -1.1030,
          'Jet2CutsGFPass' : -1.0585,
          'Jet2CutsFailVBFGF' : 1.0432,
          },
        '7TeV' : {
          'Jets01PassPtG10' : -1.0515,
          'Jets01FailPtG10' : 1.0953,
          'Jet2CutsVBFPass' : -1.1421,
          'Jet2CutsGFPass' : -1.1014,
          'Jet2CutsFailVBFGF' : 1.0763,
          },
        },
      'wh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      'zh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      }
    self.CMS_eff_j = {
      'gg' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsFailVBFGF' : 1.0181,
          'Jet2CutsGFPass' : 1.0202,
          'Jet2CutsVBFPass' : 1.0445,
        },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsFailVBFGF' : 1.0156,
          'Jet2CutsGFPass' : 1.0189,
          'Jet2CutsVBFPass' : 1.0430,
        },
      },
      'vbf' : {
        '8TeV' : {
          'Jets01PassPtG10' : 1.0128,
          'Jets01FailPtG10' : 1.0244,
          'Jet2CutsFailVBFGF' : 1.0155,
          'Jet2CutsGFPass' : 1.0168,
          'Jet2CutsVBFPass' : 1.0328,
        },
        '7TeV' : {
          'Jets01PassPtG10' : 1.0108,
          'Jets01FailPtG10' : 1.0237,
          'Jet2CutsFailVBFGF' : 1.0134,
          'Jet2CutsGFPass' : 1.0160,
          'Jet2CutsVBFPass' : 1.0324,
        },
      },
      'wh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
        },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
        },
      },
      'zh' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
        },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsFailVBFGF' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsVBFPass' : None,
        },
      },
    }

    # Muon efficiency uncertainty on yield; same for everything (not known for VH)
    self.CMS_eff_m = 1.016

    # The list of systematics to be applied
    self._keys = ["QCDscale_ggH","QCDscale_qqH","QCDscale_VH","pdf_gg","pdf_qqbar","lumi","pdf_gg_ACCEPT","pdf_qqbar_ACCEPT","CMS_scale_j","CMS_res_j","MCStat","UEPS","QCDscale_ggH_ACCEPT","QCDscale_qqH_ACCEPT","CMS_eff_m","CMS_eff_j","br_Hmm"]
    # The list of systematics which are correlated between energies
    self.keysEnergyCorr = ["QCDscale_ggH","QCDscale_qqH","QCDscale_VH","pdf_gg","pdf_qqbar","pdf_gg_ACCEPT","pdf_qqbar_ACCEPT","CMS_scale_j","CMS_res_j","UEPS","QCDscale_ggH_ACCEPT","QCDscale_qqH_ACCEPT","CMS_eff_m","CMS_eff_j","br_Hmm"]
    # The list of systematics which are not correlated between energies or categories
    self.keysNotCatCorr = ["MCStat"]
    # The list of systematics which are not correlated between energies, but correlated w/ categories 
    # (everything else)
    self.keysNotEnergyCorr = []
    for i in self._keys:
      if not i in self.keysEnergyCorr and not i in self.keysNotCatCorr:
        self.keysNotEnergyCorr.append(i)

    self.titleMap = {
        "QCDscale_ggH": "QCD Scale",
        "QCDscale_qqH": "QCD Scale",
        "QCDscale_VH": "QCD Scale",
        "pdf_gg": "PDF",
        "pdf_qqbar": "PDF",
        "lumi": "Luminosity",
        "pdf_gg_ACCEPT": "PDF",
        "pdf_qqbar_ACCEPT": "PDF",
        "CMS_scale_j": "Jet Energy Scale",
        "CMS_res_j": "Jet Energy Resolution",
        "MCStat": "MC Statistics",
        "UEPS": "UE/PS",
        "QCDscale_ggH_ACCEPT": "QCD Scale",
        "QCDscale_qqH_ACCEPT": "QCD Scale",
        "CMS_eff_m": "Muon Efficiency",
        "CMS_eff_j": "Jet Efficiency",
        "br_Hmm": r"\BF(\hmm{})",
    }
    # The list of systematics which are correlated between energies
  def keys(self):
    return self.data.keys() + self._keys
  def getTitle(self,key):
    if self.titleMap.has_key(key):
      return self.titleMap[key]
    else:
      return key
  def __setitem__(self,key,value):
    self.data[key] = value
  def __call__(self,nu,ds,category,mass):
    if nu in self.data:
      if ds in self.data[nu]:
        return self.data[nu][ds]
    match = re.match(r"([a-z]+)Hmumu([\d.]+)_([\d]+TeV)",ds)
    if not match:
      return None
    prodMode = match.group(1)
    energy = match.group(3)
    goodCorr =  self.goodCorr
    matchQCDScale =  re.match(r"^QCDscale_([gqHV]+)$",nu)
    if matchQCDScale and match:
      result = None
      if prodMode == "vbf" and matchQCDScale.group(1)=="qqH":
        result = self.vbf[energy].getQCDScaleLnN(mass)
      elif prodMode == "gg" and matchQCDScale.group(1)=="ggH":
        result = self.gg[energy].getQCDScaleLnN(mass)
      elif prodMode == "wh" and matchQCDScale.group(1)=="VH":
        result = self.wh[energy].getQCDScaleLnN(mass)
      elif prodMode == "zh" and matchQCDScale.group(1)=="VH":
        result = self.zh[energy].getQCDScaleLnN(mass)
      #print ("Info: finding QCDscale Unc for nuisance: {0} process: {1} prodMode: {2}, Unc: {3}".format(nu,ds,prodMode,result))
      return result
    matchPDFUnc = re.match(r"^pdf_([gqbar]+)$",nu)
    if matchPDFUnc and match:
      result = None
      if prodMode == "vbf" and matchPDFUnc.group(1)=="qqbar":
        result = self.vbf[energy].getPDFUncLnN(mass)
      elif prodMode == "gg" and matchPDFUnc.group(1)=="gg":
        result = self.gg[energy].getPDFUncLnN(mass)
      elif prodMode == "wh" and matchPDFUnc.group(1)=="qqbar":
        result = self.wh[energy].getPDFUncLnN(mass)
      elif prodMode == "zh" and matchPDFUnc.group(1)=="qqbar":
        result = self.zh[energy].getPDFUncLnN(mass)
      #print ("Info: finding pdf Unc for nuisance: {0} process: {1} prodMode: {2}, Unc: {3}, group1: {4}".format(nu,ds,prodMode,result,matchPDFUnc.group(1)))
      return result
    if nu == "br_Hmm" and match:
      return self.br.getLnN(mass)
    if nu == "lumi" and match:
      return self.lumi[energy]
    if nu == "pdf_gg_ACCEPT" and match:
      category = self.getBaseCat(category)
      if prodMode != "gg":
        return None
      return goodCorr(self.PDF[prodMode][energy][category])
    if nu == "pdf_qqbar_ACCEPT" and match:
      category = self.getBaseCat(category)
      if prodMode != "vbf":
        return None
      return goodCorr(self.PDF[prodMode][energy][category])
    if nu == "CMS_scale_j" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.CMS_scale_j[prodMode][energy][category])
    if nu == "CMS_res_j" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.CMS_res_j[prodMode][energy][category])
    if nu == "MCStat" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.MCStat[prodMode][energy][category])
    if nu == "QCDscale_ggH_ACCEPT" and match:
      if prodMode != "gg":
        return None
      category = self.getBaseCat(category)
      return goodCorr(self.QCDScale[prodMode][energy][category])
    if nu == "QCDscale_qqH_ACCEPT" and match:
      if prodMode != "vbf":
        return None
      category = self.getBaseCat(category)
      return goodCorr(self.QCDScale[prodMode][energy][category])
    if nu == "UEPS" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.UEPS[prodMode][energy][category])
    if nu == "CMS_eff_m" and match:
      if prodMode == "gg" or prodMode == "vbf":
        return goodCorr(self.CMS_eff_m)
      else:
        return None
    if nu == "CMS_eff_j" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.CMS_eff_j[prodMode][energy][category])
  def getBaseCat(self,cat):
    categoriesAllCCFF = ["BB","BO","BE","OO","OE","EE","CC","FF"]
    for i in categoriesAllCCFF:
      if cat[-2:] == i:
        return cat[:-2]
    return cat
  def goodCorr(self,value):
    if value == None or value > 0.:
      return value
    return 2.0 + value
    #return abs(value)

nuisanceMap = NuisanceMap()

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

TITLEMAP = {
  "All": "Muon Selection Only",
  "Jets01": "0,1-Jet",
  "Jets01PassPtG10": "0,1-Jet Tight",
  "Jets01FailPtG10": "0,1-Jet Loose",

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

  "Jet2" : "2-Jet",
  "Jet2CutsVBFPass":"2-Jet VBF Tight",
  "Jet2CutsGFPass":"2-Jet GF Tight",
  "Jet2CutsFailVBFGF":"2-Jet Loose",
  "Jet2SplitCutsGFSplit" : "2-Jet Combination",

  "Jets01SplitCatAll": "0,1-Jet Combination",
  "CombSplitAll" : "H#rightarrow#mu#mu Combination",
}

class PdfTitleMap(object):
  def __init__(self,data):
    self.data = data
  def __getitem__(self,key):
    orderMatch = re.match(r"([\d]+)(.+)",key)
    if orderMatch:
      keyNoOrder = orderMatch.group(2)
      order = orderMatch.group(1)
      order = helpers.getOrdinalStr(order)
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
    "MOverSq":"m/(m-p_{1})^{2}",
    "Old":"Voigtian+Exp",
    "ExpMOverSq":"Exp(p_{1}m)/(m-p_{2})^{2}",
    "ExpMOverSqP0":"#frac{Exp(-p_{1}^{2}m)}{(m-p_{2})}*(#frac{1}{m-p_{2}}+p_{3}^{2}m)",
    "ExpMOverSqP0New":"e^{-p_{1}^{2}m}/(m-p_{2})^{2}+p_{3}^{2}e^{-p_{1}^{2}m}",
    "Bernstein":"Bernstein",
    "BernsteinProd":"Bernstein",
    "Chebychev":"Chebychev",
    "Polynomial":"Polynomial",
    "SumExp":"Sum of Exponentials",
    "SumPow":"Sum of Power Functions",
    "Laurent":"Laurent",
    "ExpTimesBernstein":"Exp*Bernstein",
    "ExpTimesChebychev":"Exp*Chebychev",
    "ExpTimesPolynomial":"Exp*Polynomial",
    "MSSM":"Exp#times(Breit-Wigner+1/m^{2})",
    "VoigtPMm2":"Voigtian+1/m^{2}",
    "VoigtPExpMm2":"Voigtian+Exp/m^{2}",
})

##################################################

if __name__ == "__main__":
  if scaleHiggsBy != 1:  
    print("**** Higgs XSEC Scaled by Factor of: {} ****".format(scaleHiggsBy))
  print("Integrated Lumi for Datasets: [fb^-1]")
  sortedXsec = xsec.keys()
  sortedXsec.sort()
  for i in sortedXsec: 
    print("{0:<15} {1:.3f}".format(i,nEventsMap[i]/xsec[i]/1000.))
  print("xsec/nEvent Scale Factors for Datasets: [fb/event]")
  for i in sortedXsec: 
    print("{0:<15} {1:.3g}".format(i,xsec[i]/nEventsMap[i]*1000.))
