import re
import ROOT as root
import helpers

lumiDict={}
#lumiDict["8TeV"] = 19.39 #2012AARecovBCD
lumiDict["8TeV"] = 19.79 #2012ABCD 22Jan2013
lumiDict["7TeV"] = 5.05 #2011AB

#LUMI=0.807 #2012A
#LUMI=4.421 #2012B
#LUMI=0.495 #2012Cv1
#LUMI=6.311 #2012Cv2


#LUMI=2.311 #2011A
#LUMI=2.739 #2011B

brDict = helpers.readCSVXS("etc/br.csv")
ggHDict8 = helpers.readCSVXS("etc/ggH_8TeV.csv")
vbfHDict8 = helpers.readCSVXS("etc/vbfH_8TeV.csv")
wHDict8 = helpers.readCSVXS("etc/wH_8TeV.csv")
zHDict8 = helpers.readCSVXS("etc/zH_8TeV.csv")

ggHDict7 = helpers.readCSVXS("etc/ggH_7TeV.csv")
vbfHDict7 = helpers.readCSVXS("etc/vbfH_7TeV.csv")
wHDict7 = helpers.readCSVXS("etc/wH_7TeV.csv")
zHDict7 = helpers.readCSVXS("etc/zH_7TeV.csv")

class CrossSections:
  def __init__(self):
    self.data = {}
    self.br = helpers.readCSVXS("etc/br.csv")
    self.vbf = {}
    self.gg = {}
    self.wh = {}
    self.zh = {}
    for e in ["8TeV","7TeV"]:
      self.vbf[e] = helpers.readCSVXS("etc/vbfH_"+e+".csv")
      self.gg[e] = helpers.readCSVXS("etc/ggH_"+e+".csv")
      self.wh[e] = helpers.readCSVXS("etc/wH_"+e+".csv")
      self.zh[e] = helpers.readCSVXS("etc/zH_"+e+".csv")
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
      elif prodMode == "w":
        result = self.wh[energy][mass]
      elif prodMode == "z":
        result = self.zh[energy][mass]
      else:
        raise Exception("Higgs Production mode not recognized for: " + key)
      return result*self.br[mass]
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

nEventsMap["zHmumu125_8TeV"] = 10000
nEventsMap["wHmumu125_8TeV"] = 10000
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

nEventsMap["zHmumu125_7TeV"] = 10000
nEventsMap["wHmumu125_7TeV"] = 10000
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
"wHmumu125",
"zHmumu125"
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
legendEntries["DYJetsToLL"] = "DY+Jets"
legendEntries["DY2JetsToLL"] = "DY+Jets"
legendEntries["DY3JetsToLL"] = "DY+Jets"
legendEntries["DY4JetsToLL"] = "DY+Jets"
legendEntries["DYToMuMu"] = "DY->#mu#mu"
legendEntries["WJetsToLNu"] = "W\rightarrow\ell\nu+Jets"
legendEntries["ttbar"] = "t#bar{t}"
legendEntries["vbfHmumu125"] = "VBF H->#mu#mu"
legendEntries["ggHmumu125"] = "gg->H->#mu#mu"
legendEntries["vbfHmumu125ChangeEvents"] = "New VBF H->#mu#mu"
legendEntries["ggHmumu125ChangeEvents"] = "New ggH->#mu#mu"
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
colors["vbfHmumu125ChangeEvents"] = root.kBlue+3
colors["ggHmumu125ChangeEvents"] = root.kRed+2
colors["zHmumu125"] = root.kGreen+1
colors["wHmumu125"] = root.kOrange+7
colors["ttbar"] = root.kGreen-1
colors["DYToTauTau"] = root.kOrange+3 #brown
#colors["WW"] = root.kPink+9
#colors["WZ"] = root.kPink+9
#colors["ZZ"] = root.kPink+9
colors["WW"] = root.kPink+9
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

class NuisanceMap:
  def __init__(self):
    self.data = {}
    self.br = helpers.readCSVXS("etc/br.csv")
    self.vbf = {}
    self.gg = {}
    self.wh = {}
    self.zh = {}
    for e in ["8TeV","7TeV"]:
      self.vbf[e] = helpers.readCSVXS("etc/vbfH_"+e+".csv")
      self.gg[e] = helpers.readCSVXS("etc/ggH_"+e+".csv")
      self.wh[e] = helpers.readCSVXS("etc/wH_"+e+".csv")
      self.zh[e] = helpers.readCSVXS("etc/zH_"+e+".csv")
    self.lumi = {
        "14TeV" : 1.044,
        "8TeV" : 1.044,
        "7TeV" : 1.022,
        }
    self.JES = {
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
      'w' : {
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
      'z' : {
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
    self.JER = {
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
      'w' : {
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
      'z' : {
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
    self.PUID = {
      'gg' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0109,
          'Jet2CutsGFPass' : 1.0158,
          'Jet2CutsVBFPass' : 1.0368,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0142,
          'Jet2CutsGFPass' : 1.0174,
          'Jet2CutsVBFPass' : 1.0386,
          'Jets01FailPtG10' : None,
          'Jets01PassPtG10' : None,
          },
        },
      'vbf' : {
        '7TeV' : {
          'Jet2CutsFailVBFGF' : 1.0134,
          'Jet2CutsGFPass' : 1.0160,
          'Jet2CutsVBFPass' : 1.0324,
          'Jets01FailPtG10' : 1.0115,
          'Jets01PassPtG10' : 1.0108,
          },
        '8TeV' : {
          'Jet2CutsFailVBFGF' : 1.0155,
          'Jet2CutsGFPass' : 1.0168,
          'Jet2CutsVBFPass' : 1.0328,
          'Jets01FailPtG10' : 1.0129,
          'Jets01PassPtG10' : 1.0128,
          },
        },
      'w' : {
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
      'z' : {
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
      'w' : {
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
      'z' : {
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
      'w' : {
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
      'z' : {
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
    self.PU = {
      'gg' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : 1.0222,
          'Jet2CutsGFPass' : 1.0103,
          'Jet2CutsFailVBFGF' : 1.0112,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : None,
          'Jet2CutsVBFPass' : 1.0222,
          'Jet2CutsGFPass' : 1.0103,
          'Jet2CutsFailVBFGF' : 1.0112,
          },
        },
      'vbf' : {
        '8TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : 1.0207,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        '7TeV' : {
          'Jets01PassPtG10' : None,
          'Jets01FailPtG10' : 1.0207,
          'Jet2CutsVBFPass' : None,
          'Jet2CutsGFPass' : None,
          'Jet2CutsFailVBFGF' : None,
          },
        },
      'w' : {
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
      'z' : {
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
      'w' : {
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
      'z' : {
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
    self.UE = {
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
      'w' : {
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
      'z' : {
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

    # The list of systematics to be applied
    self._keys = ["xs_ggH","xs_vbfH","xs_wH","xs_zH","br_Hmm","lumi","PDF","JES","JER","PUID","MCStat","PU","UE","QCDScale"]
    # The list of systematics which are correlated between energies
    self.keysEnergyCorr = ["xs_ggH","xs_vbfH","xs_wH","xs_zH","br_Hmm","PDF","JER","UE","QCDScale"]
    # The list of systematics which are not correlated between energies or categories
    self.keysNotCatCorr = ["MCStat"]
    # The list of systematics which are not correlated between energies, but correlated w/ categories 
    # (everything else)
    self.keysNotEnergyCorr = []
    for i in self._keys:
      if not i in self.keysEnergyCorr and not i in self.keysNotCatCorr:
        self.keysNotEnergyCorr.append(i)
  def keys(self):
    return self.data.keys() + self._keys
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
    if re.match(r"^xs_.*",nu) and match:
      result = None
      if prodMode not in nu:
        return None
      if prodMode == "vbf":
        result = self.vbf[energy].lnN[mass]
      elif prodMode == "gg":
        result = self.gg[energy].lnN[mass]
      elif prodMode == "w":
        result = self.wh[energy].lnN[mass]
      elif prodMode == "z":
        result = self.zh[energy].lnN[mass]
      else:
        raise Exception("Higgs Production mode not recognized for: "+ds)
      return result
    if nu == "br_Hmm" and match:
      return self.br.lnN[mass]
    if nu == "lumi" and match:
      return self.lumi[energy]
    if nu == "PDF" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.PDF[prodMode][energy][category])
    if nu == "JES" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.JES[prodMode][energy][category])
    if nu == "JER" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.JER[prodMode][energy][category])
    if nu == "MCStat" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.MCStat[prodMode][energy][category])
    if nu == "PU" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.PU[prodMode][energy][category])
    if nu == "PUID" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.PUID[prodMode][energy][category])
    if nu == "QCDScale" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.QCDScale[prodMode][energy][category])
    if nu == "UE" and match:
      category = self.getBaseCat(category)
      return goodCorr(self.UE[prodMode][energy][category])
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
