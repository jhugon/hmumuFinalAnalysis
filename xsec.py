import ROOT as root

LUMI=5.228

#scaleHiggsBy = 50.0 #See in Z'
#scaleHiggsBy = 25.0 #See in Boosted Z
#scaleHiggsBy = 15.0
#scaleHiggsBy = 5.0 #See in VBF
scaleHiggsBy = 1.0
#scaleHiggsBy = 50.0

xsec = {}
xsec["vbfHmumu125"] = 3.347e-4 * scaleHiggsBy
xsec["ggHmumu125"] = 4.294e-3 * scaleHiggsBy
xsec["zHmumu125"] = 8.675e-5 * scaleHiggsBy
xsec["wHmumu125"] = 1.533e-4 * scaleHiggsBy

xsec["DYJetsToLL"] = 3503.71   ## madgraph
xsec["ttbar"] = 225.197   ## madgraph

xsec["DYToTauTau"] = 5745.25/3.0   ## powheg
xsec["WZ"] =   33.21 
xsec["ZZ"] =  17.654 

nEventsMap = {}
nEventsMap["vbfHmumu125"] = 9990
nEventsMap["ggHmumu125"] = 9998
nEventsMap["zHmumu125"] = 10000
nEventsMap["wHmumu125"] = 10000
nEventsMap["DYJetsToLL"] = 30459503 # GPs Ntuples w/ my changes CMSSW_5_3_X
nEventsMap["ttbar"] = 6416135 # GPs Ntuples w/ my changes CMSSW_5_2_X
nEventsMap["DYToTauTau"] = 3295238
nEventsMap["WZ"] =  10000283 
nEventsMap["ZZ"] = 9799908

backgroundList = [
"DYJetsToLL",
"ttbar",
"WZ",
"ZZ"
]

signalList = [
"ggHmumu125",
"vbfHmumu125",
"wHmumu125",
"zHmumu125"
]

dataList = [
"SingleMuRun2012Av1",
"SingleMuRun2012Bv1",
#"SingleMuRun2012Cv1",
#"SingleMuRun2012Cv2"

#"DoubleMuRun2012Av1",
#"DoubleMuRun2012Bv1",
#"DoubleMuRun2012Cv1",
#"DoubleMuRun2012Cv2"
]

legendEntries = {}
legendEntries["DYJetsToLL"] = "Z+Jets"
legendEntries["ttbar"] = "t#bar{t}"
legendEntries["vbfHmumu125"] = "VBF H->#mu#mu"
legendEntries["ggHmumu125"] = "ggH->#mu#mu"
legendEntries["zHmumu125"] = "ZH, H->#mu#mu"
legendEntries["wHmumu125"] = "WH, H->#mu#mu"
legendEntries["ttbar"] = "t#bar{t}"
legendEntries["DYToTauTau"] = "DY->#tau#tau"
legendEntries["WZ"] = "WZ"
legendEntries["ZZ"] = "ZZ"

colors = {}
colors["DYJetsToLL"] = root.kOrange
colors["vbfHmumu125"] = root.kBlue
colors["ggHmumu125"] = root.kRed
colors["zHmumu125"] = root.kGreen+1
colors["wHmumu125"] = root.kGreen+1
colors["ttbar"] = root.kGreen-1
colors["DYToTauTau"] = root.kOrange+3 #brown
colors["WZ"] = root.kPink+9
colors["ZZ"] = root.kPink+9

##################################################

if __name__ == "__main__":
  print("xsec/nEvent Scale Factors for Datasets:")
  if scaleHiggsBy != 1:  
    print("**** Higgs XSEC Scaled by Factor of: {} ****".format(scaleHiggsBy))
  for i in xsec: 
    print("{0:<15} {1:.3e}".format(i,xsec[i]/nEventsMap[i]))
