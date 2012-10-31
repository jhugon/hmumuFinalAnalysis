import ROOT as root

LUMI=0.807

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
xsec["WW"] =   54.838
xsec["WZ"] =   33.21 
xsec["ZZ"] =  17.654 
xsec["WJetsToLNu"] = 36257.2
xsec["QCD"] =  1.346e5

nEventsMap = {}
nEventsMap["vbfHmumu125"] = 9990
nEventsMap["ggHmumu125"] = 9998
nEventsMap["zHmumu125"] = 10000
nEventsMap["wHmumu125"] = 10000
nEventsMap["DYJetsToLL"] = 30459503
nEventsMap["ttbar"] = 6923750
nEventsMap["DYToTauTau"] = 3295238
nEventsMap["WW"] =  9900431 
nEventsMap["WZ"] = 10000283
nEventsMap["ZZ"] = 9799908
nEventsMap["WJetsToLNu"] =  57330800
nEventsMap["QCD"] =  20284602

backgroundList = [
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

dataList = [
"SingleMuRun2012Av1IsoMu"

#"SingleMuRun2012Av1",
#"SingleMuRun2012Bv1"
#"SingleMuRun2012Cv1",
#"SingleMuRun2012Cv2"

#"DoubleMuRun2012Av1",
#"DoubleMuRun2012Bv1",
#"DoubleMuRun2012Cv1",
#"DoubleMuRun2012Cv2"
]

legendEntries = {}
legendEntries["DYJetsToLL"] = "DY+Jets"
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

colors = {}
colors["DYJetsToLL"] = root.kOrange
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

##################################################

if __name__ == "__main__":
  print("xsec/nEvent Scale Factors for Datasets:")
  if scaleHiggsBy != 1:  
    print("**** Higgs XSEC Scaled by Factor of: {} ****".format(scaleHiggsBy))
  for i in xsec: 
    print("{0:<15} {1:.3e}".format(i,xsec[i]/nEventsMap[i]))
