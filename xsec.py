import ROOT as root

LUMI=5.723

scaleHiggsBy = 1.0
scaleHiggsBy = 50.0

xsec = {}
xsec["vbfHmumu125"] = 3.347e-4 * scaleHiggsBy
xsec["ggHmumu125"] = 4.294e-3 * scaleHiggsBy
xsec["zHmumu125"] = 8.675e-5 * scaleHiggsBy
xsec["wHmumu125"] = 1.533e-4 * scaleHiggsBy

xsec["DYJetsToLL"] = 3503.71   ## madgraph
xsec["ttbar"] = 225.197   ## madgraph

xsec["DYToMuMu"] = 5745.25/3.0   ## powheg
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
nEventsMap["DYJetsToLL"] = 29659503
nEventsMap["ttbar"] = 6923750
nEventsMap["DYToMuMu"] = 48719386
nEventsMap["DYToTauTau"] = 3295238
nEventsMap["WW"] =  10000431
nEventsMap["WZ"] = 9900283
nEventsMap["ZZ"] = 9799908
nEventsMap["WJetsToLNu"] =  55509905
nEventsMap["QCD"] =  21384602

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

dataList = [

"SingleMuRun2012Av1",
"SingleMuRun2012Bv1",
"SingleMuRun2012Cv1"
#"SingleMuRun2012Cv2"

#"DoubleMuRun2012Av1",
#"DoubleMuRun2012Bv1",
#"DoubleMuRun2012Cv1",
#"DoubleMuRun2012Cv2"
]

legendEntries = {}
legendEntries["DYJetsToLL"] = "DY+Jets"
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

colors = {}
colors["DYJetsToLL"] = root.kOrange
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

##################################################

if __name__ == "__main__":
  if scaleHiggsBy != 1:  
    print("**** Higgs XSEC Scaled by Factor of: {} ****".format(scaleHiggsBy))
  print("Integrated Lumi for Datasets:")
  for i in xsec: 
    print("{0:<15} {1:.3f}".format(i,nEventsMap[i]/xsec[i]/1000.))
  print("xsec/nEvent Scale Factors for Datasets:")
  for i in xsec: 
    print("{0:<15} {1:.3g}".format(i,xsec[i]/nEventsMap[i]*1000.))
