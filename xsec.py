import ROOT as root

LUMI=1.0

scaleHiggsBy = 1.0

xsec = {}
xsec["vbfHmumu"] = 3.338e-4 * scaleHiggsBy
xsec["vbfHmumu3"] = xsec["vbfHmumu"]
xsec["ggHmumu"] = 4.236e-3 * scaleHiggsBy
xsec["ZHmumu"] = 8.556e-5 * scaleHiggsBy
xsec["WHmumu"] = 1.512e-4 * scaleHiggsBy

xsec["Zmumujets"] = 1141.59  #mmumu>50GeV Sherpa 1141.59 pb  +/- 112.59 pb +/- 10%
xsec["ZmumujetsMgt100"] = 44.31  #mmumu>100GeV Sherpa 44.31 pb  +/- 4.40 pb +/- 10%

nEventsMap = {}
nEventsMap["vbfHmumu"] = 200000
nEventsMap["vbfHmumu3"] = 200000
nEventsMap["ggHmumu"] = 200000
nEventsMap["ZHmumu"] = 200000
nEventsMap["WHmumu"] = 200000
nEventsMap["Zmumujets"] = 2.4e6
nEventsMap["ZmumujetsMgt100"] = 3.2e6

backgroundList = [
"vbfHmumu3",
"ggHmumu",
"ZHmumu",
"WHmumu",
"ZmumujetsMgt100"
#"Zmumujets"
]

legendEntries = {}
legendEntries["Zmumujets"] = "Z+Jets"
legendEntries["ZmumujetsMgt100"] = "Z+Jets"
legendEntries["vbfHmumu"] = "VBF H->#mu#mu"
legendEntries["vbfHmumu3"] = "VBF H->#mu#mu"
legendEntries["ggHmumu"] = "ggH->#mu#mu"
legendEntries["ZHmumu"] = "ZH, H->#mu#mu"
legendEntries["WHmumu"] = "WH, H->#mu#mu"

colors = {}
colors["Zmumujets"] = root.kOrange
colors["ZmumujetsMgt100"] = root.kOrange
colors["vbfHmumu"] = root.kBlue+1
colors["vbfHmumu3"] = root.kBlue+1
colors["ggHmumu"] = root.kRed+1
colors["ZHmumu"] = root.kGreen+1
colors["WHmumu"] = root.kGreen+1
