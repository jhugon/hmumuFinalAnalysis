import ROOT as root

LUMI=20.0

scaleHiggsBy = 50.0 #See in Z'
scaleHiggsBy = 15.0 #See in Boosted Z
scaleHiggsBy = 5.0 #See in VBF
scaleHiggsBy = 1.0

xsec = {}
xsec["vbfHmumu"] = 3.338e-4 * scaleHiggsBy
xsec["vbfHmumu3"] = xsec["vbfHmumu"]
xsec["ggHmumu"] = 4.236e-3 * scaleHiggsBy
xsec["ZHmumu"] = 8.556e-5 * scaleHiggsBy
xsec["WHmumu"] = 1.512e-4 * scaleHiggsBy

xsec["Zmumujets"] = 1141.59  #mmumu>50GeV Sherpa 1141.59 pb  +/- 112.59 pb +/- 10%
xsec["ZmumujetsMgt100"] = 44.31  #mmumu>100GeV Sherpa 44.31 pb  +/- 4.40 pb +/- 10%
xsec["ZmumujetsM100_200FBjets"] = 37.35  #100GeV<mmumu<200GeV,etaj1*etaj2<0.0 Sherpa 37.35 pb  +/- 0.12 pb +/- 0.31%
xsec["ttbarMuMu"] = 2.64 # 2.64 pb; 227 +/- 3 pb for all decays TOP-12-007 W->munu BR=0.108 PDG

nEventsMap = {}
nEventsMap["vbfHmumu"] = 200000
nEventsMap["vbfHmumu3"] = 200000
nEventsMap["ggHmumu"] = 200000
nEventsMap["ZHmumu"] = 200000
nEventsMap["WHmumu"] = 200000
nEventsMap["Zmumujets"] = 2.4e6
nEventsMap["ZmumujetsMgt100"] = 3.2e6
nEventsMap["ZmumujetsM100_200FBjets"] = 3.2e6
nEventsMap["ttbarMuMu"] = 800000

backgroundList = [
"vbfHmumu3",
"ggHmumu",
#"ZHmumu",
#"WHmumu",
#"ZmumujetsM100_200FBjets"
"ZmumujetsMgt100"#,
#"Zmumujets",
#"ttbarMuMu"
]

legendEntries = {}
legendEntries["Zmumujets"] = "Z+Jets"
legendEntries["ZmumujetsMgt100"] = "Z+Jets"
legendEntries["ZmumujetsM100_200FBjets"] = "Z+Jets"
legendEntries["ttbarMuMu"] = "t#bar{t}"
legendEntries["vbfHmumu"] = "VBF H->#mu#mu"
legendEntries["vbfHmumu3"] = "VBF H->#mu#mu"
legendEntries["ggHmumu"] = "ggH->#mu#mu"
legendEntries["ZHmumu"] = "ZH, H->#mu#mu"
legendEntries["WHmumu"] = "WH, H->#mu#mu"

colors = {}
colors["Zmumujets"] = root.kOrange
colors["ZmumujetsMgt100"] = root.kOrange
colors["ZmumujetsM100_200FBjets"] = root.kOrange
colors["ttbarMuMu"] = root.kGreen+3
colors["vbfHmumu"] = root.kBlue+1
colors["vbfHmumu3"] = root.kBlue+1
colors["ggHmumu"] = root.kRed+1
colors["ZHmumu"] = root.kGreen+1
colors["WHmumu"] = root.kGreen+1

##################################################

eff = {}
eff["inc"] = {
  "ggHmumu": 0.635
}
eff["incMVAPreSelection"] = {
  "ggHmumu": 0.703
}
eff["vbf"] = {
  "vbfHmumu": 0.115
}
eff["vbfTight"] = {
  "vbfHmumu": 0.0664
}
eff["vbfLoose"] = {
  "vbfHmumu": 0.173
}
eff["vbfMVAPreSelection"] = {
  "vbfHmumu": 0.268
}
eff["pt50"] = {
  "ggHmumu": 0.185
}
eff["pt75"] = {
  "ggHmumu": 0.106
}
eff["incBZMuCuts"] = {
  "ggHmumu": 0.397
}
eff["pt50BZMuCuts"] = {
  "ggHmumu": 0.0789
}

back = {}
back["inc"] = 0.8
back["vbf"] = 8e-4
back["vbfTight"] = 1.6e-4
back["vbfLoose"] = 1.4e-3
back["pt50"] = 0.1
back["pt75"] = 4.8e-2
back["incBZMuCuts"] = 0.36
back["pt50BZMuCuts"] = 0.032
back["vbfMVAPreSelection"] = 2.0128e3*xsec["ZmumujetsMgt100"]
back["incMVAPreSelection"] = 0.171*xsec["ZmumujetsMgt100"]
