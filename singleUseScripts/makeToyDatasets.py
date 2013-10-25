#!/usr/bin/env python
import ROOT as root 

from singleHelpers import *
from helpers import *

root.gROOT.SetBatch(True)
#root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT
#PRINTLEVEL = root.RooFit.PrintLevel(1) #For MINUIT

############################################3

nToys = 10
mhVals = range(115,156,1)
width = 20.
outDirName = "output/"
energy = "8TeV"

nEvents = 33965 # between 100 and 200
categoryName = "Jets01PassPtG10BB"
plotTitle = "0,1-Jet Tight BB Toy"
expVal = 5.0e-3

#nEvents = 209 # between 100 and 200
#categoryName = "Jet2CutsVBFPass"
#plotTitle = "2-Jet VBF Tight Toy"
#expVal = 2.0e-2

############################################3

catVar = root.RooCategory("CMS_channel","CMS_channel")
catVar.defineType(categoryName+energy)

dimuonMassAll = root.RooRealVar("dimuonMassAll","M(#mu#mu) [GeV/c^{2}]",100,200)
dimuonMassAll.setRange("baseline",110,160)

#expParam = root.RooRealVar("expParam","expParam",-0.08)
#expPdf = root.RooExponential("expPdf","expPdf",dimuonMassAll,expParam)
#bakPdf = expPdf

InvPolMass = root.RooRealVar("InvPolMass","InvPolMass", 91.187, 30., 105.)
ExpMass = root.RooRealVar("ExpMass","ExpMass", expVal, -2., 2.)
baselinePdf = root.RooGenericPdf("bak","TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",root.RooArgList(dimuonMassAll,InvPolMass,ExpMass))
bakPdf = baselinePdf

dataWholeList = []
for iToy in range(nToys):
  dataWhole = bakPdf.generate(root.RooArgSet(dimuonMassAll),nEvents,root.RooFit.Extended())
  dataWholeList.append(dataWhole)

# Debug plots
rmpList = []
for dataWhole,iToy in zip(dataWholeList,range(nToys)):
  fr = bakPdf.fitTo(dataWhole,root.RooFit.Save(True),PRINTLEVEL)
  rmp = RooModelPlotter(dimuonMassAll,bakPdf,dataWhole,fr,
                    plotTitle,None,None,
                    caption2="Baseline-Based Toy Data",
                    caption3="<N_{Events}> = "+str(nEvents)
                    )
  rmp.draw(outDirName+"toyFit"+str(iToy+1)+"_"+categoryName+"_"+energy)
  rmpList.append(rmp)
  

for mh in mhVals:
  minMass = mh - width/2.
  maxMass = mh + width/2.
  dimuonMass = root.RooRealVar("dimuonMass","dimuonMass",minMass,maxMass)
  f = root.TFile(outDirName+"toysData_{0}_{1}_{2:.1f}".format(categoryName,energy,mh)+".root","RECREATE")
  toysDir = f.mkdir("toys")
  toysDir.cd()
  for dataWhole,iToy in zip(dataWholeList,range(nToys)):
    data_obs = root.RooDataSet("data_obs","Funky Generated Data",root.RooArgSet(dimuonMass,catVar))
    for iEvent in range(int(dataWhole.sumEntries())):
      tmpM = dataWhole.get(iEvent).getRealValue("dimuonMassAll")
      if tmpM >= minMass and tmpM <= maxMass:
        dimuonMass.setVal(tmpM)
        data_obs.add(root.RooArgSet(dimuonMass,catVar))
    data_obs.SetName("model_bData")
    data_obs.Write("toy_{0}".format(iToy+1))
  f.Close()
