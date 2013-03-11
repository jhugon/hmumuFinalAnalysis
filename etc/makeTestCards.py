
import ROOT as root

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

def makeExp(channel,name,mass,param,wImport,N=None):
  rParam = root.RooRealVar(channel+"_bakParam",channel+"_bakParam",param,-1.0,0.01)
  pdf = root.RooExponential(name+"_PDF",name+"_PDF",mass,rParam)
  countParam = root.RooRealVar(channel+"_bakCount",channel+"_bakCount",1.0)
  pdfExt = root.RooExtendPdf(name,name,pdf,countParam)
  wImport(pdfExt)

def makeGaus(channel,name,mass,mean,width,wImport):
  mParam = root.RooRealVar(channel+"_meanParam",channel+"_meanParam",mean)
  wParam = root.RooRealVar(channel+"_widthParam",channel+"_widthParam",width)
  pdf = root.RooGaussian(name,name,mass,mParam,wParam)
  wImport(pdf)

def genData(workspace,mass,N,bakName,sigName=None,sigN=None):
  data = None
  if sigName != None and sigN != None:
    pdfBak = workspace.pdf(bakName)
    pdfBakExtParam = root.RooRealVar("bakN","bakN",N)
    pdfBakExt = root.RooExtendPdf("bakExt","bakExt",pdfBak,pdfBakExtParam)
    pdfSig = workspace.pdf(sigName)
    pdfSigExtParam = root.RooRealVar("sigN","sigN",sigN)
    pdfSigExt = root.RooExtendPdf("sigExt","sigExt",pdfSig,pdfSigExtParam)

    pdf = root.RooAddPdf("pdfForGenData","pdfForGenData",root.RooArgList(pdfBakExt,pdfSigExt))
    data = pdf.generateBinned(root.RooArgSet(mass),int(N+sigN))
  else:
    pdf = workspace.pdf(bakName)
    data = pdf.generateBinned(root.RooArgSet(mass),int(N))
    #data = pdf.generate(root.RooArgSet(mass),int(N))
  data.SetName("data_obs")
  getattr(workspace,"import")(data)

def makeTestRoo(filename,injectSignalStrength):
  rFile = root.TFile(filename+".root","RECREATE")
  nBak1 = 10000
  nBak2 = 100
  nData1 = nBak1
  nData2 = nBak2
  expSignal1 = 10.0
  expSignal2 = 0.5
  nSig1 = None
  nSig2 = None
  if injectSignalStrength > 0.0:
    nSig1 = injectSignalStrength*expSignal1
    nSig2 = injectSignalStrength*expSignal2
    nData1 = nBak1+nSig1
    nData2 = nBak2+nSig2

  mass = root.RooRealVar("mass","mass",100.,200.)
  #binning = root.RooBinning(200,100,200)
  w1 = root.RooWorkspace("ch1")
  w2 = root.RooWorkspace("ch2")
  w1Import = getattr(w1,"import")
  w2Import = getattr(w2,"import")
  makeExp("ch1","bak",mass,-0.01,w1Import,nBak1)
  makeExp("ch2","bak",mass,-0.05,w2Import,nBak2)
  makeGaus("ch1","sig",mass,135.0,2.0,w1Import)
  makeGaus("ch2","sig",mass,135.0,3.0,w2Import)
  genData(w1,mass,nBak1,"bak","sig",nSig1)
  genData(w2,mass,nBak2,"bak","sig",nSig2)
  w1.Write()
  w2.Write()
  rFile.Close()
  
  textFile = open(filename+".txt",'w')
  textFile.write("imax {0}\n".format(2))
  textFile.write("jmax *\n")
  textFile.write("kmax *\n")
  textFile.write("---------------------\n")
  textFile.write("shapes * * {0}.root $CHANNEL:$PROCESS\n".format(filename))
  textFile.write("---------------------\n")
  textFile.write("bin {0} {1}\n".format("ch1","ch2"))
  textFile.write("observation {0} {1}\n".format(nData1,nData2))
  textFile.write("---------------------\n")
  textFile.write("bin {0} {0} {1} {1}\n".format("ch1","ch2"))
  textFile.write("process {0} {1} {0} {1}\n".format("sig","bak"))
  textFile.write("process {0} {1} {0} {1}\n".format(0,1))
  textFile.write("rate {0} {1} {2} {3}\n".format(expSignal1,nData1,expSignal2,nData2))
  textFile.write("---------------------\n")
  textFile.write("lumi lnN {0} - {0} -\n".format(1.044))
  textFile.write("Thry lnN {0} - {0} -\n".format(1.2))
  textFile.write("bakNorm lnU - {0} - {0}\n".format(1.8))
  textFile.close()

  textFile = open(filename+"Single.txt",'w')
  textFile.write("imax {0}\n".format(1))
  textFile.write("jmax *\n")
  textFile.write("kmax *\n")
  textFile.write("---------------------\n")
  textFile.write("shapes * * {0}.root $CHANNEL:$PROCESS\n".format(filename))
  textFile.write("---------------------\n")
  textFile.write("bin {0}\n".format("ch1"))
  textFile.write("observation {0}\n".format(nData1))
  textFile.write("---------------------\n")
  textFile.write("bin {0} {0}\n".format("ch1"))
  textFile.write("process {0} {1}\n".format("sig","bak"))
  textFile.write("process {0} {1}\n".format(0,1))
  textFile.write("rate {0} {1}\n".format(expSignal1,nData1))
  textFile.write("---------------------\n")
  textFile.write("lumi lnN {0} --\n".format(1.044))
  textFile.write("Thry lnN {0} -\n".format(1.2))
  textFile.write("bakNorm lnU - {0}\n".format(1.8))
  textFile.close()

if __name__ == "__main__":
  makeTestRoo("test",0.0)
  makeTestRoo("test5",5.0)
  makeTestRoo("test10",10.0)
  makeTestRoo("test50",50.0)
