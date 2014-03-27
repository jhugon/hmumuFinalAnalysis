Locations of Files and Resources
================================

H->mumu Final Analysis Dirs
---------------------------

:archer: /raid/raid9/jhugon/CMSSW_6_1_1/finalAnalysis

:HPC: /scratch/osg/jhugon/hmumu/stats/CMSSW_6_1_1/finalAnalysis

:lxplus: /afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_6_1_1/finalAnalysis

:uftrig01b: /data/uftrig01b/jhugon/hmumu/finalAnalysis/CMSSW_5_3_3_patch3/finalAnalysis

UF Cluser RAIDs
---------------

:raid5: gainesville

:raid6: newberry

:raid7: alachua

:raid8: melrose

:raid9: archer

H->mumu Ntuples
---------------

stage1 means Uf Analyzer, stage 2 means ready for hmumuFinalAnalysis:

florida cluster and HPC:

/cms/data/store/user/jhugon/hmumu/[stage1|stage2]

On Melrose:

/raid/raid8/jhugon/higgsSamples/[stage1|stage2]

On cern AFS (Stage2):

/afs/cern.ch/work/j/jhugon/public/hmumuNtuplesLevel2/unzipped/

On uftrig01 (Stage2):

/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/

H->mumu Datacards
-----------------

:PAS:  /data/uftrig01b/digiovan/baselinePP/m110to160_pixelLumi/hmumuFinalAnalysis/ 
:PAS-Era Fine Mass Scan:   /afs/cern.ch/user/j/jhugon/public/hmumu/datacards/fineMassScanMmumu110160PixlLumi/ 

H->mumu All Categories Bias Data
--------------------------------

/raid/raid9/jhugon/biasData/

Backup in /afs/cern.ch/user/j/jhugon/work/public/hmumuBkgBias/

2014/03/28
+++++++++++++++

MSSM v. Bernstein Bias studies with a veriety of reference orders:

0,1-Jet Tight BO: 5,6,7,8
2-Jet Tight VBF & GF: 2,3,4,5
All Else: 3,4,5,6

7 & 8 TeV 0 & 3 sig injected, Bernstein with p0 fixed to 1e-6

2014/03/11
++++++++++++++++++++

8 TeV 0 sigma bias of MSSM w.r.t. 4-Bernsteins with p0 fixed to 1e-6

2014/02/10
+++++++++++

Included 0 & 3 sigma signal injection on 7 and 8 TeV. Included 
Bernstein, SumExp, V+E, SM PAS, V+1/m^2, V+Exp/m^2 as Ref functions 
and MSSM as alt function.  Used Reference Study Orders for Bernstein and SumExp.  
Bernstein had all parameters floating. 
**7TeV bias study performed with 8TeV Bernstein reference orders**

2014/01/09
+++++++++++

Created 2014/01/09 with Bernstein Alts and Bernsteins and SumExp Refs (Reference SumExp 2nd Order for all)

/raid/raid9/jhugon/biasData/old/biasDataRawAllCategories

Created 2014/01/10 with 1-4Bernstein Alts, 2Bernstein and 1SumExp Refs

/raid/raid9/jhugon/biasData/old/2JetVBFTight1SumExpRef_2kToys

