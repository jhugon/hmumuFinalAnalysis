
#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TParticle.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TRandom3.h>

#include <vector>
#include <map>
#include <limits>
#include <sstream>
#include <iomanip>

using namespace std;
typedef map<string,unsigned> CandCounterType;

void addToCandCounter(CandCounterType &cc,string id)
{
    if (cc.count(id)>0)
    {
        cc[id] += 1;
    }
    else
    {
        cc[id] = 1;
    }
}

void countMultiCandEvents(string infilename="/cms/data/store/user/jhugon/hmumu/stage2/ggHmumu125_8TeV.root"){

  CandCounterType candCounts;
  CandCounterType candCounts_01JetTightBB;
  CandCounterType candCounts_01JetTightBO;
  CandCounterType candCounts_2JetVBFTight;
  CandCounterType candCounts_2JetGFTight;
  //cout << "unsigned max: " << numeric_limits<unsigned>::max() << endl;
  //cout << "unsigned long max: " << numeric_limits<unsigned long>::max() << endl;
  //cout << "unsigned long long max: " << numeric_limits<unsigned long long>::max() << endl;

  TChain * tree = new TChain("outtree");
  //string infilename = "/cms/data/store/user/jhugon/hmumu/stage2/wHmumu125_8TeV.root";
  tree->Add(infilename.c_str());
  tree->SetCacheSize(10000000);
  tree->AddBranchToCache("*");

  int eventInfo_run=0;
  int eventInfo_event=0;
  tree->SetBranchAddress("eventInfo_run",&eventInfo_run);
  tree->SetBranchAddress("eventInfo_event",&eventInfo_event);

  int muonLead_passPFRelIso=0;
  int muonSub_passPFRelIso=0;
  int muonLead_isHltMatched=0;
  int muonSub_isHltMatched=0;
  tree->SetBranchAddress("muonLead_passPFRelIso",&muonLead_passPFRelIso);
  tree->SetBranchAddress("muonSub_passPFRelIso",&muonSub_passPFRelIso);
  tree->SetBranchAddress("muonLead_isHltMatched",&muonLead_isHltMatched);
  tree->SetBranchAddress("muonSub_isHltMatched",&muonSub_isHltMatched);

  float muonLead_pt=0.;
  float muonSub_pt=0.;
  tree->SetBranchAddress("muonLead_pt",&muonLead_pt);
  tree->SetBranchAddress("muonSub_pt",&muonSub_pt);

  float dimuonMass=0.;
  float dimuonPt=0.;
  tree->SetBranchAddress("dimuonMass",&dimuonMass);
  tree->SetBranchAddress("dimuonPt",&dimuonPt);

  float jetLead_pt=0.;
  float jetSub_pt=0.;
  float ptMiss=0.;
  float deltaEtaJets=0.;
  float dijetMass=0.;
  tree->SetBranchAddress("jetLead_pt",&jetLead_pt);
  tree->SetBranchAddress("jetSub_pt",&jetSub_pt);
  tree->SetBranchAddress("ptMiss",&ptMiss);
  tree->SetBranchAddress("deltaEtaJets",&deltaEtaJets);
  tree->SetBranchAddress("dijetMass",&dijetMass);

  int eventType=0;
  tree->SetBranchAddress("eventType",&eventType);

  /*
  eventType
    if (jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)  //2-jet
    {
      if (deltaEtaJets>3.5 && dijetMass>650.)//2-jet VBF tight
      {
        addToCandCounter(candCounts_2JetVBFTight,eventID)
      }
      else if (dijetMass>250. && dimuonPt>50.)// 2-jet GF tight
      {
        addToCandCounter(candCounts_2JetGFTight,eventID)
      }
    }
  */




  unsigned selEntries = 0;

  unsigned maxEntries = numeric_limits<unsigned>::max();
  //maxEntries = 1000;
  unsigned nEntries = tree->GetEntries();
  unsigned reportEach=10000;
  if (nEntries/1000>reportEach)
    reportEach = nEntries/1000;

  ///////////////////////////////
  ///////////////////////////////
  ///////////////////////////////
  // Event Loop
  //

  for(unsigned i=0; i<nEntries;i++)
  {
    if(i >= maxEntries)
      break;

    tree->GetEntry(i);
    if (i % reportEach == 0) cout << "Entry: " << i << endl;

    // Presel cut string in python code:
    //"muonLead_passPFRelIso && muonSub_passPFRelIso && ( (muonLead_pt>25 && muonLead_isHltMatched) || (muonSub_pt >25 && muonSub_isHltMatched) )"
    bool isMuonPreselEntry = false;
    if ( muonLead_passPFRelIso && muonSub_passPFRelIso 
            && ( (muonLead_pt>25. && muonLead_isHltMatched) || (muonSub_pt >25. && muonSub_isHltMatched) )
       )
        isMuonPreselEntry = true;

    if (!isMuonPreselEntry) continue;
    if (dimuonMass < 110. || dimuonMass > 160.) continue;


    selEntries+=1;

    // event ID
    stringstream eventIDStream;
    eventIDStream << eventInfo_run << ":" << eventInfo_event;
    string eventID(eventIDStream.str());

    //cout << eventID << endl;

    addToCandCounter(candCounts,eventID);

    if (jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)  //2-jet
    {
      if (deltaEtaJets>3.5 && dijetMass>650.)//2-jet VBF tight
      {
        addToCandCounter(candCounts_2JetVBFTight,eventID);
      }
      else if (dijetMass>250. && dimuonPt>50.)// 2-jet GF tight
      {
        addToCandCounter(candCounts_2JetGFTight,eventID);
      }
    }
    else  //0,1-jet
    {
      if (dimuonPt>10.) //0,1-jet Tight
      {
        if ((16 & eventType) > 0)//0,1-jet Tight BB
        {
          addToCandCounter(candCounts_01JetTightBB,eventID);
        }
        else if ((32 & eventType) > 0)// 0,1-jet Tight BO
        {
          addToCandCounter(candCounts_01JetTightBO,eventID);
        }
      }
    }


  } // end of event loop

  CandCounterType::const_iterator candItr;
  unsigned nMulti = 0;
  unsigned nMulti_01JetTightBB = 0;
  unsigned nMulti_01JetTightBO = 0;
  unsigned nMulti_2JetVBFTight = 0;
  unsigned nMulti_2JetGFTight = 0;

  cout << "Filename: " << infilename << endl;
  //cout << "Number of cands per event" << endl;
  for(candItr = candCounts.begin(); candItr != candCounts.end(); ++candItr)
  {
    if (candItr->second > 1)
    {
      //cout << candItr->first << "    "<< candItr->second << endl;
      nMulti++;
    }
  }
  for(candItr = candCounts_01JetTightBB.begin(); candItr != candCounts_01JetTightBB.end(); ++candItr)
  {
    if (candItr->second > 1)
    {
      //cout << candItr->first << "    "<< candItr->second << endl;
      nMulti_01JetTightBB++;
    }
  }
  for(candItr = candCounts_01JetTightBO.begin(); candItr != candCounts_01JetTightBO.end(); ++candItr)
  {
    if (candItr->second > 1)
    {
      //cout << candItr->first << "    "<< candItr->second << endl;
      nMulti_01JetTightBO++;
    }
  }
  for(candItr = candCounts_2JetVBFTight.begin(); candItr != candCounts_2JetVBFTight.end(); ++candItr)
  {
    if (candItr->second > 1)
    {
      //cout << candItr->first << "    "<< candItr->second << endl;
      nMulti_2JetVBFTight++;
    }
  }
  for(candItr = candCounts_2JetGFTight.begin(); candItr != candCounts_2JetGFTight.end(); ++candItr)
  {
    if (candItr->second > 1)
    {
      //cout << candItr->first << "    "<< candItr->second << endl;
      nMulti_2JetGFTight++;
    }
  }
  cout << "Total number of entries: "<<selEntries<<endl;
  cout << "Total number of unique events: "<<candCounts.size()<<endl;
  cout << "Total number of multi-muon events: "<<nMulti<<endl;
  cout << "01JetTightBB:"<<endl;
  cout << "                unique events: "<<candCounts_01JetTightBB.size()<<endl;
  cout << "                multi-muon events: "<<nMulti_01JetTightBB<<endl;
  cout << "01JetTightBO:"<<endl;
  cout << "                unique events: "<<candCounts_01JetTightBO.size()<<endl;
  cout << "                multi-muon events: "<<nMulti_01JetTightBO<<endl;
  cout << "2JetVBFTight:"<<endl;
  cout << "                unique events: "<<candCounts_2JetVBFTight.size()<<endl;
  cout << "                multi-muon events: "<<nMulti_2JetVBFTight<<endl;
  cout << "2JetGFTight:"<<endl;
  cout << "                unique events: "<<candCounts_2JetGFTight.size()<<endl;
  cout << "                multi-muon events: "<<nMulti_2JetGFTight<<endl;
    
  
}
