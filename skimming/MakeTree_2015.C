#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>

#include <TH1.h>
#include <TH2D.h>

#include <TBranch.h>
#include <TCanvas.h>
#include "TClonesArray.h"
#include <TDirectory.h>
#include <TFile.h>
#include "TH1F.h"
#include <TLatex.h>
#include <TLegend.h>
#include "TLorentzVector.h"
#include <TMath.h>
#include "TRandom.h"
#include <TStyle.h>
#include <TSystem.h>
#include "TTree.h"
#include "TString.h"

// miscellaneous  
#include <fstream>
#include <map>
#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#endif



//this macro will run over all events in the Onia tree,
// and will store entries in a simpler way:
// - once dimuons are firing the last filter in the Trigger path X (e.g. DoubleMu0_HighQ)
// - and if the vertex probability is above 1% (by construction, you remove 1% of signal there)
//    - fill the entries with all ID variables and 
/*
Legend: Sign pairs:  +- sign=0; ++ sign=1, -- sign=2
Legend: muon pairs: Glb-Glb id=0, Glb-Trk = 1, Trk-Trk = 2

n --> 2^n -1
0  theTriggerNames.push_back("NoTrigger");
1  theTriggerNames.push_back("HLT_PAL1DoubleMuOpen_v1");
2  theTriggerNames.push_back("HLT_PAL1DoubleMu0_HighQ_v1");
3  theTriggerNames.push_back("HLT_PAL2DoubleMu3_v1");
4  theTriggerNames.push_back("HLT_PAMu3_v1");
5  theTriggerNames.push_back("HLT_PAMu7_v1");
6  theTriggerNames.push_back("HLT_PAMu12_v1");

 */


void MakeTree_2015(const char* inputOniaTree = "/afs/cern.ch/user/e/echapon/workspace/public/RunPrep2015/oniatree_upsi1S_5TeV_3.8T.root",
		   //miniUpsilon_Histos_Runs_210498-211631_HFvars_paapJune28.root",
		   //  const char* runNumber = "211739-211831",
		   //  const char* runNumber = "210498-211631",// whole pA run, with wrong alignment
		   //  const char* runNumber = "gt210658-211631", 
		   const char* runNumber = "MCtest", 
		   //const char* runNumber = "a",
		   float newVtxProbCut   = 0.01,// default in the input tree is already 0.01
		   float ptcut           = 0., // single muon cut pt
		   int nTriggerBit       = 2, 
		   // const char* dataSource= "upsiMiniTree_pp276tev_5p41_ptmu4_woPileup_june27",
		   //   const char* dataSource= "upsiMiniTree_pA5tev_ptmu4_octb15_chunk1_runGT210658",
		   const char* dataSource= "MCtest_upsi1S_5TeV_3.8T",
		   // const char* dataSource= "upsiMiniTree_pp7tev_dimu0v1_ptmu4",
		   //   const char* dataSource= "upsiMiniTree_aa276tevC50100_ppofficial_trktrk_ptmu4",
		   bool isAArereco       = false,// only for newly processed AA; old tree is with triger v1,v2 included
		   bool bAllTriggers     = true,
		   bool addExtraCentrality = true,
		   bool excludeWrongAlign_pa = false
		   )
{
  // gROOT->Macro("setTDRStyle_modified.C+");
 
  float mass_min = 7.0; 
  float mass_max = 14.0;
  int   nBins    = 30;
  
  TFile *f = TFile::Open(Form("%s",inputOniaTree));
 
  TTree *t = (TTree*)f->Get("myTree");
  Long64_t nentries = t->GetEntries();
  cout<<"Got the tree!"<<endl;
  cout << nentries << endl;
  const int NMAX=7;
  UInt_t eventNb;
  UInt_t runNb;
  Int_t Centrality;
  Int_t HLTriggers;
  Int_t Reco_QQ_size;
  Int_t Reco_QQ_type[NMAX];
  Int_t Reco_QQ_sign[NMAX];
  Int_t Reco_QQ_NtrkDeltaR03[NMAX];
  Int_t Reco_QQ_NtrkDeltaR04[NMAX];
  Int_t Reco_QQ_NtrkDeltaR05[NMAX];
  Int_t Reco_QQ_NtrkPt04[NMAX];
  Int_t Reco_QQ_NtrkPt03[NMAX];
  Int_t Reco_QQ_NtrkPt02[NMAX];

  Float_t Reco_QQ_ctau[NMAX];
  Float_t Reco_QQ_ctauErr[NMAX];
  Float_t Reco_QQ_ctauTrue[NMAX];
  Float_t Reco_QQ_VtxProb[NMAX];
  float Reco_QQ_dca[NMAX]; // new float, unused in upsilon
  Int_t Reco_QQ_trig[NMAX];
  Float_t zVtx[NMAX];

  int Reco_QQ_mupl_nTrkHits[NMAX];
  float Reco_QQ_mupl_normChi2_inner[NMAX];
  float Reco_QQ_mupl_normChi2_global[NMAX];
  float Reco_QQ_mupl_dxy[NMAX];
  float Reco_QQ_mupl_dxyErr[NMAX];
  float Reco_QQ_mupl_dz[NMAX];
  float Reco_QQ_mupl_dzErr[NMAX];
  Bool_t Reco_QQ_mupl_TrkMuArb[NMAX];
  Bool_t Reco_QQ_mupl_TMOneStaTight[NMAX]; // new bool, unused in upsilon(?)
  int Reco_QQ_mupl_nMuValHits[NMAX]; // new int, unused in upsilon
  int Reco_QQ_mupl_nPixWMea[NMAX]; // new int, unused in upsilon
  int Reco_QQ_mupl_nTrkWMea[NMAX]; // new int, unused in upsilon
  int Reco_QQ_mupl_StationsMatched[NMAX]; // new int, unused in Upsilon
  float Reco_QQ_mupl_pt_inner[NMAX]; // new float, unused.
  float Reco_QQ_mupl_pt_global[NMAX]; //new float, unused.
  float Reco_QQ_mupl_ptErr_inner[NMAX]; // new float, unused.
  float Reco_QQ_mupl_ptErr_global[NMAX]; //new float, unused.
 
  int Reco_QQ_mumi_nTrkHits[NMAX];
  float Reco_QQ_mumi_normChi2_inner[NMAX];
  float Reco_QQ_mumi_normChi2_global[NMAX];
  float Reco_QQ_mumi_dxy[NMAX];
  float Reco_QQ_mumi_dxyErr[NMAX];
  float Reco_QQ_mumi_dz[NMAX];
  float Reco_QQ_mumi_dzErr[NMAX];
  Bool_t Reco_QQ_mumi_TrkMuArb[NMAX];
  Bool_t Reco_QQ_mumi_TMOneStaTight[NMAX];  // new bool, unused in upsilon(?)
  int Reco_QQ_mumi_nMuValHits[NMAX]; // new int, unused in upsilon
  int Reco_QQ_mumi_nPixWMea[NMAX]; // new int, unused in upsilon
  int Reco_QQ_mumi_nTrkWMea[NMAX]; // new int, unused in upsilon
  int Reco_QQ_mumi_StationsMatched[NMAX]; // new int, unused in Upsilon
  float Reco_QQ_mumi_pt_inner[NMAX]; // new float, unused.
  float Reco_QQ_mumi_pt_global[NMAX]; //new float, unused.
  float Reco_QQ_mumi_ptErr_inner[NMAX]; // new float, unused.
  float Reco_QQ_mumi_ptErr_global[NMAX]; //new float, unused.

  TClonesArray *Reco_QQ_vtx = 0;;
  TClonesArray *Reco_QQ_4mom = 0;
  TClonesArray *Reco_QQ_mupl_4mom = 0;
  TClonesArray *Reco_QQ_mumi_4mom = 0;
  TClonesArray *Reco_mu_4mom = 0;

  Int_t Reco_mu_size;
  Int_t Reco_mu_type[NMAX];
  Int_t Reco_mu_charge[NMAX];

  int nPlusMu;
  int nMinusMu;
  float muPt[NMAX];
  float muEta[NMAX];
  float muPhi[NMAX];
  //Int_t nReco_QQ;
  Float_t invariantMass;
  Int_t QQtrig;
  Int_t QQsign;
  Int_t QQTrkDr03;
  Int_t QQTrkDr04;
  Int_t QQTrkDr05;

  Int_t QQTrkPt04;
  Int_t QQTrkPt03;
  Int_t QQTrkPt02;

  float weight;
  float weight2;
  Float_t dimuPt;
  Float_t dimuEta;
  Float_t dimuPhi;
  Float_t dimuRapidity;	
  Float_t vProb;
  Float_t muPlusPt;
  Float_t muMinusPt;
  Float_t muPlusEta;
  Float_t muMinusEta;
  Float_t muPlusPhi;
  Float_t muMinusPhi;
  Float_t RmuPlusPhi = 0;
  Float_t RmuMinusPhi = 0;

  // centrliaty extra stuff
  Int_t Npix, NpixelTracks,Ntracks;
  Float_t SumET_HF, SumET_HFplus, SumET_HFminus, SumET_HFplusEta4, SumET_HFminusEta4, SumET_EB, SumET_ET, SumET_EE, SumET_EEplus, SumET_EEminus, SumET_ZDC, SumET_ZDCplus, SumET_ZDCminus;

  // track extra stuff
  const int NTRKMAX=100000000;
  // int Reco_trk_size=0;
  // TClonesArray *Reco_trk_vtx = 0;
  // TClonesArray *Reco_trk_4mom = 0;
  float trkPt[NTRKMAX];
  float trkEta[NTRKMAX];
  float trkPhi[NTRKMAX];
  //---------------------------
  TBranch *b_Centrality; //!
  TBranch *b_eventNb;    //!
  TBranch *b_runNb;      //!
  TBranch *b_HLTriggers; //!
  TBranch *b_Reco_QQ_size; //!
  TBranch *b_Reco_QQ_trig; //!
  TBranch *b_Reco_QQ_type; //!
  TBranch *b_Reco_QQ_sign; //!
  TBranch *b_Reco_QQ_VtxProb;  //!
  TBranch *b_Reco_QQ_dca; // new float, unused in upsilon
  TBranch *b_Reco_QQ_ctau; //!
  TBranch *b_Reco_QQ_ctauErr; //!
  TBranch *b_Reco_QQ_ctauTrue; //!
  TBranch *b_Reco_QQ_4mom;     //!
  TBranch *b_Reco_QQ_vtx;      //!
  TBranch *b_Reco_QQ_mupl_4mom;//!
  TBranch *b_Reco_QQ_mumi_4mom;//!
  TBranch *b_Reco_QQ_NtrkDeltaR03;//!
  TBranch *b_Reco_QQ_NtrkDeltaR04;//!
  TBranch *b_Reco_QQ_NtrkDeltaR05;//!
  TBranch *b_Reco_QQ_NtrkPt04;//!
  TBranch *b_Reco_QQ_NtrkPt03;//!
  TBranch *b_Reco_QQ_NtrkPt02;//!
  //------------------------------
  TBranch *b_Reco_QQ_mupl_nTrkHits;           
  TBranch *b_Reco_QQ_mupl_normChi2_inner;
  TBranch *b_Reco_QQ_mupl_normChi2_global;
  TBranch *b_Reco_QQ_mupl_dxy;
  TBranch *b_Reco_QQ_mupl_dxyErr;
  TBranch *b_Reco_QQ_mupl_dz;
  TBranch *b_Reco_QQ_mupl_dzErr;
  TBranch *b_Reco_QQ_mupl_TrkMuArb;
  TBranch *b_Reco_QQ_mupl_TMOneStaTight; // new bool, unused in upsilon(?)
  TBranch *b_Reco_QQ_mupl_nMuValHits; // new int, unused in upsilon
  TBranch *b_Reco_QQ_mupl_nPixWMea; // new int, unused in upsilon
  TBranch *b_Reco_QQ_mupl_nTrkWMea; // new int, unused in upsilon
  TBranch *b_Reco_QQ_mupl_StationsMatched; // new int, unused in Upsilon
  TBranch *b_Reco_QQ_mupl_pt_inner; // new float, unused.
  TBranch *b_Reco_QQ_mupl_pt_global; //new float, unused.
  TBranch *b_Reco_QQ_mupl_ptErr_inner; // new float, unused.
  TBranch *b_Reco_QQ_mupl_ptErr_global; //new float, unused.
    
  TBranch *b_Reco_QQ_mumi_TMOneStaTight;  // new bool, unused in upsilon(?)
  TBranch *b_Reco_QQ_mumi_nMuValHits; // new int, unused in upsilon
  TBranch *b_Reco_QQ_mumi_nPixWMea; // new int, unused in upsilon
  TBranch *b_Reco_QQ_mumi_nTrkWMea; // new int, unused in upsilon
  TBranch *b_Reco_QQ_mumi_StationsMatched; // new int, unused in Upsilon
  TBranch *b_Reco_QQ_mumi_pt_inner; // new float, unused.
  TBranch *b_Reco_QQ_mumi_pt_global; //new float, unused.
  TBranch *b_Reco_QQ_mumi_ptErr_inner; // new float, unused.
  TBranch *b_Reco_QQ_mumi_ptErr_global; //new float, unused.
  TBranch *b_Reco_QQ_mumi_normChi2_inner;
  TBranch *b_Reco_QQ_mumi_normChi2_global;
  TBranch *b_Reco_QQ_mumi_dxy;
  TBranch *b_Reco_QQ_mumi_dxyErr;
  TBranch *b_Reco_QQ_mumi_dz;
  TBranch *b_Reco_QQ_mumi_dzErr;
  TBranch *b_Reco_QQ_mumi_TrkMuArb;
  TBranch *b_Reco_QQ_mumi_nTrkHits;
  //-------------  //-------------  //------------- 
  TBranch *b_Reco_mu_size;  //!
  TBranch *b_Reco_mu_type;  //!
  TBranch *b_Reco_mu_charge;//!
  TBranch *b_Reco_mu_4mom;  //!
  //-------------  //-------------  //------------- 
  TBranch *b_zVtx; //!
  TBranch *b_Npix; //!
  TBranch *b_NpixelTracks; //!
  TBranch *b_Ntracks; //!
  TBranch *b_SumET_HF; //!
  TBranch *b_SumET_HFplus; //!
  TBranch *b_SumET_HFminus; //!
  TBranch *b_SumET_HFplusEta4; //!
  TBranch *b_SumET_HFminusEta4; //!
  TBranch *b_SumET_ET; //!
  TBranch *b_SumET_EE; //!
  TBranch *b_SumET_EB; //!
  TBranch *b_SumET_EEplus; //!
  TBranch *b_SumET_EEminus; //!
  TBranch *b_SumET_ZDC; //!
  TBranch *b_SumET_ZDCplus; //!
  TBranch *b_SumET_ZDCminus; //!

  //------------------------------
  // TBranch *b_Reco_trk_size; //! 
  // TBranch *b_Reco_trk_4mom; //!
  // TBranch *b_Reco_trk_vtx;  //!

  /// new tree variables
  Int_t _mupl_nTrkHits;           
  Float_t _mupl_normChi2_inner;
  Float_t _mupl_normChi2_global;
  Float_t _mupl_dxy;
  Float_t _mupl_dxyErr;
  Float_t _mupl_dz;
  Float_t _mupl_dzErr;
  Bool_t _mupl_TrkMuArb;
  Bool_t _mupl_TMOneStaTight; // new bool, unused in upsilon(?)
  Int_t _mupl_nMuValHits; // new int, unused in upsilon
  Int_t _mupl_nPixWMea; // new int, unused in upsilon
  Int_t _mupl_nTrkWMea; // new int, unused in upsilon
  Int_t _mupl_StationsMatched; // new int, unused in Upsilon
  Float_t _mupl_pt_inner; // new float, unused.
  Float_t _mupl_pt_global; //new float, unused.
  Float_t _mupl_ptErr_inner; // new float, unused.
  Float_t _mupl_ptErr_global; //new float, unused.
    
  Bool_t _mumi_TMOneStaTight;  // new bool, unused in upsilon(?)
  Int_t _mumi_nMuValHits; // new int, unused in upsilon
  Int_t _mumi_nPixWMea; // new int, unused in upsilon
  Int_t _mumi_nTrkWMea; // new int, unused in upsilon
  Int_t _mumi_StationsMatched; // new int, unused in Upsilon
  Float_t _mumi_pt_inner; // new float, unused.
  Float_t _mumi_pt_global; //new float, unused.
  Float_t _mumi_ptErr_inner; // new float, unused.
  Float_t _mumi_ptErr_global; //new float, unused.
  Float_t _mumi_normChi2_inner;
  Float_t _mumi_normChi2_global;
  Float_t _mumi_dxy;
  Float_t _mumi_dxyErr;
  Float_t _mumi_dz;
  Float_t _mumi_dzErr;
  Bool_t _mumi_TrkMuArb;
  Int_t _mumi_nTrkHits;
  Float_t _dca;
   float _ctau; // new float, unused in upsilon
  float _ctauTrue; // new float, unused in upsilon
  float _ctauErr; // new float, unused in upsilon
  float _zVtx;
  // ###### read input TTree
  t->SetBranchAddress("Centrality",    &Centrality,   &b_Centrality);
  t->SetBranchAddress("eventNb",       &eventNb,      &b_eventNb);
  t->SetBranchAddress("runNb",         &runNb,        &b_runNb);
  t->SetBranchAddress("HLTriggers",    &HLTriggers,   &b_HLTriggers);
  t->SetBranchAddress("Reco_QQ_size",  &Reco_QQ_size, &b_Reco_QQ_size);
  t->SetBranchAddress("Reco_QQ_type",  &Reco_QQ_type, &b_Reco_QQ_type);
  t->SetBranchAddress("Reco_QQ_sign",  &Reco_QQ_sign, &b_Reco_QQ_sign);
  t->SetBranchAddress("Reco_QQ_4mom",  &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  t->SetBranchAddress("Reco_QQ_vtx",       &Reco_QQ_vtx,        &b_Reco_QQ_vtx);
  t->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  t->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  
  t->SetBranchAddress("Reco_mu_size",  &Reco_mu_size,   &b_Reco_mu_size);
  t->SetBranchAddress("Reco_mu_type",  Reco_mu_type,   &b_Reco_mu_type);
  t->SetBranchAddress("Reco_mu_charge",&Reco_mu_charge, &b_Reco_mu_charge);
  t->SetBranchAddress("Reco_mu_4mom",  &Reco_mu_4mom,   &b_Reco_mu_4mom);

  t->SetBranchAddress("Reco_QQ_NtrkDeltaR03",  &Reco_QQ_NtrkDeltaR03, &b_Reco_QQ_NtrkDeltaR03);
  t->SetBranchAddress("Reco_QQ_NtrkDeltaR04",  &Reco_QQ_NtrkDeltaR04, &b_Reco_QQ_NtrkDeltaR04);
  t->SetBranchAddress("Reco_QQ_NtrkDeltaR05",  &Reco_QQ_NtrkDeltaR05, &b_Reco_QQ_NtrkDeltaR05);

  t->SetBranchAddress("Reco_QQ_NtrkPt04",  &Reco_QQ_NtrkPt04, &b_Reco_QQ_NtrkPt04);
  t->SetBranchAddress("Reco_QQ_NtrkPt03",  &Reco_QQ_NtrkPt03, &b_Reco_QQ_NtrkPt03);
  t->SetBranchAddress("Reco_QQ_NtrkPt02",  &Reco_QQ_NtrkPt02, &b_Reco_QQ_NtrkPt02);
  // additional selection
  // id variables
  t->SetBranchAddress("Reco_QQ_mupl_nTrkHits", &Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
  t->SetBranchAddress("Reco_QQ_mupl_normChi2_inner", &Reco_QQ_mupl_normChi2_inner, &b_Reco_QQ_mupl_normChi2_inner);
  t->SetBranchAddress("Reco_QQ_mupl_normChi2_global", &Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
  t->SetBranchAddress("Reco_QQ_mupl_dxy", &Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  t->SetBranchAddress("Reco_QQ_mupl_dz", &Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  t->SetBranchAddress("Reco_QQ_mupl_dxyErr", &Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
  t->SetBranchAddress("Reco_QQ_mupl_dzErr", &Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  t->SetBranchAddress("Reco_QQ_mupl_TrkMuArb", &Reco_QQ_mupl_TrkMuArb, &b_Reco_QQ_mupl_TrkMuArb);
  t->SetBranchAddress("Reco_QQ_mupl_StationsMatched",&Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched); // new int, unused in Upsilon
  t->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", &Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);  // new bool, unused in upsilon(?)
  t->SetBranchAddress("Reco_QQ_mupl_nMuValHits", &Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);  // new int, unused in upsilon
  t->SetBranchAddress("Reco_QQ_mupl_nPixWMea",&Reco_QQ_mupl_nPixWMea,&b_Reco_QQ_mupl_nPixWMea); // new int, unused in upsilon
  t->SetBranchAddress("Reco_QQ_mupl_nTrkWMea",&Reco_QQ_mupl_nTrkWMea,&b_Reco_QQ_mupl_nTrkWMea); // new int, unused in upsilon
  t->SetBranchAddress("Reco_QQ_mupl_pt_inner",&Reco_QQ_mupl_pt_inner,&b_Reco_QQ_mupl_pt_inner); // new float, unused.
  t->SetBranchAddress("Reco_QQ_mupl_pt_global",&Reco_QQ_mupl_pt_global,&b_Reco_QQ_mupl_pt_global);//new float, unused.
  t->SetBranchAddress("Reco_QQ_mupl_ptErr_inner",&Reco_QQ_mupl_ptErr_inner,&b_Reco_QQ_mupl_ptErr_inner); // new float, unused.
  t->SetBranchAddress("Reco_QQ_mupl_ptErr_global",&Reco_QQ_mupl_ptErr_global,&b_Reco_QQ_mupl_ptErr_global);//new float, unused.
  //muminus
  t->SetBranchAddress("Reco_QQ_mumi_nTrkHits", &Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
  t->SetBranchAddress("Reco_QQ_mumi_normChi2_inner", &Reco_QQ_mumi_normChi2_inner, &b_Reco_QQ_mumi_normChi2_inner);
  t->SetBranchAddress("Reco_QQ_mumi_normChi2_global", &Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
  t->SetBranchAddress("Reco_QQ_mumi_dxy", &Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  t->SetBranchAddress("Reco_QQ_mumi_dz", &Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  t->SetBranchAddress("Reco_QQ_mumi_dxyErr", &Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
  t->SetBranchAddress("Reco_QQ_mumi_dzErr", &Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  t->SetBranchAddress("Reco_QQ_mumi_TrkMuArb", &Reco_QQ_mumi_TrkMuArb, &b_Reco_QQ_mumi_TrkMuArb);  
  t->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", &Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);  // new bool, unused in upsilon(?)
  t->SetBranchAddress("Reco_QQ_mumi_nMuValHits", &Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);  // new int, unused in upsilon
  t->SetBranchAddress("Reco_QQ_mumi_nPixWMea",&Reco_QQ_mumi_nPixWMea,&b_Reco_QQ_mumi_nPixWMea); // new int, unused in upsilon
  t->SetBranchAddress("Reco_QQ_mumi_nTrkWMea",&Reco_QQ_mumi_nTrkWMea,&b_Reco_QQ_mumi_nTrkWMea); // new int, unused in upsilon
  t->SetBranchAddress("Reco_QQ_mumi_StationsMatched",&Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched); // new int, unused in Upsilon
  t->SetBranchAddress("Reco_QQ_mumi_pt_inner",&Reco_QQ_mumi_pt_inner,&b_Reco_QQ_mumi_pt_inner); // new float, unused.
  t->SetBranchAddress("Reco_QQ_mumi_pt_global",&Reco_QQ_mumi_pt_global,&b_Reco_QQ_mumi_pt_global);//new float, unused.
  t->SetBranchAddress("Reco_QQ_mumi_ptErr_inner",&Reco_QQ_mumi_ptErr_inner,&b_Reco_QQ_mumi_ptErr_inner); // new float, unused.
  t->SetBranchAddress("Reco_QQ_mumi_ptErr_global",&Reco_QQ_mumi_ptErr_global,&b_Reco_QQ_mumi_ptErr_global);//new float, unused.
  //dimuons
  t->SetBranchAddress("Reco_QQ_dca",    &Reco_QQ_dca,   &b_Reco_QQ_dca); 
  t->SetBranchAddress("Reco_QQ_trig",  &Reco_QQ_trig, &b_Reco_QQ_trig);
  t->SetBranchAddress("Reco_QQ_VtxProb",   &Reco_QQ_VtxProb,    &b_Reco_QQ_VtxProb);
  t->SetBranchAddress("Reco_QQ_ctau",  &Reco_QQ_ctau, &b_Reco_QQ_ctau);
  t->SetBranchAddress("Reco_QQ_ctauErr",  &Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
  t->SetBranchAddress("Reco_QQ_ctauTrue",  &Reco_QQ_ctauTrue, &b_Reco_QQ_ctauTrue);
  t->SetBranchAddress("zVtx",             &zVtx,              &b_zVtx);
  // extra centrality variables
  t->SetBranchAddress("Npix",             &Npix,              &b_Npix);
  t->SetBranchAddress("NpixelTracks",     &NpixelTracks,      &b_NpixelTracks);
  t->SetBranchAddress("Ntracks",          &Ntracks,           &b_Ntracks);
  t->SetBranchAddress("SumET_HF",         &SumET_HF,          &b_SumET_HF);
  t->SetBranchAddress("SumET_HFplus",     &SumET_HFplus,      &b_SumET_HFplus);
  t->SetBranchAddress("SumET_HFminus",    &SumET_HFminus,     &b_SumET_HFminus);
  t->SetBranchAddress("SumET_HFplusEta4", &SumET_HFplusEta4,  &b_SumET_HFplusEta4);
  t->SetBranchAddress("SumET_HFminusEta4",&SumET_HFminusEta4, &b_SumET_HFminusEta4);
  t->SetBranchAddress("SumET_ET",         &SumET_ET,          &b_SumET_ET);
  t->SetBranchAddress("SumET_EE",         &SumET_EE,          &b_SumET_EE);
  t->SetBranchAddress("SumET_EB",         &SumET_EB,          &b_SumET_EB);
  t->SetBranchAddress("SumET_EEplus",     &SumET_EEplus,      &b_SumET_EEplus);
  t->SetBranchAddress("SumET_EEminus",    &SumET_EEminus,     &b_SumET_EEminus);
  t->SetBranchAddress("SumET_ZDC"    ,    &SumET_ZDC,         &b_SumET_ZDC);
  t->SetBranchAddress("SumET_ZDCplus",    &SumET_ZDCplus,     &b_SumET_ZDCplus);
  t->SetBranchAddress("SumET_ZDCminus",   &SumET_ZDCminus,    &b_SumET_ZDCminus);
  
  // extra track variable  
  // t->SetBranchAddress("Reco_trk_size",  &Reco_trk_size,   &b_Reco_trk_size);
  // t->SetBranchAddress("Reco_trk_4mom",  &Reco_trk_4mom,   &b_Reco_trk_4mom);
  // t->SetBranchAddress("Reco_trk_vtx",   &Reco_trk_vtx,    &b_Reco_trk_vtx);

  // #### define control histograms
  TH1F *h_QQ_mass   = new TH1F("h_QQ_mass","",nBins, mass_min,mass_max);  // all OS
  TH1F *h_QQ_mass_1 = new TH1F("h_QQ_mass_1","",nBins, mass_min,mass_max);// OS in acceptance
  TH1F *h_QQ_mass_2 = new TH1F("h_QQ_mass_2","",nBins, mass_min, mass_max);// SS 
  
  // ##### output file
  TFile *f1 = new TFile(Form("../dimuonTree_%s_Run%s_trigBit%d_allTriggers%d.root",dataSource,runNumber,nTriggerBit,bAllTriggers),"RECREATE");

  TTree *MuTree              = new TTree("MuTree","MuTree");
  TTree *DiMuonTree         = new TTree("DiMuonTree","DiMuonTree");
  TTree *DiMuonTree_allsign = new TTree("DiMuonTree_allsign","DiMuonTree_allsign");
  TTree *DiMuonTree_trkRot  = new TTree("DiMuonTree_trkRot","DiMuonTree_trkRot");

  MuTree->Branch("Reco_mu_size",   &Reco_mu_size,   "Reco_mu_size/I");
  MuTree->Branch("nPlusMu",        &nPlusMu,        "nPlusMu/I");
  MuTree->Branch("nMinusMu",       &nMinusMu,       "nMinusMu/I");
  MuTree->Branch("muPt",           &muPt,           "muPt[Reco_mu_size]/F");
  MuTree->Branch("muEta",          &muEta,          "muEta[Reco_mu_size]/F");
  MuTree->Branch("muPhi",          &muPhi,          "muPhi[Reco_mu_size]/F");
  MuTree->Branch("Reco_mu_type",   &Reco_mu_type,    "Reco_mu_type[Reco_mu_size]/I");
  MuTree->Branch("Reco_mu_charge", &Reco_mu_charge,  "Reco_mu_charge[Reco_mu_size]/I" );
  MuTree->Branch("eventNb",        &eventNb,        "eventNb/I");
  MuTree->Branch("runNb",          &runNb,          "runNb/I");

  //DiMuonTree->Branch("nReco_QQ", &nReco_QQ, "nReco_QQ/I");
  DiMuonTree->Branch("Centrality",    &Centrality,    "Centrality/I");
  DiMuonTree->Branch("HLTriggers",    &HLTriggers,    "HLTriggers/I");
  DiMuonTree->Branch("QQsign",        &QQsign,        "QQsign/I");
  DiMuonTree->Branch("QQTrkDr03",     &QQTrkDr03,     "QQTrkDr03/I");
  DiMuonTree->Branch("QQTrkDr04",     &QQTrkDr04,     "QQTrkDr04/I");
  DiMuonTree->Branch("QQTrkDr05",     &QQTrkDr05,     "QQTrkDr05/I");

  DiMuonTree->Branch("QQTrkPt02",     &QQTrkPt02,     "QQTrkPt02/I");
  DiMuonTree->Branch("QQTrkPt03",     &QQTrkPt03,     "QQTrkPt03/I");
  DiMuonTree->Branch("QQTrkPt04",     &QQTrkPt04,     "QQTrkPt04/I");
  
  DiMuonTree->Branch("eventNb",       &eventNb,     "eventNb/I");
  DiMuonTree->Branch("runNb",         &runNb,     "runNb/I");
  DiMuonTree->Branch("invariantMass", &invariantMass, "invariantMass/F");
  DiMuonTree->Branch("dimuPt",         &dimuPt, "dimuPt/F");
  DiMuonTree->Branch("dimuEta",        &dimuEta, "dimuEta/F");
  DiMuonTree->Branch("dimuPhi",        &dimuPhi, "dimuPhi/F");
  DiMuonTree->Branch("dimuRapidity",   &dimuRapidity, "dimuRapidity/F");
  DiMuonTree->Branch("muPlusPt",      &muPlusPt, "muPlusPt/F");
  DiMuonTree->Branch("muMinusPt",     &muMinusPt, "muMinusPt/F");
  DiMuonTree->Branch("muPlusEta",     &muPlusEta, "muPlusEta/F");
  DiMuonTree->Branch("muMinusEta",    &muMinusEta, "muMinusEta/F");
  DiMuonTree->Branch("muPlusPhi",     &muPlusPhi, "muPlusPhi/F");
  DiMuonTree->Branch("muMinusPhi",    &muMinusPhi, "muMinusPhi/F");
  // additional selection
  // id muplus
  DiMuonTree->Branch("_mupl_TrkMuArb",&_mupl_TrkMuArb,"_mupl_TrkMuArb/O");
  DiMuonTree->Branch("_mupl_TMOneStaTight",&_mupl_TMOneStaTight,"_mupl_TMOneStaTight/O");  // new bool, unused in upsilon(?)
  DiMuonTree->Branch("_mupl_nMuValHits",&_mupl_nMuValHits,"_mupl_nMuValHits/I"); // new int, unused in upsilon
  DiMuonTree->Branch("_mupl_nPixWMea",&_mupl_nPixWMea,"_mupl_nPixWMea/I"); // new int, unused in upsilon
  DiMuonTree->Branch("_mupl_nTrkWMea",&_mupl_nTrkWMea,"_mupl_nTrkWMea/I"); // new int, unused in upsilon
  DiMuonTree->Branch("_mupl_StationsMatched",&_mupl_StationsMatched,"_mupl_StationsMatched/I"); // new int, unused in Upsilon
  DiMuonTree->Branch("_mupl_pt_inner",&_mupl_pt_inner,"_mupl_pt_inner/F"); // new float, unused.
  DiMuonTree->Branch("_mupl_pt_global",&_mupl_pt_global,"_mupl_pt_global/F"); //new float, unused.
  DiMuonTree->Branch("_mupl_ptErr_inner",&_mupl_ptErr_inner,"_mupl_ptErr_inner/F"); // new float, unused.
  DiMuonTree->Branch("_mupl_ptErr_global",&_mupl_ptErr_global,"_mupl_ptErr_global/F"); //new float, unused.
  DiMuonTree->Branch("_mupl_nTrkHits",&_mupl_nTrkHits," _mupl_nTrkHits/I");
  DiMuonTree->Branch("_mupl_normChi2_inner", &_mupl_normChi2_inner," _mupl_normChi2_inner/F");
  DiMuonTree->Branch("_mupl_normChi2_global",&_mupl_normChi2_global,"_mupl_normChi2_global/F");
  DiMuonTree->Branch("_mupl_dxy",&_mupl_dxy,"_mupl_dxy/F");
  DiMuonTree->Branch("_mupl_dxyErr",&_mupl_dxyErr,"_mupl_dxyErr/F");
  DiMuonTree->Branch("_mupl_dz",&_mupl_dz,"_mupl_dz/F");
  DiMuonTree->Branch("_mupl_dzErr",&_mupl_dzErr,"_mupl_dzErr/F");
  //id muminus
  DiMuonTree->Branch("_mumi_TrkMuArb",&_mumi_TrkMuArb,"_mumi_TrkMuArb/O");
  DiMuonTree->Branch("_mumi_TMOneStaTight",&_mumi_TMOneStaTight,"_mumi_TMOneStaTight/O");  // new bool, unused in upsilon(?)
  DiMuonTree->Branch("_mumi_nMuValHits",&_mumi_nMuValHits,"_mumi_nMuValHits/I"); // new int, unused in upsilon
  DiMuonTree->Branch("_mumi_nPixWMea",&_mumi_nPixWMea,"_mumi_nPixWMea/I"); // new int, unused in upsilon
  DiMuonTree->Branch("_mumi_nTrkWMea",&_mumi_nTrkWMea,"_mumi_nTrkWMea/I"); // new int, unused in upsilon
  DiMuonTree->Branch("_mumi_StationsMatched",&_mumi_StationsMatched,"_mumi_StationsMatched/I"); // new int, unused in Upsilon
  DiMuonTree->Branch("_mumi_pt_inner",&_mumi_pt_inner,"_mumi_pt_inner/F"); // new float, unused.
  DiMuonTree->Branch("_mumi_pt_global",&_mumi_pt_global,"_mumi_pt_global/F"); //new float, unused.
  DiMuonTree->Branch("_mumi_ptErr_inner",&_mumi_ptErr_inner,"_mumi_ptErr_inner/F"); // new float, unused.
  DiMuonTree->Branch("_mumi_ptErr_global",&_mumi_ptErr_global,"_mumi_ptErr_global/F"); //new float, unused.
  DiMuonTree->Branch("_mumi_nTrkHits",&_mumi_nTrkHits," _mumi_nTrkHits/I");
  DiMuonTree->Branch("_mumi_normChi2_inner", &_mumi_normChi2_inner," _mumi_normChi2_inner/F");
  DiMuonTree->Branch("_mumi_normChi2_global",&_mumi_normChi2_global,"_mumi_normChi2_global/F");
  DiMuonTree->Branch("_mumi_dxy",&_mumi_dxy,"_mumi_dxy/F");
  DiMuonTree->Branch("_mumi_dxyErr",&_mumi_dxyErr,"_mumi_dxyErr/F");
  DiMuonTree->Branch("_mumi_dz",&_mumi_dz,"_mumi_dz/F");
  DiMuonTree->Branch("_mumi_dzErr",&_mumi_dzErr,"_mumi_dzErr/F");
 //dimuon variables
  DiMuonTree->Branch("QQtrig",        &QQtrig,        "QQtrig/I");
  DiMuonTree->Branch("_dca",&_dca,"_dca/F");// new float, unused in upsilon
  DiMuonTree->Branch("_ctau",&_ctau,"_ctau/F");// new float, unused in upsilon
  DiMuonTree->Branch("_ctauErr",&_ctauErr,"_ctauErr/F");// new float, unused in upsilon
  DiMuonTree->Branch("_ctauTrue",&_ctauTrue,"_ctauTrue/F");// new float, unused in upsilon
  DiMuonTree->Branch("_zVtx",  &_zVtx,"_zVtx/F");
  DiMuonTree->Branch("vProb",         &vProb,     "vProb/F");
  if(addExtraCentrality)
    {
      DiMuonTree->Branch("Npix",&Npix,"Npix/I");
      DiMuonTree->Branch("NpixelTracks",&NpixelTracks,"NpixelTracks/I");
      DiMuonTree->Branch("Ntracks", &Ntracks, "Ntracks/I");
      DiMuonTree->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
      DiMuonTree->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
      DiMuonTree->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
      DiMuonTree->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
      DiMuonTree->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
      DiMuonTree->Branch("SumET_ET",&SumET_ET,"SumET_ET/F");
      DiMuonTree->Branch("SumET_EE",&SumET_EE,"SumET_EE/F");
      DiMuonTree->Branch("SumET_EB",&SumET_EB,"SumET_EB/F");
      DiMuonTree->Branch("SumET_EEplus",&SumET_EEplus,"SumET_EEplus/F");
      DiMuonTree->Branch("SumET_EEminus",&SumET_EEminus,"SumET_EEminus/F");
      DiMuonTree->Branch("SumET_ZDC",&SumET_ZDC,"SumET_ZDC/F");
      DiMuonTree->Branch("SumET_ZDCplus",&SumET_ZDCplus,"SumET_ZDCplus/F");
      DiMuonTree->Branch("SumET_ZDCminus",&SumET_ZDCminus,"SumET_ZDCminus/F");
    }

  DiMuonTree_allsign->Branch("Centrality", &Centrality,    "Centrality/I");
  DiMuonTree_allsign->Branch("HLTriggers", &HLTriggers,    "HLTriggers/I");
  DiMuonTree_allsign->Branch("QQtrig",     &QQtrig,    "QQtrig/I");
  DiMuonTree_allsign->Branch("QQsign",     &QQsign,    "QQsign/I");
  DiMuonTree_allsign->Branch("QQTrkDr03",     &QQTrkDr03,     "QQTrkDr03/I");
  DiMuonTree_allsign->Branch("QQTrkDr04",     &QQTrkDr04,     "QQTrkDr04/I");
  DiMuonTree_allsign->Branch("QQTrkDr05",     &QQTrkDr05,     "QQTrkDr05/I");
  DiMuonTree_allsign->Branch("QQTrkPt04",     &QQTrkPt04,     "QQTrkPt04/I");
  DiMuonTree_allsign->Branch("QQTrkPt03",     &QQTrkPt03,     "QQTrkPt03/I");
  DiMuonTree_allsign->Branch("QQTrkPt02",     &QQTrkPt02,     "QQTrkPt02/I");
  DiMuonTree_allsign->Branch("weight",     &weight,    "weight/F");
  DiMuonTree_allsign->Branch("weight2",    &weight2,    "weight2/F");
  DiMuonTree_allsign->Branch("vProb",      &vProb,     "vProb/F");
  DiMuonTree_allsign->Branch("eventNb",    &eventNb,     "eventNb/I");
  DiMuonTree_allsign->Branch("runNb",      &runNb,     "runNb/I");
  DiMuonTree_allsign->Branch("invariantMass", &invariantMass, "invariantMass/F");
  DiMuonTree_allsign->Branch("dimuPt",      &dimuPt, "dimuPt/F");
  DiMuonTree_allsign->Branch("dimuEta",     &dimuEta, "dimuEta/F");
  DiMuonTree_allsign->Branch("dimuPhi",     &dimuPhi, "dimuPhi/F");
  DiMuonTree_allsign->Branch("dimuRapidity",&dimuRapidity, "dimuRapidity/F");
  DiMuonTree_allsign->Branch("muPlusPt",   &muPlusPt, "muPlusPt/F");
  DiMuonTree_allsign->Branch("muMinusPt",  &muMinusPt, "muMinusPt/F");
  DiMuonTree_allsign->Branch("muPlusEta",  &muPlusEta, "muPlusEta/F");
  DiMuonTree_allsign->Branch("muMinusEta", &muMinusEta, "muMinusEta/F");
  DiMuonTree_allsign->Branch("muPlusPhi",  &muPlusPhi, "muPlusPhi/F");
  DiMuonTree_allsign->Branch("muMinusPhi", &muMinusPhi, "muMinusPhi/F");
  
   if(addExtraCentrality)
    {
      DiMuonTree_allsign->Branch("Npix",&Npix,"Npix/I");
      DiMuonTree_allsign->Branch("NpixelTracks",&NpixelTracks,"NpixelTracks/I");
      DiMuonTree_allsign->Branch("Ntracks", &Ntracks, "Ntracks/I");
      DiMuonTree_allsign->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
      DiMuonTree_allsign->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
      DiMuonTree_allsign->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
      DiMuonTree_allsign->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
      DiMuonTree_allsign->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
      DiMuonTree_allsign->Branch("SumET_ET",&SumET_ET,"SumET_ET/F");
      DiMuonTree_allsign->Branch("SumET_EE",&SumET_EE,"SumET_EE/F");
      DiMuonTree_allsign->Branch("SumET_EB",&SumET_EB,"SumET_EB/F");
      DiMuonTree_allsign->Branch("SumET_EEplus",&SumET_EEplus,"SumET_EEplus/F");
      DiMuonTree_allsign->Branch("SumET_EEminus",&SumET_EEminus,"SumET_EEminus/F");
      DiMuonTree_allsign->Branch("SumET_ZDC",&SumET_ZDC,"SumET_ZDC/F");
      DiMuonTree_allsign->Branch("SumET_ZDCplus",&SumET_ZDCplus,"SumET_ZDCplus/F");
      DiMuonTree_allsign->Branch("SumET_ZDCminus",&SumET_ZDCminus,"SumET_ZDCminus/F");
    }

  DiMuonTree_trkRot->Branch("Centrality", &Centrality,    "Centrality/I");
  DiMuonTree_trkRot->Branch("HLTriggers", &HLTriggers,    "HLTriggers/I");
  DiMuonTree_trkRot->Branch("QQtrig",     &QQtrig,    "QQtrig/I");
  DiMuonTree_trkRot->Branch("QQsign",     &QQsign,    "QQsign/I");
  DiMuonTree_trkRot->Branch("QQTrkDr03",     &QQTrkDr03,     "QQTrkDr03/I");
  DiMuonTree_trkRot->Branch("QQTrkDr04",     &QQTrkDr04,     "QQTrkDr04/I");
  DiMuonTree_trkRot->Branch("QQTrkDr05",     &QQTrkDr05,     "QQTrkDr05/I");
  DiMuonTree_trkRot->Branch("QQTrkPt04",     &QQTrkPt04,     "QQTrkPt04/I");
  DiMuonTree_trkRot->Branch("QQTrkPt03",     &QQTrkPt03,     "QQTrkPt03/I");
  DiMuonTree_trkRot->Branch("QQTrkPt02",     &QQTrkPt02,     "QQTrkPt02/I");

  DiMuonTree_trkRot->Branch("weight",     &weight,    "weight/F");
  DiMuonTree_trkRot->Branch("weight2",    &weight2,    "weight2/F");
  DiMuonTree_trkRot->Branch("vProb",      &vProb,     "vProb/F");
  DiMuonTree_trkRot->Branch("eventNb",    &eventNb,     "eventNb/I");
  DiMuonTree_trkRot->Branch("runNb",      &runNb,     "runNb/I");
  DiMuonTree_trkRot->Branch("invariantMass", &invariantMass, "invariantMass/F");
  DiMuonTree_trkRot->Branch("dimuPt",      &dimuPt, "dimuPt/F");
  DiMuonTree_trkRot->Branch("dimuEta",     &dimuEta, "dimuEta/F");
  DiMuonTree_trkRot->Branch("dimuPhi",     &dimuPhi, "dimuPhi/F");
  DiMuonTree_trkRot->Branch("dimuRapidity",&dimuRapidity, "dimuRapidity/F");
  DiMuonTree_trkRot->Branch("muPlusPt",   &muPlusPt, "muPlusPt/F");
  DiMuonTree_trkRot->Branch("muMinusPt",  &muMinusPt, "muMinusPt/F");
  DiMuonTree_trkRot->Branch("muPlusEta",  &muPlusEta, "muPlusEta/F");
  DiMuonTree_trkRot->Branch("muMinusEta", &muMinusEta, "muMinusEta/F");
  DiMuonTree_trkRot->Branch("muPlusPhi",  &muPlusPhi, "muPlusPhi/F");
  DiMuonTree_trkRot->Branch("muMinusPhi", &muMinusPhi, "muMinusPhi/F");
  if(addExtraCentrality)
    {
      DiMuonTree_trkRot->Branch("Npix",&Npix,"Npix/I");
      DiMuonTree_trkRot->Branch("NpixelTracks",&NpixelTracks,"NpixelTracks/I");
      DiMuonTree_trkRot->Branch("Ntracks", &Ntracks, "Ntracks/I");
      DiMuonTree_trkRot->Branch("SumET_HF",&SumET_HF,"SumET_HF/F");
      DiMuonTree_trkRot->Branch("SumET_HFplus",&SumET_HFplus,"SumET_HFplus/F");
      DiMuonTree_trkRot->Branch("SumET_HFminus",&SumET_HFminus,"SumET_HFminus/F");
      DiMuonTree_trkRot->Branch("SumET_HFplusEta4",&SumET_HFplusEta4,"SumET_HFplusEta4/F");
      DiMuonTree_trkRot->Branch("SumET_HFminusEta4",&SumET_HFminusEta4,"SumET_HFminusEta4/F");
      DiMuonTree_trkRot->Branch("SumET_ET",&SumET_ET,"SumET_ET/F");
      DiMuonTree_trkRot->Branch("SumET_EE",&SumET_EE,"SumET_EE/F");
      DiMuonTree_trkRot->Branch("SumET_EB",&SumET_EB,"SumET_EB/F");
      DiMuonTree_trkRot->Branch("SumET_EEplus",&SumET_EEplus,"SumET_EEplus/F");
      DiMuonTree_trkRot->Branch("SumET_EEminus",&SumET_EEminus,"SumET_EEminus/F");
      DiMuonTree_trkRot->Branch("SumET_ZDC",&SumET_ZDC,"SumET_ZDC/F");
      DiMuonTree_trkRot->Branch("SumET_ZDCplus",&SumET_ZDCplus,"SumET_ZDCplus/F");
      DiMuonTree_trkRot->Branch("SumET_ZDCminus",&SumET_ZDCminus,"SumET_ZDCminus/F");
    }


  //____________________________________ loop over all events
  for (int i=0; i<nentries; i++) 
    {
      t->GetEntry(i);
      //      if(excludeWrongAlign_pa && runNb<=210658) continue;
     

      // ##### single muon tree
      // countng + and - single muons
      nPlusMu=0;
      nMinusMu=0;
      for(int iMu = 0; iMu < Reco_mu_size; iMu++)
   {
     if (Reco_mu_charge[iMu] == 1) nPlusMu++;
     else nMinusMu++;
     TLorentzVector *Reco_mu = (TLorentzVector *) Reco_mu_4mom->At(iMu);
     muPt[iMu]=Reco_mu->Pt();
     muEta[iMu]=Reco_mu->Eta();
     muPhi[iMu]=Reco_mu->Phi();
   }
      MuTree->Fill();

      // // ntrk loop
      //  for(int iTrk = 0; iTrk < Reco_trk_size; iTrk++)
      // 	{
      // 	  TLorentzVector *Reco_trk = (TLorentzVector *) Reco_trk_4mom->At(iTrk);
      // 	  trkPt[iTrk]=Reco_trk->Pt();
      // 	  trkEta[iTrk]=Reco_trk->Eta();
      // 	  trkPhi[iTrk]=Reco_trk->Phi();
      // 	}

      //----------------------------------------------------
       for(int iQQ = 0; iQQ < Reco_QQ_size; iQQ++)
   {
     vProb  = Reco_QQ_VtxProb[iQQ];
     QQsign = Reco_QQ_sign[iQQ];
     
     QQTrkDr03 = Reco_QQ_NtrkDeltaR03[iQQ];
     QQTrkDr04 = Reco_QQ_NtrkDeltaR04[iQQ];
     QQTrkDr05 = Reco_QQ_NtrkDeltaR05[iQQ];

     QQTrkPt04 = Reco_QQ_NtrkPt04[iQQ];
     QQTrkPt03 = Reco_QQ_NtrkPt03[iQQ];
     QQTrkPt02 = Reco_QQ_NtrkPt02[iQQ];

     //	  have to red a bit about the weighting Zhen calculates ...
     if (Reco_QQ_sign[iQQ] == 0) // opposite sign
       { 
         weight  = 1;
         weight2 = 1;
       }
     else // same sign
       {
         weight  = -1;
         float likesign_comb        = (float)nPlusMu*(nPlusMu-1.0)/2.0+(float)nMinusMu*(nMinusMu-1.0)/2.0; // number of combinatorial pairs C(nplus)^2, C(nminus)^2
       float unlikesign_bkgd_comb = (float)nPlusMu*nMinusMu - (nPlusMu>nMinusMu?nMinusMu:nPlusMu);
       weight2 = -1.0 * unlikesign_bkgd_comb/likesign_comb;
       }
     
     // dimuon variables
     TLorentzVector *Reco_QQ      = (TLorentzVector *) Reco_QQ_4mom->At(iQQ);
     TLorentzVector *Reco_QQ_mupl = (TLorentzVector *) Reco_QQ_mupl_4mom->At(iQQ);
     TLorentzVector *Reco_QQ_mumi = (TLorentzVector *) Reco_QQ_mumi_4mom->At(iQQ);
     invariantMass = Reco_QQ->M();
     dimuPt         = Reco_QQ->Pt();
     dimuEta        = Reco_QQ->Eta();
     dimuPhi        = Reco_QQ->Phi();
     dimuRapidity   = Reco_QQ->Rapidity();
     
     // single muon variables
     muMinusPt     = Reco_QQ_mumi->Pt();
     muMinusEta    = Reco_QQ_mumi->Eta();
     muMinusPhi    = Reco_QQ_mumi->Phi();
     muPlusPt      = Reco_QQ_mupl->Pt();
     muPlusEta     = Reco_QQ_mupl->Eta();
     muPlusPhi     = Reco_QQ_mupl->Phi();
          
     // apply extra selection
     // id muplus

     _mupl_nTrkHits       = Reco_QQ_mupl_nTrkHits[iQQ];
     _mupl_normChi2_inner = Reco_QQ_mupl_normChi2_inner[iQQ];
     _mupl_normChi2_global= Reco_QQ_mupl_normChi2_global[iQQ];
     _mupl_dxy            = Reco_QQ_mupl_dxy[iQQ];
     _mupl_dxyErr         = Reco_QQ_mupl_dxyErr[iQQ];
     _mupl_dz             = Reco_QQ_mupl_dz[iQQ];
     _mupl_dzErr          = Reco_QQ_mupl_dzErr[iQQ];
     _mupl_TrkMuArb       = Reco_QQ_mupl_TrkMuArb[iQQ];
     _mupl_TMOneStaTight  = Reco_QQ_mupl_TMOneStaTight[iQQ];
     _mupl_nMuValHits     = Reco_QQ_mupl_nMuValHits[iQQ];
     _mupl_nPixWMea       = Reco_QQ_mupl_nPixWMea[iQQ];
     _mupl_nTrkWMea       = Reco_QQ_mupl_nTrkWMea[iQQ];
     _mupl_StationsMatched = Reco_QQ_mupl_StationsMatched[iQQ]; 
     _mupl_pt_inner      = Reco_QQ_mupl_pt_inner[iQQ];
     _mupl_pt_global     = Reco_QQ_mupl_pt_global[iQQ];
     _mupl_ptErr_inner   = Reco_QQ_mupl_ptErr_inner[iQQ];
     _mupl_ptErr_global  = Reco_QQ_mupl_ptErr_global[iQQ];

     // id muminus
     
     _mumi_nTrkHits       = Reco_QQ_mumi_nTrkHits[iQQ];
     _mumi_normChi2_inner = Reco_QQ_mumi_normChi2_inner[iQQ];
     _mumi_normChi2_global= Reco_QQ_mumi_normChi2_global[iQQ];
     _mumi_dxy            = Reco_QQ_mumi_dxy[iQQ];
     _mumi_dxyErr         = Reco_QQ_mumi_dxyErr[iQQ];
     _mumi_dz             = Reco_QQ_mumi_dz[iQQ];
     _mumi_dzErr          = Reco_QQ_mumi_dzErr[iQQ];
     _mumi_TrkMuArb       = Reco_QQ_mumi_TrkMuArb[iQQ];	
     _mumi_TMOneStaTight  =Reco_QQ_mumi_TMOneStaTight[iQQ];
     _mumi_nMuValHits     = Reco_QQ_mumi_nMuValHits[iQQ];
     _mumi_nPixWMea       = Reco_QQ_mumi_nPixWMea[iQQ];
     _mumi_nTrkWMea       = Reco_QQ_mumi_nTrkWMea[iQQ];
     _mumi_StationsMatched = Reco_QQ_mumi_StationsMatched[iQQ];
     _mumi_pt_inner       = Reco_QQ_mumi_pt_inner[iQQ]; 
     _mumi_pt_global      = Reco_QQ_mumi_pt_global[iQQ];
     _mumi_ptErr_inner    = Reco_QQ_mumi_ptErr_inner[iQQ];
     _mumi_ptErr_global   = Reco_QQ_mumi_ptErr_global[iQQ];

     //dimuon variables
     
     QQtrig = Reco_QQ_trig[iQQ];
     _ctau = (invariantMass/3.0968)*Reco_QQ_ctau[iQQ];
       _ctauTrue = (invariantMass/3.0968)*Reco_QQ_ctauTrue[iQQ];
       _ctauErr = (invariantMass/3.0968)*Reco_QQ_ctauErr[iQQ];
     _zVtx=zVtx[iQQ];	     
     _dca=Reco_QQ_dca[iQQ];
    
     bool bProcess = false;
     if(!bAllTriggers) bProcess = ((HLTriggers&nTriggerBit)==nTriggerBit && (Reco_QQ_trig[iQQ]&nTriggerBit)==nTriggerBit && 
               vProb>newVtxProbCut && muMinusPt>ptcut && muPlusPt>ptcut);
     // for PbPB sample, the v1 and v2 same trigger is in nNtriggerBit=1 and nTriggerBit=2 respectivelly
     /* if(isAArereco) bProcess = (( ( (HLTriggers&nTriggerBit)==nTriggerBit && (Reco_QQ_trig[iQQ]&nTriggerBit)==nTriggerBit ) || 
                     ( (HLTriggers&(nTriggerBit+1))==(nTriggerBit+1) && (Reco_QQ_trig[iQQ]&(nTriggerBit+1))==(nTriggerBit+1) ) )
                   &&  vProb>newVtxProbCut && muMinusPt>ptcut && muPlusPt>ptcut);*/
     //special case of a new tree in which I match only by filter (I assume this is equivalent as having HLTriggers fired)
     else if(isAArereco) bProcess =  (( (Reco_QQ_trig[iQQ]&(nTriggerBit+1))==(nTriggerBit+1) || (Reco_QQ_trig[iQQ]&nTriggerBit)==nTriggerBit )&& vProb>newVtxProbCut && muMinusPt>ptcut && muPlusPt>ptcut ); 
     else bProcess = true;
     if (bProcess)
       {
         if (i%1000==0) 	      cout << i << endl;
         DiMuonTree_allsign->Fill();// all sign and all mass
         if (QQsign==0) // opposite sign
      {
        DiMuonTree->Fill();// OS and all mass
        if (Reco_QQ->M()>mass_min && Reco_QQ->M()<mass_max) 
          {
            h_QQ_mass->Fill(Reco_QQ->M());// all upsilons in 7->14
          }
        if (Reco_QQ_mupl->Pt()>=ptcut && Reco_QQ_mumi->Pt()>=ptcut && Reco_QQ->M()>mass_min && Reco_QQ->M()<mass_max) 
          {
            h_QQ_mass_1->Fill(Reco_QQ->M()); // all OS upsilons in 7->14 and pt_mu>4GeV/c
          }
      }//opposite sign
         else // same sign in the acceptance
      if (Reco_QQ_mupl->Pt()>=ptcut && Reco_QQ_mumi->Pt()>=ptcut && Reco_QQ->M()>mass_min && Reco_QQ->M()<mass_max) 
        {
          h_QQ_mass_2->Fill(Reco_QQ->M());// all SS upsilons in 7->14 and pt_mu>4GeV/c
        }
         

     //--------------------------------------------------------
     // %%%%%%%%%% track rotation: redefine some of the variables, the others remain the same
     Double_t ran = gRandom->Rndm();
     
     if(ran < 0.5 ) RmuPlusPhi = muPlusPhi + TMath::Pi();
     else  RmuMinusPhi = muMinusPhi + TMath::Pi();
       
     TLorentzVector mu1;
     mu1.SetPtEtaPhiM( muPlusPt, muPlusEta, RmuPlusPhi, 0.105);
     TLorentzVector mu2;
     mu2.SetPtEtaPhiM( muMinusPt, muMinusEta, RmuMinusPhi, 0.105);
     
     TLorentzVector dimuon;
     dimuon = mu1 + mu2;
     
     invariantMass = dimuon.M();
     dimuPt         = dimuon.Pt();
     dimuEta        = dimuon.Eta();
     dimuRapidity   = dimuon.Rapidity();

     DiMuonTree_trkRot->Fill();
     
   }// if bProcess
      }// for each QQ pair		
    }// for each event in the tree
  
  //  __________________________________________________________________
  // ###### plotting 
  TH1 * phMAxis = new TH1D("phMAxis",";m_{#mu#mu} [GeV/c^{2}];Events/(0.1 GeV/c^{2})",1,mass_min,mass_max);
  phMAxis->SetDirectory(0);
  phMAxis->SetMinimum(1);
  phMAxis->SetMaximum(4e4);

  TCanvas *c1 = new TCanvas("c1","c1");
  phMAxis->Draw();
  
  //h_QQ_mass->GetYaxis()->SetRangeUser(0,180);
  h_QQ_mass->SetMarkerColor(kRed);
  h_QQ_mass->SetMarkerStyle(22);
  h_QQ_mass->Draw("PEsame");
  
  h_QQ_mass_1->SetMarkerColor(kBlue);
  h_QQ_mass_1->SetMarkerStyle(20);
  h_QQ_mass_1->GetXaxis()->CenterTitle(kTRUE);
  h_QQ_mass_1->Draw("PEsame");
  
  //h_QQ_mass_2->GetYaxis()->SetRangeUser(0,400);
  h_QQ_mass_2->SetMarkerColor(kRed);
  h_QQ_mass_2->SetMarkerStyle(24);
  h_QQ_mass_2->GetXaxis()->CenterTitle(kTRUE);
  h_QQ_mass_2->Draw("PEsame");
  
  TLegend *legend = new TLegend(.4,.7,.8,.9);
  legend->SetTextSize(0.025);
  legend->AddEntry(h_QQ_mass,  Form("OS && M [%.1f,%.1f]GeV: %.0f",mass_min, mass_max,h_QQ_mass->Integral(1,nBins)),"P");
  legend->AddEntry(h_QQ_mass_1,Form("OS && M [%.1f,%.1f]GeV && p_{T}>%.1fGeV/c: %.0f",mass_min, mass_max,ptcut,h_QQ_mass_1->Integral(1,nBins)),"P");
  legend->AddEntry(h_QQ_mass_2,Form("SS && M [%.1f,%.1f]GeV && p_{T}>%.1fGeV/c:%.0f",mass_min, mass_max,ptcut,h_QQ_mass_2->Integral(1,nBins)),"P");

  legend->Draw();
  
  h_QQ_mass->Write();
  h_QQ_mass_1->Write();
  h_QQ_mass_2->Write();
  c1->SaveAs(Form("dimuonDistribution_%s_Run%s_trigBit%d_allTriggers%d.pdf",dataSource,runNumber,nTriggerBit,bAllTriggers));
  c1->SaveAs(Form("dimuonDistribution_%s_Run%s_trigBit%d_allTriggers%d.gif",dataSource,runNumber,nTriggerBit,bAllTriggers));
  c1->Write();
  f1->Write();
  f1->Close();
}


