#include "allFunctions.h"

using namespace RooFit;

void makeWorkspace2015(RooWorkspace& ws, const char* filename,float mass_l, float mass_h, float muonEtaMin, float muonEtaMax, float muonPtMin, float dimuRapMin, float dimuRapMax, float dimuPtMin,float dimuPtMax, int dimuCentralityStart, int dimuCentralityEnd ){
   double binw=0.05;
 
  cout<< "You're in makeWorkspace2015!" << endl;
  std::string finput(filename);
  cout << finput << endl;
  //  TFile f(finput,"read");
  TFile *f = new TFile(finput.c_str());
  TTree* theTree       = (TTree*)gROOT->FindObject("DiMuonTree"); // OS --- all mass
  cout << theTree << endl;
  //  //RooWorkspace* ws = new RooWorkspace("ws","DiMuon Mass Fitting");
  //ws.var("invariantMass");	
  RooRealVar* mass       = new RooRealVar("invariantMass","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}");	
  //  ws.import(*mass);
  RooRealVar* dimuPt      = new RooRealVar("dimuPt","p_{T}(#DiMuon)",0,60,"GeV/c");
  RooRealVar* dimuRapidity= new RooRealVar("dimuRapidity",  "dimuRapidity",-2.4, 2.4);
  RooRealVar* vProb      = new RooRealVar("vProb",  "vProb"  ,0.01,1.00);
  RooRealVar* QQsign     = new RooRealVar("QQsign",  "QQsign"  ,-1,5);
  RooRealVar* Centrality = new RooRealVar("Centrality","Centrality",0,200);
  RooRealVar* muPlusPt   = new RooRealVar("muPlusPt","muPlusPt",muonPtMin,1000);
  RooRealVar* muMinusPt  = new RooRealVar("muMinusPt","muMinusPt",muonPtMin,1000);
  RooRealVar* muPlusEta  = new RooRealVar("muPlusEta","muPlusEta",  muonEtaMin,muonEtaMax);
  RooRealVar* muMinusEta = new RooRealVar("muMinusEta","muMinusEta",muonEtaMin,muonEtaMax);
  RooDataSet* data0, *data;
  RooArgSet cols(*mass,*dimuPt,*dimuRapidity,*Centrality,*muPlusPt,*muMinusPt,*muPlusEta,*muMinusEta);
  data0 = new RooDataSet("data","data",theTree,cols); 
  data0->Print();
  // //  }
  // //data = ( RooDataSet*) data0->reduce(EventRange(0,100000));//,Cut("invariantMass<11"));
  TString cut_ap(Form("(%d<=Centrality && Centrality<%d) &&" 
           "(%.2f<muPlusEta && muPlusEta < %.2f) &&" 
           "(%.2f<muMinusEta && muMinusEta < %.2f) &&" 
           "(%.2f<dimuPt && dimuPt<%.2f) &&"
           "(abs(dimuRapidity)>%.2f && abs(dimuRapidity)<%.2f)  &&"
           "(muPlusPt > %.2f && muMinusPt > %.2f)",
           dimuCentralityStart,dimuCentralityEnd,muonEtaMin,muonEtaMax,muonEtaMin,muonEtaMax,dimuPtMin,dimuPtMax,dimuRapMin,dimuRapMax,muonPtMin,muonPtMin));
  cout << cut_ap << endl;
  data =  ( RooDataSet*) data0->reduce(Cut(cut_ap));
  ws.import(*data);
  data->Print();
  f->Close();
}
