using namespace RooFit;

void makeWorkspace(RooWorkspace& ws, int ChooseSample,float muonEtaMin, float muonEtaMax, float muonPtMin, float muonPtMax, float upsRapMin, float upsRapMax, float upsPtMin,float upsPtMax, int upsCentralityStart, int upsCentralityEnd ){
  double mass_l =  7.0;
  double mass_h = 14.0;
 
  cout<< "You're in makeWorkspace!" << endl;
  std::string finput;
  switch (ChooseSample) 
    {
    case 1: // test
      finput   = "/afs/cern.ch/user/e/echapon/workspace/public/RunPrep2015/fitting_utils/dimuonTree_MCtest_upsi1S_5TeV_3.8T_RunMCtest_trigBit2_allTriggers1.root";
      mass_l = 8.5;
      mass_h = 11.5;
      binw=0.05;
      break;
    default:
      cout<<"You don't know what you are doing! Pick one of the available datasets in the choseSampleCases[] array"<<endl;
      break;
    }
  cout << finput << endl;
  //  TFile f(finput,"read");
  TFile *f = new TFile(finput.c_str());
  TTree* theTree       = (TTree*)gROOT->FindObject("UpsilonTree"); // OS --- all mass
  //  //RooWorkspace* ws = new RooWorkspace("ws","Upsilon Mass Fitting");
  //ws.var("invariantMass");	
  RooRealVar* mass       = new RooRealVar("invariantMass","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}");	
  //  ws.import(*mass);
  RooRealVar* upsPt      = new RooRealVar("upsPt","p_{T}(#Upsilon)",0,60,"GeV/c");
  RooRealVar* upsRapidity= new RooRealVar("upsRapidity",  "upsRapidity",-2.4, 2.4);
  RooRealVar* vProb      = new RooRealVar("vProb",  "vProb"  ,0.01,1.00);
  RooRealVar* QQsign     = new RooRealVar("QQsign",  "QQsign"  ,-1,5);
  RooRealVar* Centrality = new RooRealVar("Centrality","Centrality",0,200);
  RooRealVar* muPlusPt   = new RooRealVar("muPlusPt","muPlusPt",muonPtMin,1000);
  RooRealVar* muMinusPt  = new RooRealVar("muMinusPt","muMinusPt",muonPtMin,1000);
  RooRealVar* muPlusEta  = new RooRealVar("muPlusEta","muPlusEta",  muonEtaMin,muonEtaMax);
  RooRealVar* muMinusEta = new RooRealVar("muMinusEta","muMinusEta",muonEtaMin,muonEtaMax);
  RooDataSet* data0, *data;
  RooArgSet cols(*mass,*upsPt,*upsRapidity,*Centrality,*muPlusPt,*muMinusPt,*muPlusEta,*muMinusEta);
  data0 = new RooDataSet("data","data",theTree,cols); 
  data0->Print();
  // //  }
  // //data = ( RooDataSet*) data0->reduce(EventRange(0,100000));//,Cut("invariantMass<11"));
  TString cut_ap(Form("(%d<=Centrality && Centrality<%d) && (%.2f<muPlusEta && muPlusEta < %.2f) && (%.2f<muMinusEta && muMinusEta < %.2f) && (%.2f<upsPt && upsPt<%.2f) &&(abs(upsRapidity)>%.2f && abs(upsRapidity)<%.2f)  &&((muPlusPt > %.2f && muMinusPt > %.2f) || (muPlusPt > %.2f && muMinusPt > %.2f))",upsCentralityStart,upsCentralityEnd,muonEtaMin,muonEtaMax,muonEtaMin,muonEtaMax,upsPtMin,upsPtMax,upsRapMin,upsRapMax,muonPtMin,muonPtMax,muonPtMax,muonPtMin));
  cout << cut_ap << endl;
  data =  ( RooDataSet*) data0->reduce(Cut(cut_ap));
  ws.import(*data);
  data->Print();
  f->Close();
}
