#include "allFunctions.h"

void buildModel(RooWorkspace& w,int chooseFitParams, int chooseSample,int whatBin, int signalModel, int bkgdModel, int doRap, int doPt,int doCent,int useRef,float muonPtMin, int fixFSR){
// C r e a t e   m o d e l  
  int nt=100000;
  // cout << "you're building a model for the quarkonium resonance of mass = "<< M1S <<" GeV/c^{2},"endl;
  RooRealVar *nsig1f   = new RooRealVar("N_{ #varUpsilon(1S)}","nsig1S",0,nt*10);
  RooRealVar* mass       = new RooRealVar("invariantMass","#mu#mu mass",mass_l,mass_h,"GeV/c^{2}");

  RooRealVar *nsig2f = NULL;
  RooRealVar *nsig3f = NULL;

 switch (chooseFitParams)
    {
    case 0://use the YIELDs of 2S and 3S as free parameters
      //minor modif here: 3S forced positive.
      nsig2f  = new RooRealVar("N_{ #varUpsilon(2S)}","nsig2S",   nt*0.25,-200,10*nt);
      nsig3f  = new RooRealVar("N_{ #varUpsilon(3S)}","nsig3S",   nt*0.25,-200,10*nt);
      cout << "you're fitting to extract yields, "<< endl;
      break;
    default:
      cout<<"Make a pick from chooseFitParams!!!"<<endl;
      break;
    }  

      RooRealVar  *mean = new RooRealVar("m_{ #varUpsilon(1S)}","#Upsilon mean",M1S,M1S-0.2,M1S+0.2);
      RooConstVar *rat2 = new RooConstVar("rat2", "rat2", M2S/M1S);
      RooConstVar *rat3 = new RooConstVar("rat3", "rat3", M3S/M1S);
      // scale mean and resolution by mass ratio
      RooFormulaVar *mean1S = new RooFormulaVar("mean1S","@0",RooArgList(*mean));
      RooFormulaVar *mean2S = new RooFormulaVar("mean2S","@0*@1", RooArgList(*mean,*rat2));
      RooFormulaVar *mean3S = new RooFormulaVar("mean3S","@0*@1", RooArgList(*mean,*rat3));

      // //detector resolution ?? where is this coming from?
      RooRealVar    *sigma1  = new RooRealVar("#sigma_{CB1}","#sigma_{CB1}",sigma_min[whatBin],sigma_max[whatBin]); // 
      RooFormulaVar *sigma1S = new RooFormulaVar("sigma1S","@0"   ,RooArgList(*sigma1));
      RooFormulaVar *sigma2S = new RooFormulaVar("sigma2S","@0*@1",RooArgList(*sigma1,*rat2));
      RooFormulaVar *sigma3S = new RooFormulaVar("sigma3S","@0*@1",RooArgList(*sigma1,*rat3));
      RooRealVar *alpha  = new RooRealVar("#alpha_{CB}","tail shift",alpha_min[whatBin],alpha_max[whatBin]);    // MC 5tev 1S pol2 
      RooRealVar *npow   = new RooRealVar("n_{CB}","power order",npow_min[whatBin],npow_max[whatBin]);    // MC 5tev 1S pol2 
      RooRealVar *sigmaFraction = new RooRealVar("sigmaFraction","Sigma Fraction",0.,1.);
      // scale the sigmaGaus with sigma1S*scale=sigmaGaus now.
      RooRealVar    *scaleWidth = new RooRealVar("#sigma_{CB2}/#sigma_{CB1}","scaleWidth",1.,2.5);
      RooFormulaVar *sigmaGaus = new RooFormulaVar("sigmaGaus","@0*@1", RooArgList(*sigma1,*scaleWidth));
      RooFormulaVar *sigmaGaus2 = new RooFormulaVar("sigmaGaus","@0*@1*@2", RooArgList(*sigma1,*scaleWidth,*rat2));
      RooFormulaVar *sigmaGaus3 = new RooFormulaVar("sigmaGaus","@0*@1*@2", RooArgList(*sigma1,*scaleWidth,*rat3));
      RooGaussian* gauss1 = new RooGaussian("gaus1s","gaus1s",
      					    *nsig1f,
      					    *mass,    //mean
      					    *sigmaGaus); //sigma
      // RooGaussian* gauss1b = new RooGaussian("gaus1sb","gaus1sb",
      // 					     *nsig1f,
      // 					     *m,    //mean
      // 					     *sigma1); //sigma
      switch(signalModel){    
      case 1: //crystal boule
	RooCBShape  *sig1S   = new RooCBShape ("cb1S_1", "FSR cb 1s",
					       *mass,*mean1S,*sigma1,*alpha,*npow);

	RooCBShape  *sig2S   = new RooCBShape ("cb2S_1", "FSR cb 1s",
					       *mass,*mean2S,*sigma2S,*alpha,*npow);
	RooCBShape  *sig3S   = new RooCBShape ("cb3S_1", "FSR cb 1s",
					       *mass,*mean3S,*sigma3S,*alpha,*npow);
	cout << "you're fitting each signal peak with a Crystal Ball function"<< endl;
	break;
      case 2: //Gaussein
	RooAbsPdf      *sig1S  = new RooGaussian ("g1", "gaus 1s",
						  *mass,*mean1S,*sigma1);
	cout << "you're fitting 1 signal peak with a Gaussian function"<< endl;
	break;
      case 3: //Gaussein + crystal boule
	RooCBShape  *cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
						 *mass,*mean1S,*sigma1,*alpha,*npow);
	RooAddPdf      *sig1S  = new RooAddPdf ("cbg", "cbgaus 1s",
						RooArgList(*gauss1,*cb1S_1),*sigmaFraction);
	cout << "you're fitting 1 signal peak with a sum of a Gaussian and a Crystal Ball function"<< endl;
	break;
      case 4: //crystal boules
	RooCBShape  *cb1S_1    = new RooCBShape ("cb1S_1", "FSR cb 1s",
						 *mass,*mean1S,*sigma1,*alpha,*npow);
       
	RooCBShape  *cb1S_2    = new RooCBShape ("cb1S_2", "FSR cb 1s",
						 *mass,*mean1S,*sigmaGaus,*alpha,*npow);
	RooAddPdf      *sig1S  = new RooAddPdf  ("cbcb","1S mass pdf",
						 RooArgList(*cb1S_1,*cb1S_2),*sigmaFraction);
	// /// Upsilon 2S
	RooCBShape  *cb2S_1    = new RooCBShape ("cb2S_1", "FSR cb 2s", 
						 *mass,*mean2S,*sigma2S,*alpha,*npow); 
	RooCBShape  *cb2S_2    = new RooCBShape ("cb2S_2", "FSR cb 2s", 
						 *mass,*mean2S,*sigmaGaus2,*alpha,*npow); 
	RooAddPdf      *sig2S  = new RooAddPdf  ("sig2S","2S mass pdf",
						 RooArgList(*cb2S_1,*cb2S_2),*sigmaFraction);
      
	// /// Upsilon 3S
	RooCBShape  *cb3S_1    = new RooCBShape ("cb3S_1", "FSR cb 3s", 
						 *mass,*mean3S,*sigma3S,*alpha,*npow); 
	RooCBShape  *cb3S_2    = new RooCBShape ("cb3S_2", "FSR cb 3s", 
						 *mass,*mean3S,*sigmaGaus3,*alpha,*npow); 
	RooAddPdf      *sig3S  = new RooAddPdf  ("sig3S","3S mass pdf",
						 RooArgList(*cb3S_1,*cb3S_2),*sigmaFraction); // = cb3S1*sigmaFrac + cb3S2*(1-sigmaFrac)
	cout << "you're fitting each signal peak with a Double Crystal Ball function"<< endl;
	break;
      
      case 5: //deux Gausseins
	RooAddPdf      *sig1S  = new RooAddPdf ("cb1S_1", "cbgaus 1s",
						RooArgList(*gauss1,*gauss1b),*sigmaFraction);
	cout << "you're fitting each signal peak with a Double Gaussian function"<< endl;
	break;
      }
      // bkg Chebychev
      RooRealVar *nbkgd   = new RooRealVar("n_{Bkgd}","nbkgd",0,nt);
      RooRealVar *bkg_a1  = new RooRealVar("a1_bkg", "bkg_{a1}", 0, -5, 5);
      RooRealVar *bkg_a2  = new RooRealVar("a2_Bkg", "bkg_{a2}", 0, -2, 2);
      RooRealVar *bkg_a3  = new RooRealVar("a3_Bkg", "bkg_{a3}", 0, -0.9, 2);

      //  likesign
      RooRealVar *nLikesignbkgd = new RooRealVar("NLikesignBkg","nlikesignbkgd",nt*0.75,0,10*nt);
      // *************************************************** bkgModel
      
      RooRealVar turnOn("turnOn","turnOn", turnOn_minCent[whatBin],turnOn_maxCent[whatBin]);
      RooRealVar width("width","width",width_minCent[whatBin],width_maxCent[whatBin]);// MB 2.63
      RooRealVar decay("decay","decay",decay_minCent[whatBin],decay_maxCent[whatBin]);// MB: 3.39
      if (doRap && !doPt)
	{
	  RooRealVar turnOn("turnOn","turnOn", turnOn_minRap[whatBin],turnOn_maxRap[whatBin]);
	  RooRealVar width("width","width",width_minRap[whatBin],width_maxRap[whatBin]);// MB 2.63
	  RooRealVar decay("decay","decay",decay_minRap[whatBin],decay_maxRap[whatBin]);// MB: 3.39
	}
      if (doPt && !doRap)
	{
	  RooRealVar turnOn("turnOn","turnOn", turnOn_minPt[whatBin],turnOn_maxPt[whatBin]);
	  RooRealVar width("width","width",width_minPt[whatBin],width_maxPt[whatBin]);// MB 2.63
	  RooRealVar decay("decay","decay",decay_minPt[whatBin],decay_maxPt[whatBin]);// MB: 3.39
	}
      
  
      width.setConstant(false);
      decay.setConstant(false);
      turnOn.setConstant(false);
 
  switch (useRef)// no reference
    {
    case 0: // forcing sigma and fsr to be left free.
      fixSigma1 = 0;
      fixFSR    = 0;
      break;
    case 1:
      //using data-driven estimates
      int dataRef=1;
      cout<<"You're using the debug mode based on data parameters. So you must not take this result as the central one."<<endl;
      break;
    case 2:
      
      cout << "doCent="<<doCent << endl;
      //using MC-driven estimates
      int dataRef=2;
      if(doCent)   //MB values, assumed to be the same with all centralities...
     	{
	  if(muonPtMin <4){
	    gROOT->LoadMacro("dataTable_loose.h");
	  }else if(muonPtMin > 3.5){
	    gROOT->LoadMacro("dataTable_tight.h");
	  }
	  npow->setVal(npow_rapBins[8]);
	  alpha->setVal(alpha_rapBins[8]);
	  sigma1->setVal(sigma1_rapBins[8]);
	  scaleWidth->setVal(scale_rapBins[8]);
	  sigmaFraction->setVal(pdFrac_rapBins[8]);
	  cout<< whatBin << endl;
	}
      if(doRap && !doPt)
	{
	  if(muonPtMin <4){
	    gROOT->LoadMacro("dataTable_loose.h");
	  }else if(muonPtMin > 3.5){
	    gROOT->LoadMacro("dataTable_tight.h");
	  }
	  npow->setVal(npow_rapBins[whatBin]);
	  alpha->setVal(alpha_rapBins[whatBin]);
	  sigma1->setVal(sigma1_rapBins[whatBin]);
	  scaleWidth->setVal(scale_rapBins[whatBin]);
	  sigmaFraction->setVal(pdFrac_rapBins[whatBin]);
	  cout<< whatBin << endl;
	}
      if(doPt && !doRap)
	{
	  // cout << "we're here" << endl;
	  if(muonPtMin <4){
	    gROOT->LoadMacro("dataTable_loose.h");
	  }else if(muonPtMin > 3.5){
	    gROOT->LoadMacro("dataTable_tight.h");
	  }
	  cout << " ok ... " <<endl;
	  npow->setVal(npow_ptBins[whatBin]);
	  alpha->setVal(alpha_ptBins[whatBin]);
	  sigma1->setVal(sigma1_ptBins[whatBin]);
	  scaleWidth->setVal(scale_ptBins[whatBin]);
	  sigmaFraction->setVal(pdFrac_ptBins[whatBin]);
	}
     
      cout<<"You're using MC parameters. So you may use this result as the central one, according to the LLR test outcome."<<endl;
      break;
    default: break;
    }

  //
  cout << "npow tried=" << npow->getVal(); 
  if(fixFSR==3 || fixFSR==1) cout << "   constant!" << endl;   else cout << "   floating!" << endl;
  cout << "alpha tried=" << alpha->getVal(); 
  if(fixFSR==2 || fixFSR==1) cout << "   constant!" << endl;   else cout << "   floating!" << endl;
  cout << "sigma1 tried=" << sigma1->getVal(); 
  if(fixFSR==4 || fixFSR==1) cout << "   constant!" << endl;   else cout << "   floating!" << endl;
  cout << "scale tried=" << scaleWidth->getVal(); 
  if(fixFSR==4 || fixFSR==1) cout << "   constant!" << endl;   else cout << "   floating!" << endl;
  cout << "normalisation tried=" << sigmaFraction->getVal();
  if(fixFSR==5 || fixFSR==1) cout << "   constant!" << endl;   else cout << "   floating!" << endl;
  
  switch (fixFSR) // 0: free;  1: both fixed 2: alpha fixed 3: npow fixed
    {      
    case 0:// all free
      alpha->setConstant(false);
      npow->setConstant(false);
      sigma1->setConstant(false);
      scaleWidth->setConstant(false);
      sigmaFraction->setConstant(false);
      break;
    case 1:// all fixed
      alpha->setConstant(true);
      npow ->setConstant(true);
      sigma1->setConstant(true);
      scaleWidth->setConstant(true);
      sigmaFraction->setConstant(true);
      break;
    case 2: // release alpha
      alpha->setConstant(false);
      npow ->setConstant(true);
      sigma1->setConstant(true);
      scaleWidth->setConstant(true);
      sigmaFraction->setConstant(true);
      break;
    case 3:// npow released
      alpha->setConstant(true);
      npow->setConstant(false);
      sigma1->setConstant(true);
      scaleWidth->setConstant(true);
      sigmaFraction->setConstant(true);
      break;
    case 4:// width+ sF +scale released
      alpha->setConstant(true);
      npow->setConstant(true);
      sigma1->setConstant(false);
      scaleWidth->setConstant(true);
      sigmaFraction->setConstant(true);
      break;
    case 5:// scale +sF
      alpha->setConstant(true);
      npow->setConstant(true);
      sigma1->setConstant(true);
      scaleWidth->setConstant(false);
      sigmaFraction->setConstant(true);
      break;
   case 6:// scale +sF
      alpha->setConstant(true);
      npow->setConstant(true);
      sigma1->setConstant(true);
      scaleWidth->setConstant(true);
      sigmaFraction->setConstant(false);
      break;
    default:
      cout<<"Donno this choice! Pick somehting for FSR parameters that I know"<<endl;
      break;
    }
   //thisPdf: form of the bkg pdf
  //pdf_combinedbkgd; // total bkg pdf. usually form*normalization  (so that you can do extended ML fits)
  switch (bkgdModel) 
    {
    case 1 :  //(erf*exp ) to fit the SS, then fix the shape and fit OS, in case of constrain option
      bkg_a3->setConstant(true);
      RooGenericPdf *ErrPdf = new  RooGenericPdf("ErrPdf","ErrPdf",
						 "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
					       RooArgList(*mass,turnOn,width,decay));
      RooFitResult* fit_1st  = ErrPdf->fitTo(*likesignData,Save(),NumCPU(4)) ; // likesign data
      if (doTrkRot) fit_1st  = thisPdf->fitTo(*TrkRotData,Save(),NumCPU(4)) ;
      
      if (doConstrainFit) 
	{ // allow parameters to vary within cenral value from above fit + their sigma
	  turnOn_constr = new RooGaussian("turnOn_constr","turnOn_constr",
					  turnOn,
					  RooConst(turnOn.getVal()),
					  RooConst(turnOn.getError()));
	  width_constr   = new RooGaussian("width_constr","width_constr",
					   width,					   RooConst(width.getVal()),
					   RooConst(width.getError()));
	  decay_constr    = new RooGaussian("decay_constr","decay_constr",
					    decay,
					    RooConst(decay.getVal()),
					    RooConst(decay.getError()));
	}
      else 
	{
	  turnOn.setConstant(kTRUE);
	  width.setConstant(kTRUE);
	  decay.setConstant(kTRUE);
	}
      RooRealVar *fLS =new RooRealVar("R_{SS/OS}","Empiric LS/SS ratio",0.,1.);
      RooAbsPdf  *ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
					      *mass, RooArgList(*bkg_a1,*bkg_a2));
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("bkgPdf","total combined background pdf",
						      RooArgList(*ErrPdf,*ChebPdf),
						      RooArgList(*fLS));
      
      break;
    case 2 : //us eRooKeysPdf to smooth the SS, then fit OS with pol+keys
      bkg_a3->setConstant(true);
      RooRealVar *fLS =new RooRealVar("R_{SS/OS}","Empiric LS/SS ratio",0.,1.);
      RooKeysPdf *KeysPdf        = new RooKeysPdf("KeysPdf","KeysPdf",*mass,*likesignData,
						  RooKeysPdf::MirrorBoth, 1.4);
      RooAbsPdf  *ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
					      *mass, RooArgList(*bkg_a1,*bkg_a2));
      if (doTrkRot) thisPdf     = new RooKeysPdf("thisPdf","thisPdf",*mass,*TrkRotData,
						 RooKeysPdf::MirrorBoth, 1.4);
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("bkgPdf","total combined background pdf",
						      RooArgList(*KeysPdf,*ChebPdf),
						      RooArgList(*fLS));
      break;
    case 3 : //use error function to fit the OS directly
      bkg_a3->setConstant(true);
      RooAbsPdf *pdf_combinedbkgd            = new  RooGenericPdf("bkgPdf","bkgPdf",
								  "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
								  RooArgList(*mass,turnOn,width,decay));
      break;
      
    case 4 : //use pol 2+ErfExp to fit the OS directly

      RooRealVar *fPol   = new RooRealVar("F_{pol}","fraction of polynomial distribution",0.0,1);
      RooAbsPdf  *ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
					      *mass, RooArgList(*bkg_a1,*bkg_a2));
      RooGenericPdf *ErrPdf     = new  RooGenericPdf("ErrPdf","ErrPdf",
						     "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
						     RooArgList(*mass,turnOn,width,decay));
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("bkgPdf","total combined background pdf",
						      RooArgList(*ChebPdf,*ErrPdf),
						      RooArgList(*fPol));

      break;
    case 5 : //use ( error function + polynomial 1) to fit the OS directly
      
      bkg_a3->setConstant(true);
      bkg_a2->setConstant(true);
      RooRealVar *fPol   = new RooRealVar("F_{pol}","fraction of polynomial distribution",0.0,1);
      RooAbsPdf  *ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
					      *mass, RooArgList(*bkg_a1,*bkg_a2,*bkg_a3));
      RooGenericPdf *ErrPdf     = new  RooGenericPdf("ErrPdf","ErrPdf",
						     "exp(-@0/decay)*(TMath::Erf((@0-turnOn)/width)+1)",
						   RooArgList(*mass,turnOn,width,decay));
      RooAbsPdf  *pdf_combinedbkgd   = new RooAddPdf ("bkgPdf","total combined background pdf",
						      RooArgList(*ChebPdf,*ErrPdf),
						      RooArgList(*fPol));
      break;
    case 6: // NOT WORKING
      RooRealVar *fPol   = new RooRealVar("F_{pol}","fraction of polynomial distribution",0.0,1);
      RooAbsPdf  *ChebPdf  = new RooChebychev("ChebPdf","ChebPdf",
					      *mass, RooArgList(*bkg_a1,*bkg_a2));
      RooGenericPdf *ExpPdf     = new  RooGenericPdf("ExpPdf","ExpPdf",
						     "exp(-@0/decay)",
						     RooArgList(*mass,decay));
      RooAbsPdf  *pdf_combinedbkgd  =  new RooAddPdf ("bkgPdf","total combined background pdf",
						      RooArgList(*ChebPdf,*ExpPdf),
						      RooArgList(*fPol));
      break;
    default :
      cout<<"Donno what you are talking about! Pick another fit option!"<<endl;
      break;
    }
  
  //###### the nominal fit with default pdf 
 
  // RooAbsPdf  *pdf; // nominal PDF
  if(chooseSample==8)
   {     
     // bkg_a1->setVal(0);// can be turned on at convenience
     // bkg_a1->setConstant();
     // bkg_a2->setVal(0);
     // bkg_a2->setConstant();
     // bkg_a3->setVal(0);
     // bkg_a3->setConstant();
     // RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total p.d.f.",
     // 						  RooArgList(*sig1S,*pdf_combinedbkgd),
     // 						  RooArgList(*nsig1f,*nbkgd));
        RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total p.d.f.",*sig1S,*nsig1f);
   } else if(chooseSample!=8)
    {
      // can remove the double crystal ball in pbpb: just commenting out and copying an appropriate version
      RooAbsPdf  *pdf             = new RooAddPdf ("pdf","total p.d.f.",
						   RooArgList(*sig1S,*sig2S,*sig3S,*pdf_combinedbkgd),
						   RooArgList(*nsig1f,*nsig2f,*nsig3f,*nbkgd));
      //  nsig3f->setVal(0); nsig3f->setConstant();
      
    }
  w.import(*pdf);
  w.Print();
}
