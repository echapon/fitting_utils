#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

//RooFit Stuff
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooKeysPdf.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooChi2Var.h"
//RooStats
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
//fitUpsilonYields package
#include "makeWorkspace.C"
#include "buildModel.C"
//#include "whatBin.C"
#include "fitAndDraw.C"
#include "allFunctions.h"
// Root stuff
#include "TROOT.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TLatex.h"

#include "TMath.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"

//bkgtable
#define nData  9
#include "bkgTable.h"
#include "dataTable.h"
// miscellaneous  
#include <fstream>
#include <new>
#include <iostream>




// -------- make some local picks
const bool doMinos      = 1;     //kFALSE;
using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace RooStats;
pair<double, double> ConfidencInterval(float, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf);
//
void fitUpsilonYields(int ChooseSample =1, 
		      int bkgModel=3, 
		      int fsr=1,
		      int useRef =2,//1: data-driven, 0: none fixed, 2: FSR from MC(standard case in 15-001)
		      int whatBin =0,
		      int doRap=0,
		      int doPt=0,
		      float muonPtCut1=3.5,
		      float muonPtCut2 = 4,  
		      float muonEtaMin =-2.4,
		      float muonEtaMax =2.4,
		      float upsRapStart = 0.0, // it's the absolute value (default) for y bins are symmetric in pp and pbpb.. must be over-ridden in pPb...
		      float upsRapEnd = 2.4,
		      float upsPtStart = 0,
		      float upsPtEnd = 40,
		      int upsCentralityStart = 0,
		      int upsCentralityEnd =40,
		      int signalModel = 4,
		      TString outDatDir = "txtOutput",// dats dir outfile location
		      TString outFigDir = "pdfOutput",
		      const char* ChooseSampleCase = "pbpb",
		      const char* outFilePrefix   = "yields", // yields, ratios
		      int ChooseFitParams = 0
		      ){
  cout << "welcome to FitUpsilonYields!"<<endl;
  gROOT->Macro("~/Project/ups2013/code/cm/logon.C+");
  // Create a workspace named 'w' that exports its contents to a same-name C++ namespace in CINT 'namespace w' !!!
  // THIS FEATURE is _CINT_ only, don't try to compile this code.
  RooWorkspace* _ws = new RooWorkspace("_ws",kTRUE);
  int doCent=0; if(ChooseSample==8){fsr=0;}
  if(doPt==0 && doRap==0) doCent=1; 
  //make the workspace where the job will be contained.
  makeWorkspace(*_ws,ChooseSample, muonEtaMin, muonEtaMax, muonPtCut1,muonPtCut2, upsRapStart, upsRapEnd, upsPtStart, upsPtEnd,upsCentralityStart,upsCentralityEnd);
  RooDataSet* data0_fit =(RooDataSet*) _ws->data(); // Use the name space prefix operator to access the workspace contents
  // data0_fit->Print(); // to check it worked.
  buildModel(*_ws,ChooseFitParams,ChooseSample,whatBin,signalModel,bkgModel,doRap,doPt,doCent,useRef,muonPtCut1,fsr); // builds the pdf and imports it to the rooworkspace
  // _ws->Print(); //still checkin it worked.
  //  RooFitResult* fitObject = _ws::pdf.fitTo(*data0_fit); // Fit gets into action.
  //_ws::pdf.Print("v");
  // write out the fitting params
  int centMin =2.5*upsCentralityStart;
  int centMax =2.5*upsCentralityEnd;

  TString figName_(Form("%s_%s_fsr%d",ChooseSampleCase,outFilePrefix,fsr));
  //  TString figName_(Form("%s_%s_cent%d%d_bkgModel%d_muonEta%.2f%.2f_muonPt%.2f-%.2f_dimuPt%.2f%.2f_dimuY%.2f%.2f_%d_ref%d_mass8p511p5",outFilePrefix,ChooseSampleCase,centMin,centMax,bkgModel,muonEtaMin,muonEtaMax,muonPtCut1,muonPtCut2,upsPtStart,upsPtEnd,upsRapStart,upsRapEnd,fsr,useRef));
  figName_.ReplaceAll("-","M");
  figName_.ReplaceAll(".","");
  fitAndDraw(*_ws,ChooseSample,figName_,outDatDir,outFigDir,ChooseFitParams,bkgModel,muonPtCut1,muonPtCut2,centMin,centMax,upsPtStart,upsPtEnd,upsRapStart,upsRapEnd,doRap,doPt,doCent); 
}

//_______________________________________________________________________
double computeSingle(RooRealVar& x, RooRealVar& y) 
{
  // pass the 1S(x) and xS/1S(y) ratios and calcualte the xS yield
  return x.getVal() * y.getVal();
}


double computeSingleError(RooRealVar& x, RooRealVar& y, double correlation ) 
{
  // pass the 1S(x) and xS/1S(y) ratios and calcualte the xS yield error
  
  double err2 = (x.getError()*x.getError())/(x.getVal()*x.getVal()) 
    + (y.getError()*y.getError())/(y.getVal()*y.getVal()) 
    + 2.*(x.getError()*y.getError())/(x.getVal()*y.getVal())*correlation;
  
  return fabs(computeSingle(x,y))*sqrt(err2);
}


//_______________________________________________________________________


//calculate the confidence interval with RooStats
pair<double, double> ConfidencInterval(float CI, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf)  {  
	ProfileLikelihoodCalculator pl(*data,*pdf,*fnll);
	pl.SetConfidenceLevel(CI); 
	int ci = 100*CI;
	LikelihoodInterval* interval = pl.GetInterval();
	LikelihoodIntervalPlot plot(interval);
	TCanvas c4; c4.cd(); 
	plot.SetRange(0.,3.,-.05,0.5);
	plot.Draw();
	TLatex latexCI;
	latexCI.SetNDC();
	latexCI.SetTextSize(0.035);
	latexCI.DrawLatex(0.5,1.-0.05*2,Form("R_{#frac{3S}{1S}} %d %% C.I.",ci));
	latexCI.DrawLatex(0.5,1.-0.05*3,Form("Upper limit: %f",interval->UpperLimit(*fnll)));
	TString intrvlName = fnll->GetTitle();
		// print out the iterval on the Parameter of Interest
	cout <<endl<< CI <<"\% interval on " <<fnll->GetName()<<" is : ["<<
		interval->LowerLimit(*fnll) << ", "<<
		interval->UpperLimit(*fnll) << "] "<<endl;
	pair<double, double> CnfdncIntrvl;
	CnfdncIntrvl.first  = interval->LowerLimit(*fnll);
	CnfdncIntrvl.second = interval->UpperLimit(*fnll);
	c4.SaveAs("pdfOutput/1806/HIMBProfileTest.pdf");

	return CnfdncIntrvl;
}

