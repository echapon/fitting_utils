#ifndef allFunctions_h
#define allFunctions_h

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooCBShape.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooArgList.h"

#include <iostream>
const double mass_l =  7.5;
const double mass_h = 14.0;
const double binw   = 0.1;    //bin width of the histogram
const double M1S = 9.460;   //upsilon 1S pgd mass value
const double M2S = 10.023;  //upsilon 2S pgd mass value
const double M3S = 10.355;  //upsilon 3S pgd mass value
void makeWorkspace(RooWorkspace& w, int chooseSample ,float muonEtaMin, float muonEtaMax , float muonPtMin, float muonPtMax , float upsRapStart, float upsRapEnd, float upsPtStart, float upsPtEnd, int upsCentralityStart, int upsCentralityEnd);
void buildModel(RooWorkspace& w,int chooseFitParams, int chooseSample, int whatBin, int signalModel, int bkgModel, int doRap, int doPt,int doCent,int useRef,float muonPtMin,int fsr);
//int whatBin(float upsPtMin, float upsPtMax, float upsRapMin, float upsRapMax, int centStart, int centEnd, int doRap,int doPt, int doCent);
void fitAndDraw(RooWorkspace& w, int chooseSample, TString fname, TString datDir ,  TString figDir, int fitparams, int bkgModel,float muonPt1, float muonPt2,int centMin, int centMax, float upsPtStart, float upsPtEnd, float upsRapStart, float upsRapEnd,int doRap, int doPt,int doCent);
pair<double, double> ConfidencInterval(float CI, RooRealVar *fnll, RooDataSet *data, RooAbsPdf *pdf);
double computeSingleError(RooRealVar& x, RooRealVar& y, double correlation );
double computeSingle(RooRealVar& x, RooRealVar& y);


#endif // #ifndef allFunctions_h
