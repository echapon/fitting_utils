#include "makeWorkspace2015.C"
#include "buildModelUpsi2015.C"

void fit2015() {
   RooWorkspace myws;
   makeWorkspace2015(myws, 1, -2.4, 2.4, 0, 0, 0, 2.4, 0, 50, 0, 200);

   buildModelUpsi2015(myws, 2, 3);

   RooRealVar* mass =(RooRealVar*) myws.var("invariantMass"); //
   RooDataSet* data0_fit =(RooDataSet*) myws.data("data");
   RooAbsPdf* pdf = (RooAbsPdf*) myws.pdf("pdf");
   RooFitResult* fitObject = pdf->fitTo(*data0_fit,Save(),Hesse(kTRUE),Extended(kTRUE)); // Fit
   pdf->Print();
   Double_t baseNll = fitObject->minNll();
   // RooMinuit m(*nll);
   RooRealVar *nsig1f   =(RooRealVar*) myws.var("N_{ #varUpsilon(1S)}");

   RooArgSet allvars = myws.allVars();
   int npars= allvars.getSize() ;
   int nbins = ceil((mass_h-mass_l)/binw); // all of which are 'globally' defined in the .h header
   RooPlot* frame = mass->frame(Bins(nbins),Range(mass_l,mass_h));  
   data0_fit->plotOn(frame);// data drawn first for pdf object to pick the proper normalisation
   pdf->plotOn(frame,Name("thePdf"));

   frame->Draw();
}
