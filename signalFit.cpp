/*
c++ `root-config --cflags --glibs` -lrooFit -lrooFitCore -o signalFit signalFit.cpp Functions.cc

TBA
http://root.cern.ch/root/html/src/RooCBShape.cxx.html#Uso4DD

*/

#include "TProfile2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TNtuple.h"
#include <fstream>
#include <vector>

#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooFFTConvPdf.h"
#include "RooKeysPdf.h"
#include "RooExponential.h"
#include "RooArgusBG.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooFFTConvPdf.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif 

#include "Functions.h"

using namespace RooFit ;
using namespace std ;


int main (int argc, char ** argv)
{
  gROOT->SetStyle ("Plain") ;	
  gStyle->SetOptStat ("mr") ;
  gStyle->SetOptFit (1111) ;
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.1);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.06);
  gStyle->SetStatW(0.3);
  gStyle->SetPalette(1);


// create pdfs  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  // variables definition
  RooRealVar x ("x", "x", 200, 2000) ; 

  RooRealVar mean_gaus ("mean_gaus", "mean of gaussian", 500., 200., 1500.) ; 
  RooRealVar sigma_gaus ("sigma_gaus", "width of gaussian", 20., 0., 500.) ; 
  RooGaussian pdf_gaus ("pdf_gaus", "gaussian PDF", x, mean_gaus, sigma_gaus) ; 

  RooRealVar mean_brwi ("mean_brwi", "mean of breit wigner", 500., 200., 1500.) ; 
  RooRealVar sigma_brwi ("sigma_brwi", "width of breit wigner", 20., 0., 500.) ; 
  RooBreitWigner pdf_brwi ("pdf_brwi", "breit wigner PDF", x, mean_brwi, sigma_brwi) ; 

  RooFFTConvPdf model ("model", "pdf_brwi (X) pdf_gaus", x, pdf_brwi, pdf_gaus) ;

  RooRealVar mean_cb ("mean_cb", "mean of crystal ball", 500., 200., 1500.) ; 
  RooRealVar sigma_cb ("sigma_cb", "width of crystal ball", 20., 0., 500.) ; 
  RooRealVar alpha_cb ("alpha_cb", "width of crystal ball", 20., 0., 500.) ; 
  RooRealVar power_cb ("power_cb", "width of crystal ball", 20., 0., 500.) ; 
  RooCBShape pdf_cb ("pdf_cb", "CB PDF", x, mean_cb, sigma_cb, alpha_cb, power_cb) ; 

//  // Create a ROOT TF1 function
//  TF1 * tf1_doubleCB = new TF1 ("tf1_doubleCB", crystalBallLowHigh, 0, 2000, 7) ;
//  RooAbsReal* rf_doubleCB = bindFunction (tf1_doubleCB, x) ;

  vector<double> masses ;
  masses.push_back (350) ;
  masses.push_back (500) ;
  masses.push_back (650) ;
  masses.push_back (800) ;
  masses.push_back (1000) ;
  
  vector<string> filenames ;
  filenames.push_back ("findInterference.350.root") ;
  filenames.push_back ("findInterference.500.root") ;
  filenames.push_back ("findInterference.650.root") ;
  filenames.push_back ("findInterference.800.root") ;
  filenames.push_back ("findInterference.1000.root") ;

  vector<TFile *> inputfiles ;
  vector<TH1F *> h_mg_signals ;
  vector<RooDataHist *> rdh_mg_signals ;
  vector<RooPlot *> c_frames ;
  
  for (int i = 0 ; i < filenames.size () ; ++i)
    {
      cout << "------------------------------------------------------------------------------------------\n\n" ;
      cout << "    working on mass " << masses.at (i) << "\n" ;
      cout << "\n------------------------------------------------------------------------------------------\n\n" ;
    
      //PG reading the info from the file
      
      TFile * f_dum = new TFile (filenames.at (i).c_str ()) ;
      inputfiles.push_back (f_dum) ;
      TH1F * h_dum = (TH1F *) f_dum->Get ("h_MWW_mg") ;
      h_dum->Sumw2 () ;
      for (int j = 1 ; j <= h_dum->GetNbinsX () ; ++j) if (h_dum->GetBinContent (j) - h_dum->GetBinError (j) < 0) h_dum->SetBinError (j, 0.) ;
      h_mg_signals.push_back (h_dum) ;
      TString name = "rdh_" ;
      name += masses.at (i) ;
      RooDataHist * rdh_dum = new RooDataHist (name.Data (), name.Data (), RooArgList (x), h_dum) ;
      rdh_mg_signals.push_back (rdh_dum) ;
      
      //PG fitting the thing
      
      RooPlot * fr_dum = x.frame () ;
      rdh_dum->plotOn (fr_dum, MarkerStyle (4), MarkerColor (kGreen + 2)) ;

      mean_brwi.setVal (masses.at (i)) ;
      sigma_brwi.setVal (h_dum->GetRMS ()) ;
      mean_gaus.setVal (0) ;
      sigma_gaus.setVal (20.) ;
      
      model.fitTo (*rdh_dum, SumW2Error (kTRUE), PrintLevel (-1)) ;
      model.plotOn (fr_dum, LineColor (kBlack), LineWidth (1)) ;
//      pdf_brwi.fitTo (*rdh_dum, SumW2Error (kTRUE), PrintLevel (-1)) ;
//      pdf_brwi.plotOn (fr_dum, LineColor (kBlue), LineWidth (1)) ;

      mean_cb.setVal (masses.at (i)) ;
      sigma_cb.setVal (h_dum->GetRMS ()) ;
      alpha_cb.setVal (1) ;
      power_cb.setVal (3) ;
      pdf_cb.fitTo (*rdh_dum, SumW2Error (kTRUE), PrintLevel (-1)) ;
      pdf_cb.plotOn (fr_dum, LineColor (kRed), LineWidth (1)) ;

      mean_gaus.setVal (masses.at (i)) ;
      sigma_gaus.setVal (h_dum->GetRMS ()) ;
      pdf_gaus.fitTo (*rdh_dum, SumW2Error (kTRUE), PrintLevel (-1)) ;
      pdf_gaus.plotOn (fr_dum, LineColor (kGreen), LineWidth (1)) ;
      
      TCanvas * c_dum = new TCanvas () ;
      c_dum->SetLogy () ;
      fr_dum->Draw () ;
      name = "roofit_signal_" ;
      name += masses.at (i) ;
      name += ".pdf" ;
      c_dum->Print (name, "pdf") ;
    
    }
  
  return 0 ;

} 

/*




  RooRealVar frac ("frac", "frac", 0.25, 0., 1.) ; 

  RooAddPdf tot ("tot","frac * gauss_1 + gauss_2",  
                 RooArgList (gauss_1,gauss_2) , 
                 RooArgList (frac) ,
                 kTRUE) ;

  RooDataSet * data = tot.generate (RooArgSet (x), 10000) ;

  // plot the distributions
  RooPlot *xplot = x.frame () ;
  data->plotOn (xplot) ;
  gauss_1.plotOn (xplot, LineColor (kRed)) ;
  gauss_2.plotOn (xplot, LineColor (kBlue)) ;
  tot.plotOn (xplot, LineColor (kGreen + 2)) ;
  TCanvas c1 ;
  xplot->Draw () ;
  c1.Print ("DF_distr.gif","gif") ;
  
  // first gaussian definition
  RooRealVar mean_gen ("mean_gen", "mean of gaussian _gen", 3., -10., 10.) ; 
  RooRealVar sigma_gen ("sigma_gen", "width of gaussian _gen", 0.5, 0., 5.) ; 
  RooGaussian gauss_gen ("gauss_gen", "gaussian PDF", x, mean_gen, sigma_gen) ; 

  gauss_gen.fitTo (*data, RooFit::Minos (1)) ;

  x.setRange ("left", mean_gen.getVal () - 2 * sigma_gen.getVal (), mean_gen.getVal ()) ;

  RooRealVar mean_3 ("mean_3", "mean of gaussian _3", 3., -10., 10.) ; 
  RooRealVar sigma_3 ("sigma_3", "width of gaussian _3", 0.5, 0., 5.) ; 
  RooGaussian gauss_3 ("gauss_3", "gaussian PDF", x, mean_3, sigma_3) ; 
  gauss_3.fitTo (*data, RooFit::Minos (1), Range ("left")) ;

  x.setRange ("right", mean_gen.getVal () + 0.5 * sigma_gen.getVal (), mean_gen.getVal () + 3 * sigma_gen.getVal ()) ;

  RooRealVar mean_4 ("mean_4", "mean of gaussian _4", 3., -10., 10.) ; 
  RooRealVar sigma_4 ("sigma_4", "width of gaussian _4", 0.5, 0., 5.) ; 
  RooGaussian gauss_4 ("gauss_4", "gaussian PDF", x, mean_4, sigma_4) ; 
  gauss_4.fitTo (*data, RooFit::Minos (1), Range ("right")) ;

  RooPlot *xplot_1 = x.frame () ;
  TLegend leg_bkg (0.5, 0.6, 0.9, 0.9, NULL, "brNDC") ;
  leg_bkg.SetBorderSize (0) ;
  leg_bkg.SetTextFont (42) ;
  leg_bkg.SetTextSize (0.04) ;
  leg_bkg.SetLineColor (1) ;
  leg_bkg.SetLineStyle (1) ;
  leg_bkg.SetLineWidth (1) ;
  leg_bkg.SetFillColor (0) ;
  leg_bkg.SetFillStyle (0) ;

  data->plotOn (xplot_1) ;
  leg_bkg.AddEntry (xplot_1->FindObject (xplot_1->nameOf (xplot_1->numItems () - 1)), "data", "lp") ;
  gauss_gen.plotOn (xplot_1, LineColor (kGreen + 2)) ;
  leg_bkg.AddEntry (xplot_1->FindObject (xplot_1->nameOf (xplot_1->numItems () - 1)), "gaussGen", "lp") ;
  gauss_3.plotOn (xplot_1, LineColor (kGreen)) ;
  leg_bkg.AddEntry (xplot_1->FindObject (xplot_1->nameOf (xplot_1->numItems () - 1)), "gauss 3", "lp") ;
  gauss_4.plotOn (xplot_1, LineColor (kBlue)) ;
  leg_bkg.AddEntry (xplot_1->FindObject (xplot_1->nameOf (xplot_1->numItems () - 1)), "gauss 4", "lp") ;
  xplot_1->Draw () ;
//  leg_bkg.Draw () ;
 
  for (int i = 0 ; i < xplot_1->numItems () ; ++i) 
    {
      TString obj_name = xplot_1->nameOf (i) ; 
//      if (obj_name=="") continue ;
      cout << Form ("%d. '%s'\n",i,obj_name.Data ()) ;
    }
  
  c1.Print ("DF_fit.gif","gif") ;

  return 0 ;
}

*/