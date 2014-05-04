//compile with: 
//g++ -o prova.o plot_all_parameters.cpp -c `root-config --cflags`;
//g++ -o prova.exe prova.o `root-config --libs`;

#include "Functions.cc"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TLine.h"
#include "TSystem.h"
#include "TStyle.h"

using namespace std;

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double max (double uno, double due)
{
  if (uno > due) return uno ;
  return due ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double min (double uno, double due)
{
  if (uno < due) return uno ;
  return due ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


/*** single crystall ball ***/
/*double crystalBallLowHigh (double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2

  double xx = x[0];
  double mean   = par[1] ; // mean
  double sigmaP = par[2] ; // sigma of the positive side of the gaussian
  double sigmaN = par[3] ; // sigma of the negative side of the gaussian
  double alpha  = par[4] ; // junction point on the positive side of the gaussian
  double n      = par[5] ; // power of the power law on the positive side of the gaussian

  if ((xx-mean)/sigmaP > fabs(alpha))
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    return par[0] * A * pow(B + (xx-mean)/sigmaP, -1.*n);
  }
    
  else if ((xx-mean) > 0)
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaP*sigmaP) );
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaN*sigmaN) );
  }
  
}

*/


/*** double crystall ball ***/
double doubleGausCrystalBallLowHigh (double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2

  double xx = x[0];
  double mean   = par[1] ; // mean
  double sigmaP = par[2] ; // sigma of the positive side of the gaussian
  double sigmaN = par[3] ; // sigma of the negative side of the gaussian
  double alpha  = par[4] ; // junction point on the positive side of the gaussian
  double n      = par[5] ; // power of the power law on the positive side of the gaussian
  double alpha2 = par[6] ; // junction point on the negative side of the gaussian
  double n2     = par[7] ; // power of the power law on the negative side of the gaussian

  if ((xx-mean)/sigmaP > fabs(alpha))
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    return par[0] * A * pow(B + (xx-mean)/sigmaP, -1.*n);
  }
  
  else if ((xx-mean)/sigmaN < -1.*fabs(alpha2))
  {
    double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
    double B = n2/fabs(alpha2) - fabs(alpha2);
    
    return par[0] * A * pow(B - (xx-mean)/sigmaN, -1.*n2);
  }
  
  else if ((xx-mean) > 0)
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaP*sigmaP) );
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaN*sigmaN) );
  }
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void printArray (double * array, int N)
{
  std::cout << "double * params[" << N << "] = {" ;
  for (int i = 0 ; i < N - 1 ; ++i)
    std::cout << array[i] << ", " ;
  std::cout << array[N-1] << "} ;\n" ;
  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void setParNamesdoubleGausCrystalBallLowHigh (TF1 * func)
{
  func->SetParName (1, "mean") ;
  func->SetParName (2, "sigma") ;
  func->SetParName (3, "alphaR") ;
  func->SetParName (4, "nR") ;
  func->SetParName (5, "alphaL") ;
  func->SetParName (6, "nL") ;
  return ;
}  


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t crystalBallLowHighRatio (Double_t * xx, Double_t * par) // (SBI - B) / S
{
  return crystalBallLowHigh (xx, par) / crystalBallLowHigh (xx, par + 7) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t relativeCrystalBallLowHighRatio (Double_t * xx, Double_t * par) // [(SBI - B) - S] / S
{
  return (crystalBallLowHigh (xx, par) - crystalBallLowHigh (xx, par + 7)) / crystalBallLowHigh (xx, par + 7) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t withLeftAddOn (Double_t * xx, Double_t * par)
{
   double centre       = par[0] ;
   double shift        = par[1] ;
   double scale        = par[2] ;
   double slope        = par[3] ;
   double secondCentre = par[4] ;
   double secondSlope  = par[4] ;
   double x            = xx[0] - centre ;
   
   if (x > 0) 
   return scale * x * TMath::Exp (-1 * slope * x) + shift +
          (x - secondCentre) * TMath::Exp (-1 * secondSlope * (x - secondCentre)) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t doubleSlope (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double slope_left  = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return scale * x * TMath::Exp (     slope_left  * x) + shift ;
//   else       return -1 * scale * x * TMath::Exp (-1 * slope_left  * x) + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----



Double_t doublePeakModel (Double_t * xx, Double_t * par)
{
  double scale    = par[0] ;
  double shift    = par[1] ;
  double distance = par[2] ;
  double gamma    = par[3] ;
  double x = xx[0] - shift ;

  double max = 1./3. * sqrt( -3. * gamma + 3. * distance * distance + 6. * sqrt( gamma * gamma + distance * distance * gamma + pow( distance, 4.)));
  double height = fabs ( 1. / ((max - distance) * (max - distance) + gamma) - 1 / ((max + distance) * (max + distance) + gamma)) ;
  double norm = 1. / height ;
  
  return scale * norm * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;

//  double norm = 1. / (shift * shift + gamma) - 1 / ((shift + 2 * distance) * (shift + 2 * distance) + gamma) ;
//  return scale * (1. / norm) * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;
//  return scale * norm * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;
//  return scale * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;

}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t singleSlope (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double slope_left  = par[3] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return scale * x * TMath::Exp (     slope_left  * x) + shift ;
//   else       return -1 * scale * x * TMath::Exp (-1 * slope_left  * x) + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t parabolicAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double power_left  = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return power_left * x * x + scale * x + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t cubicAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double pwr2_left   = par[4] ;
   double pwr3_left   = par[5] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return pwr3_left * x * x * x + pwr2_left * x * x + scale * x + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t sinusAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double freq_left   = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return scale * TMath::Sin (freq_left * x) / freq_left + shift ; //pg this guarantees C1 properties
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t logAndExp (Double_t * xx, Double_t * par)
{
   double centre      = par[0] ;
   double shift       = par[1] ;
   double scale       = par[2] ;
   double slope_right = par[3] ;
   double log_fact    = par[4] ;
   double x           = xx[0] - centre ;
   
   if (x > 0) return scale * x * TMath::Exp (-1 * slope_right * x) + shift ;
   else       return scale * x * TMath::Log (-1 * log_fact * x) / log_fact + shift ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int findBin (TH1F * h, double val)
{
  int i = 1 ;
  for ( ; i <= h->GetNbinsX () ; ++i)
    {
      if (h->GetBinLowEdge (i) > val) break ;
    } 
  return i - 1 ;  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double normaliseToBinWidth (TH1F * h)
{
  double scale = h->GetNbinsX () / (h->GetXaxis ()->GetXmax () - h->GetXaxis ()->GetXmin ()) ;
  h->Scale (scale) ;
  return scale ;  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


//PG parameters of the signal-and-interference first,
//PG parameters of the signal-only later
double ratio_crystalBallLowHigh (double* x, double* par)
{
  double den = crystalBallLowHigh (x, par + 7) ; // signal only
  if (den == 0) return -1. ;
  double num = crystalBallLowHigh (x, par) ;    // signal and interference
  return num / den ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


//PG parameters of the signal-and-interference first,
//PG parameters of the signal-only later
double ratio_crystalBallLowHighWithRise (double* x, double* par)
{
  double den = crystalBallLowHigh (x, par + 9) ; // signal only
  if (den == 0) return -1. ;
  double num = crystalBallLowHighWithRise (x, par) ;    // signal and interference
  return num / den ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


//PG parameters of the signal-and-interference first,
//PG parameters of the signal-only later
double diff_crystalBallLowHighWithRise (double* x, double* par)
{
  return crystalBallLowHighWithRise (x, par) - crystalBallLowHigh (x, par + 9) ;
  //PG    signal and interference                          signal only
}




// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


//PG parameters of the signal-and-interference first,
//PG parameters of the signal-only later
double diff_crystalBallLowHigh (double* x, double* par)
{
  return crystalBallLowHigh (x, par) - crystalBallLowHigh (x, par + 7) ;
  //PG    signal and interference          signal only
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


TF1 * FIT_madgraph_signal (TH1F * h_MWW_mg, double mass, double rangeScale, TString suffix, bool useLikelihood = false)
{  
  TCanvas * c4_mg = new TCanvas ("c4_mg", "c4_mg") ;

  //PG first fit: get the width for the second fit
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * gauss_mg = new TF1 ("gauss_mg", "gaus", 0, 2000) ;
  gauss_mg->SetNpx (10000) ;
  gauss_mg->SetLineWidth (1) ;
  gauss_mg->SetLineColor (kGray + 2) ;
  gauss_mg->SetParameter(1,mass);
  h_MWW_mg->Fit ("gauss_mg", "+", "", mass - 1. * h_MWW_mg->GetRMS () , mass + 1. * h_MWW_mg->GetRMS ()) ;


  TH1F * h_MWW_mg_error = (TH1F*)h_MWW_mg->Clone ("h_MWW_mg_error") ;
  (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (h_MWW_mg_error, 0.68) ;
  h_MWW_mg_error->SetMarkerSize (0) ;
  h_MWW_mg_error->SetFillColor (kAzure - 9) ;
  h_MWW_mg_error->SetFillStyle (3003) ;

  
  //PG second fit: first with chisq, if requested with likelihood also
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * func_mg_1 = new TF1 ("func_mg_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_mg_1->SetNpx (10000) ;
  func_mg_1->SetLineWidth (1) ;
  func_mg_1->SetLineColor (kBlue + 1) ;

  setParNamesdoubleGausCrystalBallLowHigh (func_mg_1) ;

  func_mg_1->SetParameter (0, h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ())) ;  // multiplicative scale
  func_mg_1->FixParameter (1, gauss_mg->GetParameter(1)) ;                                                  // mean
  func_mg_1->SetParameter (2, gauss_mg->GetParameter (2)) ;                            // gaussian sigma
  func_mg_1->SetParLimits (2, 0.01 * h_MWW_mg->GetRMS (), 20 * h_MWW_mg->GetRMS ()) ;
  func_mg_1->SetParameter (3, 1.5) ;                                                   // right junction
  func_mg_1->SetParameter (4, 2) ;                                                     // right power law order
  func_mg_1->SetParameter (5, 0.95) ;                                                  // left junction
  func_mg_1->SetParameter (6, 2.38) ;                                                  // left power law order

  //PG set the range of the fit
  int sign = 1 ;
  if (mass < 400) sign = -2 ;
  if (mass > 700) sign = -1.5 ;
  if (mass > 850) sign = -1. ;
  std::cout << "-------------------\nFITTING THE MADGRAPH SIGNAL\n\n-------------------\n" ;
  TFitResultPtr result = h_MWW_mg->Fit ("func_mg_1", "+L", "", 100, 2 * mass) ;
  //  TFitResultPtr result = h_MWW_mg->Fit ("func_mg_1", "+L", "", 0.5 * mass + sign * 50, 2 * mass) ;
  if (useLikelihood && mass < 800)
    {
      std::cout << "-------------------\nFITTING THE MADGRAPH SIGNAL W/ LIKELIHOOD\n\n-------------------\n" ;
      func_mg_1->SetParameters (func_mg_1->GetParameters ()) ; //PG not sure this is necessary
      func_mg_1->SetLineColor (kBlue + 3) ;
      h_MWW_mg->Fit ("func_mg_1", "+L", "", 0.5 * mass - 50, 2 * mass) ;
    }


  //PG plotting the result
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  double ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  double ymin = h_MWW_mg->GetBinContent (h_MWW_mg->GetMinimumBin ()) ;
  if (ymin <= 0) ymin = ymax / 10000. ; 
  TH1F * c4_mg_frame = (TH1F *) c4_mg->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_mg->SetLogy () ;
  c4_mg_frame->SetTitle (0) ;
  c4_mg_frame->SetStats (0) ;
  c4_mg_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_mg_error->Draw ("e3same") ;
  h_MWW_mg->Draw ("EPsame") ;

  //PG draw the junction points in the fit image
  
  double rightTh = fabs (func_mg_1->GetParameter (3)) * func_mg_1->GetParameter (2) + func_mg_1->GetParameter (1) ;
  double leftTh  = -1 * fabs (func_mg_1->GetParameter (5)) * func_mg_1->GetParameter (2) + func_mg_1->GetParameter (1) ;

  TLine * l_rightTh = new TLine (rightTh, 0.9 * ymin, rightTh, 1.1 * ymax) ;
  l_rightTh->SetLineColor (kRed) ;
  l_rightTh->Draw ("same") ;
  TLine * l_leftTh = new TLine (leftTh, 0.9 * ymin, leftTh, 1.1 * ymax) ;
  l_leftTh->SetLineColor (kRed) ;
  l_leftTh->Draw ("same") ;

  c4_mg->Update () ;
  c4_mg->Print (TString ("signals_mg_log") + suffix, "pdf") ;

  c4_mg->SetLogy (0) ;
  c4_mg->Print (TString ("signals_mg_lin") + suffix, "pdf") ;

  return func_mg_1 ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//


//PG no left rise
TF1 * FIT_phantom_signal (TH1F * diff, double mass, double rangeScale, TString suffix,  double cprime, bool useLikelihood = false)
{  
 TCanvas * c4_ph = new TCanvas ("c4_ph", "c4_ph") ;

 //   int cprime100 = cprime*100;   //cprime has some round-off problems not-well understood: if (cprime==x) does not work, so it is necessary to convert it to int

 // std::cout<<cprime*100<<" BBBB"<<std::endl;
 //  if (cprime == 1.)   
 //       std::cout<<"AAAAA"<<std::endl;
  
 
  //PG first fit: get the width for the second fit
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  double fact = 1. ; 
 TF1 * gauss_ph = new TF1 ("gauss_ph", "gaus", mass-fact*diff->GetRMS(),mass+fact*diff->GetRMS());// 0, 2000) ;
  gauss_ph->SetNpx (10000) ;
  gauss_ph->SetLineWidth (1) ;
  gauss_ph->SetLineColor (kGray + 2) ;

  gauss_ph->SetParameter(0,0.0001);
  gauss_ph->SetParameter(1,650);
    gauss_ph->SetParLimits(1,610,690);
  gauss_ph->SetParameter(2,30);
  gauss_ph->SetParLimits(2,0,150);
  //   diff->Fit ("gauss_ph", "+R");//, "", mass - fact * diff->GetRMS () , mass + fact * diff->GetRMS ()) ;

  //PG second fit: first with chisq, if requested with likelihood also
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * func_ph_1;
  if (mass==350)
        func_ph_1 = new TF1 ("func_ph_1", doubleGausCrystalBallLowHighPlusExp, 250, 450, 9) ;
  else
        func_ph_1 = new TF1 ("func_ph_1", doubleGausCrystalBallLowHighPlusExp, 200, 2000, 9) ;
  //      TF1 * func_ph_1 = new TF1 ("func_ph_1", crystalBallLowHigh, 0, 2000, 6) ;
  func_ph_1->SetNpx (10000) ;
  func_ph_1->SetLineWidth (1) ;
  func_ph_1->SetLineColor (kRed + 1) ;
  
  setParNamesdoubleGausCrystalBallLowHigh (func_ph_1) ;

  func_ph_1->SetParameter (0, diff->GetBinContent (diff->GetMaximumBin ())) ;  // multiplicative scale                                                

  if (mass==350) {
    func_ph_1->FixParameter(1,mass*(1.000-0.001*cprime));

    func_ph_1->FixParameter (2, 15.2*(0.013+0.386*cprime)) ;
  }

  else if (mass==650) {
  func_ph_1->FixParameter(1,mass*(0.997-0.010*cprime));
  func_ph_1->FixParameter(2,158*(0.023+0.361*cprime));
  }

  else if (mass==800) {

    if (cprime == 1.) { 
          func_ph_1->FixParameter (1, mass*(1.007-0.122*cprime+0.093*cprime*cprime));//gauss_ph->GetParameter(1)) ;          // mean 
	  func_ph_1->SetParameter(2, 304*(0.003+0.405*cprime) );
	  func_ph_1->SetParLimits(2, 50,300 );
    }
    else     {
          func_ph_1->FixParameter (1, mass*(1.007-0.122*cprime+0.093*cprime*cprime));//gauss_ph->GetParameter(1)) ;          // mean 
          func_ph_1->FixParameter (2, 304*(0.003+0.405*cprime) );//gauss_ph->GetParameter (2)) ; 
    }
  }

  else if (mass==1000) {
    func_ph_1->FixParameter (1, mass*(1.020-0.322*cprime+0.210*cprime*cprime)) ;                                                  // mean            
    if (cprime!=0.9)
      func_ph_1->FixParameter (2,647*0.38*cprime);
    else
      func_ph_1->FixParameter (2,156);
  }


  if (mass==1000) {

    if (cprime==1.) {
    func_ph_1->SetParameter (3, 5) ;                                                   // right junction 
    func_ph_1->SetParLimits(3,3,7);  
    func_ph_1->SetParameter (4, 4) ;   // right power law order                                                
    func_ph_1->SetParLimits (4,3.5,5);
    func_ph_1->SetParameter (5, 0.95) ;                                                  // left junction
    func_ph_1->SetParLimits (5, 0.1,3) ;                                                  // left junction
    func_ph_1->SetParameter (6, 1.5) ;                                                  // left power law order   
    func_ph_1->SetParLimits (6, 0.1,4) ;                                                  // left power law order   
    func_ph_1->SetParameter (7, 1.2);
    func_ph_1->SetParLimits (7, 0.5,1.5) ;                              
    func_ph_1->SetParameter (8,16);
    func_ph_1->SetParLimits (8,10,20) ;
    }
    else {
    func_ph_1->SetParameter (3, 5) ;                                                   // right junction 
    func_ph_1->SetParLimits(3,0.1,8);  
    func_ph_1->SetParameter (4, 4) ;   // right power law order                                                
    func_ph_1->SetParLimits (4,3.5,7);
    func_ph_1->SetParameter (5, 0.95) ;                                                  // left junction
    func_ph_1->SetParLimits (5, 0.1,3) ;                                                  // left junction
    func_ph_1->SetParameter (6, 1.5) ;                                                  // left power law order   
    func_ph_1->SetParLimits (6, 0.1,10) ;                                                  // left power law order   
    func_ph_1->SetParameter (7, 1.2);
    func_ph_1->SetParLimits (7, 0.5,1.5) ;                              
    func_ph_1->SetParameter (8,16);
    func_ph_1->SetParLimits (8,10,20) ;
    }
}

  else if (mass==800) {
    func_ph_1->SetParLimits (3, 0.1,2.5) ;                                                   // right junction  

    if (cprime==0.1)  func_ph_1->SetParLimits(4,0.6,2);
    else if (cprime==0.3)  func_ph_1->SetParLimits(4,0.6,2);
    else if (cprime==0.5) {  func_ph_1->SetParLimits(4,1,2);      func_ph_1->SetParLimits (3, 0.1,10) ;  }
    else if (cprime==0.7) {  func_ph_1->SetParLimits (3, 0.1,2.5) ;  func_ph_1->SetParLimits(4,0.5,1.5); }
    else if (cprime==0.9) {  func_ph_1->SetParLimits (3, 0.1,2.5) ;  func_ph_1->SetParLimits(4,0.5,2.5); }
    else         {  func_ph_1->SetParLimits (3, 1,2.5);      func_ph_1->SetParLimits (4,0.1,3.5); }

    if (cprime!=1.) {
    func_ph_1->SetParameter (5, 0.62) ;                                                  // left junction
    func_ph_1->SetParLimits (5, 0.1,2) ;                                                  // left junction
    func_ph_1->SetParameter (6, 1.5) ;                                                  // left power law order   
    func_ph_1->SetParLimits (6, 0.1,3) ;                                                  // left power law order   
    }

    if (cprime==0.1) { func_ph_1->SetParameter(5,1.3);  func_ph_1->SetParLimits(5,1,2); }
   if (cprime==0.5) { func_ph_1->SetParameter(5,0.9);  func_ph_1->SetParLimits(5,0.5,2); }
      if (cprime==1.) {  
    func_ph_1->SetParameter (5, 1.2) ;                                                  // left junction
    func_ph_1->SetParLimits (5, 0.5,2) ;                                                  // left junction
    func_ph_1->SetParLimits(6,1,3.5); }

    func_ph_1->SetParameter (7, 4);
    func_ph_1->SetParLimits (7, 0.1,5) ;                              
    func_ph_1->SetParameter (8,100);
    if (cprime==0.5 || cprime==0.3)   func_ph_1->SetParLimits (8,5,15);
    else               func_ph_1->SetParLimits (8, 5,500) ;

    if (cprime==1.)  {
      //      func_ph_1->SetParLimits(5,0.1,3);   
      //   func_ph_1->SetParLimits (6, 0.1,6) ;                                                  // left power law order   
      func_ph_1->SetParameter (7, 4);
      func_ph_1->SetParLimits (7, 0.1,5) ;                              
      func_ph_1->SetParLimits (8,5,30);
    }

}

  else {

  func_ph_1->SetParLimits (3, 0.1, 20) ;                                                   // right junction 
  if (mass==650 && cprime==0.5)    func_ph_1->SetParLimits (3, 2, 5) ; 
  if (mass==650 && cprime==0.7)    func_ph_1->SetParLimits (3, 2, 3) ; 
   func_ph_1->SetParameter (4, 2.1) ;   // right power law order                                                
  func_ph_1->SetParLimits (4, 0.1, 2.5) ;   // right power law order                                                
  if (mass==350)  func_ph_1->SetParLimits(4,1,5);
  if (mass==350 && cprime==0.1)  func_ph_1->SetParLimits(4,1,3);
  if (mass==350 && cprime==0.9)  func_ph_1->SetParLimits(4,4.5,5.5);
   if (mass==350 && cprime==1.)  func_ph_1->SetParLimits(4,0.1,5);
   if (mass==650 && cprime==0.3)  func_ph_1->SetParLimits(4,0.1,1.1);
   if (mass==650 && cprime==0.5)  func_ph_1->SetParLimits(4,0.1,2.3);


      if (mass==650 && cprime==0.7)  func_ph_1->SetParLimits(4,1.5,2.5);
     if (mass==650 && cprime==0.9)  func_ph_1->SetParLimits(4,1.,3.);
     if (mass==650 && cprime==1.0)  func_ph_1->SetParLimits(4,1.5,3.5);
  // if (mass==650 && cprime==0.9)   { func_ph_1->SetParLimits (4, 0.1, 1.5) ;  func_ph_1->SetParLimits (4, 0.1, 1.5) ; } 

  func_ph_1->SetParameter (5, 0.95) ;                                                  // left junction
  func_ph_1->SetParLimits (5, 0.5,4) ;                                                  // left junction
  func_ph_1->SetParameter (6, 1.5) ;                                                  // left power law order   
  func_ph_1->SetParLimits (6, 0.1,4) ;                                                  // left power law order   
  func_ph_1->SetParameter (7, 4);
  func_ph_1->SetParLimits (7, 0.01,5) ;                              
  func_ph_1->SetParameter (8,53);
  func_ph_1->SetParLimits (8, 10,200) ;                              

  }
    


     std::cout << "-------------------\nFITTING THE PHANTOM SIGNAL - mass: " << mass << "    cprime: " << cprime <<"\n\n-------------------\n" ;
  if (mass==350)
    diff->Fit ("func_ph_1", "+", "", 250, 450) ;
  else
    diff->Fit ("func_ph_1", "+", "", 200, 2 * mass) ;
  if (useLikelihood)
    {
      std::cout << "-------------------\nFITTING THE PHANTOM SIGNAL W/ LIKELIHOOD\n\n-------------------\n" ;
      func_ph_1->SetParameters (func_ph_1->GetParameters ()) ;
      func_ph_1->SetLineColor (kRed + 3) ;
      diff->Fit ("func_ph_1", "+L", "", 0.5 * mass - 50, 1.5 * mass) ;
    }
    
  double ymax = diff->GetBinContent (diff->GetMaximumBin ()) ;
  ymax = max (ymax, func_ph_1->GetMaximum ()) ;
  double ymin = diff->GetBinContent (diff->GetMinimumBin ()) ;
  ymin = max (ymin, -0.1 * ymax) ;
  if (ymin <= 0) ymin = ymax / 10000. ; 
  TH1F * c4_ph_frame;
  if (mass==350)
    c4_ph_frame = (TH1F *) c4_ph->DrawFrame (320, 0.9 * ymin, 380, 1.1 * ymax);//, 1.1 * ymax) ;
  else
    c4_ph_frame = (TH1F *) c4_ph->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax);//, 1.1 * ymax) ;
  c4_ph_frame->SetTitle (0) ;
  c4_ph_frame->SetStats (0) ;
  c4_ph_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  //    gauss_ph->Draw ("same") ;
  func_ph_1->Draw("same");
  diff->Draw ("EPsame") ;



  TH1F * h_MWW_ph_error = (TH1F*)diff->Clone ("h_MWW_ph_error") ;
  (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (h_MWW_ph_error, 0.68) ;
  h_MWW_ph_error->SetMarkerSize (0) ;
  h_MWW_ph_error->SetFillColor (kBlue) ;
  h_MWW_ph_error->SetFillStyle (3002) ;

    
  double rightTh_ph = fabs (func_ph_1->GetParameter (3)) * fabs (func_ph_1->GetParameter (2)) + func_ph_1->GetParameter (1) ;
  double leftTh_ph  = -1 * fabs (func_ph_1->GetParameter (5)) * fabs (func_ph_1->GetParameter (2)) + func_ph_1->GetParameter (1) ;

  TLine * l_rightTh_ph = new TLine (rightTh_ph, 0.9 * ymin, rightTh_ph, 1.1 * ymax) ;
  l_rightTh_ph->SetLineColor (kRed) ;
  l_rightTh_ph->Draw ("same") ;
  TLine * l_leftTh_ph = new TLine (leftTh_ph, 0.9 * ymin, leftTh_ph, 1.1 * ymax) ;
  l_leftTh_ph->SetLineColor (kRed) ;
  l_leftTh_ph->Draw ("same") ;
  
  h_MWW_ph_error->Draw("e3same");

  c4_ph->Print (TString ("signals_ph_lin") + suffix, "png") ;
  //  c4_ph->Print (TString ("signals_ph_lin") + suffix, "png") ;
  c4_ph->SetLogy () ;
  c4_ph->Print (TString ("signals_ph_log") + suffix, "png") ;
  // c4_ph->Print (TString ("signals_ph_log") + suffix, "png") ;

  return func_ph_1 ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double leftRise (double * x, double * par)
{
//  return (-1. / (par[0] * (x[0] - 180.)) + 1) ;
//  return (-1. / (par[0] * TMath::Sqrt (x[0] - 180.))) ;
//  return (-1. / (    par[0] * (x[0] - 180.)    )) ;
  return par[0] * ( 1. / (1. + TMath::Exp (-1. * par[1] * x[0]) ) - 1. )  ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double crystalBallLowHighWithRise (double* x, double* par)
{
  //PG the first 7 parameters are passed to the crystalBallLowHigh
//  return TMath::Log (par[7] * (x[0] - 180.)) * crystalBallLowHigh (x, par) ;
//  return leftRise (x, par + 7) * crystalBallLowHigh (x, par) ;
  return leftRise (x, par + 7) + crystalBallLowHigh (x, par) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


//PG left rise
TF1 * FIT_phantom_signal_2 (TH1F * diff, double mass, double rangeScale, TString suffix, bool useLikelihood = false)
{  
 TCanvas * c4_ph = new TCanvas ("c4_ph", "c4_ph") ;

  //PG first fit: get the width for the second fit
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * gauss_ph = new TF1 ("gauss_ph", "gaus", 0, 2000) ;
  gauss_ph->SetNpx (10000) ;
  gauss_ph->SetLineWidth (1) ;
  gauss_ph->SetLineColor (kGray + 2) ;
  gauss_ph->SetParameter (1, mass) ;
  gauss_ph->SetParameter (2, 0.5 * diff->GetRMS ()) ;
  double fact = 1. ;
  double span = min (fact * fabs (diff->GetRMS ()), 0.5 * mass) ;
  diff->Fit ("gauss_ph", "+", "", mass - span , mass + span) ;

  //PG preliminary fit of the left rise
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  TF1 * f_leftRise = new TF1 ("f_leftRise", leftRise, 181, 2000, 2) ;
  f_leftRise->SetNpx (10000) ;
  f_leftRise->SetLineWidth (1) ;
  f_leftRise->SetLineColor (kOrange + 2) ;
  f_leftRise->SetParLimits (0, 0., 10000.) ;
  double N = (1.5 - 4) / (800. - 500.) * (mass - 500.) + 4 ;
  double riseLimit = mass - N * gauss_ph->GetParameter (2) ;
  diff->Fit ("f_leftRise", "+", "", 180., riseLimit) ;

  //PG second fit: first with chisq, if requested with likelihood also
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * func_ph_1 = new TF1 ("func_ph_1", crystalBallLowHighWithRise, 0, 2000, 9) ;
  func_ph_1->SetNpx (10000) ;
  func_ph_1->SetLineWidth (1) ;
  func_ph_1->SetLineColor (kRed + 1) ;
  
  setParNamesdoubleGausCrystalBallLowHigh (func_ph_1) ;

  func_ph_1->SetParameter (0, diff->Integral ()) ;                //PG multiplicative scale
  func_ph_1->SetParameter (1, mass) ;                             //PG mean
  func_ph_1->SetParameter (2, 2 * gauss_ph->GetParameter (2)) ;   //PG gaussian sigma
  func_ph_1->SetParLimits (2, 0., mass) ;                         
  func_ph_1->SetParameter (3, 1) ;                                //PG right junction
  func_ph_1->SetParLimits (3, 0.1, 5) ;                           //PG right junction
  func_ph_1->FixParameter (4, 3) ;                                //PG right power law order            //PG NB THIS IS FIXED
  func_ph_1->FixParameter (5, 1) ;                                //PG left junction                    //PG NB THIS IS FIXED
  func_ph_1->FixParameter (6, 3) ;                                //PG left power law order             //PG NB THIS IS FIXED
  func_ph_1->FixParameter (7, f_leftRise->GetParameter (0)) ;     //PG from the first fit to the rise
  func_ph_1->FixParameter (8, f_leftRise->GetParameter (1)) ;     //PG from the first fit to the rise

  std::cout << "-------------------\nFITTING THE PHANTOM SIGNAL with left turn on \n\n-------------------\n" ;
//  diff->Fit ("func_ph_1", "", "", 0.5 * mass - 50, 2 * mass) ;
  diff->Fit ("func_ph_1", "", "", 200., 2 * mass) ;
  if (useLikelihood)
    {
      std::cout << "-------------------\nFITTING THE PHANTOM SIGNAL W/ LIKELIHOOD with left turn on\n\n-------------------\n" ;
      func_ph_1->SetParameters (func_ph_1->GetParameters ()) ;
      func_ph_1->SetLineColor (kRed + 3) ;
//      diff->Fit ("func_ph_1", "L+", "", 0.5 * mass - 50, 1.5 * mass) ;
      diff->Fit ("func_ph_1", "L+", "", 200., 1.5 * mass) ;
    }
    
  double ymax = diff->GetBinContent (diff->GetMaximumBin ()) ;
  double ymin = diff->GetBinContent (diff->GetMinimumBin ()) ;
  TH1F * c4_ph_frame = (TH1F *) c4_ph->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_ph_frame->SetTitle (0) ;
  c4_ph_frame->SetStats (0) ;
  c4_ph_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  gauss_ph->Draw ("same") ;
  diff->Draw ("EPsame") ;

  //PG limits of the gaussian core
  double rightTh_ph = fabs (func_ph_1->GetParameter (3)) * fabs (func_ph_1->GetParameter (2)) + func_ph_1->GetParameter (1) ;
  double leftTh_ph  = -1 * fabs (func_ph_1->GetParameter (5)) * fabs (func_ph_1->GetParameter (2)) + func_ph_1->GetParameter (1) ;
  TLine * l_rightTh_ph = new TLine (rightTh_ph, 0.9 * ymin, rightTh_ph, 1.1 * ymax) ;
  l_rightTh_ph->SetLineColor (kRed) ;
  l_rightTh_ph->Draw ("same") ;
  TLine * l_leftTh_ph = new TLine (leftTh_ph, 0.9 * ymin, leftTh_ph, 1.1 * ymax) ;
  l_leftTh_ph->SetLineColor (kRed) ;
  l_leftTh_ph->Draw ("same") ;

  //PG limit of the prefit for the tail
  TLine * l_riseLimit = new TLine (riseLimit, 0.9 * ymin, riseLimit, 1.1 * ymax) ;
  l_riseLimit->SetLineColor (kOrange + 2) ;
  l_riseLimit->SetLineStyle (2) ;
  gauss_ph->Draw ("same") ;
  l_riseLimit->Draw ("same") ;
  f_leftRise->Draw ("same") ;

  c4_ph->Print (TString ("signals_ph_lin_leftTail") + suffix, "pdf") ;

  return func_ph_1 ;
}


// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====


int macro_findInterferece (string filename, double mass, std::vector<double> & param, std::vector<double> & param_error, int initialRebin=1, double cprime=1., bool useLeftRise = false)
{        

  //    TVirtualFitter::SetDefaultFitter ("Minuit2") ;
    // gSystem->Load ("/afs/cern.ch/user/l/lbrianza/work/PHD/INTERF/lvJJInterference_MINE/Functions.cc") ;
  gStyle->SetPadTopMargin (0.1) ;
  gSystem->AddIncludePath("-I$ROOFITSYS/include");

  TFile * f = new TFile (filename.c_str ()) ;
  TH1F * h_MWW_phbkgsig = (TH1F *) f->Get ("h_MWW_phbkgsig") ;
  TH1F * h_MWW_phbkg    = (TH1F *) f->Get ("h_MWW_phbkg") ;
  TH1F * h_MWW_mg       = (TH1F *) f->Get ("h_MWW_mg") ;


  int reBin = initialRebin ;
  if      (mass > 310 && mass < 610) reBin *= 2 ;
  else if (mass > 610 && mass < 790) reBin *= 4 ;
  else if (mass > 790 && mass < 860 && cprime==1.) reBin *= 20 ;
  else if (mass > 790 && mass < 860) reBin *= 8 ;
  else if (mass > 860)               reBin *= 20 ;
  double rangeScale = 1.5 ; // for plotting purposes
  if (mass > 480) rangeScale = 2 ;

  int scaling = 0 ;
  
  std::cout << "\n--------------------------------------------------\n" ;
  std::cout << "rebinning factor: " << reBin << "\n" ;
  std::cout << "MG bins number before rebinning: " << h_MWW_mg->GetNbinsX () << "\n" ;
  
  h_MWW_phbkgsig->Rebin (reBin) ;
  h_MWW_phbkg   ->Rebin (reBin) ;
  h_MWW_mg      ->Rebin (reBin) ;

  std::cout << "MG bins number after rebinning: " << h_MWW_mg->GetNbinsX () << "\n" ;

  //PG normalize to the bin width
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  normaliseToBinWidth (h_MWW_phbkgsig) ;
  normaliseToBinWidth (h_MWW_phbkg) ;
  normaliseToBinWidth (h_MWW_mg) ;

  double ymax= 2000. ;
  double ymin= 0. ;

  //PG graphics
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  h_MWW_phbkg->SetStats (0) ;
  h_MWW_phbkg->SetTitle ("") ;
  h_MWW_phbkg->SetLineColor (kOrange) ;
  h_MWW_phbkg->SetLineWidth (2) ;
  h_MWW_phbkg->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_phbkgsig->SetStats (0) ;
  h_MWW_phbkgsig->SetTitle ("") ;
  h_MWW_phbkgsig->SetLineColor (kRed) ;
  h_MWW_phbkgsig->SetLineWidth (2) ;
  h_MWW_phbkgsig->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_mg->SetStats (0) ;
  h_MWW_mg->SetTitle ("") ;
  h_MWW_mg->SetLineColor (kBlue) ;
  h_MWW_mg->SetLineWidth (2) ;
  h_MWW_mg->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;

  //PG histograms operations
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  //PG SBI - B
  TH1F * tempoDiff = (TH1F *) h_MWW_phbkgsig->Clone ("tempoDiff") ;
  tempoDiff->Add (h_MWW_phbkg, -1) ;
  tempoDiff->SetLineColor (kGreen + 1) ;

  if (scaling > 0)
    {
      double maxSI = tempoDiff->GetBinContent (tempoDiff->GetMaximumBin ()) ;
      double maxS = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
      h_MWW_mg->Scale (maxSI / maxS) ;
    }

  //PG SBI - B
  TH1F * diff = (TH1F *) h_MWW_phbkgsig->Clone ("diff") ;
  diff->Add (h_MWW_phbkg, -1) ;
  diff->SetLineColor (kGreen + 1) ;

  double N;
  for (int i=0; i<(diff->GetEntries()+1); i++)
    {
      N = diff->GetBinContent(i);
      //      for (int c=0; c<N; c++)
	//	std::cout<<diff->GetBinCenter(i)<<std::endl;
      //      std::cout<<N<<std::endl;
    }

  //PG (SBI - B) / S
  TH1F * ratio = (TH1F *) diff->Clone ("ratio") ;
  ratio->Divide (h_MWW_mg) ;

  //PG (SBI - B) - S
  TH1F * delta = (TH1F *) diff->Clone ("delta") ;
  delta->SetTitle ("") ;
  delta->Add (h_MWW_mg, -1) ;
  delta->SetLineColor (kViolet + 1) ;

  //PG ((SBI - B) - S) / S
  TH1F * relDiff = (TH1F*)delta->Clone ("relDiff") ;
  relDiff->Divide (h_MWW_mg) ;

  //PG plotting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  TString suffix (filename) ;
  suffix.ReplaceAll ("findInterference.", "") ;
  suffix.ReplaceAll ("root", "png") ;

  
  //PG initial spectra ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  /*
  TCanvas * c1 = new TCanvas ("c1", "c1") ;
  ymax = h_MWW_phbkgsig->GetBinContent (h_MWW_phbkgsig->GetMaximumBin ()) ;
  ymin = h_MWW_mg->GetBinContent (findBin (h_MWW_mg, mass + exp (mass * 0.0058461 + 0.65385))) ;
  TH1F * c1_frame = (TH1F *) c1->DrawFrame (200, 0.2 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c1_frame->SetTitle (0) ;
  c1_frame->SetStats (0) ;
  c1->SetLogy () ;
  c1_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_phbkg->Draw ("histsame") ;
  h_MWW_phbkgsig->Draw ("histsame") ;
  h_MWW_mg->Draw ("histsame") ;

  c1_leg = new TLegend (0.6,0.8,0.9,0.95) ;
  c1_leg->SetFillStyle (0) ;
  c1_leg->SetBorderSize (0) ;
  c1_leg->SetTextFont (42) ;
  c1_leg->AddEntry (h_MWW_phbkg, "B","l") ;
  c1_leg->AddEntry (h_MWW_phbkgsig, "SBI","l") ;
  c1_leg->AddEntry (h_MWW_mg, "S","l") ;
  c1_leg->Draw () ;
  
  c1->Print (TString ("spectra") + suffix, "pdf") ;
  */
  //PG S only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  //  TF1 * func_mg_1 = FIT_madgraph_signal (h_MWW_mg, mass, rangeScale, suffix, true) ;

  //PG (SBI - B) only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TF1 * func_ph_1 ;
  //  if (useLeftRise == true) func_ph_1 = FIT_phantom_signal_2 (diff, mass, rangeScale, suffix, false) ;
  func_ph_1 = FIT_phantom_signal   (diff, mass, rangeScale, suffix, cprime, false) ;

  for (int i=0; i<9; i++) {
    param.push_back(func_ph_1->GetParameter(i));
    param_error.push_back(func_ph_1->GetParError(i));
    //  std::std::cout<<param.at(i)<<" "<<param_error.at(i)<<std::endl;
  }



  return(0);
  //PG (SBI - B) - S only   i.e   INTERFERENCE FIT 
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  /*
  std::cout << "-------------------\nFITTING THE INTERFERENCE\n\n-------------------\n" ;

  TCanvas * c3 = new TCanvas ("c3", "c3") ;
  double ymax = delta->GetBinContent (delta->GetMaximumBin ()) ;
  double ymin = delta->GetBinContent (delta->GetMinimumBin ()) ;
  if (ymin > 0) ymin *= 0.5 ;
  else          ymin *= 1.5 ;
  if (mass > 800 && ymin < -3 * ymax) ymin = -2 * ymax ;
  TH1F * c3_frame = (TH1F *) c3->DrawFrame (200, ymin, rangeScale * mass, 2 * ymax) ;
  c3_frame->SetTitle (0) ;
  c3_frame->SetStats (0) ;
  c3_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;

  TF1 * f_doublePeakModel = new TF1 ("f_doublePeakModel", doublePeakModel, 0, 2000, 4) ;
  f_doublePeakModel->SetNpx (10000) ;
  f_doublePeakModel->SetLineWidth (1) ;
  f_doublePeakModel->SetLineColor (kRed + 1) ;

  f_doublePeakModel->SetParName (0, "scale") ;
  f_doublePeakModel->SetParName (1, "shift") ;
  f_doublePeakModel->SetParName (2, "distance") ;
  f_doublePeakModel->SetParName (3, "gamma") ; 

  f_doublePeakModel->SetParameter (0, -0.05 * h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ())) ;
  f_doublePeakModel->SetParameter (1, mass) ; 
  f_doublePeakModel->FixParameter (2, 0.0008) ; 
  f_doublePeakModel->SetParameter (3, mass * mass * 0.1 * 0.1) ;

  double aveWidth = 0.5 * sqrt (
      func_ph_1->GetParameter (2) * func_ph_1->GetParameter (2) +
      func_mg_1->GetParameter (2) * func_mg_1->GetParameter (2)  
    ) ;
  delta->Fit ("f_doublePeakModel", "+", "same", 0.5 * mass - 50, 2 * mass) ;

  delta->Draw ("histsame") ;
  c3_leg = new TLegend (0.5,0.8,0.9,0.95) ;
  c3_leg->SetFillStyle (0) ;
  c3_leg->SetBorderSize (0) ;
  c3_leg->SetTextFont (42) ;
  c3_leg->AddEntry (delta, "(SBI - B) - S","l") ;
  
  //PG plot the result of the fits on top
  TF1 * f_difference = new TF1 ("f_difference", diff_crystalBallLowHighWithRise, 0, 2000, 16) ;
  Double_t params_difference[16] ;
  func_ph_1->GetParameters (params_difference) ;
  func_mg_1->GetParameters (params_difference + 9) ;
  f_difference->SetParameters (params_difference) ;
  f_difference->SetNpx (10000) ;
  f_difference->SetLineWidth (2) ;
  f_difference->SetLineColor (kGray + 1) ;
  f_difference->Draw ("same") ;

  c3_leg->Draw () ;

  c3->Print (TString ("interf") + suffix, "pdf") ;
  */
  //PG (SBI - B) / S only   i.e   CORRECTION FACTOR
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  // TCanvas * c5 = new TCanvas ("c5", "c5") ;
  /*  ymax = ratio->GetBinContent (ratio->GetMaximumBin ()) ;
  ymin = ratio->GetBinContent (ratio->GetMinimumBin ()) ;
  if (ymin > 0) ymin *= 0.5 ;
  else          ymin *= 1.5 ;
  if (ymin < 0) ymin = - 0.1 ;
  TH1F * c5_frame = (TH1F *) c5->DrawFrame (200, 0., rangeScale*mass, 12.) ;
  c5_frame->SetTitle (0) ;
  c5_frame->SetStats (0) ;
  c5_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  
  ratio->Draw ("hist same") ;
  //PG plot the result of the fits on top
  TF1 * f_ratio = new TF1 ("f_ratio", ratio_crystalBallLowHighWithRise, 0, 2000, 16) ;
  f_ratio->SetParameters (params_difference) ;
  f_ratio->SetNpx (10000) ;
  f_ratio->SetLineWidth (2) ;
  f_ratio->SetLineColor (kBlue + 2) ;
  f_ratio->Draw ("same") ;
  
  c5->Print (TString ("corr_factor") + suffix, "pdf") ;
  */
  //PG plotting S only, and (SBI - B) 
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  /*
  TCanvas * c4 = new TCanvas ("c4", "c4") ;
  double ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  double ymax2 = diff->GetBinContent (diff->GetMaximumBin ()) ;
  if (ymax2 > ymax) ymax = ymax2 ;
  double ymin = delta->GetBinContent (delta->GetMinimumBin ()) ;
  if (ymin > 0) ymin *= 0.9 ;
  else          ymin *= 1.5 ;
  if (ymin < -1.1 * ymax) ymin = -1.1 * ymax ;
  TH1F * c4_frame = (TH1F *) c4->DrawFrame (200, ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_frame->SetTitle (0) ;
  c4_frame->SetStats (0) ;
  c4_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  h_MWW_mg->Draw ("histsame") ;
  diff->Draw ("histsame") ;
  delta->Draw ("histsame") ;
  c4_leg = new TLegend (0.5,0.8,0.9,0.95) ;
  c4_leg->SetFillStyle (0) ;
  c4_leg->SetBorderSize (0) ;
  c4_leg->SetTextFont (42) ;
  c4_leg->AddEntry (diff, "SBI - B","l") ;
  c4_leg->AddEntry (h_MWW_mg, "S","l") ;
  c4_leg->AddEntry (delta, "(SBI - B) - S","l") ;
  
  c4_leg->Draw () ;

  c4->Print (TString ("signals") + suffix, "pdf") ;

  //PG output of the fitting function parameters
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TString outfilename = filename.c_str () ;
  outfilename.ReplaceAll ("findInterference", "results_interference") ;
  TFile outputfile (outfilename, "recreate") ;
  outputfile.cd () ;
  f_doublePeakModel->Write () ;
  func_mg_1->Write () ;
  func_ph_1->Write () ;
  outputfile.Close () ;  

  std::cout << "scaling applied: " << scaling << std::endl ;

  return 0 ;
  */

}                                                                                   


  ////////////////////////////////////////////////////////////

void call_fit (int mass, string massa,  std::vector<vector<double> > &param,   std::vector<vector<double> > &param_error)

{

  int N=9;
  int N2=6;
  double cprime[6] = {0.1,0.3,0.5,0.7,0.9,1.};
  std::string cprime_string[6] = {"0p1","0p3","0p5","0p7","0p9","SM"};

  std::vector<double> param_temp;
  std::vector<double> param_error_temp;

  std::string stringa;

  for (int c=0; c<N2; c++)
    {
      stringa.append("findInterference.");
      stringa.append(massa);
      stringa.append(".");
      stringa.append(cprime_string[c]);
      stringa.append(".1.root");

      macro_findInterferece (stringa, mass, param_temp, param_error_temp, 1, cprime[c]); //call fit on a single point (mass,cprime)

      param.push_back(param_temp);
      param_error.push_back(param_error_temp);

      //      for (int i=0; i<9; i++) {
      //	std::std::cout<<param.at(c).at(i)<<" "<<param_error.at(c).at(i)<<std::endl;
      // }
      
      param_temp.clear();
      param_error_temp.clear();
      stringa.clear();

    }

  //  return();
}

///////////////////////////////////////////////7

int produce_plot (TString name,  double param350_new[],  double paramerror350_new[],  double param650_new[],  double paramerror650_new[],  double param800_new[],  double paramerror800_new[],  double param1000_new[],  double paramerror1000_new[], double min, double max)
{
  int Ncprime = 6;

  TString *namefile = new TString (name);
  namefile->Append(".root");
  TFile* outfile = new TFile(namefile->Data(), "RECREATE");

  double cprime [6] = {0.1,0.3,0.5,0.7,0.9,1.0};
  double errcprime [6] = {0.,0.,0.,0.,0.,0.};

  std::cout<<"DAIIIIIIIIIIIIII "<<param350_new[2]<<std::endl;

   TGraphErrors *gr_D = new TGraphErrors ( Ncprime, cprime, param350_new, errcprime, paramerror350_new);

   gr_D -> SetTitle(name.Data());
   gr_D -> GetXaxis() -> SetTitle("C'^{2}");
   gr_D -> GetYaxis() -> SetTitle("value");
   gr_D -> GetXaxis () -> CenterTitle ();
   gr_D -> GetYaxis () -> CenterTitle ();
   gr_D -> GetXaxis () -> SetRangeUser (0,1);

   gr_D -> SetLineColor (1);
   gr_D -> SetMarkerColor (1);
   gr_D -> SetMarkerStyle (20);
   gr_D -> SetMarkerSize (0.9);

   TGraphErrors *gr_2D = new TGraphErrors ( Ncprime, cprime, param650_new, errcprime, paramerror650_new);


   gr_2D -> SetLineColor (2);
   gr_2D -> SetMarkerColor (2);
   gr_2D -> SetMarkerStyle (20);
   gr_2D -> SetMarkerSize (0.9);

   TGraphErrors *gr_3D = new TGraphErrors ( Ncprime, cprime, param800_new, errcprime, paramerror800_new);

   gr_3D -> SetLineColor (4);
   gr_3D -> SetMarkerColor (4);
   gr_3D -> SetMarkerStyle (20);
   gr_3D -> SetMarkerSize (0.9);

   TGraphErrors *gr_4D = new TGraphErrors ( Ncprime, cprime, param1000_new, errcprime, paramerror1000_new);

   gr_4D -> SetLineColor (6);
   gr_4D -> SetMarkerColor (6);
   gr_4D -> SetMarkerStyle (20);
   gr_4D -> SetMarkerSize (0.9);

  gr_D->Write();
  gr_2D -> Write();
  gr_3D -> Write();
  gr_4D -> Write();

   gr_D -> SetMinimum (min);
   gr_D -> SetMaximum (max);

   TCanvas *c2 = new TCanvas(name.Data(),"");
   c2->cd();
   c2->SetTitle(name.Data());
   c2->SetGrid();
   c2->SetTickx();
   c2->SetTicky();
   gr_D -> Draw("AP");
   gr_2D -> Draw ("Psame");
   gr_3D -> Draw("Psame");
   gr_4D -> Draw("Psame");

   c2->Print(name.Data(),"png");

  c2 -> Write();

  outfile->Close();



  return(0);
}



////////////////////////////////////////////////77

int main(int argc, char *argv[])
{

  gROOT->SetBatch(1);

  int Ncprime = 6;
  int Npar = 9;

  std::vector<vector<double> > param350;
  std::vector<vector<double> > paramerror350;
  std::vector<vector<double> > param650;
  std::vector<vector<double> > paramerror650;
  std::vector<vector<double> > param800;
  std::vector<vector<double> > paramerror800;
  std::vector<vector<double> > param1000;
  std::vector<vector<double> > paramerror1000;


  call_fit(350, "350", param350, paramerror350);
  call_fit(650, "650", param650, paramerror650);
  call_fit(800, "800", param800, paramerror800);
  call_fit(1000, "1000", param1000, paramerror1000);


  TString parameters [9] = {"Norm_CB","Mean_CB_over_Higgs_mass","Sigma_CB_over_Higgs_width","alphaR_CB_times_Sigma_CB","nR_CB","alphaL_CB_times_Sigma_CB","nL_CB","Norm_Exp","Slope_EXP"};
  TString parameters_normal [9] = {"Norm_CB","Mean_CB","Sigma_CB","alphaR_CB","nR_CB","alphaL_CB","nL_CB","Norm_Exp","Slope_EXP"};
  double min[9] = {0,0.8,0,0,0,0,0,0,0};
  double max[9] = {0.001,1.2,0.6,300,12,300,8,100,100};

  double param350_sorted [6] = {0.,0.,0.,0.,0.,0.};
  double paramerror350_sorted [6] = {0.,0.,0.,0.,0.,0.};
  double param650_sorted [6] = {0.,0.,0.,0.,0.,0.};
  double paramerror650_sorted [6] = {0.,0.,0.,0.,0.,0.};
  double param800_sorted [6] = {0.,0.,0.,0.,0.,0.};
  double paramerror800_sorted [6] = {0.,0.,0.,0.,0.,0.};
  double param1000_sorted [6] = {0.,0.,0.,0.,0.,0.};
  double paramerror1000_sorted [6] = {0.,0.,0.,0.,0.,0.};

  ofstream file;
  file.open ("fit_result.txt");  //write fit results on file

  file<<"350 GeV: \t\t 650 GeV \t\t 800 GeV \t\t 1000 GeV\n";
  for (int i=0; i<Npar; i++) {
    file<<"\n"<<parameters_normal[i] <<"\n";
    for (int c=0; c<Ncprime; c++) {
      file<<param350.at(c).at(i) <<" +/- "<<paramerror350.at(c).at(i)<<"\t"<<param650.at(c).at(i) <<" +/- "<<paramerror650.at(c).at(i)<<"\t"<<param800.at(c).at(i) <<" +/- "<<paramerror800.at(c).at(i)<<"\t"<<param1000.at(c).at(i) <<" +/- "<<paramerror1000.at(c).at(i)<<"\n";
    }
  }

  file.close();
  
  for (int i=0; i<Npar; i++) {
    for (int c=0; c<Ncprime; c++) {

      if (i==1) {          //if MeanCB, plot mean_CB/Higgs_mass
	param350_sorted[c]=param350.at(c).at(i)/350.;
	paramerror350_sorted[c]=paramerror350.at(c).at(i)/350.;
	param650_sorted[c]=param650.at(c).at(i)/650.;
	paramerror650_sorted[c]=paramerror650.at(c).at(i)/650.;
	param800_sorted[c]=param800.at(c).at(i)/800.;
	paramerror800_sorted[c]=paramerror800.at(c).at(i)/800.;
	param1000_sorted[c]=param1000.at(c).at(i)/1000.;
	paramerror1000_sorted[c]=paramerror1000.at(c).at(i)/1000.;
      }
      else if (i==2) {          //if SigmaCB, plot Sigma_CB/Higgs_width
	param350_sorted[c]=param350.at(c).at(i)/15.2;
	paramerror350_sorted[c]=paramerror350.at(c).at(i)/15.2;
	param650_sorted[c]=param650.at(c).at(i)/158.;
	paramerror650_sorted[c]=paramerror650.at(c).at(i)/158.;
	param800_sorted[c]=param800.at(c).at(i)/304.;
	paramerror800_sorted[c]=paramerror800.at(c).at(i)/304.;
	param1000_sorted[c]=param1000.at(c).at(i)/647.;
	paramerror1000_sorted[c]=paramerror1000.at(c).at(i)/647.;
      }
      else if (i==3 || i==5) {          //if alphaR or alphaL, plot alpha*sigmaCB
	param350_sorted[c]=param350.at(c).at(i)*param350.at(c).at(2);
	paramerror350_sorted[c]=paramerror350.at(c).at(i)*param350.at(c).at(2);
	param650_sorted[c]=param650.at(c).at(i)*param650.at(c).at(2);
	paramerror650_sorted[c]=paramerror650.at(c).at(i)*param650.at(c).at(2);
	param800_sorted[c]=param800.at(c).at(i)*param800.at(c).at(2);
	paramerror800_sorted[c]=paramerror800.at(c).at(i)*param800.at(c).at(2);
	param1000_sorted[c]=param1000.at(c).at(i)*param1000.at(c).at(2);
	paramerror1000_sorted[c]=paramerror1000.at(c).at(i)*param1000.at(c).at(2);
      }
      else  {
	param350_sorted[c]=param350.at(c).at(i);
	paramerror350_sorted[c]=paramerror350.at(c).at(i);
	param650_sorted[c]=param650.at(c).at(i);
	paramerror650_sorted[c]=paramerror650.at(c).at(i);
	param800_sorted[c]=param800.at(c).at(i);
	paramerror800_sorted[c]=paramerror800.at(c).at(i);
	param1000_sorted[c]=param1000.at(c).at(i);
	paramerror1000_sorted[c]=paramerror1000.at(c).at(i);
      }
    }
  


    TString *namefile = new TString (parameters[i]);  //save results (graph and root files)
  namefile->Append(".root");
  TFile* outfile = new TFile(namefile->Data(), "RECREATE");

  double cprime [6] = {0.1,0.3,0.5,0.7,0.9,1.0};
  double errcprime [6] = {0.,0.,0.,0.,0.,0.};

  

   TGraphErrors *gr_D = new TGraphErrors ( Ncprime, cprime, param350_sorted, errcprime, paramerror350_sorted);

   gr_D -> SetTitle(parameters[i].Data());
   gr_D -> GetXaxis() -> SetTitle("C'^{2}");
   gr_D -> GetYaxis() -> SetTitle("value");
   gr_D -> GetXaxis () -> CenterTitle ();
   gr_D -> GetYaxis () -> CenterTitle ();
   gr_D -> GetXaxis () -> SetRangeUser (0,1);

   gr_D -> SetLineColor (1);
   gr_D -> SetMarkerColor (1);
   gr_D -> SetMarkerStyle (20);
   gr_D -> SetMarkerSize (0.9);

   TGraphErrors *gr_2D = new TGraphErrors ( Ncprime, cprime, param650_sorted, errcprime, paramerror650_sorted);


   gr_2D -> SetLineColor (2);
   gr_2D -> SetMarkerColor (2);
   gr_2D -> SetMarkerStyle (20);
   gr_2D -> SetMarkerSize (0.9);

   TGraphErrors *gr_3D = new TGraphErrors ( Ncprime, cprime, param800_sorted, errcprime, paramerror800_sorted);

   gr_3D -> SetLineColor (4);
   gr_3D -> SetMarkerColor (4);
   gr_3D -> SetMarkerStyle (20);
   gr_3D -> SetMarkerSize (0.9);

   TGraphErrors *gr_4D = new TGraphErrors ( Ncprime, cprime, param1000_sorted, errcprime, paramerror1000_sorted);

   gr_4D -> SetLineColor (6);
   gr_4D -> SetMarkerColor (6);
   gr_4D -> SetMarkerStyle (20);
   gr_4D -> SetMarkerSize (0.9);

  gr_D->Write();
  gr_2D -> Write();
  gr_3D -> Write();
  gr_4D -> Write();

   gr_D -> SetMinimum (min[i]);
   gr_D -> SetMaximum (max[i]);

   TCanvas *c2 = new TCanvas(parameters[i].Data(),"");
   c2->cd();
   c2->SetTitle(parameters[i].Data());
   c2->SetGrid();
   c2->SetTickx();
   c2->SetTicky();
   gr_D -> Draw("AP");
   gr_2D -> Draw ("Psame");
   gr_3D -> Draw("Psame");
   gr_4D -> Draw("Psame");

   c2->Print(parameters[i].Data(),"png");

  c2 -> Write();

  outfile->Close();


  }

  return(0);
}
