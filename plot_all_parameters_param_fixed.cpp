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
#include "TH2D.h"
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
#include "TMath.h"
#include "TGraph2D.h"

using namespace std;

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

TF1 *func_mg_2, *f_ratio2;

double corr_sign (double *x, double *par)
{
  return (func_mg_2->EvalPar(x,par))*(f_ratio2->EvalPar(x,par)); 
}

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

double Mean (double mass, double cprime)
{
  if (mass==349||mass==350)         return mass*(0.999716+2.00238e-05*cprime);
  else if (mass==650)    return mass*(0.996825-0.0104439*cprime);
  else if (mass==800)    return mass*(0.994242-0.0305097*cprime);
  else if (mass==1000||mass==1001)    return mass*(0.983497-0.126495*cprime);  
}

double Sigma (double mass, double cprime)
{
  if (mass==349||mass==350)         return 15.2*(0.0190926+0.384806*cprime);
  else if (mass==650)    return 158*(0.0389309+0.329887*cprime);
  else if (mass==800)    return 304*(0.0232355+0.310422*cprime);
  else if (mass==1000||mass==1001)    return 647*(0.052344+0.273591*cprime);  
}

double alphaR (double mass, double cprime)
{
  if (mass==349||mass==350)         return (0.710543+5.00614*cprime)/Sigma(mass,cprime);
  else if (mass==650)    return 12.685566682*(0.710543+5.00614*cprime)/Sigma(mass,cprime);
  else if (mass==800)    return 24.048163133*(0.710543+5.00614*cprime)/Sigma(mass,cprime);
  else if (mass==1000||mass==1001)    return 55.400735572*(0.710543+5.00614*cprime)/Sigma(mass,cprime);  
}

double nR (double mass, double cprime)
{
  if (mass==349||mass==350)         return 3.32443*(1-exp(-100.*cprime))+0.2565*cprime;
  else if (mass==650)    return 1.140059648*(3.32443*(1-exp(-3.*cprime))+0.2565*cprime);
  else if (mass==800)    return 1.210214476*(3.32443*(1-exp(-2.*cprime))+0.2565*cprime);
  else if (mass==1000||mass==1001)    return 1.303137708*(3.32443*(1-exp(-2.*cprime))+0.2565*cprime);  
}

double alphaL (double mass, double cprime)
{
  if (mass==349||mass==350)         return ( 0.670827+4.60148*cprime)/Sigma(mass,cprime);
  else if (mass==650)    return 7.98268379*( 0.670827+4.60148*cprime)/Sigma(mass,cprime);
  else if (mass==800)    return 8.106704398*( 0.670827+4.60148*cprime)/Sigma(mass,cprime);
  else if (mass==1000||mass==1001)    return 22.917553807*( 0.670827+4.60148*cprime)/Sigma(mass,cprime);  
}

double nL (double mass, double cprime)
{
  if (mass==349||mass==350)         return 1.24122+2.22144*cprime;
  else    return ((350/(1.*mass))*1.24122+(1+(350/(1.*mass)))*2.22144*cprime);
}

double Tau (double mass, double cprime)
{
  return 86.0005+25.3303*cprime;
}

double R (double mass, double cprime)
{
  if (mass==349||mass==350)         return 0.00109501+ 0.0185918*cprime;
  else if (mass==650)    return 1.1*(0.0587824+ 0.0412092*cprime);
  else if (mass==800)    return 0.01*(0.0587824+ 0.0412092*cprime);
  else if (mass==1000||mass==1001)    return 0.0001*(0.0587824+ 0.0412092*cprime);  
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
  func->SetParName (7, "R") ;
  func->SetParName (8, "Tau") ;
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
  double num = doubleGausCrystalBallLowHighPlusExp (x, par) ;    // signal and interference
  return num / den ;
  // return num;
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
  //  h_MWW_mg->Fit ("gauss_mg", "+", "", mass - 1. * h_MWW_mg->GetRMS () , mass + 1. * h_MWW_mg->GetRMS ()) ;


  /*  TH1F * h_MWW_mg_error = (TH1F*)h_MWW_mg->Clone ("h_MWW_mg_error") ;
  (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (h_MWW_mg_error, 0.68) ;
  h_MWW_mg_error->SetMarkerSize (0) ;
  h_MWW_mg_error->SetFillColor (kAzure - 9) ;
  h_MWW_mg_error->SetFillStyle (3003) ;
  */
  
  //PG second fit: first with chisq, if requested with likelihood also
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * func_mg_1 = new TF1 ("func_mg_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_mg_1->SetNpx (10000) ;
  func_mg_1->SetLineWidth (1) ;
  func_mg_1->SetLineColor (kBlue + 1) ;

  setParNamesdoubleGausCrystalBallLowHigh (func_mg_1) ;

  func_mg_1->SetParameter (0, h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ())) ;  // multiplicative scale
  if (mass==350)   
  {
    func_mg_1->FixParameter(1,350.58);
    func_mg_1->FixParameter(2,6.78);
    func_mg_1->FixParameter(3,1.10);
    func_mg_1->FixParameter(4,1.57);
    func_mg_1->FixParameter(5,1.11);
    func_mg_1->FixParameter(6,2.71);
  }

  else if (mass==650)   
  {
    func_mg_1->FixParameter(1,663.58);
    func_mg_1->FixParameter(2,62.44);
    func_mg_1->FixParameter(3,0.80);
    func_mg_1->FixParameter(4,3.01);
    func_mg_1->FixParameter(5,0.92);
    func_mg_1->FixParameter(6,100);
  }

  else if (mass==800)   
  {
    func_mg_1->FixParameter(1,820.14);
    func_mg_1->FixParameter(2,117.22);
    //  func_mg_1->SetParLimits(2,116,118);
    func_mg_1->FixParameter(3,0.79);
    // func_mg_1->SetParLimits(3,0.75,0.8);
    func_mg_1->FixParameter(4,6.78);
    // func_mg_1->SetParLimits(4,6.5,6.9);
    func_mg_1->FixParameter(5,1.19);
    // func_mg_1->SetParLimits(5,1.15,1.24);
    func_mg_1->FixParameter(6,147);
    // func_mg_1->SetParLimits(6,140,150);
  }

  else if (mass==1000)   
  {
    func_mg_1->FixParameter(1,1044);
    func_mg_1->FixParameter(2,221);
    func_mg_1->FixParameter(3,0.91);
    func_mg_1->FixParameter(4,100);
    func_mg_1->FixParameter(5,1.45);
    func_mg_1->FixParameter(6,100);
  }


  //PG set the range of the fit
  int sign = 1 ;
  if (mass < 400) sign = -2 ;
  if (mass > 700) sign = -1.5 ;
  if (mass > 850) sign = -1. ;
  std::cout << "-------------------\nFITTING THE MADGRAPH SIGNAL\n\n-------------------\n" ;
  TFitResultPtr result = h_MWW_mg->Fit ("func_mg_1", "+", "", 100, 2 * mass) ;
  //  TFitResultPtr result = h_MWW_mg->Fit ("func_mg_1", "+L", "", 0.5 * mass + sign * 50, 2 * mass) ;
  if (useLikelihood && mass < 800)
    {
      std::cout << "-------------------\nFITTING THE MADGRAPH SIGNAL W/ LIKELIHOOD\n\n-------------------\n" ;
      func_mg_1->SetParameters (func_mg_1->GetParameters ()) ; //PG not sure this is necessary
      func_mg_1->SetLineColor (kBlue + 3) ;
      h_MWW_mg->Fit ("func_mg_1", "+L", "", 0.5 * mass - 50, 2 * mass) ;
    }


  TH1F * h_MWW_mg_error = (TH1F*)h_MWW_mg->Clone ("h_MWW_mg_error") ;
  (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (h_MWW_mg_error, 0.68) ;
  h_MWW_mg_error->SetMarkerSize (0) ;
  h_MWW_mg_error->SetFillColor (kBlue) ;
  h_MWW_mg_error->SetFillStyle (3002) ;


  //PG plotting the result
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  double ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  double ymin = h_MWW_mg->GetBinContent (h_MWW_mg->GetMinimumBin ()) ;
  if (ymin <= 0) ymin = ymax / 10000. ; 
  TH1F *c4_mg_frame;
  if (mass==350)
    c4_mg_frame = (TH1F *) c4_mg->DrawFrame (320, 0.9 * ymin, 380, 1.1 * ymax) ;
  else
    c4_mg_frame = (TH1F *) c4_mg->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
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
  c4_mg->Print (TString ("signals_mg_log") + suffix, "png") ;

  c4_mg->SetLogy (0) ;
  c4_mg->Print (TString ("signals_mg_lin") + suffix, "png") ;

  return func_mg_1 ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//


//PG no left rise
TF1 * FIT_phantom_signal (TH1F * diff, double mass, double rangeScale, TString suffix,  double cprime, bool useLikelihood = false)
{  
 TCanvas * c4_ph = new TCanvas ("c4_ph", "c4_ph") ;

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


  /*	double width = 0.;

	if (mass==350)  width=15.2;
	else if (mass==650)  width=158.;
	else if (mass==800)  width=304.;
	else if (mass==1000) width=647.;

          func_ph_1->SetParameter (0, 0.01*0.000433428+8.80617e-05*cprime) ;  // multiplicative scale                                                
      func_ph_1->SetParLimits (0, 0, 0.006);
      func_ph_1->SetParameter (7, 1.1*(0.0587824+ 0.0412092*cprime));
	func_ph_1->SetParLimits (7,0,0.1);

    func_ph_1->FixParameter(1,mass*(1+(0.177504-0.00103068*mass+1.90992e-06*mass*mass-1.18324e-09*mass*mass*mass)*cprime));
        func_ph_1->FixParameter (2, width*((0.00416317+4.1768e-05*mass)+(0.443076-0.000169142*mass)*cprime)) ;
	func_ph_1->FixParameter (3, (21.6235-0.107433*mass+0.000140711*mass*mass-0.25)*(0.710543+5.00614*cprime)/func_ph_1->GetParameter(2)) ;
	 func_ph_1->FixParameter (4,(0.836811+0.000466488*mass)*(3.32443*(1-exp((-268.667+0.986813*mass-0.00120366*mass*mass+4.83516e-07*mass*mass*mass)*cprime))+0.2565*cprime));
      func_ph_1->FixParameter (5, (8.5992-0.0386763*mass+5.21299e-05*mass*mass-0.44)*( 0.670827+4.60148*cprime)/func_ph_1->GetParameter(2));
      func_ph_1->FixParameter (6, ((350/(1.*mass))*1.24122+(0.620881+0.000332256*mass+1.95247e-06*mass*mass)*2.22144*cprime));
     func_ph_1->FixParameter (8, 86.0005+25.3303*cprime);


     //	   func_ph_1->FixParameter (4,1.210214476*(3.32443*(1-exp(-2.*cprime))+0.2565*cprime));

    //   func_ph_1->FixParameter (2, 158*(0.0389309+0.329887*cprime)) ;
	 //	 func_ph_1->FixParameter (3, 10.184622795*1.245560777*(0.710543+5.00614*cprime)/func_ph_1->GetParameter(2)) ;
	   //	 	 func_ph_1->FixParameter (4,1.140059648*(3.32443*(1-exp(-3.*cprime))+0.2565*cprime));
		 //  func_ph_1->FixParameter (5, 10.184622795*0.783797687*( 0.670827+4.60148*cprime)/func_ph_1->GetParameter(2));
		 //  func_ph_1->FixParameter (6, ((350/(1.*mass))*1.24122+(1+(350/(1.*mass)))*2.22144*cprime));
		 //  func_ph_1->FixParameter (8, 86.0005+25.3303*cprime);
		 */

  func_ph_1->FixParameter(1,Mean(mass,cprime));
  func_ph_1->FixParameter(2,Sigma(mass,cprime));
  func_ph_1->FixParameter(3,alphaR(mass,cprime));
  func_ph_1->FixParameter(4,nR(mass,cprime));
  func_ph_1->FixParameter(5,alphaL(mass,cprime));
  func_ph_1->FixParameter(6,nL(mass,cprime));
  func_ph_1->FixParameter(7,R(mass,cprime));
  func_ph_1->FixParameter(8,Tau(mass,cprime));

  if (mass==350)     func_ph_1->SetParameter(0,0.000433428+8.80617e-05*cprime);
  else if (mass==650) {
      func_ph_1->SetParameter (0, 0.01*0.000433428+8.80617e-05*cprime) ;  // multiplicative scale                                                
      func_ph_1->SetParLimits (0, 0, 0.00006);
  }
  else if (mass==800)     func_ph_1->SetParameter (0, 0.005250714*(0.000433428+8.80617e-05*cprime)) ;  // multiplicative scale         
  else if (mass==1000) {
      func_ph_1->SetParameter (0, 0.01*0.000433428+8.80617e-05*cprime) ;  // multiplicative scale                                                
      func_ph_1->SetParLimits (0, 0, 0.00006);
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


// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====


int macro_findInterferece (string filename, double mass, std::vector<double> & param, std::vector<double> & param_error, std::vector<double> & param_mg, std::vector<double> & param_error_mg, int initialRebin=1, double cprime=1., bool useLeftRise = false)
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
  else if (mass > 610 && mass < 790) reBin *= 8 ;
  //  else if (mass > 790 && mass < 860 && cprime==1.) reBin *= 20 ;
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

  h_MWW_mg->Scale(1./cprime);

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

  TLegend * c1_leg = new TLegend (0.6,0.8,0.9,0.95) ;
  c1_leg->SetFillStyle (0) ;
  c1_leg->SetBorderSize (0) ;
  c1_leg->SetTextFont (42) ;
  c1_leg->AddEntry (h_MWW_phbkg, "B","l") ;
  c1_leg->AddEntry (h_MWW_phbkgsig, "SBI","l") ;
  c1_leg->AddEntry (h_MWW_mg, "S","l") ;
  c1_leg->Draw () ;
  
  c1->Print (TString ("spectra") + suffix, "png") ;
  
  //PG S only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  func_mg_2 = FIT_madgraph_signal (h_MWW_mg, mass, rangeScale, suffix, false) ;

  for (int i=0; i<7; i++) {
    param_mg.push_back(func_mg_2->GetParameter(i));
    param_error_mg.push_back(func_mg_2->GetParError(i));
    //       std::cout<<param_mg.back()<<std::endl;
    //  getchar();

    //  std::std::cout<<param.at(i)<<" "<<param_error.at(i)<<std::endl;
  }



  //PG (SBI - B) only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TF1 * func_ph_1 ;
  //  if (useLeftRise == true) func_ph_1 = FIT_phantom_signal_2 (diff, mass, rangeScale, suffix, false) ;
  func_ph_1 = FIT_phantom_signal   (diff, mass, rangeScale, suffix, cprime, false) ;

  for (int i=0; i<9; i++) {
    param.push_back(func_ph_1->GetParameter(i));
    param_error.push_back(func_ph_1->GetParError(i));
    //  std::std::cout<<param.at(i)<<" "<<param_error.at(i)<<std::endl;
  }



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
      func_mg_2->GetParameter (2) * func_mg_2->GetParameter (2)  
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
  func_mg_2->GetParameters (params_difference + 9) ;
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

  //  ratio      ->Rebin (reBin) ;
  // normaliseToBinWidth (ratio) ;


  TCanvas * c5 = new TCanvas ("c5", "c5") ;
  ymax = ratio->GetBinContent (ratio->GetMaximumBin ()) ;
  ymin = ratio->GetBinContent (ratio->GetMinimumBin ()) ;
  //  if (ymin > 0) ymin *= 0.5 ;
  // else          ymin *= 1.5 ;
  //  if (ymin < 0) ymin = - 0.1 ;

  if (mass==350) ymax=500;
  else if (mass==650 || mass==800) ymax=2000;
  else ymax=8000;

  c5->SetLogy () ;

  //  TH1F * c5_frame = (TH1F *) c5->DrawFrame (200, 0.1, rangeScale*mass, 10) ;
  TH1F * c5_frame = (TH1F *) c5->DrawFrame (200, 0.001, rangeScale*mass, ymax) ;
  c5_frame->SetTitle (0) ;
  c5_frame->SetStats (1) ;
  c5_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;

  //  func_mg_2->FixParameter(2,func_mg_2->GetParameter(2)*cprime);

  Double_t params_difference[16] ;
  func_ph_1->GetParameters (params_difference) ;
  func_mg_2->GetParameters (params_difference + 9) ;
  
  ratio->Draw ("EPsame") ;
  //PG plot the result of the fits on top
  TF1 * f_ratio = new TF1 ("f_ratio", ratio_crystalBallLowHighWithRise, 0, 2000, 16) ;
  f_ratio->SetParameters (params_difference) ;
  f_ratio->SetNpx (10000) ;
  f_ratio->SetLineWidth (2) ;
  f_ratio->SetLineColor (kBlue + 2) ;
  f_ratio->Draw ("same") ;


  //How to correctly rescale the signal to the bsm case??? not sure...

  //  func_mg_2->FixParameter(2,func_mg_2->GetParameter(2)*cprime);
  //  func_mg_2->FixParameter(0,func_mg_2->GetParameter(0)/cprime);
  //  func_mg_2->FixParameter(3,func_mg_2->GetParameter(3)*cprime); ??
  //  func_mg_2->FixParameter(5,func_mg_2->GetParameter(5)*cprime); ??

  Double_t params_difference2[16] ;
  func_ph_1->GetParameters (params_difference2) ;
  func_mg_2->GetParameters (params_difference2 + 9) ;
  
  f_ratio2 = new TF1 ("f_ratio2", ratio_crystalBallLowHighWithRise, 0, 2000, 16) ;
  f_ratio2->SetParameters (params_difference2) ;
  f_ratio2->SetNpx (10000) ;
  f_ratio2->SetLineWidth (2) ;
  f_ratio2->SetLineColor (kRed + 2) ;
  //  f_ratio2->Draw ("same") ;

  c5->Update () ;
  c5->Print (TString ("corr_factor_log") + suffix, "png") ;


  TCanvas * c6 = new TCanvas ("c6", "c6") ;

  TH1F * c6_frame = (TH1F *) c6->DrawFrame (200, 0, rangeScale*mass, 10) ;
  c6_frame->SetTitle (0) ;
  c6_frame->SetStats (1011) ;
  c6_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  ratio->Draw ("EPsame") ;
  f_ratio->Draw ("same") ;

  c6->Print (TString ("corr_factor_lin") + suffix, "png") ;



  
  //  c5->Print (TString ("corr_factor") + suffix, "png") ;


  //plot S after correction for interference

  /*  
  TCanvas * c18_mg = new TCanvas ("c18_mg", "c18_mg") ;
  ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  ymin = h_MWW_mg->GetBinContent (h_MWW_mg->GetMinimumBin ()) ;
  if (ymin <= 0) ymin = ymax / 10000. ; 

  //  TH1F * c18_mg_frame = (TH1F *) c18_mg->DrawFrame (200, 0.001, rangeScale * mass, 1.) ;
    TH1F * c18_mg_frame = (TH1F *) c18_mg->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c18_mg->SetLogy () ;
  c18_mg_frame->SetTitle (0) ;
  c18_mg_frame->SetStats (0) ;
  c18_mg_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  //  h_MWW_mg_error->Draw ("e3same") ;
  h_MWW_mg->Draw ("EPsame") ;
  func_mg_2->SetLineColor(kRed);

  TF1 * corr_sign2 = new TF1 ("corr_sign", corr_sign, 0, 2000, 16) ;
  corr_sign2->SetLineColor(kGreen);
  // corr_sign2->Draw("same");

  func_mg_2->Draw ("same");
      
  c18_mg->Update () ;
  c18_mg->Print (TString ("corrected_signal_log") + suffix, "png") ;

  c18_mg->SetLogy (0) ;
  c18_mg->Print (TString ("corrected_signal_lin") + suffix, "png") ;
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
  func_mg_2->Write () ;
  func_ph_1->Write () ;
  outputfile.Close () ;  

  std::cout << "scaling applied: " << scaling << std::endl ;

  return 0 ;
  */

  return(0);

}                                                                                   


  ////////////////////////////////////////////////////////////

void call_fit (int mass, string massa,  std::vector<vector<double> > &param,   std::vector<vector<double> > &param_error, std::vector<vector<double> > &param_mg, std::vector<vector<double> > &param_error_mg)

{

  int N=9;
  int N2=6;
  double cprime[6] = {0.1,0.3,0.5,0.7,0.9,1.};
  std::string cprime_string[6] = {"0p1","0p3","0p5","0p7","0p9","SM"};

  std::vector<double> param_temp;
  std::vector<double> param_error_temp;
  std::vector<double> param_temp_mg;
  std::vector<double> param_error_temp_mg;

  std::string stringa;

  for (int c=0; c<N2; c++)
    {
      stringa.append("findInterference.");
      stringa.append(massa);
      stringa.append(".");
      stringa.append(cprime_string[c]);
      stringa.append(".1.root");

      macro_findInterferece (stringa, mass, param_temp, param_error_temp,  param_temp_mg, param_error_temp_mg, 1, cprime[c]); //call fit on a single point (mass,cprime)

      //      getchar();
      param.push_back(param_temp);
      param_error.push_back(param_error_temp);

      param_mg.push_back(param_temp_mg);
      param_error_mg.push_back(param_error_temp_mg);

      /*            for (int i=0; i<7; i++) {
	      if (i==2);
      	std::cout<<param_mg.at(c).at(i)<<" "<<param_error.at(c).at(i)<<std::endl;
	getchar();
             }
      */
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
  int Ncprime_alt = 15;

  double cprime [6] = {0.1,0.3,0.5,0.7,0.9,1.0};
  double errcprime [6] = {0.,0.,0.,0.,0.,0.};
    double cprime_alt [2] = {0.01,1.05};
  //   double cprime_alt [15] = {0.05,0.15,0.2,0.25,0.35,0.4,0.45,0.55,0.6,0.65,0.75,0.8,0.85,0.95,1.05};


  std::vector<vector<double> > param350;
  std::vector<vector<double> > paramerror350;
  std::vector<vector<double> > param650;
  std::vector<vector<double> > paramerror650;
  std::vector<vector<double> > param800;
  std::vector<vector<double> > paramerror800;
  std::vector<vector<double> > param1000;
  std::vector<vector<double> > paramerror1000;

  std::vector<vector<double> > param350_mg;
  std::vector<vector<double> > paramerror350_mg;
  std::vector<vector<double> > param650_mg;
  std::vector<vector<double> > paramerror650_mg;
  std::vector<vector<double> > param800_mg;
  std::vector<vector<double> > paramerror800_mg;
  std::vector<vector<double> > param1000_mg;
  std::vector<vector<double> > paramerror1000_mg;


  call_fit(350, "350", param350, paramerror350, param350_mg, paramerror350_mg);
  call_fit(650, "650", param650, paramerror650, param650_mg, paramerror650_mg);
  call_fit(800, "800", param800, paramerror800, param800_mg, paramerror800_mg);
  call_fit(1000, "1000", param1000, paramerror1000, param1000_mg, paramerror1000_mg);


  TString parameters [9] = {"Norm","Mean_CB_over_Higgs_mass","Sigma_CB_over_Higgs_width","alphaR_CB_times_Sigma_CB","nR_CB","alphaL_CB_times_Sigma_CB","nL_CB","R","Tau"};
    // TString parameters [9] = {"Norm","Mean_CB_over_Higgs_mass","Sigma_CB_over_Higgs_width","alphaR_CB_times_Sigma_CB","nR_CB","alphaL_CB_times_Sigma_CB","nL_CB","R","Tau"};
  TString parameters_normal [9] = {"Norm","Mean_CB","Sigma_CB","alphaR_CB","nR_CB","alphaL_CB","nL_CB","R","Tau"};
  double min[9] = {-16,0.8,0,0,0,0,0,0,0};
  double max[9] = {-6,1.2,0.6,300,12,300,8,100,100};

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

                  if (i==0) {          //if Norm, plot LogNorm
	param350_sorted[c]=TMath::Log(param350.at(c).at(i));
	paramerror350_sorted[c]=TMath::Log(paramerror350.at(c).at(i));
	param650_sorted[c]=TMath::Log(param650.at(c).at(i));
	paramerror650_sorted[c]=TMath::Log(paramerror650.at(c).at(i));
	param800_sorted[c]=TMath::Log(param800.at(c).at(i));
	paramerror800_sorted[c]=TMath::Log(paramerror800.at(c).at(i));
	param1000_sorted[c]=TMath::Log(param1000.at(c).at(i));
	paramerror1000_sorted[c]=TMath::Log(paramerror1000.at(c).at(i));
      }
      
		  /*	                if (i==0) {          //if Norm, plot LogNorm
	param350_sorted[c]=(param350.at(c).at(i));
	paramerror350_sorted[c]=(paramerror350.at(c).at(i));
	param650_sorted[c]=(param650.at(c).at(i));
	paramerror650_sorted[c]=(paramerror650.at(c).at(i));
	param800_sorted[c]=(param800.at(c).at(i));
	paramerror800_sorted[c]=(paramerror800.at(c).at(i));
	param1000_sorted[c]=(param1000.at(c).at(i));
	paramerror1000_sorted[c]=(paramerror1000.at(c).at(i));
      }
		  */
      else if (i==1) {          //if MeanCB, plot mean_CB/Higgs_mass
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

  
   TGraphErrors *gr_D = new TGraphErrors ( Ncprime, cprime, param350_sorted, errcprime, paramerror350_sorted);

   gr_D -> SetTitle(parameters[i].Data());
   gr_D -> GetXaxis() -> SetTitle("C'^{2}");
   gr_D -> GetYaxis() -> SetTitle("value");
   gr_D -> GetXaxis () -> CenterTitle ();
   gr_D -> GetYaxis () -> CenterTitle ();
   //   gr_D -> GetXaxis () -> SetRangeUser (0.001,1);

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

  ofstream filebis[9];

  for (int i=0; i<Npar; i++) {

    TString *namefile = new TString (parameters_normal[i]);  
    namefile->Append("_SI.txt");

    filebis[i].open (namefile->Data());

    for (int c=0; c<Ncprime_alt; c++) {
      if (parameters_normal[i].Contains("Mean_CB")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<Mean(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<Mean(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<Mean(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<Mean(1001,cprime_alt[c])<<endl;
	}
      else if (parameters_normal[i].Contains("Sigma_CB")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<Sigma(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<Sigma(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<Sigma(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<Sigma(1001,cprime_alt[c])<<endl;
	}
      else if (parameters_normal[i].Contains("alphaR_CB")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<alphaR(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<alphaR(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<alphaR(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<alphaR(1001,cprime_alt[c])<<endl;
	}
      else if (parameters_normal[i].Contains("alphaL_CB")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<alphaL(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<alphaL(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<alphaL(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<alphaL(1001,cprime_alt[c])<<endl;
	}
      else if (parameters_normal[i].Contains("nL_CB")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<nL(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<nL(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<nL(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<nL(1001,cprime_alt[c])<<endl;
	}
      else if (parameters_normal[i].Contains("nR_CB")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<nR(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<nR(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<nR(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<nR(1001,cprime_alt[c])<<endl;
	}
      else if (parameters_normal[i].Contains("R")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<R(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<R(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<R(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<R(1001,cprime_alt[c])<<endl;
	}
      else if (parameters_normal[i].Contains("Tau")) {
	  filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<Tau(349,cprime_alt[c])<<endl;
	  filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<Tau(650,cprime_alt[c])<<endl;
	  filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<Tau(800,cprime_alt[c])<<endl;
	  filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<Tau(1001,cprime_alt[c])<<endl;
	}
      else {
        int c2=0;
	if (c==0)  c2=0;
	else if (c==1)  c2=5;
	filebis[i]<<349<<" "<<cprime_alt[c]<<" "<<TMath::Log(param350.at(c2).at(i))<<endl;
	filebis[i]<<650<<" "<<cprime_alt[c]<<" "<<TMath::Log(param650.at(c2).at(i))<<endl;
	filebis[i]<<800<<" "<<cprime_alt[c]<<" "<<TMath::Log(param800.at(c2).at(i))<<endl;
	filebis[i]<<1001<<" "<<cprime_alt[c]<<" "<<TMath::Log(param1000.at(c2).at(i))<<endl;
      }
         }
    
    for (int c=0; c<Ncprime; c++) {
            if (parameters_normal[i].Contains("Norm")) {
	  filebis[i]<<350<<" "<<cprime[c]<<" "<<TMath::Log(param350.at(c).at(i))<<endl;
	  filebis[i]<<650<<" "<<cprime[c]<<" "<<TMath::Log(param650.at(c).at(i))<<endl;
	  filebis[i]<<800<<" "<<cprime[c]<<" "<<TMath::Log(param800.at(c).at(i))<<endl;
	  filebis[i]<<1000<<" "<<cprime[c]<<" "<<TMath::Log(param1000.at(c).at(i))<<endl;
	}
      else {
      	filebis[i]<<350<<" "<<cprime[c]<<" "<<param350.at(c).at(i)<<endl;
	filebis[i]<<650<<" "<<cprime[c]<<" "<<param650.at(c).at(i)<<endl;
	filebis[i]<<800<<" "<<cprime[c]<<" "<<param800.at(c).at(i)<<endl;
	filebis[i]<<1000<<" "<<cprime[c]<<" "<<param1000.at(c).at(i)<<endl;
		}

      }

    filebis[i]<<-1;
    filebis[i].close();

  }


  ofstream filetris[7];

  for (int i=0; i<(Npar-2); i++) {

    TString *namefile = new TString (parameters_normal[i]);  
    namefile->Append("_S.txt");

    filetris[i].open (namefile->Data()); 

    for (int c=0; c<Ncprime; c++) {
        if (parameters_normal[i].Contains("Norm")) {
	  filetris[i]<<350<<" "<<cprime[c]<<" "<<TMath::Log(param350_mg.at(c).at(i))<<endl;
	  filetris[i]<<650<<" "<<cprime[c]<<" "<<TMath::Log(param650_mg.at(c).at(i))<<endl;
	  filetris[i]<<800<<" "<<cprime[c]<<" "<<TMath::Log(param800_mg.at(c).at(i))<<endl;
	  filetris[i]<<1000<<" "<<cprime[c]<<" "<<TMath::Log(param1000_mg.at(c).at(i))<<endl;
	}
      else {
        filetris[i]<<350<<" "<<cprime[c]<<" "<<param350_mg.at(c).at(i)<<endl;
      filetris[i]<<650<<" "<<cprime[c]<<" "<<param650_mg.at(c).at(i)<<endl;
      filetris[i]<<800<<" "<<cprime[c]<<" "<<param800_mg.at(c).at(i)<<endl;
      filetris[i]<<1000<<" "<<cprime[c]<<" "<<param1000_mg.at(c).at(i)<<endl;
      	}
    }

    int c2=0;
    for (int c=0; c<Ncprime_alt; c++) {
        if (c==0)  c2=0;
      else if (c==1)  c2=5;
      if (parameters_normal[i].Contains("Norm")) {
	  filetris[i]<<349<<" "<<cprime_alt[c]<<" "<<TMath::Log(param350_mg.at(c2).at(i))<<endl;
	  filetris[i]<<650<<" "<<cprime_alt[c]<<" "<<TMath::Log(param650_mg.at(c2).at(i))<<endl;
	  filetris[i]<<800<<" "<<cprime_alt[c]<<" "<<TMath::Log(param800_mg.at(c2).at(i))<<endl;
	  filetris[i]<<1001<<" "<<cprime_alt[c]<<" "<<TMath::Log(param1000_mg.at(c2).at(i))<<endl;
	}
      else {
         filetris[i]<<349<<" "<<cprime_alt[c]<<" "<<param350_mg.at(c2).at(i)<<endl;
      filetris[i]<<650<<" "<<cprime_alt[c]<<" "<<param650_mg.at(c2).at(i)<<endl;
      filetris[i]<<800<<" "<<cprime_alt[c]<<" "<<param800_mg.at(c2).at(i)<<endl;
      filetris[i]<<1001<<" "<<cprime_alt[c]<<" "<<param1000_mg.at(c2).at(i)<<endl;
      	}
    }

    filetris[i]<<-1;
    filetris[i].close();

  }



  ////open parameter files and create TGraph2D

  TString *namefile_SI = new TString ("file_for_interpolation.root");
  TFile* outfile_SI = new TFile(namefile_SI->Data(), "RECREATE");

  ifstream file_SI[9];
  TGraph2D *graph_SI[9];
  TH2D* histo_SI[9];
  double mass,c,value;
  /*
  Double_t fill_param1[16];
  Double_t fill_param2[16];
  Double_t fill_param3[16];
  Double_t fill_param4[16];
  Double_t fill_param5[16];
  Double_t fill_param6[16];
  Double_t fill_param7[16];
  Double_t fill_param8[16];
  Double_t fill_param9[16];
  Double_t fill_param10[16];
  Double_t fill_param11[16];
  Double_t fill_param12[16];
  Double_t fill_param13[16];
  Double_t fill_param14[16];
  Double_t fill_param15[16];
  Double_t fill_param16[16];
  Double_t fill_param17[16];
  Double_t fill_param18[16];
  Double_t fill_param19[16];
  Double_t fill_param20[16];
  Double_t fill_param21[16];
  Double_t fill_param22[16];
  Double_t fill_param23[16];
  Double_t fill_param24[16];
  */
  for (int i=0; i<Npar; i++) {

    TString *namefile = new TString (parameters_normal[i]);  
    namefile->Append("_SI.txt");

    TString *namefile2 = new TString (parameters_normal[i]);  
    namefile2->Append("_SI_histo");

    cout<<namefile->Data()<<endl;
    file_SI[i].open (namefile->Data());  
    graph_SI[i] = new TGraph2D(namefile->Data(),namefile->Data());
    //    graph_SI[i]->SetNpx(4);
    //  graph_SI[i]->SetNpy(21);
    histo_SI[i] = new TH2D(namefile2->Data(),namefile2->Data(),4,340,1010,21,0,1.1);
    //    graph_SI[i]->SetHistogram(histo_SI[i]);
    
    int count=0;
    while (!file_SI[i].eof()) {
      file_SI[i]>>mass>>c>>value;
      if (mass==-1) continue;
      graph_SI[i]->SetPoint(count,mass,c,value);
      /* if (parameters_normal[i].Contains("Norm")){       value=exp(value); }
           if (mass==350 && c==0.1) {  fill_param1[i]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.3) {  fill_param2[i]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.5) {  fill_param3[i]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.7) {  fill_param4[i]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.9) {  fill_param5[i]=value;  cout<<value<<endl; }
      if (mass==350 && c==1.0) {  fill_param6[i]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.1) {  fill_param7[i]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.3) {  fill_param8[i]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.5) {  fill_param9[i]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.7) {  fill_param10[i]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.9) {  fill_param11[i]=value;  cout<<value<<endl; }
      if (mass==650 && c==1.0) {  fill_param12[i]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.1) {  fill_param13[i]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.3) {  fill_param14[i]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.5) {  fill_param15[i]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.7) {  fill_param16[i]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.9) {  fill_param17[i]=value;  cout<<value<<endl; }
      if (mass==800 && c==1.0) {  fill_param18[i]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.1) {  fill_param19[i]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.3) {  fill_param20[i]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.5) {  fill_param21[i]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.7) {  fill_param22[i]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.9) {  fill_param23[i]=value;  cout<<value<<endl; }
      if (mass==1000 && c==1.0) {  fill_param24[i]=value;  cout<<value<<endl; }
      */
      count++;
    }

    graph_SI[i]->Write(namefile->Data());
    //    histo_SI[i]->Write(namefile2->Data());
    file_SI[i].close();
  }


  ifstream file_S[7];
  TGraph2D *graph_S[7];
  TH2D *histo_S[9];

  for (int i=0; i<Npar-2; i++) {

    TString *namefile = new TString (parameters_normal[i]);  
    namefile->Append("_S.txt");

    TString *namefile2 = new TString (parameters_normal[i]);  
    namefile2->Append("_S_histo");

    file_S[i].open (namefile->Data()); 
    graph_S[i] = new TGraph2D(namefile->Data(),namefile->Data());
    //    graph_S[i]->SetNpx(4);
    //  graph_S[i]->SetNpy(21);
    histo_S[i] = new TH2D(namefile2->Data(),namefile2->Data(),4,340,1010,21,0,1.1);
    //    graph_S[i]->SetHistogram(histo_S[i]);

    int count=0;
    while (!file_S[i].eof()) {
      file_S[i]>>mass>>c>>value;
      if (mass==-1) continue;
      graph_S[i]->SetPoint(count,mass,c,value);
      /*       if (parameters_normal[i].Contains("Norm"))   value=exp(value);
      if (mass==350 && c==0.1) {  fill_param1[i+9]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.3) {  fill_param2[i+9]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.5) {  fill_param3[i+9]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.7) {  fill_param4[i+9]=value;  cout<<value<<endl; }
      if (mass==350 && c==0.9) {  fill_param5[i+9]=value;  cout<<value<<endl; }
      if (mass==350 && c==1.0) {  fill_param6[i+9]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.1) {  fill_param7[i+9]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.3) {  fill_param8[i+9]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.5) {  fill_param9[i+9]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.7) {  fill_param10[i+9]=value;  cout<<value<<endl; }
      if (mass==650 && c==0.9) {  fill_param11[i+9]=value;  cout<<value<<endl; }
      if (mass==650 && c==1.0) {  fill_param12[i+9]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.1) {  fill_param13[i+9]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.3) {  fill_param14[i+9]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.5) {  fill_param15[i+9]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.7) {  fill_param16[i+9]=value;  cout<<value<<endl; }
      if (mass==800 && c==0.9) {  fill_param17[i+9]=value;  cout<<value<<endl; }
      if (mass==800 && c==1.0) {  fill_param18[i+9]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.1) {  fill_param19[i+9]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.3) {  fill_param20[i+9]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.5) {  fill_param21[i+9]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.7) {  fill_param22[i+9]=value;  cout<<value<<endl; }
      if (mass==1000 && c==0.9) {  fill_param23[i+9]=value;  cout<<value<<endl; }
      if (mass==1000 && c==1.0) {  fill_param24[i+9]=value;  cout<<value<<endl; }
      */
      count++;
    }

    graph_S[i]->Write(namefile->Data());
    //    histo_S[i]->Write(namefile2->Data());
    file_S[i].close();
  }

    outfile_SI->Close();

  // for cross-check
    /*  TF1* f_ratio23 = new TF1 ("f_ratio23", ratio_crystalBallLowHighWithRise, 0, 2000, 16) ;

  f_ratio23->SetParameters (fill_param1) ;
  std::cout<<f_ratio23->Eval(350)<<std::endl;
  f_ratio23->SetParameters (fill_param2) ;
  std::cout<<f_ratio23->Eval(350)<<std::endl;
  f_ratio23->SetParameters (fill_param3) ;
  std::cout<<f_ratio23->Eval(350)<<std::endl;
  f_ratio23->SetParameters (fill_param4) ;
  std::cout<<f_ratio23->Eval(350)<<std::endl;
  f_ratio23->SetParameters (fill_param5) ;
  std::cout<<f_ratio23->Eval(350)<<std::endl;
  f_ratio23->SetParameters (fill_param6) ;
  std::cout<<f_ratio23->Eval(350)<<std::endl;

  f_ratio23->SetParameters (fill_param7) ;
  std::cout<<f_ratio23->Eval(650)<<std::endl;
  f_ratio23->SetParameters (fill_param8) ;
  std::cout<<f_ratio23->Eval(650)<<std::endl;
  f_ratio23->SetParameters (fill_param9) ;
  std::cout<<f_ratio23->Eval(650)<<std::endl;
  f_ratio23->SetParameters (fill_param10) ;
  std::cout<<f_ratio23->Eval(650)<<std::endl;
  f_ratio23->SetParameters (fill_param11) ;
  std::cout<<f_ratio23->Eval(650)<<std::endl;
  f_ratio23->SetParameters (fill_param12) ;
  std::cout<<f_ratio23->Eval(650)<<std::endl;

  f_ratio23->SetParameters (fill_param13) ;
  std::cout<<f_ratio23->Eval(800)<<std::endl;
  f_ratio23->SetParameters (fill_param14) ;
  std::cout<<f_ratio23->Eval(800)<<std::endl;
  f_ratio23->SetParameters (fill_param15) ;
  std::cout<<f_ratio23->Eval(800)<<std::endl;
  f_ratio23->SetParameters (fill_param16) ;
  std::cout<<f_ratio23->Eval(800)<<std::endl;
  f_ratio23->SetParameters (fill_param17) ;
  std::cout<<f_ratio23->Eval(800)<<std::endl;
  f_ratio23->SetParameters (fill_param18) ;
  std::cout<<f_ratio23->Eval(800)<<std::endl;

  f_ratio23->SetParameters (fill_param19) ;
  std::cout<<f_ratio23->Eval(1000)<<std::endl;
  f_ratio23->SetParameters (fill_param20) ;
  std::cout<<f_ratio23->Eval(1000)<<std::endl;
  f_ratio23->SetParameters (fill_param21) ;
  std::cout<<f_ratio23->Eval(1000)<<std::endl;
  f_ratio23->SetParameters (fill_param22) ;
  std::cout<<f_ratio23->Eval(1000)<<std::endl;
  f_ratio23->SetParameters (fill_param23) ;
  std::cout<<f_ratio23->Eval(1000)<<std::endl;
  f_ratio23->SetParameters (fill_param24) ;
  std::cout<<f_ratio23->Eval(1000)<<std::endl;
  
    */

return(0);

}

