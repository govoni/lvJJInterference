/*
r00t -q macro_findInterferece.C\(\"findInterference.350.root\",350,1\)
r00t -q macro_findInterferece.C\(\"findInterference.500.root\",500,2\)
r00t -q macro_findInterferece.C\(\"findInterference.650.root\",650,4\)
r00t -q macro_findInterferece.C\(\"findInterference.800.root\",800,4\)
r00t -q macro_findInterferece.C\(\"findInterference.1000.root\",1000,2\)
*/


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double max (double uno, double due)
{
  if (uno > due) return uno ;
  return due ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


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
  cout << "double * params[" << N << "] = {" ;
  for (int i = 0 ; i < N - 1 ; ++i)
    cout << array[i] << ", " ;
  cout << array[N-1] << "} ;\n" ;
  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void setParNamesdoubleGausCrystalBallLowHigh (TF1 * func)
{
  func->SetParName (1, "mean") ;
  func->SetParName (2, "sigma") ;
  func->SetParName (3, "alphaR") ;
  func->SetParName (4, "nL") ;
  func->SetParName (5, "alphaL") ;
  func->SetParName (6, "nR") ;
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
double diff_crystalBallLowHigh (double* x, double* par)
{
  return crystalBallLowHigh (x, par) - crystalBallLowHigh (x, par + 7) ;
  //PG    signal and interference          signal only
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


TF1 * FIT_madgraph_signal (TH1F * h_MWW_mg, float mass, float rangeScale, TString suffix, bool useLikelihood = false)
{  
  TCanvas * c4_mg = new TCanvas ("c4_mg", "c4_mg") ;

  //PG first fit: get the width for the second fit
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * gauss_mg = new TF1 ("gauss_mg", "gaus", 0, 2000) ;
  gauss_mg->SetNpx (10000) ;
  gauss_mg->SetLineWidth (1) ;
  gauss_mg->SetLineColor (kGray + 2) ;
  h_MWW_mg->Fit ("gauss_mg", "Q+", "", mass - 1. * h_MWW_mg->GetRMS () , mass + 1. * h_MWW_mg->GetRMS ()) ;
  
  //PG second fit: first with chisq, if requested with likelihood also
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * func_mg_1 = new TF1 ("func_mg_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_mg_1->SetNpx (10000) ;
  func_mg_1->SetLineWidth (1) ;
  func_mg_1->SetLineColor (kBlue + 1) ;

  setParNamesdoubleGausCrystalBallLowHigh (func_mg_1) ;

  func_mg_1->SetParameter (0, h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ())) ;  // multiplicative scale
  func_mg_1->SetParameter (1, mass) ;                                                  // mean
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
  cout << "-------------------\nFITTING THE MADGRAPH SIGNAL\n\n-------------------\n" ;
  TFitResultPtr result = h_MWW_mg->Fit ("func_mg_1", "+", "", 0.5 * mass + sign * 50, 2 * mass) ;
  if (useLikelihood && mass < 800)
    {
      cout << "-------------------\nFITTING THE MADGRAPH SIGNAL W/ LIKELIHOOD\n\n-------------------\n" ;
      func_mg_1->SetParameters (func_mg_1->GetParameters ()) ; //PG not sure this is necessary
      func_mg_1->SetLineColor (kBlue + 3) ;
      h_MWW_mg->Fit ("func_mg_1", "+L", "", 0.5 * mass - 50, 2 * mass) ;
    }

//  TH1F * h_MWW_mg_error = h_MWW_mg->Clone ("h_MWW_mg_error") ;
//  (TVirtualFitter::GetFitter ())->GetConfidenceIntervals (h_MWW_mg_error, 0.68) ;
//  h_MWW_mg_error->SetMarkerSize (0) ;
//  h_MWW_mg_error->SetFillColor (kAzure - 9) ;
//  h_MWW_mg_error->SetFillStyle (3001) ;

  //PG plotting the result
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  ymin = h_MWW_mg->GetBinContent (h_MWW_mg->GetMinimumBin ()) ;
  if (ymin <= 0) ymin = ymax / 10000. ; 
  TH1F * c4_mg_frame = (TH1F *) c4_mg->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_mg->SetLogy () ;
  c4_mg_frame->SetTitle (0) ;
  c4_mg_frame->SetStats (0) ;
  c4_mg_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
//  h_MWW_mg_error->Draw ("same") ;
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


TF1 * FIT_phantom_signal (TH1F * diff, float mass, float rangeScale, TString suffix, bool useLikelihood = false)
{  
 TCanvas * c4_ph = new TCanvas ("c4_ph", "c4_ph") ;

  //PG first fit: get the width for the second fit
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * gauss_ph = new TF1 ("gauss_ph", "gaus", 0, 2000) ;
  gauss_ph->SetNpx (10000) ;
  gauss_ph->SetLineWidth (1) ;
  gauss_ph->SetLineColor (kGray + 2) ;
  double fact = 1. ;
  diff->Fit ("gauss_ph", "Q+", "", mass - fact * diff->GetRMS () , mass + fact * diff->GetRMS ()) ;

  //PG second fit: first with chisq, if requested with likelihood also
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
  TF1 * func_ph_1 = new TF1 ("func_ph_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_ph_1->SetNpx (10000) ;
  func_ph_1->SetLineWidth (1) ;
  func_ph_1->SetLineColor (kRed + 1) ;
  
  setParNamesdoubleGausCrystalBallLowHigh (func_ph_1) ;

  func_ph_1->SetParameter (0, diff->Integral ()) ;              // multiplicative scale
  func_ph_1->SetParameter (1, mass) ;                           // mean
  func_ph_1->SetParameter (2, 2 * gauss_ph->GetParameter (2)) ; // gaussian sigma
  func_ph_1->SetParLimits (2, 0., mass) ;
  func_ph_1->SetParameter (3, 1) ;                              // right junction
  func_ph_1->SetParLimits (3, 0.1, 5) ;                         // right junction
  func_ph_1->FixParameter (4, 3) ;                              // right power law order            //PG NB THIS IS FIXED
//  func_ph_1->SetParLimits (4, 0.1, 3) ;                         // left junction
  func_ph_1->SetParameter (5, 1) ;                              // left junction
  func_ph_1->FixParameter (5, 1) ;                              // left junction                    //PG NB THIS IS FIXED
//  func_ph_1->SetParLimits (5, 0.1, 5) ;                         // left junction
  func_ph_1->FixParameter (6, 3) ;                              // left power law order             //PG NB THIS IS FIXED
//  func_ph_1->SetParLimits (6, 0.1, 3) ;                         // left junction

  cout << "-------------------\nFITTING THE PHANTOM SIGNAL\n\n-------------------\n" ;
  diff->Fit ("func_ph_1", "", "", 0.5 * mass - 50, 2 * mass) ;
  if (useLikelihood)
    {
      cout << "-------------------\nFITTING THE PHANTOM SIGNAL W/ LIKELIHOOD\n\n-------------------\n" ;
      func_ph_1->SetParameters (func_ph_1->GetParameters ()) ;
      func_ph_1->SetLineColor (kRed + 3) ;
      diff->Fit ("func_ph_1", "+L", "", 0.5 * mass - 50, 1.5 * mass) ;
    }
    
  float ymax = diff->GetBinContent (diff->GetMaximumBin ()) ;
  ymax = max (ymax, func_ph_1->GetMaximum ()) ;
  float ymin = diff->GetBinContent (diff->GetMinimumBin ()) ;
  ymin = max (ymin, -0.1 * ymax) ;
  if (ymin <= 0) ymin = ymax / 10000. ; 
  TH1F * c4_ph_frame = (TH1F *) c4_ph->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_ph_frame->SetTitle (0) ;
  c4_ph_frame->SetStats (0) ;
  c4_ph_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  gauss_ph->Draw ("same") ;
  diff->Draw ("EPsame") ;

  double rightTh_ph = fabs (func_ph_1->GetParameter (3)) * fabs (func_ph_1->GetParameter (2)) + func_ph_1->GetParameter (1) ;
  double leftTh_ph  = -1 * fabs (func_ph_1->GetParameter (5)) * fabs (func_ph_1->GetParameter (2)) + func_ph_1->GetParameter (1) ;

  TLine * l_rightTh_ph = new TLine (rightTh_ph, 0.9 * ymin, rightTh_ph, 1.1 * ymax) ;
  l_rightTh_ph->SetLineColor (kRed) ;
  l_rightTh_ph->Draw ("same") ;
  TLine * l_leftTh_ph = new TLine (leftTh_ph, 0.9 * ymin, leftTh_ph, 1.1 * ymax) ;
  l_leftTh_ph->SetLineColor (kRed) ;
  l_leftTh_ph->Draw ("same") ;

  c4_ph->Print (TString ("signals_ph_lin") + suffix, "pdf") ;
  c4_ph->SetLogy () ;
  c4_ph->Print (TString ("signals_ph_log") + suffix, "pdf") ;

  return func_ph_1 ;
}


// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====


int macro_findInterferece (string filename, double mass, int initialRebin = 1)                                                        
{        

  TVirtualFitter::SetDefaultFitter ("Minuit2") ;
  gSystem->Load ("Functions.cc") ;
  gStyle->SetPadTopMargin (0.1) ;

  TFile * f = new TFile (filename.c_str ()) ;
  TH1F * h_MWW_phbkgsig = (TH1F *) f->Get ("h_MWW_phbkgsig") ;
  TH1F * h_MWW_phbkg    = (TH1F *) f->Get ("h_MWW_phbkg") ;
  TH1F * h_MWW_mg       = (TH1F *) f->Get ("h_MWW_mg") ;

  int reBin = initialRebin ;
  if (mass > 480) reBin *= 2 ;
  if (mass > 810) reBin *= 6 ;
  double rangeScale = 1.5 ; // for plotting purposes
  if (mass > 480) rangeScale = 2 ;
  
  int scaling = 0 ;
  
  cout << "rebinning factor: " << reBin << "\n" ;
  cout << "MG bins number before rebinning: " << h_MWW_mg->GetNbinsX () << "\n" ;
  
  h_MWW_phbkgsig->Rebin (reBin) ;
  h_MWW_phbkg   ->Rebin (reBin) ;
  h_MWW_mg      ->Rebin (reBin) ;

  cout << "MG bins number after rebinning: " << h_MWW_mg->GetNbinsX () << "\n" ;

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

  //PG (SBI - B) / S
  TH1F * ratio = (TH1F *) diff->Clone ("ratio") ;
  ratio->Divide (h_MWW_mg) ;

  //PG (SBI - B) - S
  TH1F * delta = (TH1F *) diff->Clone ("delta") ;
  delta->SetTitle ("") ;
  delta->Add (h_MWW_mg, -1) ;
  delta->SetLineColor (kViolet + 1) ;

  //PG ((SBI - B) - S) / S
  TH1F * relDiff = delta->Clone ("relDiff") ;
  relDiff->Divide (h_MWW_mg) ;

  //PG plotting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  TString suffix (filename) ;
  suffix.ReplaceAll ("findInterference.", "") ;
  suffix.ReplaceAll ("root", "pdf") ;
  
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

  c1_leg = new TLegend (0.6,0.8,0.9,0.95) ;
  c1_leg->SetFillStyle (0) ;
  c1_leg->SetBorderSize (0) ;
  c1_leg->SetTextFont (42) ;
  c1_leg->AddEntry (h_MWW_phbkg, "B","l") ;
  c1_leg->AddEntry (h_MWW_phbkgsig, "SBI","l") ;
  c1_leg->AddEntry (h_MWW_mg, "S","l") ;
  c1_leg->Draw () ;
  
  c1->Print (TString ("spectra") + suffix, "pdf") ;
  
  //PG S only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TF1 * func_mg_1 = FIT_madgraph_signal (h_MWW_mg, mass, rangeScale, suffix, true) ;

  //PG (SBI - B) only ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TF1 * func_ph_1 = FIT_phantom_signal (diff, mass, rangeScale, suffix, false) ;

  //PG (SBI - B) - S only   i.e   INTERFERENCE FIT 
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  cout << "-------------------\nFITTING THE INTERFERENCE\n\n-------------------\n" ;

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
  TF1 * f_difference = new TF1 ("f_difference", diff_crystalBallLowHigh, 0, 2000, 14) ;
  Double_t params_difference[14] ;
  func_ph_1->GetParameters (params_difference) ;
  func_mg_1->GetParameters (params_difference + 7) ;
  f_difference->SetParameters (params_difference) ;
  f_difference->SetNpx (10000) ;
  f_difference->SetLineWidth (2) ;
  f_difference->SetLineColor (kGray + 1) ;
  f_difference->Draw ("same") ;

  c3_leg->Draw () ;

  c3->Print (TString ("interf") + suffix, "pdf") ;

  //PG (SBI - B) / S only   i.e   CORRECTION FACTOR
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TCanvas * c5 = new TCanvas ("c5", "c5") ;
  ymax = ratio->GetBinContent (ratio->GetMaximumBin ()) ;
  ymin = ratio->GetBinContent (ratio->GetMinimumBin ()) ;
  if (ymin > 0) ymin *= 0.5 ;
  else          ymin *= 1.5 ;
//  if (mass > 800 && ymin < -3 * ymax) ymin = -2 * ymax ;
  if (ymin < 0) ymin = - 0.1 ;
  TH1F * c5_frame = (TH1F *) c5->DrawFrame (200, 0., rangeScale * mass, 6.) ;
  c5_frame->SetTitle (0) ;
  c5_frame->SetStats (0) ;
  c5_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  
  ratio->Draw ("hist same") ;
  //PG plot the result of the fits on top
  TF1 * f_ratio = new TF1 ("f_ratio", ratio_crystalBallLowHigh, 0, 2000, 14) ;
  f_ratio->SetParameters (params_difference) ;
  f_ratio->SetNpx (10000) ;
  f_ratio->SetLineWidth (2) ;
  f_ratio->SetLineColor (kBlue + 2) ;
  f_ratio->Draw ("same") ;
  
  c5->Print (TString ("corr_factor") + suffix, "pdf") ;

  //PG plotting S only, and (SBI - B) 
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

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

  cout << "scaling applied: " << scaling << endl ;

  return 0 ;


}                                                                                   
                                                                                    
