/*
r00t -q macro_findInterferece.C\(\"findInterference.350.root\",350\)
r00t -q macro_findInterferece.C\(\"findInterference.500.root\",500\)
r00t -q macro_findInterferece.C\(\"findInterference.650.root\",650\)
r00t -q macro_findInterferece.C\(\"findInterference.800.root\",800\)
r00t -q macro_findInterferece.C\(\"findInterference.1000.root\",1000\)
*/


/*** right-hand side tail fit ***/
double expOfPowerLaw (double* x, double* par)
{
  double mean   = par[1] ; // mean
  double gamma  = par[2] ; // power law at the exponent
  double xx = x[0] - mean ;

  return par [0] * exp (-1 * pow (xx, gamma)) ;
}


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

  double norm = 1. / (shift * shift + gamma) - 1 / ((shift + 2 * distance) * (shift + 2 * distance) + gamma) ;
  return scale * (1. / norm) * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;
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


int macro_signalFit (string filename, double mass)                                                        
{        

//  TVirtualFitter::SetDefaultFitter ("Minuit2") ;
  gSystem->Load ("Functions.cc") ;
  gStyle->SetPadTopMargin (0.1) ;

  TFile * f = new TFile (filename.c_str ()) ;
  TH1F * h_MWW_mg       = (TH1F *) f->Get ("h_MWW_mg") ;

  int reBin = 1 ;
  if (mass > 480) reBin = 2 ;
  if (mass > 810) reBin = 6 ;
  double rangeScale = 1.5 ;
  if (mass > 480) rangeScale = 2 ;
  
  int scaling = 1 ;
  
  h_MWW_mg      ->Rebin (reBin) ;

  normaliseToBinWidth (h_MWW_mg) ;

  double ymax= 2000. ;
  double ymin= 0. ;

  //PG graphics
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  h_MWW_mg->SetStats (0) ;
  h_MWW_mg->SetTitle ("") ;
  h_MWW_mg->SetLineColor (kBlue) ;
  h_MWW_mg->SetLineWidth (2) ;
  h_MWW_mg->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;

  //PG plotting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  TString suffix = "." ; 
  suffix += mass ;
  suffix += ".pdf" ;
  
  //PG double crystal ball 
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  TCanvas * c4_mg = new TCanvas ("c4_mg", "c4_mg") ;
  ymax = h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ()) ;
  ymin = h_MWW_mg->GetBinContent (h_MWW_mg->GetMinimumBin ()) ;
  if (ymin <= 0) ymin = ymax / 10000. ; 

  TF1 * func_mg_1 = new TF1 ("func_mg_1", crystalBallLowHigh, 0, 2000, 7) ;
  func_mg_1->SetNpx (10000) ;
  func_mg_1->SetLineWidth (1) ;
  func_mg_1->SetLineColor (kBlue + 1) ;

  setParNamesdoubleGausCrystalBallLowHigh (func_mg_1) ;
  cout << "PREFIT RMS : " << h_MWW_mg->GetRMS () << endl ;

  func_mg_1->SetParameter (0, h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ())) ;  // multiplicative scale
  func_mg_1->SetParameter (1, mass) ;                                                  // mean
  func_mg_1->SetParameter (2, h_MWW_mg->GetRMS ()) ;                                   // gaussian sigma
  func_mg_1->SetParLimits (2, 0.1 * h_MWW_mg->GetRMS (), 20 * h_MWW_mg->GetRMS ()) ;
  func_mg_1->SetParameter (3, 1.5) ;                                                   // right junction
//  func_mg_1->SetParLimits (3, 0.1, 5) ;                                              // right junction
  func_mg_1->SetParameter (4, 2) ;                                                     // right power law order
  func_mg_1->SetParameter (5, 0.95) ;                                                   // left junction
//  func_mg_1->FixParameter (5, 0.95) ;                                                   // left junction
//  func_mg_1->SeaParLimits (5, 0.1, 5) ;                                              // left junction
  func_mg_1->SetParameter (6, 2.38) ;                                                  // left power law order

  int sign = 1 ;
  if (mass < 400) sign = -2 ;

  cout << "-------------------\nFITTING W/ DOUBLE CB\n\n-------------------\n" ;
  h_MWW_mg->Fit ("func_mg_1", "", "", 0.5 * mass + sign * 50, 2 * mass) ;
  cout << "CHI2 / NDOF = " << func_mg_1->GetChisquare () /func_mg_1->GetNDF () << endl ;
  func_mg_1->SetParameters (func_mg_1->GetParameters ()) ;
  func_mg_1->SetLineColor (kBlue + 3) ;
  cout << "-------------------\nFITTING W/ DOUBLE CB W/ LIKELIHOOD\n\n-------------------\n" ;
  h_MWW_mg->Fit ("func_mg_1", "L", "", 0.5 * mass - 50, 2 * mass) ;
  cout << "CHI2 / NDOF = " << func_mg_1->GetChisquare () /func_mg_1->GetNDF () << endl ;


  double rightTh_1 = fabs (func_mg_1->GetParameter (3)) * func_mg_1->GetParameter (2) + func_mg_1->GetParameter (1) ;
  double leftTh_1  = -1 * fabs (func_mg_1->GetParameter (5)) * func_mg_1->GetParameter (2) + func_mg_1->GetParameter (1) ;

  TLine * l_rightTh_1 = new TLine (rightTh_1, 0.9 * ymin, rightTh_1, 1.1 * ymax) ;
  l_rightTh_1->SetLineColor (kBlue + 1) ;
  l_rightTh_1->SetLineStyle (2) ;
  TLine * l_leftTh_1 = new TLine (leftTh_1, 0.9 * ymin, leftTh_1, 1.1 * ymax) ;
  l_leftTh_1->SetLineColor (kBlue + 1) ;
  l_leftTh_1->SetLineStyle (2) ;

//  //PG super gaus cum cauda 
//  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//  
//  TF1 * func_mg_2 = new TF1 ("func_mg_2", doubleSuperGausCumCauda, 0, 2000, 7) ;
//  func_mg_2->SetNpx (10000) ;
//  func_mg_2->SetLineWidth (1) ;
//  func_mg_2->SetLineColor (kRed + 1) ;
//  func_mg_2->SetParName (0, "N") ;
//  func_mg_2->SetParName (1, "mean") ;
//  func_mg_2->SetParName (2, "Nsigma") ;
//  func_mg_2->SetParName (3, "alphaR") ;
//  func_mg_2->SetParName (4, "alphaL") ;
//
////  func_mg_2->SetParameter (0, h_MWW_mg->GetBinContent (h_MWW_mg->GetMaximumBin ())) ;                  // multiplicative scale
//  func_mg_2->SetParameter (0, 1.) ;                  // multiplicative scale
//  func_mg_2->SetParameter (1, mass) ;                        // mean
//  func_mg_2->FixParameter (2, func_mg_1->GetParameter (2)) ; // gaussian sigma
////  func_mg_2->SetParameter (2, h_MWW_mg->GetRMS ()) ;       // gaussian sigma
////  func_mg_2->SetParLimits (2, 0.5 * h_MWW_mg->GetRMS (), 20 * h_MWW_mg->GetRMS ()) ;
//  func_mg_2->FixParameter (3, func_mg_1->GetParameter (3)) ;                   // right junction
//  func_mg_2->FixParameter (4,  func_mg_1->GetParameter (5)) ;                   // left junction
////  func_mg_2->SetParameter (3, 2) ;                   // right junction
////  func_mg_2->SetParameter (4, 2) ;                   // left junction
//
//  int sign = 1 ;
//  if (mass < 400) sign = -2 ;
//  cout << "-------------------\nFITTING W/ SUPER GAUSS CC\n\n-------------------\n" ;
//  h_MWW_mg->Fit ("func_mg_2", "", "", 0.5 * mass + sign * 50, 2 * mass) ;
//  cout << "CHI2 / NDOF = " << func_mg_2->GetChisquare () /func_mg_2->GetNDF () << endl ;
//  func_mg_2->SetParameters (func_mg_2->GetParameters ()) ;
//  func_mg_2->FixParameter (3, func_mg_1->GetParameter (3)) ;                   // right junction
//  func_mg_2->FixParameter (4,  func_mg_1->GetParameter (5)) ;                   // left junction
//  func_mg_2->SetLineColor (kRed + 3) ;
//  cout << "-------------------\nFITTING W/ SUPER GAUSS CC W/ LIKELIHOOD\n\n-------------------\n" ;
//  h_MWW_mg->Fit ("func_mg_2", "L", "", 0.5 * mass - 50, 2 * mass) ;
//  cout << "CHI2 / NDOF = " << func_mg_2->GetChisquare () /func_mg_2->GetNDF () << endl ;
//
//  double rightTh_2 = fabs (func_mg_2->GetParameter (3)) * func_mg_2->GetParameter (2) + func_mg_2->GetParameter (1) ;
//  double leftTh_2  = -1 * fabs (func_mg_2->GetParameter (4)) * func_mg_2->GetParameter (2) + func_mg_2->GetParameter (1) ;
//
//  TLine * l_rightTh_2 = new TLine (rightTh_2, 0.9 * ymin, rightTh_2, 1.1 * ymax) ;
//  l_rightTh_2->SetLineColor (kRed + 1) ;
//  l_rightTh_2->SetLineStyle (2) ;
//  TLine * l_leftTh_2 = new TLine (leftTh_2, 0.9 * ymin, leftTh_2, 1.1 * ymax) ;
//  l_leftTh_2->SetLineColor (kRed + 1) ;
//  l_leftTh_2->SetLineStyle (2) ;
//
//
  //PG tail fitting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

//  cout << "-------------------\nFITTING THE TAIL W/ EXP OF POWER LAW\n\n-------------------\n" ;
//
//  double start_3 = rightTh_1 ; //+ 0.5 * func_mg_2->GetParameter (2) ;
//  TF1 * func_mg_3 = new TF1 ("func_mg_3", expOfPowerLaw, 0, 2000, 7) ;
//  func_mg_3->SetNpx (10000) ;
//  func_mg_3->SetLineWidth (1) ;
//  func_mg_3->SetLineColor (kOrange + 1) ;
//  func_mg_3->SetParName (0, "N") ;
//  func_mg_3->SetParName (1, "mean") ;
//  func_mg_3->SetParName (2, "gamma") ;
//  func_mg_3->SetParameter (0, h_MWW_mg->Integral (start_3, 2 * mass)) ;            // multiplicative scale
//  func_mg_3->FixParameter (1, mass) ;                        // mean
//  func_mg_3->SetParameter (2, 1) ; // gaussian sigma
//
//  h_MWW_mg->Fit ("func_mg_3", "L", "", start_3, 2 * mass) ;
//

  //PG plotting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  
  TH1F * c4_mg_frame = (TH1F *) c4_mg->DrawFrame (200, 0.9 * ymin, rangeScale * mass, 1.1 * ymax) ;
  c4_mg->SetLogy () ;
  c4_mg_frame->SetTitle (0) ;
  c4_mg_frame->SetStats (0) ;
  c4_mg_frame->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;

  h_MWW_mg->Draw ("EPsame") ;
  func_mg_1->Draw ("same") ;
//  func_mg_2->Draw ("same") ;
//  func_mg_3->Draw ("same") ;
  
  
//  char testo[10] ;
//  sprintf (testo, "exp pol. %.2f", func_mg_3->GetParameter (2)) ;

  leg = new TLegend (0.35,0.26,0.6,0.4) ;
  leg->SetTextFont (42) ;
  leg->SetBorderSize (0) ;
  leg->SetFillStyle (0) ;
  leg->AddEntry (func_mg_1, "double CB", "l") ;
//  leg->AddEntry (func_mg_2, "exp tails", "l") ;
//  leg->AddEntry (func_mg_3, testo, "l") ;
//  leg->AddEntry (l_rightTh_2, "boundaries", "l") ;
  leg->Draw () ;

  l_rightTh_1->Draw ("same") ;
  l_leftTh_1->Draw ("same") ;

//  l_rightTh_2->Draw ("same") ;
//  l_leftTh_2->Draw ("same") ;

  c4_mg->Update () ;
  c4_mg->Print (TString ("fitSignals_mg") + suffix, "pdf") ;



  //PG output of the fitting function parameters

  std::ofstream outfile;

  outfile.open ("graphs_sig.txt", std::ios_base::app) ;
  outfile << "  \n  // ----> MASS " << mass << " ---- ---- ---- \n\n" ;

  outfile << "  // MG signal only parametrisation:\n" ;
  outfile << "  tg_sig_par0->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (0) << ") ;\n" ;
  outfile << "  tg_sig_par1->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (1) << ") ;\n" ;
  outfile << "  tg_sig_par2->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (2) << ") ;\n" ;
  outfile << "  tg_sig_par3->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (3) << ") ;\n" ;
  outfile << "  tg_sig_par4->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (4) << ") ;\n" ;
  outfile << "  tg_sig_par5->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (5) << ") ;\n" ;
  outfile << "  tg_sig_par6->SetPoint (i, " << mass << ", " << func_mg_1->GetParameter (6) << ") ;\n" ;
  outfile << "  TF1 * func_sig_" << mass << " = new TF1 (\"func_sig_" << mass << "\",crystalBallLowHigh, 200, 2000, 7) ;\n" ; 
  outfile << "  double params_sig_" << mass << "[7] = {" << func_mg_1->GetParameter (0) << ", " << func_mg_1->GetParameter (1) << ", " << func_mg_1->GetParameter (2) << ", " << func_mg_1->GetParameter (3) << ", " << func_mg_1->GetParameter (4) << ", " << func_mg_1->GetParameter (5) << ", " << func_mg_1->GetParameter (6)  << " } ;\n" ;
  outfile << "  func_sig_" << mass << "->SetParameters (params_sig_" << mass << ") ;\n" ; 
  outfile << "  i++ ;\n" ;

  outfile.close () ;


  return 0 ;


}                                                                                   
                                                                                    
