// r00t -q graphs.C > fits.txt


/*** double crystall ball ***/
double crystalBallLowHigh (double* x, double* par)
{
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2
  
  double xx = x[0];
  double mean = par[1];
  double sigma = par[2];
  double alpha = par[3];
  double n = par[4];
  double alpha2 = par[5];
  double n2 = par[6];

  if( (xx-mean)/sigma > fabs(alpha) )
  {
    double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
    double B = n/fabs(alpha) - fabs(alpha);
    
    return par[0] * A * pow(B + (xx-mean)/sigma, -1.*n);
  }
  
  else if( (xx-mean)/sigma < -1.*fabs(alpha2) )
  {
    double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
    double B = n2/fabs(alpha2) - fabs(alpha2);
    
    return par[0] * A * pow(B - (xx-mean)/sigma, -1.*n2);
  }
  
  else
  {
    return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
  }
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


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


Double_t min (Double_t uno, Double_t due)
{
  if (uno < due) return uno ;
  return due ;
} 


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t max (Double_t uno, Double_t due)
{
  if (uno > due) return uno ;
  return due ;
} 


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


Double_t intereferenceOverSignal (Double_t * xx, Double_t * par) 
{
  if (crystalBallLowHigh (xx, par + 4) < 0.0000001) return 0 ;
  return doublePeakModel (xx, par) / crystalBallLowHigh (xx, par + 4) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int graphs () 
{
  TGraph * tg_par0 = new TGraph (5) ;
  TGraph * tg_par1 = new TGraph (5) ;
  TGraph * tg_par2 = new TGraph (5) ;
  TGraph * tg_par3 = new TGraph (5) ;
  
  TGraph * tg_sig_par0 = new TGraph (5) ;
  TGraph * tg_sig_par1 = new TGraph (5) ;
  TGraph * tg_sig_par2 = new TGraph (5) ;
  TGraph * tg_sig_par3 = new TGraph (5) ;
  TGraph * tg_sig_par4 = new TGraph (5) ;
  TGraph * tg_sig_par5 = new TGraph (5) ;
  TGraph * tg_sig_par6 = new TGraph (5) ;

  TGraph * tg_sAi_par0 = new TGraph (5) ;
  TGraph * tg_sAi_par1 = new TGraph (5) ;
  TGraph * tg_sAi_par2 = new TGraph (5) ;
  TGraph * tg_sAi_par3 = new TGraph (5) ;
  TGraph * tg_sAi_par4 = new TGraph (5) ;
  TGraph * tg_sAi_par5 = new TGraph (5) ;
  TGraph * tg_sAi_par6 = new TGraph (5) ;

  int i = 0 ; 
  
  
  //PG here follow the fitting function parameters, inserted from graph.txt
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  
  
  
  // ----> MASS 350 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 350, -3.71341e-05) ;
  tg_par1->SetPoint (i, 350, 353.091) ;
  tg_par2->SetPoint (i, 350, 0.0008) ;
  tg_par3->SetPoint (i, 350, 1074.59) ;
  TF1 * func_350 = new TF1 ("func_350",doublePeakModel, 200, 2000, 4) ;
  double params_350[4] = {-3.71341e-05, 353.091, 0.0008, 1074.59 } ;
  func_350->SetParameters (params_350) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 350, 0.000952886) ;
  tg_sig_par1->SetPoint (i, 350, 350.384) ;
  tg_sig_par2->SetPoint (i, 350, 7.70264) ;
  tg_sig_par3->SetPoint (i, 350, 1.16314) ;
  tg_sig_par4->SetPoint (i, 350, 1.53064) ;
  tg_sig_par5->SetPoint (i, 350, 1.27642) ;
  tg_sig_par6->SetPoint (i, 350, 2.17316) ;
  TF1 * func_sig_350 = new TF1 ("func_sig_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_350[7] = {0.000952886, 350.384, 7.70264, 1.16314, 1.53064, 1.27642, 2.17316 } ;
  func_sig_350->SetParameters (params_sig_350) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 350, 0.000978532) ;
  tg_sAi_par1->SetPoint (i, 350, 349.72) ;
  tg_sAi_par2->SetPoint (i, 350, -7.33274) ;
  tg_sAi_par3->SetPoint (i, 350, 1.03633) ;
  tg_sAi_par4->SetPoint (i, 350, 3) ;
  tg_sAi_par5->SetPoint (i, 350, 1) ;
  tg_sAi_par6->SetPoint (i, 350, 3) ;
  TF1 * func_sAi_350 = new TF1 ("func_sAi_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_350[7] = {0.000978532, 349.72, -7.33274, 1.03633, 3, 1, 3 } ;
  func_sAi_350->SetParameters (params_sAi_350) ;
  i++ ;
  
  // ----> MASS 500 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 500, -9.23198e-06) ;
  tg_par1->SetPoint (i, 500, 504.435) ;
  tg_par2->SetPoint (i, 500, 0.0008) ;
  tg_par3->SetPoint (i, 500, 11291.6) ;
  TF1 * func_500 = new TF1 ("func_500",doublePeakModel, 200, 2000, 4) ;
  double params_500[4] = {-9.23198e-06, 504.435, 0.0008, 11291.6 } ;
  func_500->SetParameters (params_500) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 500, 8.37317e-05) ;
  tg_sig_par1->SetPoint (i, 500, 504.033) ;
  tg_sig_par2->SetPoint (i, 500, 28.0345) ;
  tg_sig_par3->SetPoint (i, 500, 0.814726) ;
  tg_sig_par4->SetPoint (i, 500, 2.42493) ;
  tg_sig_par5->SetPoint (i, 500, 1.01086) ;
  tg_sig_par6->SetPoint (i, 500, 4.83438) ;
  TF1 * func_sig_500 = new TF1 ("func_sig_500",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_500[7] = {8.37317e-05, 504.033, 28.0345, 0.814726, 2.42493, 1.01086, 4.83438 } ;
  func_sig_500->SetParameters (params_sig_500) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 500, 8.21694e-05) ;
  tg_sAi_par1->SetPoint (i, 500, 498.063) ;
  tg_sAi_par2->SetPoint (i, 500, -30.2295) ;
  tg_sAi_par3->SetPoint (i, 500, 1.22531) ;
  tg_sAi_par4->SetPoint (i, 500, 3) ;
  tg_sAi_par5->SetPoint (i, 500, 1) ;
  tg_sAi_par6->SetPoint (i, 500, 3) ;
  TF1 * func_sAi_500 = new TF1 ("func_sAi_500",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_500[7] = {8.21694e-05, 498.063, -30.2295, 1.22531, 3, 1, 3 } ;
  func_sAi_500->SetParameters (params_sAi_500) ;
  i++ ;
  
  // ----> MASS 650 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 650, -3.42623e-06) ;
  tg_par1->SetPoint (i, 650, 660.125) ;
  tg_par2->SetPoint (i, 650, 0.0008) ;
  tg_par3->SetPoint (i, 650, 21723.7) ;
  TF1 * func_650 = new TF1 ("func_650",doublePeakModel, 200, 2000, 4) ;
  double params_650[4] = {-3.42623e-06, 660.125, 0.0008, 21723.7 } ;
  func_650->SetParameters (params_650) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 650, 1.76845e-05) ;
  tg_sig_par1->SetPoint (i, 650, 661.821) ;
  tg_sig_par2->SetPoint (i, 650, 62.0886) ;
  tg_sig_par3->SetPoint (i, 650, 0.792082) ;
  tg_sig_par4->SetPoint (i, 650, 6.19973) ;
  tg_sig_par5->SetPoint (i, 650, 0.954558) ;
  tg_sig_par6->SetPoint (i, 650, 16.8303) ;
  TF1 * func_sig_650 = new TF1 ("func_sig_650",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_650[7] = {1.76845e-05, 661.821, 62.0886, 0.792082, 6.19973, 0.954558, 16.8303 } ;
  func_sig_650->SetParameters (params_sig_650) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 650, 1.75996e-05) ;
  tg_sAi_par1->SetPoint (i, 650, 641.578) ;
  tg_sAi_par2->SetPoint (i, 650, -65.6889) ;
  tg_sAi_par3->SetPoint (i, 650, 1.6513) ;
  tg_sAi_par4->SetPoint (i, 650, 3) ;
  tg_sAi_par5->SetPoint (i, 650, 1) ;
  tg_sAi_par6->SetPoint (i, 650, 3) ;
  TF1 * func_sAi_650 = new TF1 ("func_sAi_650",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_650[7] = {1.75996e-05, 641.578, -65.6889, 1.6513, 3, 1, 3 } ;
  func_sAi_650->SetParameters (params_sAi_650) ;
  i++ ;
  
  // ----> MASS 800 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 800, -1.21044e-06) ;
  tg_par1->SetPoint (i, 800, 805.215) ;
  tg_par2->SetPoint (i, 800, 0.0008) ;
  tg_par3->SetPoint (i, 800, 31300.9) ;
  TF1 * func_800 = new TF1 ("func_800",doublePeakModel, 200, 2000, 4) ;
  double params_800[4] = {-1.21044e-06, 805.215, 0.0008, 31300.9 } ;
  func_800->SetParameters (params_800) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 800, 4.2978e-06) ;
  tg_sig_par1->SetPoint (i, 800, 806.433) ;
  tg_sig_par2->SetPoint (i, 800, 108.112) ;
  tg_sig_par3->SetPoint (i, 800, 1.02477) ;
  tg_sig_par4->SetPoint (i, 800, 8.1546) ;
  tg_sig_par5->SetPoint (i, 800, 0.977263) ;
  tg_sig_par6->SetPoint (i, 800, 141.628) ;
  TF1 * func_sig_800 = new TF1 ("func_sig_800",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_800[7] = {4.2978e-06, 806.433, 108.112, 1.02477, 8.1546, 0.977263, 141.628 } ;
  func_sig_800->SetParameters (params_sig_800) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 800, 4.90625e-06) ;
  tg_sAi_par1->SetPoint (i, 800, 743.86) ;
  tg_sAi_par2->SetPoint (i, 800, -116.156) ;
  tg_sAi_par3->SetPoint (i, 800, 4.85572) ;
  tg_sAi_par4->SetPoint (i, 800, 3) ;
  tg_sAi_par5->SetPoint (i, 800, 1) ;
  tg_sAi_par6->SetPoint (i, 800, 3) ;
  TF1 * func_sAi_800 = new TF1 ("func_sAi_800",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_800[7] = {4.90625e-06, 743.86, -116.156, 4.85572, 3, 1, 3 } ;
  func_sAi_800->SetParameters (params_sAi_800) ;
  i++ ;
  
  // ----> MASS 1000 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 1000, -1.97628e-07) ;
  tg_par1->SetPoint (i, 1000, 954.874) ;
  tg_par2->SetPoint (i, 1000, 0.0008) ;
  tg_par3->SetPoint (i, 1000, 66307) ;
  TF1 * func_1000 = new TF1 ("func_1000",doublePeakModel, 200, 2000, 4) ;
  double params_1000[4] = {-1.97628e-07, 954.874, 0.0008, 66307 } ;
  func_1000->SetParameters (params_1000) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 1000, 5.40097e-07) ;
  tg_sig_par1->SetPoint (i, 1000, 934.543) ;
  tg_sig_par2->SetPoint (i, 1000, 223.329) ;
  tg_sig_par3->SetPoint (i, 1000, 1.29696) ;
  tg_sig_par4->SetPoint (i, 1000, 94.1604) ;
  tg_sig_par5->SetPoint (i, 1000, -252619) ;
  tg_sig_par6->SetPoint (i, 1000, 75.0455) ;
  TF1 * func_sig_1000 = new TF1 ("func_sig_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_1000[7] = {5.40097e-07, 934.543, 223.329, 1.29696, 94.1604, -252619, 75.0455 } ;
  func_sig_1000->SetParameters (params_sig_1000) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 1000, 9.86385e-07) ;
  tg_sAi_par1->SetPoint (i, 1000, 745.657) ;
  tg_sAi_par2->SetPoint (i, 1000, -226.263) ;
  tg_sAi_par3->SetPoint (i, 1000, 4.94692) ;
  tg_sAi_par4->SetPoint (i, 1000, 3) ;
  tg_sAi_par5->SetPoint (i, 1000, 1) ;
  tg_sAi_par6->SetPoint (i, 1000, 3) ;
  TF1 * func_sAi_1000 = new TF1 ("func_sAi_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_1000[7] = {9.86385e-07, 745.657, -226.263, 4.94692, 3, 1, 3 } ;
  func_sAi_1000->SetParameters (params_sAi_1000) ;
  i++ ;


  //PG plotting
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  //PG plot and fit the trend of each single parameter of the interference functions
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  // interference
  tg_par0->SetTitle ("") ;
  tg_par1->SetTitle ("") ;
  tg_par2->SetTitle ("") ;
  tg_par3->SetTitle ("") ;
  
  tg_par0->GetYaxis ()->SetTitle ("scale") ;
  tg_par1->GetYaxis ()->SetTitle ("shift") ;
  tg_par2->GetYaxis ()->SetTitle ("distance") ;
  tg_par3->GetYaxis ()->SetTitle ("gamma") ;
  
  tg_par0->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_par1->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_par2->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_par3->GetXaxis ()->SetTitle ("m_{WW}") ;

  // madgrapgh signal
  tg_sig_par0->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par1->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par2->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par3->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par4->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par5->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sig_par6->GetXaxis ()->SetTitle ("m_{WW}") ;

  tg_sig_par0->GetYaxis ()->SetTitle ("scale") ;
  tg_sig_par1->GetYaxis ()->SetTitle ("mean") ;
  tg_sig_par2->GetYaxis ()->SetTitle ("gaussian sigma") ;
  tg_sig_par3->GetYaxis ()->SetTitle ("right alpha") ;
  tg_sig_par4->GetYaxis ()->SetTitle ("right power law") ;
  tg_sig_par5->GetYaxis ()->SetTitle ("left alpha") ;
  tg_sig_par6->GetYaxis ()->SetTitle ("left power law") ;

  //phantom signal + interference
  tg_sAi_par0->SetTitle ("") ;
  tg_sAi_par1->SetTitle ("") ;
  tg_sAi_par2->SetTitle ("") ;
  tg_sAi_par3->SetTitle ("") ;
  tg_sAi_par4->SetTitle ("") ;
  tg_sAi_par5->SetTitle ("") ;
  tg_sAi_par6->SetTitle ("") ;

  tg_sAi_par0->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par1->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par2->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par3->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par4->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par5->GetXaxis ()->SetTitle ("m_{WW}") ;
  tg_sAi_par6->GetXaxis ()->SetTitle ("m_{WW}") ;

  tg_sAi_par0->GetYaxis ()->SetTitle ("scale") ;
  tg_sAi_par1->GetYaxis ()->SetTitle ("mean") ;
  tg_sAi_par2->GetYaxis ()->SetTitle ("gaussian sAima") ;
  tg_sAi_par3->GetYaxis ()->SetTitle ("right alpha") ;
  tg_sAi_par4->GetYaxis ()->SetTitle ("right power law") ;
  tg_sAi_par5->GetYaxis ()->SetTitle ("left alpha") ;
  tg_sAi_par6->GetYaxis ()->SetTitle ("left power law") ;

  //PG PARAMETERS OF THE INTERFERENCE
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  TCanvas * c_par = new TCanvas () ;

  c_par->Divide (2,2) ;
  i = 0 ;
  double x, y ;
  
  cout << "---> INTERF PARAM 0\n" ;
  c_par->cd (++i) ; 

//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + [1] * TMath::Power (1 - TMath::Exp (- [2] * (x * x - [3] * [3])), 2)") ;
//  fit_param0->FixParameter (0, -0.000005) ; // vertical offset
//  fit_param0->FixParameter (1, 0.0000033) ; // difference betw. min and plateau
//  fit_param0->FixParameter (2, 0.00001) ;
//  fit_param0->FixParameter (3, 460) ; // position of the minimum
//  fit_param0->SetLineWidth (1) ;
//  fit_param0->SetLineColor (kBlue + 1) ;
//  tg_par0->Fit (fit_param0, "", "", 0, 700) ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + [1] * x + [2] * x * x + [3] * x * x * x + [4] * x * x * x * x") ;
//  fit_param0->SetParameter (0,   0.00013466) ;
//  fit_param0->SetParameter (1, -8.60154e-07) ;
//  fit_param0->SetParameter (2,  1.89307e-09) ;
//  fit_param0->SetParameter (3,  -1.7511e-12) ;
//  fit_param0->SetParameter (4,  5.79317e-16) ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + [1] * exp (x * [2])") ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + [1] * exp ([2] * (x - [3]))") ;
//  fit_param0->SetParameter (0,   0.) ;
//  fit_param0->SetParameter (1, -0.00005) ;
//  fit_param0->SetParameter (1, -0.00005) ;
//  fit_param0->SetParameter (2,  -1) ;
//  fit_param0->SetParameter (3,  350) ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + [1] * x + [2] * x * x + [3] * x * x * x") ;
//  fit_param0->SetParameter (0,   0.00013466) ;
//  fit_param0->SetParameter (1, -8.60154e-07) ;
//  fit_param0->SetParameter (2,  1.89307e-09) ;
//  fit_param0->SetParameter (3,  -1.7511e-12) ;
//  fit_param0->SetParameter (4,  5.79317e-16) ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + 1. / ([2] * TMath::Sqrt ( abs (x - [1])))") ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + 1. / ([2] * (x - [1]))") ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + 1. / ([2] * (x - [1]) + [3] * (x - [1]) * (x - [1]))") ;
//  fit_param0->SetParameter (0,   0.) ;
//  fit_param0->SetParameter (1, 300) ;
//  fit_param0->SetParameter (2, -10000.) ;
//  fit_param0->SetParameter (2, -1000.) ;
//  fit_param0->SetParameter (3, -10.) ;
//  fit_param0->SetParameter (4,  5.79317e-16) ;
//  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + 1 / ([1] * (x - [1]) * (x - [1]) * (x - [1]))") ;

//  fit_param0->SetLineWidth (1) ;
//  fit_param0->SetLineColor (kBlue + 1) ;

//  tg_par0->Fit (fit_param0) ;

  tg_par0->Draw ("AL*") ; 

  cout << "---> INTERF PARAM 1\n" ;
  c_par->cd (++i) ; tg_par1->Draw ("AL*") ; 
  TF1 * fit_param1 = new TF1 ("fit_param1", "[0] + [1] * x + [2] * x * x") ;
  fit_param1->SetLineWidth (1) ;
  fit_param1->SetLineColor (kBlue + 1) ;
  tg_par1->Fit (fit_param1) ;
  cout << "---> INTERF PARAM 2\n" ;
  tg_par2->GetPoint (1, x, y) ;
  cout << y << endl ;
  c_par->cd (++i) ; tg_par2->Draw ("AL*") ; 
  
  cout << "---> INTERF PARAM 3\n" ;
  c_par->cd (++i) ; tg_par3->Draw ("AL*") ; 
//  TF1 * sigmoid = new TF1 ("sigmoid", "[2] + [0] / (1.0 + TMath::Exp(-[1] * (x - [3])))") ;
//  sigmoid->SetParameter (0, 25000) ;
//  sigmoid->SetParameter (1, 0.005) ;
//  sigmoid->SetParameter (2, 30000) ;
//  sigmoid->SetParameter (3, 700) ;
//  sigmoid->SetLineWidth (1) ;
//  sigmoid->SetLineColor (kBlue + 1) ;
  TF1 * fit_param3 = new TF1 ("fit_param3", "[0] + [1] * x + [2] * x * x") ;
  fit_param3->SetLineWidth (1) ;
  fit_param3->SetLineColor (kBlue + 1) ;
  tg_par3->Fit (fit_param3) ;
 
  c_par->Print ("params_interf.pdf", "pdf") ;

  //PG PARAMETERS OF THE MADGRAPH SIGNAL
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


  //PG get the initial params of the first param fit
  TGraph * log_tg_sig_par0 = new TGraph (5) ;
  for (int i = 0 ; i < 5 ; ++i) 
    {
      double x, y ;
      tg_sig_par0->GetPoint (i, x, y) ;
      log_tg_sig_par0->SetPoint (i, x, TMath::Log (y)) ;
    }
  TF1 * mypol3 = new TF1 ("mypol3", "[0] + [2] * (x - [1]) + [3] * (x - [1]) * (x - [1]) + [4] * (x - [1]) * (x - [1]) * (x - [1])", 200, 2000) ;
  cout << "-------------------\nmypol2\n\n-------------------\n" ;
  log_tg_sig_par0->Fit ("mypol3") ;
//  TCanvas * tempo_log = new TCanvas ("tempo_log", "tempo_log", 600, 600) ;
//  log_tg_sig_par0->Draw ("AL*") ;


//  TCanvas * tempo = new TCanvas ("tempo", "tempo", 600, 600) ;
//  tempo->DrawFrame (300, 0.000001, 1100, 0.0012) ;
//  tempo->SetLogy () ;
//  tg_sig_par0->Draw ("L*") ; 
//  cout << "-------------------\nexpo\n\n-------------------\n" ;
//  tg_sig_par0->Fit ("expo") ;
//  tg_sig_par0->GetFunction ("expo")->SetLineWidth (1) ;
//  tg_sig_par0->GetFunction ("expo")->SetLineColor (kBlue + 1) ;
//  TF1 * expo_1 = new TF1 ("expo_1","TMath::Exp ([0] + [2] * (x - [1]))", 200, 2000) ;
//  expo_1->SetLineColor (kRed) ;
//  expo_1->SetLineWidth (1) ;
//  expo_1->SetParameter (0, tg_sig_par0->GetFunction ("expo")->GetParameter (0)) ;
//  expo_1->SetParameter (1, 0.) ;
//  expo_1->SetParameter (2, tg_sig_par0->GetFunction ("expo")->GetParameter (1)) ;
//  cout << "-------------------\nexpo_1\n\n-------------------\n" ;
//  tg_sig_par0->Fit (expo_1, "+") ;
//  TF1 * expo_2 = new TF1 ("expo_2","TMath::Exp ([0] + [2] * (x - [1]) + [3] * (x - [1]) * (x - [1]))", 200, 2000) ;
//  expo_2->SetLineColor (kViolet + 1) ;
//  expo_2->SetLineWidth (1) ;
//  expo_2->SetParameter (0, mypol2->GetParameter (0)) ;
//  expo_2->SetParameter (1, mypol2->GetParameter (1)) ;
//  expo_2->SetParameter (2, mypol2->GetParameter (2)) ;
//  expo_2->SetParameter (3, mypol2->GetParameter (3)) ;
//  cout << "-------------------\nexpo_2\n\n-------------------\n" ;
//  tg_sig_par0->Fit (expo_2, "+") ;
//  TF1 * expo_3 = new TF1 ("expo_3","TMath::Exp ([0] + [2] * (x - [1]) + [3] * (x - [1]) * (x - [1]) + [4] * (x - [1]) * (x - [1]) * (x - [1]))", 200, 2000) ;
//  expo_3->SetLineColor (kViolet + 1) ;
//  expo_3->SetLineWidth (1) ;
//  expo_3->SetParameter (0, mypol3->GetParameter (0)) ;
//  expo_3->SetParameter (1, mypol3->GetParameter (1)) ;
//  expo_3->SetParameter (2, mypol3->GetParameter (2)) ;
//  expo_3->SetParameter (3, mypol3->GetParameter (3)) ;
//  expo_3->SetParameter (4, mypol3->GetParameter (4)) ;
//  cout << "-------------------\nexpo_3\n\n-------------------\n" ;
//  tg_sig_par0->Fit (expo_3, "+") ;
//  TF1 * powlaw_1 = new TF1 ("powlaw_1","[0] + [1] / (x - [2])", 200, 2000) ;
//  powlaw_1->SetLineColor (kGreen + 2) ;
//  powlaw_1->SetLineWidth (1) ;
//  powlaw_1->SetParameter (0, 0.) ;
//  powlaw_1->SetParameter (1, 0.05) ;
//  powlaw_1->SetParameter (2, 300) ;
//  cout << "-------------------\npowlaw_1\n\n-------------------\n" ;
//  tg_sig_par0->Fit (powlaw_1, "+") ;


  TCanvas * c_sig_par = new TCanvas ("c_sig_par", "c_sig_par", 4000, 600) ;
  c_sig_par->Divide (4,2) ;
  i = 0 ;
  cout << "---> SIGNAL PARAM 0\n" ;
  c_sig_par->cd (++i) ; tg_sig_par0->Draw ("AL*") ; 
  TF1 * expo_3 = new TF1 ("expo_3","TMath::Exp ([0] + [2] * (x - [1]) + [3] * (x - [1]) * (x - [1]) + [4] * (x - [1]) * (x - [1]) * (x - [1]))", 200, 2000) ;
  expo_3->SetLineColor (kBlue + 1) ;
  expo_3->SetLineWidth (1) ;
  expo_3->SetParameter (0, mypol3->GetParameter (0)) ;
  expo_3->SetParameter (1, mypol3->GetParameter (1)) ;
  expo_3->SetParameter (2, mypol3->GetParameter (2)) ;
  expo_3->SetParameter (3, mypol3->GetParameter (3)) ;
  expo_3->SetParameter (4, mypol3->GetParameter (4)) ;
  cout << "-------------------\nexpo_3\n\n-------------------\n" ;
  tg_sig_par0->Fit (expo_3, "+") ;
  cout << "---> SIGNAL PARAM 1\n" ;
  c_sig_par->cd (++i) ; tg_sig_par1->Draw ("AL*") ; 
  tg_sig_par1->Fit ("pol2") ;
  tg_sig_par1->GetFunction ("pol2")->SetLineWidth (1) ;
  tg_sig_par1->GetFunction ("pol2")->SetLineColor (kBlue + 1) ;
  cout << "---> SIGNAL PARAM 2\n" ;
  c_sig_par->cd (++i) ; tg_sig_par2->Draw ("AL*") ;
  tg_sig_par2->Fit ("pol2") ;
  tg_sig_par2->GetFunction ("pol2")->SetLineWidth (1) ;
  tg_sig_par2->GetFunction ("pol2")->SetLineColor (kBlue + 1) ;
  c_sig_par->cd (++i) ; tg_sig_par3->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par4->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par5->Draw ("AL*") ;
  c_sig_par->cd (++i) ; tg_sig_par6->Draw ("AL*") ;

  c_sig_par->Print ("params_signal.pdf", "pdf") ;
  
  //PG PARAMETERS OF THE PHANTOM SIGNAL + INTERFERENCE
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  TCanvas * c_sAi_par = new TCanvas ("c_sAi_par", "c_sAi_par", 4000, 600) ;
  c_sAi_par->Divide (4,2) ;
  i = 0 ;
  c_sAi_par->cd (++i) ; tg_sAi_par0->SetLineColor (kRed) ; tg_sAi_par0->Draw ("AL*") ; 
  c_sAi_par->cd (++i) ; tg_sAi_par1->SetLineColor (kRed) ; tg_sAi_par1->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par2->SetLineColor (kRed) ; tg_sAi_par2->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par3->SetLineColor (kRed) ; tg_sAi_par3->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par4->SetLineColor (kRed) ; tg_sAi_par4->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par5->SetLineColor (kRed) ; tg_sAi_par5->Draw ("AL*") ;
  c_sAi_par->cd (++i) ; tg_sAi_par6->SetLineColor (kRed) ; tg_sAi_par6->Draw ("AL*") ;

  c_sAi_par->Print ("params_phantom.pdf", "pdf") ;

  //PG plot the interference functions
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  // f_IOS_MASS    = interfrence over signal
  // func_MASS     = interference
  // func_sig_MASS = madgraph signal
  // func_sAi_MASS = phantom signal and intereference

  TCanvas * c_results = new TCanvas ("c_results", "c_results", 5000, 600) ;
  c_results->Divide (5,2) ;

  i = 1 ;
  double mass = 350 ;
  double rangeScale = 1.5 ;
  TH1F * bkg ;

  // 350 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 350 ;
  rangeScale = 1.5 ;
  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_350_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_350_pars[j] = params_350[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_350_pars[j + 4] = params_sig_350[j] ;
  TF1 * f_IOS_350 = new TF1 ("f_IOS_350", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_350->SetParameters (IOS_350_pars) ;

  c_results->cd (i) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * func_350->GetMinimum (), rangeScale * mass, 1.1 * max (func_sig_350->GetMaximum (), func_sAi_350->GetMaximum ())) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  func_sig_350->SetLineWidth (1) ;
  func_sig_350->SetNpx (10000) ;
  func_sig_350->Draw ("same") ;
  func_sAi_350->SetLineWidth (1) ;
  func_sAi_350->SetLineColor (kRed) ;
  func_sAi_350->SetNpx (10000) ;
  func_sAi_350->Draw ("same") ;
  func_350->SetLineWidth (1) ;
  func_350->SetNpx (10000) ;
  func_350->Draw ("same") ;
  c_results->cd (5 + i++) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * f_IOS_350->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_350->GetMaximum ()) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interf / signal") ;
  f_IOS_350->SetLineWidth (1) ;
  f_IOS_350->SetNpx (10000) ;
  f_IOS_350->Draw ("same") ;
 
  // 500 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 500 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_500_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_500_pars[j] = params_500[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_500_pars[j + 4] = params_sig_500[j] ;
  TF1 * f_IOS_500 = new TF1 ("f_IOS_500", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_500->SetParameters (IOS_500_pars) ;

  c_results->cd (i) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * func_500->GetMinimum (), rangeScale * mass, 1.1 * max (func_sig_500->GetMaximum (), func_sAi_500->GetMaximum ())) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  func_sig_500->SetLineWidth (1) ;
  func_sig_500->SetNpx (10000) ;
  func_sig_500->Draw ("same") ;
  func_sAi_500->SetLineWidth (1) ;
  func_sAi_500->SetLineColor (kRed) ;
  func_sAi_500->SetNpx (10000) ;
  func_sAi_500->Draw ("same") ;
  func_500->SetLineWidth (1) ;
  func_500->SetNpx (10000) ;
  func_500->Draw ("same") ;
  c_results->cd (5 + i++) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * f_IOS_500->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_500->GetMaximum ()) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interf / signal") ;
  f_IOS_500->SetLineWidth (1) ;
  f_IOS_500->SetNpx (10000) ;
  f_IOS_500->Draw ("same") ;

  // 650 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 650 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_650_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_650_pars[j] = params_650[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_650_pars[j + 4] = params_sig_650[j] ;
  TF1 * f_IOS_650 = new TF1 ("f_IOS_650", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_650->SetParameters (IOS_650_pars) ;

  c_results->cd (i) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * func_650->GetMinimum (), rangeScale * mass, 1.1 * max (func_sig_650->GetMaximum (), func_sAi_650->GetMaximum ())) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  func_sig_650->SetLineWidth (1) ;
  func_sig_650->SetNpx (10000) ;
  func_sig_650->Draw ("same") ;
  func_sAi_650->SetLineWidth (1) ;
  func_sAi_650->SetLineColor (kRed) ;
  func_sAi_650->SetNpx (10000) ;
  func_sAi_650->Draw ("same") ;
  func_650->SetLineWidth (1) ;
  func_650->SetNpx (10000) ;
  func_650->Draw ("same") ;
  c_results->cd (5 + i++) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * f_IOS_650->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_650->GetMaximum ()) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interf / signal") ;
  f_IOS_650->SetLineWidth (1) ;
  f_IOS_650->SetNpx (10000) ;
  f_IOS_650->Draw ("same") ;
 
  // 800 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 800 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_800_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_800_pars[j] = params_800[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_800_pars[j + 4] = params_sig_800[j] ;
  TF1 * f_IOS_800 = new TF1 ("f_IOS_800", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_800->SetParameters (IOS_800_pars) ;

  c_results->cd (i) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * func_800->GetMinimum (), rangeScale * mass, 1.1 * max (func_sig_800->GetMaximum (), func_sAi_800->GetMaximum ())) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  func_sig_800->SetLineWidth (1) ;
  func_sig_800->SetNpx (10000) ;
  func_sig_800->Draw ("same") ;
  func_sAi_800->SetLineWidth (1) ;
  func_sAi_800->SetLineColor (kRed) ;
  func_sAi_800->SetNpx (10000) ;
  func_sAi_800->Draw ("same") ;
  func_800->SetLineWidth (1) ;
  func_800->SetNpx (10000) ;
  func_800->Draw ("same") ;
  c_results->cd (5 + i++) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * f_IOS_800->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_800->GetMaximum ()) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interf / signal") ;
  f_IOS_800->SetLineWidth (1) ;
  f_IOS_800->SetNpx (10000) ;
  f_IOS_800->Draw ("same") ;

  // 1000 gev ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

  mass = 1000 ;
  rangeScale = 1.5 ;
//  if (mass[i] > 480) rangeScale = 2 ;

  Double_t IOS_1000_pars[11] ;
  for (int j = 0 ; j < 4 ; ++j) IOS_1000_pars[j] = params_1000[j] ;
  for (int j = 0 ; j < 7 ; ++j) IOS_1000_pars[j + 4] = params_sig_1000[j] ;
  TF1 * f_IOS_1000 = new TF1 ("f_IOS_1000", intereferenceOverSignal, 200, 2000, 11) ;
  f_IOS_1000->SetParameters (IOS_1000_pars) ;

  c_results->cd (i) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * func_1000->GetMinimum (), rangeScale * mass, 1.1 * max (func_sig_1000->GetMaximum (), func_sAi_1000->GetMaximum ())) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  func_sig_1000->SetLineWidth (1) ;
  func_sig_1000->SetNpx (10000) ;
  func_sig_1000->Draw ("same") ;
  func_sAi_1000->SetLineWidth (1) ;
  func_sAi_1000->SetLineColor (kRed) ;
  func_sAi_1000->SetNpx (10000) ;
  func_sAi_1000->Draw ("same") ;
  func_1000->SetLineWidth (1) ;
  func_1000->SetNpx (10000) ;
  func_1000->Draw ("same") ;
  c_results->cd (5 + i++) ;
  bkg = (TH1F *) gPad->DrawFrame (200, 1.1 * f_IOS_1000->GetMinimum (), rangeScale * mass, 1.1 * f_IOS_1000->GetMaximum ()) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interf / signal") ;
  f_IOS_1000->SetLineWidth (1) ;
  f_IOS_1000->SetNpx (10000) ;
  f_IOS_1000->Draw ("same") ;
 
  c_results->Print ("params_masses.pdf", "pdf") ;
  
  //PG put all the interefence shapes together
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  TCanvas * c_merge = new TCanvas () ;
  bkg = (TH1F *) c_merge->DrawFrame (200, 1.1 * func_350->GetMinimum (), 1500, 1.1 * func_350->GetMaximum ()) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interference") ;
  func_350->Draw ("same") ;
  func_500->Draw ("same") ;
  func_650->Draw ("same") ;
  func_800->Draw ("same") ;
  func_1000->Draw ("same") ;
  c_merge->Print ("overlapped_interf.pdf", "pdf") ;


  //PG write the functions for the interpolation
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  int k = 0 ;
//  cout << "TF1 * fit_param0 = new TF1 (\"fit_param0\", \"[0] + [1] * x + [2] * x * x + [3] * x * x * x + [4] * x * x * x * x\") ;\n" ;
//  cout << "fit_param0->FixParameter (" << k++ << ", " << fit_param0->GetParameter (0) << ") ;\n" ;
//  cout << "fit_param0->FixParameter (" << k++ << ", " << fit_param0->GetParameter (1) << ") ;\n" ;
//  cout << "fit_param0->FixParameter (" << k++ << ", " << fit_param0->GetParameter (2) << ") ;\n" ;
//  cout << "fit_param0->FixParameter (" << k++ << ", " << fit_param0->GetParameter (3) << ") ;\n" ;
//  cout << "fit_param0->FixParameter (" << k++ << ", " << fit_param0->GetParameter (4) << ") ;\n" ;
  
  k = 0 ;
  cout << "TF1 * fit_param1 = new TF1 (\"fit_param1\", \"[0] + [1] * x + [2] * x * x\") ;                                        \n" ;
  cout << "fit_param1->FixParameter (" << k++ << ", " << fit_param1->GetParameter (0) << ") ;\n" ;
  cout << "fit_param1->FixParameter (" << k++ << ", " << fit_param1->GetParameter (1) << ") ;\n" ;
  cout << "fit_param1->FixParameter (" << k++ << ", " << fit_param1->GetParameter (2) << ") ;\n" ;
  
  k = 0 ;
  cout << "TF1 * fit_param3 = new TF1 (\"fit_param3\", \"[0] + [1] * x + [2] * x * x\") ;                                        \n" ;
  cout << "fit_param3->FixParameter (" << k++ << ", " << fit_param3->GetParameter (0) << ") ;\n" ;
  cout << "fit_param3->FixParameter (" << k++ << ", " << fit_param3->GetParameter (1) << ") ;\n" ;
  cout << "fit_param3->FixParameter (" << k++ << ", " << fit_param3->GetParameter (2) << ") ;\n" ;

  cout << "TF1 * expo_3 = new TF1 (\"expo_3\",\"TMath::Exp ([0] + [2] * (x - [1]) + [3] * (x - [1]) * (x - [1]) + [4] * (x - [1]) * (x - [1]) * (x - [1]))\", 200, 2000) ;\n" ;
  cout << "expo_3->SetLineColor (kBlue + 1) ;\n" ;
  cout << "expo_3->SetLineWidth (1) ;\n" ;
  cout << "expo_3->SetParameter (0, " << expo_3->GetParameter (0) << ") ;\n" ;
  cout << "expo_3->SetParameter (1, " << expo_3->GetParameter (1) << ") ;\n" ;
  cout << "expo_3->SetParameter (2, " << expo_3->GetParameter (2) << ") ;\n" ;
  cout << "expo_3->SetParameter (3, " << expo_3->GetParameter (3) << ") ;\n" ;
  cout << "expo_3->SetParameter (4, " << expo_3->GetParameter (4) << ") ;\n" ;
  

}  
  
