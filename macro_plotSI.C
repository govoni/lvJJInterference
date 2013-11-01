
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


int macro_plotSI ()
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
  tg_sig_par2->SetPoint (i, 350, 7.70265) ;
  tg_sig_par3->SetPoint (i, 350, 1.16314) ;
  tg_sig_par4->SetPoint (i, 350, 1.53064) ;
  tg_sig_par5->SetPoint (i, 350, 1.27642) ;
  tg_sig_par6->SetPoint (i, 350, 2.17316) ;
  TF1 * func_sig_350 = new TF1 ("func_sig_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_350[7] = {0.000952886, 350.384, 7.70265, 1.16314, 1.53064, 1.27642, 2.17316 } ;
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
  tg_sAi_par0->SetPoint (i, 1000, 9.86367e-07) ;
  tg_sAi_par1->SetPoint (i, 1000, 745.645) ;
  tg_sAi_par2->SetPoint (i, 1000, -226.281) ;
  tg_sAi_par3->SetPoint (i, 1000, 4.04549) ;
  tg_sAi_par4->SetPoint (i, 1000, 3) ;
  tg_sAi_par5->SetPoint (i, 1000, 1) ;
  tg_sAi_par6->SetPoint (i, 1000, 3) ;
  TF1 * func_sAi_1000 = new TF1 ("func_sAi_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_1000[7] = {9.86367e-07, 745.645, -226.281, 4.04549, 3, 1, 3 } ;
  func_sAi_1000->SetParameters (params_sAi_1000) ;
  i++ ;

  //PG with the interpolation of the parameters
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


  TCanvas * c_merge = new TCanvas () ;
  bkg = (TH1F *) c_merge->DrawFrame (200, 0.0000001, 1500, 0.002) ;
  c_merge->SetLogy () ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("signal and interference") ;

  TLine * l_leftLim = new TLine (350, 0.0000001, 350, 0.002) ;
  l_leftLim->SetLineColor (kGray + 1) ;
  l_leftLim->SetLineStyle (2) ;
  l_leftLim->SetLineWidth (1) ;

  int legred = 0 ; int legblue = 0 ; int legblack = 0 ; int leggray = 0 ;
  leg = new TLegend (0.6,0.6,0.85,0.8) ;
  leg->SetTextFont (42) ;
  leg->SetBorderSize (0) ;
  leg->SetFillStyle (0) ;

  for (double mass = 300 ; mass < 1000 ; mass += 25)
    {
      TF1 * func = new TF1 ("func",crystalBallLowHigh, 200, 2000, 7) ;
      double params[7] = {tg_sAi_par0->Eval (mass), 
                          tg_sAi_par1->Eval (mass), 
                          tg_sAi_par2->Eval (mass), 
                          tg_sAi_par3->Eval (mass), 
                          tg_sAi_par4->Eval (mass), 
                          tg_sAi_par5->Eval (mass), 
                          tg_sAi_par6->Eval (mass)} ;
      func->SetParameters (params) ;
      func->SetLineWidth (1) ;
      func->SetNpx (10000) ;
      if (int (mass) % 100 == 0) 
        {
          func->SetLineColor (kRed) ;
          if (legred == 0)
            {
              leg->AddEntry (func, "N x 100 GeV mass", "l") ;
              legred = 1 ;
            }
        }  
      else
        {
          func->SetLineColor (kGray + 2) ;
          if (leggray == 0)
            {
              leg->AddEntry (func, "N x 25 GeV mass", "l") ;
              leggray = 1 ;
            }           
        }  
      if (mass == 350 ||
          mass == 500 ||
          mass == 650 ||
          mass == 800 ||
          mass == 1000)          
        {  
          if (legblue == 0)
            {
              leg->AddEntry (func, "extrap. samples", "l") ;
              legblue = 1 ;
            }           
          func->SetLineColor (kBlue) ;
        }
      func->Draw ("same") ;
    }

  leg->AddEntry (func_sig_1000, "origi. samples", "l") ;

  func_sAi_350->SetLineWidth (1) ;
  func_sAi_350->SetNpx (10000) ;
  func_sAi_500->SetLineWidth (1) ;
  func_sAi_500->SetNpx (10000) ;
  func_sAi_650->SetLineWidth (1) ;
  func_sAi_650->SetNpx (10000) ;
  func_sAi_800->SetLineWidth (1) ;
  func_sAi_800->SetNpx (10000) ;
  func_sAi_1000->SetLineWidth (1) ;
  func_sAi_1000->SetNpx (10000) ;

  func_350-> SetLineStyle (2) ;
  func_500-> SetLineStyle (2) ;
  func_650-> SetLineStyle (2) ;
  func_800-> SetLineStyle (2) ;
  func_1000->SetLineStyle (2) ;

  func_sig_350->Draw ("same") ;
  func_sig_500->Draw ("same") ;
  func_sig_650->Draw ("same") ;
  func_sig_800->Draw ("same") ;
  func_sig_1000->Draw ("same") ;

  leg->Draw () ;
  l_leftLim->Draw ("same") ;
  c_merge->Print ("sAi_interpol.pdf", "pdf") ;
  c_merge->SetLogy (0) ;
  c_merge->Print ("sAi_interpol_lin.pdf", "pdf") ;

  //PG with the fit on the first parameter
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TF1 * expo_2 = new TF1 ("expo_2","TMath::Exp ([0] + [2] * (x - [1]) + [3] * (x - [1]) * (x - [1]))", 200, 2000) ;
  expo_2->SetLineColor (kRed + 1) ;
  expo_2->SetLineWidth (1) ;
  expo_2->SetParameter (0,    -0.680127) ;
  expo_2->SetParameter (1,    -0.336481) ;
  expo_2->SetParameter (2,    -0.020906) ;
  expo_2->SetParameter (3,  7.83171e-06) ;

  TF1 * expo_3 = new TF1 ("expo_3","TMath::Exp ([0] + [2] * (x - [1]) + [3] * (x - [1]) * (x - [1]) + [4] * (x - [1]) * (x - [1]) * (x - [1]))", 200, 2000) ;
  expo_3->SetLineColor (kBlue + 1) ;
  expo_3->SetLineWidth (1) ;
  expo_3->SetParameter (0,      4.38191) ;
  expo_3->SetParameter (1,     -5.32204) ;
  expo_3->SetParameter (2,   -0.0463521) ;
  expo_3->SetParameter (3,  4.77379e-05) ;
  expo_3->SetParameter (4, -1.95508e-08) ;

  bkg->Draw () ;
  c_merge->SetLogy (1) ;

  for (double mass = 300 ; mass < 1000 ; mass += 25)
    {
      TF1 * func = new TF1 ("func",crystalBallLowHigh, 200, 2000, 7) ;
      double params[7] = {expo_3->Eval (mass), tg_sAi_par1->Eval (mass), tg_sAi_par2->Eval (mass), 
                          tg_sAi_par3->Eval (mass), tg_sAi_par4->Eval (mass), tg_sAi_par5->Eval (mass), 
                          tg_sAi_par6->Eval (mass)} ;
      func->SetParameters (params) ;
      func->SetLineWidth (1) ;
      func->SetNpx (10000) ;
      if (int (mass) % 100 == 0) 
        {
          func->SetLineColor (kRed) ;
        }  
      else
        {
          func->SetLineColor (kGray + 2) ;
        }  
      if (mass == 350 ||
          mass == 500 ||
          mass == 650 ||
          mass == 800 ||
          mass == 1000)          
        {  
          func->SetLineColor (kBlue) ;
        }
      func->Draw ("same") ;
    }
  func_sAi_350->Draw ("same") ;
  func_sAi_500->Draw ("same") ;
  func_sAi_650->Draw ("same") ;
  func_sAi_800->Draw ("same") ;
  func_sAi_1000->Draw ("same") ;
  
  leg->Draw () ;
  l_leftLim->Draw ("same") ;
  c_merge->Print ("sAi_hybrid.pdf", "pdf") ;
  c_merge->SetLogy (0) ;
  c_merge->Print ("sAi_hybrid_lin.pdf", "pdf") ;
  

}
