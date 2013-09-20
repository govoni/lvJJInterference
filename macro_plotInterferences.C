
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


int macro_plotInterferences ()
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
  tg_par0->SetPoint (i, 350, -3.75686e-05) ;
  tg_par1->SetPoint (i, 350, 350.889) ;
  tg_par2->SetPoint (i, 350, 0.0008) ;
  tg_par3->SetPoint (i, 350, 1057.25) ;
  TF1 * func_350 = new TF1 ("func_350",doublePeakModel, 200, 2000, 4) ;
  double params_350[4] = {-3.75686e-05, 350.889, 0.0008, 1057.25 } ;
  func_350->SetParameters (params_350) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 350, 0.000990502) ;
  tg_sig_par1->SetPoint (i, 350, 350.346) ;
  tg_sig_par2->SetPoint (i, 350, 7.64082) ;
  tg_sig_par3->SetPoint (i, 350, 1.15547) ;
  tg_sig_par4->SetPoint (i, 350, 1.55098) ;
  tg_sig_par5->SetPoint (i, 350, 1.26826) ;
  tg_sig_par6->SetPoint (i, 350, 2.15497) ;
  TF1 * func_sig_350 = new TF1 ("func_sig_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_350[7] = {0.000990502, 350.346, 7.64082, 1.15547, 1.55098, 1.26826, 2.15497 } ;
  func_sig_350->SetParameters (params_sig_350) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 350, 0.000984922) ;
  tg_sAi_par1->SetPoint (i, 350, 349.732) ;
  tg_sAi_par2->SetPoint (i, 350, 7.29879) ;
  tg_sAi_par3->SetPoint (i, 350, 1.0354) ;
  tg_sAi_par4->SetPoint (i, 350, 3) ;
  tg_sAi_par5->SetPoint (i, 350, 1) ;
  tg_sAi_par6->SetPoint (i, 350, 3) ;
  TF1 * func_sAi_350 = new TF1 ("func_sAi_350",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_350[7] = {0.000984922, 349.732, 7.29879, 1.0354, 3, 1, 3 } ;
  func_sAi_350->SetParameters (params_sAi_350) ;
  i++ ;
  
  // ----> MASS 500 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 500, -1.01505e-05) ;
  tg_par1->SetPoint (i, 500, 494.806) ;
  tg_par2->SetPoint (i, 500, 0.0008) ;
  tg_par3->SetPoint (i, 500, 11820.7) ;
  TF1 * func_500 = new TF1 ("func_500",doublePeakModel, 200, 2000, 4) ;
  double params_500[4] = {-1.01505e-05, 494.806, 0.0008, 11820.7 } ;
  func_500->SetParameters (params_500) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 500, 8.89251e-05) ;
  tg_sig_par1->SetPoint (i, 500, 504.033) ;
  tg_sig_par2->SetPoint (i, 500, 28.0346) ;
  tg_sig_par3->SetPoint (i, 500, 0.814726) ;
  tg_sig_par4->SetPoint (i, 500, 2.42493) ;
  tg_sig_par5->SetPoint (i, 500, 1.01087) ;
  tg_sig_par6->SetPoint (i, 500, 4.83435) ;
  TF1 * func_sig_500 = new TF1 ("func_sig_500",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_500[7] = {8.89251e-05, 504.033, 28.0346, 0.814726, 2.42493, 1.01087, 4.83435 } ;
  func_sig_500->SetParameters (params_sig_500) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 500, 8.25563e-05) ;
  tg_sAi_par1->SetPoint (i, 500, 498.203) ;
  tg_sAi_par2->SetPoint (i, 500, -29.9421) ;
  tg_sAi_par3->SetPoint (i, 500, 1.22454) ;
  tg_sAi_par4->SetPoint (i, 500, 3) ;
  tg_sAi_par5->SetPoint (i, 500, 1) ;
  tg_sAi_par6->SetPoint (i, 500, 3) ;
  TF1 * func_sAi_500 = new TF1 ("func_sAi_500",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_500[7] = {8.25563e-05, 498.203, -29.9421, 1.22454, 3, 1, 3 } ;
  func_sAi_500->SetParameters (params_sAi_500) ;
  i++ ;
  
  // ----> MASS 650 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 650, -3.65649e-06) ;
  tg_par1->SetPoint (i, 650, 649.011) ;
  tg_par2->SetPoint (i, 650, 0.0008) ;
  tg_par3->SetPoint (i, 650, 23478.6) ;
  TF1 * func_650 = new TF1 ("func_650",doublePeakModel, 200, 2000, 4) ;
  double params_650[4] = {-3.65649e-06, 649.011, 0.0008, 23478.6 } ;
  func_650->SetParameters (params_650) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 650, 1.87391e-05) ;
  tg_sig_par1->SetPoint (i, 650, 661.553) ;
  tg_sig_par2->SetPoint (i, 650, 62.2655) ;
  tg_sig_par3->SetPoint (i, 650, 0.797603) ;
  tg_sig_par4->SetPoint (i, 650, 6.13542) ;
  tg_sig_par5->SetPoint (i, 650, 0.955861) ;
  tg_sig_par6->SetPoint (i, 650, 15.3988) ;
  TF1 * func_sig_650 = new TF1 ("func_sig_650",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_650[7] = {1.87391e-05, 661.553, 62.2655, 0.797603, 6.13542, 0.955861, 15.3988 } ;
  func_sig_650->SetParameters (params_sig_650) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 650, 1.76301e-05) ;
  tg_sAi_par1->SetPoint (i, 650, 642.004) ;
  tg_sAi_par2->SetPoint (i, 650, -65.4949) ;
  tg_sAi_par3->SetPoint (i, 650, 1.67731) ;
  tg_sAi_par4->SetPoint (i, 650, 3) ;
  tg_sAi_par5->SetPoint (i, 650, 1) ;
  tg_sAi_par6->SetPoint (i, 650, 3) ;
  TF1 * func_sAi_650 = new TF1 ("func_sAi_650",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_650[7] = {1.76301e-05, 642.004, -65.4949, 1.67731, 3, 1, 3 } ;
  func_sAi_650->SetParameters (params_sAi_650) ;
  i++ ;
  
  // ----> MASS 800 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 800, -1.71241e-06) ;
  tg_par1->SetPoint (i, 800, 763.576) ;
  tg_par2->SetPoint (i, 800, 0.0008) ;
  tg_par3->SetPoint (i, 800, 36336.1) ;
  TF1 * func_800 = new TF1 ("func_800",doublePeakModel, 200, 2000, 4) ;
  double params_800[4] = {-1.71241e-06, 763.576, 0.0008, 36336.1 } ;
  func_800->SetParameters (params_800) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 800, 5.43552e-06) ;
  tg_sig_par1->SetPoint (i, 800, 806.479) ;
  tg_sig_par2->SetPoint (i, 800, 108.152) ;
  tg_sig_par3->SetPoint (i, 800, 1.0251) ;
  tg_sig_par4->SetPoint (i, 800, 8.15215) ;
  tg_sig_par5->SetPoint (i, 800, 0.951711) ;
  tg_sig_par6->SetPoint (i, 800, 141.733) ;
  TF1 * func_sig_800 = new TF1 ("func_sig_800",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_800[7] = {5.43552e-06, 806.479, 108.152, 1.0251, 8.15215, 0.951711, 141.733 } ;
  func_sig_800->SetParameters (params_sig_800) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 800, 4.91246e-06) ;
  tg_sAi_par1->SetPoint (i, 800, 744.754) ;
  tg_sAi_par2->SetPoint (i, 800, -115.026) ;
  tg_sAi_par3->SetPoint (i, 800, 4.92705) ;
  tg_sAi_par4->SetPoint (i, 800, 3) ;
  tg_sAi_par5->SetPoint (i, 800, 1) ;
  tg_sAi_par6->SetPoint (i, 800, 3) ;
  TF1 * func_sAi_800 = new TF1 ("func_sAi_800",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_800[7] = {4.91246e-06, 744.754, -115.026, 4.92705, 3, 1, 3 } ;
  func_sAi_800->SetParameters (params_sAi_800) ;
  i++ ;
  
  // ----> MASS 1000 ---- ---- ---- 

  // interference parametrisation:
  tg_par0->SetPoint (i, 1000, -7.22506e-07) ;
  tg_par1->SetPoint (i, 1000, 817.131) ;
  tg_par2->SetPoint (i, 1000, 0.0008) ;
  tg_par3->SetPoint (i, 1000, 85380.2) ;
  TF1 * func_1000 = new TF1 ("func_1000",doublePeakModel, 200, 2000, 4) ;
  double params_1000[4] = {-7.22506e-07, 817.131, 0.0008, 85380.2 } ;
  func_1000->SetParameters (params_1000) ;
  // MG signal only parametrisation:
  tg_sig_par0->SetPoint (i, 1000, 1.26711e-06) ;
  tg_sig_par1->SetPoint (i, 1000, 933.347) ;
  tg_sig_par2->SetPoint (i, 1000, 224.458) ;
  tg_sig_par3->SetPoint (i, 1000, 1.30561) ;
  tg_sig_par4->SetPoint (i, 1000, 92.1374) ;
  tg_sig_par5->SetPoint (i, 1000, -252354) ;
  tg_sig_par6->SetPoint (i, 1000, 77.2872) ;
  TF1 * func_sig_1000 = new TF1 ("func_sig_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sig_1000[7] = {1.26711e-06, 933.347, 224.458, 1.30561, 92.1374, -252354, 77.2872 } ;
  func_sig_1000->SetParameters (params_sig_1000) ;
  // PG SBI - B  parametrisation:
  tg_sAi_par0->SetPoint (i, 1000, 9.78941e-07) ;
  tg_sAi_par1->SetPoint (i, 1000, 753.335) ;
  tg_sAi_par2->SetPoint (i, 1000, -222.197) ;
  tg_sAi_par3->SetPoint (i, 1000, 4.01578) ;
  tg_sAi_par4->SetPoint (i, 1000, 3) ;
  tg_sAi_par5->SetPoint (i, 1000, 1) ;
  tg_sAi_par6->SetPoint (i, 1000, 3) ;
  TF1 * func_sAi_1000 = new TF1 ("func_sAi_1000",crystalBallLowHigh, 200, 2000, 7) ;
  double params_sAi_1000[7] = {9.78941e-07, 753.335, -222.197, 4.01578, 3, 1, 3 } ;
  func_sAi_1000->SetParameters (params_sAi_1000) ;
  i++ ;


  //PG with the fitted function
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

  TF1 * fit_param1 = new TF1 ("fit_param1", "[0] + [1] * x + [2] * x * x") ;                                        
  fit_param1->FixParameter (0, -201.868) ;
  fit_param1->FixParameter (1, 1.83565) ;
  fit_param1->FixParameter (2, -0.000810224) ;
  TF1 * fit_param3 = new TF1 ("fit_param3", "[0] + [1] * x + [2] * x * x") ;                                        
  fit_param3->FixParameter (0, 24345.2) ;
  fit_param3->FixParameter (1, -123.834) ;
  fit_param3->FixParameter (2, 0.182752) ;

  TCanvas * c_merge = new TCanvas () ;
  bkg = (TH1F *) c_merge->DrawFrame (200, -0.00005, 1500, 0.00005) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interference") ;

  int legred = 0 ; int legblue = 0 ; int legblack = 0 ; int leggray = 0 ;
  leg = new TLegend (0.6,0.6,0.85,0.8) ;
  leg->SetTextFont (42) ;
  leg->SetBorderSize (0) ;
  leg->SetFillStyle (0) ;

  for (double mass = 300 ; mass < 1100 ; mass += 25)
    {
      TF1 * func = new TF1 ("func", doublePeakModel, 200, 2000, 4) ;
      double params[4] = {tg_par0->Eval (mass), fit_param1->Eval (mass), 0.0008, fit_param3->Eval (mass) } ;
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
          func->SetLineColor (kGray) ;
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

  double params_350[4] = {-3.75686e-05, 350.889, 0.0008, 1057.25 } ;
  double params_500[4] = {-1.01505e-05, 494.806, 0.0008, 11820.7 } ;
  double params_650[4] = {-3.65649e-06, 649.011, 0.0008, 23478.6 } ;
  double params_800[4] = {-1.71241e-06, 763.576, 0.0008, 36336.1 } ;
  double params_1000[4] = {-7.22506e-07, 817.131, 0.0008, 85380.2 } ;

  TF1 * func_350 = new TF1 ("func_350",doublePeakModel, 200, 2000, 4) ;
  func_350->SetParameters (params_350) ;
  TF1 * func_500 = new TF1 ("func_500",doublePeakModel, 200, 2000, 4) ;
  func_500->SetParameters (params_500) ;
  TF1 * func_650 = new TF1 ("func_650",doublePeakModel, 200, 2000, 4) ;
  func_650->SetParameters (params_650) ;
  TF1 * func_800 = new TF1 ("func_800",doublePeakModel, 200, 2000, 4) ;
  func_800->SetParameters (params_800) ;
  TF1 * func_1000 = new TF1 ("func_1000",doublePeakModel, 200, 2000, 4) ;
  func_1000->SetParameters (params_1000) ;

  leg->AddEntry (func_1000, "origi. samples", "l") ;


  func_350->SetLineWidth (1) ;
  func_350->SetNpx (10000) ;
  func_500->SetLineWidth (1) ;
  func_500->SetNpx (10000) ;
  func_650->SetLineWidth (1) ;
  func_650->SetNpx (10000) ;
  func_800->SetLineWidth (1) ;
  func_800->SetNpx (10000) ;
  func_1000->SetLineWidth (1) ;
  func_1000->SetNpx (10000) ;

  func_350-> SetLineStyle (2) ;
  func_500-> SetLineStyle (2) ;
  func_650-> SetLineStyle (2) ;
  func_800-> SetLineStyle (2) ;
  func_1000->SetLineStyle (2) ;


  func_350->Draw ("same") ;
  func_500->Draw ("same") ;
  func_650->Draw ("same") ;
  func_800->Draw ("same") ;
  func_1000->Draw ("same") ;

  leg->Draw () ;
  c_merge->Print ("interf_fit.pdf", "pdf") ;

  //PG with linear interpolation
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  

  bkg->Draw () ;
  for (double mass = 300 ; mass < 1100 ; mass += 25)
    {
      TF1 * func = new TF1 ("func", doublePeakModel, 200, 2000, 4) ;
      double params[4] = {tg_par0->Eval (mass), tg_par1->Eval (mass), tg_par2->Eval (mass), tg_par3->Eval (mass) } ;
      func->SetParameters (params) ;
      func->SetLineWidth (1) ;
      func->SetNpx (10000) ;
      if (int (mass) % 100 == 0) func->SetLineColor (kRed) ;
      else                       func->SetLineColor (kGray) ;
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
  
  func_350->Draw ("same") ;
  func_500->Draw ("same") ;
  func_650->Draw ("same") ;
  func_800->Draw ("same") ;
  func_1000->Draw ("same") ;

  leg->Draw () ;
  c_merge->Print ("interf_interpol.pdf", "pdf") ;
  
  //PG try a hybrid combination done by heart
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  
  bkg->Draw () ;
  for (double mass = 300 ; mass < 1100 ; mass += 25)
    {
      TF1 * func = new TF1 ("func", doublePeakModel, 200, 2000, 4) ;
//      double params[4] = {fit_param0->Eval (mass), fit_param1->Eval (mass), 0.0008, fit_param3->Eval (mass) } ;
//      double params[4] = {tg_par0->Eval (mass), tg_par1->Eval (mass), tg_par2->Eval (mass), tg_par3->Eval (mass) } ;
      double params[4] = {tg_par0->Eval (mass), fit_param1->Eval (mass), 0.0008, tg_par3->Eval (mass) } ;
      func->SetParameters (params) ;
      func->SetLineWidth (1) ;
      func->SetNpx (10000) ;
      if (int (mass) % 100 == 0) func->SetLineColor (kRed) ;
      else                       func->SetLineColor (kGray) ;
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
  
  func_350->Draw ("same") ;
  func_500->Draw ("same") ;
  func_650->Draw ("same") ;
  func_800->Draw ("same") ;
  func_1000->Draw ("same") ;

  leg->Draw () ;
  c_merge->Print ("interf_hybrid.pdf", "pdf") ;
  

}
