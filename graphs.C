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


struct intCont
{
  
  TF1 * f_sig ;
  TF1 * f_sAi ;
  float m_mass ;
  
  intCont (TString filename, float mass)
    {
      m_mass = mass ;
      TFile input (filename) ;
      TString name = "f_sig_" ;
      name += m_mass ;
      f_sig = (TF1 *) ((TF1 *) input.Get ("func_mg_1"))->Clone (name) ;
      f_sAi = (TF1 *) ((TF1 *) input.Get ("func_ph_1"))->Clone (name) ;
      input.Close () ;
    } 
  void setsigGraphs (TGraph ** local_tg_sig_par, int index)
    {
      for (int k = 0 ; k < 7 ; ++k)
        local_tg_sig_par[k]->SetPoint (index, m_mass, f_sig->GetParameter (k)) ; 
      return ;
    }
  void setsAiGraphs (TGraph ** local_tg_sAi_par, int index)
    {
      for (int k = 0 ; k < 7 ; ++k)
        local_tg_sAi_par[k]->SetPoint (index, m_mass, f_sAi->GetParameter (k)) ; 
      return ;
    }

} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


TGraph * makeLog (TGraph * orig)
  {
    TGraph * dummy = new TGraph (orig->GetN ()) ;
    for (int k = 0 ; k < orig->GetN () ; ++k) 
      {
        double x, y ;
        orig->GetPoint (k, x, y) ;
        dummy->SetPoint (k, x, TMath::Log (y)) ;
      }
    return dummy ;
  }


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void drawShapes (TGraph ** tg_sample_par)
  {
    TGraph * log_tg_sample_par0 = makeLog (tg_sample_par[0]) ;
    //PG loop on candidate masses
    for (double mass = 300 ; mass < 1100 ; mass += 25)
      {
        TF1 * func = new TF1 ("func",crystalBallLowHigh, 200, 2000, 7) ;
        func->SetParameter (0, TMath::Exp (log_tg_sample_par0->Eval (mass))) ;
        for (int iPar = 1 ; iPar < 7 ; ++iPar)
          func->SetParameter (iPar, tg_sample_par[iPar]->Eval (mass)) ;
        func->SetLineWidth (1) ;
        func->SetLineColor (kBlue + 1) ;
        func->SetNpx (10000) ;
        func->Draw ("same") ;
      } //PG loop on candidate masses
    return ;
  }


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void showInterpolation (TGraph ** tg_sample_par, TString sample)
  {

    TCanvas * c_merge = new TCanvas () ;
    bkg = (TH1F *) c_merge->DrawFrame (200, 0.0000001, 1500, 0.002) ;
    c_merge->SetLogy () ;
    bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
    bkg->GetYaxis ()->SetTitle ("signal") ;

    drawShapes (tg_sample_par) ;
    
    c_merge->Print (TString ("interpol_log_") + sample + TString (".pdf"), "pdf") ;
    c_merge->SetLogy (0) ;
    c_merge->Print (TString ("interpol_lin_") + sample + TString (".pdf"), "pdf") ;
    return ;
  }


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void testInterpolation (TGraph ** tg_sample_par, TString sample)
  {

    //PG testing sets for the test
    TGraph ** tg_test_par = new TGraph * [7] ; // [parameter][mass]
    for (int k = 0 ; k < 7 ; ++k) tg_test_par[k] = new TGraph (4) ;
    
    TCanvas * c_merge = new TCanvas () ;
    bkg = (TH1F *) c_merge->DrawFrame (200, 0.0000001, 1500, 0.002) ;
    c_merge->SetLogy () ;
    bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
    bkg->GetYaxis ()->SetTitle ("signal") ;

    //PG loop on the point to be removed
    for (int iTest = 0 ; iTest < 5 ; ++iTest)
      {
        cout << " ==== iTest " << iTest << endl ;
        int iActualMass = 0 ;
        //PG fill the temporary graphs
        for (int iMass = 0 ; iMass < 5 ; ++iMass)
          {
            if (iMass == iTest) continue ;
            double x = 0 ; 
            double y = 0 ;
            for (int k = 0 ; k < 7 ; ++k) 
              {
                tg_sample_par[k]->GetPoint (iMass, x, y) ;      
                tg_test_par[k]->SetPoint (iActualMass, x, y) ;
              }
            ++iActualMass ;
          } //PG fill the temporary graphs
        
        // plot the interpolation
        bkg->Draw () ;
        drawShapes (tg_test_par) ;

        // overlap the missing point
        double dummy = 0 ;
        double mass = 0 ;
        tg_sample_par[0]->GetPoint (iTest, mass, dummy) ;
        TF1 * func = new TF1 ("func",crystalBallLowHigh, 200, 2000, 7) ;
        TGraph * log_tg_this_sample_par0 = makeLog (tg_sample_par[0]) ;
        func->SetParameter (0, TMath::Exp (log_tg_this_sample_par0->Eval (mass))) ;
        for (int iPar = 1 ; iPar < 7 ; ++iPar)
          func->SetParameter (iPar, tg_sample_par[iPar]->Eval (mass)) ;
        func->SetLineColor (kRed) ;
        func->Draw ("same") ;
        
        // output        
        c_merge->SetLogy (1) ;
        TString name = "test_interpol_log_" ;
        name += sample ; name += "_" ; name += iTest ; name += ".pdf" ;
        c_merge->Print (name, "pdf") ;
        name = "test_interpol_lin_" ;
        name += sample ; name += "_" ; name += iTest ; name += ".pdf" ;
        c_merge->SetLogy (0) ;
        c_merge->Print (name, "pdf") ;
      
      } //PG loop on the point to be removed

    return ;
  }


// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====


int graphs () 
{
  //PG prepare the graphs
  //PG ---- ---- ---- ---- ---- ---- ---- ----

  TGraph ** tg_sig_par = new TGraph * [7] ; // [parameter][mass]
  for (int k = 0 ; k < 7 ; ++k) tg_sig_par[k] = new TGraph (5) ;
  TGraph ** tg_sAi_par = new TGraph * [7] ;
  for (int k = 0 ; k < 7 ; ++k) tg_sAi_par[k] = new TGraph (5) ;

  //PG fill the graphs from the input files
  //PG ---- ---- ---- ---- ---- ---- ---- ----

  int i = 0 ; 

  intCont ic_350  ("results_interference.350.root", 350) ;      
  ic_350.setsigGraphs (tg_sig_par, i) ;
  ic_350.setsAiGraphs (tg_sAi_par, i++) ;
  intCont ic_500  ("results_interference.500.root", 500) ;      
  ic_500.setsigGraphs (tg_sig_par, i) ;  
  ic_500.setsAiGraphs (tg_sAi_par, i++) ;  
  intCont ic_650  ("results_interference.650.root", 650) ;      
  ic_650.setsigGraphs (tg_sig_par, i) ;  
  ic_650.setsAiGraphs (tg_sAi_par, i++) ;  
  intCont ic_800  ("results_interference.800.root", 800) ;      
  ic_800.setsigGraphs (tg_sig_par, i) ;  
  ic_800.setsAiGraphs (tg_sAi_par, i++) ;  
  intCont ic_1000 ("results_interference.1000.root", 1000) ;   
  ic_1000.setsigGraphs (tg_sig_par, i) ;   
  ic_1000.setsAiGraphs (tg_sAi_par, i++) ;   

  //PG plot the interpolation for the signal
  //PG ---- ---- ---- ---- ---- ---- ---- ----
  showInterpolation (tg_sig_par, "sig") ;
  showInterpolation (tg_sAi_par, "sAi") ;

  testInterpolation (tg_sig_par, "sig") ;
  testInterpolation (tg_sAi_par, "sAi") ;


  return 0 ;  

}  
  
