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
      f_sAi = (TF1 *) ((TF1 *) input.Get ("func_ph_2"))->Clone (name) ;
      input.Close () ;
    } 
  void setsigGraphs (TGraph ** local_tg_sig_par, int index, int Npar = 7)
    {
      for (int k = 0 ; k < Npar ; ++k)
        local_tg_sig_par[k]->SetPoint (index, m_mass, f_sig->GetParameter (k)) ; 
      return ;
    }
  void setsAiGraphs (TGraph ** local_tg_sAi_par, int index, int Npar = 9)
    {
      for (int k = 0 ; k < Npar ; ++k)
        local_tg_sAi_par[k]->SetPoint (index, m_mass, f_sAi->GetParameter (k)) ; 
      return ;
    }

} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


struct corrections
{
  TGraph ** tg_sig_par ;
  TGraph ** tg_sAi_par ;
  string m_path ;
  string m_scale ;
  TF1 * f_sig ;
  TF1 * f_sAi ;
  
  corrections (string path, string scale = "1") : m_path (path), m_scale (scale)
    {
      tg_sig_par = new TGraph * [7] ; // [parameter][mass]
      for (int k = 0 ; k < 7 ; ++k) tg_sig_par[k] = new TGraph (5) ;
      tg_sAi_par = new TGraph * [9] ;
      for (int k = 0 ; k < 9 ; ++k) tg_sAi_par[k] = new TGraph (5) ;

      //PG fill the graphs from the input files
      //PG ---- ---- ---- ---- ---- ---- ---- ----
    
      int i = 0 ; 
    
      intCont ic_350  (m_path + string ("/results_interference.350.") + m_scale + string (".root"), 350) ;      
      ic_350.setsigGraphs (tg_sig_par, i) ;
      ic_350.setsAiGraphs (tg_sAi_par, i++) ;
      intCont ic_500  (m_path + string ("/results_interference.500.") + m_scale + string (".root"), 500) ;      
      ic_500.setsigGraphs (tg_sig_par, i) ;  
      ic_500.setsAiGraphs (tg_sAi_par, i++) ;  
      intCont ic_650  (m_path + string ("/results_interference.650.") + m_scale + string (".root"), 650) ;      
      ic_650.setsigGraphs (tg_sig_par, i) ;  
      ic_650.setsAiGraphs (tg_sAi_par, i++) ;  
      intCont ic_800  (m_path + string ("/results_interference.800.") + m_scale + string (".root"), 800) ;      
      ic_800.setsigGraphs (tg_sig_par, i) ;  
      ic_800.setsAiGraphs (tg_sAi_par, i++) ;  
      intCont ic_1000 (m_path + string ("/results_interference.1000.") + m_scale + string (".root"), 1000) ;   
      ic_1000.setsigGraphs (tg_sig_par, i) ;   
      ic_1000.setsAiGraphs (tg_sAi_par, i++) ;   
    }


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


  TF1 * initiateFunction_sig (TGraph ** tg_sample_par, TString fname, double mass)
    {
      TGraph * log_tg_sample_par0 = makeLog (tg_sample_par[0]) ;
      TF1 * f_sample = new TF1 (fname, crystalBallLowHigh, 200, 2000, 7) ;
      f_sample->SetParameter (0, TMath::Exp (log_tg_sample_par0->Eval (mass))) ;
      for (int iPar = 1 ; iPar < 7 ; ++iPar)
        {
          f_sample->SetParameter (iPar, tg_sample_par[iPar]->Eval (mass)) ;
        }
      return f_sample ;
    }


  TF1 * initiateFunction_sAi (TGraph ** tg_sample_par, TString fname, double mass)
    {
      TGraph * log_tg_sample_par0 = makeLog (tg_sample_par[0]) ;
      TF1 * f_sample = new TF1 (fname, crystalBallLowHighWithRise, 200, 2000, 9) ;
      f_sample->SetParameter (0, TMath::Exp (log_tg_sample_par0->Eval (mass))) ;
      for (int iPar = 1 ; iPar < 9 ; ++iPar)
        {
          f_sample->SetParameter (iPar, tg_sample_par[iPar]->Eval (mass)) ;
        }
      return f_sample ;
    }


  void prepareFunctions (float mass)
    {
      f_sig = initiateFunction_sig (tg_sig_par, "f_sig", mass) ;
      f_sAi = initiateFunction_sAi (tg_sAi_par, "f_sAi", mass) ;
      return ;
    }


  float getCorrFactor (float mWW)
    {
      float val = f_sig->Eval (mWW) ;
      if (val == 0) return 0. ;
      return f_sAi->Eval (mWW) / val ;
    }

} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


TH1F * correctHisto (TString outputName, TH1F * input, corrections & corrector)
{
  TH1F * output = (TH1F *) input->Clone (outputName) ;
  for (int i = 1 ; i <= input->GetNbinsX () ; ++i)
    {
      float mWW = input->GetBinCenter (i) ;
      output->SetBinContent (i, input->GetBinContent (i) * corrector.getCorrFactor (mWW)) ;
    }
  return output ;  
}


// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
// ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====


int macro_testCorrection () 
{

  string correctionsFolder = "./" ;
//  float mass = 500. ;
  float mass = 800. ;

  //PG read the correction factors
  corrections corrTool_1 (correctionsFolder, "1") ;    corrTool_1.prepareFunctions (mass) ; // central
  corrections corrTool_H (correctionsFolder, "0.5") ;  corrTool_H.prepareFunctions (mass) ; // half
  corrections corrTool_D (correctionsFolder, "2") ;    corrTool_D.prepareFunctions (mass) ; // double
  
//  TFile f_signal ("testCorrection.500.root") ;
  TFile f_signal ("testCorrection.800.root") ;
  TH1F * h_MWW_mg = (TH1F *) f_signal.Get ("h_MWW_mg") ;
  h_MWW_mg->Rebin (2) ;

  TH1F * correct_1 = correctHisto ("corr_1", h_MWW_mg, corrTool_1) ;
  TH1F * correct_H = correctHisto ("corr_H", h_MWW_mg, corrTool_H) ;
  TH1F * correct_D = correctHisto ("corr_D", h_MWW_mg, corrTool_D) ;

  correct_H->SetLineColor (kGray + 1) ;
  correct_D->SetLineColor (kGray + 1) ;

  correct_H->SetFillColor (kGray) ;
  correct_1->SetLineColor (kGreen + 2) ;
  correct_D->SetFillColor (10) ;

  
  TCanvas * c1 = new TCanvas ("c1") ;
  TH1F * h_c1 = (TH1F *) c1->DrawFrame (200, 0., 1200, 0.01) ;
  h_c1->SetTitle (0) ;
  h_c1->SetStats (0) ;
  h_c1->GetXaxis ()->SetTitle ("m_{WW} (GeV)") ;
  correct_H->Draw ("histo same") ;
  correct_D->Draw ("histo same") ;
  h_MWW_mg->Draw ("same") ;
  correct_1->Draw ("same") ;


//  c1->Print ("reweighedSig.500.pdf", "pdf") ;
  c1->Print ("reweighedSig.800.pdf", "pdf") ;
  

  return 0 ;  

}  
  
