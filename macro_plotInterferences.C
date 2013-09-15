

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


int macro_plotInterferences ()
{
  TF1 * fit_param0 = new TF1 ("fit_param0", "[0] + [1] * x + [2] * x * x + [3] * x * x * x + [4] * x * x * x * x") ;
  fit_param0->FixParameter (0, 0.000134687) ;
  fit_param0->FixParameter (1, -8.60149e-07) ;
  fit_param0->FixParameter (2, 1.89306e-09) ;
  fit_param0->FixParameter (3, -1.75111e-12) ;
  fit_param0->FixParameter (4, 5.79309e-16) ;
  TF1 * fit_param1 = new TF1 ("fit_param1", "[0] + [1] * x + [2] * x * x") ;                                        
  fit_param1->FixParameter (0, -205.885) ;
  fit_param1->FixParameter (1, 1.85257) ;
  fit_param1->FixParameter (2, -0.000827262) ;
  TF1 * fit_param3 = new TF1 ("fit_param3", "[0] + [1] * x + [2] * x * x") ;                                        
  fit_param3->FixParameter (0, 24472.7) ;
  fit_param3->FixParameter (1, -124.226) ;
  fit_param3->FixParameter (2, 0.182908) ;

  TCanvas * c_merge = new TCanvas () ;
  bkg = (TH1F *) c_merge->DrawFrame (200, -0.0004, 1500, 0.0004) ;
  bkg->GetXaxis ()->SetTitle ("m_{WW}") ;
  bkg->GetYaxis ()->SetTitle ("interference") ;

  for (double mass = 300 ; mass < 1100 ; mass += 25)
    {
      TF1 * func = new TF1 ("func", doublePeakModel, 200, 2000, 4) ;
      double params[4] = {fit_param0->Eval (mass), fit_param1->Eval (mass), 0.0008, fit_param3->Eval (mass) } ;
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
          func->Draw ("same") ;
        }
    }

  double params_350[4] = {-7.27716e-07, 350.92, 0.0008, 1057.76 } ;
  double params_500[4] = {-4.85722e-06, 494.671, 0.0008, 11800.6 } ;
  double params_650[4] = {-2.15074e-06, 648.512, 0.0008, 23441.5 } ;
  double params_800[4] = {-1.18176e-06, 762.783, 0.0008, 36214.6 } ;
  double params_1000[4] = {-4.20059e-06, 812.788, 0.0008, 85283.8 } ;

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


}
