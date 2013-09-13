// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//PG DA FARE MEGLIO
//PG ---------------

Double_t doublePeakModel (Double_t * xx, Double_t * par)
{
  double scale    = par[0] ;
  double shift    = par[1] ;
  double distance = par[2] ;
  double gamma    = par[3] ;
  double x = xx[0] - shift ;

  double norm = 1. / (shift * shift + gamma) - 1 / ((shift + 2 * distance) * (shift + 2 * distance) + gamma) ;
  return scale * (1. / norm) * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;
  return scale * ( 1. / ((x - distance) * (x - distance) + gamma) - 1 / ((x + distance) * (x + distance) + gamma)) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


// funzioni per i parametri dell'interferenza
double getIntPar (int par, double mass)
{
  if (par == 0) return -2e-06 ;
  if (par == 1) return -106.461 + 1.4656 * mass - 0.000496861 * mass * mass ;
  if (par == 2) return -0.002 ;
  if (par == 3) return 3093.54 + 52864.9 / (1.0 + TMath::Exp (-0.0203282 * (mass - 676.986))) ;
  cerr << "getIntPar wrong par ID: " << par << endl ;
  return -1 ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


// funzioni per i parametri del segnale
double getSigPar (int par, double mass)
{
  if (par == 0) return TMath::Exp (-0.88736 - 0.011128 * mass) ;
  if (par == 1) return -115.928 + 1.45672 * mass - 0.000398279 * mass * mass ;
  if (par == 2) return 21.1992 - 0.165296 * mass + 0.000361834 * mass * mass ;
  if (par == 3) return 0.8 + (mass - 350) * (1.6 - 0.8) / (1000 - 350) ;
  if (par == 4) return -1.5 + (mass - 350) * (3.37 - 1.5) / (1000 - 350) ;
  if (par == 5) return 1.25 ;
  if (par == 6) return -2.15 + (mass - 350) * (7 - 2.15) / (1000 - 350) ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int macro_interpolInterference ()
{
  
  TCanvas * c1 = new TCanvas () ;
  c1->DrawFrame (200, -0.05, 2000, 0.05) ;

  for (int mass = 350 ; mass < 1001 ;mass += 200)
    {
      cout << mass << endl ;
      TF1 * func = new TF1 (TString ("func_") + TString (mass), doublePeakModel, 200, 2000, 4) ;
      for (int j = 0 ; j < 4 ; ++j) func->SetParameter (j, getIntPar (j, mass)) ;
      func->Draw ("same") ;
    }
}