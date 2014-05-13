////open parameter file and interpolate to find the value of the interference weigth

//compile with: 
//g++ -o prova.o read_and_interpolate.cpp -c `root-config --cflags`;
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

double ratio_crystalBallLowHighWithRise (double* x, double* par)
{
  double den = crystalBallLowHigh (x, par + 9) ; // signal only
  if (den == 0) return -1. ;
  double num = doubleGausCrystalBallLowHighPlusExp (x, par) ;    // signal and interference
  return num / den ;
}

//////////////////////////////////////

int main(int argc, char *argv[])
{

  TString parameters_normal [9] = {"Norm","Mean_CB","Sigma_CB","alphaR_CB","nR_CB","alphaL_CB","nL_CB","R","Tau"};

  double mass_chosen[5] = {600,700,800,900,1000};
  double c_chosen[6] = {0.100001,0.3,0.500001,0.7,0.900001,0.999};
  int Npar=9;

  TString *readfile = new TString ("file_for_interpolation.root"); //file with the values of the all parameters
  TFile* SI = new TFile(readfile->Data());
  Double_t fill_param[16];

  ofstream grid;
  grid.open ("final_grid.txt");  //write interference weight values on file

  TGraph2D *igraph_SI[9];
  TGraph2D *igraph_S[7];

  for (int i=0; i<Npar; i++) {

	TString *name = new TString (parameters_normal[i]);  
	name->Append("_SI.txt");
	igraph_SI[i] = (TGraph2D*)SI->Get(name->Data());
  }

  for (int i=0; i<Npar-2; i++) {

	TString *name = new TString (parameters_normal[i]);  
	name->Append("_S.txt");
	igraph_S[i] = (TGraph2D*)SI->Get(name->Data());
   }


  for (int m=0; m<5; m++) {
    for (int c=0; c<6; c++) {

      for (int i=0; i<Npar; i++) {

	cout<<m<<" "<<c<<" "<<i<<endl;
	fill_param[i]=igraph_SI[i]->Interpolate(mass_chosen[m],c_chosen[c]);
      }
      for (int i=0; i<Npar-2; i++) {
	cout<<m<<" "<<c<<" "<<i+9<<endl;
	fill_param[i+9]=igraph_S[i]->Interpolate(mass_chosen[m],c_chosen[c]);
      }

      //  for (int i=0; i<16; i++)   std::cout<<fill_param[i]<<std::endl;

      TF1* f_ratio23 = new TF1 ("f_ratio23", ratio_crystalBallLowHighWithRise, 0, 2000, 16) ;
      f_ratio23->SetParameters (fill_param) ;

      std::cout<<f_ratio23->Eval(mass_chosen[m])<<std::endl;
      grid<<mass_chosen[m]<<" "<<c_chosen[c]<<" "<<f_ratio23->Eval(mass_chosen[m])<<"\n";
    }
  }

  grid.close();
  SI->Close();
  

  return(0);

}