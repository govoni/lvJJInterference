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
#include "TH2D.h"
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
  cout<<"den: "<<den<<endl;
  if (den == 0) return -1. ;
  double num = doubleGausCrystalBallLowHighPlusExp (x, par) ;    // signal and interference
    return num / den ;
    // return num;
}

//////////////////////////////////////

int main(int argc, char *argv[])
{

  TString parameters_normal [9] = {"Norm","Mean_CB","Sigma_CB","alphaR_CB","nR_CB","alphaL_CB","nL_CB","R","Tau"};

      double mass_chosen[5] = {600,700,800,900,1000};
    double c_chosen[6] = {0.1,0.3,0.5,0.7,0.9,1.0};
  //   double mass_chosen[1] = {800};
  //  double c_chosen[1] = {1.0};

  int Npar=9;

  TString *readfile = new TString ("file_for_interpolation.root"); //file with the values of the all parameters
  TFile* SI = new TFile(readfile->Data());
  Double_t fill_param[16];

  ofstream grid;
  grid.open ("final_grid.txt");  //write interference weight values on file

  TGraph2D *igraph_SI[9];
  TGraph2D *igraph_S[7];

  TH2D *ihisto_SI[9];
  TH2D *ihisto_S[7];

  for (int i=0; i<Npar; i++) {

	TString *name = new TString (parameters_normal[i]);  
	name->Append("_SI.txt");
	TString *name2 = new TString (parameters_normal[i]);  
	name2->Append("_SI_histo");

	igraph_SI[i] = (TGraph2D*)SI->Get(name->Data());
	//	ihisto_SI[i] = (TH2D*)SI->Get(name2->Data());
	//ihisto_SI[i]->SetDirectory(0);
	//	igraph_SI[i]->GetHistogram("empty");
	//	igraph_SI[i]->SetMaxIter(500000);
  }

  for (int i=0; i<Npar-2; i++) {

	TString *name = new TString (parameters_normal[i]);  
	name->Append("_S.txt");
	TString *name2 = new TString (parameters_normal[i]);  
	name2->Append("_S_histo");

	igraph_S[i] = (TGraph2D*)SI->Get(name->Data());
	//	igraph_S[i]->SetMaxIter(500000);

	//	ihisto_S[i] = (TH2D*)SI->Get(name2->Data());
	//	ihisto_S[i]->SetDirectory(0);
	//		igraph_S[i]->GetHistogram("empty");
   }


  for (int m=0; m<5; m++) {
    for (int c=0; c<6; c++) {

      for (int i=0; i<Npar; i++) {

	if (parameters_normal[i].Contains("Norm")) 
	    fill_param[i]=exp(igraph_SI[i]->Interpolate(mass_chosen[m],c_chosen[c]));
	else
	    fill_param[i]=igraph_SI[i]->Interpolate(mass_chosen[m],c_chosen[c]);
      }
      for (int i=0; i<Npar-2; i++) {
	//	cout<<m<<" "<<c<<" "<<i+9<<endl;
	if (parameters_normal[i].Contains("Norm"))
	    fill_param[i+9]=exp(igraph_S[i]->Interpolate(mass_chosen[m],c_chosen[c]));
	else if (parameters_normal[i].Contains("nL"))
	  if (igraph_S[i]->Interpolate(mass_chosen[m],c_chosen[c])>100)
	    fill_param[i+9]=100;
	  else
	    fill_param[i+9]=igraph_S[i]->Interpolate(mass_chosen[m],c_chosen[c]);

	else
	    fill_param[i+9]=igraph_S[i]->Interpolate(mass_chosen[m],c_chosen[c]);
      }
             if (m==0 && c==0) {
      	        for (int i=0; i<16; i++)   
      	  std::cout<<fill_param[i]<<std::endl;
            }

      TF1* f_ratio23 = new TF1 ("f_ratio23", ratio_crystalBallLowHighWithRise, 0, 2000, 16) ;
      f_ratio23->SetParameters (fill_param) ;

      //std::cout<<mass_chosen[m]<<" "<<c_chosen[c]<<" a 1000: "<<f_ratio23->Eval(588.131164551)<<std::endl;


	        grid<<mass_chosen[m]<<" "<<c_chosen[c]<<" "<<f_ratio23->Eval(mass_chosen[m])<<"\n";
    }
  }

  grid.close();
  SI->Close();
  

  return(0);

}
