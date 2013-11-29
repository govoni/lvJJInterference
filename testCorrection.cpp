// c++ -o testCorrection `root-config --glibs --cflags` `lhapdf-config --cppflags  --ldflags`  -lGenVector -lm testCorrection.cpp


#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cassert>
#include <fstream>

#include "TString.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
// #include "TLorentzVector.h"
// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math ;
using namespace std ;

#include "findInterferenceTools.h"


struct histos
{
  TH1F *  m_h_MWW ;
  TH1F *  m_h_scale ;
  TString m_name ;
  double  m_XS ;
  
  histos (TString name, double XS) : m_name (name), m_XS (XS)
    {
      m_h_MWW = new TH1F (TString ("h_MWW_") + name, 
                          TString ("h_MWW_") + name, 800., 200., 1800.) ;
      m_h_MWW->Sumw2 () ;
      m_h_scale = new TH1F (TString ("h_scale_") + name, 
                          TString ("h_scale_") + name, 100, 0., 1000.) ;
      m_h_scale->Sumw2 () ;
    }
 
  void norm (double total = 0)
    {
      double factor = m_XS / m_h_MWW->GetEntries () ;
      if (total != 0) factor = m_XS / total ;
      m_h_MWW->Scale (factor) ;
      m_h_scale->Scale (factor) ;
    }
  
  ~histos ()
    {
      delete m_h_MWW ;
      delete m_h_scale ;
    }  
    
  void save (TFile & outfile) 
    {
      outfile.cd () ;
      m_h_MWW->Write () ;
      m_h_scale->Write () ;
    }
  

} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


double 
fillHistos (LHEF::Reader & reader, histos & Histos, double XS, double referenceScale = 0, int max = -1)
{
  double totalCount = 0. ;
  int events = 0 ;
   
  //PG loop over input events
  while (reader.readEvent ()) 
    {
//      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      if (events % 10000 == 0) cout << "reading " << events << " event" << endl ;
          
      vector<lorentzVector> v_f_Hs ;
      vector<lorentzVector> v_f_quarks ;
      
      double x[2] = {0., 0.} ;
      int flavour[2] = {0, 0} ;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
          lorentzVector dummy = buildLP (reader.hepeup, iPart) ;

           // incoming particles   
           int iInc = 0 ;     
           if (reader.hepeup.ISTUP.at (iPart) == -1) 
             {
               x[iInc] = dummy.P () / 4000. ;
               flavour[iInc++] = reader.hepeup.IDUP.at (iPart) ;
             } // incoming particles         

          // outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
           // quarks
           if (abs (reader.hepeup.IDUP.at (iPart)) < 7) 
             {
               v_f_quarks.push_back (dummy) ;        
             } // quarks
           else if (abs (reader.hepeup.IDUP.at (iPart)) == 25)
             {
               v_f_Hs.push_back (dummy) ;
             }
         } // outgoing particles
       } // loop over particles in the event

//      if (totalCount < 10) cout << "PARTICLES " <<  v_f_leptons.size () << "\t" << v_f_neutrinos.size () << "\t" << v_f_quarks.size () << "\n" ;

      double weight = 1. ;
      totalCount += weight ;

      int warnNum = 0 ;
      if (v_f_quarks.size () < 2)
        {
          cout << "warning, not enough quarks" << endl ;
          ++warnNum ;
        }
      if (v_f_Hs.size () < 1)
        {
          cout << "warning, not enough Higgs bosons" << endl ;
          ++warnNum ;
        }
      if (warnNum > 0) continue ;


      //PG apply all the production cuts from phantom and madgraph to the sample
      //PG ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
      
      int cont = 0 ;
      for (int i = 0 ; i < 2 ; ++i) 
        if (v_f_quarks.at (i).Pt () > 20 && 
            v_f_quarks.at (i).E () > 20 && 
            fabs (v_f_quarks.at (i).Eta ()) < 6.5) cont += 1 ;
      if (cont < 2) continue ;
      
      pair<int, int> detaIndices = findPairWithLargestDeta (v_f_quarks) ;
      if (v_f_quarks.at (detaIndices.second).Eta () - v_f_quarks.at (detaIndices.first).Eta () < 2) continue ;
      lorentzVector largestPair = v_f_quarks.at (detaIndices.second) + v_f_quarks.at (detaIndices.first) ;
      if (largestPair.M () < 100) continue ; //PG selection applied in phantom

      Histos.m_h_MWW->Fill (v_f_Hs[0].M (), weight) ;
      ++events ;
      if (max > 0 && max < events) 
        {
          cout << max << " events reached, exiting" << endl ;
          break ;
        }

    } //PG loop over input events

  Histos.norm (totalCount) ;

  return totalCount ;
  
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


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


struct corrections
{
  TGraph ** tg_sig_par ;
  TGraph ** tg_sAi_par ;

  string m_path ;
  corrections (string path) : m_path (path)
    {
      tg_sig_par = new TGraph * [7] ; // [parameter][mass]
      tg_sAi_par = new TGraph * [7] ;
      for (int k = 0 ; k < 7 ; ++k) tg_sig_par[k] = new TGraph (5) ;
      for (int k = 0 ; k < 7 ; ++k) tg_sAi_par[k] = new TGraph (5) ;

      //PG fill the graphs from the input files
      //PG ---- ---- ---- ---- ---- ---- ---- ----
    
      int i = 0 ; 
    
      intCont ic_350  (m_path + string ("/results_interference.350.root"), 350) ;      
      ic_350.setsigGraphs (tg_sig_par, i) ;
      ic_350.setsAiGraphs (tg_sAi_par, i++) ;
      intCont ic_500  (m_path + string ("/results_interference.500.root"), 500) ;      
      ic_500.setsigGraphs (tg_sig_par, i) ;  
      ic_500.setsAiGraphs (tg_sAi_par, i++) ;  
      intCont ic_650  (m_path + string ("/results_interference.650.root"), 650) ;      
      ic_650.setsigGraphs (tg_sig_par, i) ;  
      ic_650.setsAiGraphs (tg_sAi_par, i++) ;  
      intCont ic_800  (m_path + string ("/results_interference.800.root"), 800) ;      
      ic_800.setsigGraphs (tg_sig_par, i) ;  
      ic_800.setsAiGraphs (tg_sAi_par, i++) ;  
      intCont ic_1000 (m_path + string ("/results_interference.1000.root"), 1000) ;   
      ic_1000.setsigGraphs (tg_sig_par, i) ;   
      ic_1000.setsAiGraphs (tg_sAi_par, i++) ;   


    }


} ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int main (int argc, char ** argv) 
{


  if (argc < 5) 
    {
      cout << "usage " << argv[0] << " mass signalFile signalXS correctionsLocation" << endl ;
      exit (1) ;
    }

  double mass = atof (argv[1]) ;
  string signalFile = argv[2] ;
  double XS_sig = atof (argv[3]) ;
  string correctionsFolder = argv[4] ;

  //PG read the signal sample
  
  std::ifstream ifs_sig (signalFile.c_str ()) ;
  LHEF::Reader reader_sig (ifs_sig) ;
  histos H_sig ("mg", XS_sig) ;
  double entries_sig = fillHistos (reader_sig, H_sig, XS_sig, mass) ;
  TFile outfile ("testCorrection.root", "recreate") ;
  H_sig.save (outfile) ;
  outfile.Close () ;

  cout << "signal events : " << entries_sig << endl ;

  //PG read the correction factors
  corrections corrTool (correctionsFolder) ;
  

  return 0 ;
}
