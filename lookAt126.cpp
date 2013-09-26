// c++ -o lookAt126 `root-config --glibs --cflags` `lhapdf-config --cppflags  --ldflags` -lm lookAt126.cpp

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cassert>

#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"
// CINT does not understand some files included by LorentzVector
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "LHAPDF/LHAPDF.h"


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
                          TString ("h_MWW_") + name, 3400., 100., 1800.) ;
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
fillHistos (LHEF::Reader & reader, histos & Histos, double XS, double referenceScale = -1, int max = -1)
{
  double totalCount = 0. ;
  int events = 0 ;

  int dynamicScale = 0 ;
  if (referenceScale < 0) dynamicScale = 1 ;
   
  //PG loop over input events
  while (reader.readEvent ()) 
    {
//      if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

      if (events % 10000 == 0) cout << "reading " << events << " event" << endl ;
          
      vector<TLorentzVector> v_f_Ws ;
      vector<TLorentzVector> v_f_quarks ;
      vector<TLorentzVector> v_f_leptons ;
      vector<TLorentzVector> v_f_neutrinos ;
      
      double x[2] = {0., 0.} ;
      int flavour[2] = {0, 0} ;

      // loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); ++iPart) 
        {
          TLorentzVector dummy = buildP (reader.hepeup, iPart) ;

           // incoming particles        
           if (reader.hepeup.ISTUP.at (iPart) == -1) 
             {
               x[iPart] = dummy.P () / 4000. ;
               flavour[iPart] = reader.hepeup.IDUP.at (iPart) ;
             } // incoming particles         

          // intermediate particles          
          if (reader.hepeup.ISTUP.at (iPart) == 2)
            {
              if (abs (reader.hepeup.IDUP.at (iPart)) == 24) 
                {
                  v_f_Ws.push_back (dummy) ;
                }              
            } // intermediate particles
            
          // outgoing particles          
          if (reader.hepeup.ISTUP.at (iPart) == 1)
            {
           // quarks
           if (abs (reader.hepeup.IDUP.at (iPart)) < 7) 
             {
               v_f_quarks.push_back (dummy) ;        
             } // quarks
           else if (abs (reader.hepeup.IDUP.at (iPart)) == 11 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 13 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 15)
             {
               v_f_leptons.push_back (dummy) ;
             }
           else if (abs (reader.hepeup.IDUP.at (iPart)) == 12 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 14 ||
                    abs (reader.hepeup.IDUP.at (iPart)) == 16)
             {
               v_f_neutrinos.push_back (dummy) ;        
             }
         } // outgoing particles
       } // loop over particles in the event

//      if (totalCount < 10) cout << "PARTICLES " <<  v_f_leptons.size () << "\t" << v_f_neutrinos.size () << "\t" << v_f_quarks.size () << "\n" ;

      int warnNum = 0 ;
      if (v_f_quarks.size () < 4)
        {
          cout << "warning, not enough quarks" << endl ;
          ++warnNum ;
        }
      if (v_f_leptons.size () < 1)
        {
          cout << "warning, not enough leptons" << endl ;
          ++warnNum ;
        }
      if (v_f_neutrinos.size () < 1)
        {
          cout << "warning, not enough neutrinos" << endl ;
          ++warnNum ;
        }
      if (warnNum > 0) continue ;

      double weight = 1. ;
      float scale = reader.hepeup.SCALUP ;
      if (dynamicScale > 0)
        {
          referenceScale = 80.385 * 80.385 + 
              (v_f_quarks.at (0).Pt () * v_f_quarks.at (0).Pt () +
               v_f_quarks.at (1).Pt () * v_f_quarks.at (1).Pt () +
               v_f_quarks.at (2).Pt () * v_f_quarks.at (2).Pt () +
               v_f_quarks.at (3).Pt () * v_f_quarks.at (3).Pt () +
               v_f_leptons.at (0).Pt () * v_f_leptons.at (0).Pt () +
               v_f_neutrinos.at (0).Pt () * v_f_neutrinos.at (0).Pt ()) / 6. ;
          referenceScale = sqrt (referenceScale) ;
        }
      
      if (referenceScale != 0 )
        weight = LHAPDF::xfx (x[0], referenceScale, flavour[0]) * LHAPDF::xfx (x[1], referenceScale, flavour[1]) /
                 (LHAPDF::xfx (x[0], scale, flavour[0]) * LHAPDF::xfx (x[1], scale, flavour[1])) ;
      totalCount += weight ;

      //PG the scale:
      Histos.m_h_scale->Fill (referenceScale) ;

      pair<int, int> detaIndices = findPairWithLargestDeta (v_f_quarks) ;
      if (v_f_quarks.at (detaIndices.second).Eta () - v_f_quarks.at (detaIndices.first).Eta () < 2) continue ;
      TLorentzVector largestPair = v_f_quarks.at (detaIndices.second) + v_f_quarks.at (detaIndices.first) ;
      if (largestPair.M () < 100) continue ; //PG selection applied in phantom

      //PG do I need this cut?! FIXME
      int cont = 0 ;
      for (int iJ = 0 ; iJ < 4 ; ++iJ)
        for (int iJ2 = iJ + 1 ; iJ2 < 4 ; ++iJ2)
          if (v_f_quarks.at (iJ).DeltaR (v_f_quarks.at (iJ2)) < 0.4) cont = 1 ;
      if (cont == 1) continue ;

      //PG the first two are the VBF jets, the following ones the W jets
      sort (v_f_quarks.rbegin (), v_f_quarks.rend (), ptsort ()) ;  
      
//      pair<int, int> Wpair (2, 3) ;
      pair<int, int> Wpair = findPairWithWMass (v_f_quarks) ;

      if (Wpair.first > 3 || Wpair.second > 3)
        {
          cout << "warning, wrong quarks in W determination\n" ;
        }

      TLorentzVector total = (v_f_leptons.at (0) + v_f_neutrinos.at (0)) + 
                             (v_f_quarks.at (Wpair.first) + v_f_quarks.at (Wpair.second)) ;

      Histos.m_h_MWW->Fill (total.M (), weight) ;
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


int main (int argc, char ** argv) 
{

  const int SUBSET = 0 ;
  const string NAME = "cteq6ll" ; //"cteq6l1"

  LHAPDF::initPDFSet (NAME, LHAPDF::LHPDF, SUBSET) ;
  const int NUMBER = LHAPDF::numberPDF () ;

  LHAPDF::initPDF (0) ;

  int maxEventsPerSample = -1 ;
  if (argc >= 2)
    {
      int dummy = atoi (argv[1]) ;
      if (dummy > 0) maxEventsPerSample = dummy ;
    }

  int dynamicScale = 0 ; 
  if (argc == 3)
    {
      int dummy = atoi (argv[2]) ;
      if (dummy > 0) dynamicScale = 1 ;
    }

  cout << "use dynamic scale: " << dynamicScale << endl ;
  double mass = 126 ;

  //PG choose the samples
  //PG ---- ---- ---- ---- ---- ---- ---- ---- ----

  // pietro's sample
  string filename_phbkg = "/Users/govoni/data/lvjj_samples/interference/phantom/total.126.lhe" ;
  double XS_phbkg = 0.07756069 * 2 ; // 7.7560687011E-002 // pb the factor 2 accounts for muons, electrons
  string filename_mg = "/Users/govoni/data/lvjj_samples/interference/madgraph/madgraph_mH126.lhe" ;
  double XS_mg =  0.0077079 * 2 ;
                    
  //PG messages
  
  cout << "\nworking with mass : " << mass << endl ;
  cout << "sample PH     :\t" << filename_phbkg << endl ;
  cout << "sample MG     :\t" << filename_mg << endl ;

  //PG ---- phantom ---- background only

  std::ifstream ifs_phbkg (filename_phbkg.c_str ()) ;
  LHEF::Reader reader_phbkg (ifs_phbkg) ;
  histos H_phbkg ("phbkg", XS_phbkg) ;

  std::ifstream ifs_mg (filename_mg.c_str ()) ;
  LHEF::Reader reader_mg (ifs_mg) ;
  histos H_mg ("mg", XS_mg) ;
  
  double scale = mass ;
  if (dynamicScale == 1) scale = -1 ;
  double entries_phbkg = fillHistos (reader_phbkg, H_phbkg, XS_phbkg, scale, maxEventsPerSample) ;
  cout << "total number of PH events : " << entries_phbkg << endl ;
  double entries_mg    = fillHistos (reader_mg   , H_mg   , XS_mg   , scale, maxEventsPerSample) ;
  cout << "total number of MG events : " << entries_mg << endl ;


  //PG saving the histograms

  TString name = "lookAt126.root" ;
  TFile f (name, "recreate") ;
  H_phbkg.save (f) ;
  H_mg.save (f) ;
  f.Close () ;


  return 0 ;
}
