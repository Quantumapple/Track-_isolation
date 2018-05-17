#define sw_treeMaker_cxx
#include "sw_treeMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>
#include <string>
using namespace std;

void sw_treeMaker::Loop(int nevent = 0, int debug = 0)
{
   if (fChain == 0) return;
   
   class track
   {
       public:
       float DR;
       float pos_x, pos_y, pos_z;
       track() { DR = 0.; pos_x = 0.; pos_y = 0.; pos_z = 0.; }
       track(float a, float b, float c, float d)
       {
	   DR = a; pos_x = b; pos_y = c; pos_z = d;
       }
       static bool comp(const track &t1, const track &t2)
       {
	   return ( t1.DR < t2.DR );
       }
   };

   TFile *file = new TFile("data.root","recreate");

   TTree* out_tree = new TTree("t","t");

   vector<float> bsl1_l1l2_dphi; out_tree->Branch("bsl1_l1l2_dphi",&bsl1_l1l2_dphi);
   vector<float> genpT;          out_tree->Branch("genpT",&genpT);
   
   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   if(nevent) nentries = nevent;
   cout << " nentries: " <<  nentries << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      cout << "Process: (" << jentry << "/" << nentries << ")" << endl;

      bsl1_l1l2_dphi.clear();
      genpT.clear();
      
      vector<track> v1; 
      vector<track> v2; 

      v1.clear();
      v2.clear();

      for(Int_t i = 0; i < bHitN; i++)
      {
	  if( fabs(bHitGz->at(i)) > 28.5 ) continue;
	  Float_t R = sqrt(pow(bHitGx->at(i),2)+pow(bHitGy->at(i),2));
	  TVector3 pixel;
	  pixel.SetXYZ( bHitGx->at(i) - simVx->at(0), bHitGy->at(i) - simVy->at(0), bHitGz->at(i) - simVz->at(0));
	  Float_t pixelPhi = pixel.Phi();
	  Float_t pixelEta = pixel.PseudoRapidity();
	  Float_t deltaPhi = pixelPhi - genPartPhi->at(0);
	  if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
	  if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
	  Float_t deltaEta = pixelEta - genPartEta->at(0);
	  Float_t deltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta,2));

	  if( R < 5.5 )            v1.push_back(track(deltaR,bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
	  if( R > 5.5 && R < 8.5 ) v2.push_back(track(deltaR,bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
      }

      sort(v1.begin(), v1.end(), track::comp); 
      sort(v2.begin(), v2.end(), track::comp); 
      
      file->cd();
      
      if( v1.size() != 0 && v2.size() != 0 ) // PVL1 - L1L2
      {
	  TVector3 pixel1;
	  pixel1.SetXYZ( v1[0].pos_x - simVx->at(0), v1[0].pos_y - simVy->at(0), v1[0].pos_z - simVz->at(0) );
	  
	  TVector3 pixel2;
	  pixel2.SetXYZ( v2[0].pos_x - v1[0].pos_x, v2[0].pos_y - v1[0].pos_y, v2[0].pos_z - v1[0].pos_z );

	  Float_t phi1 = pixel1.Phi();
	  Float_t phi2 = pixel2.Phi();
	  Float_t deltaPhi = phi1 - phi2;
	  if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
	  if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl1_l1l2_dphi.push_back(deltaPhi);
          //genpT.push_back(propgenPartPt->at(0));
          genpT.push_back(genPartPt->at(0));
      }

      out_tree->Fill(); 
 
       

   } // end of event loop

   file->Write();
}

