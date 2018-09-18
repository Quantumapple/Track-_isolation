#define sw_treeMaker_cxx
#include "sw_treeMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>
#include <string>
using namespace std;

void sw_treeMaker::Loop()
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

   vector<float> bsl1_l1l3_dphi; out_tree->Branch("bsl1_l1l3_dphi",&bsl1_l1l3_dphi);
   vector<float> bsl1_l1l4_dphi; out_tree->Branch("bsl1_l1l4_dphi",&bsl1_l1l4_dphi);
   vector<float> bsl2_l2l3_dphi; out_tree->Branch("bsl2_l2l3_dphi",&bsl2_l2l3_dphi);
   vector<float> bsl2_l2l4_dphi; out_tree->Branch("bsl2_l2l4_dphi",&bsl2_l2l4_dphi);
   vector<float> bsl3_l3l4_dphi; out_tree->Branch("bsl3_l3l4_dphi",&bsl3_l3l4_dphi);
   
   vector<float> genpT_case1;    out_tree->Branch("genpT_case1",&genpT_case1);
   vector<float> genpT_case2;    out_tree->Branch("genpT_case2",&genpT_case2);
   vector<float> genpT_case3;    out_tree->Branch("genpT_case3",&genpT_case3);
   vector<float> genpT_case4;    out_tree->Branch("genpT_case4",&genpT_case4);
   vector<float> genpT_case5;    out_tree->Branch("genpT_case5",&genpT_case5);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      cout << "Process: (" << jentry << "/" << nentries << ")" << endl;

      bsl1_l1l3_dphi.clear(); 
      bsl1_l1l4_dphi.clear();
      bsl2_l2l3_dphi.clear();
      bsl2_l2l4_dphi.clear();
      bsl3_l3l4_dphi.clear();
   
      genpT_case1.clear();
      genpT_case2.clear();
      genpT_case3.clear();
      genpT_case4.clear();
      genpT_case5.clear();
      
      vector<track> v1; 
      vector<track> v2; 
      vector<track> v3; 
      vector<track> v4; 

      v1.clear();
      v2.clear();
      v3.clear();
      v4.clear();

      for(Int_t i = 0; i < bRecHitN; i++)
      {
          Float_t R = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
          TVector3 pixel;
          pixel.SetXYZ( bRecHitGx->at(i) - simVx->at(0), bRecHitGy->at(i) - simVy->at(0), bRecHitGz->at(i) - simVz->at(0));
          Float_t pixelPhi = pixel.Phi();
          Float_t pixelEta = pixel.PseudoRapidity();
          Float_t deltaPhi = pixelPhi - propgenElPartPhi->at(0);
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
          Float_t deltaEta = pixelEta - propgenElPartEta->at(0);
          Float_t deltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta,2));

          if( R < 5.5 )            v1.push_back(track(deltaR,bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          if( R > 5.5 && R < 8.5 ) v2.push_back(track(deltaR,bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          if( R > 8.5 && R < 13. ) v3.push_back(track(deltaR,bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          if( R > 13. && R < 18. ) v4.push_back(track(deltaR,bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
      }

      sort(v1.begin(), v1.end(), track::comp); 
      sort(v2.begin(), v2.end(), track::comp); 
      sort(v3.begin(), v3.end(), track::comp); 
      sort(v4.begin(), v4.end(), track::comp); 
      
      file->cd();
      
      if( v1.size() != 0 && v3.size() != 0 ) // PVL1 - L1L3
      {
          TVector3 pixel1;
          pixel1.SetXYZ( v1[0].pos_x - simVx->at(0), v1[0].pos_y - simVy->at(0), v1[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( v3[0].pos_x - v1[0].pos_x, v3[0].pos_y - v1[0].pos_y, v3[0].pos_z - v1[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl1_l1l3_dphi.push_back(deltaPhi);
          genpT_case1.push_back(genPartPt->at(0));
      }
      
      if( v1.size() != 0 && v4.size() != 0 ) // PVL1 - L1L3
      {
          TVector3 pixel1;
          pixel1.SetXYZ( v1[0].pos_x - simVx->at(0), v1[0].pos_y - simVy->at(0), v1[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( v4[0].pos_x - v1[0].pos_x, v4[0].pos_y - v1[0].pos_y, v4[0].pos_z - v1[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl1_l1l4_dphi.push_back(deltaPhi);
          genpT_case2.push_back(genPartPt->at(0));
      }
      
      if( v2.size() != 0 && v3.size() != 0 ) // PVL2 - L2L3
      {
          TVector3 pixel1;
          pixel1.SetXYZ( v2[0].pos_x - simVx->at(0), v2[0].pos_y - simVy->at(0), v2[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( v3[0].pos_x - v2[0].pos_x, v3[0].pos_y - v2[0].pos_y, v3[0].pos_z - v2[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl2_l2l3_dphi.push_back(deltaPhi);
          genpT_case3.push_back(genPartPt->at(0));
      }
      
      if( v2.size() != 0 && v4.size() != 0 ) // PVL2 - L2L4
      {
          TVector3 pixel1;
          pixel1.SetXYZ( v2[0].pos_x - simVx->at(0), v2[0].pos_y - simVy->at(0), v2[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( v4[0].pos_x - v2[0].pos_x, v4[0].pos_y - v2[0].pos_y, v4[0].pos_z - v2[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl2_l2l4_dphi.push_back(deltaPhi);
          genpT_case4.push_back(genPartPt->at(0));
      }
      
      if( v3.size() != 0 && v4.size() != 0 ) // PVL3 - L3L4
      {
          TVector3 pixel1;
          pixel1.SetXYZ( v3[0].pos_x - simVx->at(0), v3[0].pos_y - simVy->at(0), v3[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( v4[0].pos_x - v3[0].pos_x, v4[0].pos_y - v3[0].pos_y, v4[0].pos_z - v3[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl3_l3l4_dphi.push_back(deltaPhi);
          genpT_case5.push_back(genPartPt->at(0));
      }
      
       
      out_tree->Fill(); 
   
   } // event loop

   file->Write();   

}
