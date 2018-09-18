#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>
#include <string>

using namespace std;

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   //TFile *file = new TFile("data.root","recreate");
   TFile *file = new TFile("data_bsd2_d2d3.root","recreate");

   TTree* out_tree = new TTree("t","t");

   //vector<float> bsl1_l1d1_dphi;  out_tree->Branch("bsl1_l1d1_dphi",&bsl1_l1d1_dphi);
   //vector<float> bsl1_l1d2_dphi;  out_tree->Branch("bsl1_l1d2_dphi",&bsl1_l1d2_dphi);
   //vector<float> bsl2_l2d1_dphi;  out_tree->Branch("bsl2_l2d1_dphi",&bsl2_l2d1_dphi);
   //vector<float> bsl2_l2d2_dphi;  out_tree->Branch("bsl2_l2d2_dphi",&bsl2_l2d2_dphi);
   //vector<float> bsl3_l3d1_dphi;  out_tree->Branch("bsl3_l3d1_dphi",&bsl3_l3d1_dphi);
   //vector<float> bsd1_d1d2_dphi;  out_tree->Branch("bsd1_d1d2_dphi",&bsd1_d1d2_dphi);
   //vector<float> bsd1_d1d3_dphi;  out_tree->Branch("bsd1_d1d3_dphi",&bsd1_d1d3_dphi);
   //vector<float> bsd1_d1d4_dphi;  out_tree->Branch("bsd1_d1d4_dphi",&bsd1_d1d4_dphi);
   vector<float> bsd2_d2d3_dphi;  out_tree->Branch("bsd2_d2d3_dphi",&bsd2_d2d3_dphi);
   //vector<float> bsd2_d2d4_dphi;  out_tree->Branch("bsd2_d2d4_dphi",&bsd2_d2d4_dphi);
   //vector<float> bsd2_d2d5_dphi;  out_tree->Branch("bsd2_d2d5_dphi",&bsd2_d2d5_dphi);
   //vector<float> bsd3_d3d4_dphi;  out_tree->Branch("bsd3_d3d4_dphi",&bsd3_d3d4_dphi);
   //vector<float> bsd3_d3d5_dphi;  out_tree->Branch("bsd3_d3d5_dphi",&bsd3_d3d5_dphi);
   //vector<float> bsd3_d3d6_dphi;  out_tree->Branch("bsd3_d3d6_dphi",&bsd3_d3d6_dphi);
   //vector<float> bsd4_d4d6_dphi;  out_tree->Branch("bsd4_d4d6_dphi",&bsd4_d4d6_dphi);
   
   //vector<float> genpT_bsl1_l1d1; out_tree->Branch("genpT_bsl1_l1d1",&genpT_bsl1_l1d1);
   //vector<float> genpT_bsl1_l1d2; out_tree->Branch("genpT_bsl1_l1d2",&genpT_bsl1_l1d2);
   //vector<float> genpT_bsl2_l2d1; out_tree->Branch("genpT_bsl2_l2d1",&genpT_bsl2_l2d1);
   //vector<float> genpT_bsl2_l2d2; out_tree->Branch("genpT_bsl2_l2d2",&genpT_bsl2_l2d2);
   //vector<float> genpT_bsl3_l3d1; out_tree->Branch("genpT_bsl3_l3d1",&genpT_bsl3_l3d1);
   //vector<float> genpT_bsd1_d1d2; out_tree->Branch("genpT_bsd1_d1d2",&genpT_bsd1_d1d2);
   //vector<float> genpT_bsd1_d1d3; out_tree->Branch("genpT_bsd1_d1d3",&genpT_bsd1_d1d3);
   //vector<float> genpT_bsd1_d1d4; out_tree->Branch("genpT_bsd1_d1d4",&genpT_bsd1_d1d4);
   vector<float> genpT_bsd2_d2d3; out_tree->Branch("genpT_bsd2_d2d3",&genpT_bsd2_d2d3);
   //vector<float> genpT_bsd2_d2d4; out_tree->Branch("genpT_bsd2_d2d4",&genpT_bsd2_d2d4);
   //vector<float> genpT_bsd2_d2d5; out_tree->Branch("genpT_bsd2_d2d5",&genpT_bsd2_d2d5);
   //vector<float> genpT_bsd3_d3d4; out_tree->Branch("genpT_bsd3_d3d4",&genpT_bsd3_d3d4);
   //vector<float> genpT_bsd3_d3d5; out_tree->Branch("genpT_bsd3_d3d5",&genpT_bsd3_d3d5);
   //vector<float> genpT_bsd3_d3d6; out_tree->Branch("genpT_bsd3_d3d6",&genpT_bsd3_d3d6);
   //vector<float> genpT_bsd4_d4d6; out_tree->Branch("genpT_bsd4_d4d6",&genpT_bsd4_d4d6);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<100000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      cout << "Process: " << jentry+1 << endl;

      //bsl1_l1d1_dphi.clear();  
      //bsl1_l1d2_dphi.clear();  
      //bsl2_l2d1_dphi.clear();  
      //bsl2_l2d2_dphi.clear();  
      //bsl3_l3d1_dphi.clear();  
      //bsd1_d1d2_dphi.clear();  
      //bsd1_d1d3_dphi.clear();  
      //bsd1_d1d4_dphi.clear();  
      bsd2_d2d3_dphi.clear();  
      //bsd2_d2d4_dphi.clear();  
      //bsd2_d2d5_dphi.clear();  
      //bsd3_d3d4_dphi.clear();  
      //bsd3_d3d5_dphi.clear();  
      //bsd3_d3d6_dphi.clear();  
      //bsd4_d4d6_dphi.clear();  
   
      //genpT_bsl1_l1d1.clear(); 
      //genpT_bsl1_l1d2.clear(); 
      //genpT_bsl2_l2d1.clear(); 
      //genpT_bsl2_l2d2.clear(); 
      //genpT_bsl3_l3d1.clear(); 
      //genpT_bsd1_d1d2.clear(); 
      //genpT_bsd1_d1d3.clear(); 
      //genpT_bsd1_d1d4.clear(); 
      genpT_bsd2_d2d3.clear(); 
      //genpT_bsd2_d2d4.clear(); 
      //genpT_bsd2_d2d5.clear(); 
      //genpT_bsd3_d3d4.clear(); 
      //genpT_bsd3_d3d5.clear(); 
      //genpT_bsd3_d3d6.clear(); 
      //genpT_bsd4_d4d6.clear(); 
      
      Float_t propGenPhi  = propgenElPartPhi->at(0);
      Float_t propGenEta  = propgenElPartEta->at(0);
      
      bLayer1.clear(); 
      bLayer2.clear(); 
      bLayer3.clear();
      bLayer4.clear();
      fDisk1.clear(); 
      fDisk2.clear(); 
      fDisk3.clear(); 
      fDisk4.clear(); 
      fDisk5.clear(); 
      fDisk6.clear();
      
      // Store pixel clusters w.r.t. eta range
      StorePixelHits(propGenPhi, propGenEta); 

      sort(bLayer1.begin(), bLayer1.end(), track::comp);
      sort(bLayer2.begin(), bLayer2.end(), track::comp);
      sort(bLayer3.begin(), bLayer3.end(), track::comp);
      sort(bLayer4.begin(), bLayer4.end(), track::comp);

      sort(fDisk1.begin(), fDisk1.end(), track::comp);
      sort(fDisk2.begin(), fDisk2.end(), track::comp);
      sort(fDisk3.begin(), fDisk3.end(), track::comp);
      sort(fDisk4.begin(), fDisk4.end(), track::comp);
      sort(fDisk5.begin(), fDisk5.end(), track::comp);
      sort(fDisk6.begin(), fDisk6.end(), track::comp);

      file->cd();
     
      /*
      if( bLayer1.size() != 0 && fDisk1.size() != 0 ) // PVL1 - L1D1
      {
          TVector3 pixel1;
          pixel1.SetXYZ( bLayer1[0].pos_x - simVx->at(0), bLayer1[0].pos_y - simVy->at(0), bLayer1[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk1[0].pos_x - bLayer1[0].pos_x, fDisk1[0].pos_y - bLayer1[0].pos_y, fDisk1[0].pos_z - bLayer1[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl1_l1d1_dphi.push_back(deltaPhi);
          genpT_bsl1_l1d1.push_back(genPartPt->at(0));
      }


      if( bLayer1.size() != 0 && fDisk2.size() != 0 ) // PVL1 - L1D2
      {
          TVector3 pixel1;
          pixel1.SetXYZ( bLayer1[0].pos_x - simVx->at(0), bLayer1[0].pos_y - simVy->at(0), bLayer1[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk2[0].pos_x - bLayer1[0].pos_x, fDisk2[0].pos_y - bLayer1[0].pos_y, fDisk2[0].pos_z - bLayer1[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl1_l1d2_dphi.push_back(deltaPhi);
          genpT_bsl1_l1d2.push_back(genPartPt->at(0));
      }

      if( bLayer2.size() != 0 && fDisk1.size() != 0 ) // PVL2 - L2D1
      {
          TVector3 pixel1;
          pixel1.SetXYZ( bLayer2[0].pos_x - simVx->at(0), bLayer2[0].pos_y - simVy->at(0), bLayer2[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk1[0].pos_x - bLayer2[0].pos_x, fDisk1[0].pos_y - bLayer2[0].pos_y, fDisk1[0].pos_z - bLayer2[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl2_l2d1_dphi.push_back(deltaPhi);
          genpT_bsl2_l2d1.push_back(genPartPt->at(0));
      }

      if( bLayer2.size() != 0 && fDisk2.size() != 0 ) // PVL2 - L2D2
      {
          TVector3 pixel1;
          pixel1.SetXYZ( bLayer2[0].pos_x - simVx->at(0), bLayer2[0].pos_y - simVy->at(0), bLayer2[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk2[0].pos_x - bLayer2[0].pos_x, fDisk2[0].pos_y - bLayer2[0].pos_y, fDisk2[0].pos_z - bLayer2[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl2_l2d2_dphi.push_back(deltaPhi);
          genpT_bsl2_l2d2.push_back(genPartPt->at(0));
      }
      

      if( bLayer3.size() != 0 && fDisk1.size() != 0 ) // PVL3 - L3D1
      {
          TVector3 pixel1;
          pixel1.SetXYZ( bLayer3[0].pos_x - simVx->at(0), bLayer3[0].pos_y - simVy->at(0), bLayer3[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk1[0].pos_x - bLayer3[0].pos_x, fDisk1[0].pos_y - bLayer3[0].pos_y, fDisk1[0].pos_z - bLayer3[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsl3_l3d1_dphi.push_back(deltaPhi);
          genpT_bsl3_l3d1.push_back(genPartPt->at(0));
      }

      
      if( fDisk1.size() != 0 && fDisk2.size() != 0 ) // PVD1 - D1D2
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk1[0].pos_x - simVx->at(0), fDisk1[0].pos_y - simVy->at(0), fDisk1[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk2[0].pos_x - fDisk1[0].pos_x, fDisk2[0].pos_y - fDisk1[0].pos_y, fDisk2[0].pos_z - fDisk1[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd1_d1d2_dphi.push_back(deltaPhi);
          genpT_bsd1_d1d2.push_back(genPartPt->at(0));
      }

      if( fDisk1.size() != 0 && fDisk3.size() != 0 ) // PVD1 - D1D3
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk1[0].pos_x - simVx->at(0), fDisk1[0].pos_y - simVy->at(0), fDisk1[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk3[0].pos_x - fDisk1[0].pos_x, fDisk3[0].pos_y - fDisk1[0].pos_y, fDisk3[0].pos_z - fDisk1[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd1_d1d3_dphi.push_back(deltaPhi);
          genpT_bsd1_d1d3.push_back(genPartPt->at(0));
      }
      
      if( fDisk1.size() != 0 && fDisk4.size() != 0 ) // PVD1 - D1D4
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk1[0].pos_x - simVx->at(0), fDisk1[0].pos_y - simVy->at(0), fDisk1[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk4[0].pos_x - fDisk1[0].pos_x, fDisk4[0].pos_y - fDisk1[0].pos_y, fDisk4[0].pos_z - fDisk1[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd1_d1d4_dphi.push_back(deltaPhi);
          genpT_bsd1_d1d4.push_back(genPartPt->at(0));
      }
      */
      
      if( fDisk2.size() != 0 && fDisk3.size() != 0 ) // PVD2 - D2D3
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk2[0].pos_x - simVx->at(0), fDisk2[0].pos_y - simVy->at(0), fDisk2[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk3[0].pos_x - fDisk2[0].pos_x, fDisk3[0].pos_y - fDisk2[0].pos_y, fDisk3[0].pos_z - fDisk2[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd2_d2d3_dphi.push_back(deltaPhi);
          genpT_bsd2_d2d3.push_back(genPartPt->at(0));
      }

      /*
      if( fDisk2.size() != 0 && fDisk4.size() != 0 ) // PVD2 - D2D4
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk2[0].pos_x - simVx->at(0), fDisk2[0].pos_y - simVy->at(0), fDisk2[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk4[0].pos_x - fDisk2[0].pos_x, fDisk4[0].pos_y - fDisk2[0].pos_y, fDisk4[0].pos_z - fDisk2[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd2_d2d4_dphi.push_back(deltaPhi);
          genpT_bsd2_d2d4.push_back(genPartPt->at(0));
      }

      if( fDisk2.size() != 0 && fDisk5.size() != 0 ) // PVD2 - D2D5
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk2[0].pos_x - simVx->at(0), fDisk2[0].pos_y - simVy->at(0), fDisk2[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk5[0].pos_x - fDisk2[0].pos_x, fDisk5[0].pos_y - fDisk2[0].pos_y, fDisk5[0].pos_z - fDisk2[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd2_d2d5_dphi.push_back(deltaPhi);
          genpT_bsd2_d2d5.push_back(genPartPt->at(0));
      }

      if( fDisk3.size() != 0 && fDisk4.size() != 0 ) // PVD3 - D3D4
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk3[0].pos_x - simVx->at(0), fDisk3[0].pos_y - simVy->at(0), fDisk3[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk4[0].pos_x - fDisk3[0].pos_x, fDisk4[0].pos_y - fDisk3[0].pos_y, fDisk4[0].pos_z - fDisk3[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd3_d3d4_dphi.push_back(deltaPhi);
          genpT_bsd3_d3d4.push_back(genPartPt->at(0));
      }
      
      if( fDisk3.size() != 0 && fDisk5.size() != 0 ) // PVD3 - D3D5
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk3[0].pos_x - simVx->at(0), fDisk3[0].pos_y - simVy->at(0), fDisk3[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk5[0].pos_x - fDisk3[0].pos_x, fDisk5[0].pos_y - fDisk3[0].pos_y, fDisk5[0].pos_z - fDisk3[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd3_d3d5_dphi.push_back(deltaPhi);
          genpT_bsd3_d3d5.push_back(genPartPt->at(0));
      }

      if( fDisk3.size() != 0 && fDisk6.size() != 0 ) // PVD3 - D3D6
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk3[0].pos_x - simVx->at(0), fDisk3[0].pos_y - simVy->at(0), fDisk3[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk6[0].pos_x - fDisk3[0].pos_x, fDisk6[0].pos_y - fDisk3[0].pos_y, fDisk6[0].pos_z - fDisk3[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd3_d3d6_dphi.push_back(deltaPhi);
          genpT_bsd3_d3d6.push_back(genPartPt->at(0));
      }

      if( fDisk4.size() != 0 && fDisk6.size() != 0 ) // PVD4 - D4D6
      {
          TVector3 pixel1;
          pixel1.SetXYZ( fDisk4[0].pos_x - simVx->at(0), fDisk4[0].pos_y - simVy->at(0), fDisk4[0].pos_z - simVz->at(0) );

          TVector3 pixel2;
          pixel2.SetXYZ( fDisk6[0].pos_x - fDisk4[0].pos_x, fDisk6[0].pos_y - fDisk4[0].pos_y, fDisk6[0].pos_z - fDisk4[0].pos_z );

          Float_t phi1 = pixel1.Phi();
          Float_t phi2 = pixel2.Phi();
          Float_t deltaPhi = phi1 - phi2;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

          bsd4_d4d6_dphi.push_back(deltaPhi);
          genpT_bsd4_d4d6.push_back(genPartPt->at(0));
      }
      */
      
      out_tree->Fill(); 
   
   } // event loop
   
   file->Write();   
   
}
