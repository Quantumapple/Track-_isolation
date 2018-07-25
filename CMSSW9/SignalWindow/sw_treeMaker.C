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


   bool debug = false;

   TFile *file = new TFile("R1_ROI.root","recreate");

   TTree* out_tree = new TTree("t","t");

   float abs_pterr; out_tree->Branch("abs_pterr", &abs_pterr, "abs_pterr/F");
   float closest_egEt; out_tree->Branch("closest_egEt", &closest_egEt, "closest_egEt/F");
   float closest_egEta; out_tree->Branch("closest_egEta", &closest_egEta, "closest_egEta/F");
   float closest_egPhi; out_tree->Branch("closest_egPhi", &closest_egPhi, "closest_egPhi/F");
   float closest_eg_dr; out_tree->Branch("closest_eg_dr", &closest_eg_dr, "closest_egEt/F");

   vector<float>* EM_L1_dPhi= new std::vector<float>; out_tree->Branch("EM_L1_dPhi", &EM_L1_dPhi);
   vector<float>* EM_L1_R= new std::vector<float>; out_tree->Branch("EM_L1_R", &EM_L1_R);
   vector<float>* EM_L2_dPhi= new std::vector<float>; out_tree->Branch("EM_L2_dPhi", &EM_L2_dPhi);
   vector<float>* EM_L2_R= new std::vector<float>; out_tree->Branch("EM_L2_R", &EM_L2_R);
   vector<float>* EM_L3_dPhi= new std::vector<float>; out_tree->Branch("EM_L3_dPhi", &EM_L3_dPhi);
   vector<float>* EM_L3_R= new std::vector<float>; out_tree->Branch("EM_L3_R", &EM_L3_R);
   vector<float>* EM_L4_dPhi= new std::vector<float>; out_tree->Branch("EM_L4_dPhi", &EM_L4_dPhi);
   vector<float>* EM_L4_R= new std::vector<float>; out_tree->Branch("EM_L4_R", &EM_L4_R);

   vector<float>* EM_pixel_dPhi= new std::vector<float>; out_tree->Branch("EM_pixel_dPhi", &EM_pixel_dPhi);

   vector<float>* EM_D1_dPhi= new std::vector<float>; out_tree->Branch("EM_D1_dPhi", &EM_D1_dPhi);
   vector<float>* EM_D2_dPhi= new std::vector<float>; out_tree->Branch("EM_D2_dPhi", &EM_D2_dPhi);
   vector<float>* EM_D3_dPhi= new std::vector<float>; out_tree->Branch("EM_D3_dPhi", &EM_D3_dPhi);
   vector<float>* EM_D4_dPhi= new std::vector<float>; out_tree->Branch("EM_D4_dPhi", &EM_D4_dPhi);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      EM_L1_dPhi->clear();
      EM_L2_dPhi->clear();
      EM_L3_dPhi->clear();
      EM_L4_dPhi->clear();

      EM_D1_dPhi->clear();
      EM_D2_dPhi->clear();
      EM_D3_dPhi->clear();
      EM_D4_dPhi->clear();

      EM_pixel_dPhi->clear();

      EM_L1_R->clear();
      EM_L2_R->clear();
      EM_L3_R->clear();
      EM_L4_R->clear();

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //find closest egamma object to the gen electron
      int egN_ = egCrysClusterEt->size();
      int cl3d_N_ = cl3d_pt->size();
      if(egN_==0 && cl3d_N_==0) continue; // skip if there is no eg objects

      float closest_dr = 9999.;
      int closest_eg = 0;
      int egCount = 0;

      for(int i=0; i < egN_;i++){

         float dPhi = propgenElPartPhi->at(0)-egCrysClusterPhi->at(i);
         if( dPhi > float(M_PI) ) dPhi -= float(2*M_PI);
         else if( dPhi <= -float(M_PI) ) dPhi += float(2*M_PI);

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
         if(egCrysClusterEt->at(i) < 5) continue;
         egCount++;
         if(current_dr < closest_dr){
           closest_dr = current_dr;
           closest_eg = i;
         }
      }// end of loop to find the closest egamma to gen electron 

      // HGCAL 3D cluster
      float closest_cl3d_dr = 9999.;
      int closest_cl3d = 0;
      int cl3d_Count = 0;

      for(int i=0; i < cl3d_N_;i++){

         float dPhi = propgenElPartPhi->at(0)-cl3d_phi->at(i);
         if( dPhi > float(M_PI) ) dPhi -= float(2*M_PI);
         else if( dPhi <= -float(M_PI) ) dPhi += float(2*M_PI);

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-cl3d_eta->at(i),2));
         if(cl3d_pt->at(i) < 5) continue;
         cl3d_Count++;
         if(current_dr < closest_cl3d_dr){
           closest_cl3d_dr = current_dr;
           closest_cl3d = i;
         }
      }// end of loop to find the closest egamma to gen electron 
      
      //if(egCount!=1) continue;
      if(egCount == 0 && cl3d_Count == 0) continue;
     
      TVector3 emvector;

      if( closest_dr < closest_cl3d_dr ){ 
         if(debug) cout << jentry << "th event " << " closest eg dR: " << closest_dr << " gen pt: " << genPartPt->at(0) << " eg et: " << egEt->at(closest_eg) << endl;
         abs_pterr = fabs(genPartPt->at(0)-egCrysClusterEt->at(closest_eg))/genPartPt->at(0);
         closest_egEt = egCrysClusterEt->at(closest_eg);
         closest_egEta = egCrysClusterEta->at(closest_eg);
         closest_egPhi = egCrysClusterPhi->at(closest_eg);
         closest_eg_dr = closest_dr;
         emvector.SetXYZ(egCrysClusterGx->at(closest_eg),egCrysClusterGy->at(closest_eg), egCrysClusterGz->at(closest_eg));
      }
      else{
         abs_pterr = fabs(genPartPt->at(0)-cl3d_pt->at(closest_cl3d))/genPartPt->at(0);
         closest_egEt =  cl3d_pt->at(closest_cl3d);
         closest_egEta = cl3d_eta->at(closest_cl3d);
         closest_egPhi = cl3d_phi->at(closest_cl3d);
         closest_eg_dr = closest_cl3d_dr;
         emvector.SetXYZ(cl3d_x->at(closest_cl3d),cl3d_y->at(closest_cl3d), cl3d_z->at(closest_cl3d));
      }

      std::vector<TVector3> first_layer_hits;
      std::vector<TVector3> second_layer_hits;
      std::vector<TVector3> third_layer_hits;
      std::vector<TVector3> fourth_layer_hits;

      std::vector<TVector3> first_disk_hits;
      std::vector<TVector3> second_disk_hits;
      std::vector<TVector3> third_disk_hits;
      std::vector<TVector3> fourth_disk_hits;
      std::vector<int> hitted_layers;

      int layers[9] = {};  // initialize as 0 for each event, save number of hits for each pixel layer 
      layers[0] = 1; // beam spot

      file->cd();

      int bpix_size = bRecHitGx->size();
      for(int a=0; a<bpix_size; a++){
         TVector3 current_hit;
         current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
      
         float temp_dphi = current_hit.Phi() - closest_egPhi;
         if( temp_dphi > float(M_PI) ) temp_dphi -= float(2*M_PI);
         else if( temp_dphi <= -float(M_PI) ) temp_dphi += float(2*M_PI);

         if( bRecHitLayer->at(a) == 1 ){
            //if(temp_dphi < L1_Dphi_cut1 && temp_dphi > L1_Dphi_cut2){
            layers[1]++;
            first_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            //}

         } 
         if( bRecHitLayer->at(a) == 2 ){
            //if(temp_dphi < L2_Dphi_cut1 && temp_dphi > L2_Dphi_cut2){
            layers[2]++;
            second_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            //}
         } 
         if( bRecHitLayer->at(a) == 3 ){
            //if(temp_dphi < L3_Dphi_cut1 && temp_dphi > L3_Dphi_cut2){
            layers[3]++;
            third_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            //}
         } 
         if( bRecHitLayer->at(a) == 4 ){
            //if(temp_dphi < L4_Dphi_cut1 && temp_dphi > L4_Dphi_cut2){
            layers[4]++;
            fourth_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            //}
         }
      }

      int fpix_size = fRecHitGx->size();
      for(int a=0; a<fpix_size; a++){
         TVector3 current_hit;
         current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );

         float temp_dphi = current_hit.Phi() - closest_egPhi;
         if( temp_dphi > float(M_PI) ) temp_dphi -= float(2*M_PI);
         else if( temp_dphi <= -float(M_PI) ) temp_dphi += float(2*M_PI);

         if( fRecHitDisk->at(a) == 1 ){
            //if(temp_dphi < L1_Dphi_cut1 && temp_dphi > L1_Dphi_cut2){
            layers[5]++;
            first_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            //}

         }
         if( fRecHitDisk->at(a) == 2 ){
            //if(temp_dphi < L2_Dphi_cut1 && temp_dphi > L2_Dphi_cut2){
            layers[6]++;
            second_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            //}
         }
         if( fRecHitDisk->at(a) == 3 ){
            //if(temp_dphi < L3_Dphi_cut1 && temp_dphi > L3_Dphi_cut2){
            layers[7]++;
            third_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            //}
         }
         if( fRecHitDisk->at(a) == 4 ){
            //if(temp_dphi < L4_Dphi_cut1 && temp_dphi > L4_Dphi_cut2){
            layers[8]++;
            fourth_disk_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            //}
         }
      }


     int three_outof_four = 0; // initialize as 0 for each event
     hitted_layers.push_back(0); // 0 mean beam spot i.e., (0,0,0)
     for( int i=1; i < 9; i++){
        if( layers[i] != 0 ){
          hitted_layers.push_back(i); //check if hits on each barrel or disk exists
          three_outof_four++;
        }
     }

     if( three_outof_four >= 1 ){ // save if there is at least one pixel hit
          for( std::vector<int>::iterator first_hit = hitted_layers.begin()+1; first_hit != hitted_layers.end(); first_hit++){ // hitted_layers.begin()+1: to avoid using beam spot
             for( int k=0; k < layers[*first_hit]; k++){
                  double dPhi = 0.;
                  double R = 0.;

                  TVector3 pixel_vector;

                  if( *first_hit == 1 ) pixel_vector = first_layer_hits[k];
                  if( *first_hit == 2 ) pixel_vector = second_layer_hits[k];
                  if( *first_hit == 3 ) pixel_vector = third_layer_hits[k];
                  if( *first_hit == 4 ) pixel_vector = fourth_layer_hits[k];
                  if( *first_hit == 5 ) pixel_vector = first_disk_hits[k];
                  if( *first_hit == 6 ) pixel_vector = second_disk_hits[k];
                  if( *first_hit == 7 ) pixel_vector = third_disk_hits[k];
                  if( *first_hit == 8 ) pixel_vector = fourth_disk_hits[k];

                  dPhi = deltaPhi(pixel_vector.Phi(), closest_egPhi);
                  R = emvector.Perp() - pixel_vector.Perp();

                  if( *first_hit == 1  ) {EM_L1_dPhi->push_back(dPhi); EM_L1_R->push_back(R);}
                  if( *first_hit == 2  ) {EM_L2_dPhi->push_back(dPhi); EM_L2_R->push_back(R);} 
                  if( *first_hit == 3  ) {EM_L3_dPhi->push_back(dPhi); EM_L3_R->push_back(R);}
                  if( *first_hit == 4  ) {EM_L4_dPhi->push_back(dPhi); EM_L4_R->push_back(R);}
                  if( *first_hit == 5  ) {EM_D1_dPhi->push_back(dPhi); }
                  if( *first_hit == 6  ) {EM_D2_dPhi->push_back(dPhi); } 
                  if( *first_hit == 7  ) {EM_D3_dPhi->push_back(dPhi); }
                  if( *first_hit == 8  ) {EM_D4_dPhi->push_back(dPhi); }

                  EM_pixel_dPhi->push_back(dPhi);
             }
           }
     }

      out_tree->Fill(); 
   }

   file->Write();
}
