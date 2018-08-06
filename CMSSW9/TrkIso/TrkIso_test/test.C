#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <TLorentzVector.h>

using namespace std;

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
  
   float dr_cut = 0.1;

   bit1 = 0x1;
   bit2 = 0x1;

   debug = false;

   //const double EM_PiX_dphi_width_[27] = {0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.50};
   //const double EM_PiX_deta_width_[27] = {0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027, 0.029, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.50};

   //const double PiX_PiX_dphi_width_[27] = {0.0005, 0.0007, 0.0009, 0.0011, 0.0013, 0.0015, 0.0017, 0.0019, 0.0021, 0.0023, 0.0025, 0.0027, 0.0029, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.02, 0.05};
   //const double PiX_PiX_deta_width_[27] = {0.0005, 0.0007, 0.0009, 0.0011, 0.0013, 0.0015, 0.0017, 0.0019, 0.0021, 0.0023, 0.0025, 0.0027, 0.0029, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.02, 0.05};


   const double EM_PiX_dphi_width_[9] = {0.025, 0.04, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05};
   const double EM_PiX_deta_width_[9] = {0.015, 0.03, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05};

   const double PiX_PiX_dphi_width_[9] = {0.0025, 0.0035, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005};
   const double PiX_PiX_deta_width_[9] = {0.0045, 0.0055, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005};

   // Histogram for number of tracks
   TH1F *h_track = new TH1F("h1"," ; Number of tracks; ",10,0,10);
   TH1F *h_pT = new TH1F("h2"," ; Isolation; ",100,0.,1.);

  //nentries = 100;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;
      FillCutFlow("NoCut", 1.);

      matchedEgEt  = -999.;
      matchedEgEta = -999.;
      matchedEgPhi = -999.;
      fired        = 0;;

      ntnEg2 = 0;
      ntEgEt.clear();
      ntEgEta.clear();
      ntEgPhi.clear();

      PiXTRKbit.clear();
      trigger_bit_width.clear();
      pix_comb.clear();

      ntCl_match.clear();
      isTrack_match.clear();
      chi2.clear();
      track_dr.clear();
      withoutEM_match.clear();
      withEM_match.clear();

      ntfirstPix.clear();
      ntsecondPix.clear();
      ntthirdPix.clear();
      ntfourthPix.clear();

      nt_genPhi = propgenElPartPhi->at(0);
      nt_genEta = propgenElPartEta->at(0);
      nt_genPt = propgenElPartPt->at(0);

      nt_lastSimtkpt = lastSimtkpt;
      nt_initialSimtkpt = initialSimtkpt;
       
      // Track isolation variables initialize
      /*
      ntnEg3 = 0;

      iso_ntEgEt.clear(); 
      iso_ntEgEta.clear(); 
      iso_ntEgPhi.clear(); 
   
      PiXTRKbit_iso.clear();
      pix_comb_iso.clear();
      trigger_bit_width_iso.clear();

      iso_ntCl_match.clear();
      iso_only_match.clear();
      */
 
      float tempDR = 999.;
      int   indx = -1;
      int   indx_closestEg = -1;

      //find closest egamma object to the gen electron
      EgN=egCrysClusterEt->size();
      int cl3d_N_ = cl3d_pt->size();

      float closest_dr = 9999.;
      int closest_eg = 0;
      int egCount = 0;

      // loop over barrel egamma objects
      for(int i=0; i < EgN;i++){

         float dPhi = deltaPhi(propgenElPartPhi->at(0), egCrysClusterPhi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
         if(egCrysClusterEt->at(i) < 10) continue;
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

         // simplge Egamma ID for HGCAL
         //if(cl3d_coreshowerlength->at(i) < 3 || cl3d_coreshowerlength->at(i) > 18) continue;
         //if(cl3d_srrtot->at(i) < 0.002 || cl3d_srrtot->at(i) > 0.005) continue;
         //if(cl3d_maxlayer->at(i) < 8 || cl3d_maxlayer->at(i) > 20) continue;
         //if(cl3d_firstlayer->at(i) > 5) continue;

         if(cl3d_egid->at(i) != 1) continue;

         float dPhi = deltaPhi(propgenElPartPhi->at(0), cl3d_phi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-cl3d_eta->at(i),2));
         if(cl3d_pt->at(i) < 10) continue;
         cl3d_Count++;
         if(current_dr < closest_cl3d_dr){
           closest_cl3d_dr = current_dr;
           closest_cl3d = i;
         }
      }// end of loop to find the closest egamma to gen electron 

      pix_comb_ = 0x0;

      nPix123_segments = 0;
      nPix124_segments = 0;
      nPix134_segments = 0;
      nPix234_segments = 0;

     // find egamma objects passing pixtrk signal windows
     if((closest_cl3d_dr < dr_cut && closest_cl3d_dr != 9999.)|| (closest_dr < dr_cut && closest_dr != 9999.)){

      debug = false;

      indx++; //remove this variable

      if( closest_dr < closest_cl3d_dr ){
        EgEt =egCrysClusterEt ->at(closest_eg);
        EgEta=egCrysClusterEta->at(closest_eg);
        EgPhi=egCrysClusterPhi->at(closest_eg);

        isTrack_match.push_back(isTrackMatched->at(closest_eg));
        chi2.push_back(trackHighestPtCutChi2Chi2->at(closest_eg));
        track_dr.push_back(trackmatchingdR->at(closest_eg));

        float EgGx = egCrysClusterGx->at(closest_eg);
        float EgGy = egCrysClusterGy->at(closest_eg);
        float EgGz = egCrysClusterGz->at(closest_eg);
        emvector.SetXYZ(EgGx,EgGy,EgGz);

        tempDR = closest_dr; 
        matchedEgEta = EgEta;
        matchedEgPhi = EgPhi;
        matchedEgEt  = EgEt;
        indx_closestEg = indx;
      }
      else{

          EgEt =cl3d_pt->at(closest_cl3d);
          EgEta=cl3d_eta->at(closest_cl3d);
          EgPhi=cl3d_phi->at(closest_cl3d);

          isTrack_match.push_back(hgcal_isTrackMatched->at(closest_cl3d));
          chi2.push_back(hgcal_trackHighestPtCutChi2Chi2->at(closest_cl3d));
          track_dr.push_back(hgcal_trackmatchingdR->at(closest_cl3d));

          float EgGx = cl3d_x->at(closest_cl3d);
          float EgGy = cl3d_y->at(closest_cl3d);
          float EgGz = (float)cl3d_z->at(closest_cl3d);
          emvector.SetXYZ(EgGx,EgGy,EgGz);

          tempDR = closest_cl3d_dr; 
          matchedEgEta = EgEta;
          matchedEgPhi = EgPhi;
          matchedEgEt  = EgEt;
          indx_closestEg = indx;
      }


      if( fabs(EgEta) <= 0.8 ) eta_region =1;
      //if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
      //if( fabs(EgEta) <= 1.8 && fabs(EgEta) > 1.4 ) eta_region =3;
      //if( fabs(EgEta) <= 2.7 && fabs(EgEta) > 1.8 ) eta_region =4;
      //if( fabs(EgEta) <= 2.9 && fabs(EgEta) > 2.7 ) eta_region =5;
      //if( fabs(EgEta) <= 3.0 && fabs(EgEta) > 2.9 ) eta_region =6;

      //if( fabs(EgEta) > 3. ) continue;
      if( fabs(EgEta) > 0.8 ) continue;

      Bool_t flag123 = false;
      Bool_t flag124 = false;
      Bool_t flag134 = false;
      Bool_t flag234 = false;

      Float_t recoPV = 0.;

      Float_t zp14 = -99.;
      Float_t zp13 = -99.;
      Float_t zp24 = -99.;

      ntnEg2++;
      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);
      ntEgPhi.push_back(EgPhi);
      
      // set regin of interest
      SetROI(eta_region);

      // initialize pixel hit variables
      first_layer_hits.clear();
      second_layer_hits.clear();
      third_layer_hits.clear();
      fourth_layer_hits.clear();

      first_layer_hits_Ele_or_Pos.clear();
      second_layer_hits_Ele_or_Pos.clear();
      third_layer_hits_Ele_or_Pos.clear();
      fourth_layer_hits_Ele_or_Pos.clear();
      hitted_layers.clear();

      
      layers[0] = 1; // beam spot
      layers[1] = 0; layers[2] = 0; layers[3] = 0; layers[4] = 0;
      r = 0;

      StorePixelHit(eta_region); // save pixel hits in Region of Interest for the given eta region

      // check which pixel has hits
       for( int i=1; i < 5; i++){ 
          if( layers[i] != 0 ){ 
            hitted_layers.push_back(i); 
          }
       }

       int global_index_width = 0;
       trigger_bit_width_ = 0x0;
       // set pixtrk signal boundary
       for(int nth_eg_pix_deta = 0; nth_eg_pix_deta < 9; nth_eg_pix_deta++){
       if( nth_eg_pix_deta != 0) continue;
       
       if(eta_region == 1) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);
       else if(eta_region == 2) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
       else if(eta_region == 6) SetSingalBoundary(5, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
       else SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);

       //SetSingalBoundary(1);

       // PixTRK algorithm 
       PixTrkPassed = false;
       withoutEM_count_Ele = 0, withEM_count_Ele = 0;

       fourth_layer_missing = 0;
       third_layer_missing = 0;
       second_layer_missing = 0;
       first_layer_missing = 0;

       // loop over every 3 out of 4 pixel combination 
       for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
          for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
              for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){

                 
                 // loop over every pixel hits in the given pixel combination
                 for( int k=0; k < layers[*first_hit]; k++){
                    for( int i=0; i < layers[*second_hit]; i++){
                        _pass_Ele = 0, _pass_Pos = 0;
                        L012_pass_Ele = 0, L012_pass_Pos = 0;
                        L013_pass_Ele = 0, L013_pass_Pos = 0;
                        L023_pass_Ele = 0, L023_pass_Pos = 0;

                        if( *first_hit == 1 && *second_hit == 2 )
                          TriggeringWith_1st2ndPixel(k,i);

                        if( *first_hit == 1 && *second_hit == 3 )
                          TriggeringWith_1st3rdPixel(k,i);

                        if( *first_hit == 2 && *second_hit == 3 )
                          TriggeringWith_2nd3rdPixel(k,i);

                        // skip only if both _pass_Ele and _pass_Pos are 0 i.e., both electron and positron signal window are not satisfied
                        if( !_pass_Ele && !_pass_Pos ) continue;

                        for( int j=0; j < layers[*third_hit]; j++){
                            all_cut_pass_Ele = 0, all_cut_pass_Pos = 0;
                            withoutEM_pass_Ele = 0, withEM_pass_Ele = 0;

                            L012_pass_Ele = 0, L012_pass_Pos = 0;
                            L013_pass_Ele = 0, L013_pass_Pos = 0;
                            L014_pass_Ele = 0, L014_pass_Pos = 0;
                            L023_pass_Ele = 0, L023_pass_Pos = 0;
                            L024_pass_Ele = 0, L024_pass_Pos = 0;
                            L034_pass_Ele = 0, L034_pass_Pos = 0;
                            L123_pass_Ele = 0, L123_pass_Pos = 0;
                            L124_pass_Ele = 0, L124_pass_Pos = 0;
                            L134_pass_Ele = 0, L134_pass_Pos = 0;
                            L234_pass_Ele = 0, L234_pass_Pos = 0;

                            L12_EM_Ele = 0, L12_EM_Pos = 0;
                            L13_EM_Ele = 0, L13_EM_Pos = 0;
                            L14_EM_Ele = 0, L14_EM_Pos = 0;
                            L23_EM_Ele = 0, L23_EM_Pos = 0;
                            L24_EM_Ele = 0, L24_EM_Pos = 0;
                            L34_EM_Ele = 0, L34_EM_Pos = 0;

            	            dPhi = StandaloneDPhi( *first_hit, *second_hit, *third_hit, k, i, j );
                            dEta = StandaloneDEta( *first_hit, *second_hit, *third_hit, k, i, j );

                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second hit is in the second pixel layer and the third hit is in the third layer. 
                                TriggeringWithout_4thPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                         all_cut_pass_Ele = 1; 
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele)
                                         withoutEM_pass_Ele = 1;
                                      if(L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                         withEM_pass_Ele = 1;
                                }
 
                                if( L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos &&
            			    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 1);
                                   nPix123_segments++;
                                }
                                // Save vertex coordinate and check whether pass or not
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                    flag123 = true;
                                    Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    Float_t R3 = sqrt(pow(third_layer_hits[j].X(),2)+pow(third_layer_hits[j].Y(),2));
                                    Float_t Z1 = first_layer_hits[k].Z();
                                    Float_t Z3 = third_layer_hits[j].Z();
                                    zp13 = (R3*Z1 - R1*Z3)/(R3-R1);
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit loop
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  fourth_layer_missing = 1;
                   		 }
                               }
                              }
                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second is in the second layer and the third hit is in the fourth layer.
                                TriggeringWithout_3rdPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele && L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos &&
            			    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 2);
                                   nPix124_segments++;
                                }

                                // Save vertex coordinate and check whether pass or not
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                    flag124 = true;
                                    Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    Float_t Z1 = first_layer_hits[k].Z();
                                    Float_t Z4 = fourth_layer_hits[j].Z();
                                    zp14 = (R4*Z1 - R1*Z4)/(R4-R1);
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
            	     	          third_layer_missing = 1;
                   		 }
                               }
                              }
                              if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                TriggeringWithout_2ndPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele && L13_EM_Ele && L14_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                }

                                if( L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos &&
            			    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 3);
                                   nPix134_segments++;
                                }
                                // Save vertex coordinate and check whether pass or not
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                    flag134 = true;
                                    Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    Float_t Z1 = first_layer_hits[k].Z();
                                    Float_t Z4 = fourth_layer_hits[j].Z();
                                    zp14 = (R4*Z1 - R1*Z4)/(R4-R1);
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  second_layer_missing = 1; 
                   		 }
                               }
                              }
                              if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                TriggeringWithout_1stPixel(k, i, j);

                                if( (second_layer_hits_Ele_or_Pos[k] == 1 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele && L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos &&
            			    (second_layer_hits_Ele_or_Pos[k] == 2 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 
                   
                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 4);
                                   nPix234_segments++;
                                } 
                                // Save vertex coordinate and check whether pass or not
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ) {
                                    flag234 = true;
                                    Float_t R2 = sqrt(pow(second_layer_hits[k].X(),2)+pow(second_layer_hits[k].Y(),2));
                                    Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    Float_t Z2 = second_layer_hits[k].Z();
                                    Float_t Z4 = fourth_layer_hits[j].Z();
                                    zp24 = (R4*Z2 - R2*Z4)/(R4-R2);
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  first_layer_missing = 1;
                   		 } 
                               }
                              }

                         if( all_cut_pass_Ele == 1 ) { PixTrkPassed = true;}
                         if( withoutEM_pass_Ele == 1 ) withoutEM_count_Ele = 1;
                         if( withEM_pass_Ele == 1 ) withEM_count_Ele = 1;
                       } // loop for third layer hits
                   } // loop for second layer hits       
                 } // loop for first layer hits

          }          
        }
      }

     if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 1.3 && PixTrkPassed == false) {
       for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
          for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
                 
                 // loop over every pixel hits in the given pixel combination
                 for( int k=0; k < layers[*first_hit]; k++){
                    for( int i=0; i < layers[*second_hit]; i++){
                        _pass_Ele = 0, _pass_Pos = 0;

                        L012_pass_Ele = 0, L012_pass_Pos = 0;
                        L013_pass_Ele = 0, L013_pass_Pos = 0;
                        L014_pass_Ele = 0, L014_pass_Pos = 0;
                        L023_pass_Ele = 0, L023_pass_Pos = 0;
                        L024_pass_Ele = 0, L024_pass_Pos = 0;
                        L034_pass_Ele = 0, L034_pass_Pos = 0;
                        L123_pass_Ele = 0, L123_pass_Pos = 0;
                        L124_pass_Ele = 0, L124_pass_Pos = 0;
                        L134_pass_Ele = 0, L134_pass_Pos = 0;
                        L234_pass_Ele = 0, L234_pass_Pos = 0;

                        L12_EM_Ele = 0, L12_EM_Pos = 0;
                        L13_EM_Ele = 0, L13_EM_Pos = 0;
                        L14_EM_Ele = 0, L14_EM_Pos = 0;
                        L23_EM_Ele = 0, L23_EM_Pos = 0;
                        L24_EM_Ele = 0, L24_EM_Pos = 0;
                        L34_EM_Ele = 0, L34_EM_Pos = 0;

                        // skip only if both _pass_Ele and _pass_Pos are 0 i.e., both electron and positron signal window are not satisfied

                        if( *first_hit == 1 && *second_hit == 2 ){ // for efficiency counting  !!caution of the position of this codition

                          TriggeringWith_1st2ndPixel_v2(k,i); 

                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
            		    k = layers[*first_hit];
                            i = layers[*second_hit];
                   	   }
                         }
                        }

                        if( *first_hit == 1 && *second_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_1st3rdPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                        if( *first_hit == 1 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_1st4thPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                        if( *first_hit == 2 && *second_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_2nd3rdPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                        if( *first_hit == 2 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition

                          TriggeringWith_2nd4thPixel_v2(k,i);

                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         }
                        }

                        if( *first_hit == 3 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_3rd4thPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                       if( _pass_Ele == 1 ) { PixTrkPassed = true;}
                   } // loop for second layer hits       
                 } // loop for first layer hits

        }
      }
     }
  //    cout << "global_index_width: " << global_index_width << " PixTrkPassed: " << PixTrkPassed << " (bit1 << global_index_width) " << (bit1 << global_index_width) <<  endl;
      if( PixTrkPassed ){ 
          //trigger_bit_width_ = trigger_bit_width_| (bit1 << global_index_width);
          trigger_bit_width_ = trigger_bit_width_| (bit1 << nth_eg_pix_deta);
      }

      /////////////////////////////////

     global_index_width++;
     }

     //PiXTRKbit.push_back(PiXTRKbit_);
     //trigger_bit_width.push_back(trigger_bit_width_);

     // Start track isolation algorithm
     
       
     // Limited to barrel region
     
     if( flag123 == false && flag124 == false && flag134 == false && flag234 == false ) 
     {
         //iso_ntCl_match.push_back(false);
         continue;  
     }
     if( EgEta > 0.8 ) continue; 
     
     if( flag124 || flag134 ) recoPV = zp14;
     if( !flag124 && !flag134 && flag123 ) recoPV = zp13;
     if( !flag124 && !flag134 && !flag123 && flag234 ) recoPV = zp24;

     if( recoPV == -99. ) continue;
     L1.clear(); 
     L2.clear();
     L3.clear(); 
     L4.clear();
    
     for(Int_t i = 0; i < bRecHitN; i++)
     {
         Float_t R = sqrt(pow(bRecHitGx->at(i), 2)+ pow(bRecHitGy->at(i), 2));
         TVector3 pixel;         
         pixel.SetXYZ(bRecHitGx->at(i), bRecHitGy->at(i), bRecHitGz->at(i) - recoPV);
         Float_t pixelPhi = pixel.Phi();
         Float_t pixelEta = pixel.Eta();
         Float_t deltaPhi = EgPhi - pixelPhi;
         if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
         if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
         Float_t deltaEta = EgEta - pixelEta;
         Float_t deltaR = sqrt(pow(deltaPhi, 2) + pow(deltaEta, 2));
         if( fabs(deltaPhi) > 0.3 ) continue;

         if( R < 5. && deltaR < 0.3 )             L1.push_back(TVector3(bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
         if( R > 5. && R < 9. && deltaR < 0.3 )   L2.push_back(TVector3(bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
         if( R > 9. && R < 14. && deltaR < 0.3 )  L3.push_back(TVector3(bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
         if( R > 14. && R < 18. && deltaR < 0.3 ) L4.push_back(TVector3(bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));

     }
     
     L123.clear();
     L124.clear();
     L134.clear();
     L234.clear();

     // Call TrkIso_Combi function to store pixel clusters combinations
     TrkIso_L123(recoPV);
     TrkIso_L124(recoPV);
     TrkIso_L134(recoPV);
     TrkIso_L234(recoPV);

     // Erase duplication in each combination
     if( L123.size() >= 2 ) 
     {
         sort(L123.begin(), L123.end(), track::comp3);
         L123.erase(unique(L123.begin(), L123.end(), track::uni3),L123.end());
         sort(L123.begin(), L123.end(), track::comp2);
         L123.erase(unique(L123.begin(), L123.end(), track::uni2),L123.end());
         sort(L123.begin(), L123.end(), track::comp1);
         L123.erase(unique(L123.begin(), L123.end(), track::uni1),L123.end());
     }
     if( L124.size() >= 2 ) 
     {
         sort(L124.begin(), L124.end(), track::comp3);
         L124.erase(unique(L124.begin(), L124.end(), track::uni3),L124.end());
         sort(L124.begin(), L124.end(), track::comp2);
         L124.erase(unique(L124.begin(), L124.end(), track::uni2),L124.end());
         sort(L124.begin(), L124.end(), track::comp1);
         L124.erase(unique(L124.begin(), L124.end(), track::uni1),L124.end());
     }
     if( L134.size() >= 2 ) 
     {
         sort(L134.begin(), L134.end(), track::comp3);
         L134.erase(unique(L134.begin(), L134.end(), track::uni3),L134.end());
         sort(L134.begin(), L134.end(), track::comp2);
         L134.erase(unique(L134.begin(), L134.end(), track::uni2),L134.end());
         sort(L134.begin(), L134.end(), track::comp1);
         L134.erase(unique(L134.begin(), L134.end(), track::uni1),L134.end());
     }
     if( L234.size() >= 2 ) 
     {
         sort(L234.begin(), L234.end(), track::comp3);
         L234.erase(unique(L234.begin(), L234.end(), track::uni3),L234.end());
         sort(L234.begin(), L234.end(), track::comp2);
         L234.erase(unique(L234.begin(), L234.end(), track::uni2),L234.end());
         sort(L234.begin(), L234.end(), track::comp1);
         L234.erase(unique(L234.begin(), L234.end(), track::uni1),L234.end());
     }

     // Make vector to contain all pixel clusters combinations from different layer combinations and erase duplication
     vector<track> all;
     all.clear();

     for(vector<track>::iterator a1 = L123.begin(); a1 != L123.end(); ++a1) 
         all.push_back(track((*a1).pos_x3, (*a1).pos_x2, (*a1).pos_x1, (*a1).pos_y3, (*a1).pos_y2, (*a1).pos_y1, (*a1).pos_z3, (*a1).pos_z2, (*a1).pos_z1, 1 ));
     for(vector<track>::iterator a2 = L124.begin(); a2 != L124.end(); ++a2) 
         all.push_back(track((*a2).pos_x3, (*a2).pos_x2, (*a2).pos_x1, (*a2).pos_y3, (*a2).pos_y2, (*a2).pos_y1, (*a2).pos_z3, (*a2).pos_z2, (*a2).pos_z1, 2 ));
     for(vector<track>::iterator a3 = L134.begin(); a3 != L134.end(); ++a3) 
         all.push_back(track((*a3).pos_x3, (*a3).pos_x2, (*a3).pos_x1, (*a3).pos_y3, (*a3).pos_y2, (*a3).pos_y1, (*a3).pos_z3, (*a3).pos_z2, (*a3).pos_z1, 3 ));
     for(vector<track>::iterator a4 = L234.begin(); a4 != L234.end(); ++a4) 
         all.push_back(track((*a4).pos_x3, (*a4).pos_x2, (*a4).pos_x1, (*a4).pos_y3, (*a4).pos_y2, (*a4).pos_y1, (*a4).pos_z3, (*a4).pos_z2, (*a4).pos_z1, 4 ));

     sort(all.begin(), all.end(), track::comp3);
     all.erase(unique(all.begin(), all.end(), track::uni3),all.end());
     sort(all.begin(), all.end(), track::comp2);
     all.erase(unique(all.begin(), all.end(), track::uni2),all.end());
     sort(all.begin(), all.end(), track::comp1);
     all.erase(unique(all.begin(), all.end(), track::uni1),all.end());

     cout << "Number of tracks: " << all.size() << endl;
     cout << "=============================" << endl << endl;
     if( EgEt > 20. ) 
     {
         h_track->Fill(all.size()); 
         vector<Float_t> pT_vector; 
         if( all.size() <= 1 ) h_pT->Fill(0.01);
         if( all.size() >= 2 )
         {
             //for(vector<track>::iterator cur = all.begin(); cur != all.end(); ++cur)
             for(Int_t cur = 0; cur < all.size(); cur++)
             {
                 pT_vector.clear();
                 if( all[cur].index == 1 )
                 {
                     Float_t p[5] = {};
                     p[0] = -0.0865347;
                     p[1] = 0.0055108;
                     p[2] = -1.05129;
                     p[3] = -0.0268743;
                     p[4] = 3.3316;
                     
                     TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x1, all[cur].pos_y1, all[cur].pos_z1 - recoPV );
                     TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x1, all[cur].pos_y3 - all[cur].pos_y1, all[cur].pos_z3 - all[cur].pos_z1 );
                    
                     Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                     Float_t deltaPhi = phi1 - phi2;
                     if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
                     if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
                     Float_t x = fabs(deltaPhi);
                     Float_t recopT = p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);
                     pT_vector.push_back(recopT);
                 }
                 if( all[cur].index == 2 )
                 {
                     Float_t p[5] = {};
                     p[0] = -0.228916;
                     p[1] = 0.0059613;
                     p[2] = -1.1059;
                     p[3] = -0.0653335;
                     p[4] = 3.52806;

                     TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x2, all[cur].pos_y2, all[cur].pos_z2 - recoPV );
                     TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x2, all[cur].pos_y3 - all[cur].pos_y2, all[cur].pos_z3 - all[cur].pos_z2 );
                    
                     Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                     Float_t deltaPhi = phi1 - phi2;
                     if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
                     if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
                     Float_t x = fabs(deltaPhi);
                     Float_t recopT = p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);
                     pT_vector.push_back(recopT);
                 }
                 if( all[cur].index == 3 )
                 {
                     Float_t p[5] = {};
                     p[0] = 0.194839; 
                     p[1] = 0.00561084; 
                     p[2] = -1.24477; 
                     p[3] = -0.113004; 
                     p[4] = 3.23352; 

                     TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x2, all[cur].pos_y2, all[cur].pos_z2 - recoPV );
                     TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x2, all[cur].pos_y3 - all[cur].pos_y2, all[cur].pos_z3 - all[cur].pos_z2 );
                    
                     Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                     Float_t deltaPhi = phi1 - phi2;
                     if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
                     if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
                     Float_t x = fabs(deltaPhi);
                     Float_t recopT = p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);
                     pT_vector.push_back(recopT);
                 }
                 if( all[cur].index == 4 )
                 {
                     Float_t p[5] = {};
                     p[0] = -0.228916;
                     p[1] = 0.0059613;
                     p[2] = -1.1059;
                     p[3] = -0.0653335;
                     p[4] = 3.52806;

                     TVector3 pixel1; pixel1.SetXYZ( all[cur].pos_x1, all[cur].pos_y1, all[cur].pos_z1 - recoPV );
                     TVector3 pixel2; pixel2.SetXYZ( all[cur].pos_x3 - all[cur].pos_x2, all[cur].pos_y3 - all[cur].pos_y2, all[cur].pos_z3 - all[cur].pos_z2 );
                    
                     Float_t phi1 = pixel1.Phi(); Float_t phi2 = pixel2.Phi();
                     Float_t deltaPhi = phi1 - phi2;
                     if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
                     if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
                     Float_t x = fabs(deltaPhi);
                     Float_t recopT = p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);
                     pT_vector.push_back(recopT);
                 }
             }

             sort(pT_vector.begin(), pT_vector.end());
             Float_t denomi = 0.; Float_t nomi = 0.;
             for(Int_t k = 0; k < pT_vector.size(); k++) denomi += pT_vector.at(k);
             for(Int_t k = 0; k < pT_vector.size()-1; k++) nomi += pT_vector.at(k);
             h_pT->Fill(nomi/denomi);
         }
     }
     
  } // end of egamma loop    

  trigger_bit_width.push_back(trigger_bit_width_);
  pix_comb.push_back(pix_comb_);

  pixtrk_tree->Fill();
 } // end of entries loop 
   outfile->Write();
}

