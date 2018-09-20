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

   int bit1 = 0x1;

   //const double EM_PiX_dphi_width_[9] = {0.025, 0.04, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; // as of current 3 Sep
   //const double EM_PiX_deta_width_[9] = {0.015, 0.03, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; // as of current 3 Sep

   const double EM_PiX_dphi_width_[9] = {0.02, 0.04, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; // 
   const double EM_PiX_deta_width_[9] = {0.01, 0.03, 0.031, 0.033, 0.035, 0.037, 0.039, 0.04, 0.05}; // 

   //const double PiX_PiX_dphi_width_[9] = {0.0025, 0.0035, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005}; // as of current 3 Sep
   //const double PiX_PiX_deta_width_[9] = {0.0045, 0.0055, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005}; // as of current 3 Sep

   const double PiX_PiX_dphi_width_[9] = {0.002, 0.0035, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005};
   const double PiX_PiX_deta_width_[9] = {0.004, 0.0055, 0.0031, 0.0033, 0.0035, 0.0037, 0.0039, 0.004, 0.005};
   //nentries = 2000;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
   //for (Long64_t jentry=0; jentry<10;jentry++) { //nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;
      FillCutFlow("NoCut", 1.);
        
      //EgN=egN;
      EgN=egCrysEt->size();;

      pass_egobjects_check = 0;
      ntnEg2 = 0;
      ntEgEt.clear();
      ntEgEta.clear();
      ntclusterIsEG.clear();
      ntCl_match.clear();
      isTrack_match.clear();
      chi2.clear();
      track_dr.clear();
      withoutEM_match.clear();
      withEM_match.clear();
      pass_Ele.clear();
      pass_Pos.clear();
      pass_ElePos.clear();

      trigger_bit_width.clear();
      trigger_4Hits_bit_width.clear();

      ntL1TkEgEt.clear();
      ntL1TkEgEta.clear();
      ntL1TkEgPhi.clear();

      ntL1TkLooseEgEt.clear();
      ntL1TkLooseEgEta.clear();
      ntL1TkLooseEgPhi.clear();

      ntCl_match_wo4thPix.clear();
      ntCl_match_wo3thPix.clear();
      ntCl_match_wo2thPix.clear();
      ntCl_match_wo1thPix.clear();

      Npass_woEM_wo4thPix.clear();
      Npass_woEM_wo3thPix.clear();
      Npass_woEM_wo2thPix.clear();
      Npass_woEM_wo1thPix.clear();

      Npass_wEM_wo4thPix.clear();
      Npass_wEM_wo3thPix.clear();
      Npass_wEM_wo2thPix.clear();
      Npass_wEM_wo1thPix.clear();

      ntfirstPix.clear();
      ntsecondPix.clear();
      ntthirdPix.clear();
      ntfourthPix.clear();

      ntL1TkEleEt.clear();
      ntL1TkEleEta.clear();
      ntL1TkElePhi.clear();

      ntL1TkEleIsoEt.clear();
      ntL1TkEleIsoEta.clear();
      ntL1TkEleIsoPhi.clear();
    
      dphi_L12EM_wo4thPix.clear();
      dphi_L12EM_wo3thPix.clear();
      dphi_L12EM_wo2thPix.clear();
      dphi_L12EM_wo1thPix.clear();

      deta_L12EM_wo4thPix.clear();
      deta_L12EM_wo3thPix.clear();
      deta_L12EM_wo2thPix.clear();
      deta_L12EM_wo1thPix.clear();

      all_cut_pass_eg = 0;
      event_nominator = 0;
      event_denominator = 0;
   
      // cout flow
      float EgEtCut  = 0;
      for(int k=0; k<EgN; k++) {
        EgEt =egCrysEt ->at(k);
        if(EgEt > 8 ) EgEtCut = 1;
      }
      //if(EgEtCut == 0) continue;
      if(EgEtCut == 1) FillCutFlow("MinEtCut", 1.);

      //save L1TkElectron
      for(int i = 0; i < tkEGEt->size(); i++){
         ntL1TkEgEt.push_back(tkEGEt->at(i));
         ntL1TkEgEta.push_back(tkEGEta->at(i));
         ntL1TkEgPhi.push_back(tkEGPhi->at(i));
      }

      //save L1TkElectron
      for(int i = 0; i < tkEGLooseEt->size(); i++){
         ntL1TkLooseEgEt.push_back(tkEGLooseEt->at(i));
         ntL1TkLooseEgEta.push_back(tkEGLooseEta->at(i));
         ntL1TkLooseEgPhi.push_back(tkEGLoosePhi->at(i));
      }

      int EtaCutFlow = 0;
      for(int k=0; k<EgN; k++) {
        EgEta=egCrysEta->at(k);
        EgEt =egCrysEt ->at(k);
        //if(fabs(EgEta) < 1.3 && EgEt > 10) EtaCutFlow = 1; // for first η region
        //if(fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 && EgEt > 10) EtaCutFlow = 1;
        if(fabs(EgEta) < 1.7 && EgEt > 8) EtaCutFlow = 1;
      }
      //if(EtaCutFlow == 0) continue;
      if(EtaCutFlow == 1) FillCutFlow("EtaCut", 1.);


     // find egamma objects passing pixtrk signal windows
     for( int q=0; q<EgN; q++){ 
      EgEt =egCrysEt ->at(q);
      EgEta=egCrysEta->at(q);
      EgPhi=egCrysPhi->at(q);

      float EgGx = egCrysGx->at(q);
      float EgGy = egCrysGy->at(q);
      float EgGz = egCrysGz->at(q);
      emvector.SetXYZ(EgGx,EgGy,EgGz);

      if(EgEt < 8) continue;

      eta_region = 0; // initialize variable 
      if( fabs(EgEta) <= 0.8 ) eta_region =1;
      if( fabs(EgEta) <= 1.15 && fabs(EgEta) > 0.8 ) eta_region =2;
      if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 1.15 ) eta_region =3;
      if( fabs(EgEta) <= 1.7 && fabs(EgEta) > 1.4 ) eta_region =4;
      if( fabs(EgEta) <= 2.25 && fabs(EgEta) > 1.7 ) eta_region =5;
      if( fabs(EgEta) <= 2.7 && fabs(EgEta) > 2.25 ) eta_region =6;
      if( fabs(EgEta) <= 3.0 && fabs(EgEta) > 2.7 ) eta_region =7;

      if( fabs(EgEta) > 3. ) continue;
      if( eta_region < 2 || eta_region > 3 ) continue;

      pass_egobjects_check = 1;
      ntnEg2++;
      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);

      isTrack_match.push_back(isTrackMatched->at(q));
      chi2.push_back(trackHighestPtCutChi2Chi2->at(q));
      track_dr.push_back(trackmatchingdR->at(q));
      
      // set regin of interest
      // for η < 1.3, Δφ < 0.05 
      //if(eta_region==1)SetROI(2); 
      //else SetROI(eta_region);
      //SetROI(1);

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

      bool isFourHits = true;
      // check which pixel has hits
       for( int i=1; i < 5; i++){
          if( layers[i] != 0 ){
            hitted_layers.push_back(i);
          }
          else{
              isFourHits = false;
          }
       }
   
       // set pixtrk signal boundary
       //if(eta_region==1)SetSingalBoundary(2);
       //else SetSingalBoundary(eta_region);
       //SetSingalBoundary(1);

       int global_index_width = 0;
       trigger_bit_width_ = 0x0;
       PiXTRKbit_4Hits_ = 0x0; 

       for(int nth_eg_pix_deta = 0; nth_eg_pix_deta < 9; nth_eg_pix_deta++){
       if(nth_eg_pix_deta != 0) continue;

       //SetSingalBoundary(eta_region);
       if(eta_region == 1) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);
       else if(eta_region == 2) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
       else if(eta_region == 6) SetSingalBoundary(5, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
       else SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);

       // PixTRK algorithm 
       pass_count = 0;
       pass_count_wo4thPix = 0, pass_count_wo3thPix = 0, pass_count_wo2thPix = 0, pass_count_wo1thPix = 0;
       woEM_pass_count_wo4thPix = 0, woEM_pass_count_wo3thPix = 0, woEM_pass_count_wo2thPix = 0, woEM_pass_count_wo1thPix = 0;
       wEM_pass_count_wo4thPix = 0, wEM_pass_count_wo3thPix = 0, wEM_pass_count_wo2thPix = 0, wEM_pass_count_wo1thPix = 0;
       withoutEM_count_Ele = 0, withEM_count_Ele = 0;

       fourth_layer_missing = 0;
       third_layer_missing = 0;
       second_layer_missing = 0;
       first_layer_missing = 0;

      bool PixTrkFourHits = false;

/*
      // require 4 pixel hits
      if(isFourHits){

         for( int k=0; k < layers[1]; k++){
            for( int i=0; i < layers[2]; i++){

               _pass_Ele = 0, _pass_Pos = 0;
               TriggeringWith_1st2ndPixel(k,i);
               if( !_pass_Ele && !_pass_Pos ) continue;

               for( int j=0; j < layers[3]; j++){
                  bool wo4thHits = 0;

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

                  dPhi = StandaloneDPhi( 1, 2, 3, k, i, j );
                  dEta = StandaloneDEta( 1, 2, 3, k, i, j );

                  TriggeringWithout_4thPixel(k, i, j);

                  if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                       (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                       (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3)){
                         if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                            wo4thHits = true;
                  }
                  if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                      (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                      (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3)){
                      if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) wo4thHits = true;
                  }

                  if( wo4thHits == 0) continue;


                  for( int h=0; h < layers[4]; h++){

                     bool check124 = false;
                     bool check134 = false;
                     bool check234 = false;
                     TriggeringWithout_3rdPixel(k, i, h); // 124

                     if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                         (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                         (fourth_layer_hits_Ele_or_Pos[h] == 1 || fourth_layer_hits_Ele_or_Pos[h] ==3)){
                         if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele && L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) check124 = true;
                     }
                     if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                         (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                         (fourth_layer_hits_Ele_or_Pos[h] == 2 || fourth_layer_hits_Ele_or_Pos[h] ==3)){
                         if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) check124 = true;
                     }
                     if(!check124) continue;

                     TriggeringWithout_2ndPixel(k, j, h); // 134

                     if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                         (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3) &&
                         (fourth_layer_hits_Ele_or_Pos[h] == 1 || fourth_layer_hits_Ele_or_Pos[h] ==3)){
                         if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele && L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) check134 = true;
                     }
                     if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
                         (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3) &&
                         (fourth_layer_hits_Ele_or_Pos[h] == 2 || fourth_layer_hits_Ele_or_Pos[h] ==3)){
                         if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) check134 = true;
                     }  
                     if(!check134) continue;

                     TriggeringWithout_1stPixel(i, j, h);

                     if( (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                         (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3) &&
                         (fourth_layer_hits_Ele_or_Pos[h] == 1 || fourth_layer_hits_Ele_or_Pos[h] ==3)){
                         if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele && L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) check234 = true;
                     }
                     if( (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
                         (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3) &&
                         (fourth_layer_hits_Ele_or_Pos[h] == 2 || fourth_layer_hits_Ele_or_Pos[h] ==3)){
                         if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) check234 = true;
                     }   
                     if(!check234) continue;

                     PixTrkFourHits = true;                  
                  }// loop for 4th layer
               }// loop for 3rd layer
            }
         }
      }
*/
       // loop over every 3 out of 4 pixel combination 
       for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
          for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
              for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){
                 pass_count_EleorPos = 0;
                 pass_count_Ele = 0;
                 pass_count_Pos = 0;              
                 
                 // loop over every pixel hits in the given pixel combination
                 for( int k=0; k < layers[*first_hit]; k++){
                    for( int i=0; i < layers[*second_hit]; i++){
                        _pass_Ele = 0, _pass_Pos = 0;
                        L012_pass_Ele = 0, L012_pass_Pos = 0;
                        L013_pass_Ele = 0, L013_pass_Pos = 0;
                        L023_pass_Ele = 0, L023_pass_Pos = 0;
                      
                        // check L012 dphi, L013 dphi, L023 dphi
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
                            withoutEM_pass_Pos = 0, withEM_pass_Pos = 0;

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
                            dR = sqrt( pow(dEta,2) + pow(dPhi,2) );

                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second hit is in the second pixel layer and the third hit is in the third layer. 
                                TriggeringWithout_4thPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                      //if(L012_pass_Ele && L123_pass_Ele && L12_EM_Ele && L23_EM_Ele)
                                         all_cut_pass_Ele = 1; 
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele)
                                         withoutEM_pass_Ele = 1;
                                      if(L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                         withEM_pass_Ele = 1;
                                }
 
                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L012_pass_Pos && L123_pass_Pos && L12_EM_Pos && L23_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) withEM_pass_Pos = 1;
                                 }

                               if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo4thPix = 1;
                                }

                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo4thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo4thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit loop
                                  cout << "without fourth layer: " << (k+1) * (i+1) * ( j+1) << endl;
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
                                    //if(L012_pass_Ele && L124_pass_Ele && L12_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L012_pass_Pos && L124_pass_Pos && L12_EM_Pos && L24_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) withEM_pass_Pos = 1;
                                }

                                if( all_cut_pass_Ele || all_cut_pass_Pos ){
                                  pass_count_wo3thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo3thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo3thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without third layer: " << (k+1) * (i+1) * ( j+1) << endl;
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
                                    //if(L013_pass_Ele && L134_pass_Ele && L13_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                }

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L013_pass_Pos && L134_pass_Pos && L13_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                 }

                                if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo2thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo2thPix++;
                                }

                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo2thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without second layer: " << (k+1) * (i+1) * ( j+1) << endl;
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
                                    //if(L023_pass_Ele && L234_pass_Ele && L23_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (second_layer_hits_Ele_or_Pos[k] == 2 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L023_pass_Pos && L234_pass_Pos && L23_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                }


                                if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo1thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo1thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo1thPix++;
                                }
                    
                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without first layer: " << (k+1) * (i+1) * ( j+1) << endl;
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  first_layer_missing = 1;
                   		 } 
                               }
                              }

                         if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1) {pass_count_EleorPos++; pass_count = 1;}
                         if( all_cut_pass_Ele == 1 ) {pass_count_Ele++;}
                         if( all_cut_pass_Pos == 1 ) pass_count_Pos++;
                         if( withoutEM_pass_Ele == 1 ) withoutEM_count_Ele = 1;
                         if( withEM_pass_Ele == 1 ) withEM_count_Ele = 1;
                       } // loop for third layer hits
                   } // loop for second layer hits       
                 } // loop for first layer hits

                 pass_Ele.push_back(pass_count_Ele);
                 pass_Pos.push_back(pass_count_Pos);
                 pass_ElePos.push_back(pass_count_EleorPos);
          }          
        }
      }
/*
     if( (fabs(EgEta) <= 1.4 && fabs(EgEta) > 1.3) && pass_count == false) {
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
                         
                        if( *first_hit == 3 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                         
                          TriggeringWith_3rd4thPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                       if( _pass_Ele == 1 || _pass_Pos == 1) { pass_count = 1;}
                   } // loop for second layer hits       
                 } // loop for first layer hits

        }
      }
     }  
*/
      if(PixTrkFourHits){
         PiXTRKbit_4Hits_ = PiXTRKbit_4Hits_| (bit1 << nth_eg_pix_deta);
      } 
 
      if( pass_count ){
          trigger_bit_width_ = trigger_bit_width_| (bit1 << nth_eg_pix_deta);
      }

      global_index_width++; // remove this variable
      }

      trigger_bit_width.push_back(trigger_bit_width_);
      trigger_4Hits_bit_width.push_back(PiXTRKbit_4Hits_);

      // remove 
      if( pass_count ){ 
          ntCl_match.push_back(true);
          all_cut_pass_eg = 1; 
      }
      else ntCl_match.push_back(false);

      if( pass_count_wo4thPix ){
          ntCl_match_wo4thPix.push_back(true);
      }
      else ntCl_match_wo4thPix.push_back(false);

      Npass_woEM_wo4thPix.push_back(woEM_pass_count_wo4thPix);
      Npass_wEM_wo4thPix.push_back(wEM_pass_count_wo4thPix);

      if( pass_count_wo3thPix ){
          ntCl_match_wo3thPix.push_back(true);
      }
      else ntCl_match_wo3thPix.push_back(false);

      Npass_woEM_wo3thPix.push_back(woEM_pass_count_wo3thPix);
      Npass_wEM_wo3thPix.push_back(wEM_pass_count_wo3thPix);

      if( pass_count_wo2thPix ){
          ntCl_match_wo2thPix.push_back(true);
      }
      else ntCl_match_wo2thPix.push_back(false);

      Npass_woEM_wo2thPix.push_back(woEM_pass_count_wo2thPix);
      Npass_wEM_wo2thPix.push_back(wEM_pass_count_wo2thPix);

      if( pass_count_wo1thPix ){
          ntCl_match_wo1thPix.push_back(true);
      }
      else ntCl_match_wo1thPix.push_back(false);

      Npass_woEM_wo1thPix.push_back(woEM_pass_count_wo1thPix);
      Npass_wEM_wo1thPix.push_back(wEM_pass_count_wo1thPix);

      if(withoutEM_count_Ele){
         withoutEM_match.push_back(true);
      }
      else withoutEM_match.push_back(false);

      if(withEM_count_Ele){
         withEM_match.push_back(true);
      }
      else withEM_match.push_back(false);
      // remove since this was just for only one signal window width

         /////////////////////////////////
         //
  } // end of egamma loop    

     // find egamma objects (HGCAL) passing pixtrk signal windows
     for( int q=0; q<cl3d_n; q++){ 

      //if(cl3d_coreshowerlength->at(q) < 3 || cl3d_coreshowerlength->at(q) > 18) continue;
      //if(cl3d_srrtot->at(q) < 0.002 || cl3d_srrtot->at(q) > 0.005) continue;
      //if(cl3d_maxlayer->at(q) < 8 || cl3d_maxlayer->at(q) > 20) continue;
      //if(cl3d_firstlayer->at(q) > 5) continue;

      if(cl3d_egid->at(q) != 1) continue;

      EgEt =cl3d_pt ->at(q);
      EgEta=cl3d_eta->at(q);
      EgPhi=cl3d_phi->at(q);

      float EgGx = cl3d_x->at(q);
      float EgGy = cl3d_y->at(q);
      float EgGz = cl3d_z->at(q);
      emvector.SetXYZ(EgGx,EgGy,EgGz);

      if(EgEt < 8) continue;

      if( fabs(EgEta) <= 0.8 ) eta_region =1;
      if( fabs(EgEta) <= 1.15 && fabs(EgEta) > 0.8 ) eta_region =2;
      if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 1.15 ) eta_region =3;
      if( fabs(EgEta) <= 1.7 && fabs(EgEta) > 1.4 ) eta_region =4;
      if( fabs(EgEta) <= 2.25 && fabs(EgEta) > 1.7 ) eta_region =5;
      if( fabs(EgEta) <= 2.7 && fabs(EgEta) > 2.25 ) eta_region =6;
      if( fabs(EgEta) <= 3.0 && fabs(EgEta) > 2.7 ) eta_region =7;

      if( eta_region < 2 || eta_region > 3 ) continue;
      if( fabs(EgEta) > 3. ) continue;
      eta_region = 0; // initialize variable 

      pass_egobjects_check = 1;
      ntnEg2++;
      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);

      isTrack_match.push_back(hgcal_isTrackMatched->at(q));
      chi2.push_back(hgcal_trackHighestPtCutChi2Chi2->at(q));
      track_dr.push_back(hgcal_trackmatchingdR->at(q));
      
      // set regin of interest
      // for η < 1.3, Δφ < 0.05 
      //if(eta_region==1)SetROI(2); 
      //else SetROI(eta_region);
      //SetROI(1);

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

      bool isFourHits = true;
      // check which pixel has hits
       for( int i=1; i < 5; i++){
          if( layers[i] != 0 ){
            hitted_layers.push_back(i);
          }
          else{
              isFourHits = false;
          }
       }
   
       // set pixtrk signal boundary
       //
       int global_index_width = 0;
       trigger_bit_width_ = 0x0;
       PiXTRKbit_4Hits_ = 0x0; 

       for(int nth_eg_pix_deta = 0; nth_eg_pix_deta < 9; nth_eg_pix_deta++){
       if(nth_eg_pix_deta != 0) continue;

       //SetSingalBoundary(eta_region);
       if(eta_region == 1) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);
       else if(eta_region == 2) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
       else if(eta_region == 6) SetSingalBoundary(5, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);
       else SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);

       // PixTRK algorithm 
       pass_count = 0;
       pass_count_wo4thPix = 0, pass_count_wo3thPix = 0, pass_count_wo2thPix = 0, pass_count_wo1thPix = 0;
       woEM_pass_count_wo4thPix = 0, woEM_pass_count_wo3thPix = 0, woEM_pass_count_wo2thPix = 0, woEM_pass_count_wo1thPix = 0;
       wEM_pass_count_wo4thPix = 0, wEM_pass_count_wo3thPix = 0, wEM_pass_count_wo2thPix = 0, wEM_pass_count_wo1thPix = 0;
       withoutEM_count_Ele = 0, withEM_count_Ele = 0;

       fourth_layer_missing = 0;
       third_layer_missing = 0;
       second_layer_missing = 0;
       first_layer_missing = 0;

       // loop over every 3 out of 4 pixel combination 
       for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
          for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
              for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){
                 pass_count_EleorPos = 0;
                 pass_count_Ele = 0;
                 pass_count_Pos = 0;              
                 
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
                            withoutEM_pass_Pos = 0, withEM_pass_Pos = 0;

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
                            dR = sqrt( pow(dEta,2) + pow(dPhi,2) );

                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second hit is in the second pixel layer and the third hit is in the third layer. 
                                TriggeringWithout_4thPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                      //if(L012_pass_Ele && L123_pass_Ele && L12_EM_Ele && L23_EM_Ele)
                                         all_cut_pass_Ele = 1; 
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele)
                                         withoutEM_pass_Ele = 1;
                                      if(L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                         withEM_pass_Ele = 1;
                                }
 
                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L012_pass_Pos && L123_pass_Pos && L12_EM_Pos && L23_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) withEM_pass_Pos = 1;
                                 }

                               if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo4thPix = 1;
                                }

                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo4thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo4thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit loop
                                  cout << "without fourth layer: " << (k+1) * (i+1) * ( j+1) << endl;
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
                                    //if(L012_pass_Ele && L124_pass_Ele && L12_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L012_pass_Pos && L124_pass_Pos && L12_EM_Pos && L24_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) withEM_pass_Pos = 1;
                                }

                                if( all_cut_pass_Ele || all_cut_pass_Pos ){
                                  pass_count_wo3thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo3thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo3thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without third layer: " << (k+1) * (i+1) * ( j+1) << endl;
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
                                    //if(L013_pass_Ele && L134_pass_Ele && L13_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                }

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L013_pass_Pos && L134_pass_Pos && L13_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                 }

                                if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo2thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo2thPix++;
                                }

                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo2thPix++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without second layer: " << (k+1) * (i+1) * ( j+1) << endl;
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
                                    //if(L023_pass_Ele && L234_pass_Ele && L23_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (second_layer_hits_Ele_or_Pos[k] == 2 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    //if(L023_pass_Pos && L234_pass_Pos && L23_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                }


                                if( all_cut_pass_Ele || all_cut_pass_Pos){
                                  pass_count_wo1thPix = 1;
                                }
                                if( withoutEM_pass_Ele || withoutEM_pass_Pos) {
                                  woEM_pass_count_wo1thPix++;
                                }
                                if( withEM_pass_Ele || withEM_pass_Pos) {
                                  wEM_pass_count_wo1thPix++;
                                }
                    
                                if(skip){
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 ){ // if pass exit the for loops
                                  cout << "without first layer: " << (k+1) * (i+1) * ( j+1) << endl;
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  first_layer_missing = 1;
                   		 } 
                               }
                              }

                         if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1) {pass_count_EleorPos++; pass_count = 1;}
                         if( all_cut_pass_Ele == 1 ) {pass_count_Ele++;}
                         if( all_cut_pass_Pos == 1 ) pass_count_Pos++;
                         if( withoutEM_pass_Ele == 1 ) withoutEM_count_Ele = 1;
                         if( withEM_pass_Ele == 1 ) withEM_count_Ele = 1;
                       } // loop for third layer hits
                   } // loop for second layer hits       
                 } // loop for first layer hits

                 pass_Ele.push_back(pass_count_Ele);
                 pass_Pos.push_back(pass_count_Pos);
                 pass_ElePos.push_back(pass_count_EleorPos); // FIXME: save information on pixel combination for example to see if an event is passed if only L123 is passed or all L123, 124, 134, 234 is passed etc. bit can be used to save these information.
          }          
        }
      }
/*
     if( (fabs(EgEta) <= 1.4 && fabs(EgEta) > 1.3) && pass_count == false) {
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

                        if( *first_hit == 3 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_3rd4thPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                       if( _pass_Ele == 1 || _pass_Pos == 1) { pass_count = 1;}
                   } // loop for second layer hits       
                 } // loop for first layer hits

        }
      }
     }
*/
      if( pass_count ){
          trigger_bit_width_ = trigger_bit_width_| (bit1 << nth_eg_pix_deta);
      }

      global_index_width++;
      }

      trigger_bit_width.push_back(trigger_bit_width_);
      trigger_4Hits_bit_width.push_back(PiXTRKbit_4Hits_);

      if( pass_count ){ 
          ntCl_match.push_back(true);
          all_cut_pass_eg = 1; 
      }
      else ntCl_match.push_back(false);

      if( pass_count_wo4thPix ){
          ntCl_match_wo4thPix.push_back(true);
      }
      else ntCl_match_wo4thPix.push_back(false);

      Npass_woEM_wo4thPix.push_back(woEM_pass_count_wo4thPix);
      Npass_wEM_wo4thPix.push_back(wEM_pass_count_wo4thPix);

      if( pass_count_wo3thPix ){
          ntCl_match_wo3thPix.push_back(true);
      }
      else ntCl_match_wo3thPix.push_back(false);

      Npass_woEM_wo3thPix.push_back(woEM_pass_count_wo3thPix);
      Npass_wEM_wo3thPix.push_back(wEM_pass_count_wo3thPix);

      if( pass_count_wo2thPix ){
          ntCl_match_wo2thPix.push_back(true);
      }
      else ntCl_match_wo2thPix.push_back(false);

      Npass_woEM_wo2thPix.push_back(woEM_pass_count_wo2thPix);
      Npass_wEM_wo2thPix.push_back(wEM_pass_count_wo2thPix);

      if( pass_count_wo1thPix ){
          ntCl_match_wo1thPix.push_back(true);
      }
      else ntCl_match_wo1thPix.push_back(false);

      Npass_woEM_wo1thPix.push_back(woEM_pass_count_wo1thPix);
      Npass_wEM_wo1thPix.push_back(wEM_pass_count_wo1thPix);

      if(withoutEM_count_Ele){
         withoutEM_match.push_back(true);
      }
      else withoutEM_match.push_back(false);

      if(withEM_count_Ele){
         withEM_match.push_back(true);
      }
      else withEM_match.push_back(false);

         /////////////////////////////////
         //
  } // end of egamma loop  

  if(pass_egobjects_check){ event_denominator = 1; FillCutFlow("EvtCut", 1.);}
  if(all_cut_pass_eg) event_nominator = 1; 
     pixtrk_tree->Fill();
 } // end of entries loop 
   file3->Write();
}

