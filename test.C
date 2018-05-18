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
   
   class track
   {
       public:
             float pos_x1; float pos_y1; float pos_z1;
             float pos_x2; float pos_y2; float pos_z2;
             float pos_x3; float pos_y3; float pos_z3;
	     int index;

	     track() 
	     { 
		 pos_x3 = 0; pos_y2 = 0; pos_z1 = 0; 
		 pos_x3 = 0; pos_y2 = 0; pos_z1 = 0; 
		 pos_x3 = 0; pos_y2 = 0; pos_z1 = 0; 
	     }
	     track(float a, float b, float c, float d, float e, float f, float g, float h, float i, int k) 
	     { 
		 pos_x3 = a; pos_x2 = b; pos_x1 = c; 
		 pos_y3 = d; pos_y2 = e; pos_y1 = f; 
		 pos_z3 = g; pos_z2 = h; pos_z1 = i; 
		 index = k;
	     }
	     
	     static bool comp3(const track &t1, const track &t2)
	     {
		 return ( t1.pos_x3 < t2.pos_x3 );
	     }
	     static bool uni3(const track &t1, const track &t2)
	     {
		 return( t1.pos_x3 == t2.pos_x3 );
	     }
	     
	     static bool comp2(const track &t1, const track &t2)
	     {
		 return ( t1.pos_x2 < t2.pos_x2 );
	     }
	     static bool uni2(const track &t1, const track &t2)
	     {
		 return( t1.pos_x2 == t2.pos_x2 );
	     }
	     
	     static bool comp1(const track &t1, const track &t2)
	     {
		 return ( t1.pos_x1 < t2.pos_x1 );
	     }
	     static bool uni1(const track &t1, const track &t2)
	     {
		 return( t1.pos_x1 == t2.pos_x1 );
	     }
   };

   TH1F *h1 = new TH1F("h1","E_{T} > 10GeV",30,0,30);
   TH1F *h2 = new TH1F("h2","E_{T} > 20GeV",30,0,30);
   TH1F *h3 = new TH1F("h3","1.3 < #eta < 1.6 && Et > 10 GeV",30,0,30);
   TH1F *h4 = new TH1F("h4","1.3 < #eta < 1.6 && Et > 20 GeV",30,0,30);


//  nentries = 100;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
   //for (Long64_t jentry=0; jentry<15;jentry++) { //nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;
      FillCutFlow("NoCut", 1.);
        
      EgN=egN;

      pass_egobjects_check = 0;
      ntnEg2 = 0;
      ntnEg3 = 0;
      ntEgEt.clear();
      ntEgEta.clear();
      ntEgPhi.clear();
      ntEgEt_iso.clear();
      ntEgEta_iso.clear();
      ntEgPhi_iso.clear();
      ntCl_match.clear();
      ntCl_iso_match.clear();
      only_iso_match.clear();
      withoutEM_match.clear();
      withEM_match.clear();
      pass_Ele.clear();
      pass_Pos.clear();
      pass_ElePos.clear();

      ntfirstPix.clear();
      ntsecondPix.clear();
      ntthirdPix.clear();
      ntfourthPix.clear();

      all_cut_pass_eg = 0;
      event_nominator = 0;
      event_denominator = 0;
   
      // cout flow
      float EgEtCut  = 0;
      for(int k=0; k<EgN; k++) {
        EgEt =egEt ->at(k);
        if(EgEt > 10) EgEtCut = 1;
      }
      //if(EgEtCut == 0) continue;
      //FillCutFlow("MinEtCut", 1.);
      if(EgEtCut == 1) FillCutFlow("MinEtCut", 1.);

      int EtaCutFlow = 0;
      for(int k=0; k<EgN; k++) {
        EgEta=egEta->at(k);
        EgEt =egEt ->at(k);
        if(fabs(EgEta) < 1.3 && EgEt > 10) EtaCutFlow = 1; // for first η region
        //if(fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) < 2.5 && EgEt > 10) EtaCutFlow = 1;
      }
      //if(EtaCutFlow == 0) continue;
      //FillCutFlow("EtaCut", 1.);
      if(EtaCutFlow == 1) FillCutFlow("EtaCut", 1.);

      int DRCutFlow = 0;
      for(int k=0; k<EgN; k++) {
        EgEt =egEt ->at(k);
        EgEta=egEta->at(k);
        EgPhi=egPhi->at(k);

        float delta_phi = 0., delta_eta = 0;
        delta_phi = propgenPartPhi->at(0)  - EgPhi;
        delta_eta = propgenPartEta->at(0)  - EgEta;

        //if( sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < dr_cut && fabs(EgEta) < 1.3 && EgEt > 10 ) DRCutFlow = 1;
        //if( sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 && EgEt > 10 ) DRCutFlow = 1;
        //if( sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 && EgEt > 10 ) DRCutFlow = 1;
        //if( sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 && EgEt > 10 ) DRCutFlow = 1;
        if( sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) < 2.5 && EgEt > 10 ) DRCutFlow = 1;
      }
      //if(DRCutFlow == 0) continue;
      //FillCutFlow("DRCut", 1.);
      if(DRCutFlow == 1) FillCutFlow("DRCut", 1.);
/*
      float PtErrorCutFlow = 0;
      for(int k=0; k<EgN; k++) {
        EgEt =egEt ->at(k);
        EgEta=egEta->at(k);
        EgPhi=egPhi->at(k);

        float delta_phi = 0., delta_eta = 0;
        delta_phi = propgenPartPhi->at(0)  - EgPhi;
        delta_eta = propgenPartEta->at(0)  - EgEta;

        //if( fabs(propgenPartPt->at(0) - EgEt)/propgenPartPt->at(0) < .5 && sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < dr_cut && fabs(EgEta) < 1.3 && EgEt > 10 ) PtErrorCutFlow = 1;
        //if( fabs(propgenPartPt->at(0) - EgEt)/propgenPartPt->at(0) < .5 && sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 && EgEt > 10 ) PtErrorCutFlow = 1;
        //if( fabs(propgenPartPt->at(0) - EgEt)/propgenPartPt->at(0) < .5 && sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 && EgEt > 10 ) PtErrorCutFlow = 1;
        //if( fabs(propgenPartPt->at(0) - EgEt)/propgenPartPt->at(0) < .5 && sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 && EgEt > 10 ) PtErrorCutFlow = 1;
        if( fabs(propgenPartPt->at(0) - EgEt)/propgenPartPt->at(0) < .5 && sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) < 0.1 && fabs(EgEta) < 2.5 && EgEt > 10 ) PtErrorCutFlow = 1;
      }
      if(PtErrorCutFlow == 0) continue;
      FillCutFlow("PtErrCut", 1.);
*/
      nt_genPhi = propgenPartPhi->at(0);
      nt_genEta = propgenPartEta->at(0);
      nt_genPt = propgenPartPt->at(0);


     // find egamma objects passing pixtrk signal windows
     for( int q=0; q<EgN; q++){ 
      cout << "Start " << q << "th L1 egamma" << endl;
      
      EgEt =egEt ->at(q);
      EgEta=egEta->at(q);
      EgPhi=egPhi->at(q);

      float EgGx = egGx->at(q);
      float EgGy = egGy->at(q);
      float EgGz = egGz->at(q);
      emvector.SetXYZ(EgGx,EgGy,EgGz);

      if(EgEt < 10 ) continue;

      if( fabs(EgEta) < 1.3 ) eta_region =1;
      if( fabs(EgEta) < 1.6 && fabs(EgEta) > 1.3 ) eta_region =2;
      //if( fabs(EgEta) < 1.9 && fabs(EgEta) > 1.6 ) eta_region =3;
      //if( fabs(EgEta) < 2.5 && fabs(EgEta) > 1.9 ) eta_region =4;
      //if( fabs(EgEta) < 2.8 && fabs(EgEta) > 2.5 ) eta_region =5;
      //if( fabs(EgEta) < 3.0 && fabs(EgEta) > 2.8 ) eta_region =6;
      //if( eta_region != 1 ) continue;
      //if( fabs(EgEta) > 2.5 ) continue;
      //if( fabs(EgEta) > 3.0 ) continue;
      if( fabs(EgEta) > 1.6 ) continue;

      Int_t flag123 = 0;
      Int_t flag124 = 0;
      Int_t flag134 = 0;
      Int_t flag234 = 0;

      Float_t recoPV = 0.;

      Float_t zp14 = 0.;
      Float_t zp13 = 0.;
      Float_t zp24 = 0.;

      Bool_t PixTRK_check = false;

      float delta_phi = 0., delta_eta = 0;
      delta_phi = propgenPartPhi->at(0) - EgPhi;
      delta_eta = propgenPartEta->at(0) - EgEta;
 
      if( sqrt( pow(delta_phi,2) + pow(delta_eta,2) ) > dr_cut ) continue;
      if( fabs(propgenPartPt->at(0) - EgEt)/propgenPartPt->at(0) > .5 ) continue;

      pass_egobjects_check = 1;
      ntnEg2++;
      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);
      ntEgPhi.push_back(EgPhi);
      
      // set regin of interest
      // for η < 1.3, Δφ < 0.05 
      //
      if( eta_region == 1) SetROI(2);
      else SetROI(eta_region); 
      //SetROI(4); 

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
   
       // set pixtrk signal boundary
       if( eta_region == 1 ) SetSingalBoundary(2);
       else SetSingalBoundary(eta_region);
       //SetSingalBoundary(4);

       // PixTRK algorithm 
       pass_count = 0;
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

                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 )
                                {
                                    //if( all_cut_pass_Ele && !all_cut_pass_Pos ) cout << "  L123 Electron passed" << endl; 
                                    //if( !all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L123 Positron passed" << endl; 
                                    //if( all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L123 Ele & Pos passed" << endl; 
                                    flag123++;
                                    Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    Float_t R3 = sqrt(pow(third_layer_hits[j].X(),2)+pow(third_layer_hits[j].Y(),2));
                                    Float_t Z1 = first_layer_hits[k].Z();
                                    Float_t Z3 = third_layer_hits[j].Z();

                                    zp13 = (R3*Z1 - R1*Z3) / (R3 - R1);
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
                                
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 )
                                {
                                    //if( all_cut_pass_Ele && !all_cut_pass_Pos ) cout << "  L124 Electron passed" << endl; 
                                    //if( !all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L124 Positron passed" << endl; 
                                    //if( all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L124 Ele & Pos passed" << endl; 
                                    flag124++;
                                    Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    Float_t Z1 = first_layer_hits[k].Z();
                                    Float_t Z4 = fourth_layer_hits[j].Z();

                                    zp14 = (R4*Z1 - R1*Z4) / (R4 - R1);
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
                                
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 )
                                {
                                    //if( all_cut_pass_Ele && !all_cut_pass_Pos ) cout << "  L134 Electron passed" << endl; 
                                    //if( !all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L134 Positron passed" << endl; 
                                    //if( all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L134 Ele & Pos passed" << endl; 
                                    flag134++;
                                    Float_t R1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    Float_t Z1 = first_layer_hits[k].Z();
                                    Float_t Z4 = fourth_layer_hits[j].Z();

                                    zp14 = (R4*Z1 - R1*Z4) / (R4 - R1);
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
                                
                                if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1 )
                                {
                                    if( all_cut_pass_Ele && !all_cut_pass_Pos ) cout << "  L234 Electron passed" << endl; 
                                    if( !all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L234 Positron passed" << endl; 
                                    if( all_cut_pass_Ele && all_cut_pass_Pos ) cout << "  L234 Ele & Pos passed" << endl; 
                                    flag234++;
                                    Float_t R2 = sqrt(pow(second_layer_hits[k].X(),2)+pow(second_layer_hits[k].Y(),2));
                                    Float_t R4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    Float_t Z2 = second_layer_hits[k].Z();
                                    Float_t Z4 = fourth_layer_hits[j].Z();

                                    zp24 = (R4*Z2 - R2*Z4) / (R4 - R2);
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

                         if( all_cut_pass_Ele == 1 || all_cut_pass_Pos == 1) {pass_count_EleorPos++;}
                         if( all_cut_pass_Ele == 1 ) {pass_count_Ele++; pass_count = 1;}
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
      if( pass_count ){ 
          ntCl_match.push_back(true);
          all_cut_pass_eg = 1;
	  PixTRK_check = true;
      }
      else ntCl_match.push_back(false);

      if(withoutEM_count_Ele){
         withoutEM_match.push_back(true);
      }
      else withoutEM_match.push_back(false);

      if(withEM_count_Ele){
         withEM_match.push_back(true);
      }
      else withEM_match.push_back(false);


      //////////////// Track isolation algorithm ////////////////

      //cout << "Pass L123: " << flag123 << ", pass L124: " << flag124 << ", pass L134: " << flag134 << ", pass L234: " << flag234 << endl;

      //if( EgEt < 10 || fabs(EgEta) > 1.3 ) continue;
      if( flag123 == 0 && flag124 == 0 && flag134 == 0 && flag234 == 0 ) 
      {
          cout << "   This " << q << "th egamma passed Eta and Et cut but there are no matched pixel combinations" << endl;
          cout << endl;
          ntCl_iso_match.push_back(false);
          continue;
      }

      Bool_t isL1234 = false;
      Bool_t isL123D1 = false;
      Bool_t isL12D12 = false;
      Bool_t isL1D123 = false;

      if( eta_region = 1 ) isL1234 = true; 
      if( eta_region = 2 ) isL123D1 = true;
      if( eta_region = 3 ) isL12D12 = true;
      if( eta_region = 4 ) isL1D123 = true;

      cout << "    Let's start track isolation" << endl;
      //cout << "        EgEt: " << EgEt << ", EgEta: " << EgEta << endl;

      ntnEg3++;
      ntEgEt_iso.push_back(EgEt);
      ntEgEta_iso.push_back(EgEta);
      ntEgPhi_iso.push_back(EgPhi);

      if( flag124 || flag134 ) { recoPV = zp14; } /// L124 or L134 reconstruct 
      if( !flag124 && !flag134 && flag123 ) { recoPV = zp13; } /// L123 reconstruct
      if( !flag124 && !flag134 && !flag123 && flag234 ) { recoPV = zp24; } /// L234 reconstruct

      //cout << "        Gen pv: " << simVz->at(0) << ", recoPV: " << recoPV << endl;

      vector<TVector3> L1;
      vector<TVector3> L2;
      vector<TVector3> L3;
      vector<TVector3> L4;

      vector<TVector3> D1_L1;
      vector<TVector3> D1_L2;
      vector<TVector3> D1_L3;
      vector<TVector3> D1_L4;
      
      L1.clear();
      L2.clear();
      L3.clear();
      L4.clear();
      
      D1_L1.clear();
      D1_L2.clear();
      D1_L3.clear();
      D1_L4.clear();

      for(Int_t i = 0; i < bHitN; i++)
      {
          if(bHitGz->at(i) > 28.5) continue;
          Float_t R = sqrt(pow(bHitGx->at(i), 2)+ pow(bHitGy->at(i), 2));
          TVector3 pixel;	     
          pixel.SetXYZ(bHitGx->at(i), bHitGy->at(i), bHitGz->at(i) - recoPV);
          Float_t pixelPhi = pixel.Phi();
          Float_t pixelEta = pixel.Eta();
          Float_t deltaPhi = EgPhi - pixelPhi;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
          Float_t deltaEta = EgEta - pixelEta;
          Float_t deltaR = sqrt(pow(deltaPhi, 2) + pow(deltaEta, 2));
          if( fabs(deltaPhi) > 0.3 ) continue;
          
          if( R < 5.5 && deltaR < 0.3 )            L1.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
          if( R > 5.5 && R < 8.5 && deltaR < 0.3 ) L2.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
          if( R > 8.5 && R < 13 && deltaR < 0.3 )  L3.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
          if( R > 13 && R < 18 && deltaR < 0.3 )   L4.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
      }

      for(Int_t j = 0; j < fHitN; j++)
      {
          Float_t R = sqrt(pow(fHitGx->at(j),2)+pow(fHitGy->at(j),2));
          Float_t Z = fabs(fHitGz->at(j));
          TVector3 diskPixel;
          diskPixel.SetXYZ( fHitGx->at(j), fHitGy->at(j), fHitGz->at(j) - recoPV );
          Float_t diskPhi = diskPixel.Phi();
          Float_t diskEta = diskPixel.Eta();
          Float_t deltaPhi = EgPhi - diskPhi;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
          Float_t deltaEta = EgEta - diskEta;
          Float_t deltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta,2));
          if( fabs(deltaPhi) > 0.3 ) continue;

          if( Z > 28. && Z < 36. && deltaR < 0.3  D1.push_back(TVector3(fHitGx->at(j),fHitGy->at(j),fHitGz->at(j)));
      }

      // Vectors are declared for comparing x-positions between 3 out of 4 combinations //
      vector<track> L123_region1;
      L123_region1.clear();
      vector<track> L124_region1;
      L124_region1.clear();
      vector<track> L134_region1;
      L134_region1.clear();
      vector<track> L234_region1;
      L234_region1.clear();

      if( isL1234 )
      {
      /////////////////////// 1 2 3 ///////////////////////
      // {{{
      
      Float_t dz_L13L12 = 0.048;
      Float_t dz_L13L13 = 0.037;
      Float_t dz_L13L23 = 0.06;
      
      Float_t dEta_L123cut1 = 0.00486; // L12 - L23
      Float_t dEta_L123cut2 = 0.00312; // L12 - L13
      Float_t dEta_L123cut3 = 0.00300; // L13 - L23
      
      Float_t dEta_L123PVcut1 = 0.00162; // PVL3 - PVL2
      Float_t dEta_L123PVcut2 = 0.00312; // PVL3 - PVL1

      for(std::vector<TVector3>::iterator a = L1.begin(); a != L1.end(); ++a)
      {
          TVector3 PVL1;
	  PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL1 = PVL1.Phi();
	  Float_t etaPVL1 = PVL1.PseudoRapidity();
	  
	  Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z1 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L2.begin(); b != L2.end(); ++b)
	  {
	      TVector3 L1L2;
	      L1L2.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL1L2 = L1L2.Phi();
	      Float_t etaL1L2 = L1L2.PseudoRapidity();
	      
	      TVector3 PVL2;
	      PVL2.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL2 = PVL2.PseudoRapidity();
	      
	      Float_t R2 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z2 = (*b).Z();
	      Float_t pz12 = ( R2*Z1 - R1*Z2 ) / ( R2 - R1 );
	      Float_t dz12 = pz12 - recoPV;
	      
	      if( fabs(dz12) < dz_L13L12 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L3.begin(); c != L3.end(); ++c)
		  {
		      TVector3 L1L3;
		      L1L3.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL1L3 = L1L3.PseudoRapidity();
		      
		      TVector3 L2L3;
		      L2L3.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL2L3 = L2L3.Phi();
		      Float_t etaL2L3 = L2L3.PseudoRapidity();
		      
		      TVector3 PVL3;
		      PVL3.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL3 = PVL3.PseudoRapidity();
		      
		      Float_t R3 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z3 = (*c).Z();
		      
		      Float_t pz13 = ( R3*Z1 - R1*Z3 ) / ( R3 - R1 );
		      Float_t dz13 = pz13 - recoPV;
		      
		      Float_t pz23 = ( R3*Z2 - R2*Z3 ) / ( R3 - R2 );
		      Float_t dz23 = pz23 - recoPV;

		      Float_t dPhi1 = phiPVL1 - phiL1L2;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL1L2 - phiL2L3;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz13) < dz_L13L13 && fabs(dz23) < dz_L13L23 ) // Delta z cut
			  if( fabs(etaL2L3-etaL1L2) < dEta_L123cut1 && fabs(etaL1L3-etaL1L2) < dEta_L123cut2 && fabs(etaL2L3-etaL1L3) < dEta_L123cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL3 - etaPVL2) < dEta_L123PVcut1 && fabs(etaPVL3 - etaPVL1) < dEta_L123PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.05 && cut < 0.08 )
				      {
                          L123_region1.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 1));
				      }
		  }
	      }
	  }
      }
      
      // }}}
      
      /////////////////////// 1 2 4 ///////////////////////
      // {{{ 
      
      Float_t dz_L14L12 = 0.047;
      Float_t dz_L14L14 = 0.023;
      Float_t dz_L14L24 = 0.031;
      
      Float_t dEta_L124cut1 = 0.00448; // L12 - L24
      Float_t dEta_L124cut2 = 0.00324; // L12 - L14
      Float_t dEta_L124cut3 = 0.00138; // L14 - L24
      
      Float_t dEta_L124PVcut1 = 0.00142; // PVL4 - PVL2
      Float_t dEta_L124PVcut2 = 0.00300; // PVL4 - PVL1

      for(std::vector<TVector3>::iterator a = L1.begin(); a != L1.end(); ++a)
      {
          TVector3 PVL1;
	  PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL1 = PVL1.Phi();
	  Float_t etaPVL1 = PVL1.PseudoRapidity();
	  
	  Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z1 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L2.begin(); b != L2.end(); ++b)
	  {
	      TVector3 L1L2;
	      L1L2.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL1L2 = L1L2.Phi();
	      Float_t etaL1L2 = L1L2.PseudoRapidity();
	      
	      TVector3 PVL2;
	      PVL2.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL2 = PVL2.PseudoRapidity();
	      
	      Float_t R2 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z2 = (*b).Z();
	      Float_t pz12 = ( R2*Z1 - R1*Z2 ) / ( R2 - R1 );
	      Float_t dz12 = pz12 - recoPV;
	      
	      if( fabs(dz12) < dz_L14L12 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L4.begin(); c != L4.end(); ++c)
		  {
		      TVector3 L1L4;
		      L1L4.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL1L4 = L1L4.PseudoRapidity();
		      
		      TVector3 L2L4;
		      L2L4.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL2L4 = L2L4.Phi();
		      Float_t etaL2L4 = L2L4.PseudoRapidity();
		      
		      TVector3 PVL4;
		      PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL4 = PVL4.PseudoRapidity();
		      
		      Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z4 = (*c).Z();
		      
		      Float_t pz14 = ( R4*Z1 - R1*Z4 ) / ( R4 - R1 );
		      Float_t dz14 = pz14 - recoPV;
		      
		      Float_t pz24 = ( R4*Z2 - R2*Z4 ) / ( R4 - R2 );
		      Float_t dz24 = pz24 - recoPV;

		      Float_t dPhi1 = phiPVL1 - phiL1L2;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL1L2 - phiL2L4;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz14) < dz_L14L14 && fabs(dz24) < dz_L14L24 ) // Delta z cut
			  if( fabs(etaL2L4-etaL1L2) < dEta_L124cut1 && fabs(etaL1L4-etaL1L2) < dEta_L124cut2 && fabs(etaL2L4-etaL1L4) < dEta_L124cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL4 - etaPVL2) < dEta_L124PVcut1 && fabs(etaPVL4 - etaPVL1) < dEta_L124PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.04 && cut < 0.1 )
				      {
                          L124_region1.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 2));
				      }
		  }
	      }
	  }
      }
      // }}}
      
      /////////////////////// 1 3 4 ///////////////////////
      // {{{ 
      
      Float_t dz_L14L13 = 0.03;
      Float_t dz_L14L34 = 0.056;
      
      Float_t dEta_L134cut1 = 0.00404; // L13 - L34
      Float_t dEta_L134cut2 = 0.00152; // L13 - L14
      Float_t dEta_L134cut3 = 0.00252; // L14 - L34
      
      Float_t dEta_L134PVcut1 = 0.00100; // PVL4 - PVL3
      Float_t dEta_L134PVcut2 = 0.00300; // PVL4 - PVL1
      
      for(std::vector<TVector3>::iterator a = L1.begin(); a != L1.end(); ++a)
      {
          TVector3 PVL1;
	  PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL1 = PVL1.Phi();
	  Float_t etaPVL1 = PVL1.PseudoRapidity();
	  
	  Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z1 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L3.begin(); b != L3.end(); ++b)
	  {
	      TVector3 L1L3;
	      L1L3.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL1L3 = L1L3.Phi();
	      Float_t etaL1L3 = L1L3.PseudoRapidity();
	      
	      TVector3 PVL3;
	      PVL3.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL3 = PVL3.PseudoRapidity();
	      
	      Float_t R3 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z3 = (*b).Z();
	      Float_t pz13 = ( R3*Z1 - R1*Z3 ) / ( R3 - R1 );
	      Float_t dz13 = pz13 - recoPV;
	      
	      if( fabs(dz13) < dz_L14L13 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L4.begin(); c != L4.end(); ++c)
		  {
		      TVector3 L1L4;
		      L1L4.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL1L4 = L1L4.PseudoRapidity();
		      
		      TVector3 L3L4;
		      L3L4.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL3L4 = L3L4.Phi();
		      Float_t etaL3L4 = L3L4.PseudoRapidity();
		      
		      TVector3 PVL4;
		      PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL4 = PVL4.PseudoRapidity();
		      
		      Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z4 = (*c).Z();
		      
		      Float_t pz14 = ( R4*Z1 - R1*Z4 ) / ( R4 - R1 );
		      Float_t dz14 = pz14 - recoPV;
		      
		      Float_t pz34 = ( R4*Z3 - R3*Z4 ) / ( R4 - R3 );
		      Float_t dz34 = pz34 - recoPV;

		      Float_t dPhi1 = phiPVL1 - phiL1L3;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL1L3 - phiL3L4;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz14) < dz_L14L14 && fabs(dz34) < dz_L14L34 ) // Delta z cut
			  if( fabs(etaL3L4-etaL1L3) < dEta_L134cut1 && fabs(etaL1L4-etaL1L3) < dEta_L134cut2 && fabs(etaL3L4-etaL1L4) < dEta_L134cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL4 - etaPVL3) < dEta_L134PVcut1 && fabs(etaPVL4 - etaPVL1) < dEta_L134PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.04 && cut < 0.08 )
				      {
                          L134_region1.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 3));
				      }
		  }
	      }
	  }
      }
      // }}}

      /////////////////////// 2 3 4 ///////////////////////
      // {{{
      
      Float_t dz_L24L23 = 0.06;
      Float_t dz_L24L24 = 0.038;
      Float_t dz_L24L34 = 0.058;
      
      Float_t dEta_L234cut1 = 0.00417; // L23 - L34
      Float_t dEta_L234cut2 = 0.00296; // L23 - L24
      Float_t dEta_L234cut3 = 0.00238; // L24 - L34
      
      Float_t dEta_L234PVcut1 = 0.00100; // PVL4 - PVL3
      Float_t dEta_L234PVcut2 = 0.00142; // PVL4 - PVL2

      for(std::vector<TVector3>::iterator a = L2.begin(); a != L2.end(); ++a)
      {
          TVector3 PVL2;
	  PVL2.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL2 = PVL2.Phi();
	  Float_t etaPVL2 = PVL2.PseudoRapidity();
	  
	  Float_t R2 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z2 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L3.begin(); b != L3.end(); ++b)
	  {
	      TVector3 L2L3;
	      L2L3.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL2L3 = L2L3.Phi();
	      Float_t etaL2L3 = L2L3.PseudoRapidity();
	      
	      TVector3 PVL3;
	      PVL3.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL3 = PVL3.PseudoRapidity();
	      
	      Float_t R3 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z3 = (*b).Z();
	      Float_t pz23 = ( R3*Z2 - R2*Z3 ) / ( R3 - R2 );
	      Float_t dz23 = pz23 - recoPV;
	      
	      if( fabs(dz23) < dz_L24L23 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L4.begin(); c != L4.end(); ++c)
		  {
		      TVector3 L2L4;
		      L2L4.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL2L4 = L2L4.PseudoRapidity();
		      
		      TVector3 L3L4;
		      L3L4.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL3L4 = L3L4.Phi();
		      Float_t etaL3L4 = L3L4.PseudoRapidity();
		      
		      TVector3 PVL4;
		      PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL4 = PVL4.PseudoRapidity();
		      
		      Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z4 = (*c).Z();
		      
		      Float_t pz24 = ( R4*Z2 - R2*Z4 ) / ( R4 - R2 );
		      Float_t dz24 = pz24 - recoPV;
		      
		      Float_t pz34 = ( R4*Z3 - R3*Z4 ) / ( R4 - R3 );
		      Float_t dz34 = pz34 - recoPV;

		      Float_t dPhi1 = phiPVL2 - phiL2L3;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL2L3 - phiL3L4;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz24) < dz_L24L24 && fabs(dz34) < dz_L24L34 ) // Delta z cut
			  if( fabs(etaL3L4-etaL2L3) < dEta_L234cut1 && fabs(etaL2L4-etaL2L3) < dEta_L234cut2 && fabs(etaL3L4-etaL2L4) < dEta_L234cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL4 - etaPVL3) < dEta_L234PVcut1 && fabs(etaPVL4 - etaPVL2) < dEta_L234PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.06 && cut < 0.06 )
				      {
                          L234_region1.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 4));
				      }
		  }
	      }
	  }
      }
      // }}} 
      }
      
      vector<track> L123_region2;
      L123_region2.clear();
      vector<track> L12D1_region2;
      L12D1_region2.clear();
      vector<track> L13D1_region2;
      L13D1_region2.clear();
      vector<track> L23D1_region2;
      L23D1_region2.clear();

      if( isL123D1 )
      {
      /////////////////////// 1 2 3 ///////////////////////
      // {{{
      
      Float_t dz_L13L12 = 0.048;
      Float_t dz_L13L13 = 0.037;
      Float_t dz_L13L23 = 0.06;
      
      Float_t dEta_L123cut1 = 0.00486; // L12 - L23
      Float_t dEta_L123cut2 = 0.00312; // L12 - L13
      Float_t dEta_L123cut3 = 0.00300; // L13 - L23
      
      Float_t dEta_L123PVcut1 = 0.00162; // PVL3 - PVL2
      Float_t dEta_L123PVcut2 = 0.00312; // PVL3 - PVL1

      for(std::vector<TVector3>::iterator a = L1.begin(); a != L1.end(); ++a)
      {
          TVector3 PVL1;
	  PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL1 = PVL1.Phi();
	  Float_t etaPVL1 = PVL1.PseudoRapidity();
	  
	  Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z1 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L2.begin(); b != L2.end(); ++b)
	  {
	      TVector3 L1L2;
	      L1L2.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL1L2 = L1L2.Phi();
	      Float_t etaL1L2 = L1L2.PseudoRapidity();
	      
	      TVector3 PVL2;
	      PVL2.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL2 = PVL2.PseudoRapidity();
	      
	      Float_t R2 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z2 = (*b).Z();
	      Float_t pz12 = ( R2*Z1 - R1*Z2 ) / ( R2 - R1 );
	      Float_t dz12 = pz12 - recoPV;
	      
	      if( fabs(dz12) < dz_L13L12 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L3.begin(); c != L3.end(); ++c)
		  {
		      TVector3 L1L3;
		      L1L3.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL1L3 = L1L3.PseudoRapidity();
		      
		      TVector3 L2L3;
		      L2L3.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL2L3 = L2L3.Phi();
		      Float_t etaL2L3 = L2L3.PseudoRapidity();
		      
		      TVector3 PVL3;
		      PVL3.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL3 = PVL3.PseudoRapidity();
		      
		      Float_t R3 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z3 = (*c).Z();
		      
		      Float_t pz13 = ( R3*Z1 - R1*Z3 ) / ( R3 - R1 );
		      Float_t dz13 = pz13 - recoPV;
		      
		      Float_t pz23 = ( R3*Z2 - R2*Z3 ) / ( R3 - R2 );
		      Float_t dz23 = pz23 - recoPV;

		      Float_t dPhi1 = phiPVL1 - phiL1L2;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL1L2 - phiL2L3;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz13) < dz_L13L13 && fabs(dz23) < dz_L13L23 ) // Delta z cut
			  if( fabs(etaL2L3-etaL1L2) < dEta_L123cut1 && fabs(etaL1L3-etaL1L2) < dEta_L123cut2 && fabs(etaL2L3-etaL1L3) < dEta_L123cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL3 - etaPVL2) < dEta_L123PVcut1 && fabs(etaPVL3 - etaPVL1) < dEta_L123PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.05 && cut < 0.08 )
				      {
                          L123_region2.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 1));
				      }
		  }
	      }
	  }
      }
      
      // }}}
      
      /////////////////////// 1 2 D1 ///////////////////////
      // {{{ 
      
      Float_t dz_L14L12 = 0.047;
      Float_t dz_L14L14 = 0.023;
      Float_t dz_L14L24 = 0.031;
      
      Float_t dEta_L124cut1 = 0.00448; // L12 - L24
      Float_t dEta_L124cut2 = 0.00324; // L12 - L14
      Float_t dEta_L124cut3 = 0.00138; // L14 - L24
      
      Float_t dEta_L124PVcut1 = 0.00142; // PVL4 - PVL2
      Float_t dEta_L124PVcut2 = 0.00300; // PVL4 - PVL1

      for(std::vector<TVector3>::iterator a = L1.begin(); a != L1.end(); ++a)
      {
          TVector3 PVL1;
	  PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL1 = PVL1.Phi();
	  Float_t etaPVL1 = PVL1.PseudoRapidity();
	  
	  Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z1 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L2.begin(); b != L2.end(); ++b)
	  {
	      TVector3 L1L2;
	      L1L2.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL1L2 = L1L2.Phi();
	      Float_t etaL1L2 = L1L2.PseudoRapidity();
	      
	      TVector3 PVL2;
	      PVL2.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL2 = PVL2.PseudoRapidity();
	      
	      Float_t R2 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z2 = (*b).Z();
	      Float_t pz12 = ( R2*Z1 - R1*Z2 ) / ( R2 - R1 );
	      Float_t dz12 = pz12 - recoPV;
	      
	      if( fabs(dz12) < dz_L14L12 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L4.begin(); c != L4.end(); ++c)
		  {
		      TVector3 L1L4;
		      L1L4.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL1L4 = L1L4.PseudoRapidity();
		      
		      TVector3 L2L4;
		      L2L4.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL2L4 = L2L4.Phi();
		      Float_t etaL2L4 = L2L4.PseudoRapidity();
		      
		      TVector3 PVL4;
		      PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL4 = PVL4.PseudoRapidity();
		      
		      Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z4 = (*c).Z();
		      
		      Float_t pz14 = ( R4*Z1 - R1*Z4 ) / ( R4 - R1 );
		      Float_t dz14 = pz14 - recoPV;
		      
		      Float_t pz24 = ( R4*Z2 - R2*Z4 ) / ( R4 - R2 );
		      Float_t dz24 = pz24 - recoPV;

		      Float_t dPhi1 = phiPVL1 - phiL1L2;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL1L2 - phiL2L4;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz14) < dz_L14L14 && fabs(dz24) < dz_L14L24 ) // Delta z cut
			  if( fabs(etaL2L4-etaL1L2) < dEta_L124cut1 && fabs(etaL1L4-etaL1L2) < dEta_L124cut2 && fabs(etaL2L4-etaL1L4) < dEta_L124cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL4 - etaPVL2) < dEta_L124PVcut1 && fabs(etaPVL4 - etaPVL1) < dEta_L124PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.04 && cut < 0.1 )
				      {
                          L124_region1.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 2));
				      }
		  }
	      }
	  }
      }
      // }}}
      
      /////////////////////// 1 3 D1 ///////////////////////
      // {{{ 
      
      Float_t dz_L14L13 = 0.03;
      Float_t dz_L14L34 = 0.056;
      
      Float_t dEta_L134cut1 = 0.00404; // L13 - L34
      Float_t dEta_L134cut2 = 0.00152; // L13 - L14
      Float_t dEta_L134cut3 = 0.00252; // L14 - L34
      
      Float_t dEta_L134PVcut1 = 0.00100; // PVL4 - PVL3
      Float_t dEta_L134PVcut2 = 0.00300; // PVL4 - PVL1
      
      for(std::vector<TVector3>::iterator a = L1.begin(); a != L1.end(); ++a)
      {
          TVector3 PVL1;
	  PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL1 = PVL1.Phi();
	  Float_t etaPVL1 = PVL1.PseudoRapidity();
	  
	  Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z1 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L3.begin(); b != L3.end(); ++b)
	  {
	      TVector3 L1L3;
	      L1L3.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL1L3 = L1L3.Phi();
	      Float_t etaL1L3 = L1L3.PseudoRapidity();
	      
	      TVector3 PVL3;
	      PVL3.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL3 = PVL3.PseudoRapidity();
	      
	      Float_t R3 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z3 = (*b).Z();
	      Float_t pz13 = ( R3*Z1 - R1*Z3 ) / ( R3 - R1 );
	      Float_t dz13 = pz13 - recoPV;
	      
	      if( fabs(dz13) < dz_L14L13 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L4.begin(); c != L4.end(); ++c)
		  {
		      TVector3 L1L4;
		      L1L4.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL1L4 = L1L4.PseudoRapidity();
		      
		      TVector3 L3L4;
		      L3L4.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL3L4 = L3L4.Phi();
		      Float_t etaL3L4 = L3L4.PseudoRapidity();
		      
		      TVector3 PVL4;
		      PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL4 = PVL4.PseudoRapidity();
		      
		      Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z4 = (*c).Z();
		      
		      Float_t pz14 = ( R4*Z1 - R1*Z4 ) / ( R4 - R1 );
		      Float_t dz14 = pz14 - recoPV;
		      
		      Float_t pz34 = ( R4*Z3 - R3*Z4 ) / ( R4 - R3 );
		      Float_t dz34 = pz34 - recoPV;

		      Float_t dPhi1 = phiPVL1 - phiL1L3;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL1L3 - phiL3L4;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz14) < dz_L14L14 && fabs(dz34) < dz_L14L34 ) // Delta z cut
			  if( fabs(etaL3L4-etaL1L3) < dEta_L134cut1 && fabs(etaL1L4-etaL1L3) < dEta_L134cut2 && fabs(etaL3L4-etaL1L4) < dEta_L134cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL4 - etaPVL3) < dEta_L134PVcut1 && fabs(etaPVL4 - etaPVL1) < dEta_L134PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.04 && cut < 0.08 )
				      {
                          L134_region1.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 3));
				      }
		  }
	      }
	  }
      }
      // }}}

      /////////////////////// 2 3 D1 ///////////////////////
      // {{{
      
      Float_t dz_L24L23 = 0.06;
      Float_t dz_L24L24 = 0.038;
      Float_t dz_L24L34 = 0.058;
      
      Float_t dEta_L234cut1 = 0.00417; // L23 - L34
      Float_t dEta_L234cut2 = 0.00296; // L23 - L24
      Float_t dEta_L234cut3 = 0.00238; // L24 - L34
      
      Float_t dEta_L234PVcut1 = 0.00100; // PVL4 - PVL3
      Float_t dEta_L234PVcut2 = 0.00142; // PVL4 - PVL2

      for(std::vector<TVector3>::iterator a = L2.begin(); a != L2.end(); ++a)
      {
          TVector3 PVL2;
	  PVL2.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
	  Float_t phiPVL2 = PVL2.Phi();
	  Float_t etaPVL2 = PVL2.PseudoRapidity();
	  
	  Float_t R2 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
	  Float_t Z2 = (*a).Z();

	  for(std::vector<TVector3>::iterator b = L3.begin(); b != L3.end(); ++b)
	  {
	      TVector3 L2L3;
	      L2L3.SetXYZ( (*b).X() - (*a).X(), (*b).Y() - (*a).Y(), (*b).Z() - (*a).Z() );
	      Float_t phiL2L3 = L2L3.Phi();
	      Float_t etaL2L3 = L2L3.PseudoRapidity();
	      
	      TVector3 PVL3;
	      PVL3.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
	      Float_t etaPVL3 = PVL3.PseudoRapidity();
	      
	      Float_t R3 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
	      Float_t Z3 = (*b).Z();
	      Float_t pz23 = ( R3*Z2 - R2*Z3 ) / ( R3 - R2 );
	      Float_t dz23 = pz23 - recoPV;
	      
	      if( fabs(dz23) < dz_L24L23 ) 
	      {
		  for(std::vector<TVector3>::iterator c = L4.begin(); c != L4.end(); ++c)
		  {
		      TVector3 L2L4;
		      L2L4.SetXYZ( (*c).X() - (*a).X(), (*c).Y() - (*a).Y(), (*c).Z() - (*a).Z() );
		      Float_t etaL2L4 = L2L4.PseudoRapidity();
		      
		      TVector3 L3L4;
		      L3L4.SetXYZ( (*c).X() - (*b).X(), (*c).Y() - (*b).Y(), (*c).Z() - (*b).Z() );
		      Float_t phiL3L4 = L3L4.Phi();
		      Float_t etaL3L4 = L3L4.PseudoRapidity();
		      
		      TVector3 PVL4;
		      PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
		      Float_t etaPVL4 = PVL4.PseudoRapidity();
		      
		      Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
		      Float_t Z4 = (*c).Z();
		      
		      Float_t pz24 = ( R4*Z2 - R2*Z4 ) / ( R4 - R2 );
		      Float_t dz24 = pz24 - recoPV;
		      
		      Float_t pz34 = ( R4*Z3 - R3*Z4 ) / ( R4 - R3 );
		      Float_t dz34 = pz34 - recoPV;

		      Float_t dPhi1 = phiPVL2 - phiL2L3;
		      if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
		      if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

		      Float_t dPhi2 = phiL2L3 - phiL3L4;
		      if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
		      if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

		      Float_t cut = dPhi1 - dPhi2;

		      if( fabs(dz24) < dz_L24L24 && fabs(dz34) < dz_L24L34 ) // Delta z cut
			  if( fabs(etaL3L4-etaL2L3) < dEta_L234cut1 && fabs(etaL2L4-etaL2L3) < dEta_L234cut2 && fabs(etaL3L4-etaL2L4) < dEta_L234cut3  ) // Delta eta(pixel,pixel) cut
			      if( fabs(etaPVL4 - etaPVL3) < dEta_L234PVcut1 && fabs(etaPVL4 - etaPVL2) < dEta_L234PVcut2 ) // Delta eta(PV, pixel) cut
				  if( ( dPhi1 < 0 && dPhi2 < 0 ) || ( dPhi1 > 0 && dPhi2 > 0 ) ) // Delta phi sign cut
				      if( cut > -0.06 && cut < 0.06 )
				      {
                          L234_region1.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 4));
				      }
		  }
	      }
	  }
      }
      // }}} 
      }


      //cout << "L123: " << L123.size() << ", L124: " << L124.size() << ", L134: " << L134.size() << ", L234: " << L234.size() << endl;

      if( L123_region1.size() >= 2 ) 
      {
	  sort(L123_region1.begin(), L123_region1.end(), track::comp3);
	  L123_region1.erase(unique(L123_region1.begin(), L123_region1.end(), track::uni3),L123_region1.end());
	  sort(L123_region1.begin(), L123_region1.end(), track::comp2);
	  L123_region1.erase(unique(L123_region1.begin(), L123_region1.end(), track::uni2),L123_region1.end());
	  sort(L123_region1.begin(), L123_region1.end(), track::comp1);
	  L123_region1.erase(unique(L123_region1.begin(), L123_region1.end(), track::uni1),L123_region1.end());
      }
      if( L124_region1.size() >= 2 ) 
      {
	  sort(L124_region1.begin(), L124_region1.end(), track::comp3);
	  L124_region1.erase(unique(L124_region1.begin(), L124_region1.end(), track::uni3),L124_region1.end());
	  sort(L124_region1.begin(), L124_region1.end(), track::comp2);
	  L124_region1.erase(unique(L124_region1.begin(), L124_region1.end(), track::uni2),L124_region1.end());
	  sort(L124_region1.begin(), L124_region1.end(), track::comp1);
	  L124_region1.erase(unique(L124_region1.begin(), L124_region1.end(), track::uni1),L124_region1.end());
      }
      if( L134_region1.size() >= 2 ) 
      {
	  sort(L134_region1.begin(), L134_region1.end(), track::comp3);
	  L134_region1.erase(unique(L134_region1.begin(), L134_region1.end(), track::uni3),L134_region1.end());
	  sort(L134_region1.begin(), L134_region1.end(), track::comp2);
	  L134_region1.erase(unique(L134_region1.begin(), L134_region1.end(), track::uni2),L134_region1.end());
	  sort(L134_region1.begin(), L134_region1.end(), track::comp1);
	  L134_region1.erase(unique(L134_region1.begin(), L134_region1.end(), track::uni1),L134_region1.end());
      }
      if( L234_region1.size() >= 2 ) 
      {
	  sort(L234_region1.begin(), L234_region1.end(), track::comp3);
	  L234_region1.erase(unique(L234_region1.begin(), L234_region1.end(), track::uni3),L234_region1.end());
	  sort(L234_region1.begin(), L234_region1.end(), track::comp2);
	  L234_region1.erase(unique(L234_region1.begin(), L234_region1.end(), track::uni2),L234_region1.end());
	  sort(L234_region1.begin(), L234_region1.end(), track::comp1);
	  L234_region1.erase(unique(L234_region1.begin(), L234_region1.end(), track::uni1),L234_region1.end());
      }

      vector<track> all_region1;
      all_region1.clear();

      for(vector<track>::iterator a1 = L123_region1.begin(); a1 != L123_region1.end(); ++a1) all_region1.push_back(track((*a1).pos_x3, (*a1).pos_x2, (*a1).pos_x1, (*a1).pos_y3, (*a1).pos_y2, (*a1).pos_y1, (*a1).pos_z3, (*a1).pos_z2, (*a1).pos_z1, 1 ));
      for(vector<track>::iterator a2 = L124_region1.begin(); a2 != L124_region1.end(); ++a2) all_region1.push_back(track((*a2).pos_x3, (*a2).pos_x2, (*a2).pos_x1, (*a2).pos_y3, (*a2).pos_y2, (*a2).pos_y1, (*a2).pos_z3, (*a2).pos_z2, (*a2).pos_z1, 2 ));
      for(vector<track>::iterator a3 = L134_region1.begin(); a3 != L134_region1.end(); ++a3) all_region1.push_back(track((*a3).pos_x3, (*a3).pos_x2, (*a3).pos_x1, (*a3).pos_y3, (*a3).pos_y2, (*a3).pos_y1, (*a3).pos_z3, (*a3).pos_z2, (*a3).pos_z1, 3 ));
      for(vector<track>::iterator a4 = L234_region1.begin(); a4 != L234_region1.end(); ++a4) all_region1.push_back(track((*a4).pos_x3, (*a4).pos_x2, (*a4).pos_x1, (*a4).pos_y3, (*a4).pos_y2, (*a4).pos_y1, (*a4).pos_z3, (*a4).pos_z2, (*a4).pos_z1, 4 ));

      sort(all_region1.begin(), all_region1.end(), track::comp3);
      all_region1.erase(unique(all_region1.begin(), all_region1.end(), track::uni3),all_region1.end());
      sort(all_region1.begin(), all_region1.end(), track::comp2);
      all_region1.erase(unique(all_region1.begin(), all_region1.end(), track::uni2),all_region1.end());
      sort(all_region1.begin(), all_region1.end(), track::comp1);
      all_region1.erase(unique(all_region1.begin(), all_region1.end(), track::uni1),all_region1.end());

      //cout << "Number of tracks: " << all.size() << endl;
      //cout << "============================" << endl << endl;

      /*
      for(vector<track>::iterator a1 = all.begin(); a1 !=  all.end(); ++a1)
      {
	  cout << "Which combination: " << (*a1).index << endl;
	  cout << "3rd X: " << (*a1).pos_x3 << ", Y: " << (*a1).pos_y3 << ", Z: " << (*a1).pos_z3 << endl;
	  cout << "2nd X: " << (*a1).pos_x2 << ", Y: " << (*a1).pos_y2 << ", Z: " << (*a1).pos_z2 << endl;
	  cout << "1st X: " << (*a1).pos_x1 << ", Y: " << (*a1).pos_y1 << ", Z: " << (*a1).pos_z1 << endl;
          cout << endl;

	  Float_t R1 = sqrt(pow((*a1).pos_x1,2)+pow((*a1).pos_y1,2));
	  Float_t Z1 = (*a1).pos_z1;

	  Float_t R2 = sqrt(pow((*a1).pos_x2,2)+pow((*a1).pos_y2,2));
	  Float_t Z2 = (*a1).pos_z2;
		      
	  Float_t R3 = sqrt(pow((*a1).pos_x3,2)+pow((*a1).pos_y3,2));
	  Float_t Z3 = (*a1).pos_z3;
	      
	  Float_t pz12 = ( R2*Z1 - R1*Z2 ) / ( R2 - R1 );
	  Float_t dz12 = pz12 - recoPV;
	  
	  Float_t pz13 = ( R3*Z1 - R1*Z3 ) / ( R3 - R1 );
	  Float_t dz13 = pz13 - recoPV;
	  
	  Float_t pz23 = ( R3*Z2 - R2*Z3 ) / ( R3 - R2 );
	  Float_t dz23 = pz23 - recoPV;
	  
	  TVector3 seg1;
	  seg1.SetXYZ( (*a1).pos_x3 - (*a1).pos_x2, (*a1).pos_y3 - (*a1).pos_y2, (*a1).pos_z3 - (*a1).pos_z2);
	  Float_t phi1 = seg1.Phi();
	  Float_t eta1 = seg1.PseudoRapidity();
		      
	  TVector3 seg2;
	  seg2.SetXYZ( (*a1).pos_x2 - (*a1).pos_x1, (*a1).pos_y2 - (*a1).pos_y1, (*a1).pos_z2 - (*a1).pos_z1);
	  Float_t phi2 = seg2.Phi();
	  Float_t eta2 = seg2.PseudoRapidity();
	  
	  TVector3 seg3;
	  seg3.SetXYZ( (*a1).pos_x1, (*a1).pos_y1, (*a1).pos_z1 - recoPV );
	  Float_t phi3 = seg3.Phi();
	  Float_t eta3 = seg3.PseudoRapidity();
	  
	  TVector3 seg4;
	  seg4.SetXYZ( (*a1).pos_x3 - (*a1).pos_x1, (*a1).pos_y3 - (*a1).pos_y1, (*a1).pos_z3 - (*a1).pos_z1);
	  Float_t eta4 = seg4.PseudoRapidity();
	  
	  TVector3 seg5;
	  seg5.SetXYZ( (*a1).pos_x2, (*a1).pos_y2, (*a1).pos_z2 - recoPV );
	  Float_t eta5 = seg5.PseudoRapidity();
	  
	  TVector3 seg6;
	  seg6.SetXYZ( (*a1).pos_x3, (*a1).pos_y3, (*a1).pos_z3 - recoPV );
	  Float_t eta6 = seg6.PseudoRapidity();
	  
	  Float_t dPhi1 = phi1 - phi2;
	  if( dPhi1 >= TMath::Pi() ) dPhi1 -= 2.*TMath::Pi();
	  if( dPhi1 < -TMath::Pi() ) dPhi1 += 2.*TMath::Pi();

	  Float_t dPhi2 = phi2 - phi3;
	  if( dPhi2 >= TMath::Pi() ) dPhi2 -= 2.*TMath::Pi();
	  if( dPhi2 < -TMath::Pi() ) dPhi2 += 2.*TMath::Pi();

	  Float_t cut = dPhi1 - dPhi2;

	  cout << "Gen Phi: " << genPartPhi->at(0) << ", gen Eta: " << genPartEta->at(0) << ", gen pT: " << genPartPt->at(0) << endl;
	  cout << endl;
	  
	  if( (*a1).index == 1 ) 
	  {
	      cout << "======================= pixel information (L123)  ============================" << endl;
	      cout << "   Phi1: " << phi1 << ", Phi2: " << phi2 << ", Phi3: " << phi3 << endl;
	      cout << "   Eta1: " << eta1 << ", Eta2: " << eta2 << ", Eta3: " << eta4 << endl;
	      cout << endl;
	      cout << "   PV Eta1: " << eta4 << ", Eta2: " << eta5 << ", Eta3: " << eta3 << endl;
	      cout << endl; 
	      cout << "======================= delta information ============================" << endl;
	      cout << "   dZ1 value(1st, 2nd)(dZ1 cut: " << dz_L13L12 << ")= " << fabs(dz12) <<
	      ",   dZ2 value(1st, 3rd)(dZ2 cut: " << dz_L13L13 << ")= " << fabs(dz13) <<
	      ",   dZ3 value(2nd, 3rd)(dZ3 cut: " << dz_L13L23 << ")= " << fabs(dz23) << endl;
	      cout << "   Delta eta(pixel,pixel) dEta1(cut: " << dEta_L123cut1 << ")= " << fabs(eta1-eta2) <<   
	      ",   Delta eta(pixel,pixel) dEta2(cut: " << dEta_L123cut2 << ")= " << fabs(eta2-eta4) <<   
	      ",   Delta eta(pixel,pixel) dEta3(cut: " << dEta_L123cut3 << ")= " << fabs(eta1-eta4) << endl;  
              cout << "   Delta eta(PV, pixel) dEta1(cut: " << dEta_L123PVcut1 << ")= " << fabs(eta6-eta5) <<  
              ",   Delta eta(PV, pixel) dEta2(cut: " << dEta_L123PVcut2 << ")= " << fabs(eta6-eta3) << endl;  
	      cout << "   dPhi1: " << dPhi1 << ", dPhi2: " << dPhi2 << 
	      ", ddPhi: " << dPhi1 - dPhi2 << ", -0.05 < cut < 0.08" << endl;
	      cout << "===================================================" << endl;
	      cout << endl;
	  }
	  
	  if( (*a1).index == 2 ) 
	  {
	      cout << "======================= pixel information (L124) ============================" << endl;
	      cout << "   Phi1: " << phi1 << ", Phi2: " << phi2 << ", Phi3: " << phi3 << endl;
	      cout << "   Eta1: " << eta1 << ", Eta2: " << eta2 << ", Eta3: " << eta4 << endl;
	      cout << endl;
	      cout << "   PV Eta1: " << eta4 << ", Eta2: " << eta5 << ", Eta3: " << eta3 << endl;
	      cout << endl; 
	      cout << "======================= delta information ============================" << endl;
	      cout << "   dZ1 value(1st, 2nd)(dZ1 cut: " << dz_L14L12 << ")= " << fabs(dz12) <<
	      ",   dZ2 value(1st, 3rd)(dZ2 cut: " << dz_L14L14 << ")= " << fabs(dz13) <<
	      ",   dZ3 value(2nd, 3rd)(dZ3 cut: " << dz_L14L24 << ")= " << fabs(dz23) << endl;
	      cout << "   Delta eta(pixel,pixel) dEta1(cut: " << dEta_L124cut1 << ")= " << fabs(eta1-eta2) <<   
	      ",   Delta eta(pixel,pixel) dEta2(cut: " << dEta_L124cut2 << ")= " << fabs(eta2-eta4) <<   
	      ",   Delta eta(pixel,pixel) dEta3(cut: " << dEta_L124cut3 << ")= " << fabs(eta1-eta4) << endl;  
              cout << "   Delta eta(PV, pixel) dEta1(cut: " << dEta_L124PVcut1 << ")= " << fabs(eta6-eta5) <<  
              ",   Delta eta(PV, pixel) dEta2(cut: " << dEta_L124PVcut2 << ")= " << fabs(eta6-eta3) << endl; 
	      cout << "   dPhi1: " << dPhi1 << ", dPhi2: " << dPhi2 << 
	      ", ddPhi: " << dPhi1 - dPhi2 << ", -0.05 < cut < 0.08" << endl;
	      cout << "===================================================" << endl;
	      cout << endl;
	  }
	  
	  if( (*a1).index == 3 ) 
	  {
	      cout << "======================= pixel information (L134) ============================" << endl;
	      cout << "   Phi1: " << phi1 << ", Phi2: " << phi2 << ", Phi3: " << phi3 << endl;
	      cout << "   Eta1: " << eta1 << ", Eta2: " << eta2 << ", Eta3: " << eta4 << endl;
	      cout << endl;
	      cout << "   PV Eta1: " << eta4 << ", Eta2: " << eta5 << ", Eta3: " << eta3 << endl;
	      cout << endl; 
	      cout << "======================= delta information ============================" << endl;
	      cout << "   dZ1 value(1st, 2nd)(dZ1 cut: " << dz_L14L13 << ")= " << fabs(dz12) <<
	      ",   dZ2 value(1st, 3rd)(dZ2 cut: " << dz_L14L14 << ")= " << fabs(dz13) <<
	      ",   dZ3 value(2nd, 3rd)(dZ3 cut: " << dz_L14L34 << ")= " << fabs(dz23) << endl;
	      cout << "   Delta eta(pixel,pixel) dEta1(cut: " << dEta_L134cut1 << ")= " << fabs(eta1-eta2) <<   
	      ",   Delta eta(pixel,pixel) dEta2(cut: " << dEta_L134cut2 << ")= " << fabs(eta2-eta4) <<   
	      ",   Delta eta(pixel,pixel) dEta3(cut: " << dEta_L134cut3 << ")= " << fabs(eta1-eta4) << endl;  
              cout << "   Delta eta(PV, pixel) dEta1(cut: " << dEta_L134PVcut1 << ")= " << fabs(eta6-eta5) <<  
              ",   Delta eta(PV, pixel) dEta2(cut: " << dEta_L134PVcut2 << ")= " << fabs(eta6-eta3) << endl; 
	      cout << "   dPhi1: " << dPhi1 << ", dPhi2: " << dPhi2 << 
	      ", ddPhi: " << dPhi1 - dPhi2 << ", -0.05 < cut < 0.08" << endl;
	      cout << "===================================================" << endl;
	      cout << endl;
	  }
	  
	  if( (*a1).index == 4 ) 
	  {
	      cout << "======================= pixel information (L234) ============================" << endl;
	      cout << "   Phi1: " << phi1 << ", Phi2: " << phi2 << ", Phi3: " << phi3 << endl;
	      cout << "   Eta1: " << eta1 << ", Eta2: " << eta2 << ", Eta3: " << eta4 << endl;
	      cout << endl;
	      cout << "   PV Eta1: " << eta4 << ", Eta2: " << eta5 << ", Eta3: " << eta3 << endl;
	      cout << endl; 
	      cout << "======================= delta information ============================" << endl;
	      cout << "   dZ1 value(1st, 2nd)(dZ1 cut: " << dz_L24L23 << ")= " << fabs(dz12) <<
	      ",   dZ2 value(1st, 3rd)(dZ2 cut: " << dz_L24L24 << ")= " << fabs(dz13) <<
	      ",   dZ3 value(2nd, 3rd)(dZ3 cut: " << dz_L24L34 << ")= " << fabs(dz23) << endl;
	      cout << "   Delta eta(pixel,pixel) dEta1(cut: " << dEta_L234cut1 << ")= " << fabs(eta1-eta2) <<   
	      ",   Delta eta(pixel,pixel) dEta2(cut: " << dEta_L234cut2 << ")= " << fabs(eta2-eta4) <<   
	      ",   Delta eta(pixel,pixel) dEta3(cut: " << dEta_L234cut3 << ")= " << fabs(eta1-eta4) << endl;  
              cout << "   Delta eta(PV, pixel) dEta1(cut: " << dEta_L234PVcut1 << ")= " << fabs(eta6-eta5) <<  
              ",   Delta eta(PV, pixel) dEta2(cut: " << dEta_L234PVcut2 << ")= " << fabs(eta6-eta3) << endl; 
	      cout << "   dPhi1: " << dPhi1 << ", dPhi2: " << dPhi2 << 
	      ", ddPhi: " << dPhi1 - dPhi2 << ", -0.05 < cut < 0.08" << endl;
	      cout << "===================================================" << endl;
	      cout << endl;
	  }
      }
      */

      //cout << "        This e/gamma is (0: false, 1: true): " << PixTRK_check << ", and number of tracks: " << all.size() << endl;

      //if( all.size() <= 2 ) iso_count++;


      if( all_region1.size() <= 2 ) 
      {
              //cout << "           This " << q << "th egamma passed Eta, Et cut and track isolation" << endl;
              ntCl_iso_match.push_back(true);
      }
      else ntCl_iso_match.push_back(false);

      if( all.size_region1() <= 2 ) only_iso_match.push_back(true);
      else only_iso_match.push_back(false);

      h1->Fill( all_region1.size() );
      if( EgEt > 20. ) h2->Fill( all_region1.size() );


  } // end of egamma loop    
  if(pass_egobjects_check){ event_denominator = 1; FillCutFlow("EvtCut", 1.); }
  if(all_cut_pass_eg) event_nominator = 1; 
     pixtrk_tree->Fill();
 } // end of entries loop 
   file3->Write();
}
