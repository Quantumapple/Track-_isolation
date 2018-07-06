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

   useDR == 1;

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
  
   TH1F *L123_case1 = new TH1F("L123_case1","L123 with Bsl1-l1l3",1000,0,250);
   TH1F *L123_case2 = new TH1F("L123_case2","L123 with Bsl2-l2l3",1000,0,250);

   TH1F *L124_case1 = new TH1F("L124_case1","L124 with Bsl1-l1l4",1000,0,250);
   TH1F *L124_case2 = new TH1F("L124_case2","L124 with Bsl2-l2l4",1000,0,250);
   
   TH1F *L134_case1 = new TH1F("L134_case1","L134 with Bsl1-l1l4",1000,0,250);
   TH1F *L134_case2 = new TH1F("L134_case2","L134 with Bsl3-l3l4",1000,0,250);

   TH1F *L234_case1 = new TH1F("L234_case1","L234 with Bsl2-l2l4",1000,0,250);
   TH1F *L234_case2 = new TH1F("L234_case2","L234 with Bsl3-l3l4",1000,0,250);

   L123_case1->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L123_case1->GetXaxis()->CenterTitle(true);
   L123_case1->GetXaxis()->SetTitleSize(0.05);
   L123_case1->GetYaxis()->SetTitle("Number of tracks");
   L123_case1->GetYaxis()->CenterTitle(true);
   L123_case1->GetYaxis()->SetTitleSize(0.05);
   
   L123_case2->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L123_case2->GetXaxis()->CenterTitle(true);
   L123_case2->GetXaxis()->SetTitleSize(0.05);
   L123_case2->GetYaxis()->SetTitle("Number of tracks");
   L123_case2->GetYaxis()->CenterTitle(true);
   L123_case2->GetYaxis()->SetTitleSize(0.05);
   
   L124_case1->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L124_case1->GetXaxis()->CenterTitle(true);
   L124_case1->GetXaxis()->SetTitleSize(0.05);
   L124_case1->GetYaxis()->SetTitle("Number of tracks");
   L124_case1->GetYaxis()->CenterTitle(true);
   L124_case1->GetYaxis()->SetTitleSize(0.05);
   
   L124_case2->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L124_case2->GetXaxis()->CenterTitle(true);
   L124_case2->GetXaxis()->SetTitleSize(0.05);
   L124_case2->GetYaxis()->SetTitle("Number of tracks");
   L124_case2->GetYaxis()->CenterTitle(true);
   L124_case2->GetYaxis()->SetTitleSize(0.05);
   
   L134_case1->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L134_case1->GetXaxis()->CenterTitle(true);
   L134_case1->GetXaxis()->SetTitleSize(0.05);
   L134_case1->GetYaxis()->SetTitle("Number of tracks");
   L134_case1->GetYaxis()->CenterTitle(true);
   L134_case1->GetYaxis()->SetTitleSize(0.05);
   
   L134_case2->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L134_case2->GetXaxis()->CenterTitle(true);
   L134_case2->GetXaxis()->SetTitleSize(0.05);
   L134_case2->GetYaxis()->SetTitle("Number of tracks");
   L134_case2->GetYaxis()->CenterTitle(true);
   L134_case2->GetYaxis()->SetTitleSize(0.05);
   
   L234_case1->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L234_case1->GetXaxis()->CenterTitle(true);
   L234_case1->GetXaxis()->SetTitleSize(0.05);
   L234_case1->GetYaxis()->SetTitle("Number of tracks");
   L234_case1->GetYaxis()->CenterTitle(true);
   L234_case1->GetYaxis()->SetTitleSize(0.05);
   
   L234_case2->GetXaxis()->SetTitle("Reconstructed P_{T} (GeV)");
   L234_case2->GetXaxis()->CenterTitle(true);
   L234_case2->GetXaxis()->SetTitleSize(0.05);
   L234_case2->GetYaxis()->SetTitle("Number of tracks");
   L234_case2->GetYaxis()->CenterTitle(true);
   L234_case2->GetYaxis()->SetTitleSize(0.05);

   ofstream fout1; fout1.open("L123.log");
   ofstream fout2; fout2.open("L124.log");
   ofstream fout3; fout3.open("L134.log");
   ofstream fout4; fout4.open("L234.log");
  
  //nentries = 100;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
   //for (Long64_t jentry=0; jentry<10;jentry++) { //nentries
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
      ntEgEt_iso.clear();
      ntEgEta_iso.clear();
      ntEgPhi_iso.clear();
      ntclusterIsEG.clear();
      ntCl_match.clear();
      ntCl_iso_match.clear();
      only_iso_match.clear();
      withoutEM_match.clear();
      withEM_match.clear();
      pass_Ele.clear();
      pass_Pos.clear();
      pass_ElePos.clear();

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
      
      recon_pT.clear();

      all_cut_pass_eg = 0;
      event_nominator = 0;
      event_denominator = 0;

      int flag123;
      int flag124;
      int flag134;
      int flag234;

      float recoPV;
   
      // cout flow
      float EgEtCut  = 0;
      for(int k=0; k<EgN; k++) {
        EgEt =egEt ->at(k);
        if(EgEt > 10. ) EgEtCut = 1;
      }
      if(EgEtCut == 0) continue;
      FillCutFlow("MinEtCut", 1.);

      int EtaCutFlow = 0;
      for(int k=0; k<EgN; k++) {
        EgEta=egEta->at(k);
        EgEt =egEt ->at(k);
        //if(fabs(EgEta) < 1.3 && EgEt > 10) EtaCutFlow = 1; // for first η region
        //if(fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 && EgEt > 10) EtaCutFlow = 1;
        //if(fabs(EgEta) < 2.5 && EgEt > 9.5) EtaCutFlow = 1;
        if(fabs(EgEta) < 1.3 && EgEt > 10. ) EtaCutFlow = 1;
      }
      if(EtaCutFlow == 0) continue;
      FillCutFlow("EtaCut", 1.);

      //save L1TkElectron
      L1TkEleN = l1tkegN;
      for(int i = 0; i < l1tkegN; i++){
         ntL1TkEleEt.push_back(l1tkegEt->at(i));
	 ntL1TkEleEta.push_back(l1tkegEta->at(i));
         ntL1TkElePhi.push_back(l1tkegPhi->at(i));
      }

      //save L1TkElectron_Iso
      L1TkEleIsoN = l1tkegIsoN;
      for(int i = 0; i < l1tkegIsoN; i++){
         ntL1TkEleIsoEt.push_back(l1tkegIsoEt->at(i));
	 ntL1TkEleIsoEta.push_back(l1tkegIsoEta->at(i));
         ntL1TkEleIsoPhi.push_back(l1tkegIsoPhi->at(i));
      }


     // find egamma objects passing pixtrk signal windows
     for( int q=0; q<EgN; q++){ 

	 flag123 = 0;
	 flag124 = 0;
	 flag134 = 0;
	 flag234 = 0;
	 recoPV = 0;

	 float zp14;
	 float zp13;
	 float zp24;

      EgEt =egEt ->at(q);
      EgEta=egEta->at(q);
      EgPhi=egPhi->at(q);

      float EgGx = egGx->at(q);
      float EgGy = egGy->at(q);
      float EgGz = egGz->at(q);
      emvector.SetXYZ(EgGx,EgGy,EgGz);


      if(EgEt < 10.) continue;


      eta_region = 0; // initialize variable 
      if( fabs(EgEta) < 1.3 ) eta_region =1;
      //if( fabs(EgEta) < 1.6 && fabs(EgEta) > 1.3 ) eta_region =2;
      //if( fabs(EgEta) < 1.9 && fabs(EgEta) > 1.6 ) eta_region =3;
      //if( fabs(EgEta) < 2.5 && fabs(EgEta) > 1.9 ) eta_region =4;
      //if( fabs(EgEta) < 2.8 && fabs(EgEta) > 2.5 ) eta_region =5;
      //if( eta_region != 4 ) continue;
      //if( fabs(EgEta) > 2.5 ) continue;
      if( fabs(EgEta) > 1.3 ) continue;


      pass_egobjects_check = 1;
      ntnEg2++;
      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);

      //ntclusterIsEG.push_back(clusterIsEG->at(q));
      //ntclusterIsEG.push_back(seedIsEG->at(q));
     

      // set regin of interest
      // for η < 1.3, Δφ < 0.05 
      if(eta_region==1)SetROI(2); 
      else SetROI(eta_region);
      

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
       if(eta_region==1)SetSingalBoundary(2);
       else SetSingalBoundary(eta_region);
      

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
                                    if(L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L13_EM_Pos && L23_EM_Pos) withEM_pass_Pos = 1;
                                 }

                                if(all_cut_pass_Ele || all_cut_pass_Pos  )
                                {
                                    flag123++;

                                    float r1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    float z1 = first_layer_hits[k].Z();
                                    float r3 = sqrt(pow(third_layer_hits[j].X(),2)+pow(third_layer_hits[j].Y(),2));
                                    float z3 = third_layer_hits[j].Z();

                                    zp13 = (z1*r3-r1*z3)/(r3-r1);
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
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L12_EM_Pos && L14_EM_Pos && L24_EM_Pos) withEM_pass_Pos = 1;
                                }

                                if(all_cut_pass_Ele || all_cut_pass_Pos){
                                    flag124++;

                                    float r1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    float z1 = first_layer_hits[k].Z();
                                    float r4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    float z4 = fourth_layer_hits[j].Z();

                                    zp14 = (z1*r4-r1*z4)/(r4-r1);
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
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                }

                                if( (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L13_EM_Pos && L14_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                 }

                                if(all_cut_pass_Ele || all_cut_pass_Pos){
                                    flag124++;

                                    float r1 = sqrt(pow(first_layer_hits[k].X(),2)+pow(first_layer_hits[k].Y(),2));
                                    float z1 = first_layer_hits[k].Z();
                                    float r4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    float z4 = fourth_layer_hits[j].Z();

                                    zp14 = (z1*r4-r1*z4)/(r4-r1);
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
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                } 

                                if( (second_layer_hits_Ele_or_Pos[k] == 2 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) all_cut_pass_Pos = 1; 
                                    if(L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos) withoutEM_pass_Pos = 1;
                                    if(L23_EM_Pos && L24_EM_Pos && L34_EM_Pos) withEM_pass_Pos = 1;
                                }


                                if(all_cut_pass_Ele || all_cut_pass_Pos){
                                    flag234++;

                                    float r2 = sqrt(pow(second_layer_hits[k].X(),2)+pow(second_layer_hits[k].Y(),2));
                                    float z2 = second_layer_hits[k].Z();
                                    float r4 = sqrt(pow(fourth_layer_hits[j].X(),2)+pow(fourth_layer_hits[j].Y(),2));
                                    float z4 = fourth_layer_hits[j].Z();

                                    zp24 = (z2*r4-r2*z4)/(r4-r2);
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



      if( flag123 == 0 && flag124 == 0 && flag134 == 0 && flag234 == 0 )
      {
          ntCl_iso_match.push_back(false);
          continue;
      }
      if( flag124 || flag134 ){ recoPV = zp14; }
      if( !flag124 && !flag134 && flag123 ){ recoPV = zp13; }
      if( !flag124 && !flag134 && !flag123 && flag234 ){ recoPV = zp24; }


      /************ find recoPV is finished *****************/


      ntnEg3++;
      ntEgEt_iso.push_back(EgEt);
      ntEgEta_iso.push_back(EgEta);
      ntEgPhi_iso.push_back(EgPhi);

      cout << "Let's start track isolation" << endl;
      cout << "Pass L123: " << flag123 << ", pass L124: " << flag124 << ", pass L134: " << flag134 << ", pass L234: " << flag234 << endl;

      vector<TVector3> L1;
      vector<TVector3> L2;
      vector<TVector3> L3;
      vector<TVector3> L4;
      
      L1.clear();
      L2.clear();
      L3.clear();
      L4.clear();
      
      for( int i = 0; i < bHitN; i++){
         if(bHitGz->at(i) > 28.5) continue;
         float R = sqrt(pow(bHitGx->at(i),2)+pow(bHitGy->at(i),2));
         TVector3 bHit;	     
         bHit.SetXYZ(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)-recoPV);
         float phi = bHit.Phi();
         float eta = bHit.Eta();
         float deta = EgEta - eta;
         float dphi = EgPhi - phi;
         if( dphi < -TMath::Pi() ) dphi = dphi + 2*TMath::Pi();
         if( dphi >= TMath::Pi() ) dphi = dphi - 2*TMath::Pi();

         if( fabs(dphi) > 0.3 ) continue;

         float dR = sqrt(pow(dphi,2)+pow(deta,2));
         
         if( R < 5.5 && dR < 0.3 )            L1.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
         if( R > 5.5 && R < 8.5 && dR < 0.3 ) L2.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
         if( R > 8.5 && R < 13 && dR < 0.3 )  L3.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));
         if( R > 13 && R < 18 && dR < 0.3 )   L4.push_back(TVector3(bHitGx->at(i),bHitGy->at(i),bHitGz->at(i)));

      }//b loop

      vector<track> L123;
      L123.clear();
      vector<track> L124;
      L124.clear();
      vector<track> L134;
      L134.clear();
      vector<track> L234;
      L234.clear();
      // Above vector is declared for comparing x-positions between 3 out of 4 combinations //


      /////////////////////// 1 2 3 ///////////////////////
      // {{{
      
      Float_t dz_L13L12 = 0.048*5/3.;
      Float_t dz_L13L13 = 0.037*5/3.;
      Float_t dz_L13L23 = 0.060*5/3.;
      
      Float_t dEta_L123cut1 = 0.00486*5/3.; // L12 - L23
      Float_t dEta_L123cut2 = 0.00312*5/3.; // L12 - L13
      Float_t dEta_L123cut3 = 0.00300*5/3.; // L13 - L23
      
      Float_t dEta_L123PVcut1 = 0.00162*5/3.; // PVL3 - PVL2
      Float_t dEta_L123PVcut2 = 0.00312*5/3.; // PVL3 - PVL1

      Float_t ddphi_L123left = -0.05*5/3.;
      Float_t ddphi_L123right = 0.08*5/3.;

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
				      if( cut > ddphi_L123left && cut < ddphi_L123right )
				      {
					  L123.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 1));
				      }
		  }
	      }
	  }
      }
      
      // }}}
      
      /////////////////////// 1 2 4 ///////////////////////
           
      // {{{ 
      
      Float_t dz_L14L12 = 0.047*5/3.;
      Float_t dz_L14L14 = 0.023*5/3.;
      Float_t dz_L14L24 = 0.031*5/3.;
      
      Float_t dEta_L124cut1 = 0.00448*5/3.; // L12 - L24
      Float_t dEta_L124cut2 = 0.00324*5/3.; // L12 - L14
      Float_t dEta_L124cut3 = 0.00138*5/3.; // L14 - L24
      
      Float_t dEta_L124PVcut1 = 0.00142*5/3.; // PVL4 - PVL2
      Float_t dEta_L124PVcut2 = 0.00300*5/3.; // PVL4 - PVL1
      
      Float_t ddphi_L124left = -0.04*5/3.;
      Float_t ddphi_L124right = 0.10*5/3.;

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
				      if( cut > ddphi_L124left && cut < ddphi_L124right )
				      {
					  L124.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 2));
				      }
		  }
	      }
	  }
      }
      // }}}
      
      /////////////////////// 1 3 4 ///////////////////////
       
      // {{{ 
      
      Float_t dz_L14L13 = 0.030*5/3.;
      Float_t dz_L14L34 = 0.056*5/3.;
      
      Float_t dEta_L134cut1 = 0.00404*5/3.; // L13 - L34
      Float_t dEta_L134cut2 = 0.00152*5/3.; // L13 - L14
      Float_t dEta_L134cut3 = 0.00252*5/3.; // L14 - L34
      
      Float_t dEta_L134PVcut1 = 0.00100*5/3.; // PVL4 - PVL3
      Float_t dEta_L134PVcut2 = 0.00300*5/3.; // PVL4 - PVL1
      
      Float_t ddphi_L134left = -0.04*5/3.;
      Float_t ddphi_L134right = 0.08*5/3.;
      
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
				      if( cut > ddphi_L134left && cut < ddphi_L134right )
				      {
					  L134.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 3));
				      }
		  }
	      }
	  }
      }
      // }}}

      /////////////////////// 2 3 4 ///////////////////////
      
      // {{{
      
      Float_t dz_L24L23 = 0.060*5/3.;
      Float_t dz_L24L24 = 0.038*5/3.;
      Float_t dz_L24L34 = 0.058*5/3.;
      
      Float_t dEta_L234cut1 = 0.00417*5/3.; // L23 - L34
      Float_t dEta_L234cut2 = 0.00296*5/3.; // L23 - L24
      Float_t dEta_L234cut3 = 0.00238*5/3.; // L24 - L34
      
      Float_t dEta_L234PVcut1 = 0.00100*5/3.; // PVL4 - PVL3
      Float_t dEta_L234PVcut2 = 0.00142*5/3.; // PVL4 - PVL2
      
      Float_t ddphi_L234left = -0.06*5/3.;
      Float_t ddphi_L234right = 0.06*5/3.;

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
				      if( cut > ddphi_L234left && cut < ddphi_L234right )
				      {
					  L234.push_back(track((*c).X(), (*b).X(), (*a).X(), (*c).Y(), (*b).Y(), (*a).Y(), (*c).Z(), (*b).Z(), (*a).Z(), 4));
				      }
		  }
	      }
	  }
      }
      // }}} 
      
      /////////////////////// 3 out of 4 loop is end ///////////////////////

      /*
      cout << "L123" << endl;
      cout << "Size: " << L123.size() << endl;
      for(Int_t i = 0; i < L123.size(); i++)
      {
	  cout << "3rd X: " << L123[i].pos_x3 << ", Y: " << L123[i].pos_y3 << ", Z: " << L123[i].pos_z3 << endl;
	  cout << "2nd X: " << L123[i].pos_x2 << ", Y: " << L123[i].pos_y2 << ", Z: " << L123[i].pos_z2 << endl;
	  cout << "1st X: " << L123[i].pos_x1 << ", Y: " << L123[i].pos_y1 << ", Z: " << L123[i].pos_z1 << endl;
          cout << endl;
      }
      
      cout << "L124" << endl;
      cout << "Size: " << L124.size() << endl;
      for(Int_t i = 0; i < L124.size(); i++)
      {
	  cout << "3rd X: " << L124[i].pos_x3 << ", Y: " << L124[i].pos_y3 << ", Z: " << L124[i].pos_z3 << endl;
	  cout << "2nd X: " << L124[i].pos_x2 << ", Y: " << L124[i].pos_y2 << ", Z: " << L124[i].pos_z2 << endl;
	  cout << "1st X: " << L124[i].pos_x1 << ", Y: " << L124[i].pos_y1 << ", Z: " << L124[i].pos_z1 << endl;
          cout << endl;
      }
      
      cout << "L134" << endl;
      cout << "Size: " << L134.size() << endl;
      for(Int_t i = 0; i < L134.size(); i++)
      {
	  cout << "3rd X: " << L134[i].pos_x3 << ", Y: " << L134[i].pos_y3 << ", Z: " << L134[i].pos_z3 << endl;
	  cout << "2nd X: " << L134[i].pos_x2 << ", Y: " << L134[i].pos_y2 << ", Z: " << L134[i].pos_z2 << endl;
	  cout << "1st X: " << L134[i].pos_x1 << ", Y: " << L134[i].pos_y1 << ", Z: " << L134[i].pos_z1 << endl;
          cout << endl;
      }
      
      cout << "L234" << endl;
      cout << "Size: " << L234.size() << endl;
      for(Int_t i = 0; i < L234.size(); i++)
      {
	  cout << "3rd X: " << L234[i].pos_x3 << ", Y: " << L234[i].pos_y3 << ", Z: " << L234[i].pos_z3 << endl;
	  cout << "2nd X: " << L234[i].pos_x2 << ", Y: " << L234[i].pos_y2 << ", Z: " << L234[i].pos_z2 << endl;
	  cout << "1st X: " << L234[i].pos_x1 << ", Y: " << L234[i].pos_y1 << ", Z: " << L234[i].pos_z1 << endl;
          cout << endl;
      }
      */


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
      
      /*
      cout << "After erase" << endl << endl;
      cout << "L123" << endl;
      cout << "Size: " << L123.size() << endl;
      for(Int_t i = 0; i < L123.size(); i++)
      {
	  cout << "3rd X: " << L123[i].pos_x3 << ", Y: " << L123[i].pos_y3 << ", Z: " << L123[i].pos_z3 << endl;
	  cout << "2nd X: " << L123[i].pos_x2 << ", Y: " << L123[i].pos_y2 << ", Z: " << L123[i].pos_z2 << endl;
	  cout << "1st X: " << L123[i].pos_x1 << ", Y: " << L123[i].pos_y1 << ", Z: " << L123[i].pos_z1 << endl;
          cout << endl;
      }
      
      cout << "L124" << endl;
      cout << "Size: " << L124.size() << endl;
      for(Int_t i = 0; i < L124.size(); i++)
      {
	  cout << "3rd X: " << L124[i].pos_x3 << ", Y: " << L124[i].pos_y3 << ", Z: " << L124[i].pos_z3 << endl;
	  cout << "2nd X: " << L124[i].pos_x2 << ", Y: " << L124[i].pos_y2 << ", Z: " << L124[i].pos_z2 << endl;
	  cout << "1st X: " << L124[i].pos_x1 << ", Y: " << L124[i].pos_y1 << ", Z: " << L124[i].pos_z1 << endl;
          cout << endl;
      }
      
      cout << "L134" << endl;
      cout << "Size: " << L134.size() << endl;
      for(Int_t i = 0; i < L134.size(); i++)
      {
	  cout << "3rd X: " << L134[i].pos_x3 << ", Y: " << L134[i].pos_y3 << ", Z: " << L134[i].pos_z3 << endl;
	  cout << "2nd X: " << L134[i].pos_x2 << ", Y: " << L134[i].pos_y2 << ", Z: " << L134[i].pos_z2 << endl;
	  cout << "1st X: " << L134[i].pos_x1 << ", Y: " << L134[i].pos_y1 << ", Z: " << L134[i].pos_z1 << endl;
          cout << endl;
      }
      
      cout << "L234" << endl;
      cout << "Size: " << L234.size() << endl;
      for(Int_t i = 0; i < L234.size(); i++)
      {
	  cout << "3rd X: " << L234[i].pos_x3 << ", Y: " << L234[i].pos_y3 << ", Z: " << L234[i].pos_z3 << endl;
	  cout << "2nd X: " << L234[i].pos_x2 << ", Y: " << L234[i].pos_y2 << ", Z: " << L234[i].pos_z2 << endl;
	  cout << "1st X: " << L234[i].pos_x1 << ", Y: " << L234[i].pos_y1 << ", Z: " << L234[i].pos_z1 << endl;
          cout << endl;
      }
      */

      vector<track> all;
      all.clear();

      for(vector<track>::iterator a1 = L123.begin(); a1 != L123.end(); ++a1) all.push_back(track((*a1).pos_x3, (*a1).pos_x2, (*a1).pos_x1, (*a1).pos_y3, (*a1).pos_y2, (*a1).pos_y1, (*a1).pos_z3, (*a1).pos_z2, (*a1).pos_z1, 1 ));
      for(vector<track>::iterator a2 = L124.begin(); a2 != L124.end(); ++a2) all.push_back(track((*a2).pos_x3, (*a2).pos_x2, (*a2).pos_x1, (*a2).pos_y3, (*a2).pos_y2, (*a2).pos_y1, (*a2).pos_z3, (*a2).pos_z2, (*a2).pos_z1, 2 ));
      for(vector<track>::iterator a3 = L134.begin(); a3 != L134.end(); ++a3) all.push_back(track((*a3).pos_x3, (*a3).pos_x2, (*a3).pos_x1, (*a3).pos_y3, (*a3).pos_y2, (*a3).pos_y1, (*a3).pos_z3, (*a3).pos_z2, (*a3).pos_z1, 3 ));
      for(vector<track>::iterator a4 = L234.begin(); a4 != L234.end(); ++a4) all.push_back(track((*a4).pos_x3, (*a4).pos_x2, (*a4).pos_x1, (*a4).pos_y3, (*a4).pos_y2, (*a4).pos_y1, (*a4).pos_z3, (*a4).pos_z2, (*a4).pos_z1, 4 ));

      sort(all.begin(), all.end(), track::comp3);
      all.erase(unique(all.begin(), all.end(), track::uni3),all.end());
      sort(all.begin(), all.end(), track::comp2);
      all.erase(unique(all.begin(), all.end(), track::uni2),all.end());
      sort(all.begin(), all.end(), track::comp1);
      all.erase(unique(all.begin(), all.end(), track::uni1),all.end());

      cout << "Number of tracks: " << all.size() << endl;
      cout << "============================" << endl << endl;
      
      /*
      for(Int_t i = 0; i < all.size(); i++)
      {
	  cout << "Which combination: " << all[i].index << endl;
	  cout << "3rd X: " << all[i].pos_x3 << ", Y: " << all[i].pos_y3 << ", Z: " << all[i].pos_z3 << endl;
	  cout << "2nd X: " << all[i].pos_x2 << ", Y: " << all[i].pos_y2 << ", Z: " << all[i].pos_z2 << endl;
	  cout << "1st X: " << all[i].pos_x1 << ", Y: " << all[i].pos_y1 << ", Z: " << all[i].pos_z1 << endl;
          cout << endl;

	  Float_t R1 = sqrt(pow(all[i].pos_x1,2)+pow(all[i].pos_y1,2));
	  Float_t Z1 = all[i].pos_z1;

	  Float_t R2 = sqrt(pow(all[i].pos_x2,2)+pow(all[i].pos_y2,2));
	  Float_t Z2 = all[i].pos_z2;
		      
	  Float_t R3 = sqrt(pow(all[i].pos_x3,2)+pow(all[i].pos_y3,2));
	  Float_t Z3 = all[i].pos_z3;
	      
	  Float_t pz12 = ( R2*Z1 - R1*Z2 ) / ( R2 - R1 );
	  Float_t dz12 = pz12 - recoPV;
	  
	  Float_t pz13 = ( R3*Z1 - R1*Z3 ) / ( R3 - R1 );
	  Float_t dz13 = pz13 - recoPV;
	  
	  Float_t pz23 = ( R3*Z2 - R2*Z3 ) / ( R3 - R2 );
	  Float_t dz23 = pz23 - recoPV;
	  
	  TVector3 seg1;
	  seg1.SetXYZ( all[i].pos_x3 - all[i].pos_x2, all[i].pos_y3 - all[i].pos_y2, all[i].pos_z3 - all[i].pos_z2);
	  Float_t phi1 = seg1.Phi();
	  Float_t eta1 = seg1.PseudoRapidity();
		      
	  TVector3 seg2;
	  seg2.SetXYZ( all[i].pos_x2 - all[i].pos_x1, all[i].pos_y2 - all[i].pos_y1, all[i].pos_z2 - all[i].pos_z1);
	  Float_t phi2 = seg2.Phi();
	  Float_t eta2 = seg2.PseudoRapidity();
	  
	  TVector3 seg3;
	  seg3.SetXYZ( all[i].pos_x1, all[i].pos_y1, all[i].pos_z1 - recoPV );
	  Float_t phi3 = seg3.Phi();
	  Float_t eta3 = seg3.PseudoRapidity();
	  
	  TVector3 seg4;
	  seg4.SetXYZ( all[i].pos_x3 - all[i].pos_x1, all[i].pos_y3 - all[i].pos_y1, all[i].pos_z3 - all[i].pos_z1);
	  Float_t eta4 = seg4.PseudoRapidity();
	  
	  TVector3 seg5;
	  seg5.SetXYZ( all[i].pos_x2, all[i].pos_y2, all[i].pos_z2 - recoPV );
	  Float_t eta5 = seg5.PseudoRapidity();
	  
	  TVector3 seg6;
	  seg6.SetXYZ( all[i].pos_x3, all[i].pos_y3, all[i].pos_z3 - recoPV );
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
	  
	  if( all[i].index == 1 ) 
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
	  
	  if( all[i].index == 2 ) 
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
	  
	  if( all[i].index == 3 ) 
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
	  
	  if( all[i].index == 4 ) 
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


      } */

      if( all.size() <= 2 ) ntCl_iso_match.push_back(true);
      else ntCl_iso_match.push_back(false);

      if( all.size() <= 2 ) only_iso_match.push_back(true);
      else only_iso_match.push_back(false);

      h1->Fill( all.size() );
      if( EgEt > 20. ) h2->Fill( all.size() );
      
      if( all.size() <= 1 ) recon_pT.push_back(0.01);

      if( all.size() >= 2 )
      {
          for(vector<track>::iterator cur = all.begin(); cur != all.end(); ++cur)
          {
              Float_t recopT_L123 = -1.;
              Float_t recopT_L124 = -1.;
              Float_t recopT_L134 = -1.;
              Float_t recopT_L234 = -1.;

              if( (*cur).index == 1 ) // L123 case
              {
                  // Bsl1 - l1l3
                  Float_t p1[5] = {};
                  p1[0] = 0.0944605;
                  p1[1] = 0.00390068;
                  p1[2] = -1.02422;
                  p1[3] = -1.87021e-05;
                  p1[4] = 3.59429;

                  TVector3 pixel1;
                  pixel1.SetXYZ( (*cur).pos_x1, (*cur).pos_y1, (*cur).pos_z1 - recoPV );

                  TVector3 pixel2;
                  pixel2.SetXYZ( (*cur).pos_x3 - (*cur).pos_x1, (*cur).pos_y3 - (*cur).pos_y1, (*cur).pos_z3 - (*cur).pos_z1 );

                  Float_t phi1 = pixel1.Phi();
                  Float_t phi2 = pixel2.Phi();
                  Float_t deltaPhi1 = phi1 - phi2;
                  if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                  if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                  Float_t x1 = fabs(deltaPhi1);
                  Float_t reco_pT1 = p1[0]*pow(x1,0) + p1[1]*pow(x1,p1[2])*exp(-pow(x1,p1[3])+p1[4]);

                  L123_case1->Fill(reco_pT1);
                  
                  recopT_L123 = reco_pT1;

                  // Bsl2 - l2l3
                  Float_t p2[5] = {};
                  p2[0] = 0.183807;
                  p2[1] = 0.00390363;
                  p2[2] = -1.03385;
                  p2[3] = -1.82917e-05;
                  p2[4] = 3.5336;

                  TVector3 pixel3;
                  pixel3.SetXYZ( (*cur).pos_x2, (*cur).pos_y2, (*cur).pos_z2 - recoPV );

                  TVector3 pixel4;
                  pixel4.SetXYZ( (*cur).pos_x3 - (*cur).pos_x2, (*cur).pos_y3 - (*cur).pos_y2, (*cur).pos_z3 - (*cur).pos_z2 );

                  Float_t phi3 = pixel3.Phi();
                  Float_t phi4 = pixel4.Phi();
                  Float_t deltaPhi2 = phi3 - phi4;
                  if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                  if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                  Float_t x2 = fabs(deltaPhi2);
                  Float_t reco_pT2 = p2[0]*pow(x2,0) + p2[1]*pow(x2,p2[2])*exp(-pow(x2,p2[3])+p2[4]);

                  L123_case2->Fill(reco_pT2);

                  fout1 << "Process: " << jentry << endl;
                  fout1 << "   Number of tracks: " << all.size() << endl;
                  fout1 << "     Gen pT: " << propgenPartPt->at(0) << endl;
                  fout1 << "        Reco pT 1st case: " << reco_pT1 << ", deltaPhi1: " << fabs(deltaPhi1) << endl;
                  fout1 << "        Reco pT 2nd case: " << reco_pT2 << ", deltaPhi2: " << fabs(deltaPhi2)<< endl;
                  fout1 << "=========================================================" << endl << endl;
              }

              if( (*cur).index == 2 ) // L124 case
              {
                  // Bsl1 - l1l4
                  Float_t p1[5] = {};
                  p1[0] = 0.181996;
                  p1[1] = 0.00313492;
                  p1[2] = -1.02863;
                  p1[3] = -4.31965e-05;
                  p1[4] = 4.17543;

                  TVector3 pixel1;
                  pixel1.SetXYZ( (*cur).pos_x1, (*cur).pos_y1, (*cur).pos_z1 - recoPV );

                  TVector3 pixel2;
                  pixel2.SetXYZ( (*cur).pos_x3 - (*cur).pos_x1, (*cur).pos_y3 - (*cur).pos_y1, (*cur).pos_z3 - (*cur).pos_z1 );

                  Float_t phi1 = pixel1.Phi();
                  Float_t phi2 = pixel2.Phi();
                  Float_t deltaPhi1 = phi1 - phi2;
                  if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                  if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                  Float_t x1 = fabs(deltaPhi1);
                  Float_t reco_pT1 = p1[0]*pow(x1,0) + p1[1]*pow(x1,p1[2])*exp(-pow(x1,p1[3])+p1[4]);

                  L124_case1->Fill(reco_pT1);

                  // Bsl2 - l2l4
                  Float_t p2[5] = {};
                  p2[0] = 0.239752;
                  p2[1] = 0.00344332;
                  p2[2] = -1.26978;
                  p2[3] = -0.114929;
                  p2[4] = 3.65731;

                  TVector3 pixel3;
                  pixel3.SetXYZ( (*cur).pos_x2, (*cur).pos_y2, (*cur).pos_z2 - recoPV );

                  TVector3 pixel4;
                  pixel4.SetXYZ( (*cur).pos_x3 - (*cur).pos_x2, (*cur).pos_y3 - (*cur).pos_y2, (*cur).pos_z3 - (*cur).pos_z2 );

                  Float_t phi3 = pixel3.Phi();
                  Float_t phi4 = pixel4.Phi();
                  Float_t deltaPhi2 = phi3 - phi4;
                  if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                  if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                  Float_t x2 = fabs(deltaPhi2);
                  Float_t reco_pT2 = p2[0]*pow(x2,0) + p2[1]*pow(x2,p2[2])*exp(-pow(x2,p2[3])+p2[4]);

                  L124_case2->Fill(reco_pT2);
                  
                  recopT_L124 = reco_pT2;

                  fout2 << "Process: " << jentry << endl;
                  fout2 << "   Number of tracks: " << all.size() << endl;
                  fout2 << "     Gen pT: " << propgenPartPt->at(0) << endl;
                  fout2 << "        Reco pT 1st case: " << reco_pT1 << ", deltaPhi1: " << fabs(deltaPhi1) << endl;
                  fout2 << "        Reco pT 2nd case: " << reco_pT2 << ", deltaPhi2: " << fabs(deltaPhi2)<< endl;
                  fout2 << "=========================================================" << endl << endl;
              }

              if( (*cur).index == 3 ) // L134 case
              {
                  // Bsl1 - l1l4
                  Float_t p1[5] = {};
                  p1[0] = 0.181996;
                  p1[1] = 0.00313492;
                  p1[2] = -1.02863;
                  p1[3] = -4.31965e-05;
                  p1[4] = 4.17543;

                  TVector3 pixel1;
                  pixel1.SetXYZ( (*cur).pos_x1, (*cur).pos_y1, (*cur).pos_z1 - recoPV );

                  TVector3 pixel2;
                  pixel2.SetXYZ( (*cur).pos_x3 - (*cur).pos_x1, (*cur).pos_y3 - (*cur).pos_y1, (*cur).pos_z3 - (*cur).pos_z1 );

                  Float_t phi1 = pixel1.Phi();
                  Float_t phi2 = pixel2.Phi();
                  Float_t deltaPhi1 = phi1 - phi2;
                  if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                  if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                  Float_t x1 = fabs(deltaPhi1);
                  Float_t reco_pT1 = p1[0]*pow(x1,0) + p1[1]*pow(x1,p1[2])*exp(-pow(x1,p1[3])+p1[4]);

                  L134_case1->Fill(reco_pT1);

                  // Bsl3 - l3l4
                  Float_t p2[5] = {};
                  p2[0] = 0.808079;
                  p2[1] = 0.00343442;
                  p2[2] = -1.60152;
                  p2[3] = -0.177588;
                  p2[4] = 2.53451;

                  TVector3 pixel3;
                  pixel3.SetXYZ( (*cur).pos_x2, (*cur).pos_y2, (*cur).pos_z2 - recoPV );

                  TVector3 pixel4;
                  pixel4.SetXYZ( (*cur).pos_x3 - (*cur).pos_x2, (*cur).pos_y3 - (*cur).pos_y2, (*cur).pos_z3 - (*cur).pos_z2 );

                  Float_t phi3 = pixel3.Phi();
                  Float_t phi4 = pixel4.Phi();
                  Float_t deltaPhi2 = phi3 - phi4;
                  if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                  if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                  Float_t x2 = fabs(deltaPhi2);
                  Float_t reco_pT2 = p2[0]*pow(x2,0) + p2[1]*pow(x2,p2[2])*exp(-pow(x2,p2[3])+p2[4]);

                  L134_case2->Fill(reco_pT2);
                  
                  recopT_L134 = reco_pT2;

                  fout3 << "Process: " << jentry << endl;
                  fout3 << "   Number of tracks: " << all.size() << endl;
                  fout3 << "     Gen pT: " << propgenPartPt->at(0) << endl;
                  fout3 << "        Reco pT 1st case: " << reco_pT1 << ", deltaPhi1: " << fabs(deltaPhi1) << endl;
                  fout3 << "        Reco pT 2nd case: " << reco_pT2 << ", deltaPhi2: " << fabs(deltaPhi2)<< endl;
                  fout3 << "=========================================================" << endl << endl;
              }

              if( (*cur).index == 4 ) // L234 case
              {
                  // Bsl2 - l2l4
                  Float_t p1[5] = {};
                  p1[0] = 0.239752;
                  p1[1] = 0.00344332;
                  p1[2] = -1.26978;
                  p1[3] = -0.114929;
                  p1[4] = 3.65731;

                  TVector3 pixel1;
                  pixel1.SetXYZ( (*cur).pos_x1, (*cur).pos_y1, (*cur).pos_z1 - recoPV );

                  TVector3 pixel2;
                  pixel2.SetXYZ( (*cur).pos_x3 - (*cur).pos_x1, (*cur).pos_y3 - (*cur).pos_y1, (*cur).pos_z3 - (*cur).pos_z1 );

                  Float_t phi1 = pixel1.Phi();
                  Float_t phi2 = pixel2.Phi();
                  Float_t deltaPhi1 = phi1 - phi2;
                  if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                  if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                  Float_t x1 = fabs(deltaPhi1);
                  Float_t reco_pT1 = p1[0]*pow(x1,0) + p1[1]*pow(x1,p1[2])*exp(-pow(x1,p1[3])+p1[4]);

                  L234_case1->Fill(reco_pT1);
                  
                  recopT_L234 = reco_pT1;

                  // Bsl3 - l3l4
                  Float_t p2[5] = {};
                  p2[0] = 0.808079;
                  p2[1] = 0.00343442;
                  p2[2] = -1.60152;
                  p2[3] = -0.177588;
                  p2[4] = 2.53451;

                  TVector3 pixel3;
                  pixel3.SetXYZ( (*cur).pos_x2, (*cur).pos_y2, (*cur).pos_z2 - recoPV );

                  TVector3 pixel4;
                  pixel4.SetXYZ( (*cur).pos_x3 - (*cur).pos_x2, (*cur).pos_y3 - (*cur).pos_y2, (*cur).pos_z3 - (*cur).pos_z2 );

                  Float_t phi3 = pixel3.Phi();
                  Float_t phi4 = pixel4.Phi();
                  Float_t deltaPhi2 = phi3 - phi4;
                  if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                  if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                  Float_t x2 = fabs(deltaPhi2);
                  Float_t reco_pT2 = p2[0]*pow(x2,0) + p2[1]*pow(x2,p2[2])*exp(-pow(x2,p2[3])+p2[4]);

                  L234_case2->Fill(reco_pT2);

                  fout3 << "Process: " << jentry << endl;
                  fout3 << "   Number of tracks: " << all.size() << endl;
                  fout3 << "     Gen pT: " << propgenPartPt->at(0) << endl;
                  fout3 << "        Reco pT 1st case: " << reco_pT1 << ", deltaPhi1: " << fabs(deltaPhi1) << endl;
                  fout3 << "        Reco pT 2nd case: " << reco_pT2 << ", deltaPhi2: " << fabs(deltaPhi2)<< endl;
                  fout3 << "=========================================================" << endl << endl;
              }
              
              if( recopT_L123 != -1. ) recon_pT.push_back(recopT_L123);
              if( recopT_L124 != -1. ) recon_pT.push_back(recopT_L124);
              if( recopT_L134 != -1. ) recon_pT.push_back(recopT_L134);
              if( recopT_L234 != -1. ) recon_pT.push_back(recopT_L234);
          }
      }

  } // end of egamma loop    
  if(pass_egobjects_check){ event_denominator = 1; FillCutFlow("EvtCut", 1.);}
  if(all_cut_pass_eg) event_nominator = 1; 
     pixtrk_tree->Fill();
 } // end of entries loop 
   
   fout1.close();
   fout2.close();
   fout3.close();
   fout4.close();
   
   
   file3->Write();
}

