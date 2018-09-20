#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Float_t std_genVz;

   TString file_name;
   
   Bool_t eta_flag1 = false;
   Bool_t eta_flag2 = false;
   Bool_t eta_flag3 = false;
   Bool_t eta_flag4 = false;
   Bool_t eta_flag5 = false;
   Bool_t eta_flag6 = false;
   Bool_t eta_flag7 = false;

   TH1F *a1 = new TH1F("FiCombi_case1","; #Delta z; Number of combinations",400,-0.2,0.2);
   TH1F *a2 = new TH1F("FiCombi_case2","; #Delta z; Number of combinations",400,-0.2,0.2);
   TH1F *a3 = new TH1F("FiCombi_case3","; #Delta z; Number of combinations",400,-0.2,0.2);
   
   TH1F *b1 = new TH1F("SeCombi_case1","; #Delta z; Number of combinations",400,-0.2,0.2);
   TH1F *b2 = new TH1F("SeCombi_case2","; #Delta z; Number of combinations",400,-0.2,0.2);
   TH1F *b3 = new TH1F("SeCombi_case3","; #Delta z; Number of combinations",400,-0.2,0.2);
   
   TH1F *c1 = new TH1F("ThCombi_case1","; #Delta z; Number of combinations",400,-0.2,0.2);
   TH1F *c2 = new TH1F("ThCombi_case2","; #Delta z; Number of combinations",400,-0.2,0.2);
   TH1F *c3 = new TH1F("ThCombi_case3","; #Delta z; Number of combinations",400,-0.2,0.2);
   
   TH1F *d1 = new TH1F("FoCombi_case1","; #Delta z; Number of combinations",400,-0.2,0.2);
   TH1F *d2 = new TH1F("FoCombi_case2","; #Delta z; Number of combinations",400,-0.2,0.2);

   Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      Float_t genPt       = genPartPt->at(0);
      Float_t propGenPt   = propgenElPartPt->at(0);
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

      SaveFi_case1.clear();  
      SaveFi_case2.clear();  
      SaveFi_case3.clear();  

      SaveSe_case1.clear();  
      SaveSe_case2.clear();  
      SaveSe_case3.clear();  

      SaveTh_case1.clear();  
      SaveTh_case2.clear();  
      SaveTh_case3.clear();  

      SaveFo_case1.clear();  
      SaveFo_case2.clear();  

      if( (jentry+1) %2 != 0 ) std_genVz = -999.;

      if( genPt < 10. ) continue;
      
      Int_t eta_region = 0;
      //if( fabs(propGenEta) < 0.8  ) eta_region = 1;                              // L1234
      //if( fabs(propGenEta) > 0.8 && fabs(propGenEta) < 1.15 ) eta_region = 2;    // L1234 
      //if( fabs(propGenEta) > 1.15 && fabs(propGenEta) < 1.4 ) eta_region = 3;    // L123D1
      //if( fabs(propGenEta) > 1.4 && fabs(propGenEta) < 1.7 ) eta_region = 4;     // L12D12
      //if( fabs(propGenEta) > 1.7 && fabs(propGenEta) < 2.25 ) eta_region = 5;    // L1D123
      //if( fabs(propGenEta) > 2.25 && fabs(propGenEta) < 2.7 ) eta_region = 6;    // D2345
      if( fabs(propGenEta) > 2.7 && fabs(propGenEta) < 3.0 ) eta_region = 7;     // D3456
      if( fabs(propGenEta) > 3.0 ) continue;

      if( eta_region == 1 || eta_region == 2 ) if( bRecHitN == 0 ) continue;
      if( eta_region == 3 || eta_region == 4 || eta_region == 5 ) if( bRecHitN == 0 || fRecHitN == 0 ) continue;
      if( eta_region > 5 ) if( fRecHitN == 0 ) continue;
      
      cout << "Process: " << jentry+1 << endl;
      
      if( eta_region == 1 ) { file_name = "dZEtaRegion1.root"; eta_flag1 = true; }   
      if( eta_region == 2 ) { file_name = "dZEtaRegion2.root"; eta_flag2 = true; }
      if( eta_region == 3 ) { file_name = "dZEtaRegion3.root"; eta_flag3 = true; }
      if( eta_region == 4 ) { file_name = "dZEtaRegion4.root"; eta_flag4 = true; }
      if( eta_region == 5 ) { file_name = "dZEtaRegion5.root"; eta_flag5 = true; }
      if( eta_region == 6 ) { file_name = "dZEtaRegion6.root"; eta_flag6 = true; }
      if( eta_region == 7 ) { file_name = "dZEtaRegion7.root"; eta_flag7 = true; }
      
      Float_t bR[4] = {}, bZ[4] = {};
      Float_t fR[6] = {}, fZ[6] = {};

      // Store pixel clusters w.r.t. eta range
      StorePixelHits(eta_region, propGenPhi, propGenEta); 

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

      if( bLayer1.size() != 0 ) { bR[0] = sqrt(pow(bLayer1[0].pos_x,2)+pow(bLayer1[0].pos_y,2)); bZ[0] = bLayer1[0].pos_z; }
      if( bLayer2.size() != 0 ) { bR[1] = sqrt(pow(bLayer2[0].pos_x,2)+pow(bLayer2[0].pos_y,2)); bZ[1] = bLayer2[0].pos_z; }
      if( bLayer3.size() != 0 ) { bR[2] = sqrt(pow(bLayer3[0].pos_x,2)+pow(bLayer3[0].pos_y,2)); bZ[2] = bLayer3[0].pos_z; }
      if( bLayer4.size() != 0 ) { bR[3] = sqrt(pow(bLayer4[0].pos_x,2)+pow(bLayer4[0].pos_y,2)); bZ[3] = bLayer4[0].pos_z; }
     
      if( fDisk1.size() != 0 ) { fR[0] = sqrt(pow(fDisk1[0].pos_x,2)+pow(fDisk1[0].pos_y,2)); fZ[0] = fDisk1[0].pos_z; }
      if( fDisk2.size() != 0 ) { fR[1] = sqrt(pow(fDisk2[0].pos_x,2)+pow(fDisk2[0].pos_y,2)); fZ[1] = fDisk2[0].pos_z; }
      if( fDisk3.size() != 0 ) { fR[2] = sqrt(pow(fDisk3[0].pos_x,2)+pow(fDisk3[0].pos_y,2)); fZ[2] = fDisk3[0].pos_z; }
      if( fDisk4.size() != 0 ) { fR[3] = sqrt(pow(fDisk4[0].pos_x,2)+pow(fDisk4[0].pos_y,2)); fZ[3] = fDisk4[0].pos_z; }
      if( fDisk5.size() != 0 ) { fR[4] = sqrt(pow(fDisk5[0].pos_x,2)+pow(fDisk5[0].pos_y,2)); fZ[4] = fDisk5[0].pos_z; }
      if( fDisk6.size() != 0 ) { fR[5] = sqrt(pow(fDisk6[0].pos_x,2)+pow(fDisk6[0].pos_y,2)); fZ[5] = fDisk6[0].pos_z; }

      if( (jentry+1) %2 != 0 )
      { 
          if( eta_region == 1 || eta_region == 2 ) // 0 ~ 0.8, 0.8 ~ 1.15
          {
              if( bR[0] && bZ[0] && bR[2] && bZ[2] ) std_Vz[0] = (bR[2]*bZ[0]-bR[0]*bZ[2])/(bR[2]-bR[0]); else std_Vz[0] = 0.; // L13
              if( bR[1] && bZ[1] && bR[3] && bZ[3] ) std_Vz[1] = (bR[3]*bZ[1]-bR[1]*bZ[3])/(bR[3]-bR[1]); else std_Vz[1] = 0.; // L24
              if( bR[0] && bZ[0] && bR[3] && bZ[3] ) std_Vz[2] = (bR[3]*bZ[0]-bR[0]*bZ[3])/(bR[3]-bR[0]); else std_Vz[2] = 0.; // L14
          }
          
          if( eta_region == 3 ) // 1.15 ~ 1.4
          {
              if( bR[0] && bZ[0] && bR[2] && bZ[2] ) std_Vz[0] = (bR[2]*bZ[0]-bR[0]*bZ[2])/(bR[2]-bR[0]); else std_Vz[0] = 0.; // L1L3
              if( bR[1] && bZ[1] && fR[0] && fZ[0] ) std_Vz[1] = (fR[0]*bZ[1]-bR[1]*fZ[0])/(fR[0]-bR[1]); else std_Vz[1] = 0.; // L2D1
              if( bR[2] && bZ[2] && fR[0] && fZ[0] ) std_Vz[2] = (fR[0]*bZ[2]-bR[2]*fZ[0])/(fR[0]-bR[2]); else std_Vz[2] = 0.; // L3D1
          }
          
          if( eta_region == 4 ) // 1.4 ~ 1.7
          {
              if( bR[1] && bZ[1] && fR[0] && fZ[0] ) std_Vz[1] = (fR[0]*bZ[1]-bR[1]*fZ[0])/(fR[0]-bR[1]); else std_Vz[0] = 0.;  // L2D1
              if( bR[0] && bZ[0] && fR[1] && fZ[1] ) std_Vz[0] = (fR[1]*bZ[0]-bR[0]*fZ[1])/(fR[1]-bR[0]); else std_Vz[1] = 0.;  // L1D2
              if( bR[1] && bZ[1] && fR[1] && fZ[1] ) std_Vz[2] = (fR[1]*bZ[1]-bR[1]*fZ[1])/(fR[1]-bR[1]); else std_Vz[2] = 0.;  // L2D2
          }
          
          if( eta_region == 5 ) // 1.7 ~ 2.25
          {
              if( bR[0] && bZ[0] && fR[1] && fZ[1] ) std_Vz[0] = (fR[1]*bZ[0]-bR[0]*fZ[1])/(fR[1]-bR[0]); else std_Vz[1] = 0.; // L1D2
              if( fR[0] && fZ[0] && fR[2] && fZ[2] ) std_Vz[1] = (fR[2]*fZ[0]-fR[0]*fZ[2])/(fR[2]-fR[0]); else std_Vz[0] = 0.; // D1D3
              if( bR[0] && bZ[0] && fR[2] && fZ[2] ) std_Vz[2] = (fR[2]*bZ[0]-bR[0]*fZ[2])/(fR[2]-bR[0]); else std_Vz[2] = 0.; // L1D3
          }
          
          if( eta_region == 6 ) // 2.25 ~ 2.7
          {
              if( fR[1] && fZ[1] && fR[3] && fZ[3] ) std_Vz[0] = (fR[3]*fZ[1]-fR[1]*fZ[3])/(fR[3]-fR[1]); else std_Vz[0] = 0.; // D2D4
              if( fR[2] && fZ[2] && fR[4] && fZ[4] ) std_Vz[1] = (fR[4]*fZ[2]-fR[2]*fZ[4])/(fR[4]-fR[2]); else std_Vz[1] = 0.; // D3D5
              if( fR[1] && fZ[1] && fR[4] && fZ[4] ) std_Vz[2] = (fR[4]*fZ[1]-fR[1]*fZ[4])/(fR[4]-fR[1]); else std_Vz[2] = 0.; // D2D5
          }

          if( eta_region == 7 ) // 2.7 ~ 3.0
          {
              if( fR[2] && fZ[2] && fR[4] && fZ[4] ) std_Vz[0] = (fR[4]*fZ[2]-fR[2]*fZ[4])/(fR[4]-fR[2]); else std_Vz[0] = 0.; // D3D5
              if( fR[3] && fZ[3] && fR[5] && fZ[5] ) std_Vz[1] = (fR[5]*fZ[3]-fR[3]*fZ[5])/(fR[5]-fR[3]); else std_Vz[1] = 0.; // D4D6
              if( fR[2] && fZ[2] && fR[5] && fZ[5] ) std_Vz[2] = (fR[5]*fZ[2]-fR[2]*fZ[5])/(fR[5]-fR[2]); else std_Vz[2] = 0.; // D3D6
          }
          
         
          for(Int_t i = 0; i < 3; i++)
          {
              //cout << "This event saves reco vertex and gen vertex" << endl;
              //cout << "   Reconstructed vertex: " << i+1 << "th " << std_Vz[i] << endl;
              //cout << "   Generator level vertex: " << simVz->at(0) << endl;
              if( std_Vz[i] != 0. ) std_genVz = simVz->at(0);
          }
          //std_genVz = simVz->at(0);
      }

      cout << " ================================= " << endl << endl;
      if( std_genVz == -999. ) continue;

      if( (jentry+1) %2 == 0 ) 
      {
          MeasureDeltaZ(eta_region, std_genVz, std_Vz); 
          
          cout << "  First combination" << endl;
          if( SaveFi_case1.size() != 0 ) for(auto it : SaveFi_case1) { cout << " Delta Z: " << it << endl; a1->Fill(it); }  
          if( SaveFi_case2.size() != 0 ) for(auto it : SaveFi_case2) { cout << "  Delta Z: " << it << endl; a2->Fill(it); }
          if( SaveFi_case3.size() != 0 ) for(auto it : SaveFi_case3) { cout << "   Delta Z: " << it << endl; a3->Fill(it); }
          
          cout << "  Second combination" << endl;
          if( SaveSe_case1.size() != 0 ) for(auto it : SaveSe_case1) { cout << " Delta Z: " << it << endl; b1->Fill(it); }   
          if( SaveSe_case2.size() != 0 ) for(auto it : SaveSe_case2) { cout << "  Delta Z: " << it << endl; b2->Fill(it); }
          if( SaveSe_case3.size() != 0 ) for(auto it : SaveSe_case3) { cout << "   Delta Z: " << it << endl; b3->Fill(it); }
          
          cout << "  Third combination" << endl;
          if( SaveTh_case1.size() != 0 ) for(auto it : SaveTh_case1) { cout << " Delta Z: " << it << endl; c1->Fill(it); }
          if( SaveTh_case2.size() != 0 ) for(auto it : SaveTh_case2) { cout << "  Delta Z: " << it << endl; c2->Fill(it); }
          if( SaveTh_case3.size() != 0 ) for(auto it : SaveTh_case3) { cout << "   Delta Z: " << it << endl; c3->Fill(it); }
          
          cout << "  Fourth combination" << endl;
          if( SaveFo_case1.size() != 0 ) for(auto it : SaveFo_case1) { cout << " Delta Z: " << it << endl; d1->Fill(it); }
          if( SaveFo_case2.size() != 0 ) for(auto it : SaveFo_case2) { cout << "  Delta Z: " << it << endl; d2->Fill(it); }
                                                                     
          for(Int_t i = 0; i < 3; i++) std_Vz[i] = 0.;
      }

      cout << endl << endl;
      

   } // event loop
 
   TFile *output = new TFile(file_name,"RECREATE");
   
   if( eta_flag1 || eta_flag2 )
   {
      a1->SetTitle("L_{1}L_{3}-L_{1}L_{2}");
      a2->SetTitle("L_{1}L_{3}-L_{1}L_{3}");
      a3->SetTitle("L_{1}L_{3}-L_{2}L_{3}");
      
      b1->SetTitle("L_{2}L_{4}-L_{2}L_{4}");
      b2->SetTitle("L_{2}L_{4}-L_{2}L_{4}");
      b3->SetTitle("L_{2}L_{4}-L_{3}L_{4}");
      
      c1->SetTitle("L_{1}L_{4}-L_{1}L_{2}");
      c2->SetTitle("L_{1}L_{4}-L_{1}L_{4}");
      c3->SetTitle("L_{1}L_{4}-L_{2}L_{4}");
      
      d1->SetTitle("L_{1}L_{4}-L_{1}L_{4}");
      d2->SetTitle("L_{1}L_{4}-L_{3}L_{4}");
   }

   if( eta_flag3 )
   {
      a1->SetTitle("L_{1}L_{3}-L_{1}L_{2}");
      a2->SetTitle("L_{1}L_{3}-L_{1}L_{3}");
      a3->SetTitle("L_{1}L_{3}-L_{2}L_{3}");
      
      b1->SetTitle("L_{2}D_{1}-L_{1}L_{2}");
      b2->SetTitle("L_{2}D_{1}-L_{1}D_{1}");
      b3->SetTitle("L_{2}D_{1}-L_{2}D_{1}");
      
      c1->SetTitle("L_{3}D_{1}-L_{1}L_{3}");
      c2->SetTitle("L_{3}D_{1}-L_{1}D_{1}");
      c3->SetTitle("L_{3}D_{1}-L_{3}D_{1}");
      
      d1->SetTitle("L_{3}D_{1}-L_{2}L_{3}");
      d2->SetTitle("L_{3}D_{1}-L_{2}D_{1}");
   }
   
   if( eta_flag4 )
   {
      a1->SetTitle("L_{1}D_{2}-L_{1}D_{1}");
      a2->SetTitle("L_{1}D_{2}-L_{1}D_{2}");
      a3->SetTitle("L_{1}D_{2}-D_{1}D_{2}");
      
      b1->SetTitle("L_{2}D_{1}-L_{1}L_{2}");
      b2->SetTitle("L_{2}D_{1}-L_{1}D_{1}");
      b3->SetTitle("L_{2}D_{1}-L_{2}D_{1}");
      
      c1->SetTitle("L_{2}D_{2}-L_{1}L_{2}");
      c2->SetTitle("L_{2}D_{2}-L_{1}D_{2}");
      c3->SetTitle("L_{2}D_{2}-L_{2}D_{2}");
      
      d1->SetTitle("L_{2}D_{2}-L_{2}D_{1}");
      d2->SetTitle("L_{2}D_{2}-D_{1}D_{2}");
   }

   if( eta_flag5 )
   {
      a1->SetTitle("L_{1}D_{2}-L_{1}D_{1}");
      a2->SetTitle("L_{1}D_{2}-L_{1}D_{2}");
      a3->SetTitle("L_{1}D_{2}-D_{1}D_{2}");
      
      b1->SetTitle("D_{1}D_{3}-D_{1}D_{2}");
      b2->SetTitle("D_{1}D_{3}-D_{1}D_{3}");
      b3->SetTitle("D_{1}D_{3}-D_{2}D_{3}");
      
      c1->SetTitle("L_{1}D_{3}-L_{1}D_{1}");
      c2->SetTitle("L_{1}D_{3}-L_{1}D_{3}");
      c3->SetTitle("L_{1}D_{3}-D_{1}D_{3}");
      
      d1->SetTitle("L_{1}D_{3}-L_{1}D_{2}");
      d2->SetTitle("L_{1}D_{3}-D_{2}D_{3}");
   }
   
   if( eta_flag6 )
   {
      a1->SetTitle("D_{2}D_{4}-D_{2}D_{3}");
      a2->SetTitle("D_{2}D_{4}-D_{2}D_{4}");
      a3->SetTitle("D_{2}D_{4}-D_{3}D_{4}");
      
      b1->SetTitle("D_{3}D_{5}-D_{3}D_{4}");
      b2->SetTitle("D_{3}D_{5}-D_{3}D_{5}");
      b3->SetTitle("D_{3}D_{5}-D_{4}D_{5}");
      
      c1->SetTitle("D_{2}D_{5}-D_{2}D_{3}");
      c2->SetTitle("D_{2}D_{5}-D_{2}D_{5}");
      c3->SetTitle("D_{2}D_{5}-D_{3}D_{5}");
      
      d1->SetTitle("D_{2}D_{5}-D_{2}D_{4}");
      d2->SetTitle("D_{2}D_{5}-D_{4}D_{5}");
   }
   
   if( eta_flag7 )
   {
      a1->SetTitle("D_{3}D_{5}-D_{3}D_{4}");
      a2->SetTitle("D_{3}D_{5}-D_{3}D_{5}");
      a3->SetTitle("D_{3}D_{5}-D_{4}D_{5}");
      
      b1->SetTitle("D_{4}D_{6}-D_{4}D_{5}");
      b2->SetTitle("D_{4}D_{6}-D_{4}D_{6}");
      b3->SetTitle("D_{4}D_{6}-D_{5}D_{6}");
      
      c1->SetTitle("D_{3}D_{6}-D_{3}D_{4}");
      c2->SetTitle("D_{3}D_{6}-D_{3}D_{6}");
      c3->SetTitle("D_{3}D_{6}-D_{4}D_{6}");
      
      d1->SetTitle("D_{3}D_{6}-D_{3}D_{5}");
      d2->SetTitle("D_{3}D_{6}-D_{5}D_{6}");
   }
   
   a1->GetXaxis()->CenterTitle(true);
   a1->GetXaxis()->SetTitleSize(0.05);
   a1->GetYaxis()->SetTitleSize(0.04);
   a2->GetXaxis()->CenterTitle(true);
   a2->GetXaxis()->SetTitleSize(0.05);
   a2->GetYaxis()->SetTitleSize(0.04);
   a3->GetXaxis()->CenterTitle(true);
   a3->GetXaxis()->SetTitleSize(0.05);
   a3->GetYaxis()->SetTitleSize(0.04);
   
   b1->GetXaxis()->CenterTitle(true);
   b1->GetXaxis()->SetTitleSize(0.05);
   b1->GetYaxis()->SetTitleSize(0.04);
   b2->GetXaxis()->CenterTitle(true);
   b2->GetXaxis()->SetTitleSize(0.05);
   b2->GetYaxis()->SetTitleSize(0.04);
   b3->GetXaxis()->CenterTitle(true);
   b3->GetXaxis()->SetTitleSize(0.05);
   b3->GetYaxis()->SetTitleSize(0.04);
   
   c1->GetXaxis()->CenterTitle(true);
   c1->GetXaxis()->SetTitleSize(0.05);
   c1->GetYaxis()->SetTitleSize(0.04);
   c2->GetXaxis()->CenterTitle(true);
   c2->GetXaxis()->SetTitleSize(0.05);
   c2->GetYaxis()->SetTitleSize(0.04);
   c3->GetXaxis()->CenterTitle(true);
   c3->GetXaxis()->SetTitleSize(0.05);
   c3->GetYaxis()->SetTitleSize(0.04);
   
   d1->GetXaxis()->CenterTitle(true);
   d1->GetXaxis()->SetTitleSize(0.05);
   d1->GetYaxis()->SetTitleSize(0.04);
   d2->GetXaxis()->CenterTitle(true);
   d2->GetXaxis()->SetTitleSize(0.05);
   d2->GetYaxis()->SetTitleSize(0.04);
   
   a1->Write();
   a2->Write();
   a3->Write();
   b1->Write();
   b2->Write();
   b3->Write();
   c1->Write();
   c2->Write();
   c3->Write();
   d1->Write();
   d2->Write();
   
   output->Close();

}
