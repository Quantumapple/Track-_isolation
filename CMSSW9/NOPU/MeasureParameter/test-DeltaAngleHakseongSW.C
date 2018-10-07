#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;

   TString file_name;
   
   // Histograms for DeltaEta(pixel, pixel)
   TH1F *h1 = new TH1F("dEtaCase1pp1","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *h2 = new TH1F("dEtaCase1pp2","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *h3 = new TH1F("dEtaCase1pp3","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   
   TH1F *a1 = new TH1F("dEtaCase2pp1","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *a2 = new TH1F("dEtaCase2pp2","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *a3 = new TH1F("dEtaCase2pp3","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   
   TH1F *b1 = new TH1F("dEtaCase3pp1","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *b2 = new TH1F("dEtaCase3pp2","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *b3 = new TH1F("dEtaCase3pp3","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);

   TH1F *c1 = new TH1F("dEtaCase4pp1","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *c2 = new TH1F("dEtaCase4pp2","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *c3 = new TH1F("dEtaCase4pp3","; #Delta#eta (pixel,pixel); Number of combinations",1000,-0.01,0.01);
   
   // Histograms for DeltaEta(pv, pixel)
   TH1F *hh1 = new TH1F("dEta1STPV1","; #Delta#eta (PV,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *hh2 = new TH1F("dEta2NDPV2","; #Delta#eta (PV,pixel); Number of combinations",1000,-0.01,0.01);
   
   TH1F *aa1 = new TH1F("dEta3RDPV1","; #Delta#eta (PV,pixel); Number of combinations",1000,-0.01,0.01);
   TH1F *aa2 = new TH1F("dEta4THPV2","; #Delta#eta (PV,pixel); Number of combinations",1000,-0.01,0.01);
   
   TH1F *bb1 = new TH1F("dEta5THPV1","; #Delta#eta (PV,pixel); Number of combinations",1000,-0.01,0.01);
   
   // Histograms for DeltaPhi difference
   TH1F *hhh1 = new TH1F("dPhi1ST","; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);
   TH1F *aaa1 = new TH1F("dPhi2ND","; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);
   TH1F *bbb1 = new TH1F("dPhi3RD","; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);
   TH1F *ccc1 = new TH1F("dPhi4TH","; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);

   Bool_t eta_flag1 = false;
   Bool_t eta_flag2 = false;
   Bool_t eta_flag3 = false;
   Bool_t eta_flag4 = false;
   Bool_t eta_flag5 = false;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<1000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      bLayer1.clear();
      bLayer2.clear();
      bLayer3.clear();
      bLayer4.clear();

      fDisk1.clear();
      fDisk2.clear();
      fDisk3.clear();
      fDisk4.clear();
      fDisk5.clear();
   
      saveParaCase1.clear();
      saveParaCase2.clear();
      saveParaCase3.clear();
      saveParaCase4.clear();
   
      Float_t genPt       = genPartPt->at(0);
      Float_t propGenPt   = propgenElPartPt->at(0);
      Float_t propGenPhi  = propgenElPartPhi->at(0);
      Float_t propGenEta  = propgenElPartEta->at(0);
      
      if( genPt < 10. ) continue;

      Int_t eta_region = 0;
      //if( fabs(propGenEta) < 0.8  ) eta_region = 1;                           // L1234
      //if( fabs(propGenEta) > 0.8 && fabs(propGenEta) < 1.4 ) eta_region = 2;  // L123D1 
      //if( fabs(propGenEta) > 1.4 && fabs(propGenEta) < 1.8 ) eta_region = 3;  // L12D12
      //if( fabs(propGenEta) > 1.8 && fabs(propGenEta) < 2.4 ) eta_region = 4;  // L1D123
      if( fabs(propGenEta) > 2.4 && fabs(propGenEta) < 3.0 ) eta_region = 5;  // D2345
      if( fabs(propGenEta) > 3.0 ) continue;

      if( eta_region == 1 ) if( bRecHitN == 0 ) continue;
      if( eta_region >= 2 && eta_region == 4 ) if( bRecHitN == 0 || fRecHitN == 0 ) continue;
      if( eta_region == 5 ) if( fRecHitN == 0 ) continue;
      
      cout << "Process: " << jentry << endl;

      if( eta_region == 1 ) { file_name = "dAngleEtaRegion1.root"; eta_flag1 = true; }   
      if( eta_region == 2 ) { file_name = "dAngleEtaRegion2.root"; eta_flag2 = true; }
      if( eta_region == 3 ) { file_name = "dAngleEtaRegion3.root"; eta_flag3 = true; }
      if( eta_region == 4 ) { file_name = "dAngleEtaRegion4.root"; eta_flag4 = true; }
      if( eta_region == 5 ) { file_name = "dAngleEtaRegion5.root"; eta_flag5 = true; }
    
      // Store pixel clusters w.r.t eta range
      StorePixelHits(eta_region, propGenPhi, propGenEta);
     
      // Sort vectors to get closest pixel clusters
      sort(bLayer1.begin(), bLayer1.end(), track::comp);
      sort(bLayer2.begin(), bLayer2.end(), track::comp);
      sort(bLayer3.begin(), bLayer3.end(), track::comp);
      sort(bLayer4.begin(), bLayer4.end(), track::comp);
      sort(fDisk1.begin(), fDisk1.end(), track::comp);
      sort(fDisk2.begin(), fDisk2.end(), track::comp);
      sort(fDisk3.begin(), fDisk3.end(), track::comp);
      sort(fDisk4.begin(), fDisk4.end(), track::comp);
      sort(fDisk5.begin(), fDisk5.end(), track::comp);

      // Calculate and save dEta(pixel,pixel), dEta(PV, pixel), dPhi difference
      CalculateParameters(eta_region);

      // Fill Histograms
      if( saveParaCase1.size() != 0 )
      {
          for(Int_t i = 0; i < saveParaCase1.size(); i++)
          {
              if( i%6 == 0 ) h1->Fill(saveParaCase1.at(i)); 
              if( i%6 == 1 ) h2->Fill(saveParaCase1.at(i)); 
              if( i%6 == 2 ) h3->Fill(saveParaCase1.at(i)); 
              if( i%6 == 3 ) hh1->Fill(saveParaCase1.at(i));
              if( i%6 == 4 ) hh2->Fill(saveParaCase1.at(i));
              if( i%6 == 5 ) hhh1->Fill(saveParaCase1.at(i));
          }
      }
      
      if( saveParaCase2.size() != 0 )
      {
          for(Int_t i = 0; i < saveParaCase2.size(); i++)
          {
              if( i%6 == 0 ) a1->Fill(saveParaCase2.at(i)); 
              if( i%6 == 1 ) a2->Fill(saveParaCase2.at(i)); 
              if( i%6 == 2 ) a3->Fill(saveParaCase2.at(i)); 
              if( i%6 == 3 ) aa1->Fill(saveParaCase2.at(i));
              if( i%6 == 4 ) aa2->Fill(saveParaCase2.at(i));
              if( i%6 == 5 ) aaa1->Fill(saveParaCase2.at(i));
          }
      }
      
      if( saveParaCase3.size() != 0 )
      {
          for(Int_t i = 0; i < saveParaCase3.size(); i++)
          {
              if( i%5 == 0) b1->Fill(saveParaCase3.at(i)); 
              if( i%5 == 1) b2->Fill(saveParaCase3.at(i)); 
              if( i%5 == 2) b3->Fill(saveParaCase3.at(i)); 
              if( i%5 == 3) bb1->Fill(saveParaCase3.at(i));
              if( i%5 == 4) bbb1->Fill(saveParaCase3.at(i));
          }
      }
      
      if( saveParaCase4.size() != 0 )
      {
          for(Int_t i = 0; i < saveParaCase4.size(); i++)
          {
              if( i%4 == 0 ) c1->Fill(saveParaCase4.at(i)); 
              if( i%4 == 1 ) c2->Fill(saveParaCase4.at(i)); 
              if( i%4 == 2 ) c3->Fill(saveParaCase4.at(i)); 
              if( i%4 == 3 ) ccc1->Fill(saveParaCase4.at(i));
          }
      }
      

   } // event loop
   
   TFile *output = new TFile(file_name,"RECREATE");
  
   if( eta_flag1 )
   {
      h1->SetTitle("L_{1}L_{2}-L_{1}L_{3}");
      h2->SetTitle("L_{1}L_{2}-L_{2}L_{3}");
      h3->SetTitle("L_{1}L_{3}-L_{2}L_{3}");
      hh1->SetTitle("PVL3-PVL2");
      hh2->SetTitle("PVL3-PVL1");
      hhh1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2L3)]");
      
      a1->SetTitle("L_{1}L_{2}-L_{1}L_{4}"); 
      a2->SetTitle("L_{1}L_{2}-L_{2}L_{4}"); 
      a3->SetTitle("L_{1}L_{4}-L_{2}L_{4}"); 
      aa1->SetTitle("PVL4-PVL2");
      aa2->SetTitle("PVL4-PVL1");
      aaa1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2L4)]");
      
      b1->SetTitle("L_{1}L_{3}-L_{1}L_{4}");
      b2->SetTitle("L_{1}L_{3}-L_{3}L_{4}");
      b3->SetTitle("L_{1}L_{4}-L_{3}L_{4}");
      bb1->SetTitle("PVL4-PVL3");
      bbb1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L3)]-#Delta#phi[#phi(L1L3)-#phi(L3L4)]");
      
      c1->SetTitle("L_{2}L_{3}-L_{2}L_{4}");
      c2->SetTitle("L_{2}L_{3}-L_{3}L_{4}");
      c3->SetTitle("L_{2}L_{4}-L_{3}L_{4}");
      ccc1->SetTitle("#Delta#phi[#phi(PVL2)-#phi(L2L3)]-#Delta#phi[#phi(L2L3)-#phi(L3L4)]");
   }

   if( eta_flag2 )
   {
      h1->SetTitle("L_{1}L_{2}-L_{1}L_{3}");
      h2->SetTitle("L_{1}L_{2}-L_{2}L_{3}");
      h3->SetTitle("L_{1}L_{3}-L_{2}L_{3}");
      hh1->SetTitle("PVL3-PVL2");
      hh2->SetTitle("PVL3-PVL1");
      hhh1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2L3)]");
      
      a1->SetTitle("L_{1}L_{2}-L_{1}D_{1}"); 
      a2->SetTitle("L_{1}L_{2}-L_{2}D_{1}"); 
      a3->SetTitle("L_{1}D_{1}-L_{2}D_{1}"); 
      aa1->SetTitle("PVD1-PVL2");
      aa2->SetTitle("PVD1-PVL1");
      aaa1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2D1)]");
      
      b1->SetTitle("L_{1}L_{3}-L_{1}D_{1}");
      b2->SetTitle("L_{1}L_{3}-L_{3}D_{1}");
      b3->SetTitle("L_{1}D_{1}-L_{3}D_{1}");
      bb1->SetTitle("PVD1-PVL3");
      bbb1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L3)]-#Delta#phi[#phi(L1L3)-#phi(L3D1)]");
      
      c1->SetTitle("L_{2}L_{3}-L_{2}D_{1}");
      c2->SetTitle("L_{2}L_{3}-L_{3}D_{1}");
      c3->SetTitle("L_{2}D_{1}-L_{3}D_{1}");
      ccc1->SetTitle("#Delta#phi[#phi(PVL2)-#phi(L2L3)]-#Delta#phi[#phi(L2L3)-#phi(L3D1)]");
   }
   
   if( eta_flag3 )
   {
      h1->SetTitle("L_{1}L_{2}-L_{1}D_{1}");
      h2->SetTitle("L_{1}L_{2}-L_{2}D_{1}");
      h3->SetTitle("L_{1}D_{1}-L_{2}D_{1}");
      hh1->SetTitle("PVD1-PVL2");
      hh2->SetTitle("PVD1-PVL1");
      hhh1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2D1)]");
      
      a1->SetTitle("L_{1}L_{2}-L_{1}D_{2}"); 
      a2->SetTitle("L_{1}L_{2}-L_{2}D_{2}"); 
      a3->SetTitle("L_{1}D_{2}-L_{2}D_{2}"); 
      aa1->SetTitle("PVD2-PVL2");
      aa2->SetTitle("PVD2-PVL1");
      aaa1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2D2)]");
      
      b1->SetTitle("L_{1}D_{1}-L_{1}D_{2}");
      b2->SetTitle("L_{1}D_{1}-D_{1}D_{2}");
      b3->SetTitle("L_{1}D_{2}-D_{1}D_{2}");
      bb1->SetTitle("PVD2-PVD1");
      bbb1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1D1)]-#Delta#phi[#phi(L1D1)-#phi(D1D2)]");
      
      c1->SetTitle("L_{2}D_{1}-L_{2}D_{2}");
      c2->SetTitle("L_{2}D_{1}-D_{1}D_{2}");
      c3->SetTitle("L_{2}D_{2}-D_{1}D_{2}");
      ccc1->SetTitle("#Delta#phi[#phi(PVL2)-#phi(L2D1)]-#Delta#phi[#phi(L2D1)-#phi(D1D2)]");
   }

   if( eta_flag4 )
   {
      h1->SetTitle("L_{1}D_{1}-L_{1}D_{2}");
      h2->SetTitle("L_{1}D_{1}-D_{1}D_{2}");
      h3->SetTitle("L_{1}D_{2}-D_{1}D_{2}");
      hh1->SetTitle("PVD2-PVD1");
      hh2->SetTitle("PVD2-PVL1");
      hhh1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1D1)]-#Delta#phi[#phi(L1D1)-#phi(D1D2)]");
      
      a1->SetTitle("L_{1}D_{1}-L_{1}D_{3}"); 
      a2->SetTitle("L_{1}D_{1}-D_{1}D_{3}"); 
      a3->SetTitle("L_{1}D_{3}-D_{1}D_{3}"); 
      aa1->SetTitle("PVD3-PVD1");
      aa2->SetTitle("PVD3-PVL1");
      aaa1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1D1)]-#Delta#phi[#phi(L1D1)-#phi(D1D3)]");
      
      b1->SetTitle("L_{1}D_{2}-L_{1}D_{3}");
      b2->SetTitle("L_{1}D_{2}-D_{2}D_{3}");
      b3->SetTitle("L_{1}D_{3}-D_{2}D_{3}");
      bb1->SetTitle("PVD3-PVD2");
      bbb1->SetTitle("#Delta#phi[#phi(PVL1)-#phi(L1D2)]-#Delta#phi[#phi(L1D2)-#phi(D2D3)]");
      
      c1->SetTitle("D_{1}D_{2}-D_{1}D_{3}");
      c2->SetTitle("D_{1}D_{2}-D_{2}D_{3}");
      c3->SetTitle("D_{1}D_{3}-D_{2}D_{3}");
      ccc1->SetTitle("#Delta#phi[#phi(PVD1)-#phi(D1D2)]-#Delta#phi[#phi(D1D2)-#phi(D2D3)]");
   }
   
   if( eta_flag5 )
   {
      h1->SetTitle("D_{2}D_{3}-D_{2}D_{4}");
      h2->SetTitle("D_{2}D_{3}-D_{3}D_{4}");
      h3->SetTitle("D_{2}D_{4}-D_{3}D_{4}");
      hh1->SetTitle("PVD4-PVD3");
      hh2->SetTitle("PVD4-PVD2");
      hhh1->SetTitle("#Delta#phi[#phi(PVD2)-#phi(D2D3)]-#Delta#phi[#phi(D2D3)-#phi(D3D4)]");
      
      a1->SetTitle("D_{2}D_{3}-D_{2}D_{5}"); 
      a2->SetTitle("D_{2}D_{3}-D_{3}D_{5}"); 
      a3->SetTitle("D_{2}D_{5}-D_{3}D_{5}"); 
      aa1->SetTitle("PVD5-PVD3");
      aa2->SetTitle("PVD5-PVD2");
      aaa1->SetTitle("#Delta#phi[#phi(PVD2)-#phi(D2D3)]-#Delta#phi[#phi(D2D3)-#phi(D3D5)]");
      
      b1->SetTitle("D_{2}D_{4}-D_{2}D_{5}");
      b2->SetTitle("D_{2}D_{4}-D_{4}D_{5}");
      b3->SetTitle("D_{2}D_{5}-D_{4}D_{5}");
      bb1->SetTitle("PVD5-PVD4");
      bbb1->SetTitle("#Delta#phi[#phi(PVD2)-#phi(D2D4)]-#Delta#phi[#phi(D2D4)-#phi(D4D5)]");
      
      c1->SetTitle("D_{3}D_{4}-D_{3}D_{5}");
      c2->SetTitle("D_{3}D_{4}-D_{4}D_{5}");
      c3->SetTitle("D_{3}D_{5}-D_{4}D_{5}");
      ccc1->SetTitle("#Delta#phi[#phi(PVD3)-#phi(D3D4)]-#Delta#phi[#phi(D3D4)-#phi(D4D5)]");
   }

   h1->GetXaxis()->CenterTitle(true);
   h1->GetXaxis()->SetTitleSize(0.05);
   h1->GetYaxis()->SetTitleSize(0.04);
   h2->GetXaxis()->CenterTitle(true);
   h2->GetXaxis()->SetTitleSize(0.05);
   h2->GetYaxis()->SetTitleSize(0.04);
   h3->GetXaxis()->CenterTitle(true);
   h3->GetXaxis()->SetTitleSize(0.05);
   h3->GetYaxis()->SetTitleSize(0.04);
   hh1->GetXaxis()->CenterTitle(true);
   hh1->GetXaxis()->SetTitleSize(0.05);
   hh1->GetYaxis()->SetTitleSize(0.04);
   hh2->GetXaxis()->CenterTitle(true);
   hh2->GetXaxis()->SetTitleSize(0.05);
   hh2->GetYaxis()->SetTitleSize(0.04);
   hhh1->GetXaxis()->CenterTitle(true);
   hhh1->GetXaxis()->SetTitleSize(0.05);
   hhh1->GetYaxis()->SetTitleSize(0.04);

   h1->Write();
   h2->Write();
   h3->Write();
   hh1->Write();
   hh2->Write();
   hhh1->Write();
   
   a1->GetXaxis()->CenterTitle(true);
   a1->GetXaxis()->SetTitleSize(0.05);
   a1->GetYaxis()->SetTitleSize(0.04);
   a2->GetXaxis()->CenterTitle(true);
   a2->GetXaxis()->SetTitleSize(0.05);
   a2->GetYaxis()->SetTitleSize(0.04);
   a3->GetXaxis()->CenterTitle(true);
   a3->GetXaxis()->SetTitleSize(0.05);
   a3->GetYaxis()->SetTitleSize(0.04);
   aa1->GetXaxis()->CenterTitle(true);
   aa1->GetXaxis()->SetTitleSize(0.05);
   aa1->GetYaxis()->SetTitleSize(0.04);
   aa2->GetXaxis()->CenterTitle(true);
   aa2->GetXaxis()->SetTitleSize(0.05);
   aa2->GetYaxis()->SetTitleSize(0.04);
   aaa1->GetXaxis()->CenterTitle(true);
   aaa1->GetXaxis()->SetTitleSize(0.05);
   aaa1->GetYaxis()->SetTitleSize(0.04);
   
   a1->Write();
   a2->Write();
   a3->Write();
   aa1->Write();
   aa2->Write();
   aaa1->Write();
   
   b1->GetXaxis()->CenterTitle(true);
   b1->GetXaxis()->SetTitleSize(0.05);
   b1->GetYaxis()->SetTitleSize(0.04);
   b2->GetXaxis()->CenterTitle(true);
   b2->GetXaxis()->SetTitleSize(0.05);
   b2->GetYaxis()->SetTitleSize(0.04);
   b3->GetXaxis()->CenterTitle(true);
   b3->GetXaxis()->SetTitleSize(0.05);
   b3->GetYaxis()->SetTitleSize(0.04);
   bb1->GetXaxis()->CenterTitle(true);
   bb1->GetXaxis()->SetTitleSize(0.05);
   bb1->GetYaxis()->SetTitleSize(0.04);
   bbb1->GetXaxis()->CenterTitle(true);
   bbb1->GetXaxis()->SetTitleSize(0.05);
   bbb1->GetYaxis()->SetTitleSize(0.04);
   
   b1->Write();
   b2->Write();
   b3->Write();
   bb1->Write();
   bbb1->Write();
   
   c1->GetXaxis()->CenterTitle(true);
   c1->GetXaxis()->SetTitleSize(0.05);
   c1->GetYaxis()->SetTitleSize(0.04);
   c2->GetXaxis()->CenterTitle(true);
   c2->GetXaxis()->SetTitleSize(0.05);
   c2->GetYaxis()->SetTitleSize(0.04);
   c3->GetXaxis()->CenterTitle(true);
   c3->GetXaxis()->SetTitleSize(0.05);
   c3->GetYaxis()->SetTitleSize(0.04);
   ccc1->GetXaxis()->CenterTitle(true);
   ccc1->GetXaxis()->SetTitleSize(0.05);
   ccc1->GetYaxis()->SetTitleSize(0.04);
   
   c1->Write();
   c2->Write();
   c3->Write();
   ccc1->Write();

   output->Close();
}
