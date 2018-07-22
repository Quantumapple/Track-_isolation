#define para_cxx
#include "para.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void para::Loop()
{
//   In a ROOT session, you can do:
//      root> .L para.C
//      root> para t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *output = new TFile("paras.root","RECREATE");
   
   // Histograms for DeltaEta(pixel, pixel)
   TH1F *h1 = new TH1F("dEtaL123pp1","L_{1}L_{2}-L_{1}L_{3}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *h2 = new TH1F("dEtaL123pp2","L_{1}L_{2}-L_{2}L_{3}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *h3 = new TH1F("dEtaL123pp3","L_{1}L_{3}-L_{2}L_{3}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   
   TH1F *a1 = new TH1F("dEtaL124pp1","L_{1}L_{2}-L_{1}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *a2 = new TH1F("dEtaL124pp2","L_{1}L_{2}-L_{2}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *a3 = new TH1F("dEtaL124pp3","L_{1}L_{4}-L_{2}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   
   TH1F *b1 = new TH1F("dEtaL134pp1","L_{1}L_{3}-L_{1}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *b2 = new TH1F("dEtaL134pp2","L_{1}L_{3}-L_{3}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *b3 = new TH1F("dEtaL134pp3","L_{1}L_{4}-L_{3}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);

   TH1F *c1 = new TH1F("dEtaL234pp1","L_{2}L_{3}-L_{2}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *c2 = new TH1F("dEtaL234pp2","L_{2}L_{3}-L_{3}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *c3 = new TH1F("dEtaL234pp3","L_{2}L_{3}-L_{3}L_{4}; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   
   // Histograms for DeltaEta(pv, pixel)
   TH1F *hh1 = new TH1F("dEtaL123PV1","PVL3-PVL2; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *hh2 = new TH1F("dEtaL123PV2","PVL3-PVL1; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   
   TH1F *aa1 = new TH1F("dEtaL124PV1","PVL4-PVL2; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   TH1F *aa2 = new TH1F("dEtaL124PV2","PVL4-PVL1; #Delta#eta; Number of combinations",1000,-0.01,0.01);
   
   TH1F *bb1 = new TH1F("dEtaL134PV1","PVL4-PVL3; #Delta#eta; Number of combinations",1000,-0.01,0.01);

   // Histograms for DeltaPhi difference
   TH1F *hhh1 = new TH1F("dPhiL123","#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2L3)]; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);
   TH1F *aaa1 = new TH1F("dPhiL124","#Delta#phi[#phi(PVL1)-#phi(L1L2)]-#Delta#phi[#phi(L1L2)-#phi(L2L4)]; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);
   TH1F *bbb1 = new TH1F("dPhiL134","#Delta#phi[#phi(PVL1)-#phi(L1L3)]-#Delta#phi[#phi(L1L3)-#phi(L3L4)]; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);
   TH1F *ccc1 = new TH1F("dPhiL234","#Delta#phi[#phi(PVL2)-#phi(L2L3)]-#Delta#phi[#phi(L2L3)-#phi(L3L4)]; #Delta#phi_{i}-#Delta#phi_{j}; Number of combinations",1000,-0.01,0.01);

   // Class for pixel hits 
   class track
   {
       public:
       float DR;
       float pos_x, pos_y, pos_z;
       track() { DR = 0.; pos_x = 0.; pos_y = 0.; pos_z = 0.; }
       track(float a, float b, float c, float d) { DR = a; pos_x = b; pos_y = c; pos_z = d; }
       static bool comp(const track &t1, const track &t2)
       {
	   return ( t1.DR < t2.DR );
       }
   };

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      Float_t genPt       = genPartPt->at(0);
      Float_t propGenPhi  = propgenElPartPhi->at(0);
      Float_t propGenEta  = propgenElPartEta->at(0);

      cout << "Process: " << jentry+1 << "/" << nentries << endl;

      if( genPt < 10. ) continue;

      vector<track> v1;
      vector<track> v2;
      vector<track> v3;
      vector<track> v4;

      v1.clear();
      v2.clear();
      v3.clear();
      v4.clear();

      if( bRecHitN == 0 ) continue;
      for(Int_t i = 0; i < bRecHitN; i++)
      {
          Float_t R = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
          TVector3 pixel;
          pixel.SetXYZ(bRecHitGx->at(i) - simVx->at(0), bRecHitGy->at(i) - simVy->at(0), bRecHitGz->at(i) - simVz->at(0));
          Float_t pixelPhi = pixel.Phi();
          Float_t pixelEta = pixel.PseudoRapidity();
          Float_t deltaPhi = pixelPhi - propGenPhi;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
          Float_t deltaEta = pixelEta - propGenEta;
          Float_t deltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta,2));

          if( R < 5. )             v1.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          if( R > 5. && R < 9. )   v2.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          if( R > 9. && R < 14. )  v3.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          if( R > 14.&& R < 18. )  v4.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
      }

      sort(v1.begin(), v1.end(), track::comp);
      sort(v2.begin(), v2.end(), track::comp);
      sort(v3.begin(), v3.end(), track::comp);
      sort(v4.begin(), v4.end(), track::comp);

      //// -------- L123 --------- //// 
      if( v1.size() != 0 && v2.size() != 0 && v3.size() != 0 )
      {
          for(std::vector<track>::iterator it1 = v1.begin(); it1 != v1.end(); ++it1)
          {
              TVector3 pv1; // PVL1
              pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
              Float_t pv1_eta = pv1.Eta();
              Float_t pv1_phi = pv1.Phi();
              
              for(std::vector<track>::iterator it2 = v2.begin(); it2 != v2.end(); ++it2)
              {
                  TVector3 pp1; // L12
                  pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                  Float_t pp1_phi = pp1.Phi();
                  Float_t pp1_eta = pp1.Eta();

                  TVector3 pv2; // PVL2
                  pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                  Float_t pv2_eta = pv2.Eta();

                  for(std::vector<track>::iterator it3 = v3.begin(); it3 != v3.end(); ++it3 )
                  {
                      TVector3 pp2; // L23
                      pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                      Float_t pp2_phi = pp2.Phi();
                      Float_t pp2_eta = pp2.Eta();
                      
                      TVector3 pp3; // L13
                      pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                      Float_t pp3_eta = pp3.Eta();

                      TVector3 pv3; // PVL3
                      pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                      Float_t pv3_eta = pv3.Eta();
                  
                      Float_t dEtaL1223 = pp1_eta - pp2_eta;
                      Float_t dEtaL1213 = pp1_eta - pp3_eta;
                      Float_t dEtaL1323 = pp3_eta - pp2_eta;
                      Float_t dEtaPVL13 = pv3_eta - pv1_eta;
                      Float_t dEtaPVL23 = pv3_eta - pv2_eta;
                 
                      Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                      if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                      if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                      Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2L3
                      if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                      if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                      Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                      h1->Fill(dEtaL1213); 
                      h2->Fill(dEtaL1223); 
                      h3->Fill(dEtaL1323); 
                      hh1->Fill(dEtaPVL23);
                      hh2->Fill(dEtaPVL13);
                      hhh1->Fill(ddPhi);
                  }
              }
          }
      }
      
      //// -------- L124 --------- //// 
      if( v1.size() != 0 && v2.size() != 0 && v4.size() != 0 )
      {
          for(std::vector<track>::iterator it1 = v1.begin(); it1 != v1.end(); ++it1)
          {
              TVector3 pv1; // PVL1
              pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
              Float_t pv1_eta = pv1.Eta();
              Float_t pv1_phi = pv1.Phi();
              
              for(std::vector<track>::iterator it2 = v2.begin(); it2 != v2.end(); ++it2)
              {
                  TVector3 pp1; // L12
                  pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                  Float_t pp1_phi = pp1.Phi();
                  Float_t pp1_eta = pp1.Eta();

                  TVector3 pv2; // PVL2
                  pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                  Float_t pv2_eta = pv2.Eta();

                  for(std::vector<track>::iterator it3 = v4.begin(); it3 != v4.end(); ++it3 )
                  {
                      TVector3 pp2; // L24
                      pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                      Float_t pp2_phi = pp2.Phi();
                      Float_t pp2_eta = pp2.Eta();
                      
                      TVector3 pp3; // L14
                      pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                      Float_t pp3_eta = pp3.Eta();

                      TVector3 pv3; // PVL4
                      pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                      Float_t pv3_eta = pv3.Eta();
                  
                      Float_t dEtaL1224 = pp1_eta - pp2_eta;
                      Float_t dEtaL1214 = pp1_eta - pp3_eta;
                      Float_t dEtaL1424 = pp3_eta - pp2_eta;
                      Float_t dEtaPVL14 = pv3_eta - pv1_eta;
                      Float_t dEtaPVL24 = pv3_eta - pv2_eta;
                 
                      Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                      if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                      if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                      Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2L4
                      if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                      if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                      Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                      a1->Fill(dEtaL1214); 
                      a2->Fill(dEtaL1224); 
                      a3->Fill(dEtaL1424); 
                      aa1->Fill(dEtaPVL24);
                      aa2->Fill(dEtaPVL14);
                      aaa1->Fill(ddPhi);
                  }
              }
          }
      }
      
      //// -------- L134 --------- //// 
      if( v1.size() != 0 && v3.size() != 0 && v4.size() != 0 )
      {
          for(std::vector<track>::iterator it1 = v1.begin(); it1 != v1.end(); ++it1)
          {
              TVector3 pv1; // PVL1
              pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
              Float_t pv1_phi = pv1.Phi();
              for(std::vector<track>::iterator it2 = v3.begin(); it2 != v3.end(); ++it2)
              {
                  TVector3 pp1; // L13
                  pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                  Float_t pp1_phi = pp1.Phi();
                  Float_t pp1_eta = pp1.Eta();

                  TVector3 pv2; // PVL3
                  pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                  Float_t pv2_eta = pv2.Eta();

                  for(std::vector<track>::iterator it3 = v4.begin(); it3 != v4.end(); ++it3 )
                  {
                      TVector3 pp2; // L34
                      pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                      Float_t pp2_phi = pp2.Phi();
                      Float_t pp2_eta = pp2.Eta();
                      
                      TVector3 pp3; // L14
                      pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                      Float_t pp3_eta = pp3.Eta();

                      TVector3 pv3; // PVL4
                      pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                      Float_t pv3_eta = pv3.Eta();
                  
                      Float_t dEtaL1334 = pp1_eta - pp2_eta;
                      Float_t dEtaL1314 = pp1_eta - pp3_eta;
                      Float_t dEtaL1434 = pp3_eta - pp2_eta;
                      Float_t dEtaPVL34 = pv3_eta - pv2_eta;
                 
                      Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L3
                      if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                      if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                      Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L3 - L3L4
                      if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                      if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                      Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                      b1->Fill(dEtaL1314); 
                      b2->Fill(dEtaL1334); 
                      b3->Fill(dEtaL1434); 
                      bb1->Fill(dEtaPVL34);
                      bbb1->Fill(ddPhi);
                  }
              }
          }
      }
      
      //// -------- L234 --------- //// 
      if( v2.size() != 0 && v3.size() != 0 && v4.size() != 0 )
      {
          for(std::vector<track>::iterator it1 = v2.begin(); it1 != v2.end(); ++it1)
          {
              TVector3 pv1; // PVL2
              pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
              Float_t pv1_phi = pv1.Phi();
              for(std::vector<track>::iterator it2 = v3.begin(); it2 != v3.end(); ++it2)
              {
                  TVector3 pp1; // L23
                  pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                  Float_t pp1_phi = pp1.Phi();
                  Float_t pp1_eta = pp1.Eta();

                  for(std::vector<track>::iterator it3 = v4.begin(); it3 != v4.end(); ++it3 )
                  {
                      TVector3 pp2; // L34
                      pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                      Float_t pp2_phi = pp2.Phi();
                      Float_t pp2_eta = pp2.Eta();
                      
                      TVector3 pp3; // L24
                      pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                      Float_t pp3_eta = pp3.Eta();

                      Float_t dEtaL2334 = pp1_eta - pp2_eta;
                      Float_t dEtaL2324 = pp1_eta - pp3_eta;
                      Float_t dEtaL2434 = pp3_eta - pp2_eta;
                 
                      Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL2 - L2L3
                      if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                      if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                      Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L2L3 - L3L4
                      if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                      if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                      Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                      c1->Fill(dEtaL2324); 
                      c2->Fill(dEtaL2334); 
                      c3->Fill(dEtaL2434); 
                      ccc1->Fill(ddPhi);
                  }
              }
          }
      }

   
   } // event loop

   h1->Write();
   h2->Write();
   h3->Write();
   hh1->Write();
   hh2->Write();
   hhh1->Write();
   
   a1->Write();
   a2->Write();
   a3->Write();
   aa1->Write();
   aa2->Write();
   aaa1->Write();
   
   b1->Write();
   b2->Write();
   b3->Write();
   bb1->Write();
   bbb1->Write();
   
   c1->Write();
   c2->Write();
   c3->Write();
   ccc1->Write();

   output->Close();

}
