#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void test::Loop()
{
//   In a ROOT session, you can do:
//      root> .L test.C
//      root> test t
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

   TFile *output = new TFile("pars.root","RECREATE");

   TH2F *rz = new TH2F("rz","rz plane; z; r",2000,-30,30,2000,0,20);

   TH1F *h1 = new TH1F("h1","odd L_{1}L_{4} - even L_{1}L_{2}",2000,-1,1);
   TH1F *h2 = new TH1F("h2","odd L_{1}L_{4} - even L_{1}L_{3}",2000,-1,1);
   TH1F *h3 = new TH1F("h3","odd L_{1}L_{4} - even L_{1}L_{4}",2000,-1,1);
   TH1F *h4 = new TH1F("h4","odd L_{1}L_{4} - even L_{2}L_{3}",2000,-1,1);
   TH1F *h5 = new TH1F("h5","odd L_{1}L_{4} - even L_{2}L_{4}",2000,-1,1);
   TH1F *h6 = new TH1F("h6","odd L_{1}L_{4} - even L_{3}L_{4}",2000,-1,1);
   
   TH1F *a1 = new TH1F("a1","odd L_{1}L_{3} - even L_{1}L_{2}",2000,-1,1);
   TH1F *a2 = new TH1F("a2","odd L_{1}L_{3} - even L_{1}L_{3}",2000,-1,1);
   TH1F *a3 = new TH1F("a3","odd L_{1}L_{3} - even L_{1}L_{4}",2000,-1,1);
   TH1F *a4 = new TH1F("a4","odd L_{1}L_{3} - even L_{2}L_{3}",2000,-1,1);
   TH1F *a5 = new TH1F("a5","odd L_{1}L_{3} - even L_{2}L_{4}",2000,-1,1);
   TH1F *a6 = new TH1F("a6","odd L_{1}L_{3} - even L_{3}L_{4}",2000,-1,1);
   
   TH1F *b1 = new TH1F("b1","odd L_{2}L_{4} - even L_{1}L_{2}",2000,-1,1);
   TH1F *b2 = new TH1F("b2","odd L_{2}L_{4} - even L_{1}L_{3}",2000,-1,1);
   TH1F *b3 = new TH1F("b3","odd L_{2}L_{4} - even L_{1}L_{4}",2000,-1,1);
   TH1F *b4 = new TH1F("b4","odd L_{2}L_{4} - even L_{2}L_{3}",2000,-1,1);
   TH1F *b5 = new TH1F("b5","odd L_{2}L_{4} - even L_{2}L_{4}",2000,-1,1);
   TH1F *b6 = new TH1F("b6","odd L_{2}L_{4} - even L_{3}L_{4}",2000,-1,1);
  
   Float_t odd_zp14 = -999.;
   Float_t odd_zp13 = -999.;
   Float_t odd_zp24 = -999.;
   Float_t odd_gen = -999.;

   Bool_t odd_flag1;
   Bool_t odd_flag2;
   Bool_t odd_flag3;
   Int_t event = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    
      // Drawing pixel barrels in rz plane
      for(Int_t i = 0; i < bRecHitN; i++)
      {
          Float_t R = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
          rz->Fill(bRecHitGz->at(0),R);
      }
      
      Float_t genPt       = genPartPt->at(0);
      //Float_t genPhi      = genPartPhi->at(0);
      //Float_t genEta      = genPartEta->at(0);
      Float_t propGenPt   = propgenElPartPt->at(0);
      Float_t propGenPhi  = propgenElPartPhi->at(0);
      Float_t propGenEta  = propgenElPartEta->at(0);

      if( genPt < 10. ) continue;
      
      cout << "Process: " << jentry << "/" << nentries << endl;
      event++;
      if( event % 2 != 0 ){ odd_flag1 = false; odd_flag2 = false; odd_flag3 = false;} // odd event number

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

      Float_t R1 = 0.; 
      Float_t Z1 = 0.;
      Float_t R2 = 0.;
      Float_t Z2 = 0.;
      Float_t R3 = 0.;
      Float_t Z3 = 0.;
      Float_t R4 = 0.;
      Float_t Z4 = 0.;

      if( v1.size() != 0 ) R1 = sqrt(pow(v1[0].pos_x,2)+pow(v1[0].pos_y,2));
      if( v1.size() != 0 ) Z1 = v1[0].pos_z;
      if( v2.size() != 0 ) R2 = sqrt(pow(v2[0].pos_x,2)+pow(v2[0].pos_y,2));
      if( v2.size() != 0 ) Z2 = v2[0].pos_z;
      if( v3.size() != 0 ) R3 = sqrt(pow(v3[0].pos_x,2)+pow(v3[0].pos_y,2));
      if( v3.size() != 0 ) Z3 = v3[0].pos_z;
      if( v4.size() != 0 ) R4 = sqrt(pow(v4[0].pos_x,2)+pow(v4[0].pos_y,2));
      if( v4.size() != 0 ) Z4 = v4[0].pos_z;
    
      Float_t zp12 = -999.;
      Float_t zp13 = -999.;
      Float_t zp14 = -999.;
      Float_t zp23 = -999.;
      Float_t zp24 = -999.;
      Float_t zp34 = -999.;

      if( R1 && R2 && Z1 && Z2 ) zp12 = (R2*Z1-R1*Z2)/(R2-R1);
      if( R1 && R3 && Z1 && Z3 ) zp13 = (R3*Z1-R1*Z3)/(R3-R1);
      if( R1 && R4 && Z1 && Z4 ) zp14 = (R4*Z1-R1*Z4)/(R4-R1);
      if( R2 && R3 && Z2 && Z3 ) zp23 = (R3*Z2-R2*Z3)/(R3-R2);
      if( R2 && R4 && Z2 && Z4 ) zp24 = (R4*Z2-R2*Z4)/(R4-R2);
      if( R3 && R4 && Z3 && Z4 ) zp34 = (R4*Z3-R3*Z4)/(R4-R3);

      if( event % 2 != 0 ) // odd event number
      {
          if( zp14 != -999. ) { odd_flag1 = true; odd_zp14 = zp14; cout << "Odd14: " << jentry << endl;}
          if( zp13 != -999. ) { odd_flag2 = true; odd_zp13 = zp13; cout << "Odd13: " << jentry << endl;}
          if( zp24 != -999. ) { odd_flag3 = true; odd_zp24 = zp24; cout << "Odd24: " << jentry << endl;}
          if( odd_flag1 || odd_flag2 || odd_flag3 ) { cout << "Odd_gen: " << jentry << endl; odd_gen = simVz->at(0); }
          cout << "=============================" << endl << endl;
      }


      if( event % 2 == 0 ) // even event number
      {
          if( zp12 != -999. ) 
          {
              cout << "even zp12" << endl;
              if( odd_flag1 ) // zp14
              {
                  Float_t diff_z = fabs(odd_zp14 - zp12);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  h1->Fill(diff_z-gen_diff);
              }
              if( odd_flag2 ) // zp13
              {
                  Float_t diff_z = fabs(odd_zp13 - zp12);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  a1->Fill(diff_z-gen_diff);
              }
              if( odd_flag3 ) // zp24
              {
                  Float_t diff_z = fabs(odd_zp24 - zp12);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  b1->Fill(diff_z-gen_diff);
              }
          }

          if( zp13 != -999. ) 
          {  
              cout << "even zp13" << endl;
              if( odd_flag1 ) // zp14
              {
                  Float_t diff_z = fabs(odd_zp14 - zp13);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  h2->Fill(diff_z-gen_diff);
              }
              if( odd_flag2 ) // zp13
              {
                  Float_t diff_z = fabs(odd_zp13 - zp13);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  a2->Fill(diff_z-gen_diff);
              }
              if( odd_flag3 ) // zp24
              {
                  Float_t diff_z = fabs(odd_zp24 - zp13);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  b2->Fill(diff_z-gen_diff);
              }
          }

          if( zp14 != -999. ) 
          {  
              cout << "even zp14" << endl;
              if( odd_flag1 ) // zp14
              {
                  Float_t diff_z = fabs(odd_zp14 - zp14);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  h3->Fill(diff_z-gen_diff);
              }
              if( odd_flag2 ) // zp13
              {
                  Float_t diff_z = fabs(odd_zp13 - zp14);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  a3->Fill(diff_z-gen_diff);
              }
              if( odd_flag3 ) // zp24
              {
                  Float_t diff_z = fabs(odd_zp24 - zp14);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  b3->Fill(diff_z-gen_diff);
              }
          }

          if( zp23 != -999. ) 
          {  
              cout << "even zp23" << endl;
              if( odd_flag1 ) // zp14
              {
                  Float_t diff_z = fabs(odd_zp14 - zp23);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  h4->Fill(diff_z-gen_diff);
              }
              if( odd_flag2 ) // zp13
              {
                  Float_t diff_z = fabs(odd_zp13 - zp23);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  a4->Fill(diff_z-gen_diff);
              }
              if( odd_flag3 ) // zp24
              {
                  Float_t diff_z = fabs(odd_zp24 - zp23);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  b4->Fill(diff_z-gen_diff);
              }
          }

          if( zp24 != -999. ) 
          {  
              cout << "even zp24" << endl;
              if( odd_flag1 ) // zp14
              {
                  Float_t diff_z = fabs(odd_zp14 - zp24);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  h5->Fill(diff_z-gen_diff);
              }
              if( odd_flag2 ) // zp13
              {
                  Float_t diff_z = fabs(odd_zp13 - zp24);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  a5->Fill(diff_z-gen_diff);
              }
              if( odd_flag3 ) // zp24
              {
                  Float_t diff_z = fabs(odd_zp24 - zp24);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  b5->Fill(diff_z-gen_diff);
              }
          }

          if( zp34 != -999. ) 
          {  
              cout << "even zp34" << endl;
              if( odd_flag1 ) // zp14
              {
                  Float_t diff_z = fabs(odd_zp14 - zp34);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  h6->Fill(diff_z-gen_diff);
              }
              if( odd_flag2 ) // zp13
              {
                  Float_t diff_z = fabs(odd_zp13 - zp34);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  a6->Fill(diff_z-gen_diff);
              }
              if( odd_flag3 ) // zp24
              {
                  Float_t diff_z = fabs(odd_zp24 - zp34);
                  Float_t gen_diff = fabs(odd_gen - simVz->at(0));
                  b6->Fill(diff_z-gen_diff);
              }
          }
          cout << endl;
      }


      
      /*
      // Check only for barrel region
      if( fabs(propGenEta) > 0.8 ) continue;

      if( bRecHitN == 0 ) continue;

      Float_t layer1R = -999.; Float_t layer1Z = -999.;
      Float_t layer2R = -999.; Float_t layer2Z = -999.;
      Float_t layer3R = -999.; Float_t layer3Z = -999.;
      Float_t layer4R = -999.; Float_t layer4Z = -999.;

      Float_t clodRL1 = 999.;
      Float_t clodRL2 = 999.;
      Float_t clodRL3 = 999.;
      Float_t clodRL4 = 999.;
      
      // Divide pixel clusters by its radius
      for(Int_t i = 0; i < bRecHitN; i++)
      {
          Float_t R = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
          TVector3 pixel;
          pixel.SetXYZ( bRecHitGx->at(i) - simVx->at(0), bRecHitGy->at(i) - simVy->at(0), bRecHitGz->at(i) - simVz->at(0));
          Float_t pixelPhi = pixel.Phi();
          Float_t pixelEta = pixel.Eta();
          Float_t deltaPhi = propGenPhi-pixelPhi;
          if( deltaPhi < -TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi >= TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
          Float_t deltaEta = propGenEta-pixelEta;
          Float_t deltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta,2));
          if( R < 4. )              { if( deltaR < clodRL1 ) clodRL1 = deltaR; layer1R = R; layer1Z = bRecHitGz->at(0); }
          if( R > 4. && R < 8. )    { if( deltaR < clodRL2 ) clodRL2 = deltaR; layer2R = R; layer2Z = bRecHitGz->at(0); }
          if( R > 8. && R < 13.)    { if( deltaR < clodRL3 ) clodRL3 = deltaR; layer3R = R; layer3Z = bRecHitGz->at(0); }
          if( R > 13.&& R < 17. )   { if( deltaR < clodRL4 ) clodRL4 = deltaR; layer4R = R; layer4Z = bRecHitGz->at(0); }
      }

      if( layer1R != -999. && layer2R != -999.) dz12->Fill( ( layer2R*layer1Z - layer1R*layer2Z ) / ( layer2R - layer1R ) - simVz->at(0) );
      if( layer1R != -999. && layer3R != -999.) dz13->Fill( ( layer3R*layer1Z - layer1R*layer3Z ) / ( layer3R - layer1R ) - simVz->at(0) );
      if( layer1R != -999. && layer4R != -999.) dz14->Fill( ( layer4R*layer1Z - layer1R*layer4Z ) / ( layer4R - layer1R ) - simVz->at(0) );
      if( layer2R != -999. && layer3R != -999.) dz23->Fill( ( layer3R*layer2Z - layer2R*layer3Z ) / ( layer3R - layer2R ) - simVz->at(0) );
      if( layer2R != -999. && layer4R != -999.) dz24->Fill( ( layer4R*layer2Z - layer2R*layer4Z ) / ( layer4R - layer2R ) - simVz->at(0) );
      if( layer3R != -999. && layer4R != -999.) dz34->Fill( ( layer4R*layer3Z - layer3R*layer4Z ) / ( layer4R - layer3R ) - simVz->at(0) );
      */
     
      
   


   } // event loop

   h1->SetXTitle("#Deltaz[cm]");
   h1->GetXaxis()->CenterTitle(true);
   h1->GetXaxis()->SetTitleSize(0.05);
   h1->Write();
   h2->SetXTitle("#Deltaz[cm]");
   h2->GetXaxis()->CenterTitle(true);
   h2->GetXaxis()->SetTitleSize(0.05);
   h2->Write();
   h3->SetXTitle("#Deltaz[cm]");
   h3->GetXaxis()->CenterTitle(true);
   h3->GetXaxis()->SetTitleSize(0.05);
   h3->Write();
   h4->SetXTitle("#Deltaz[cm]");
   h4->GetXaxis()->CenterTitle(true);
   h4->GetXaxis()->SetTitleSize(0.05);
   h4->Write();
   h5->SetXTitle("#Deltaz[cm]");
   h5->GetXaxis()->CenterTitle(true);
   h5->GetXaxis()->SetTitleSize(0.05);
   h5->Write();
   h6->SetXTitle("#Deltaz[cm]");
   h6->GetXaxis()->CenterTitle(true);
   h6->GetXaxis()->SetTitleSize(0.05);
   h6->Write();

   a1->SetXTitle("#Deltaz[cm]");
   a1->GetXaxis()->CenterTitle(true);
   a1->GetXaxis()->SetTitleSize(0.05);
   a1->Write();
   a2->SetXTitle("#Deltaz[cm]");
   a2->GetXaxis()->CenterTitle(true);
   a2->GetXaxis()->SetTitleSize(0.05);
   a2->Write();
   a3->SetXTitle("#Deltaz[cm]");
   a3->GetXaxis()->CenterTitle(true);
   a3->GetXaxis()->SetTitleSize(0.05);
   a3->Write();
   a4->SetXTitle("#Deltaz[cm]");
   a4->GetXaxis()->CenterTitle(true);
   a4->GetXaxis()->SetTitleSize(0.05);
   a4->Write();
   a5->SetXTitle("#Deltaz[cm]");
   a5->GetXaxis()->CenterTitle(true);
   a5->GetXaxis()->SetTitleSize(0.05);
   a5->Write();
   a6->SetXTitle("#Deltaz[cm]");
   a6->GetXaxis()->CenterTitle(true);
   a6->GetXaxis()->SetTitleSize(0.05);
   a6->Write();

   b1->SetXTitle("#Deltaz[cm]");
   b1->GetXaxis()->CenterTitle(true);
   b1->GetXaxis()->SetTitleSize(0.05);
   b1->Write();
   b2->SetXTitle("#Deltaz[cm]");
   b2->GetXaxis()->CenterTitle(true);
   b2->GetXaxis()->SetTitleSize(0.05);
   b2->Write();
   b3->SetXTitle("#Deltaz[cm]");
   b3->GetXaxis()->CenterTitle(true);
   b3->GetXaxis()->SetTitleSize(0.05);
   b3->Write();
   b4->SetXTitle("#Deltaz[cm]");
   b4->GetXaxis()->CenterTitle(true);
   b4->GetXaxis()->SetTitleSize(0.05);
   b4->Write();
   b5->SetXTitle("#Deltaz[cm]");
   b5->GetXaxis()->CenterTitle(true);
   b5->GetXaxis()->SetTitleSize(0.05);
   b5->Write();
   b6->SetXTitle("#Deltaz[cm]");
   b6->GetXaxis()->CenterTitle(true);
   b6->GetXaxis()->SetTitleSize(0.05);
   b6->Write();
   output->Close();
}
