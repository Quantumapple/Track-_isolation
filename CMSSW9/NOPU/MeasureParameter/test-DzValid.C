#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

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

   TFile *output = new TFile("validation.root","RECREATE");
   TH1F *a1 = new TH1F("a1","L1-L3; #Delta z; Number of combinations",2000,-1,1);
   
   TH1F *h1 = new TH1F("h1","L1-D1; #Delta z; Number of combinations",2000,-1,1);
   TH1F *h2 = new TH1F("h2","L2-D1; #Delta z; Number of combinations",2000,-1,1);
   TH1F *h3 = new TH1F("h3","L3-D1; #Delta z; Number of combinations",2000,-1,1);
   
   TH1F *h4 = new TH1F("h4","L1-D2; #Delta z; Number of combinations",2000,-1,1);
   TH1F *h5 = new TH1F("h5","L2-D2; #Delta z; Number of combinations",2000,-1,1);

   TH1F *h6 = new TH1F("h6","L1-D3; #Delta z; Number of combinations",2000,-1,1);
   
   TH1F *disk1A = new TH1F("disk1A","D1-D2; #Delta z; Number of combinations",2000,-1,1);
   TH1F *disk1B = new TH1F("disk1B","D1-D3; #Delta z; Number of combinations",2000,-1,1);
   TH1F *disk1C = new TH1F("disk1C","D1-D4; #Delta z; Number of combinations",2000,-1,1);
   
   TH1F *disk2A = new TH1F("disk2A","D2-D3; #Delta z; Number of combinations",2000,-1,1);
   TH1F *disk2B = new TH1F("disk2B","D2-D4; #Delta z; Number of combinations",2000,-1,1);
   TH1F *disk2C = new TH1F("disk2C","D2-D5; #Delta z; Number of combinations",2000,-1,1);
   
   TH1F *disk3A = new TH1F("disk3A","D3-D4; #Delta z; Number of combinations",2000,-1,1);
   TH1F *disk3B = new TH1F("disk3B","D3-D5; #Delta z; Number of combinations",2000,-1,1);
   TH1F *disk3C = new TH1F("disk3C","D3-D6; #Delta z; Number of combinations",2000,-1,1);
   
   TH1F *disk4A = new TH1F("disk4A","D4-D5; #Delta z; Number of combinations",2000,-1,1);
   TH1F *disk4B = new TH1F("disk4B","D4-D6; #Delta z; Number of combinations",2000,-1,1);
   
   TH1F *disk5A = new TH1F("disk5A","D5-D6; #Delta z; Number of combinations",2000,-1,1);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<100;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      Float_t genPt = genPartPt->at(0); 
      if( genPt < 10. ) continue;
      
      cout << "Process: " << jentry << endl;
      
      if( bRecHitN != 0 && fRecHitN != 0 ) 
      {
          for(Int_t i = 0; i < bRecHitN; i++)
          {
              for(Int_t j = 0; j < fRecHitN; j++)
              {
                  if( bRecHitLayer->at(i) == 1 && fRecHitDisk->at(j) == 1 )
                  {
                      Float_t bR = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
                      Float_t bZ = bRecHitGz->at(i);
                      Float_t fR = sqrt(pow(fRecHitGx->at(j),2)+pow(fRecHitGy->at(j),2));
                      Float_t fZ = fRecHitGz->at(j);
                      Float_t recoZ = (fR*bZ - bR*fZ) / (fR - bR);
                      Float_t deltaZ = simVz->at(0) - recoZ;
                      h1->Fill(deltaZ);
                  }

                  if( bRecHitLayer->at(i) == 2 && fRecHitDisk->at(j) == 1 )
                  {
                      Float_t bR = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
                      Float_t bZ = bRecHitGz->at(i);
                      Float_t fR = sqrt(pow(fRecHitGx->at(j),2)+pow(fRecHitGy->at(j),2));
                      Float_t fZ = fRecHitGz->at(j);
                      Float_t recoZ = (fR*bZ - bR*fZ) / (fR - bR);
                      Float_t deltaZ = simVz->at(0) - recoZ;
                      h2->Fill(deltaZ);
                  }

                  if( bRecHitLayer->at(i) == 3 && fRecHitDisk->at(j) == 1 )
                  {
                      Float_t bR = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
                      Float_t bZ = bRecHitGz->at(i);
                      Float_t fR = sqrt(pow(fRecHitGx->at(j),2)+pow(fRecHitGy->at(j),2));
                      Float_t fZ = fRecHitGz->at(j);
                      Float_t recoZ = (fR*bZ - bR*fZ) / (fR - bR);
                      Float_t deltaZ = simVz->at(0) - recoZ;
                      h3->Fill(deltaZ);
                  }
                  
                  if( bRecHitLayer->at(i) == 1 && fRecHitDisk->at(j) == 2 )
                  {
                      Float_t bR = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
                      Float_t bZ = bRecHitGz->at(i);
                      Float_t fR = sqrt(pow(fRecHitGx->at(j),2)+pow(fRecHitGy->at(j),2));
                      Float_t fZ = fRecHitGz->at(j);
                      Float_t recoZ = (fR*bZ - bR*fZ) / (fR - bR);
                      Float_t deltaZ = simVz->at(0) - recoZ;
                      h4->Fill(deltaZ);
                  }

                  if( bRecHitLayer->at(i) == 2 && fRecHitDisk->at(j) == 2 )
                  {
                      Float_t bR = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
                      Float_t bZ = bRecHitGz->at(i);
                      Float_t fR = sqrt(pow(fRecHitGx->at(j),2)+pow(fRecHitGy->at(j),2));
                      Float_t fZ = fRecHitGz->at(j);
                      Float_t recoZ = (fR*bZ - bR*fZ) / (fR - bR);
                      Float_t deltaZ = simVz->at(0) - recoZ;
                      h5->Fill(deltaZ);
                  }
                  
                  if( bRecHitLayer->at(i) == 1 && fRecHitDisk->at(j) == 3 )
                  {
                      Float_t bR = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));
                      Float_t bZ = bRecHitGz->at(i);
                      Float_t fR = sqrt(pow(fRecHitGx->at(j),2)+pow(fRecHitGy->at(j),2));
                      Float_t fZ = fRecHitGz->at(j);
                      Float_t recoZ = (fR*bZ - bR*fZ) / (fR - bR);
                      Float_t deltaZ = simVz->at(0) - recoZ;
                      h6->Fill(deltaZ);
                  }
              }
          }
      }

      vector<TVector3> fDisk1;
      vector<TVector3> fDisk2;
      vector<TVector3> fDisk3;
      vector<TVector3> fDisk4;
      vector<TVector3> fDisk5;
      vector<TVector3> fDisk6;

      fDisk1.clear();
      fDisk2.clear();
      fDisk3.clear();
      fDisk4.clear();
      fDisk5.clear();
      fDisk6.clear();

      if( fRecHitN != 0 )
      {
          for(Int_t j = 0; j < fRecHitN; j++)
          {
              if( fRecHitDisk->at(j) == 1 )
              {
                  fDisk1.push_back(TVector3(fRecHitGx->at(j), fRecHitGy->at(j), fRecHitGz->at(j)));
              }
              if( fRecHitDisk->at(j) == 2 )
              {
                  fDisk2.push_back(TVector3(fRecHitGx->at(j), fRecHitGy->at(j), fRecHitGz->at(j)));
              }
              if( fRecHitDisk->at(j) == 3 )
              {
                  fDisk3.push_back(TVector3(fRecHitGx->at(j), fRecHitGy->at(j), fRecHitGz->at(j)));
              }
              if( fRecHitDisk->at(j) == 4 )
              {
                  fDisk4.push_back(TVector3(fRecHitGx->at(j), fRecHitGy->at(j), fRecHitGz->at(j)));
              }
              if( fRecHitDisk->at(j) == 5 )
              {
                  fDisk5.push_back(TVector3(fRecHitGx->at(j), fRecHitGy->at(j), fRecHitGz->at(j)));
              }
              if( fRecHitDisk->at(j) == 6 )
              {
                  fDisk6.push_back(TVector3(fRecHitGx->at(j), fRecHitGy->at(j), fRecHitGz->at(j)));
              }
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk1.begin(); it1 != fDisk1.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk2.begin(); it2 != fDisk2.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk1A->Fill(deltaZ);
          }
      }

      for(vector<TVector3>::iterator it1 = fDisk1.begin(); it1 != fDisk1.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk3.begin(); it2 != fDisk3.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk1B->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk1.begin(); it1 != fDisk1.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk4.begin(); it2 != fDisk4.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk1C->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk2.begin(); it1 != fDisk2.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk3.begin(); it2 != fDisk3.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk2A->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk2.begin(); it1 != fDisk2.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk4.begin(); it2 != fDisk4.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk2B->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk2.begin(); it1 != fDisk2.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk5.begin(); it2 != fDisk5.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk2C->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk3.begin(); it1 != fDisk3.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk4.begin(); it2 != fDisk4.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk3A->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk3.begin(); it1 != fDisk3.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk5.begin(); it2 != fDisk5.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk3B->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk3.begin(); it1 != fDisk3.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk6.begin(); it2 != fDisk6.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk3C->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk4.begin(); it1 != fDisk4.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk5.begin(); it2 != fDisk5.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk4A->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk4.begin(); it1 != fDisk4.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk6.begin(); it2 != fDisk6.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk4B->Fill(deltaZ);
          }
      }
      
      for(vector<TVector3>::iterator it1 = fDisk5.begin(); it1 != fDisk5.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = fDisk6.begin(); it2 != fDisk6.end(); ++it2)
          {
              Float_t fR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t fZ1 = (*it1).Z();
              Float_t fR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t fZ2 = (*it2).Z();
              Float_t recoZ = (fR2*fZ1 - fR1*fZ2)/(fR2 - fR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              disk5A->Fill(deltaZ);
          }
      }
      
      vector<TVector3> bLayer1;
      vector<TVector3> bLayer2;

      bLayer1.clear();
      bLayer2.clear();

      if( bRecHitN != 0 )
      {
          for(Int_t j = 0; j < bRecHitN; j++)
          {
              if( bRecHitLayer->at(j) == 1 )
              {
                  bLayer1.push_back(TVector3(bRecHitGx->at(j), bRecHitGy->at(j), bRecHitGz->at(j)));
              }
              if( bRecHitLayer->at(j) == 2 )
              {
                  bLayer2.push_back(TVector3(bRecHitGx->at(j), bRecHitGy->at(j), bRecHitGz->at(j)));
              }
          }
      }
      
      for(vector<TVector3>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
      {
          for(vector<TVector3>::iterator it2 = bLayer2.begin(); it2 != bLayer2.end(); ++it2)
          {
              Float_t bR1 = sqrt(pow((*it1).X(),2)+pow((*it1).Y(),2));
              Float_t bZ1 = (*it1).Z();
              Float_t bR2 = sqrt(pow((*it2).X(),2)+pow((*it2).Y(),2));
              Float_t bZ2 = (*it2).Z();
              Float_t recoZ = (bR2*bZ1 - bR1*bZ2)/(bR2 - bR1);
              Float_t deltaZ = simVz->at(0) - recoZ;
              a1->Fill(deltaZ);
          }
      }

   } // event loop
 
   a1->Write();
   h1->Write();
   h2->Write();
   h3->Write();
   h4->Write();
   h5->Write();
   h6->Write();

   disk1A->Write();
   disk1B->Write();
   disk1C->Write();
   disk2A->Write();
   disk2B->Write();
   disk2C->Write();
   disk3A->Write();
   disk3B->Write();
   disk3C->Write();
   disk4A->Write();
   disk4B->Write();
   disk5A->Write();

}
