#define bkg_cxx
#include "bkg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void bkg::Loop()
{
//   In a ROOT session, you can do:
//      root> .L bkg.C
//      root> bkg t
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

   vector<float> v1;
   TH1F *h1 = new TH1F("h1","",100,0,1);
   TH1F *a1 = new TH1F("a1","",100,0,1);
   TH2F *h2 = new TH2F("h2","",15,0,15,50,0,1);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      Float_t nomi = 0.;
      Float_t denomi = 0.;

      v1.clear();
      if( !ntnEg3 ) continue;
      if( recon_pT->size() == 0 ) continue;
      
      Int_t eg_size = ntEgEt_iso->size();
      
      for(Int_t j = 0; j < eg_size; j++)
      {
          if( ntEgEt_iso->at(j) < 20. ) continue;
          if( recon_pT->size() == 1 )
          {
              h1->Fill(recon_pT->at(0));
              h2->Fill(v1.size(), recon_pT->at(0));
          }

          if( recon_pT->size() >= 2 )
          {
              for(Int_t i = 0; i < recon_pT->size(); i++) v1.push_back(recon_pT->at(i));

              sort(v1.begin(), v1.end());

              for(Int_t j = 0; j < v1.size();   j++) denomi += v1.at(j);
              for(Int_t k = 0; k < v1.size()-1; k++) nomi += v1.at(k);

              Float_t ratio = nomi/denomi;
              h1->Fill(ratio);
              h2->Fill(v1.size(), ratio);
              a1->Fill(ratio);
          }
      }
      
   }


   TFile *output = new TFile("plot_bkg.root","RECREATE");
   h1->Write();
   h2->Write();
   a1->Write();
   output->Close();


}
