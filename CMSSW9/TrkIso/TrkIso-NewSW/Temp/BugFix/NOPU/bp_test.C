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
   for (Long64_t jentry=0; jentry<50000;jentry++) {
   //for (Long64_t jentry=285; jentry<286;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      Float_t genPt = propgenElPartPt->at(0);

      if( genPt > 38. || genPt < 37. ) continue;

      cout << jentry << endl;
      for(Int_t i = 0; i < bRecHitN; i++)
      {
          Int_t flag = bRecHitLayer->at(i);
          cout << i << ", " << flag << endl;
      }
      cout << endl;

   } // event loop
  
   
}
