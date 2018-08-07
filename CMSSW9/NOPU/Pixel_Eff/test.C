#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TVector3.h>

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   int nbins = 16;
   float x1 = 2.0;
   float x2 = 2.8;
   float nth_bin_left;
   float nth_bin_right;

   TH1F *bkg_plot = new TH1F("h1"," ; |#eta|; Pixel efficiency",nbins,x1,x2); 
   TH1F *L1 = new TH1F("L1","L1",nbins,x1,x2);
   TH1F *D1 = new TH1F("D1","D1",nbins,x1,x2);
   TH1F *D2 = new TH1F("D2","D2",nbins,x1,x2);
   TH1F *D3 = new TH1F("D3","D3",nbins,x1,x2);
 
   Float_t denominator = 0.;
   Float_t nominator[4] = {};
 
   for(int i = 0; i < nbins; i++){
	   nth_bin_left = x1 + (i*0.05);
	   nth_bin_right = x1 + ((i+1)*0.05);
       
       for (Long64_t jentry=0; jentry<nentries;jentry++) {
           Long64_t ientry = LoadTree(jentry);
           if (ientry < 0) break;
           nb = fChain->GetEntry(jentry);   nbytes += nb;
           // if (Cut(ientry) < 0) continue;
      
           float genPhi = genPartPhi->at(0);
           float genEta = genPartEta->at(0);
           float ME0Eta;
           float ME0Phi;
      
           float dr = 9999.;
           for(int j = 0; j < me0SegNum; j++){
               TVector3 me0SegPos;
               me0SegPos.SetXYZ(me0SegPosX->at(j),me0SegPosY->at(j),me0SegPosZ->at(j));
               Float_t me0Eta = me0SegPos.Eta();
               Float_t me0Phi = me0SegPos.Phi();
               Float_t dPhi = deltaPhi(genPhi, me0Phi);
               Float_t deltaEta = genEta - me0Eta;
               Float_t tmp_dr = sqrt(pow(dPhi,2)+pow(deltaEta,2));

               if( tmp_dr < dr )
               {
                   dr = tmp_dr;
                   ME0Eta = me0Eta;
                   ME0Phi = me0Phi;
               }
           } // ME0 loop

           // Ignore if we can't find the closest muon to gen direction
           if( dr == 9999. )
           {
               //cout << "Something wrong" << endl;
               continue;
           }

           // Ignore outer range muon
           if( fabs(ME0Eta) > nth_bin_right || fabs(ME0Eta) < nth_bin_left ) continue;

           denominator++;

           Bool_t isFirst = false;
           for(Int_t k = 0; k < bRecHitN; k++)
           {
               Bool_t isBarrelPhi = false;
               Float_t Rad = sqrt(pow(bRecHitGx->at(k),2)+pow(bRecHitGy->at(k),2));
               TVector3 pixel;
               pixel.SetXYZ(bRecHitGx->at(k),bRecHitGy->at(k),bRecHitGz->at(k));
               Float_t pixelPhi = pixel.Phi();
               Float_t dPhi = deltaPhi(ME0Phi, pixelPhi);
               if(fabs(dPhi) < 0.1 ) isBarrelPhi = true;

               if( isBarrelPhi ) 
                   if( Rad < 4. ) isFirst = true;
               
           } // barrel pixel loop

           if( isFirst ) nominator[0]++; 
           
           Bool_t isDiskFirst = false, isDiskSecond = false, isDiskThird = false;
           for(int k = 0; k < fRecHitN; k++){
               Bool_t isDiskPhi = false;
               float diskZ = fRecHitGz->at(k);
               TVector3 fRecHit;
               fRecHit.SetXYZ(fRecHitGx->at(k),fRecHitGy->at(k),fRecHitGz->at(k));
               Float_t pixelPhi = fRecHit.Phi();
               Float_t dPhi = deltaPhi(ME0Phi, pixelPhi);
               if( fabs(dPhi) < 0.1 ) isDiskPhi = true;

               if( isDiskPhi )
               {
                   if( fabs(diskZ) > 22 && fabs(diskZ) < 28 ) isDiskFirst = true; 
                   if( fabs(diskZ) > 29 && fabs(diskZ) < 34 ) isDiskSecond = true;
                   if( fabs(diskZ) > 37 && fabs(diskZ) < 43 ) isDiskThird = true;
               }
           } // disk pixel loop

           if( isDiskFirst ) nominator[1]++;
           if( isDiskSecond ) nominator[2]++;
           if( isDiskThird ) nominator[3]++;
                  

       } // event loop

       Float_t eff_la1 = nominator[0]/denominator;
       Float_t eff_di1 = nominator[1]/denominator;
       Float_t eff_di2 = nominator[2]/denominator;
       Float_t eff_di3 = nominator[3]/denominator;

       if( denominator == 0. ) { 
           L1->SetBinContent(i+1, 0.);
           L1->SetBinError(i+1, 0.); 
           D1->SetBinContent(i+1, 0.);
           D1->SetBinError(i+1, 0.); 
           D2->SetBinContent(i+1, 0.);
           D2->SetBinError(i+1, 0.); 
           D3->SetBinContent(i+1, 0.);
           D3->SetBinError(i+1, 0.); 
       }
       else {
           L1->SetBinContent(i+1, eff_la1);
           L1->SetBinError(i+1, sqrt( eff_la1 * (1-eff_la1) / denominator) ); 
           D1->SetBinContent(i+1, eff_di1);
           D1->SetBinError(i+1, sqrt( eff_di1 * (1-eff_di1) / denominator) ); 
           D2->SetBinContent(i+1, eff_di2);
           D2->SetBinError(i+1, sqrt( eff_di2 * (1-eff_di1) / denominator) ); 
           D3->SetBinContent(i+1, eff_di3);
           D3->SetBinError(i+1, sqrt( eff_di3 * (1-eff_di1) / denominator) ); 
       }

       cout << "Eta bin from " << nth_bin_left << " to " << nth_bin_right << endl;
       cout << "At 1st layer: " << nominator[0]/denominator << endl;
       cout << "At 1st disk: " << nominator[1]/denominator << endl;
       cout << "At 2nd disk: " << nominator[2]/denominator << endl;
       cout << "At 3rd disk: " << nominator[3]/denominator << endl;
       cout << endl;
       cout << "Move to next bin" << endl;
       cout << endl;
   
       denominator = 0.;
       nominator[0] = 0.; nominator[1] = 0.; nominator[2] = 0.; nominator[3] = 0.;

   } // eta bin loop

  TCanvas *c1 = new TCanvas("c1","",1366,768);
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);

  bkg_plot->GetXaxis()->SetTitleSize(0.05);
  bkg_plot->GetXaxis()->CenterTitle(true);
  bkg_plot->GetXaxis()->SetNdivisions(516);
  bkg_plot->GetXaxis()->SetLabelSize(0.02);
  bkg_plot->GetYaxis()->SetTitleSize(0.05);
  bkg_plot->GetYaxis()->SetNdivisions(506);
  bkg_plot->GetYaxis()->SetLabelSize(0.02);
  bkg_plot->GetYaxis()->SetRangeUser(0, 1.02);
  bkg_plot->Draw();

  L1->SetMarkerStyle(20);
  L1->SetMarkerColor(2);
  L1->SetMarkerSize(1.4);
  L1->SetLineWidth(1);
  L1->SetLineColor(2);
  L1->Draw("HIST PL same");
  
  D1->SetMarkerStyle(21);
  D1->SetMarkerColor(3);
  D1->SetMarkerSize(1.4);
  D1->SetLineWidth(1);
  D1->SetLineColor(3);
  D1->Draw("HIST PL same");
  
  D2->SetMarkerStyle(22);
  D2->SetMarkerColor(4);
  D2->SetMarkerSize(1.4);
  D2->SetLineWidth(1);
  D2->SetLineColor(4);
  D2->Draw("HIST PL same");
  
  D3->SetMarkerStyle(29);
  D3->SetMarkerColor(kYellow-1);
  D3->SetMarkerSize(1.4);
  D3->SetLineWidth(1);
  D3->SetLineColor(kYellow-1);
  D3->Draw("HIST PL same");

  gStyle->SetOptStat(0);
  
  TLegend *lgd1 = new TLegend(0.2, 0.2, 0.55, 0.45);
  lgd1->AddEntry(L1, "1st Barrel" ,"lp");
  lgd1->AddEntry(D1, "1st Disk" ,"lp");
  lgd1->AddEntry(D2, "2nd Disk" ,"lp");
  lgd1->AddEntry(D3, "3rd Disk" ,"lp");
  lgd1->Draw();
  
  c1->SaveAs("eff.png");

}
