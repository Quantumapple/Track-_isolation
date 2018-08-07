#define sw_treeMaker_cxx
#include "sw_treeMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void sw_treeMaker::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   Int_t nbins = 15; 
   Float_t x1 = 0.0; Float_t x2 = 1.5;
  
   TH1F *bkg_plot = new TH1F("h1"," ; |#eta|; Pixel efficiency",nbins,x1,x2);
   TH1F *eff_la1 = new TH1F("eff_la1","1st Barrel",nbins,x1,x2);
   TH1F *eff_la2 = new TH1F("eff_la2","2nd Barrel",nbins,x1,x2);
   TH1F *eff_la3 = new TH1F("eff_la3","3rd Barrel",nbins,x1,x2);
   TH1F *eff_la4 = new TH1F("eff_la4","4th Barrel",nbins,x1,x2);

   Float_t denominator = 0.;
   Float_t nominator[4] = {};

   for(Int_t bin = 0; bin < nbins; bin++)
   {
       Float_t bin_left = (Float_t)bin*0.1;
       Float_t bin_right = bin_left + 0.1;
       cout << "Range from: " << bin_left << " to " << bin_right << endl;

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
           //for (Long64_t jentry=0; jentry<15;jentry++) {
           Long64_t ientry = LoadTree(jentry);
           if (ientry < 0) break;
           nb = fChain->GetEntry(jentry);   nbytes += nb;
           // if (Cut(ientry) < 0) continue;

           Float_t genPhi = propgenElPartPhi->at(0);     
           Float_t genEta = propgenElPartEta->at(0);

           Float_t EgEta = - 99.; Float_t EgPhi = -99.; Float_t EgEt = 0.;
           Float_t clo_dr = 9999.;
           for(Int_t i = 0; i < egCrysN; i++)
           {
               Float_t dPhi = deltaPhi(genPhi, egCrysPhi->at(i));
               Float_t tmp_dr = sqrt(pow(dPhi,2)+pow(genEta - egCrysEta->at(i),2));
               if( tmp_dr < clo_dr )
               {
                   clo_dr = tmp_dr;
                   EgPhi = egCrysPhi->at(i);
                   EgEta = egCrysEta->at(i);
                   EgEt  = egCrysEt->at(i);
               }
           }

           // Ignore when event cannot find closest egamma to gen electron
           if( dr == 999. ) 
           {
               cout << " This event seems something wrong" << endl;
               cout << "  Process: " << jentry << endl;
               continue;
           }

           // Ignore low Et egamma
           if( EgEt < 10. ) continue;
           // Ignore outer range egamma
           if( fabs(EgEta) > bin_right || fabs(EgEta) < bin_left ) continue;
      

           // Increase denominator because of presense of egamma
           denominator++;

           cout << " Process: " << jentry << endl;
           cout << "  The closest egamma position and energy information" << endl;
           cout << "   ------------------ Position ------------------" << endl;
           cout << "     genPhi: " << genPhi << ", EgPhi: " << EgPhi << endl;
           cout << "     genEta: " << genEta << ", EgEta: " << EgEta << endl;
           cout << endl;
           cout << "   ------------------ Energy ------------------" << endl;
           cout << "     EgEt: " << EgEt << endl;
           cout << endl;


           Bool_t isFirst = false, isSecond = false, isThird = false, isFourth = false;
           for(Int_t j = 0; j < bRecHitN; j++)
           {
               Bool_t isBarrelPhi = false;
               Float_t Rad = sqrt(pow(bRecHitGx->at(j),2)+pow(bRecHitGy->at(j),2));
               TVector3 pixel;
               pixel.SetXYZ( bRecHitGx->at(j), bRecHitGy->at(j), bRecHitGz->at(j) );
               Float_t pixelPhi = pixel.Phi();
               Float_t dPhi = deltaPhi(EgPhi, pixelPhi);
               if( fabs(dPhi) < 0.1 ) isBarrelPhi = true;

               if( isBarrelPhi)
               {
                   if( Rad < 5. ) isFirst = true;
                   if( Rad > 5.  && Rad < 9.  ) isSecond = true;
                   if( Rad > 9.  && Rad < 14. ) isThird = true;
                   if( Rad > 14. && Rad < 18. ) isFourth = true;
               }
           } // barrel pixel loop

           if( isFirst ) nominator[0]++;
           if( isSecond ) nominator[1]++;
           if( isThird ) nominator[2]++;
           if( isFourth ) nominator[3]++;

       } // event loop
   
       cout << endl;
       cout << "In this range, number of egamma is: " << denominator << endl;
       cout << "  At 1st layer: " << nominator[0] << endl;
       cout << "  At 2nd layer: " << nominator[1] << endl;
       cout << "  At 3rd layer: " << nominator[2] << endl;
       cout << "  At 4th layer: " << nominator[3] << endl;
       cout << endl;
       cout << "So the efficiency of pixel detectors are: " << endl;
       cout << "  At 1st layer: " << nominator[0]/denominator << endl;
       cout << "  At 2nd layer: " << nominator[1]/denominator << endl;
       cout << "  At 3rd layer: " << nominator[2]/denominator << endl;
       cout << "  At 4th layer: " << nominator[3]/denominator << endl;
       cout << endl;
       
       cout << "############# Move to next eta bin ###########" << endl << endl;
 
       Float_t ratio1 = nominator[0]/denominator;
       Float_t ratio2 = nominator[1]/denominator;
       Float_t ratio3 = nominator[2]/denominator;
       Float_t ratio4 = nominator[3]/denominator;

       if( denominator == 0 ) {
           eff_la1->SetBinContent(bin+1, 0.);
           eff_la1->SetBinError(bin+1, 0.);
           eff_la2->SetBinContent(bin+1, 0.);
           eff_la2->SetBinError(bin+1, 0.);
           eff_la3->SetBinContent(bin+1, 0.);
           eff_la3->SetBinError(bin+1, 0.);
           eff_la4->SetBinContent(bin+1, 0.);
           eff_la4->SetBinError(bin+1, 0.);
       }
       else {
           eff_la1->SetBinContent(bin+1, ratio1);
           eff_la1->SetBinError(bin+1, sqrt( ratio1 * (1-ratio1)/denominator));
           eff_la2->SetBinContent(bin+1, ratio2);
           eff_la2->SetBinError(bin+1, sqrt( ratio2 * (1-ratio2)/denominator));
           eff_la3->SetBinContent(bin+1, ratio3);
           eff_la3->SetBinError(bin+1, sqrt( ratio3 * (1-ratio3)/denominator));
           eff_la4->SetBinContent(bin+1, ratio4);
           eff_la4->SetBinError(bin+1, sqrt( ratio4 * (1-ratio4)/denominator));
       }

       // Initialize for next eta bin
       denominator = 0.;
       nominator[0] = 0.; nominator[1] = 0.; nominator[2] = 0.; nominator[3] = 0.;
        
  } // eta bin loop

  TCanvas *c1 = new TCanvas("c1","",1366,768);
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);

  bkg_plot->GetXaxis()->SetTitleSize(0.05);
  bkg_plot->GetXaxis()->CenterTitle(true);
  bkg_plot->GetXaxis()->SetNdivisions(515);
  bkg_plot->GetXaxis()->SetLabelSize(0.02);
  bkg_plot->GetYaxis()->SetTitleSize(0.05);
  bkg_plot->GetYaxis()->SetNdivisions(506);
  bkg_plot->GetYaxis()->SetLabelSize(0.02);
  bkg_plot->GetYaxis()->SetRangeUser(0, 1.02);
  bkg_plot->Draw();

  eff_la1->SetMarkerStyle(20);
  eff_la1->SetMarkerColor(2);
  eff_la1->SetMarkerSize(1.4);
  eff_la1->SetLineWidth(1);
  eff_la1->SetLineColor(2);
  eff_la1->Draw("HIST PL same");
  
  eff_la2->SetMarkerStyle(21);
  eff_la2->SetMarkerColor(3);
  eff_la2->SetMarkerSize(1.4);
  eff_la2->SetLineWidth(1);
  eff_la2->SetLineColor(3);
  eff_la2->Draw("HIST PL same");
  
  eff_la3->SetMarkerStyle(22);
  eff_la3->SetMarkerColor(4);
  eff_la3->SetMarkerSize(1.4);
  eff_la3->SetLineWidth(1);
  eff_la3->SetLineColor(4);
  eff_la3->Draw("HIST PL same");
  
  eff_la4->SetMarkerStyle(29);
  eff_la4->SetMarkerColor(kYellow-1);
  eff_la4->SetMarkerSize(1.4);
  eff_la4->SetLineWidth(1);
  eff_la4->SetLineColor(kYellow-1);
  eff_la4->Draw("HIST PL same");

  //gPad->BuildLegend();
  gStyle->SetOptStat(0);
  
  TLegend *lgd1 = new TLegend(0.2, 0.2, 0.55, 0.45);
  lgd1->AddEntry(eff_la1, "1st Barrel" ,"lp");
  lgd1->AddEntry(eff_la2, "2nd Barrel" ,"lp");
  lgd1->AddEntry(eff_la3, "3rd Barrel" ,"lp");
  lgd1->AddEntry(eff_la4, "4th Barrel" ,"lp");
  lgd1->Draw();
  
  c1->SaveAs("eff.png");
}
