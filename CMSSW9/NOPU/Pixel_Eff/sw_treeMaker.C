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

   Int_t nbins = 30; 
   Float_t x1 = 0.0; Float_t x2 = 3.0;
   Float_t dr_cut = 0.1;
  
   TH1F *bkg_plot = new TH1F("h1"," ; |#eta|; Pixel efficiency",nbins,x1,x2);
   TH1F *eff_la1 = new TH1F("eff_la1","1st Barrel",nbins,x1,x2);
   TH1F *eff_la2 = new TH1F("eff_la2","2nd Barrel",nbins,x1,x2);
   TH1F *eff_la3 = new TH1F("eff_la3","3rd Barrel",nbins,x1,x2);
   TH1F *eff_la4 = new TH1F("eff_la4","4th Barrel",nbins,x1,x2);
   TH1F *eff_di1 = new TH1F("eff_di1","1st Disk",nbins,x1,x2);
   TH1F *eff_di2 = new TH1F("eff_di2","2nd Disk",nbins,x1,x2);
   TH1F *eff_di3 = new TH1F("eff_di3","3rd Disk",nbins,x1,x2);
   TH1F *eff_di4 = new TH1F("eff_di4","4th Disk",nbins,x1,x2);
   TH1F *eff_di5 = new TH1F("eff_di5","5th Disk",nbins,x1,x2);
   TH1F *eff_di6 = new TH1F("eff_di6","6th Disk",nbins,x1,x2);

   Float_t denominator = 0.;
   Float_t Bnominator[4] = {};
   Float_t Fnominator[6] = {};

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

           Int_t index = 0;
           Float_t genPhi = propgenElPartPhi->at(0);     
           Float_t genEta = propgenElPartEta->at(0);
           Float_t EgEta = - 99.; Float_t EgPhi = -99.; Float_t EgEt = 0.;

           // ECAL loop
           Int_t closest_eg = -1;
           Float_t clo_dr = 9999.;
           for(Int_t i = 0; i < egCrysN; i++)
           {
               Float_t dPhi = deltaPhi(genPhi, egCrysPhi->at(i));
               Float_t tmp_dr = sqrt(pow(dPhi,2)+pow(genEta - egCrysEta->at(i),2));
               if( egCrysClusterEt->at(i) < 10. ) continue;
               if( tmp_dr < clo_dr )
               {
                   clo_dr = tmp_dr;
                   closest_eg = i;
               }
           } // ECAL loop end

           // HGCAL loop
           Float_t clo_cl3d_dr = 9999.;
           Int_t closest_cl3d = -1;
           for(Int_t i = 0; i < cl3d_n; i++)
           {
               if(cl3d_coreshowerlength->at(i) < 3 || cl3d_coreshowerlength->at(i) > 18 ) continue;
               if(cl3d_srrtot->at(i) < 0.002 || cl3d_srrtot->at(i) > 0.005 ) continue;
               if(cl3d_maxlayer->at(i) < 8 || cl3d_maxlayer->at(i) > 20 ) continue;
               if(cl3d_firstlayer->at(i) > 5 ) continue;

               Float_t dPhi = deltaPhi(genPhi, cl3d_phi->at(i));

               Float_t tmp_dr = sqrt(pow(genPhi - cl3d_phi->at(i),2)+pow(genEta - cl3d_eta->at(i),2));
               if( cl3d_pt->at(i) < 10. ) continue;
               if( tmp_dr < clo_cl3d_dr )
               {
                   clo_cl3d_dr = tmp_dr;
                   closest_cl3d = i;
               }
               
           }
           // HGCAL loop end

           // If there is no egamma both ECAL and HGCAL ignore it
           if( clo_dr == 9999. && clo_cl3d_dr == 9999. ) continue;

           // Determine which egamma is closer to gen electron direction
           if( (clo_dr != 9999. && clo_cl3d_dr != 9999.) || (clo_dr != 9999. && clo_cl3d_dr == 9999.) || (clo_dr == 9999. && clo_cl3d_dr != 9999.) ) 
           {
               if( clo_dr < clo_cl3d_dr )
               {
                   EgEt  = egCrysClusterEt->at(closest_eg);
                   EgEta = egCrysClusterEta->at(closest_eg);
                   EgPhi = egCrysClusterPhi->at(closest_eg);
                   index = 1;
               }
               else
               {
                   EgEt  = cl3d_pt->at(closest_cl3d);
                   EgEta = cl3d_eta->at(closest_cl3d);
                   EgPhi = cl3d_phi->at(closest_cl3d);
                   index = 2;
               }
           }

           // Ignore outer range egamma
           if( fabs(EgEta) > bin_right || fabs(EgEta) < bin_left ) continue;
      

           // Increase denominator because of presense of egamma
           denominator++;
 
           /*
           cout << " Process: " << jentry << endl;
           if( index == 1 ) cout << "------------------ ECAL egamma selected !! ------------------" << endl;
           if( index == 2 ) cout << "------------------ HGCAL egamma selected !! ------------------" << endl;
           cout << "  The closest egamma position and energy information" << endl;
           cout << "   ------------------ Position ------------------" << endl;
           cout << "     genPhi: " << genPhi << ", EgPhi: " << EgPhi << endl;
           cout << "     genEta: " << genEta << ", EgEta: " << EgEta << endl;
           cout << endl;
           cout << "   ------------------ Energy ------------------" << endl;
           cout << "     EgEt: " << EgEt << endl;
           cout << endl;
           */
           
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

           if( isFirst ) Bnominator[0]++;
           if( isSecond ) Bnominator[1]++;
           if( isThird ) Bnominator[2]++;
           if( isFourth ) Bnominator[3]++;

           Bool_t isDFirst = false, isDSecond = false, isDThird = false, isDFourth = false, isDFifth = false, isDSixth = false;
           for(Int_t j = 0; j < fRecHitN; j++)
           {
               Bool_t isDiskPhi = false;
               Float_t dist = fRecHitGz->at(j);
               TVector3 pixel;
               pixel.SetXYZ( fRecHitGx->at(j), fRecHitGy->at(j), fRecHitGz->at(j) );
               Float_t pixelPhi = pixel.Phi();
               Float_t deltaPhi = EgPhi - pixelPhi;
               if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
               if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
               if( fabs(deltaPhi) < 0.1 ) isDiskPhi = true;
               
               if( isDiskPhi )
               { 
                   if( fabs(dist) > 22 && fabs(dist) < 28 ) isDFirst = true;   
                   if( fabs(dist) > 29 && fabs(dist) < 34 ) isDSecond = true; 
                   if( fabs(dist) > 37 && fabs(dist) < 43 ) isDThird = true; 
                   if( fabs(dist) > 49 && fabs(dist) < 55 ) isDFourth = true; 
                   if( fabs(dist) > 63 && fabs(dist) < 70 ) isDFifth = true;
                   if( fabs(dist) > 80 && fabs(dist) < 90 ) isDSixth = true;
               }
           
           } // Disk pixel loop
           index = 0;
           
           if( isDFirst  ) Fnominator[0]++;
           if( isDSecond ) Fnominator[1]++;
           if( isDThird  ) Fnominator[2]++;
           if( isDFourth ) Fnominator[3]++;
           if( isDFifth  ) Fnominator[4]++;
           if( isDSixth  ) Fnominator[5]++;

       } // event loop
       
       cout << endl;
       cout << " From " << bin_left << " to " << bin_right << endl;
       cout << "  In this range, number of egamma is: " << denominator << endl;
       cout << "  At 1st layer: " << Bnominator[0] << endl;
       cout << "  At 2nd layer: " << Bnominator[1] << endl;
       cout << "  At 3rd layer: " << Bnominator[2] << endl;
       cout << "  At 4th layer: " << Bnominator[3] << endl;
       cout << "  At 1st Disk:  " << Fnominator[0] << endl;
       cout << "  At 2nd Disk:  " << Fnominator[1] << endl;
       cout << "  At 3rd Disk:  " << Fnominator[2] << endl;
       cout << "  At 4th Disk:  " << Fnominator[3] << endl;
       cout << "  At 5th Disk:  " << Fnominator[4] << endl;
       cout << "  At 6th Disk:  " << Fnominator[5] << endl;
       cout << endl;
       cout << "############# Move to next eta bin ###########" << endl;
       cout << endl;

       Float_t ratio1 = Bnominator[0]/denominator;
       Float_t ratio2 = Bnominator[1]/denominator;
       Float_t ratio3 = Bnominator[2]/denominator;
       Float_t ratio4 = Bnominator[3]/denominator;
       
       Float_t fratio1 = Fnominator[0]/denominator;
       Float_t fratio2 = Fnominator[1]/denominator;
       Float_t fratio3 = Fnominator[2]/denominator;
       Float_t fratio4 = Fnominator[3]/denominator;
       Float_t fratio5 = Fnominator[4]/denominator;
       Float_t fratio6 = Fnominator[5]/denominator;

       if( denominator == 0 ) {
           eff_la1->SetBinContent(bin+1, 0.);
           eff_la1->SetBinError(bin+1, 0.);
           eff_la2->SetBinContent(bin+1, 0.);
           eff_la2->SetBinError(bin+1, 0.);
           eff_la3->SetBinContent(bin+1, 0.);
           eff_la3->SetBinError(bin+1, 0.);
           eff_la4->SetBinContent(bin+1, 0.);
           eff_la4->SetBinError(bin+1, 0.);

           eff_di1->SetBinContent(bin+1, 0.);
           eff_di1->SetBinError(bin+1, 0.);
           eff_di2->SetBinContent(bin+1, 0.);
           eff_di2->SetBinError(bin+1, 0.);
           eff_di3->SetBinContent(bin+1, 0.);
           eff_di3->SetBinError(bin+1, 0.);
           eff_di4->SetBinContent(bin+1, 0.);
           eff_di4->SetBinError(bin+1, 0.);
           eff_di5->SetBinContent(bin+1, 0.);
           eff_di5->SetBinError(bin+1, 0.);
           eff_di6->SetBinContent(bin+1, 0.);
           eff_di6->SetBinError(bin+1, 0.);
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
           
           eff_di1->SetBinContent(bin+1, fratio1);
           eff_di1->SetBinError(bin+1, sqrt( fratio1 * (1-fratio1)/denominator));
           eff_di2->SetBinContent(bin+1, fratio2);
           eff_di2->SetBinError(bin+1, sqrt( fratio2 * (1-fratio2)/denominator));
           eff_di3->SetBinContent(bin+1, fratio3);
           eff_di3->SetBinError(bin+1, sqrt( fratio3 * (1-fratio3)/denominator));
           eff_di4->SetBinContent(bin+1, fratio4);
           eff_di4->SetBinError(bin+1, sqrt( fratio4 * (1-fratio4)/denominator));
           eff_di5->SetBinContent(bin+1, fratio5);
           eff_di5->SetBinError(bin+1, sqrt( fratio5 * (1-fratio4)/denominator));
           eff_di6->SetBinContent(bin+1, fratio6);
           eff_di6->SetBinError(bin+1, sqrt( fratio6 * (1-fratio4)/denominator));
       }

       // Initialize for next eta bin
       denominator = 0.;
       Bnominator[0] = 0.; Bnominator[1] = 0.; Bnominator[2] = 0.; Bnominator[3] = 0.;
       Fnominator[0] = 0.; Fnominator[1] = 0.; Fnominator[2] = 0.; Fnominator[3] = 0.; Fnominator[4] = 0.; Fnominator[5] = 0.;

  } // eta bin loop

  TFile *output = new TFile("eff_plot.root","RECREATE");     
  
  //TCanvas *c1 = new TCanvas("c1","",1366,768);
  //c1->SetLeftMargin(0.12);
  //c1->SetBottomMargin(0.12);

  bkg_plot->GetXaxis()->SetTitleSize(0.05);
  bkg_plot->GetXaxis()->CenterTitle(true);
  bkg_plot->GetXaxis()->SetNdivisions(515);
  bkg_plot->GetXaxis()->SetLabelSize(0.02);
  bkg_plot->GetYaxis()->SetTitleSize(0.05);
  bkg_plot->GetYaxis()->SetNdivisions(506);
  bkg_plot->GetYaxis()->SetLabelSize(0.02);
  bkg_plot->GetYaxis()->SetRangeUser(0, 1.02);
  bkg_plot->Write();

  eff_la1->SetMarkerStyle(20);
  eff_la1->SetMarkerColor(2);
  eff_la1->SetMarkerSize(1.4);
  eff_la1->SetLineWidth(1);
  eff_la1->SetLineColor(2);
  //eff_la1->Draw("HIST PL same");
  eff_la1->Write();
  
  eff_la2->SetMarkerStyle(21);
  eff_la2->SetMarkerColor(3);
  eff_la2->SetMarkerSize(1.4);
  eff_la2->SetLineWidth(1);
  eff_la2->SetLineColor(3);
  //eff_la2->Draw("HIST PL same");
  eff_la2->Write();
  
  eff_la3->SetMarkerStyle(22);
  eff_la3->SetMarkerColor(4);
  eff_la3->SetMarkerSize(1.4);
  eff_la3->SetLineWidth(1);
  eff_la3->SetLineColor(4);
  //eff_la3->Draw("HIST PL same");
  eff_la3->Write();
  
  eff_la4->SetMarkerStyle(29);
  eff_la4->SetMarkerColor(kYellow-1);
  eff_la4->SetMarkerSize(1.4);
  eff_la4->SetLineWidth(1);
  eff_la4->SetLineColor(kYellow-1);
  //eff_la4->Draw("HIST PL same");
  eff_la4->Write();
  
  eff_di1->SetMarkerStyle(33);
  eff_di1->SetMarkerColor(kRed-2);
  eff_di1->SetMarkerSize(1.4);
  eff_di1->SetLineWidth(1);
  eff_di1->SetLineColor(kRed-2);
  //eff_di1->Draw("HIST PL same");
  eff_di1->Write();
  
  eff_di2->SetMarkerStyle(34);
  eff_di2->SetMarkerColor(kBlue-2);
  eff_di2->SetMarkerSize(1.4);
  eff_di2->SetLineWidth(1);
  eff_di2->SetLineColor(kBlue-2);
  //eff_di2->Draw("HIST PL same");
  eff_di2->Write();
  
  eff_di3->SetMarkerStyle(29);
  eff_di3->SetMarkerColor(kGreen-2);
  eff_di3->SetMarkerSize(1.4);
  eff_di3->SetLineWidth(1);
  eff_di3->SetLineColor(kGreen-2);
  //eff_di3->Draw("HIST PL same");
  eff_di3->Write();
  
  eff_di4->SetMarkerStyle(22);
  eff_di4->SetMarkerColor(6);
  eff_di4->SetMarkerSize(1.4);
  eff_di4->SetLineWidth(1);
  eff_di4->SetLineColor(6);
  //eff_di4->Draw("HIST PL same");
  eff_di4->Write();
  
  eff_di5->SetMarkerStyle(21);
  eff_di5->SetMarkerColor(kCyan+1);
  eff_di5->SetMarkerSize(1.4);
  eff_di5->SetLineWidth(1);
  eff_di5->SetLineColor(kCyan+1);
  //eff_di5->Draw("HIST PL same");
  eff_di5->Write();
  
  eff_di6->SetMarkerStyle(20);
  eff_di6->SetMarkerColor(kOrange);
  eff_di6->SetMarkerSize(1.4);
  eff_di6->SetLineWidth(1);
  eff_di6->SetLineColor(kOrange);
  //eff_di6->Draw("HIST PL same");
  eff_di6->Write();

  //gPad->BuildLegend();
  gStyle->SetOptStat(0);
  
  //TLegend *lgd1 = new TLegend(0.2, 0.2, 0.55, 0.45);
  //lgd1->AddEntry(eff_la1, "1st Barrel" ,"lp");
  //lgd1->AddEntry(eff_la2, "2nd Barrel" ,"lp");
  //lgd1->AddEntry(eff_la3, "3rd Barrel" ,"lp");
  //lgd1->AddEntry(eff_la4, "4th Barrel" ,"lp");
  //lgd1->Draw();
  
  //c1->SaveAs("eff.png");
  
}
