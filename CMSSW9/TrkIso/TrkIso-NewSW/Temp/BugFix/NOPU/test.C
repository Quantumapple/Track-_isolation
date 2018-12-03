#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *output = new TFile("new_draw.root","RECREATE");

   Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=285; jentry<286;jentry++) {
   //for (Long64_t jentry=404; jentry<405;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      //Float_t genPhi = genPartPhi->at(1);
      //Float_t genEta = genPartEta->at(1);
      Float_t genPhi = propgenElPartPhi->at(0);
      Float_t genEta = propgenElPartEta->at(0);

      Float_t pixR1 = 0, pixZ1 = 0.;
      Float_t pixR2 = 0, pixZ2 = 0.;
      Float_t pixR3 = 0, pixZ3 = 0.;
      Float_t pixR4 = 0, pixZ4 = 0.;
      
      Float_t pixX1 = 0, pixY1 = 0.;
      Float_t pixX2 = 0, pixY2 = 0.;
      Float_t pixX3 = 0, pixY3 = 0.;
      Float_t pixX4 = 0, pixY4 = 0.;

      for(Int_t i = 0; i < bRecHitN; i++) {
          Int_t flag = bRecHitLayer->at(i);
          TVector3 pixel;
          pixel.SetXYZ( bRecHitGx->at(i) - simVx->at(1), bRecHitGy->at(i) - simVy->at(1), bRecHitGz->at(i) - simVz->at(1) );
          Float_t pixelPhi = pixel.Phi();
          Float_t pixelEta = pixel.Eta();
          Float_t dPhi = deltaPhi(genPhi, pixelPhi);
          if( fabs(dPhi) > 0.3 ) continue;
          Float_t dEta = genEta - pixelEta;
          Float_t DR = sqrt(pow(dPhi,2)+pow(dEta,2));
          if( DR > 0.3 ) continue;
          Float_t R = sqrt(pow(bRecHitGx->at(i),2)+pow(bRecHitGy->at(i),2));

          switch(flag) {
              case 1:
                  pixR1 = R; pixZ1 = bRecHitGz->at(i);
                  pixX1 = bRecHitGx->at(i); pixY1 = bRecHitGy->at(i);
                  break;

              case 2:
                  pixR2 = R; pixZ2 = bRecHitGz->at(i);
                  pixX2 = bRecHitGx->at(i); pixY2 = bRecHitGy->at(i);
                  break;
              
              case 3:
                  pixR3 = R; pixZ3 = bRecHitGz->at(i);
                  pixX3 = bRecHitGx->at(i); pixY3 = bRecHitGy->at(i);
                  break;
              
              case 4:
                  pixR4 = R; pixZ4 = bRecHitGz->at(i);
                  pixX4 = bRecHitGx->at(i); pixY4 = bRecHitGy->at(i);
                  break;
          
              default:
                  break;
          }
      
      } // pixel loop
 
      Float_t simVr = sqrt(pow(simVx->at(1),2)+pow(simVy->at(1),2));

      Int_t egCrysN = egCrysE->size();
      Int_t egCrysClusterN = egCrysClusterE->size();

      cout << "Gen information" << endl;
      cout << genEta << endl;
      cout << genPhi << endl;
      cout << endl;

      Float_t L1_EgCrysR = sqrt(pow(egCrysGx->at(0),2)+pow(egCrysGy->at(0),2));
      Float_t L1_EgCrysTheta = atan(L1_EgCrysR/egCrysGz->at(0));
      if( L1_EgCrysTheta < 0 ) L1_EgCrysTheta += TMath::Pi();
      Float_t L1_EgCrysEta = -log(tan(L1_EgCrysTheta/2.0));
     
      TVector3 EgCrysV1;
      EgCrysV1.SetXYZ( egCrysGx->at(0), egCrysGy->at(0), egCrysGz->at(0) );

      TVector3 EgCrysV2;
      EgCrysV2.SetXYZ( egCrysGx->at(0) - simVx->at(1), egCrysGy->at(0) - simVy->at(1), egCrysGz->at(0) - simVz->at(1) );
      
      cout << "EgCrys" << endl;
      cout << "From eta branch : " << egCrysEta->at(0) << endl;
      cout << "From calculation : " << L1_EgCrysEta << endl;
      cout << "From TVector3 w/o BS info: " << EgCrysV1.Eta() << endl;
      cout << "From TVector3 w BS info: " << EgCrysV2.Eta() << endl;
      
      cout << endl; 
      Float_t L1_EgCrysClusterR = sqrt(pow(egCrysClusterGx->at(0),2)+pow(egCrysClusterGy->at(0),2));
      Float_t L1_EgCrysClusterTheta = atan(L1_EgCrysClusterR/egCrysClusterGz->at(0));
      if( L1_EgCrysClusterTheta < 0 ) L1_EgCrysClusterTheta += TMath::Pi();
      Float_t L1_EgCrysClusterEta = -log(tan(L1_EgCrysClusterTheta/2.0));
      Float_t L1_RZlength = sqrt(pow(L1_EgCrysClusterR,2)+pow(egCrysClusterGz->at(0),2));

      TVector3 EgCrysClusterV1;
      EgCrysClusterV1.SetXYZ( egCrysClusterGx->at(0), egCrysClusterGy->at(0), egCrysClusterGz->at(0) );

      TVector3 EgCrysClusterV2;
      EgCrysClusterV2.SetXYZ( egCrysClusterGx->at(0) - simVx->at(1), egCrysClusterGy->at(0) - simVy->at(1), egCrysClusterGz->at(0) - simVz->at(1) );
      
      cout << "EgCrysCluster" << endl;
      cout << "From eta branch : " << egCrysClusterEta->at(0) << endl;
      cout << "From calculation : " << L1_EgCrysClusterEta << endl;
      cout << "From TVector3 w/o BS info: " << EgCrysClusterV1.Eta() << endl;
      cout << "From TVector3 w BS info: " << EgCrysClusterV2.Eta() << endl;

      TGraph *pixelRZ = new TGraph(); 
      pixelRZ->SetPoint(0, simVz->at(1), simVr);
      pixelRZ->SetPoint(1, pixZ1, pixR1);
      pixelRZ->SetPoint(2, pixZ2, pixR2);
      pixelRZ->SetPoint(3, pixZ3, pixR3);
      pixelRZ->SetPoint(4, pixZ4, pixR4);

      pixelRZ->SetMarkerStyle(29);
      pixelRZ->SetMarkerSize(2.5);
      pixelRZ->SetMarkerColor(kAzure);
      pixelRZ->SetLineWidth(2);
      pixelRZ->SetLineColor(kAzure);
      
      TGraph *pixelXY = new TGraph(); 
      pixelXY->SetPoint(0, simVx->at(1), simVy->at(1));
      pixelXY->SetPoint(1, pixX1, pixY1);
      pixelXY->SetPoint(2, pixX2, pixY2);
      pixelXY->SetPoint(3, pixX3, pixY3);
      pixelXY->SetPoint(4, pixX4, pixY4);

      pixelXY->SetMarkerStyle(29);
      pixelXY->SetMarkerSize(2.5);
      pixelXY->SetMarkerColor(kAzure);
      pixelXY->SetLineWidth(2);
      pixelXY->SetLineColor(kAzure);
  
      TGraph *propGenXY = new TGraph();
      propGenXY->SetPoint(0, simVx->at(1), simVy->at(1));
      propGenXY->SetPoint(1, L1_EgCrysClusterR*cos(genPhi), L1_EgCrysClusterR*sin(genPhi));

      propGenXY->SetMarkerStyle(20);
      propGenXY->SetMarkerSize(1.5);
      propGenXY->SetMarkerColor(kRed);
      propGenXY->SetLineWidth(2);
      propGenXY->SetLineColor(kRed);

      Float_t propGenTheta = 2.*atan(exp(-genEta));
      TGraph *propGenRZ = new TGraph();
      propGenRZ->SetPoint(0, simVz->at(1), simVr);
      propGenRZ->SetPoint(1, L1_RZlength*cos(propGenTheta), L1_RZlength*sin(propGenTheta));

      propGenRZ->SetMarkerStyle(20);
      propGenRZ->SetMarkerSize(1.5);
      propGenRZ->SetMarkerColor(kRed);
      propGenRZ->SetLineWidth(2);
      propGenRZ->SetLineColor(kRed);

      TGraph *GenXY = new TGraph();
      GenXY->SetPoint(0, simVx->at(1), simVy->at(1));
      GenXY->SetPoint(1, L1_EgCrysClusterR*cos(genPartPhi->at(1)), L1_EgCrysClusterR*sin(genPartPhi->at(1)));

      GenXY->SetMarkerStyle(20);
      GenXY->SetMarkerSize(1.5);
      GenXY->SetMarkerColor(kOrange);
      GenXY->SetLineWidth(2);
      GenXY->SetLineColor(kOrange);

      Float_t GenTheta = 2.*atan(exp(-genPartEta->at(1)));
      TGraph *GenRZ = new TGraph();
      GenRZ->SetPoint(0, simVz->at(1), simVr);
      GenRZ->SetPoint(1, L1_RZlength*cos(GenTheta), L1_RZlength*sin(GenTheta));

      GenRZ->SetMarkerStyle(20);
      GenRZ->SetMarkerSize(1.5);
      GenRZ->SetMarkerColor(kOrange);
      GenRZ->SetLineWidth(2);
      GenRZ->SetLineColor(kOrange);

      TGraph *L1EgXY = new TGraph();
      L1EgXY->SetPoint(0, simVx->at(1), simVy->at(1));
      L1EgXY->SetPoint(1, L1_EgCrysClusterR*cos(egCrysClusterPhi->at(0)), L1_EgCrysClusterR*sin(egCrysClusterPhi->at(0)));

      L1EgXY->SetMarkerStyle(22);
      L1EgXY->SetMarkerSize(1.5);
      L1EgXY->SetMarkerColor(kViolet);
      L1EgXY->SetLineWidth(2);
      L1EgXY->SetLineColor(kViolet);

      Float_t L1_Theta = 2.*atan(exp(-egCrysClusterEta->at(0)));
      TGraph *L1EgRZ = new TGraph();
      L1EgRZ->SetPoint(0, simVz->at(1), simVr);
      L1EgRZ->SetPoint(1, L1_RZlength*cos(L1_Theta), L1_RZlength*sin(L1_Theta));

      L1EgRZ->SetMarkerStyle(22);
      L1EgRZ->SetMarkerSize(1.5);
      L1EgRZ->SetMarkerColor(kViolet);
      L1EgRZ->SetLineWidth(2);
      L1EgRZ->SetLineColor(kViolet);

      //pixelXY->Draw("P same");
      //propGenXY->Draw("LP same");
      //GenXY->Draw("LP same");
      
      pixelXY->Write("pixelXY");
      propGenXY->Write("propGenXY");
      GenXY->Write("GenXY");
      L1EgXY->Write("L1EgXY");

      //TLegend *lgd1 = new TLegend(0.3, 0.3, 0.6, 0.7);
      //lgd1->AddEntry(pixelXY, "Pixel clusters", "p");
      //lgd1->AddEntry(propGenXY, "Propgen direction", "lp");
      //lgd1->AddEntry(GenXY, "Gen direction", "lp");
      //lgd1->Draw();
      //lgd1->Write("LegendBoxXY");

      //pixelRZ->Draw("P same");
      //propGenRZ->Draw("LP same");
      //GenRZ->Draw("LP same");
      
      pixelRZ->Write("pixelRZ");
      propGenRZ->Write("propGenRZ");
      GenRZ->Write("GenRZ");
      L1EgRZ->Write("L1EgRZ");

      //TLegend *lgd2 = new TLegend(0.3, 0.3, 0.6, 0.7);
      //lgd2->AddEntry(pixelRZ, "Pixel clusters", "p");
      //lgd2->AddEntry(propGenRZ, "Propgen direction", "lp");
      //lgd2->AddEntry(GenRZ, "Gen direction", "lp");
      //lgd2->Draw();
      //lgd2->Write("LegendBoxRZ");
      

   } // event loop

   output->Close();
}
