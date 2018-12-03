#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <TLorentzVector.h>
#include <TGraph.h>
using namespace std;

void test::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    float dr_cut = 0.1;
    
    TH1F *aPix1 = new TH1F("allPix1","",250,0,15000);
    TH1F *aPix2 = new TH1F("allPix2","",250,0,15000);
    TH1F *aPix3 = new TH1F("allPix3","",250,0,15000);
    TH1F *aPix4 = new TH1F("allPix4","",250,0,15000);

    TH1F *L1PixdP1 = new TH1F("L1PixdP1","",50,0,500);
    TH1F *L1PixdP2 = new TH1F("L1PixdP2","",50,0,500);
    TH1F *L1PixdP3 = new TH1F("L1PixdP3","",50,0,500);
    TH1F *L1PixdP4 = new TH1F("L1PixdP4","",50,0,500);

    TH1F *L1PixdR1 = new TH1F("L1PixdR1","",25,0,50);
    TH1F *L1PixdR2 = new TH1F("L1PixdR2","",25,0,50);
    TH1F *L1PixdR3 = new TH1F("L1PixdR3","",25,0,50);
    TH1F *L1PixdR4 = new TH1F("L1PixdR4","",25,0,50);
    
    TH1F *PixTrkdP1 = new TH1F("PixTrkdP1","",150,0,1500);
    TH1F *PixTrkdP2 = new TH1F("PixTrkdP2","",150,0,1500);
    TH1F *PixTrkdP3 = new TH1F("PixTrkdP3","",150,0,1500);
    TH1F *PixTrkdP4 = new TH1F("PixTrkdP4","",150,0,1500);

    TH1F *PixTrkdR1 = new TH1F("PixTrkdR1","",50,0,250);
    TH1F *PixTrkdR2 = new TH1F("PixTrkdR2","",50,0,250);
    TH1F *PixTrkdR3 = new TH1F("PixTrkdR3","",50,0,250);
    TH1F *PixTrkdR4 = new TH1F("PixTrkdR4","",50,0,250);
    
    int regionFlag = 0;

    //for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
    for (Long64_t jentry=0; jentry<500;jentry++) { //nentries
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;

        first_hits.clear();
        second_hits.clear();
        third_hits.clear();
        fourth_hits.clear();

        float genPt = propgenElPartPt->at(0);  

        //find closest egamma object to the gen electron
        EgN=egCrysClusterEt->size();
        int cl3d_N_ = cl3d_pt->size();

        float closest_dr = 9999.;
        int closest_eg = 0;

        Int_t TrkIsodP[4] = {};
        Int_t TrkIsodR[4] = {};
        Int_t L1PixdP[4] = {};
        Int_t L1PixdR[4] = {};


        cout << "Check intialization" << endl;
        cout << " L1 pixel: " << L1PixdP[0] << ", " << L1PixdP[1] << ", " << L1PixdP[2] << ", " << L1PixdP[3] << endl;
        cout << " TrkIso:  " << TrkIsodP[0] << ", " << L1PixdP[1] << ", " << L1PixdP[2] << ", " << L1PixdP[3] << endl;

        // loop over barrel egamma objects
        for(int i=0; i < EgN;i++){

            float dPhi = deltaPhi(propgenElPartPhi->at(0), egCrysClusterPhi->at(i));

            float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
            if(egCrysClusterEt->at(i) < 10) continue;
            if(current_dr < closest_dr){
                closest_dr = current_dr;
                closest_eg = i;
            }
        }// end of loop to find the closest egamma to gen electron 

        // HGCAL 3D cluster
        float closest_cl3d_dr = 9999.;
        int closest_cl3d = 0;
        int cl3d_Count = 0;

        for(int i=0; i < cl3d_N_;i++){

            if(cl3d_egid->at(i) != 1) continue;

            float dPhi = deltaPhi(propgenElPartPhi->at(0), cl3d_phi->at(i));

            float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-cl3d_eta->at(i),2));
            if(cl3d_pt->at(i) < 10) continue;
            cl3d_Count++;
            if(current_dr < closest_cl3d_dr){
                closest_cl3d_dr = current_dr;
                closest_cl3d = i;
            }
        }// end of loop to find the closest egamma to gen electron 

        // find egamma objects passing pixtrk signal windows
        if((closest_cl3d_dr < dr_cut && closest_cl3d_dr != 9999.)|| (closest_dr < dr_cut && closest_dr != 9999.)){

            if( closest_dr < closest_cl3d_dr ){
                EgEt =egCrysClusterEt ->at(closest_eg);
                EgEta=egCrysClusterEta->at(closest_eg);
                EgPhi=egCrysClusterPhi->at(closest_eg);

                EgGx = egCrysClusterGx->at(closest_eg);
                EgGy = egCrysClusterGy->at(closest_eg);
                EgGz = egCrysClusterGz->at(closest_eg);
                emvector.SetXYZ(EgGx,EgGy,EgGz);

            }
            else{

                EgEt =cl3d_pt->at(closest_cl3d);
                EgEta=cl3d_eta->at(closest_cl3d);
                EgPhi=cl3d_phi->at(closest_cl3d);

                EgGx = cl3d_x->at(closest_cl3d);
                EgGy = cl3d_y->at(closest_cl3d);
                EgGz = (float)cl3d_z->at(closest_cl3d);
                emvector.SetXYZ(EgGx,EgGy,EgGz);
            }

            eta_region = 0; // intialize
            if( fabs(EgEta) <= 0.8 ) eta_region =1;
            if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
            if( fabs(EgEta) <= 1.7 && fabs(EgEta) > 1.4 ) eta_region =3;
            if( fabs(EgEta) <= 2.1 && fabs(EgEta) > 1.7 ) eta_region =4;
            if( fabs(EgEta) <= 2.7 && fabs(EgEta) > 2.1 ) eta_region =5;
            if( fabs(EgEta) <= 3.0 && fabs(EgEta) > 2.7 ) eta_region =6;
            if( fabs(EgEta) > 3.0 ) continue;
            if( eta_region != 1 ) continue;

            regionFlag = eta_region;
            StorePixelHit(eta_region);

            // Fill # of total pixel clusters
            aPix1->Fill(first_hits.size());
            aPix2->Fill(second_hits.size());
            aPix3->Fill(third_hits.size());
            aPix4->Fill(fourth_hits.size());
             
            for(Int_t  i = 0; i < first_hits.size(); i++ )
            {
                TVector3 pixel;
                pixel.SetXYZ( first_hits[i].X() - simVx->at(1), first_hits[i].Y() - simVy->at(1), first_hits[i].Z() - simVz->at(1) );
                Float_t pixPhi = pixel.Phi();
                Float_t pixEta = pixel.Eta();
                Float_t dphi = deltaPhi(EgPhi, pixPhi);
                Float_t deta = EgEta - pixEta;
                Float_t dR = radius(dphi, deta);

                if( fabs(dphi) < 0.3 ) TrkIsodP[0]++;
                if( fabs(dphi) < 0.1 ) L1PixdP[0]++;
                if( dR < 0.3 ) TrkIsodR[0]++;
                if( dR < 0.1 ) L1PixdR[0]++;
            }
            
            for(Int_t  i = 0; i < second_hits.size(); i++ )
            {
                TVector3 pixel;
                pixel.SetXYZ( second_hits[i].X() - simVx->at(1), second_hits[i].Y() - simVy->at(1), second_hits[i].Z() - simVz->at(1) );
                Float_t pixPhi = pixel.Phi();
                Float_t pixEta = pixel.Eta();
                Float_t dphi = deltaPhi(EgPhi, pixPhi);
                Float_t deta = EgEta - pixEta;
                Float_t dR = radius(dphi, deta);

                if( fabs(dphi) < 0.3 ) TrkIsodP[1]++;
                if( fabs(dphi) < 0.1 ) L1PixdP[1]++;
                if( dR < 0.3 ) TrkIsodR[1]++;
                if( dR < 0.1 ) L1PixdR[1]++;
            }
            
            for(Int_t  i = 0; i < third_hits.size(); i++ )
            {
                TVector3 pixel;
                pixel.SetXYZ( third_hits[i].X() - simVx->at(1), third_hits[i].Y() - simVy->at(1), third_hits[i].Z() - simVz->at(1) );
                Float_t pixPhi = pixel.Phi();
                Float_t pixEta = pixel.Eta();
                Float_t dphi = deltaPhi(EgPhi, pixPhi);
                Float_t deta = EgEta - pixEta;
                Float_t dR = radius(dphi, deta);

                if( fabs(dphi) < 0.3 ) TrkIsodP[2]++;
                if( fabs(dphi) < 0.1 ) L1PixdP[2]++;
                if( dR < 0.3 ) TrkIsodR[2]++;
                if( dR < 0.1 ) L1PixdR[2]++;
            }

            for(Int_t  i = 0; i < fourth_hits.size(); i++ )
            {
                TVector3 pixel;
                pixel.SetXYZ( fourth_hits[i].X() - simVx->at(1), fourth_hits[i].Y() - simVy->at(1), fourth_hits[i].Z() - simVz->at(1) );
                Float_t pixPhi = pixel.Phi();
                Float_t pixEta = pixel.Eta();
                Float_t dphi = deltaPhi(EgPhi, pixPhi);
                Float_t deta = EgEta - pixEta;
                Float_t dR = radius(dphi, deta);

                if( fabs(dphi) < 0.3 ) TrkIsodP[3]++;
                if( fabs(dphi) < 0.1 ) L1PixdP[3]++;
                if( dR < 0.3 ) TrkIsodR[3]++;
                if( dR < 0.1 ) L1PixdR[3]++;
            }


            L1PixdP1->Fill(L1PixdP[0]);   
            L1PixdP2->Fill(L1PixdP[1]);
            L1PixdP3->Fill(L1PixdP[2]);
            L1PixdP4->Fill(L1PixdP[3]);

            L1PixdR1->Fill(L1PixdR[0]);
            L1PixdR2->Fill(L1PixdR[1]);
            L1PixdR3->Fill(L1PixdR[2]);
            L1PixdR4->Fill(L1PixdR[3]);

            PixTrkdP1->Fill(TrkIsodP[0]);  
            PixTrkdP2->Fill(TrkIsodP[1]);
            PixTrkdP3->Fill(TrkIsodP[2]);
            PixTrkdP4->Fill(TrkIsodP[3]);

            PixTrkdR1->Fill(TrkIsodR[0]);
            PixTrkdR2->Fill(TrkIsodR[1]);
            PixTrkdR3->Fill(TrkIsodR[2]);
            PixTrkdR4->Fill(TrkIsodR[3]);

            /*
            cout << endl;
            cout << "# of pixel clusters w/o any restriction" << endl;
            cout << " 1st layer: " << AllPixCount[0] << endl;
            cout << " 2nd layer: " << AllPixCount[1] << endl;
            cout << " 3rd layer: " << AllPixCount[2] << endl;
            cout << " 4th layer: " << AllPixCount[3] << endl;
            cout << endl;
            cout << "=====================================" << endl;
            cout << endl;
            cout << "# of pixel cluster in fabs(dphi) < 0.3 (for pixel track isolation)" << endl;
            cout << " 1st layer: " << TrkIsoPixCountdP[0] << endl;
            cout << " 2nd layer: " << TrkIsoPixCountdP[1] << endl;
            cout << " 3rd layer: " << TrkIsoPixCountdP[2] << endl;
            cout << " 4th layer: " << TrkIsoPixCountdP[3] << endl;
            cout << endl;
            cout << "=====================================" << endl;
            cout << endl;
            cout << "# of pixel cluster in delta R < 0.3 (for pixel track isolation)" << endl;
            cout << " 1st layer: " << TrkIsoPixCountdR[0] << endl;
            cout << " 2nd layer: " << TrkIsoPixCountdR[1] << endl;
            cout << " 3rd layer: " << TrkIsoPixCountdR[2] << endl;
            cout << " 4th layer: " << TrkIsoPixCountdR[3] << endl;
            cout << endl;
            cout << "=====================================" << endl;
            cout << endl;
            cout << "# of pixel cluster in fabs(dphi) < 0.1 (for level-1 pixel trigger)" << endl;
            cout << " 1st layer: " << L1PixTrkCountdP[0] << endl;
            cout << " 2nd layer: " << L1PixTrkCountdP[1] << endl;
            cout << " 3rd layer: " << L1PixTrkCountdP[2] << endl;
            cout << " 4th layer: " << L1PixTrkCountdP[3] << endl;
            cout << endl;
            cout << "=====================================" << endl;
            cout << endl;
            cout << "# of pixel cluster in delta R < 0.1 (for level-1 pixel trigger)" << endl;
            cout << " 1st layer: " << L1PixTrkCountdR[0] << endl;
            cout << " 2nd layer: " << L1PixTrkCountdR[1] << endl;
            cout << " 3rd layer: " << L1PixTrkCountdR[2] << endl;
            cout << " 4th layer: " << L1PixTrkCountdR[3] << endl;
            cout << endl;
            */
        } // if function

        //pixtrk_tree->Fill();
    } // end of entries loop 

    if( regionFlag == 1 )
    {
        aPix1->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix2->SetTitle("# of total pixel clusters on the 2nd pixel barrel");
        aPix3->SetTitle("# of total pixel clusters on the 3rd pixel barrel");
        aPix4->SetTitle("# of total pixel clusters on the 4th pixel barrel");
        
        L1PixdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        L1PixdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel barrel");
        L1PixdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel barrel");
        L1PixdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 4th pixel barrel");

        L1PixdR1->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        L1PixdR2->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel barrel");
        L1PixdR3->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel barrel");
        L1PixdR4->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 4th pixel barrel");

        PixTrkdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        PixTrkdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel barrel");
        PixTrkdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel barrel");
        PixTrkdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 4th pixel barrel");

        PixTrkdR1->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        PixTrkdR2->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel barrel");
        PixTrkdR3->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel barrel");
        PixTrkdR4->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 4th pixel barrel");
    }

    else if( regionFlag == 2 )
    {
        aPix1->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix2->SetTitle("# of total pixel clusters on the 2nd pixel barrel");
        aPix3->SetTitle("# of total pixel clusters on the 3rd pixel barrel");
        aPix4->SetTitle("# of total pixel clusters on the 1st pixel disk");

        L1PixdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        L1PixdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel barrel");
        L1PixdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel barrel");
        L1PixdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");

        L1PixdR1->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        L1PixdR2->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel barrel");
        L1PixdR3->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel barrel");
        L1PixdR4->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");

        PixTrkdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        PixTrkdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel barrel");
        PixTrkdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel barrel");
        PixTrkdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");

        PixTrkdR1->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        PixTrkdR2->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel barrel");
        PixTrkdR3->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel barrel");
        PixTrkdR4->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 4th pixel disk");
    }

    else if( regionFlag == 3 )
    {
        aPix1->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix2->SetTitle("# of total pixel clusters on the 2nd pixel barrel");
        aPix3->SetTitle("# of total pixel clusters on the 1st pixel disk");
        aPix4->SetTitle("# of total pixel clusters on the 2nd pixel disk");

        L1PixdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        L1PixdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel barrel");
        L1PixdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        L1PixdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");

        L1PixdR1->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        L1PixdR2->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel barrel");
        L1PixdR3->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        L1PixdR4->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");

        PixTrkdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        PixTrkdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel barrel");
        PixTrkdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");
        PixTrkdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");

        PixTrkdR1->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        PixTrkdR2->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel barrel");
        PixTrkdR3->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
        PixTrkdR4->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
    }

    else if( regionFlag == 4 )
    {
        aPix1->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix2->SetTitle("# of total pixel clusters on the 1st pixel disk");
        aPix3->SetTitle("# of total pixel clusters on the 2nd pixel disk");
        aPix4->SetTitle("# of total pixel clusters on the 3rd pixel disk");

        L1PixdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        L1PixdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        L1PixdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");
        L1PixdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel disk");

        L1PixdR1->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        L1PixdR2->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        L1PixdR3->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");
        L1PixdR4->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel disk");

        PixTrkdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        PixTrkdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");
        PixTrkdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");
        PixTrkdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel disk");

        PixTrkdR1->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        PixTrkdR2->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
        PixTrkdR3->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
        PixTrkdR4->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel disk");
    }

    else if( regionFlag == 5 )
    {
        aPix1->SetTitle("# of total pixel clusters on the 1st pixel disk");
        aPix2->SetTitle("# of total pixel clusters on the 2nd pixel disk");
        aPix3->SetTitle("# of total pixel clusters on the 3rd pixel disk");
        aPix4->SetTitle("# of total pixel clusters on the 4th pixel disk");

        L1PixdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        L1PixdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");
        L1PixdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel disk");
        L1PixdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 4th pixel disk");

        L1PixdR1->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        L1PixdR2->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");
        L1PixdR3->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel disk");
        L1PixdR4->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 4th pixel disk");

        PixTrkdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");
        PixTrkdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");
        PixTrkdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel disk");
        PixTrkdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 4th pixel disk");

        PixTrkdR1->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
        PixTrkdR2->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
        PixTrkdR3->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel disk");
        PixTrkdR4->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 4th pixel disk");
    }

    else if( regionFlag == 6 )
    {
        aPix1->SetTitle("# of total pixel clusters on the 1st pixel disk");
        aPix2->SetTitle("# of total pixel clusters on the 2nd pixel disk");
        aPix3->SetTitle("# of total pixel clusters on the 3rd pixel disk");
        aPix4->SetTitle("# of total pixel clusters on the 4th pixel disk");

        L1PixdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        L1PixdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");
        L1PixdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel disk");
        L1PixdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 4th pixel disk");

        L1PixdR1->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        L1PixdR2->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");
        L1PixdR3->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel disk");
        L1PixdR4->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 4th pixel disk");

        PixTrkdP1->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");
        PixTrkdP2->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");
        PixTrkdP3->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel disk");
        PixTrkdP4->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 4th pixel disk");

        PixTrkdR1->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
        PixTrkdR2->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
        PixTrkdR3->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel disk");
        PixTrkdR4->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 4th pixel disk");
    }
    
    aPix1->Write();
    aPix2->Write();
    aPix3->Write();
    aPix4->Write();

    L1PixdP1->Write();   
    L1PixdP2->Write();
    L1PixdP3->Write();
    L1PixdP4->Write();

    L1PixdR1->Write();
    L1PixdR2->Write();
    L1PixdR3->Write();
    L1PixdR4->Write();

    PixTrkdP1->Write();
    PixTrkdP2->Write();
    PixTrkdP3->Write();
    PixTrkdP4->Write();

    PixTrkdR1->Write();
    PixTrkdR2->Write();
    PixTrkdR3->Write();
    PixTrkdR4->Write();

    outfile->Close();
    //outfile->Write();
}

