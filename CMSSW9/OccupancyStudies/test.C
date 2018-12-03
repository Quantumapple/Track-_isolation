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
    
    TH1F *aPix[4];
    aPix[0] = new TH1F("allPix1","",400,0,8000);
    aPix[1] = new TH1F("allPix2","",400,0,8000);
    aPix[2] = new TH1F("allPix3","",400,0,8000);
    aPix[3] = new TH1F("allPix4","",400,0,8000);
   
    TH1F *dPHI03[4];
    dPHI03[0] = new TH1F("dPhiFir1","",150,0,1500);
    dPHI03[1] = new TH1F("dPhiFir2","",150,0,1500);
    dPHI03[2] = new TH1F("dPhiFir3","",150,0,1500);
    dPHI03[3] = new TH1F("dPhiFir4","",150,0,1500);
    
    TH1F *dPHI02[4];
    dPHI02[0] = new TH1F("dPhiSec1","",250,0,500);
    dPHI02[1] = new TH1F("dPhiSec2","",250,0,500);
    dPHI02[2] = new TH1F("dPhiSec3","",250,0,500);
    dPHI02[3] = new TH1F("dPhiSec4","",250,0,500);

    TH1F *dPHI01[4];
    dPHI01[0] = new TH1F("dPhiThi1","",125,0,250);
    dPHI01[1] = new TH1F("dPhiThi2","",125,0,250);
    dPHI01[2] = new TH1F("dPhiThi3","",125,0,250);
    dPHI01[3] = new TH1F("dPhiThi4","",125,0,250);
    
    TH1F *DR03[4];
    DR03[0] = new TH1F("dRFir1","",50,0,250);
    DR03[1] = new TH1F("dRFir2","",50,0,250);
    DR03[2] = new TH1F("dRFir3","",50,0,250);
    DR03[3] = new TH1F("dRFir4","",50,0,250);
    
    TH1F *DR02[4];
    DR02[0] = new TH1F("dRSec1","",75,0,125);
    DR02[1] = new TH1F("dRSec2","",75,0,125);
    DR02[2] = new TH1F("dRSec3","",75,0,125);
    DR02[3] = new TH1F("dRSec4","",75,0,125);

    TH1F *DR01[4];
    DR01[0] = new TH1F("dRThi1","",25,0,50);
    DR01[1] = new TH1F("dRThi2","",25,0,50);
    DR01[2] = new TH1F("dRThi3","",25,0,50);
    DR01[3] = new TH1F("dRThi4","",25,0,50);
    
    int RegionFlag = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
    //for (Long64_t jentry=0; jentry<500;jentry++) { //nentries
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;

        first_hits.clear();
        second_hits.clear();
        third_hits.clear();
        fourth_hits.clear();

        //find closest egamma object to the gen electron
        EgN=egCrysClusterEt->size();
        int cl3d_N_ = cl3d_pt->size();

        float closest_dr = 9999.;
        int closest_eg = 0;

        //Int_t TrkIsodP[4] = {};
        //Int_t TrkIsodR[4] = {};
        //Int_t L1PixdP[4] = {};
        //Int_t L1PixdR[4] = {};
        
        Int_t dphi01[4] = {};
        Int_t dphi02[4] = {};
        Int_t dphi03[4] = {};

        Int_t dR01[4] = {};
        Int_t dR02[4] = {};
        Int_t dR03[4] = {};

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
            if( eta_region != 6 ) continue;

            RegionFlag = eta_region;
            StorePixelHit(eta_region);

            // Fill # of total pixel clusters
            aPix[0]->Fill(first_hits.size());
            aPix[1]->Fill(second_hits.size());
            aPix[2]->Fill(third_hits.size());
            aPix[3]->Fill(fourth_hits.size());
             
            for(Int_t  i = 0; i < first_hits.size(); i++ )
            {
                TVector3 pixel;
                pixel.SetXYZ( first_hits[i].X() - simVx->at(1), first_hits[i].Y() - simVy->at(1), first_hits[i].Z() - simVz->at(1) );
                Float_t pixPhi = pixel.Phi();
                Float_t pixEta = pixel.Eta();
                Float_t dphi = deltaPhi(EgPhi, pixPhi);
                Float_t deta = EgEta - pixEta;
                Float_t dR = radius(dphi, deta);

                if( fabs(dphi) > 0.3 ) continue;

                if( fabs(dphi) < 0.3 ) dphi03[0]++; 
                if( fabs(dphi) < 0.2 ) dphi02[0]++;
                if( fabs(dphi) < 0.1 ) dphi01[0]++;
                if( dR < 0.3 ) dR03[0]++; 
                if( dR < 0.2 ) dR02[0]++;
                if( dR < 0.1 ) dR01[0]++;
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

                if( fabs(dphi) > 0.3 ) continue;
                
                if( fabs(dphi) < 0.3 ) dphi03[1]++; 
                if( fabs(dphi) < 0.2 ) dphi02[1]++;
                if( fabs(dphi) < 0.1 ) dphi01[1]++;
                if( dR < 0.3 ) dR03[1]++; 
                if( dR < 0.2 ) dR02[1]++;
                if( dR < 0.1 ) dR01[1]++;
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

                if( fabs(dphi) > 0.3 ) continue;
                
                if( fabs(dphi) < 0.3 ) dphi03[2]++; 
                if( fabs(dphi) < 0.2 ) dphi02[2]++;
                if( fabs(dphi) < 0.1 ) dphi01[2]++;
                if( dR < 0.3 ) dR03[2]++; 
                if( dR < 0.2 ) dR02[2]++;
                if( dR < 0.1 ) dR01[2]++;
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

                if( fabs(dphi) > 0.3 ) continue;
                
                if( fabs(dphi) < 0.3 ) dphi03[3]++; 
                if( fabs(dphi) < 0.2 ) dphi02[3]++;
                if( fabs(dphi) < 0.1 ) dphi01[3]++;
                if( dR < 0.3 ) dR03[3]++; 
                if( dR < 0.2 ) dR02[3]++;
                if( dR < 0.1 ) dR01[3]++;
            }
            
            dPHI01[0]->Fill(dphi01[0]); dPHI01[1]->Fill(dphi01[1]); dPHI01[2]->Fill(dphi01[2]); dPHI01[3]->Fill(dphi01[3]);
            dPHI02[0]->Fill(dphi02[0]); dPHI02[1]->Fill(dphi02[1]); dPHI02[2]->Fill(dphi02[2]); dPHI02[3]->Fill(dphi02[3]);
            dPHI03[0]->Fill(dphi03[0]); dPHI03[1]->Fill(dphi03[1]); dPHI03[2]->Fill(dphi03[2]); dPHI03[3]->Fill(dphi03[3]);
            
            DR01[0]->Fill(dR01[0]); DR01[1]->Fill(dR01[1]); DR01[2]->Fill(dR01[2]); DR01[3]->Fill(dR01[3]);
            DR02[0]->Fill(dR02[0]); DR02[1]->Fill(dR02[1]); DR02[2]->Fill(dR02[2]); DR02[3]->Fill(dR02[3]);
            DR03[0]->Fill(dR03[0]); DR03[1]->Fill(dR03[1]); DR03[2]->Fill(dR03[2]); DR03[3]->Fill(dR03[3]);


        } // if function

    } // end of entries loop 

    if( RegionFlag == 1 )
    {
        aPix[0]->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix[1]->SetTitle("# of total pixel clusters on the 2nd pixel barrel");
        aPix[2]->SetTitle("# of total pixel clusters on the 3rd pixel barrel");
        aPix[3]->SetTitle("# of total pixel clusters on the 4th pixel barrel");
        
        dPHI01[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        dPHI01[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel barrel");
        dPHI01[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel barrel");
        dPHI01[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 4th pixel barrel");
        
        dPHI02[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel barrel");   
        dPHI02[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 2nd pixel barrel");
        dPHI02[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 3rd pixel barrel");
        dPHI02[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 4th pixel barrel");
        
        dPHI03[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        dPHI03[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel barrel");
        dPHI03[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel barrel");
        dPHI03[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 4th pixel barrel");

        DR01[0]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        DR01[1]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel barrel");
        DR01[2]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel barrel");
        DR01[3]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 4th pixel barrel");
        
        DR02[0]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel barrel");
        DR02[1]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 2nd pixel barrel");
        DR02[2]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 3rd pixel barrel");
        DR02[3]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 4th pixel barrel");
        
        DR03[0]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        DR03[1]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel barrel");
        DR03[2]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel barrel");
        DR03[3]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 4th pixel barrel");
    }

    else if( RegionFlag == 2 )
    {
        aPix[0]->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix[1]->SetTitle("# of total pixel clusters on the 2nd pixel barrel");
        aPix[2]->SetTitle("# of total pixel clusters on the 3rd pixel barrel");
        aPix[3]->SetTitle("# of total pixel clusters on the 1st pixel disk");
        
        dPHI01[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        dPHI01[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel barrel");
        dPHI01[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel barrel");
        dPHI01[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        
        dPHI02[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel barrel");   
        dPHI02[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 2nd pixel barrel");
        dPHI02[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 3rd pixel barrel");
        dPHI02[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel disk");
        
        dPHI03[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        dPHI03[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel barrel");
        dPHI03[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel barrel");
        dPHI03[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");

        DR01[0]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        DR01[1]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel barrel");
        DR01[2]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel barrel");
        DR01[3]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        
        DR02[0]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel barrel");
        DR02[1]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 2nd pixel barrel");
        DR02[2]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 3rd pixel barrel");
        DR02[3]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel disk");
        
        DR03[0]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        DR03[1]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel barrel");
        DR03[2]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel barrel");
        DR03[3]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
    }

    else if( RegionFlag == 3 )
    {
        aPix[0]->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix[1]->SetTitle("# of total pixel clusters on the 2nd pixel barrel");
        aPix[2]->SetTitle("# of total pixel clusters on the 1st pixel disk");
        aPix[3]->SetTitle("# of total pixel clusters on the 2nd pixel disk");
        
        dPHI01[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        dPHI01[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel barrel");
        dPHI01[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        dPHI01[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");
        
        dPHI02[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel barrel");   
        dPHI02[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 2nd pixel barrel");
        dPHI02[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel disk");
        dPHI02[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 2nd pixel disk");
        
        dPHI03[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        dPHI03[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel barrel");
        dPHI03[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");
        dPHI03[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");

        DR01[0]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        DR01[1]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel barrel");
        DR01[2]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        DR01[3]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");
        
        DR02[0]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel barrel");
        DR02[1]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 2nd pixel barrel");
        DR02[2]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel disk");
        DR02[3]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 2nd pixel disk");
        
        DR03[0]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        DR03[1]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel barrel");
        DR03[2]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
        DR03[3]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
    }

    else if( RegionFlag == 4 )
    {
        aPix[0]->SetTitle("# of total pixel clusters on the 1st pixel barrel");
        aPix[1]->SetTitle("# of total pixel clusters on the 1st pixel disk");
        aPix[2]->SetTitle("# of total pixel clusters on the 2nd pixel disk");
        aPix[3]->SetTitle("# of total pixel clusters on the 3rd pixel disk");
        
        dPHI01[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel barrel");   
        dPHI01[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        dPHI01[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");
        dPHI01[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel disk");
        
        dPHI02[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel barrel");   
        dPHI02[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel disk");
        dPHI02[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 2nd pixel disk");
        dPHI02[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 3rd pixel disk");
        
        dPHI03[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel barrel");
        dPHI03[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");
        dPHI03[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");
        dPHI03[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel disk");

        DR01[0]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel barrel");
        DR01[1]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        DR01[2]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");
        DR01[3]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel disk");
        
        DR02[0]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel barrel");
        DR02[1]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel disk");
        DR02[2]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 2nd pixel disk");
        DR02[3]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 3rd pixel disk");
        
        DR03[0]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel barrel");
        DR03[1]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
        DR03[2]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
        DR03[3]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel disk");
    }

    else if( RegionFlag == 5 )
    {
        aPix[0]->SetTitle("# of total pixel clusters on the 1st pixel disk");
        aPix[1]->SetTitle("# of total pixel clusters on the 2nd pixel disk");
        aPix[2]->SetTitle("# of total pixel clusters on the 3rd pixel disk");
        aPix[3]->SetTitle("# of total pixel clusters on the 4th pixel disk");
        
        dPHI01[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 1st pixel disk");
        dPHI01[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");
        dPHI01[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel disk");
        dPHI01[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 4th pixel disk");
        
        dPHI02[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 1st pixel disk");
        dPHI02[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 2nd pixel disk");
        dPHI02[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 3rd pixel disk");
        dPHI02[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 4th pixel disk");
        
        dPHI03[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 1st pixel disk");
        dPHI03[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");
        dPHI03[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel disk");
        dPHI03[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 4th pixel disk");

        DR01[0]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 1st pixel disk");
        DR01[1]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");
        DR01[2]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel disk");
        DR01[3]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 4th pixel disk");
        
        DR02[0]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 1st pixel disk");
        DR02[1]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 2nd pixel disk");
        DR02[2]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 3rd pixel disk");
        DR02[3]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 4th pixel disk");
        
        DR03[0]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 1st pixel disk");
        DR03[1]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
        DR03[2]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel disk");
        DR03[3]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 4th pixel disk");
    }

    else if( RegionFlag == 6 )
    {
        aPix[0]->SetTitle("# of total pixel clusters on the 2nd pixel disk");
        aPix[1]->SetTitle("# of total pixel clusters on the 3rd pixel disk");
        aPix[2]->SetTitle("# of total pixel clusters on the 4th pixel disk");
        aPix[3]->SetTitle("# of total pixel clusters on the 5th pixel disk");
        
        dPHI01[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 2nd pixel disk");
        dPHI01[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 3rd pixel disk");
        dPHI01[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 4th pixel disk");
        dPHI01[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.1 on the 5th pixel disk");
        
        dPHI02[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 2nd pixel disk");
        dPHI02[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 3rd pixel disk");
        dPHI02[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 4th pixel disk");
        dPHI02[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.2 on the 5th pixel disk");
        
        dPHI03[0]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 2nd pixel disk");
        dPHI03[1]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 3rd pixel disk");
        dPHI03[2]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 4th pixel disk");
        dPHI03[3]->SetTitle("# of pixel clusters in |#Delta#phi| < 0.3 on the 5th pixel disk");

        DR01[0]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 2nd pixel disk");
        DR01[1]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 3rd pixel disk");
        DR01[2]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 4th pixel disk");
        DR01[3]->SetTitle("# of pixel clusters in #DeltaR < 0.1 on the 5th pixel disk");
                                                                          
        DR02[0]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 2nd pixel disk");
        DR02[1]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 3rd pixel disk");
        DR02[2]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 4th pixel disk");
        DR02[3]->SetTitle("# of pixel clusters in #DeltaR < 0.2 on the 5th pixel disk");
                                                                          
        DR03[0]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 2nd pixel disk");
        DR03[1]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 3rd pixel disk");
        DR03[2]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 4th pixel disk");
        DR03[3]->SetTitle("# of pixel clusters in #DeltaR < 0.3 on the 5th pixel disk");
    }
    
    aPix[0]->Write();
    aPix[1]->Write();
    aPix[2]->Write();
    aPix[3]->Write();

    dPHI01[0]->Write();   
    dPHI01[1]->Write();
    dPHI01[2]->Write();
    dPHI01[3]->Write();

    dPHI02[0]->Write();   
    dPHI02[1]->Write();
    dPHI02[2]->Write();
    dPHI02[3]->Write();
    
    dPHI03[0]->Write();   
    dPHI03[1]->Write();
    dPHI03[2]->Write();
    dPHI03[3]->Write();
    
    DR01[0]->Write();   
    DR01[1]->Write();
    DR01[2]->Write();
    DR01[3]->Write();

    DR02[0]->Write();   
    DR02[1]->Write();
    DR02[2]->Write();
    DR02[3]->Write();
    
    DR03[0]->Write();   
    DR03[1]->Write();
    DR03[2]->Write();
    DR03[3]->Write();

    outfile->Close();
    //outfile->Write();
}

