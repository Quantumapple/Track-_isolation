#define sw_treeMaker_cxx
#include "sw_treeMaker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
#include <iostream>
#include <string>
using namespace std;


void sw_treeMaker::Loop(int nevent = 0, int debug = 0)
{
   if (fChain == 0) return;

   TFile *file = new TFile("ROI.root","recreate");

   TTree* out_tree = new TTree("t","t");

   float abs_pterr; out_tree->Branch("abs_pterr", &abs_pterr, "abs_pterr/F");
   float closest_egEt; out_tree->Branch("closest_egEt", &closest_egEt, "closest_egEt/F");
   float closest_egEta; out_tree->Branch("closest_egEta", &closest_egEta, "closest_egEta/F");
   float closest_egPhi; out_tree->Branch("closest_egPhi", &closest_egPhi, "closest_egPhi/F");
   float closest_eg_dr; out_tree->Branch("closest_eg_dr", &closest_eg_dr, "closest_egEt/F");

   vector<float>* bpix_l1_dphi_roi = new std::vector<float>; out_tree->Branch("bpix_l1_dphi_roi", &bpix_l1_dphi_roi);
   vector<float>* bpix_l2_dphi_roi = new std::vector<float>; out_tree->Branch("bpix_l2_dphi_roi", &bpix_l2_dphi_roi);
   vector<float>* bpix_l3_dphi_roi = new std::vector<float>; out_tree->Branch("bpix_l3_dphi_roi", &bpix_l3_dphi_roi);
   vector<float>* bpix_l4_dphi_roi = new std::vector<float>; out_tree->Branch("bpix_l4_dphi_roi", &bpix_l4_dphi_roi);

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   if(nevent) nentries = nevent;
   cout << " nentries: " <<  nentries << endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      bpix_l1_dphi_roi->clear();
      bpix_l2_dphi_roi->clear();
      bpix_l3_dphi_roi->clear();
      bpix_l4_dphi_roi->clear();

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //find closest egamma object to the gen electron
      int egN_ = egN;
      if(!egN_) continue;
      if(genPartPt->at(0) < 5) continue;

      float closest_dr = 9999.;
      int closest_eg = 0;
     
      int egCount = 0; 
      for(int i=0; i < egN_;i++){

         float dPhi = propgenPartPhi->at(0)-egPhi->at(i);
         if( dPhi > float(M_PI) ) dPhi -= float(2*M_PI);
         else if( dPhi <= -float(M_PI) ) dPhi += float(2*M_PI);
 
         float current_dr = sqrt(pow(dPhi,2)+pow(propgenPartEta->at(0)-egEta->at(i),2));
         if(egEt->at(i) < 5) continue;
         egCount++;
         if(current_dr < closest_dr){
           closest_dr = current_dr;
           closest_eg = i;
         } 
      }// end of loop to find the closest egamma to gen electron 
      if(egCount!=1) continue;
      if(debug) cout << jentry << "th event " << " closest eg dR: " << closest_dr << " gen pt: " << genPartPt->at(0) << " eg et: " << egEt->at(closest_eg) << endl;
      abs_pterr = fabs(genPartPt->at(0)-egEt->at(closest_eg))/genPartPt->at(0);
      closest_egEt = egEt->at(closest_eg);
      closest_egEta = egEta->at(closest_eg);
      closest_egPhi = egPhi->at(closest_eg);
      closest_eg_dr = closest_dr;

      file->cd();

      for(int a=0; a<bHitN; a++){
         TVector3 current_hit;
         current_hit.SetXYZ( bHitGx->at(a), bHitGy->at(a), bHitGz->at(a) );

         float temp_dphi = current_hit.Phi() - closest_egPhi;
         if( temp_dphi > float(M_PI) ) temp_dphi -= float(2*M_PI);
         else if( temp_dphi <= -float(M_PI) ) temp_dphi += float(2*M_PI);

         if(bHitLayer->at(a) == 1){
            bpix_l1_dphi_roi->push_back(temp_dphi);
         } 
         if(bHitLayer->at(a) == 2){
            bpix_l2_dphi_roi->push_back(temp_dphi);
         } 
         if(bHitLayer->at(a) == 3){
            bpix_l3_dphi_roi->push_back(temp_dphi);
         } 
         if(bHitLayer->at(a) == 4){
            bpix_l4_dphi_roi->push_back(temp_dphi);
         } 
      }
       
      out_tree->Fill(); 
 
       

   } // end of event loop

   file->Write();
}

