#define l1rate_cxx
#include "l1rate.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void l1rate::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   float genN = 282700.0; 
   int eta_r = 0;

   int nbins = 50; float x1 = 5.0 ; float x2 = 55.0 ;
   TH1F* hEG = new TH1F("hEG",";E_{T} threshold (GeV); Rate (kHz)",nbins,x1,x2);
   TH1F* hPXEG = new TH1F("hPXEG",";E_{T} threshold (GeV); Rate (kHz)",nbins,x1,x2);
   TH1F* hPXEG_iso = new TH1F("hPXEG_iso",";E_{T} threshold (GeV); Rate (kHz)",nbins,x1,x2);

   TH1F *h1 = new TH1F("h1","h1",nbins,x1,x2);
   TH1F *h2 = new TH1F("h2","h2",nbins,x1,x2);
   TH1F *h3 = new TH1F("h3","h3",nbins,x1,x2);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      if(!ntnEg2) continue;

      float max_eget = -999.;
      for(int i=0; i < ntnEg2; i++){
         if(eta_r == 0 && ntEgEt->at(i) > max_eget) max_eget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_eget && (fabs(ntEgEta->at(i)) < 1.3) ) max_eget = ntEgEt->at(i);
           
      }// eg loop

      float max_pxeget = -999.;
      for(int i=0; i < ntnEg2; i++){
         if(eta_r == 0 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i)) max_pxeget = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i) && (fabs(ntEgEta->at(i)) < 1.3)) max_pxeget = ntEgEt->at(i);
      
      }// pxeg loop
      
      float max_pxeget_iso = -999.;
      for(int i=0; i < ntnEg2; i++){
         if(eta_r == 0 && ntEgEt->at(i) > max_pxeget_iso && ntCl_iso_match->at(i)) max_pxeget_iso = ntEgEt->at(i);
         //if(eta_r == 1 && ntEgEt->at(i) > max_pxeget && ntCl_match->at(i) && (fabs(ntEgEta->at(i)) < 1.3)) max_pxeget = ntEgEt->at(i);
      
      }// pxeg loop
   
      hEG->Fill(max_eget);
      if(max_pxeget > 0) hPXEG->Fill(max_pxeget);
      if(max_pxeget_iso > 0) hPXEG_iso->Fill(max_pxeget_iso);

   
   } // event loop

   h1->Add(hPXEG_iso);

   for (int i=0; i<= nbins; i++) {
       hEG->SetBinContent(i, hEG -> Integral(i, nbins+1) );
       hEG->SetBinError(i, sqrt( hEG -> GetBinContent(i) ) ); 

       hPXEG->SetBinContent(i, hPXEG -> Integral(i, nbins+1) );
       hPXEG->SetBinError(i, sqrt( hPXEG -> GetBinContent(i) ) ); 
       
       hPXEG_iso->SetBinContent(i, hPXEG_iso -> Integral(i, nbins+1) );
       hPXEG_iso->SetBinError(i, sqrt( hPXEG_iso -> GetBinContent(i) ) ); 
   }

   h2->Add(hPXEG_iso);
   h3->Add(hPXEG);

   hEG -> Scale(30000./genN);
   hPXEG -> Scale(30000./genN);
   hPXEG_iso -> Scale(30000./genN);

   TCanvas *c1 = new TCanvas("c1","c1",700,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(2); // axis width, default is 1
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.12);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.2);
   c1->SetLogy();
   c1->SetGrid();
   c1->SetTicky(1);
   c1->SetTickx(1);

   hEG->SetTitle("");
   hEG->GetXaxis()->CenterTitle(true);
   hEG->GetXaxis()->SetTitleOffset(1.);
   hEG->GetXaxis()->SetTitleSize(0.055);
   hEG->GetXaxis()->SetNdivisions(505);
   hEG->GetYaxis()->SetNdivisions(546);
   hEG->GetXaxis()->SetLabelSize(0.05);
   hEG->GetYaxis()->SetLabelSize(0.05);
   hEG->GetXaxis()->SetRangeUser(10.,55.);
   hEG->GetYaxis()->SetRangeUser(1,9000);
   hEG->GetYaxis()->SetTitleOffset(1.2);
   hEG->GetYaxis()->SetTitleSize(0.055);

   hEG->SetMarkerColor(2);
   hEG->SetLineColor(2);
   hEG->SetLineWidth(2);
   hEG->SetMarkerStyle(20);
   hEG->SetMarkerSize(1.);
   hEG->Draw("HIST e");

   hPXEG->SetMarkerColor(4);
   hPXEG->SetLineColor(4);
   hPXEG->SetLineWidth(2);
   hPXEG->SetMarkerStyle(29);
   hPXEG->SetMarkerSize(1.8);
   hPXEG->Draw("HIST same e");
   
   hPXEG_iso->SetMarkerColor(8);
   hPXEG_iso->SetLineColor(8);
   hPXEG_iso->SetLineWidth(2);
   hPXEG_iso->SetMarkerStyle(29);
   hPXEG_iso->SetMarkerSize(1.8);
   hPXEG_iso->Draw("HIST same e");
   
   TLegend *Lgd = new TLegend(0.47, 0.7, 0.93, 0.85);
   Lgd->SetFillColor(0);
   Lgd->SetTextFont(42);
   Lgd->SetTextSize(0.04);
   Lgd->SetBorderSize(1);
   Lgd->AddEntry(hEG,"L1 Electron-Gamma","lp");
   Lgd->AddEntry(hPXEG,"L1 Pixel Detector","lp");
   Lgd->AddEntry(hPXEG_iso,"PixTRK + TrkIso","lp");
   Lgd->Draw();

   TLatex t(12,9600,"CMS Preliminary Simulation, Phase 2, <PU>=200");
   t.SetTextSize(0.035);
   t.Draw();

   TString eta_str;
   if(eta_r == 0) eta_str = "#lbar#bf{#eta}#lbar < 0.8";

   TLatex eta_range(43,4500,eta_str);
   eta_range.SetTextSize(0.05);
   eta_range.Draw();

   TString l1rate = "l1rate";
   if(eta_r == 0) l1rate = l1rate + "_all.png";

   c1->Print(l1rate);

   TFile *output = new TFile("file.root","RECREATE");
   h1->Write();
   h2->Write();
   h3->Write();
   hEG->Write();
   hPXEG->Write();
   hPXEG_iso->Write();
   output->Close();

}
