#define eff_plot_cxx
#include "eff_plot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void eff_plot::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   int nbins = 60; float x1 = -3.0; float x2 = 3.0;
   
   TH1F* hEG_denom = new TH1F("hEG_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hEG_nom = new TH1F("hEG_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hEG_nom->Sumw2();
   hEG_denom->Sumw2();

   TH1F* hPix_denom = new TH1F("hPix_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hPix_nom = new TH1F("hPix_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hPix_nom->Sumw2();
   hPix_denom->Sumw2();
   
   TH1F* hPix_iso_denom = new TH1F("hPix_iso_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hPix_iso_nom = new TH1F("hPix_iso_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hPix_iso_nom->Sumw2();
   hPix_iso_denom->Sumw2();
   
   int width_bit[27] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000, 0x10000, 0x20000, 0x40000, 0x80000,
                        0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, 0x4000000};
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if(nt_genPt < 20.) continue;
      //if( fabs(nt_genEta) > 0.8 ) continue;
   
      int eg_size = ntEgEt->size();
      int gen_matched_eg = -1;
      float dr = 999.;
   
      // find the closest L1 egamma to gen electron
      for(int i = 0; i < eg_size; i++){
         float temp_dr = sqrt(pow(nt_genPhi-ntEgPhi->at(i),2)+pow(nt_genEta-ntEgEta->at(i),2));
         if( temp_dr < dr){
           dr = temp_dr;
           gen_matched_eg = i;
         }
      }// eg loop
  
      float pt_err = 0.; 
      if(gen_matched_eg != -1) {pt_err = fabs(nt_genPt-ntEgEt->at(gen_matched_eg))/nt_genPt; }
      if( gen_matched_eg != -1 && pt_err < 0.5){
        if((trigger_bit_width->at(gen_matched_eg)&width_bit[0])==width_bit[0]) hPix_nom->Fill(nt_genEta, 1.);
        if((trigger_bit_width_iso->at(gen_matched_eg)&width_bit[0])==width_bit[0]) hPix_iso_nom->Fill(nt_genEta, 1.);
        hEG_nom->Fill(nt_genEta, 1.);
      }

      hEG_denom->Fill(nt_genEta, 1.);
      hPix_denom->Fill(nt_genEta, 1.);
      hPix_iso_denom->Fill(nt_genEta, 1.);
   } // event loop
   
   TGraphAsymmErrors* hEG = new TGraphAsymmErrors(hEG_nom, hEG_denom,"B");
   TGraphAsymmErrors* hPix = new TGraphAsymmErrors(hPix_nom, hPix_denom,"B");
   TGraphAsymmErrors* hPixIso = new TGraphAsymmErrors(hPix_iso_nom, hPix_iso_denom,"B");
   
   TCanvas *c1 = new TCanvas("c1","c1",800,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(1); // axis width, default is 1
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.12);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.15);
   c1->SetGrid();
   c1->SetTicky(1);
   c1->SetTickx(1);
   c1->cd();

   hPix->GetXaxis()->SetTitle("#eta_{gen} ");
   hPix->GetXaxis()->SetTitleOffset(1.);
   hPix->GetXaxis()->SetTitleSize(0.055);
   hPix->GetXaxis()->SetNdivisions(510);
   hPix->GetYaxis()->SetNdivisions(506);
   hPix->GetXaxis()->SetLabelSize(0.05);
   hPix->GetYaxis()->SetLabelSize(0.05);
   hPix->GetXaxis()->SetRangeUser(-3.0, 3.0);
   hPix->GetYaxis()->SetRangeUser(0.5, 1.1);
   hPix->GetYaxis()->SetTitle("Efficiency");
   hPix->GetYaxis()->SetTitleOffset(1.2);
   hPix->GetYaxis()->SetTitleSize(0.055);

   hPix->SetMarkerColor(4);
   hPix->SetLineColor(4);
   hPix->SetLineWidth(1);
   hPix->SetMarkerStyle(29);
   hPix->SetMarkerSize(1.8);
   hPix->Draw("ape");

   hEG->SetMarkerColor(2);
   hEG->SetLineColor(2);
   hEG->SetLineWidth(1);
   hEG->SetMarkerStyle(20);
   hEG->SetMarkerSize(1.);
   hEG->Draw("pe same");
   
   hPixIso->SetMarkerColor(8);
   hPixIso->SetLineColor(8);
   hPixIso->SetLineWidth(1);
   hPixIso->SetMarkerStyle(29);
   hPixIso->SetMarkerSize(1.8);
   hPixIso->Draw("pe same");

   TLegend *Lgd = new TLegend(0.25, 0.9, 0.95, 0.95);

   Lgd-> SetNColumns(4);
   Lgd->SetFillColor(0);
   Lgd->SetTextFont(42);
   Lgd->SetTextSize(0.025);
   Lgd->SetBorderSize(0);
   Lgd->SetFillStyle(0);
   Lgd->AddEntry(hEG,"Phase-2 L1 EG","lp");
   Lgd->AddEntry(hPix,"Pixel matching","lp");
   Lgd->AddEntry(hPixIso,"Pixel matching + Isolation","lp");
   Lgd->Draw();

   TLatex t(-1.8,1.11,"CMS Preliminary Simulation, Phase 2, <PU>=200");
   t.SetTextSize(0.035);
   t.Draw();

   TLatex pt_cut(2.,0.2,"p_{T}^{gen} > 20 GeV");
   pt_cut.SetTextSize(0.035);
   pt_cut.Draw();

   c1->Print("Eff-iso-v2.png");
   c1->Print("Eff-iso-v2.pdf");
}
