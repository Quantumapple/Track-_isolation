#define roofit_input_cxx
#include "roofit_input.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TLatex.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TPaveStats.h>
#include <TPaletteAxis.h>
#include <TGraphErrors.h>

void roofit_input::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //TH2F* pmDPhi1_dist = new TH2F("","", 90,10,100,100,-0.15,0.15);
   TH2F* pmDPhi1_dist = new TH2F("","", 90,10,100,100,-0.2,0.2);
   const int binSize = 90;
   float x[binSize], median[binSize];
   float x_err[binSize], median_err[binSize];

   //int nbins = 90; float x1 = 10.; float x2 = 100.;

   TString directory[4] = {"EG#minusL1", "EG#minusL2","EG#minusL3","EG#minusL4"};

   ofstream fit_result;
   char fit_parameter[50];
   sprintf(fit_parameter, "./dphi.txt");

   fit_result.open(fit_parameter);
   
   //for( int nth_sw = 0; nth_sw < 4; nth_sw++){
      //if(nth_sw != 5) continue;
   
      float low_et = 0.;
      float high_et = 0.;
   
      for( int nth = 0; nth < 90; nth++){
         low_et = 10. + nth;
         high_et = low_et + 1.;
      
         vector<float> pmDPhi1_;
      
         Long64_t nbytes = 0, nb = 0;
         for (Long64_t jentry=0; jentry<nentries;jentry++) {
            Long64_t ientry = LoadTree(jentry);
            if (ientry < 0) break;
            nb = fChain->GetEntry(jentry);   nbytes += nb;
            // if (Cut(ientry) < 0) continue;
            
            // test for bpix1
            if(closest_egEt > low_et && closest_egEt < high_et && fabs(closest_egEta) < 1.4 && closest_eg_dr < 0.1){
            int size_em_phi_hit = 0;
            size_em_phi_hit = EM_pixel_dPhi->size();

              for(int i=0; i < size_em_phi_hit; i++){
                    if(fabs(EM_pixel_dPhi->at(i)) < 0.15) pmDPhi1_.push_back(EM_pixel_dPhi->at(i));
                    pmDPhi1_dist->Fill(closest_egEt, EM_pixel_dPhi->at(i));
              }
            }
      
         }// end of event loop
   
        std::sort (pmDPhi1_.begin(), pmDPhi1_.end()); 
   
        int size = pmDPhi1_.size();
        if(size == 0) continue;
   
        median[nth] = getMedian(pmDPhi1_);; 
        x[nth] = low_et; 
        x_err[nth] = 0.;
        median_err[nth] = getMedianErr(pmDPhi1_);
      
      }// end of et loop

      TCanvas *c1 = new TCanvas("c1","c1",800,600);
      gStyle->SetOptStat(0);
      gStyle->SetPalette(1);
      c1->SetRightMargin(0.1);
      c1->SetLeftMargin(0.15);
      c1->SetBottomMargin(0.15);
      c1->SetGridy();
      c1->SetGridx();
      c1->SetLogz();

      pmDPhi1_dist->Draw("COLZ");
      c1->Update();
      TPaletteAxis *palette = (TPaletteAxis*)pmDPhi1_dist->GetListOfFunctions()->FindObject("palette");
      float x2_palette = palette->GetX2NDC();
      float x1_palette = palette->GetX1NDC();
      palette->SetX2NDC( x2_palette - (x2_palette-x1_palette) * 0.15 );
      palette->SetX1NDC( x1_palette + (x2_palette-x1_palette) * 0.15 );

      pmDPhi1_dist->GetYaxis()->SetDecimals(3);
      pmDPhi1_dist->GetXaxis()->SetTitle("L1 EG E_{T} (crystal) [GeV]");
      pmDPhi1_dist->GetXaxis()->SetTitleSize(0.05);
      pmDPhi1_dist->GetYaxis()->SetTitle("#Delta#phi [rad.]");
      pmDPhi1_dist->GetYaxis()->SetTitleSize(0.05);
      pmDPhi1_dist->GetYaxis()->SetTitleOffset(1.35);

      TGraphErrors* median_points = new TGraphErrors(binSize,x,median,x_err,median_err);
      median_points->SetMarkerStyle(24);
      median_points->SetMarkerSize(0.3);
      median_points->Draw("same peZ");

      TGraph* median_gr = new TGraph(binSize,x,median);

      TLatex eta_range;
      eta_range.SetTextSize(0.04);
      eta_range.DrawLatexNDC(0.18,0.85, "|#eta| < 1.4"); 

      TF1 *fit_median_func = new TF1("func_median","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 100.);

      median_gr->Fit(fit_median_func,"0W");
      fit_median_func->Draw("lsame");
      fit_median_func->SetLineColor(kBlack);
      fit_median_func->SetLineStyle(1);

      TLegend* leg = new TLegend(0.45, 0.1, 0.85, 0.3,"","brNDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.04);
      leg->SetBorderSize(0);
      //leg->Draw();

      // just to comapre to old ROI windows
      TF1 *fit_func_old_up = new TF1("func_old_up","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 100.);
      TF1 *fit_func_old_down = new TF1("func_old_down","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10, 100);

      fit_func_old_up->SetParameter(0,-0.00616578);
      fit_func_old_up->SetParameter(1,0.0130247);
      fit_func_old_up->SetParameter(2,1.06586);
      fit_func_old_up->SetParameter(3,0.291827);
      fit_func_old_up->SetParameter(4,0.418001);

      fit_func_old_down->SetParameter(0,-0.0476494);
      fit_func_old_down->SetParameter(1, -0.51066);
      fit_func_old_down->SetParameter(2,-0.175453);
      fit_func_old_down->SetParameter(3, 0.305778);
      fit_func_old_down->SetParameter(4, 0.798656);

      fit_func_old_up->SetLineColor(kRed);
      fit_func_old_down->SetLineColor(kRed);

      //fit_func_old_up->Draw("lsame");
      //fit_func_old_down->Draw("lsame");

      fit_result << endl;
      for( int i=0; i < 5; i++){
          fit_result << "p[" << i << "] = " << fit_median_func->GetParameter(i) << ";" << endl;
      }
      fit_result << endl;

      c1->SaveAs("./EG_Pixel.pdf");
      pmDPhi1_dist->Reset();
      delete c1;
   //}// loop for each signal windows

   fit_result.close();

}
