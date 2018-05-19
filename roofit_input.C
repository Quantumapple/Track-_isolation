#define roofit_input_cxx
#include "roofit_input.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooChebychev.h"
#include "RooLandau.h"

#include <iostream>
#include <fstream>

using namespace RooFit ;

void roofit_input::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
  
   Int_t count = 0;
   
   Int_t nbins = 150; Float_t x1 = 0.; Float_t x2 = 150.;
   TH1F* hist_bsl1_l1l2_a = new TH1F("hist_bsl1_l1l2_a",";P_{T} (GeV); #Delta#phi", 18, 0.2, 2.0);
   TH1F* hist_bsl1_l1l2_b = new TH1F("hist_bsl1_l1l2_b",";P_{T} (GeV); #Delta#phi", 23, 2.0, 25.);
   TH1F* hist_bsl1_l1l2_c = new TH1F("hist_bsl1_l1l2_c",";P_{T} (GeV); #Delta#phi", 10, 25., 50.);
   TH1F* hist_bsl1_l1l2_d = new TH1F("hist_bsl1_l1l2_d",";P_{T} (GeV); #Delta#phi", 20, 50., 150.);

   ofstream Bsl1_l1l2_points;
   Bsl1_l1l2_points.open("./Bsl1_l1l2_points.h");

   Bsl1_l1l2_points << "#ifndef Bsl1_l1l2_points_h" << endl;
   Bsl1_l1l2_points << "#define Bsl1_l1l2_points_h" << endl;
   Bsl1_l1l2_points << endl;

   Bsl1_l1l2_points << "double gen_Pt[70] ={}, point[70] = {}; " << endl;
   Bsl1_l1l2_points << "int points = 70;" << endl;

   Bsl1_l1l2_points << "void set_arrary(){" << endl;
   Bsl1_l1l2_points << endl;

   Float_t low_pt = 0.;
   Float_t high_pt = 1.;
   
   TFile *output = new TFile("fit.root","RECREATE");
   TH1F *H[75];
   Char_t histname[15];

   Float_t min = -1.;
   Float_t max = 1.;

   TF1 *func = new TF1("func","[0]*TMath::Gaus(x,[1],[2],1)",min,max);

   for(Int_t nth = 2; nth < 4; nth++)
   {
       cout << "Process: (" << count+1 << "/70)" << endl;
       TString file_ = "./Fix_results/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 0.+ nth*0.1;
       high_pt = low_pt + 0.1;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",250,-0.3,0.05);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   int size = bsl1_l1l2_dphi->size();
	   int pT_size = genpT->size(); 

	   if( size == 0 ) continue;

	   //cout << "Size of dphi vector: " << size << ", size of pT vector: " << pT_size << endl;
	   //cout << "dPhi: " << bsl1_l1l2_dphi->at(0) << ", pT: " << genpT->at(0) << endl;

           Float_t gen_pT = genpT->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;
	   //cout << "generator level pT: " << gen_pT << endl;
           
	   H[count]->Fill(bsl1_l1l2_dphi->at(0));

       }// end of event loop

       func->SetParameters(70, H[count]->GetMean(), 0.01);
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       Bsl1_l1l2_points << "gen_Pt[" << count << "] = " << count << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << ";" << endl;
       hist_bsl1_l1l2_a->SetBinContent(count+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       c1->Write();

       c1->Clear();
       delete c1;
  
       count++;
       cout << endl; 

   } // end of nth loop ( 0.2 ~ 0.3 GeV )
   
   for(Int_t nth = 4; nth < 10; nth++)
   {
       cout << "Process: (" << count+1 << "/70)" << endl;
       TString file_ = "./Fix_results/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 0.+ nth*0.1;
       high_pt = low_pt + 0.1;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",250,-0.11,0.03);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   int size = bsl1_l1l2_dphi->size();
	   int pT_size = genpT->size(); 

	   if( size == 0 ) continue;

	   //cout << "Size of dphi vector: " << size << ", size of pT vector: " << pT_size << endl;
	   //cout << "dPhi: " << bsl1_l1l2_dphi->at(0) << ", pT: " << genpT->at(0) << endl;

           Float_t gen_pT = genpT->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;
	   //cout << "generator level pT: " << gen_pT << endl;
           
	   H[count]->Fill(bsl1_l1l2_dphi->at(0));

       }// end of event loop

       func->SetParameters(70, H[count]->GetMean(), 0.001);
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       Bsl1_l1l2_points << "gen_Pt[" << count << "] = " << count << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << ";" << endl;
       hist_bsl1_l1l2_a->SetBinContent(count+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       c1->Write();

       c1->Clear();
       delete c1;
  
       count++;
       cout << endl; 

   } // end of nth loop ( 0.4 ~ 1 GeV )
  
   for(Int_t nth = 10; nth < 20; nth++)
   {
       cout << "Process: (" << count+1 << "/70)" << endl;
       TString file_ = "./Fix_results/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 0.+ nth*0.1;
       high_pt = low_pt + 0.1;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",250,-0.04,0.);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   int size = bsl1_l1l2_dphi->size();
	   int pT_size = genpT->size(); 

	   if( size == 0 ) continue;

	   //cout << "Size of dphi vector: " << size << ", size of pT vector: " << pT_size << endl;
	   //cout << "dPhi: " << bsl1_l1l2_dphi->at(0) << ", pT: " << genpT->at(0) << endl;

           Float_t gen_pT = genpT->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;
	   //cout << "generator level pT: " << gen_pT << endl;
           
	   H[count]->Fill(bsl1_l1l2_dphi->at(0));

       }// end of event loop

       func->SetParameters(70, H[count]->GetMean(), 0.001);
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       Bsl1_l1l2_points << "gen_Pt[" << count << "] = " << count << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << ";" << endl;
       hist_bsl1_l1l2_a->SetBinContent(count+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       c1->Write();

       c1->Clear();
       delete c1;
  
       count++;
       cout << endl; 

   } // end of nth loop ( 1 ~ 2 GeV )

   for(Int_t nth = 2; nth < 25; nth++)
   {
       cout << "Process: (" << count+1 << "/70)" << endl;
       TString file_ = "./Fix_results/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = nth;
       high_pt = low_pt + 1.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",250,-0.04,0.04);
       
       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   int size = bsl1_l1l2_dphi->size();
	   int pT_size = genpT->size(); 

	   if( size == 0 ) continue;

	   //cout << "Size of dphi vector: " << size << ", size of pT vector: " << pT_size << endl;
	   //cout << "dPhi: " << bsl1_l1l2_dphi->at(0) << ", pT: " << genpT->at(0) << endl;

           Float_t gen_pT = genpT->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;
	   //cout << "generator level pT: " << gen_pT << endl;

	   H[count]->Fill(bsl1_l1l2_dphi->at(0));

       }// end of event loop
       cout << endl;

       func->SetParameters(70, H[count]->GetMean(), 0.001);
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       Bsl1_l1l2_points << "gen_Pt[" << count << "] = " << count << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << ";" << endl;
       hist_bsl1_l1l2_b->SetBinContent(nth-1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       c1->Write();

       c1->Clear();
       delete c1;
  
       count++;
       cout << endl; 

   } // end of nth loop ( 2 ~ 25 GeV )
  
   
   for(Int_t nth = 0; nth < 10; nth++)
   {
       cout << "Process: (" << count+1 << "/70)" << endl;
       TString file_ = "./Fix_results/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";
       
       low_pt = 25. + nth*2.5;
       high_pt = low_pt + 2.5;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",250,-0.012,0.012);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
       //for (Long64_t jentry=0; jentry<10;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   int size = bsl1_l1l2_dphi->size();
	   int pT_size = genpT->size(); 

	   if( size == 0 ) continue;

	   //cout << "Size of dphi vector: " << size << ", size of pT vector: " << pT_size << endl;
	   //cout << "dPhi: " << bsl1_l1l2_dphi->at(0) << ", pT: " << genpT->at(0) << endl;

           Float_t gen_pT = genpT->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;
	   //cout << "generator level pT: " << gen_pT << endl;

	   H[count]->Fill(bsl1_l1l2_dphi->at(0));

       }// end of event loop
       cout << endl;
       
       func->SetParameters(70, H[count]->GetMean(), 0.001);
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       Bsl1_l1l2_points << "gen_Pt[" << count << "] = " << count << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << ";" << endl;
       hist_bsl1_l1l2_c->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       c1->Write();

       c1->Clear();
       delete c1;
  
       count++;
       cout << endl; 

   } // end of nth loop ( 25 ~ 50 GeV )

   for(Int_t nth = 0; nth < 20; nth++)
   {
       cout << "Process: (" << count+1 << "/70)" << endl;
       TString file_ = "./Fix_results/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";
       
       low_pt = 50. + nth*5.;
       high_pt = low_pt + 5.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",250,-0.005,0.005);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
       //for (Long64_t jentry=0; jentry<10;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   int size = bsl1_l1l2_dphi->size();
	   int pT_size = genpT->size(); 

	   if( size == 0 ) continue;

	   //cout << "Size of dphi vector: " << size << ", size of pT vector: " << pT_size << endl;
	   //cout << "dPhi: " << bsl1_l1l2_dphi->at(0) << ", pT: " << genpT->at(0) << endl;

           Float_t gen_pT = genpT->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;
	   //cout << "generator level pT: " << gen_pT << endl;

	   H[count]->Fill(bsl1_l1l2_dphi->at(0));

       }// end of event loop
       cout << endl;
       
       func->SetParameters(70, H[count]->GetMean(), 0.001);
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       Bsl1_l1l2_points << "gen_Pt[" << count << "] = " << count << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << ";" << endl;
       hist_bsl1_l1l2_d->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       c1->Write();

       c1->Clear();
       delete c1;
  
       count++;
       cout << endl; 

   } // end of nth loop ( 50 ~ 150 GeV )
   
   TCanvas *c2 = new TCanvas("c2","c2",700,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(2); // axis width, default is 1
   c2->SetTopMargin(0.05);
   c2->SetBottomMargin(0.12);
   c2->SetRightMargin(0.03);
   c2->SetLeftMargin(0.2);
   c2->SetGrid();
   c2->SetTicky(1);
   c2->SetTickx(1);
   c2->SetLogy();

   TH1F* h1= new TH1F("h1",";P_{T} (GeV); #Delta#phi", nbins, x1, x2);
   h1->SetTitle("");
   h1->GetXaxis()->SetTitleOffset(1.3);
   h1->GetXaxis()->SetTitleSize(0.045);
   h1->GetXaxis()->SetNdivisions(505);
   h1->GetYaxis()->SetNdivisions(506);
   h1->GetXaxis()->SetLabelSize(0.05);
   h1->GetYaxis()->SetLabelSize(0.05);
   h1->GetXaxis()->SetRangeUser(x1,x2);
   h1->GetYaxis()->SetRangeUser(0.0001, 1.);
   h1->GetYaxis()->SetTitleOffset(1.5);
   h1->GetYaxis()->SetTitleSize(0.050);
   h1->Draw();

   hist_bsl1_l1l2_a->SetMarkerColor(1);
   hist_bsl1_l1l2_a->SetMarkerStyle(20);
   hist_bsl1_l1l2_a->SetMarkerSize(1.);
   hist_bsl1_l1l2_a->Draw("hist p same");
  
   hist_bsl1_l1l2_b->SetMarkerColor(1);
   hist_bsl1_l1l2_b->SetMarkerStyle(20);
   hist_bsl1_l1l2_b->SetMarkerSize(1.);
   hist_bsl1_l1l2_b->Draw("hist p same");
   
   hist_bsl1_l1l2_c->SetMarkerColor(1);
   hist_bsl1_l1l2_c->SetMarkerStyle(20);
   hist_bsl1_l1l2_c->SetMarkerSize(1.);
   hist_bsl1_l1l2_c->Draw("hist p same");
  
   hist_bsl1_l1l2_d->SetMarkerColor(1);
   hist_bsl1_l1l2_d->SetMarkerStyle(20);
   hist_bsl1_l1l2_d->SetMarkerSize(1.);
   hist_bsl1_l1l2_d->Draw("hist p same");

   Bsl1_l1l2_points << "}" << endl;
   Bsl1_l1l2_points << "#endif" << endl;
   Bsl1_l1l2_points.close();

   h1->Write();
   hist_bsl1_l1l2_a->Write();
   hist_bsl1_l1l2_b->Write();
   hist_bsl1_l1l2_c->Write();
   hist_bsl1_l1l2_d->Write();
   c2->Print("Fix_results/Bsl1_l1l2_dphi.pdf");
   

}
