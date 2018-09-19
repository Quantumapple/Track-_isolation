#define gaussFit_cxx
#include "gaussFit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void gaussFit::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   TH1F *hist_a = new TH1F("hist_a",";P_{T} (GeV); #Delta#phi", 10, 2., 3.);     // 0.1 GeV
   TH1F *hist_b = new TH1F("hist_b",";P_{T} (GeV); #Delta#phi", 15, 3., 6.);     // 0.2 GeV
   TH1F *hist_c = new TH1F("hist_c",";P_{T} (GeV); #Delta#phi", 4, 6., 8.);      // 0.5 GeV
   TH1F *hist_d = new TH1F("hist_d",";P_{T} (GeV); #Delta#phi", 12, 8., 20.);    // 1 GeV
   TH1F *hist_e = new TH1F("hist_e",";P_{T} (GeV); #Delta#phi", 15, 20., 50.);   // 2 GeV
   TH1F *hist_f = new TH1F("hist_f",";P_{T} (GeV); #Delta#phi", 10, 50., 100.);  // 5 GeV
   
   ofstream points;
   TString file_name = "fitting_points.h";
   Int_t count = 0;

   points.open(file_name.Data());

   points << "#ifndef fitting_points_h" << endl;
   points << "#define fitting_points_h" << endl;
   points << endl;

   points << "double gen_Pt[66] ={}, point[66] = {}, error[66] = {};" << endl;
   points << "int points = 66;" << endl;

   points << "void set_arrary(){" << endl;
   points << endl;
   
   Float_t low_pt = 2.;
   Float_t high_pt = 3.;
   Float_t min = -1.;
   Float_t max = 1.;
   
   TH1F *H[70];
   Char_t histname[15];
   TF1 *func = new TF1("func","[0]*TMath::Gaus(x,[1],[2],1)",min,max);
   
   for(Int_t nth = 0; nth < 10; nth++)
   {
       cout << "Process: (" << count+1 << "/66)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 2.+ nth*0.1;
       high_pt = low_pt + 0.1;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",25,-0.1,0.1);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   Int_t size = genpT_bsd2_d2d3->size();
	   if( size == 0 ) continue;

       Float_t gen_pT = genpT_bsd2_d2d3->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = bsd2_d2d3_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.06, H[count]->GetMean()+0.06);
       func->SetParameters(70, H[count]->GetMean(), 0.005);
       
       func->SetParLimits(0.,1e-4,10);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << ";    point[" << count << "] = " << fabs(func->GetParameter(1)) << ";     error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_a->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       //c1->Write();

       c1->Clear();
       delete c1;

       count++;
       cout << endl; 

   } // end of nth loop ( 2 ~ 3 GeV )
    
    
   for(Int_t nth = 0; nth < 15; nth++)
   {
       cout << "Process: (" << count+1 << "/66)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 3.+ nth*0.2;
       high_pt = low_pt + 0.2;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",50,-0.1,0.1);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   Int_t size = genpT_bsd2_d2d3->size();
	   if( size == 0 ) continue;

       Float_t gen_pT = genpT_bsd2_d2d3->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = bsd2_d2d3_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.04, H[count]->GetMean()+0.04);
       func->SetParameters(70, H[count]->GetMean(), 0.005);
       
       func->SetParLimits(0.,1e-4,10);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << ";    point[" << count << "] = " << fabs(func->GetParameter(1)) << ";     error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_b->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       //c1->Write();

       c1->Clear();
       delete c1;

       count++;
       cout << endl; 

   } // end of nth loop ( 3 ~ 6 GeV )
   
   for(Int_t nth = 0; nth < 4; nth++)
   {
       cout << "Process: (" << count+1 << "/66)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 6.+ nth*0.5;
       high_pt = low_pt + 0.5;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",50,-0.05,0.05);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;
	   
       Int_t size = genpT_bsd2_d2d3->size();
	   if( size == 0 ) continue;

       Float_t gen_pT = genpT_bsd2_d2d3->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = bsd2_d2d3_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.02, H[count]->GetMean()+0.02);
       func->SetParameters(70, H[count]->GetMean(), 0.01);
       
       func->SetParLimits(0.,1e-4,50);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << ";    point[" << count << "] = " << fabs(func->GetParameter(1)) << ";     error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_c->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       //c1->Write();

       c1->Clear();
       delete c1;

       count++;
       cout << endl; 

   } // end of nth loop ( 6 ~ 8 GeV )
   
   for(Int_t nth = 0; nth < 12; nth++)
   {
       cout << "Process: (" << count+1 << "/66)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 8.+ nth;
       high_pt = low_pt + 1.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",50,-0.025,0.025);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   Int_t size = genpT_bsd2_d2d3->size();
	   if( size == 0 ) continue;

       Float_t gen_pT = genpT_bsd2_d2d3->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = bsd2_d2d3_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParameters(70, H[count]->GetMean(), 0.001);
       
       func->SetParLimits(0.,1e-4,50);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << ";    point[" << count << "] = " << fabs(func->GetParameter(1)) << ";     error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_d->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       //c1->Write();

       c1->Clear();
       delete c1;

       count++;
       cout << endl; 

   } // end of nth loop ( 8 ~ 20 GeV )
   
   for(Int_t nth = 0; nth < 15; nth++)
   {
       cout << "Process: (" << count+1 << "/66)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 20.+ nth*2;
       high_pt = low_pt + 2.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",100,-0.025,0.025);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   Int_t size = genpT_bsd2_d2d3->size();
	   if( size == 0 ) continue;

       Float_t gen_pT = genpT_bsd2_d2d3->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = bsd2_d2d3_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParameters(70, H[count]->GetMean(), 0.001);
       
       func->SetParLimits(0.,1e-4,50);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << ";    point[" << count << "] = " << fabs(func->GetParameter(1)) << ";     error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_e->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       //c1->Write();

       c1->Clear();
       delete c1;

       count++;
       cout << endl; 

   } // end of nth loop ( 20 ~ 50 GeV )
   
   for(Int_t nth = 0; nth < 10; nth++)
   {
       cout << "Process: (" << count+1 << "/66)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 50.+ nth*5;
       high_pt = low_pt + 5.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",100,-0.01,0.01);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   Int_t size = genpT_bsd2_d2d3->size();
	   if( size == 0 ) continue;

       Float_t gen_pT = genpT_bsd2_d2d3->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = bsd2_d2d3_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.005, H[count]->GetMean()+0.005);
       func->SetParameters(70, H[count]->GetMean(), 0.001);
       
       func->SetParLimits(0.,1e-4,50);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << ";    point[" << count << "] = " << fabs(func->GetParameter(1)) << ";     error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_f->SetBinContent(nth+1, fabs(func->GetParameter(1)));
       gStyle->SetOptFit(1111);

       func->SetLineColor(kRed);
       func->SetLineWidth(2);

       TCanvas* c1 = new TCanvas("c1","",800,800) ;
       c1->cd(1); 
       gPad->SetLeftMargin(0.15);  
       H[count]->Draw();
       func->Draw("same");

       c1->SaveAs(file_);
       //c1->Write();

       c1->Clear();
       delete c1;

       count++;
       cout << endl; 

   } // end of nth loop ( 50 ~ 100 GeV )

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

   Int_t nbins = 105; Float_t x1 = 0.; Float_t x2 = 105.;
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

   hist_a->SetMarkerColor(1);
   hist_a->SetMarkerStyle(20);
   hist_a->SetMarkerSize(1.);
   hist_a->Draw("hist p same");

   hist_b->SetMarkerColor(1);
   hist_b->SetMarkerStyle(20);
   hist_b->SetMarkerSize(1.);
   hist_b->Draw("hist p same");

   hist_c->SetMarkerColor(1);
   hist_c->SetMarkerStyle(20);
   hist_c->SetMarkerSize(1.);
   hist_c->Draw("hist p same");

   hist_d->SetMarkerColor(1);
   hist_d->SetMarkerStyle(20);
   hist_d->SetMarkerSize(1.);
   hist_d->Draw("hist p same");
   
   hist_e->SetMarkerColor(1);
   hist_e->SetMarkerStyle(20);
   hist_e->SetMarkerSize(1.);
   hist_e->Draw("hist p same");
   
   hist_f->SetMarkerColor(1);
   hist_f->SetMarkerStyle(20);
   hist_f->SetMarkerSize(1.);
   hist_f->Draw("hist p same");

   points << "}" << endl;
   points << "#endif" << endl;
   points.close();

   c2->Print("Fit_points.pdf");

   TFile *output = new TFile("fit_test.root", "RECREATE");
   h1->Write();
   hist_a->Write();
   hist_b->Write();
   hist_c->Write();
   hist_d->Write();
   hist_e->Write();
   hist_f->Write();
}
