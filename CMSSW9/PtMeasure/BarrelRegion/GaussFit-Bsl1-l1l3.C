#define GaussFit_cxx
#include "GaussFit.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void GaussFit::Loop(int nth_case)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   TH1F *hist_bsl1_l1l2_a = new TH1F("hist_bsl1_l1l2_a",";P_{T} (GeV); #Delta#phi", 2, 2, 4);       // 1 GeV
   TH1F *hist_bsl1_l1l2_b = new TH1F("hist_bsl1_l1l2_b",";P_{T} (GeV); #Delta#phi", 6, 4., 10.);    // 1 GeV
   TH1F *hist_bsl1_l1l2_c = new TH1F("hist_bsl1_l1l2_c",";P_{T} (GeV); #Delta#phi", 10, 10., 30.);   // 2 GeV
   TH1F *hist_bsl1_l1l2_d = new TH1F("hist_bsl1_l1l2_d",";P_{T} (GeV); #Delta#phi", 10, 30., 60.);   // 3 GeV
   TH1F *hist_bsl1_l1l2_e = new TH1F("hist_bsl1_l1l2_e",";P_{T} (GeV); #Delta#phi", 28, 60., 200.);   // 5 GeV

   ofstream points;
   TString file_name = "fitting_points.h";
   Int_t count = 0;

   points.open(file_name.Data());

   points << "#ifndef fitting_points_h" << endl;
   points << "#define fitting_points_h" << endl;
   points << endl;

   points << "double gen_Pt[56] ={}, point[56] = {}, error[56] = {};" << endl;
   points << "int points = 56;" << endl;

   points << "void set_arrary(){" << endl;
   points << endl;

   Float_t low_pt = 0.;
   Float_t high_pt = 1.;
   Float_t min = -1.;
   Float_t max = 1.;

   TH1F *H[70];
   Char_t histname[15];
   TF1 *func = new TF1("func","[0]*TMath::Gaus(x,[1],[2],1)",min,max);
   
   for(Int_t nth = 0; nth < 2; nth++)
   {
       cout << "Process: (" << count+1 << "/56)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 2.+ nth;
       high_pt = low_pt + 1.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",500,-0.1,0.1);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   Int_t size = 0;
	   if( nth_case == 1 ) size = genpT_case1->size();
	   if( nth_case == 2 ) size = genpT_case2->size();
	   if( nth_case == 3 ) size = genpT_case3->size();
	   if( nth_case == 4 ) size = genpT_case4->size();
	   if( nth_case == 5 ) size = genpT_case5->size();

	   if( size == 0 ) continue;

       Float_t gen_pT = -1.;
	   if( nth_case == 1 ) gen_pT = genpT_case1->at(0);
	   if( nth_case == 2 ) gen_pT = genpT_case2->at(0);
	   if( nth_case == 3 ) gen_pT = genpT_case3->at(0);
	   if( nth_case == 4 ) gen_pT = genpT_case4->at(0);
	   if( nth_case == 5 ) gen_pT = genpT_case5->at(0);

	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = -999.;
	   if( nth_case == 1 ) dphi = bsl1_l1l3_dphi->at(0);
	   if( nth_case == 2 ) dphi = bsl1_l1l4_dphi->at(0);
	   if( nth_case == 3 ) dphi = bsl2_l2l3_dphi->at(0);
	   if( nth_case == 4 ) dphi = bsl2_l2l4_dphi->at(0);
	   if( nth_case == 5 ) dphi = bsl3_l3l4_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.02, H[count]->GetMean()+0.02);
       func->SetParameters(70, H[count]->GetMean(), 0.001);
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << "; error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_bsl1_l1l2_a->SetBinContent(nth+1, fabs(func->GetParameter(1)));
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

   } // end of nth loop ( 2 ~ 4 GeV )

   for(Int_t nth = 0; nth < 6; nth++)
   {
       cout << "Process: (" << count+1 << "/56)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 4.+ nth;
       high_pt = low_pt + 1.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",500,-0.08,0.08);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   Int_t size = 0;
	   if( nth_case == 1 )size = genpT_case1->size();
	   if( nth_case == 2 )size = genpT_case2->size();
	   if( nth_case == 3 )size = genpT_case3->size();
	   if( nth_case == 4 )size = genpT_case4->size();
	   if( nth_case == 5 )size = genpT_case5->size();

	   if( size == 0 ) continue;

       Float_t gen_pT = -1.;
	   if( nth_case == 1 ) gen_pT = genpT_case1->at(0);
	   if( nth_case == 2 ) gen_pT = genpT_case2->at(0);
	   if( nth_case == 3 ) gen_pT = genpT_case3->at(0);
	   if( nth_case == 4 ) gen_pT = genpT_case4->at(0);
	   if( nth_case == 5 ) gen_pT = genpT_case5->at(0);

	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = -999.;
	   if( nth_case == 1 ) dphi = bsl1_l1l3_dphi->at(0);
	   if( nth_case == 2 ) dphi = bsl1_l1l4_dphi->at(0);
	   if( nth_case == 3 ) dphi = bsl2_l2l3_dphi->at(0);
	   if( nth_case == 4 ) dphi = bsl2_l2l4_dphi->at(0);
	   if( nth_case == 5 ) dphi = bsl3_l3l4_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.008, H[count]->GetMean()+0.008);
       func->SetParameters(70, H[count]->GetMean(), 0.001);

       func->SetParLimits(0.,1e-4,10);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);

       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << "; error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_bsl1_l1l2_b->SetBinContent(nth+1, fabs(func->GetParameter(1)));
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

   } // end of nth loop ( 4 ~ 10 GeV )

   for(Int_t nth = 0; nth < 4; nth++)
   {
       cout << "Process: (" << count+1 << "/56)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 10.+ nth*2;
       high_pt = low_pt + 2.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",500,-0.06,0.06);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;
	   
	   Int_t size = 0;
	   if( nth_case == 1 )size = genpT_case1->size();
	   if( nth_case == 2 )size = genpT_case2->size();
	   if( nth_case == 3 )size = genpT_case3->size();
	   if( nth_case == 4 )size = genpT_case4->size();
	   if( nth_case == 5 )size = genpT_case5->size();

	   if( size == 0 ) continue;

       Float_t gen_pT = -1.;
	   if( nth_case == 1 ) gen_pT = genpT_case1->at(0);
	   if( nth_case == 2 ) gen_pT = genpT_case2->at(0);
	   if( nth_case == 3 ) gen_pT = genpT_case3->at(0);
	   if( nth_case == 4 ) gen_pT = genpT_case4->at(0);
	   if( nth_case == 5 ) gen_pT = genpT_case5->at(0);

	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = -999.;
	   if( nth_case == 1 ) dphi = bsl1_l1l3_dphi->at(0);
	   if( nth_case == 2 ) dphi = bsl1_l1l4_dphi->at(0);
	   if( nth_case == 3 ) dphi = bsl2_l2l3_dphi->at(0);
	   if( nth_case == 4 ) dphi = bsl2_l2l4_dphi->at(0);
	   if( nth_case == 5 ) dphi = bsl3_l3l4_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.005, H[count]->GetMean()+0.005);
       func->SetParameters(70, H[count]->GetMean(), 0.0001);
       
       func->SetParLimits(0.,1e-4,10);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << "; error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
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

   } // end of nth loop ( 10 ~ 18 GeV )

   for(Int_t nth = 0; nth < 6; nth++)
   {
       cout << "Process: (" << count+1 << "/56)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 18.+ nth*2;
       high_pt = low_pt + 2.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",1000,-0.06,0.06);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;
	   
	   Int_t size = 0;
	   if( nth_case == 1 )size = genpT_case1->size();
	   if( nth_case == 2 )size = genpT_case2->size();
	   if( nth_case == 3 )size = genpT_case3->size();
	   if( nth_case == 4 )size = genpT_case4->size();
	   if( nth_case == 5 )size = genpT_case5->size();

	   if( size == 0 ) continue;

       Float_t gen_pT = -1.;
	   if( nth_case == 1 ) gen_pT = genpT_case1->at(0);
	   if( nth_case == 2 ) gen_pT = genpT_case2->at(0);
	   if( nth_case == 3 ) gen_pT = genpT_case3->at(0);
	   if( nth_case == 4 ) gen_pT = genpT_case4->at(0);
	   if( nth_case == 5 ) gen_pT = genpT_case5->at(0);

	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = -999.;
	   if( nth_case == 1 ) dphi = bsl1_l1l3_dphi->at(0);
	   if( nth_case == 2 ) dphi = bsl1_l1l4_dphi->at(0);
	   if( nth_case == 3 ) dphi = bsl2_l2l3_dphi->at(0);
	   if( nth_case == 4 ) dphi = bsl2_l2l4_dphi->at(0);
	   if( nth_case == 5 ) dphi = bsl3_l3l4_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.003, H[count]->GetMean()+0.003);
       func->SetParameters(70, H[count]->GetMean(), 0.0001);
       
       func->SetParLimits(0.,1e-4,10);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << "; error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_bsl1_l1l2_c->SetBinContent(nth+4, fabs(func->GetParameter(1)));
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

   } // end of nth loop ( 18 ~ 30 GeV )
  
      
   for(Int_t nth = 0; nth < 10; nth++)
   {
       cout << "Process: (" << count+1 << "/56)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 30. + nth*3;
       high_pt = low_pt + 3.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",1000,-0.05,0.05);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;
	   
	   Int_t size = 0;
	   if( nth_case == 1 )size = genpT_case1->size();
	   if( nth_case == 2 )size = genpT_case2->size();
	   if( nth_case == 3 )size = genpT_case3->size();
	   if( nth_case == 4 )size = genpT_case4->size();
	   if( nth_case == 5 )size = genpT_case5->size();

	   if( size == 0 ) continue;

       Float_t gen_pT = -1.;
	   if( nth_case == 1 ) gen_pT = genpT_case1->at(0);
	   if( nth_case == 2 ) gen_pT = genpT_case2->at(0);
	   if( nth_case == 3 ) gen_pT = genpT_case3->at(0);
	   if( nth_case == 4 ) gen_pT = genpT_case4->at(0);
	   if( nth_case == 5 ) gen_pT = genpT_case5->at(0);

	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = -999.;
	   if( nth_case == 1 ) dphi = bsl1_l1l3_dphi->at(0);
	   if( nth_case == 2 ) dphi = bsl1_l1l4_dphi->at(0);
	   if( nth_case == 3 ) dphi = bsl2_l2l3_dphi->at(0);
	   if( nth_case == 4 ) dphi = bsl2_l2l4_dphi->at(0);
	   if( nth_case == 5 ) dphi = bsl3_l3l4_dphi->at(0);

       H[count]->Fill(dphi);

       }// end of event loop
       cout << endl;

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.002, H[count]->GetMean()+0.002);
       func->SetParameters(70, H[count]->GetMean(), 0.0001);
       
       func->SetParLimits(0.,1e-4,10);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << "; error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
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

   } // end of nth loop ( 30 ~ 60 GeV )

   for(Int_t nth = 0; nth < 28; nth++)
   {
       cout << "Process: (" << count+1 << "/56)" << endl;
       TString file_ = "./Plots/phi_";
       TString ith_;

       ith_.Form("%d", count+1);
       file_ = file_ + ith_ + "th.pdf";

       low_pt = 60. + nth*5.;
       high_pt = low_pt + 5.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       sprintf(histname, "hist_%dth", count+1);
       H[count] = new TH1F(histname,"",1000,-0.05,0.05);

       for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;
	   
	   Int_t size = 0;
	   if( nth_case == 1 )size = genpT_case1->size();
	   if( nth_case == 2 )size = genpT_case2->size();
	   if( nth_case == 3 )size = genpT_case3->size();
	   if( nth_case == 4 )size = genpT_case4->size();
	   if( nth_case == 5 )size = genpT_case5->size();

	   if( size == 0 ) continue;

       Float_t gen_pT = -1.;
	   if( nth_case == 1 ) gen_pT = genpT_case1->at(0);
	   if( nth_case == 2 ) gen_pT = genpT_case2->at(0);
	   if( nth_case == 3 ) gen_pT = genpT_case3->at(0);
	   if( nth_case == 4 ) gen_pT = genpT_case4->at(0);
	   if( nth_case == 5 ) gen_pT = genpT_case5->at(0);

	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;

       Float_t dphi = -999.;
	   if( nth_case == 1 ) dphi = bsl1_l1l3_dphi->at(0);
	   if( nth_case == 2 ) dphi = bsl1_l1l4_dphi->at(0);
	   if( nth_case == 3 ) dphi = bsl2_l2l3_dphi->at(0);
	   if( nth_case == 4 ) dphi = bsl2_l2l4_dphi->at(0);
	   if( nth_case == 5 ) dphi = bsl3_l3l4_dphi->at(0);

       H[count]->Fill(dphi);
       
       }// end of event loop
       cout << endl;

       H[count]->GetXaxis()->SetRangeUser(H[count]->GetMean()-0.002, H[count]->GetMean()+0.002);
       func->SetParameters(70, H[count]->GetMean(), 0.0001);
       
       func->SetParLimits(0.,1e-4,10);
       func->SetParLimits(1,H[count]->GetMean()-0.01, H[count]->GetMean()+0.01);
       func->SetParLimits(2,1e-5,0.1);
       
       H[count]->Fit(func,"0R");
       cout << "Fit parameters: " << "Const: " << func->GetParameter(0) << ", Mean: " << func->GetParameter(1) << ", sigma: " << func->GetParameter(2) << endl;
       points << "gen_Pt[" << count << "] = " << low_pt << "; point[" << count << "] = " << fabs(func->GetParameter(1)) << "; error[" << count << "] = " << fabs(func->GetParameter(2)) << ";" << endl;
       hist_bsl1_l1l2_e->SetBinContent(nth+1, fabs(func->GetParameter(1)));
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

   } // end of nth loop ( 60 ~ 200 GeV )
  
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

   hist_bsl1_l1l2_e->SetMarkerColor(1);
   hist_bsl1_l1l2_e->SetMarkerStyle(20);
   hist_bsl1_l1l2_e->SetMarkerSize(1.);
   hist_bsl1_l1l2_e->Draw("hist p same");

   points << "}" << endl;
   points << "#endif" << endl;
   points.close();

   c2->Print("Fit_points.pdf");

   TFile *output = new TFile("fit.root", "RECREATE");
   h1->Write();
   hist_bsl1_l1l2_a->Write();
   hist_bsl1_l1l2_b->Write();
   hist_bsl1_l1l2_c->Write();
   hist_bsl1_l1l2_d->Write();
   hist_bsl1_l1l2_e->Write();
   
}
