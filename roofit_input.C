#define roofit_input_cxx
#include "roofit_input.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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

   
   Int_t nbins = 150; Float_t x1 = 0.; Float_t x2 = 150.;
   TH1F* hist_bsl1_l1l2 = new TH1F("hist_bsl1_l1l2",";P_{T} (GeV); #Delta#phi", nbins, x1, x2);
  
   ofstream Bsl1_l1l2_points;
   Bsl1_l1l2_points.open("./Bsl1_l1l2_points.h");

   Bsl1_l1l2_points << "#ifndef Bsl1_l1l2_points_h" << endl;
   Bsl1_l1l2_points << "#define Bsl1_l1l2_points_h" << endl;
   Bsl1_l1l2_points << endl;

   Bsl1_l1l2_points << "double gen_Pt[150] ={}, point[150] = {}; " << endl;
   Bsl1_l1l2_points << "int points = " << nbins << ";" << endl;

   Bsl1_l1l2_points << "void set_arrary(){" << endl;
   Bsl1_l1l2_points << endl;

   Float_t low_pt = 0.;
   Float_t high_pt = 1.;
   
   //for(Int_t nth = 0; nth < 150; nth++)
   for(Int_t nth = 0; nth < 50; nth++)
   {
       cout << "Process: (" << nth << "/150)" << endl;
       TString file_ = "./Fix_results/phi_Pt";
       TString ith_;

       ith_.Form("%d", nth);
       file_ = file_ + ith_ + ".pdf";

       low_pt = 0.+ nth;
       high_pt = low_pt + 1.;

       cout << "Low pT bound: " << low_pt << ", high pT bound: " << high_pt << endl;

       TTree *out_tree = new TTree("t","t");
       Float_t bsl1_l1l2_dphi_; 
       out_tree->Branch("bsl1_l1l2_Dphi", &bsl1_l1l2_dphi_, "bsl1_l1l2_Dphi/F");
       
       for (Long64_t jentry=0; jentry<nentries;jentry++) {
       //for (Long64_t jentry=0; jentry<10;jentry++) {
	   Long64_t ientry = LoadTree(jentry);
	   if (ientry < 0) break;
	   nb = fChain->GetEntry(jentry);   nbytes += nb;
	   // if (Cut(ientry) < 0) continue;

	   //cout << "Process: (" << jentry << "/" << nentries << ")" << endl;

	   int size = bsl1_l1l2_dphi->size();
	   int pT_size = genpT->size(); 

	   if( size == 0 ) continue;

	   //cout << "Size of dphi vector: " << size << ", size of pT vector: " << pT_size << endl;
	   //cout << "dPhi: " << bsl1_l1l2_dphi->at(0) << ", pT: " << genpT->at(0) << endl;

           Float_t gen_pT = genpT->at(0);
	   if( gen_pT > high_pt || gen_pT < low_pt ) continue;
	   cout << "generator level pT: " << gen_pT << endl;

           bsl1_l1l2_dphi_ = bsl1_l1l2_dphi->at(0);
	   out_tree->Fill();

       }// end of event loop
       cout << endl;

       float low_x = -0.1, high_x = 0.1;
       RooRealVar bsl1_l1l2_Dphi("bsl1_l1l2_Dphi","bsl1_l1l2_Dphi", low_x, high_x) ;

       RooDataSet ds1("ds1","ds1",RooArgSet(bsl1_l1l2_Dphi),Import(*out_tree)) ;

       RooPlot* frame1 = bsl1_l1l2_Dphi.frame(Title("bsl1_l1l2_Dphi")) ;
       ds1.plotOn(frame1) ;

       RooRealVar mean1("mean1","mean1", -.05, -.1, .1) ;
       RooRealVar sigma1("sigma1","sigma1", .05, .0001, 10.) ;
       RooGaussian gauss1("gauss1","gauss1",bsl1_l1l2_Dphi,mean1,sigma1) ;

       RooRealVar l1_a0("a0","a0",0.5,0.,1.) ;
       RooRealVar l1_a1("a1","a1",-0.2,0.,1.) ;
       RooChebychev l1_bkg("bkg","Background",bsl1_l1l2_Dphi,RooArgSet(l1_a0,l1_a1)) ;

       RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.98,0.95,1.) ;
       RooAddPdf sig("sig","Signal",RooArgList(gauss1, l1_bkg),sig1frac) ;

       sig.fitTo(ds1) ;
       sig.plotOn(frame1,LineColor(kCyan)) ;
       sig.plotOn(frame1, Components(gauss1), LineStyle(kDashed), LineColor(kRed)) ;
       sig.plotOn(frame1, Components(l1_bkg), LineStyle(kDashed), LineColor(kBlue)) ;

       hist_bsl1_l1l2->SetBinContent(nth+1, fabs(mean1.getVal()));

       Bsl1_l1l2_points << "gen_Pt[" << nth << "] = " << nth << "; point[" << nth << "] = " << fabs(mean1.getVal()) << ";" << endl;

       TCanvas* c = new TCanvas("Region of interest","Region of interest",800,800) ;

       c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;  

       c->SaveAs(file_);

       out_tree->Delete();
       c->Clear();
       delete c;
   
   } // end of nth loop
   
   TCanvas *c1 = new TCanvas("c1","c1",700,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(2); // axis width, default is 1
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.12);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.2);
   c1->SetGrid();
   c1->SetTicky(1);
   c1->SetTickx(1);

   hist_bsl1_l1l2->SetTitle("");
   hist_bsl1_l1l2->GetXaxis()->SetTitleOffset(1.3);
   hist_bsl1_l1l2->GetXaxis()->SetTitleSize(0.045);
   hist_bsl1_l1l2->GetXaxis()->SetNdivisions(505);
   hist_bsl1_l1l2->GetYaxis()->SetNdivisions(506);
   hist_bsl1_l1l2->GetXaxis()->SetLabelSize(0.05);
   hist_bsl1_l1l2->GetYaxis()->SetLabelSize(0.05);
   hist_bsl1_l1l2->GetXaxis()->SetRangeUser(x1,x2);
   hist_bsl1_l1l2->GetYaxis()->SetRangeUser(-0.01, 0.08);
   hist_bsl1_l1l2->GetYaxis()->SetTitleOffset(1.5);
   hist_bsl1_l1l2->GetYaxis()->SetTitleSize(0.050);

   hist_bsl1_l1l2->SetMarkerColor(1);
   hist_bsl1_l1l2->SetMarkerStyle(20);
   hist_bsl1_l1l2->SetMarkerSize(1.);
   hist_bsl1_l1l2->Draw("hist p");

   Bsl1_l1l2_points << "}" << endl;
   Bsl1_l1l2_points << "#endif" << endl;
   Bsl1_l1l2_points.close();

   c1->Print("Fix_results/Bsl1_l1l2_dphi.pdf");


}
