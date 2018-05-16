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

void roofit_input::Loop(float low_et = 10., float high_et = 11.)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
/*
   TTree* out_tree1 = new TTree("roi","roi");
   TTree* out_tree2 = new TTree("roi","roi");
   TTree* out_tree3 = new TTree("roi","roi");
   TTree* out_tree4 = new TTree("roi","roi");

   float bpix_l1_dphi_; out_tree1->Branch("bpix_l1_dphi", &bpix_l1_dphi_, "bpix_l1_dphi/F");
   float bpix_l2_dphi_; out_tree2->Branch("bpix_l2_dphi", &bpix_l2_dphi_, "bpix_l2_dphi/F");
   float bpix_l3_dphi_; out_tree3->Branch("bpix_l3_dphi", &bpix_l3_dphi_, "bpix_l3_dphi/F");
   float bpix_l4_dphi_; out_tree4->Branch("bpix_l4_dphi", &bpix_l4_dphi_, "bpix_l4_dphi/F");
*/
int nbins = 100; float x1 = 10.; float x2 = 110.;
TH1F* BPix1_ROI_up = new TH1F("BPix1_ROI_up",";E_{T} (GeV); dPhi",nbins,x1,x2);
TH1F* BPix1_ROI_down = new TH1F("BPix1_ROI_down",";E_{T} (GeV); dPhi",nbins,x1,x2);

TH1F* BPix2_ROI_up = new TH1F("BPix2_ROI_up",";E_{T} (GeV); dPhi",nbins,x1,x2);
TH1F* BPix2_ROI_down = new TH1F("BPix2_ROI_down",";E_{T} (GeV); dPhi",nbins,x1,x2);

TH1F* BPix3_ROI_up = new TH1F("BPix3_ROI_up",";E_{T} (GeV); dPhi",nbins,x1,x2);
TH1F* BPix3_ROI_down = new TH1F("BPix3_ROI_down",";E_{T} (GeV); dPhi",nbins,x1,x2);

TH1F* BPix4_ROI_up = new TH1F("BPix4_ROI_up",";E_{T} (GeV); dPhi",nbins,x1,x2);
TH1F* BPix4_ROI_down = new TH1F("BPix4_ROI_down",";E_{T} (GeV); dPhi",nbins,x1,x2);

// for the RoI of the first barrel layer
ofstream BPix1_ROI_points;
BPix1_ROI_points.open("./BPix1_ROI_points.h");

BPix1_ROI_points << "#ifndef BPix1_ROI_points_h" << endl;
BPix1_ROI_points << "#define BPix1_ROI_points_h" << endl;
BPix1_ROI_points << endl;

BPix1_ROI_points << "double EG_Et[100] ={}, up_bound[100] = {}; " << endl;
BPix1_ROI_points << "double low_bound[100] = {}; " << endl;
BPix1_ROI_points << "int points = " << nbins << ";" << endl;

BPix1_ROI_points << "void set_arrary(){" << endl;
BPix1_ROI_points << endl;

// for the RoI of the second barrel layer
ofstream BPix2_ROI_points;
BPix2_ROI_points.open("./BPix2_ROI_points.h");

BPix2_ROI_points << "#ifndef BPix2_ROI_points_h" << endl;
BPix2_ROI_points << "#define BPix2_ROI_points_h" << endl;
BPix2_ROI_points << endl;

BPix2_ROI_points << "double EG_Et[100] ={}, up_bound[100] = {}; " << endl;
BPix2_ROI_points << "double low_bound[100] = {}; " << endl;
BPix2_ROI_points << "int points = " << nbins << ";" << endl;

BPix2_ROI_points << "void set_arrary(){" << endl;
BPix2_ROI_points << endl;

// for the RoI of the first barrel layer
ofstream BPix3_ROI_points;
BPix3_ROI_points.open("./BPix3_ROI_points.h");

BPix3_ROI_points << "#ifndef BPix3_ROI_points_h" << endl;
BPix3_ROI_points << "#define BPix3_ROI_points_h" << endl;
BPix3_ROI_points << endl;

BPix3_ROI_points << "double EG_Et[100] ={}, up_bound[100] = {}; " << endl;
BPix3_ROI_points << "double low_bound[100] = {}; " << endl;
BPix3_ROI_points << "int points = " << nbins << ";" << endl;

BPix3_ROI_points << "void set_arrary(){" << endl;
BPix3_ROI_points << endl;

// for the RoI of the first barrel layer
ofstream BPix4_ROI_points;
BPix4_ROI_points.open("./BPix4_ROI_points.h");

BPix4_ROI_points << "#ifndef BPix4_ROI_points_h" << endl;
BPix4_ROI_points << "#define BPix4_ROI_points_h" << endl;
BPix4_ROI_points << endl;

BPix4_ROI_points << "double EG_Et[100] ={}, up_bound[100] = {}; " << endl;
BPix4_ROI_points << "double low_bound[100] = {}; " << endl;
BPix4_ROI_points << "int points = " << nbins << ";" << endl;

BPix4_ROI_points << "void set_arrary(){" << endl;
BPix4_ROI_points << endl;


for( int nth = 0; nth < 100; nth++){
   TString file_ = "./ROI_Et";
   TString ith_;

   ith_.Form("%d", nth + 10);
   file_ = file_ + ith_ +".pdf";

   low_et = 10. + nth;
   high_et = low_et + 1.;

   TTree* out_tree1 = new TTree("roi","roi");
   TTree* out_tree2 = new TTree("roi","roi");
   TTree* out_tree3 = new TTree("roi","roi");
   TTree* out_tree4 = new TTree("roi","roi");

   float bpix_l1_dphi_; out_tree1->Branch("bpix_l1_dphi", &bpix_l1_dphi_, "bpix_l1_dphi/F");
   float bpix_l2_dphi_; out_tree2->Branch("bpix_l2_dphi", &bpix_l2_dphi_, "bpix_l2_dphi/F");
   float bpix_l3_dphi_; out_tree3->Branch("bpix_l3_dphi", &bpix_l3_dphi_, "bpix_l3_dphi/F");
   float bpix_l4_dphi_; out_tree4->Branch("bpix_l4_dphi", &bpix_l4_dphi_, "bpix_l4_dphi/F");


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      // test for bpix1
      if(closest_egEt > low_et && closest_egEt < high_et && fabs(closest_egEta) < 1.3 && closest_eg_dr < 0.1){
        int size_bpix_l1_hit = bpix_l1_dphi_roi->size();
        int size_bpix_l2_hit = bpix_l2_dphi_roi->size();
        int size_bpix_l3_hit = bpix_l3_dphi_roi->size();
        int size_bpix_l4_hit = bpix_l4_dphi_roi->size();

        for(int i=0; i < size_bpix_l1_hit; i++){
           bpix_l1_dphi_ = bpix_l1_dphi_roi->at(i);
           out_tree1->Fill();
        }

        for(int i=0; i < size_bpix_l2_hit; i++){
           bpix_l2_dphi_ = bpix_l2_dphi_roi->at(i);
           out_tree2->Fill();
        }

        for(int i=0; i < size_bpix_l3_hit; i++){
           bpix_l3_dphi_ = bpix_l3_dphi_roi->at(i);
           out_tree3->Fill();
        }

        for(int i=0; i < size_bpix_l4_hit; i++){
           bpix_l4_dphi_ = bpix_l4_dphi_roi->at(i);
           out_tree4->Fill();
        }

      }

   }// end of event loop

   float low_x = -0.25, high_x = 0.25;
   RooRealVar bpix_l1_dphi("bpix_l1_dphi","bpix_l1_dphi", low_x, high_x) ;
   RooRealVar bpix_l2_dphi("bpix_l2_dphi","bpix_l2_dphi", low_x, high_x) ;
   RooRealVar bpix_l3_dphi("bpix_l3_dphi","bpix_l3_dphi", low_x, high_x) ;
   RooRealVar bpix_l4_dphi("bpix_l4_dphi","bpix_l4_dphi", low_x, high_x) ;

   RooDataSet ds1("ds1","ds1",RooArgSet(bpix_l1_dphi),Import(*out_tree1)) ;
   RooDataSet ds2("ds2","ds2",RooArgSet(bpix_l2_dphi),Import(*out_tree2)) ;
   RooDataSet ds3("ds3","ds3",RooArgSet(bpix_l3_dphi),Import(*out_tree3)) ;
   RooDataSet ds4("ds4","ds4",RooArgSet(bpix_l4_dphi),Import(*out_tree4)) ;

   RooPlot* frame1 = bpix_l1_dphi.frame(Title("bpix_l1_dphi")) ;
   ds1.plotOn(frame1) ;

   RooPlot* frame2 = bpix_l2_dphi.frame(Title("bpix_l2_dphi")) ;
   ds2.plotOn(frame2) ;

   RooPlot* frame3 = bpix_l3_dphi.frame(Title("bpix_l3_dphi")) ;
   ds3.plotOn(frame3) ;

   RooPlot* frame4 = bpix_l4_dphi.frame(Title("bpix_l4_dphi")) ;
   ds4.plotOn(frame4) ;

   RooRealVar mean1("mean1","mean1", -.05, -.2, .2) ;
   RooRealVar sigma1("sigma1","sigma1", .05, .0001, 10.) ;
   RooGaussian gauss1("gauss1","gauss1",bpix_l1_dphi,mean1,sigma1) ;

   RooRealVar l1_a0("a0","a0",0.5,0.,1.) ;
   RooRealVar l1_a1("a1","a1",-0.2,0.,1.) ;
   RooChebychev l1_bkg("bkg","Background",bpix_l1_dphi,RooArgSet(l1_a0,l1_a1)) ;

   RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.98,0.95,1.) ;
   RooAddPdf sig("sig","Signal",RooArgList(gauss1, l1_bkg),sig1frac) ;

   sig.fitTo(ds1) ;
   sig.plotOn(frame1,LineColor(kCyan)) ;
   sig.plotOn(frame1, Components(gauss1), LineStyle(kDashed), LineColor(kRed)) ;
   sig.plotOn(frame1, Components(l1_bkg), LineStyle(kDashed), LineColor(kBlue)) ;

   BPix1_ROI_up->SetBinContent(nth+1, mean1.getVal() + 3.*sigma1.getVal());
   BPix1_ROI_down->SetBinContent(nth+1, mean1.getVal() - 3.*sigma1.getVal());

   BPix1_ROI_points << "EG_Et[" << nth << "] = " << nth+10 << "; up_bound[" << nth << "] = " << mean1.getVal() + 3.*sigma1.getVal() << ";" 
                    << " low_bound[" << nth << "] = " << mean1.getVal() - 3.*sigma1.getVal() << ";" << endl;


   RooRealVar mean2("mean2","mean2", -.05, -.2, .2) ;
   RooRealVar sigma2("sigma2","sigma2", .05, .0001, 10.) ;
   RooGaussian gauss2("gauss2","gauss2",bpix_l2_dphi,mean2,sigma2) ;

   RooRealVar l2_a0("a0","a0",0.5,0.,1.) ;
   RooRealVar l2_a1("a1","a1",-0.2,0.,1.) ;
   RooChebychev l2_bkg("bkg","Background",bpix_l2_dphi,RooArgSet(l2_a0,l2_a1)) ;

   RooRealVar sig2frac("sig2frac","fraction of component 1 in signal",0.98,0.95,1.) ;
   RooAddPdf sig2("sig2","Signal",RooArgList(gauss2, l2_bkg),sig2frac) ;

   sig2.fitTo(ds2) ;
   sig2.plotOn(frame2,LineColor(kCyan)) ;
   sig2.plotOn(frame2, Components(gauss2), LineStyle(kDashed), LineColor(kRed)) ;
   sig2.plotOn(frame2, Components(l2_bkg), LineStyle(kDashed), LineColor(kBlue)) ;

   BPix2_ROI_up->SetBinContent(nth+1, mean2.getVal() + 3.*sigma2.getVal());
   BPix2_ROI_down->SetBinContent(nth+1, mean2.getVal() - 3.*sigma2.getVal());

   BPix2_ROI_points << "EG_Et[" << nth << "] = " << nth+10 << "; up_bound[" << nth << "] = " << mean2.getVal() + 3.*sigma2.getVal() << ";" 
                    << " low_bound[" << nth << "] = " << mean2.getVal() - 3.*sigma2.getVal() << ";" << endl;

   RooRealVar mean3("mean3","mean3", -.05, -.2, .2) ;
   RooRealVar sigma3("sigma3","sigma3", .05, .0001, 10.) ;
   RooGaussian gauss3("gauss3","gauss3",bpix_l3_dphi,mean3,sigma3) ;

   RooRealVar l3_a0("a0","a0",0.5,0.,1.) ;
   RooRealVar l3_a1("a1","a1",-0.2,0.,1.) ;
   RooChebychev l3_bkg("bkg","Background",bpix_l3_dphi,RooArgSet(l3_a0,l3_a1)) ;

   RooRealVar sig3frac("sig3frac","fraction of component 1 in signal",0.98,0.95,1.) ;
   RooAddPdf sig3("sig3","Signal",RooArgList(gauss3, l3_bkg),sig3frac) ;

   sig3.fitTo(ds3) ;
   sig3.plotOn(frame3,LineColor(kCyan)) ;
   sig3.plotOn(frame3, Components(gauss3), LineStyle(kDashed), LineColor(kRed)) ;
   sig3.plotOn(frame3, Components(l3_bkg), LineStyle(kDashed), LineColor(kBlue)) ;

   BPix3_ROI_up->SetBinContent(nth+1, mean3.getVal() + 3.*sigma3.getVal());
   BPix3_ROI_down->SetBinContent(nth+1, mean3.getVal() - 3.*sigma3.getVal());

   BPix3_ROI_points << "EG_Et[" << nth << "] = " << nth+10 << "; up_bound[" << nth << "] = " << mean3.getVal() + 3.*sigma3.getVal() << ";" 
                    << " low_bound[" << nth << "] = " << mean3.getVal() - 3.*sigma3.getVal() << ";" << endl;

   RooRealVar mean4("mean4","mean4", -.05, -.2, .2) ;
   RooRealVar sigma4("sigma4","sigma4", .05, .0001, 10.) ;
   RooGaussian gauss4("gauss4","gauss4",bpix_l4_dphi,mean4,sigma4) ;

   RooRealVar l4_a0("a0","a0",0.5,0.,1.) ;
   RooRealVar l4_a1("a1","a1",-0.2,0.,1.) ;
   RooChebychev l4_bkg("bkg","Background",bpix_l4_dphi,RooArgSet(l4_a0,l3_a1)) ;

   RooRealVar sig4frac("sig4frac","fraction of component 1 in signal",0.996,0.99,1.) ;
   RooAddPdf sig4("sig4","Signal",RooArgList(gauss4, l4_bkg),sig4frac) ;

   sig4.fitTo(ds4) ;
   sig4.plotOn(frame4,LineColor(kCyan)) ;
   sig4.plotOn(frame4, Components(gauss4), LineStyle(kDashed), LineColor(kRed)) ;
   sig4.plotOn(frame4, Components(l4_bkg), LineStyle(kDashed), LineColor(kBlue)) ;

   BPix4_ROI_up->SetBinContent(nth+1, mean4.getVal() + 3.*sigma4.getVal());
   BPix4_ROI_down->SetBinContent(nth+1, mean4.getVal() - 3.*sigma4.getVal());

   BPix4_ROI_points << "EG_Et[" << nth << "] = " << nth+10 << "; up_bound[" << nth << "] = " << mean4.getVal() + 3.*sigma4.getVal() << ";" 
                    << " low_bound[" << nth << "] = " << mean4.getVal() - 3.*sigma4.getVal() << ";" << endl;

   TCanvas* c = new TCanvas("Region of interest","Region of interest",800,800) ;
   c->Divide(2,2) ;

   c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;  
   c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;  
   c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;  
   c->cd(4) ; gPad->SetLeftMargin(0.15) ; frame4->GetYaxis()->SetTitleOffset(1.4) ; frame4->Draw() ;  
   
   c->SaveAs(file_);

   out_tree1->Delete();
   out_tree2->Delete();
   out_tree3->Delete();
   out_tree4->Delete();
   c->Clear();
   delete c;

}// end of et loop

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

   BPix1_ROI_up->SetTitle("");
   BPix1_ROI_up->GetXaxis()->SetTitleOffset(1.3);
   BPix1_ROI_up->GetXaxis()->SetTitleSize(0.045);
   BPix1_ROI_up->GetXaxis()->SetNdivisions(505);
   BPix1_ROI_up->GetYaxis()->SetNdivisions(506);
   BPix1_ROI_up->GetXaxis()->SetLabelSize(0.05);
   BPix1_ROI_up->GetYaxis()->SetLabelSize(0.05);
   BPix1_ROI_up->GetXaxis()->SetRangeUser(x1,x2);
   BPix1_ROI_up->GetYaxis()->SetRangeUser(-0.25, 0.25);
   BPix1_ROI_up->GetYaxis()->SetTitleOffset(1.5);
   BPix1_ROI_up->GetYaxis()->SetTitleSize(0.050);

   BPix1_ROI_up->SetMarkerColor(1);
   BPix1_ROI_up->SetMarkerStyle(20);
   BPix1_ROI_up->SetMarkerSize(1.);
   BPix1_ROI_up->Draw("hist p");
   BPix1_ROI_down->SetMarkerColor(2);
   BPix1_ROI_down->SetMarkerStyle(20);
   BPix1_ROI_down->SetMarkerSize(1.);
   BPix1_ROI_down->Draw("hist same p");

   BPix1_ROI_points  << "}" << endl;
   BPix1_ROI_points  << "#endif" << endl;
   BPix1_ROI_points .close();

   c1->Print("BPix1_ROI.pdf");

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
  
   c2->cd();

   BPix2_ROI_up->SetTitle("");
   BPix2_ROI_up->GetXaxis()->SetTitleOffset(1.3);
   BPix2_ROI_up->GetXaxis()->SetTitleSize(0.045);
   BPix2_ROI_up->GetXaxis()->SetNdivisions(505);
   BPix2_ROI_up->GetYaxis()->SetNdivisions(506);
   BPix2_ROI_up->GetXaxis()->SetLabelSize(0.05);
   BPix2_ROI_up->GetYaxis()->SetLabelSize(0.05);
   BPix2_ROI_up->GetXaxis()->SetRangeUser(x1,x2);
   BPix2_ROI_up->GetYaxis()->SetRangeUser(-0.25, 0.25);
   BPix2_ROI_up->GetYaxis()->SetTitleOffset(1.5);
   BPix2_ROI_up->GetYaxis()->SetTitleSize(0.050);

   BPix2_ROI_up->SetMarkerColor(1);
   BPix2_ROI_up->SetMarkerStyle(20);
   BPix2_ROI_up->SetMarkerSize(1.);
   BPix2_ROI_up->Draw("hist p");
   BPix2_ROI_down->SetMarkerColor(2);
   BPix2_ROI_down->SetMarkerStyle(20);
   BPix2_ROI_down->SetMarkerSize(1.);
   BPix2_ROI_down->Draw("hist same p");

   BPix2_ROI_points  << "}" << endl;
   BPix2_ROI_points  << "#endif" << endl;
   BPix2_ROI_points .close();

   c2->Print("BPix2_ROI.pdf");

   TCanvas *c3 = new TCanvas("c3","c3",700,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(2); // axis width, default is 1
   c3->SetTopMargin(0.05);
   c3->SetBottomMargin(0.12);
   c3->SetRightMargin(0.03);
   c3->SetLeftMargin(0.2);
   c3->SetGrid();
   c3->SetTicky(1);
   c3->SetTickx(1);
  
   c3->cd();

   BPix3_ROI_up->SetTitle("");
   BPix3_ROI_up->GetXaxis()->SetTitleOffset(1.3);
   BPix3_ROI_up->GetXaxis()->SetTitleSize(0.045);
   BPix3_ROI_up->GetXaxis()->SetNdivisions(505);
   BPix3_ROI_up->GetYaxis()->SetNdivisions(506);
   BPix3_ROI_up->GetXaxis()->SetLabelSize(0.05);
   BPix3_ROI_up->GetYaxis()->SetLabelSize(0.05);
   BPix3_ROI_up->GetXaxis()->SetRangeUser(x1,x2);
   BPix3_ROI_up->GetYaxis()->SetRangeUser(-0.25, 0.25);
   BPix3_ROI_up->GetYaxis()->SetTitleOffset(1.5);
   BPix3_ROI_up->GetYaxis()->SetTitleSize(0.050);

   BPix3_ROI_up->SetMarkerColor(1);
   BPix3_ROI_up->SetMarkerStyle(20);
   BPix3_ROI_up->SetMarkerSize(1.);
   BPix3_ROI_up->Draw("hist p");
   BPix3_ROI_down->SetMarkerColor(2);
   BPix3_ROI_down->SetMarkerStyle(20);
   BPix3_ROI_down->SetMarkerSize(1.);
   BPix3_ROI_down->Draw("hist same p");

   BPix3_ROI_points  << "}" << endl;
   BPix3_ROI_points  << "#endif" << endl;
   BPix3_ROI_points .close();

   c3->Print("BPix3_ROI.pdf");

   TCanvas *c4 = new TCanvas("c4","c4",700,700);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(2); // axis width, default is 1
   c4->SetTopMargin(0.05);
   c4->SetBottomMargin(0.12);
   c4->SetRightMargin(0.03);
   c4->SetLeftMargin(0.2);
   c4->SetGrid();
   c4->SetTicky(1);
   c4->SetTickx(1);
  
   c4->cd();

   BPix4_ROI_up->SetTitle("");
   BPix4_ROI_up->GetXaxis()->SetTitleOffset(1.3);
   BPix4_ROI_up->GetXaxis()->SetTitleSize(0.045);
   BPix4_ROI_up->GetXaxis()->SetNdivisions(505);
   BPix4_ROI_up->GetYaxis()->SetNdivisions(506);
   BPix4_ROI_up->GetXaxis()->SetLabelSize(0.05);
   BPix4_ROI_up->GetYaxis()->SetLabelSize(0.05);
   BPix4_ROI_up->GetXaxis()->SetRangeUser(x1,x2);
   BPix4_ROI_up->GetYaxis()->SetRangeUser(-0.25, 0.25);
   BPix4_ROI_up->GetYaxis()->SetTitleOffset(1.5);
   BPix4_ROI_up->GetYaxis()->SetTitleSize(0.050);

   BPix4_ROI_up->SetMarkerColor(1);
   BPix4_ROI_up->SetMarkerStyle(20);
   BPix4_ROI_up->SetMarkerSize(1.);
   BPix4_ROI_up->Draw("hist p");
   BPix4_ROI_down->SetMarkerColor(2);
   BPix4_ROI_down->SetMarkerStyle(20);
   BPix4_ROI_down->SetMarkerSize(1.);
   BPix4_ROI_down->Draw("hist same p");

   BPix4_ROI_points  << "}" << endl;
   BPix4_ROI_points  << "#endif" << endl;
   BPix4_ROI_points .close();

   c4->Print("BPix4_ROI.pdf");

}
