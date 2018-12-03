#include <iostream>
#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>

void draw()
{
    Int_t eta_r = 1;
    
    TFile *input1 = new TFile("../L1PixEle-eff/Results/dZ-Distribution/SE_PU200_Region1/output/Tree/Tree_SE_PU200.root","open");
    //TFile *input1 = new TFile("../L1PixEle-eff/Results/dZ-Distribution/SE_PU200_Region2/output/Tree/Tree_SE_PU200.root","open");
    //TFile *input1 = new TFile("../L1PixEle-eff/Results/dZ-Distribution/SE_PU200_Region3/output/Tree/Tree_SE_PU200.root","open");
    //TFile *input1 = new TFile("../L1PixEle-eff/Results/dZ-Distribution/SE_PU200_Region4/output/Tree/Tree_SE_PU200.root","open");
    //TFile *input1 = new TFile("../L1PixEle-eff/Results/dZ-Distribution/SE_PU200_Region5/output/Tree/Tree_SE_PU200.root","open");
    //TFile *input1 = new TFile("../L1PixEle-eff/Results/dZ-Distribution/SE_PU200_Region6/output/Tree/Tree_SE_PU200.root","open");
   
    TH1F *dZdist = (TH1F*)input1->Get("h3");
    dZdist->SetTitle("");

    dZdist->SetName("htemp");
    //if( eta_r == 1 ) dZdist->SetName("#Delta z dist in Region1");
    //if( eta_r == 2 ) dZdist->SetName("#Delta z dist in Region2");
    //if( eta_r == 3 ) dZdist->SetName("#Delta z dist in Region3");
    //if( eta_r == 4 ) dZdist->SetName("#Delta z dist in Region4");
    //if( eta_r == 5 ) dZdist->SetName("#Delta z dist in Region5");
    //if( eta_r == 6 ) dZdist->SetName("#Delta z dist in Region6");

    TCanvas *c1 = new TCanvas("c1","",800,700);
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.12);

    dZdist->GetXaxis()->SetTitle("z_{vtx} - z_{gen} (cm)");
    dZdist->GetXaxis()->CenterTitle(true);
    dZdist->GetXaxis()->SetTitleSize(0.05);
    
    if( eta_r == 1 ) dZdist->GetXaxis()->SetRangeUser(-0.04,0.04); 
    if( eta_r == 2 ) dZdist->GetXaxis()->SetRangeUser(-0.04,0.04); 
    if( eta_r == 3 ) dZdist->GetXaxis()->SetRangeUser(-0.06,0.06); 
    if( eta_r == 4 ) dZdist->GetXaxis()->SetRangeUser(-0.1,0.1); 
    if( eta_r == 5 ) dZdist->GetXaxis()->SetRangeUser(-0.2,0.2); 
    if( eta_r == 6 ) dZdist->GetXaxis()->SetRangeUser(-0.2,0.2); 

    dZdist->Draw();

    if( eta_r == 1 ) c1->Print("dZ-dist-Region1.png"); 
    if( eta_r == 2 ) c1->Print("dZ-dist-Region2.png"); 
    if( eta_r == 3 ) c1->Print("dZ-dist-Region3.png");
    if( eta_r == 4 ) c1->Print("dZ-dist-Region4.png");
    if( eta_r == 5 ) c1->Print("dZ-dist-Region5.png");
    if( eta_r == 6 ) c1->Print("dZ-dist-Region6.png");

    //gStyle->SetOptStat(0);


}
