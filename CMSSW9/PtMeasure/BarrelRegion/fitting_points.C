#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include "TKey.h"
#include "TMacro.h"
#include <string>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include "./fitting_points.h"

using namespace std;

void fitting_points(){

  set_arrary();

  Int_t nth = 0;;

  TString fit_name = "Bsl3_l3l4";

  double *points_x, *points_y, *points_ex; 
  points_x = new double[points]; 
  points_y = new double[points];
  points_ex = new double[points];
  
  cout << "points: " << points << endl;
  for( Int_t i = 0; i < points; i++)
  {
     points_x[i] = gen_Pt[i]; 
     points_y[i] = point[i];
     points_ex[i] = error[i];
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // From Dr.Moon //
  gStyle->SetOptStat(0);
  
  //TCanvas *cc = new TCanvas("cc","cc",800,600);
  //cc->cd();
  //cc->SetLogy();
  
  TCanvas *c1 = new TCanvas("c1","",800,800);
  c1->cd();
  
  //TH2F *axis = new TH2F("","",150,0,150,4000,-0.30,0.20);
  TH2F *axis = new TH2F("","",5000,0,0.2,155,0,155);
  char signal_title[200];
  
  axis->SetTitle(fit_name);
  axis->SetTitleSize(0.050);
  axis->GetXaxis()->CenterTitle(true);
  axis->GetXaxis()->SetTitleOffset(1.1);
  axis->GetXaxis()->SetTitleSize(0.040);
  axis->GetXaxis()->SetTitle("#Delta#phi");
  //axis->GetXaxis()->SetNdivisions(521);
  //axis->GetYaxis()->SetNdivisions(517);
  //axis->GetYaxis()->SetRangeUser(0.00001, 1.);
  axis->GetXaxis()->SetRangeUser(0.00001, 1.);
  axis->GetYaxis()->SetRangeUser(0, 155);
  axis->GetYaxis()->CenterTitle(true);
  axis->GetYaxis()->SetTitleOffset(1.2);
  axis->GetYaxis()->SetTitleSize(0.040);
  axis->GetYaxis()->SetTitle("Generator level p_{T} (GeV)");
  axis->Draw();


  //TGraphErrors *a = new TGraphErrors(points, points_y, points_x,points_ey, 0);
  TGraphErrors *a = new TGraphErrors(points, points_y, points_x, points_ex, 0);
  a->SetMarkerStyle(20);
  a->SetMarkerColor(2);
  a->SetMarkerSize(1);
  a->Draw("PE"); 
  
  TF1 *fit_func = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 0, 1);
  
  fit_func->SetParameter(1,0.001);
  
  //fit_func->SetParLimits(0,3.5e-9,3.8e-9);
  //fit_func->SetParameter(0,3.6e-9);

  //fit_func->SetParLimits(1,0.0002,0.0003);
  //fit_func->SetParameter(1,0.00003);

  //fit_func->SetParLimits(2,-0.57,-0.55);
  //fit_func->SetParameter(2,-0.56);

  //fit_func->SetParLimits(3,0.190,0.200);
  //fit_func->SetParameter(3,0.199);

  //fit_func->SetParLimits(4,5.3,5.4);
  //fit_func->SetParameter(4,5.4);

  a->Fit(fit_func,"0W"); 
  fit_func->Draw("lsame");
  fit_func->SetLineColor(1);
  
  ofstream fit_result;
  fit_result.open("parameters.txt");    

  fit_result << endl;  
  for( int i=0; i < 5; i++){
      fit_result << "p[" << i << "] = " << fit_func->GetParameter(i) << ";" << endl;
  }   
  fit_result << "p[0]*pow(x,0) + p[1]*pow(x,p[2])*exp(-pow(x,p[3])+p[4]);" << endl;
  fit_result << endl;  
  fit_result.close();

  //c1->SaveAs("./SW_linear.pdf");

  //cc->Update();
  
  //cc->SaveAs("./SW.pdf");
  
  TFile *output = new TFile("Pt_fit.root","RECREATE");
  axis->Write();
  a->Write();
  fit_func->Write();


}
