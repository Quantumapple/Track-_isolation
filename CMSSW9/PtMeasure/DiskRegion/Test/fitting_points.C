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

  TString fit_name = "Bsd4_d4d6";

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
  axis->GetYaxis()->SetRangeUser(0, 100);
  axis->GetYaxis()->CenterTitle(true);
  axis->GetYaxis()->SetTitleOffset(1.2);
  axis->GetYaxis()->SetTitleSize(0.040);
  axis->GetYaxis()->SetTitle("Generator level p_{T} (GeV)");
  axis->Draw();


  TGraphErrors *a = new TGraphErrors(points, points_y, points_x, points_ex, 0);
  a->SetMarkerStyle(20);
  a->SetMarkerColor(2);
  a->SetMarkerSize(1);
  a->Draw("PE"); 
  
  TF1 *fit_func = new TF1("func","( [0]*(1./x) + [1] )", 0, 1);
  
  a->Fit(fit_func,"0W"); 
  fit_func->Draw("lsame");
  fit_func->SetLineColor(1);
  
  ofstream fit_result;
  fit_result.open("parameters.txt");    

  fit_result << endl;  
  for( int i=0; i < 2; i++){
      fit_result << "p[" << i << "] = " << fit_func->GetParameter(i) << ";" << endl;
  }   
  fit_result << "p[0]*(1./x) + p[1];" << endl;
  fit_result << endl;  
  fit_result.close();

  TFile *output = new TFile("Pt_fit.root","RECREATE");
  axis->Write();
  a->Write();
  fit_func->Write();


}
