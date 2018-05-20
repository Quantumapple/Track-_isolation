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
#include "./Bsl1_l1l2_points.h"

using namespace std;

void Bsl1_l1l2_points(){

  set_arrary();

  Int_t nth = 0;;

  TString fit_name = "Bsl1_l1l2";

  double *points_x, *points_y; 
  points_x = new double[points]; 
  points_y = new double[points]; 
  
  cout << "points: " << points << endl;
  for( Int_t i = 0; i < points; i++)
  {
     points_x[i] = gen_Pt[i]; 
     points_y[i] = point[i];
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // From Dr.Moon //
  gStyle->SetOptStat(0);
  
  TCanvas *cc = new TCanvas("cc","cc",800,600);
  cc->cd();
  cc->SetLogy();
  
  TH2F *axis = new TH2F("","",150,0,150,4000,-0.30,0.20);
  char signal_title[200];
  
  axis->SetTitle(fit_name);
  axis->SetTitleSize(0.050);
  axis->GetXaxis()->SetTitle("Generator level p_{T} (GeV)");
  axis->GetXaxis()->CenterTitle(true);
  axis->GetXaxis()->SetTitleOffset(1.1);
  axis->GetXaxis()->SetTitleSize(0.040);
  axis->GetXaxis()->SetNdivisions(521);
  axis->GetYaxis()->SetNdivisions(517);
  axis->GetYaxis()->SetRangeUser(0.00001, 1.);
  axis->GetXaxis()->SetRangeUser(0, 155);
  axis->GetYaxis()->CenterTitle(true);
  axis->GetYaxis()->SetTitleOffset(1.2);
  axis->GetYaxis()->SetTitleSize(0.040);
  axis->GetYaxis()->SetTitle("#Delta#phi");
  axis->Draw();

  
  TGraphErrors *a = new TGraphErrors(points, points_x, points_y,0,0);
  a->SetMarkerStyle(20);
  a->SetMarkerColor(2);
  a->SetMarkerSize(1);
  a->Draw("P"); 
  
  TF1 *fit_func = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 0, 150);
  
  fit_func->SetParameter(1,0.01);
 
  a->Fit(fit_func,"0W"); 
  fit_func->Draw("lsame");
  fit_func->SetLineColor(1);
  
  /*
  double p[10]={0}, p2[10]={0};
  fit_func->GetParameters(p);
  fit_func2->GetParameters(p2);
  cout<<"\n/// Set Parameter before fit"<<endl;
  for(int k=0; k<5; k++) if(p[k]!=0)  cout<<"p["<<k<<"] = "<<p[k] <<";"<<endl;
  cout<<endl;
  for(int k=0; k<5; k++) if(p2[k]!=0) cout<<"p["<<k<<"] = "<<p2[k]<<";"<<endl;
  cout<<endl;

  ofstream fit_result;
  fit_result.open(roi_name + ".txt");    

  fit_result << endl;  
  fit_result << "if( i == " << nth_roi << " && up_down == 0){" <<endl;
  for( int i=0; i < 5; i++){
      fit_result << "p[" << i << "] = " << fit_func->GetParameter(i) << ";" << endl;
  }   
  fit_result << "return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);" << endl;
  fit_result << "}" << endl; 

  fit_result << "if( i == " << nth_roi << " && up_down == 1){" << endl;
  for( int i=0; i < 5; i++){
      fit_result << "p[" << i << "] = " << fit_func2->GetParameter(i) << ";" << endl;
  }   
  fit_result << "return p[0]*pow(eget,0) + p[1]*pow(eget,p[2])*exp(-pow(eget,p[3])+p[4]);" << endl;
  fit_result << "}" << endl; 
  fit_result << endl;  
  fit_result.close();
  */

  cc->Update();
  
  cc->SaveAs("./SW.pdf");

}
