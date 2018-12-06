#include <TCanvas.h>
#include <TH2.h>
#include <TFile.h>

void draw()
{
    Int_t region = 1;

    TFile *input;
    if( region == 1 ) input = new TFile("Region1.root","OPEN");
    if( region == 2 ) input = new TFile("Region2.root","OPEN");
    if( region == 3 ) input = new TFile("Region3.root","OPEN");
    if( region == 4 ) input = new TFile("Region4.root","OPEN");
    if( region == 5 ) input = new TFile("Region5.root","OPEN");
    if( region == 6 ) input = new TFile("Region6.root","OPEN");
   
    /*
    TH1F *a1 = (TH1F*)input->Get("dRThi1");
    TH1F *a2 = (TH1F*)input->Get("dRThi2");
    TH1F *a3 = (TH1F*)input->Get("dRThi3");
    TH1F *a4 = (TH1F*)input->Get("dRThi4");

    TH1F *h1 = (TH1F*)input->Get("dRSec1");
    TH1F *h2 = (TH1F*)input->Get("dRSec2");
    TH1F *h3 = (TH1F*)input->Get("dRSec3");
    TH1F *h4 = (TH1F*)input->Get("dRSec4");
    
    TH1F *b1 = (TH1F*)input->Get("dRFir1");
    TH1F *b2 = (TH1F*)input->Get("dRFir2");
    TH1F *b3 = (TH1F*)input->Get("dRFir3");
    TH1F *b4 = (TH1F*)input->Get("dRFir4");
    */

    TH1F *a1 = (TH1F*)input->Get("dPhiThi1");
    TH1F *a2 = (TH1F*)input->Get("dPhiThi2");
    TH1F *a3 = (TH1F*)input->Get("dPhiThi3");
    TH1F *a4 = (TH1F*)input->Get("dPhiThi4");

    TH1F *h1 = (TH1F*)input->Get("dPhiSec1");
    TH1F *h2 = (TH1F*)input->Get("dPhiSec2");
    TH1F *h3 = (TH1F*)input->Get("dPhiSec3");
    TH1F *h4 = (TH1F*)input->Get("dPhiSec4");
    
    TH1F *b1 = (TH1F*)input->Get("dPhiFir1");
    TH1F *b2 = (TH1F*)input->Get("dPhiFir2");
    TH1F *b3 = (TH1F*)input->Get("dPhiFir3");
    TH1F *b4 = (TH1F*)input->Get("dPhiFir4");

    TCanvas *c1 = new TCanvas("c1","",1200,800);
    c1->Divide(2,2);
    c1->cd(1);
    h1->Draw();
    c1->cd(2);
    h2->Draw();
    c1->cd(3);
    h3->Draw();
    c1->cd(4);
    h4->Draw();
    
    TCanvas *c2 = new TCanvas("c2","",1200,800);
    c2->Divide(2,2);
    c2->cd(1);
    b1->Draw();
    c2->cd(2);
    b2->Draw();
    c2->cd(3);
    b3->Draw();
    c2->cd(4);
    b4->Draw();

    TCanvas *c3 = new TCanvas("c3","",1200,800);
    c3->Divide(2,2);
    c3->cd(1);
    a1->Draw();
    c3->cd(2);
    a2->Draw();
    c3->cd(3);
    a3->Draw();
    c3->cd(4);
    a4->Draw();

    TString filename1;
    TString filename2;
    TString filename3;
    
    /*
    if( region == 1 ) {
        filename1 = "dR02_r1.png";   
        filename2 = "dR03_r1.png";
        filename3 = "dR01_r1.png";
    }
    if( region == 2 ) {
        filename1 = "dR02_r2.png";
        filename2 = "dR03_r2.png";
        filename3 = "dR01_r2.png";
    }
    if( region == 3 ) {
        filename1 = "dR02_r3.png";
        filename2 = "dR03_r3.png";
        filename3 = "dR01_r3.png";
    }
    if( region == 4 ) {
        filename1 = "dR02_r4.png";
        filename2 = "dR03_r4.png";
        filename3 = "dR01_r4.png";
    }
    if( region == 5 ) {
        filename1 = "dR02_r5.png";
        filename2 = "dR03_r5.png";
        filename3 = "dR01_r5.png";
    }
    if( region == 6 ) {
        filename1 = "dR02_r6.png";
        filename2 = "dR03_r6.png";
        filename3 = "dR01_r6.png";
    }
    */
    
    if( region == 1 ) {
        filename1 = "dPhi02_r1.png";   
        filename2 = "dPhi03_r1.png";
        filename3 = "dPhi01_r1.png";
    }
    if( region == 2 ) {
        filename1 = "dPhi02_r2.png";
        filename2 = "dPhi03_r2.png";
        filename3 = "dPhi01_r2.png";
    }
    if( region == 3 ) {
        filename1 = "dPhi02_r3.png";
        filename2 = "dPhi03_r3.png";
        filename3 = "dPhi01_r3.png";
    }
    if( region == 4 ) {
        filename1 = "dPhi02_r4.png";
        filename2 = "dPhi03_r4.png";
        filename3 = "dPhi01_r4.png";
    }
    if( region == 5 ) {
        filename1 = "dPhi02_r5.png";
        filename2 = "dPhi03_r5.png";
        filename3 = "dPhi01_r5.png";
    }
    if( region == 6 ) {
        filename1 = "dPhi02_r6.png";
        filename2 = "dPhi03_r6.png";
        filename3 = "dPhi01_r6.png";
    }

    c1->Print(filename1);
    c2->Print(filename2);
    c3->Print(filename3);
     
    c1->Close();
    c2->Close();
    c3->Close();
}
