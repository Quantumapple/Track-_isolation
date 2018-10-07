#include <TH1.h>
#include <TFile.h>

void fit_Region1()
{
    TFile *input = new TFile("dAngleEtaRegion1.root","OPEN");

    TH1F *hist_pp[12];
    hist_pp[0] =  (TH1F*)input->Get("dEtaCase1pp1");
    hist_pp[1] =  (TH1F*)input->Get("dEtaCase1pp2");
    hist_pp[2] =  (TH1F*)input->Get("dEtaCase1pp3");
    hist_pp[3] =  (TH1F*)input->Get("dEtaCase2pp1");
    hist_pp[4] =  (TH1F*)input->Get("dEtaCase2pp2");
    hist_pp[5] =  (TH1F*)input->Get("dEtaCase2pp3");
    hist_pp[6] =  (TH1F*)input->Get("dEtaCase3pp1");
    hist_pp[7] =  (TH1F*)input->Get("dEtaCase3pp2");
    hist_pp[8] =  (TH1F*)input->Get("dEtaCase3pp3");
    hist_pp[9] =  (TH1F*)input->Get("dEtaCase4pp1");
    hist_pp[10] = (TH1F*)input->Get("dEtaCase4pp2");
    hist_pp[11] = (TH1F*)input->Get("dEtaCase4pp3");
   
    TH1F *hist_PV[5];
    hist_PV[0] = (TH1F*)input->Get("dEta1STPV1");
    hist_PV[1] = (TH1F*)input->Get("dEta2NDPV2");
    hist_PV[2] = (TH1F*)input->Get("dEta3RDPV1");
    hist_PV[3] = (TH1F*)input->Get("dEta4THPV2");
    hist_PV[4] = (TH1F*)input->Get("dEta5THPV1");

    cout << "Delta eta(pixel, pixel)" << endl;
    // Loop for dEta(pixel, pixel)
    for(Int_t i = 0; i < 12; i++)
    {   
        cout.width(4);
        cout << i+1 << "th histogram" << endl;

        Float_t denomi;
        if( i == 5 || i == 6 ) denomi = hist_pp[i]->Integral(450,550); // fabs(0.001)
        if( i == 2 || i == 9 ) denomi = hist_pp[i]->Integral(400,600); // fabs(0.002)
        if( i == 0 || i == 3 || i == 8 || i == 11 ) denomi = hist_pp[i]->Integral(375,625); // fabs(0.0025)
        if( i == 4 || i == 7 ) denomi = hist_pp[i]->Integral(325,675); // fabs(0.0035)
        if( i == 1 || i == 10 ) denomi = hist_pp[i]->Integral(300,700); // fabs(0.004)
        
        Int_t maxbin = hist_pp[i]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_pp[i]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout.width(8);
                cout << "Maxbin: " << hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_pp[i]->GetBinCenter(maxbin+k) - hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_pp[i]->GetBinCenter(maxbin) - hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
   
    cout << endl << endl;
    cout << "Delta eta(PV, pixel)" << endl;
    // Loop for dEta(PV, pixel)
    for(Int_t j = 0; j < 5; j++)
    {
        cout.width(4);
        cout << j+1 << "th histogram" << endl;

        Float_t denomi;
        if( j == 4 ) denomi = hist_PV[j]->Integral(450,550); // fabs(0.001)
        if( j == 0 || j == 2 ) denomi = hist_PV[j]->Integral(440,560); // fabs(0.0012)
        if( j == 3 ) denomi = hist_PV[j]->Integral(400,600); // fabs(0.002)
        if( j == 1 ) denomi = hist_PV[j]->Integral(375,625); // fabs(0.0025)

        Int_t maxbin = hist_PV[j]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_PV[j]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "Maxbin: " << hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_PV[j]->GetBinCenter(maxbin+k) - hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_PV[j]->GetBinCenter(maxbin) - hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
    cout << endl;
}

void fit_Region2()
{
    TFile *input = new TFile("dAngleEtaRegion2.root","OPEN");

    TH1F *hist_pp[12];
    hist_pp[0] =  (TH1F*)input->Get("dEtaCase1pp1");
    hist_pp[1] =  (TH1F*)input->Get("dEtaCase1pp2");
    hist_pp[2] =  (TH1F*)input->Get("dEtaCase1pp3");
    hist_pp[3] =  (TH1F*)input->Get("dEtaCase2pp1");
    hist_pp[4] =  (TH1F*)input->Get("dEtaCase2pp2");
    hist_pp[5] =  (TH1F*)input->Get("dEtaCase2pp3");
    hist_pp[6] =  (TH1F*)input->Get("dEtaCase3pp1");
    hist_pp[7] =  (TH1F*)input->Get("dEtaCase3pp2");
    hist_pp[8] =  (TH1F*)input->Get("dEtaCase3pp3");
    hist_pp[9] =  (TH1F*)input->Get("dEtaCase4pp1");
    hist_pp[10] = (TH1F*)input->Get("dEtaCase4pp2");
    hist_pp[11] = (TH1F*)input->Get("dEtaCase4pp3");
   
    TH1F *hist_PV[5];
    hist_PV[0] = (TH1F*)input->Get("dEta1STPV1");
    hist_PV[1] = (TH1F*)input->Get("dEta2NDPV2");
    hist_PV[2] = (TH1F*)input->Get("dEta3RDPV1");
    hist_PV[3] = (TH1F*)input->Get("dEta4THPV2");
    hist_PV[4] = (TH1F*)input->Get("dEta5THPV1");

    cout << "Delta eta(pixel, pixel)" << endl;
    // Loop for dEta(pixel, pixel)
    for(Int_t i = 0; i < 12; i++)
    {   
        cout.width(4);
        cout << i+1 << "th histogram" << endl;

        Float_t denomi;
        if( i == 6 ) denomi = hist_pp[i]->Integral(476,525); // fabs(0.0005)
        if( i == 5 || i == 9 ) denomi = hist_pp[i]->Integral(461,540); // fabs(0.0008)
        if( i == 0 || i == 2 || i == 3 || i == 8 || i == 11 ) denomi = hist_pp[i]->Integral(451,550); // fabs(0.001)
        if( i == 4 || i == 7 ) denomi = hist_pp[i]->Integral(426,575); // fabs(0.0015)
        if( i == 1 || i == 10 ) denomi = hist_pp[i]->Integral(401,600); // fabs(0.002)
        
        Int_t maxbin = hist_pp[i]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_pp[i]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout.width(8);
                cout << "Maxbin: " << hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_pp[i]->GetBinCenter(maxbin+k) - hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_pp[i]->GetBinCenter(maxbin) - hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
   
    cout << endl << endl;
    cout << "Delta eta(PV, pixel)" << endl;
    // Loop for dEta(PV, pixel)
    for(Int_t j = 0; j < 5; j++)
    {
        cout.width(4);
        cout << j+1 << "th histogram" << endl;

        Float_t denomi;
        if( j == 4 ) denomi = hist_PV[j]->Integral(486,515); // fabs(0.0003)
        if( j == 2 ) denomi = hist_PV[j]->Integral(476,525); // fabs(0.0005)
        if( j == 0 ) denomi = hist_PV[j]->Integral(471,530); // fabs(0.0006)
        if( j == 3 ) denomi = hist_PV[j]->Integral(450,550); // fabs(0.001)
        if( j == 1 ) denomi = hist_PV[j]->Integral(441,560); // fabs(0.0012)

        Int_t maxbin = hist_PV[j]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_PV[j]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "Maxbin: " << hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_PV[j]->GetBinCenter(maxbin+k) - hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_PV[j]->GetBinCenter(maxbin) - hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
    cout << endl;
}

void fit_Region3()
{
    TFile *input = new TFile("dAngleEtaRegion3.root","OPEN");

    TH1F *hist_pp[12];
    hist_pp[0] =  (TH1F*)input->Get("dEtaCase1pp1");
    hist_pp[1] =  (TH1F*)input->Get("dEtaCase1pp2");
    hist_pp[2] =  (TH1F*)input->Get("dEtaCase1pp3");
    hist_pp[3] =  (TH1F*)input->Get("dEtaCase2pp1");
    hist_pp[4] =  (TH1F*)input->Get("dEtaCase2pp2");
    hist_pp[5] =  (TH1F*)input->Get("dEtaCase2pp3");
    hist_pp[6] =  (TH1F*)input->Get("dEtaCase3pp1");
    hist_pp[7] =  (TH1F*)input->Get("dEtaCase3pp2");
    hist_pp[8] =  (TH1F*)input->Get("dEtaCase3pp3");
    hist_pp[9] =  (TH1F*)input->Get("dEtaCase4pp1");
    hist_pp[10] = (TH1F*)input->Get("dEtaCase4pp2");
    hist_pp[11] = (TH1F*)input->Get("dEtaCase4pp3");
   
    TH1F *hist_PV[5];
    hist_PV[0] = (TH1F*)input->Get("dEta1STPV1");
    hist_PV[1] = (TH1F*)input->Get("dEta2NDPV2");
    hist_PV[2] = (TH1F*)input->Get("dEta3RDPV1");
    hist_PV[3] = (TH1F*)input->Get("dEta4THPV2");
    hist_PV[4] = (TH1F*)input->Get("dEta5THPV1");

    cout << "Delta eta(pixel, pixel)" << endl;
    // Loop for dEta(pixel, pixel)
    for(Int_t i = 0; i < 12; i++)
    {   
        cout.width(4);
        cout << i+1 << "th histogram" << endl;

        Float_t denomi;
        //if( i ==  ) denomi = hist_pp[i]->Integral(475,525); // fabs(0.0005)
        if( i == 5 ) denomi = hist_pp[i]->Integral(466,535); // fabs(0.0007)
        if( i == 0 || i == 3 ) denomi = hist_pp[i]->Integral(451,550); // fabs(0.001)
        if( i == 2 || i == 6 ) denomi = hist_pp[i]->Integral(441,560); // fabs(0.0012)
        if( i == 4 || i == 10 ) denomi = hist_pp[i]->Integral(426,575); // fabs(0.0015)
        if( i == 1 || i == 8 || i == 11 ) denomi = hist_pp[i]->Integral(401,600); // fabs(0.002)
        if( i == 7 || i == 9 ) denomi = hist_pp[i]->Integral(351,650); // fabs(0.003)
        
        Int_t maxbin = hist_pp[i]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_pp[i]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout.width(8);
                cout << "Maxbin: " << hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_pp[i]->GetBinCenter(maxbin+k) - hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_pp[i]->GetBinCenter(maxbin) - hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
   
    cout << endl << endl;
    cout << "Delta eta(PV, pixel)" << endl;
    // Loop for dEta(PV, pixel)
    for(Int_t j = 0; j < 5; j++)
    {
        cout.width(4);
        cout << j+1 << "th histogram" << endl;

        Float_t denomi;
        if( j == 0 || j == 2 ) denomi = hist_PV[j]->Integral(476,525); // fabs(0.0005)
        if( j == 4 ) denomi = hist_PV[j]->Integral(471,530); // fabs(0.0006)
        if( j == 1 ) denomi = hist_PV[j]->Integral(461,540); // fabs(0.0008)
        if( j == 3 ) denomi = hist_PV[j]->Integral(450,550); // fabs(0.001)

        Int_t maxbin = hist_PV[j]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_PV[j]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "Maxbin: " << hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_PV[j]->GetBinCenter(maxbin+k) - hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_PV[j]->GetBinCenter(maxbin) - hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
    cout << endl;
}

void fit_Region4()
{
    TFile *input = new TFile("dAngleEtaRegion4.root","OPEN");

    TH1F *hist_pp[12];
    hist_pp[0] =  (TH1F*)input->Get("dEtaCase1pp1");
    hist_pp[1] =  (TH1F*)input->Get("dEtaCase1pp2");
    hist_pp[2] =  (TH1F*)input->Get("dEtaCase1pp3");
    hist_pp[3] =  (TH1F*)input->Get("dEtaCase2pp1");
    hist_pp[4] =  (TH1F*)input->Get("dEtaCase2pp2");
    hist_pp[5] =  (TH1F*)input->Get("dEtaCase2pp3");
    hist_pp[6] =  (TH1F*)input->Get("dEtaCase3pp1");
    hist_pp[7] =  (TH1F*)input->Get("dEtaCase3pp2");
    hist_pp[8] =  (TH1F*)input->Get("dEtaCase3pp3");
    hist_pp[9] =  (TH1F*)input->Get("dEtaCase4pp1");
    hist_pp[10] = (TH1F*)input->Get("dEtaCase4pp2");
    hist_pp[11] = (TH1F*)input->Get("dEtaCase4pp3");
   
    TH1F *hist_PV[5];
    hist_PV[0] = (TH1F*)input->Get("dEta1STPV1");
    hist_PV[1] = (TH1F*)input->Get("dEta2NDPV2");
    hist_PV[2] = (TH1F*)input->Get("dEta3RDPV1");
    hist_PV[3] = (TH1F*)input->Get("dEta4THPV2");
    hist_PV[4] = (TH1F*)input->Get("dEta5THPV1");

    cout << "Delta eta(pixel, pixel)" << endl;
    // Loop for dEta(pixel, pixel)
    for(Int_t i = 0; i < 12; i++)
    {   
        cout.width(4);
        cout << i+1 << "th histogram" << endl;

        Float_t denomi;
        if( i == 5 || i == 6 ) denomi = hist_pp[i]->Integral(401,600); // fabs(0.002)
        if( i == 0 || i == 3 ) denomi = hist_pp[i]->Integral(351,650); // fabs(0.003)
        if( i == 8 || i == 11 ) denomi = hist_pp[i]->Integral(326,675); // fabs(0.0035)
        if( i == 2 || i == 4 || i == 7 || i == 9 ) denomi = hist_pp[i]->Integral(251,750); // fabs(0.005)
        if( i == 1 || i == 10 ) denomi = hist_pp[i]->Integral(101,900); // fabs(0.008)
        
        Int_t maxbin = hist_pp[i]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_pp[i]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout.width(8);
                cout << "Maxbin: " << hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_pp[i]->GetBinCenter(maxbin+k) - hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_pp[i]->GetBinCenter(maxbin) - hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
   
    cout << endl << endl;
    cout << "Delta eta(PV, pixel)" << endl;
    // Loop for dEta(PV, pixel)
    for(Int_t j = 0; j < 5; j++)
    {
        cout.width(4);
        cout << j+1 << "th histogram" << endl;

        Float_t denomi;
        if( j == 1 || j == 3 || j == 4 ) denomi = hist_PV[j]->Integral(451,550); // fabs(0.001)
        if( j == 0 || j == 2 ) denomi = hist_PV[j]->Integral(441,560); // fabs(0.0012)

        Int_t maxbin = hist_PV[j]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_PV[j]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "Maxbin: " << hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_PV[j]->GetBinCenter(maxbin+k) - hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_PV[j]->GetBinCenter(maxbin) - hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
    cout << endl;
}

void fit_Region5()
{
    TFile *input = new TFile("dAngleEtaRegion5.root","OPEN");

    TH1F *hist_pp[12];
    hist_pp[0] =  (TH1F*)input->Get("dEtaCase1pp1");
    hist_pp[1] =  (TH1F*)input->Get("dEtaCase1pp2");
    hist_pp[2] =  (TH1F*)input->Get("dEtaCase1pp3");
    hist_pp[3] =  (TH1F*)input->Get("dEtaCase2pp1");
    hist_pp[4] =  (TH1F*)input->Get("dEtaCase2pp2");
    hist_pp[5] =  (TH1F*)input->Get("dEtaCase2pp3");
    hist_pp[6] =  (TH1F*)input->Get("dEtaCase3pp1");
    hist_pp[7] =  (TH1F*)input->Get("dEtaCase3pp2");
    hist_pp[8] =  (TH1F*)input->Get("dEtaCase3pp3");
    hist_pp[9] =  (TH1F*)input->Get("dEtaCase4pp1");
    hist_pp[10] = (TH1F*)input->Get("dEtaCase4pp2");
    hist_pp[11] = (TH1F*)input->Get("dEtaCase4pp3");
   
    TH1F *hist_PV[5];
    hist_PV[0] = (TH1F*)input->Get("dEta1STPV1");
    hist_PV[1] = (TH1F*)input->Get("dEta2NDPV2");
    hist_PV[2] = (TH1F*)input->Get("dEta3RDPV1");
    hist_PV[3] = (TH1F*)input->Get("dEta4THPV2");
    hist_PV[4] = (TH1F*)input->Get("dEta5THPV1");

    cout << "Delta eta(pixel, pixel)" << endl;
    // Loop for dEta(pixel, pixel)
    for(Int_t i = 0; i < 12; i++)
    {   
        cout.width(4);
        cout << i+1 << "th histogram" << endl;
        if( i == 1 ) cout << "Pass this histogram since cut < 0.01" << endl << endl;

        Float_t denomi;
        if( i == 5 ) denomi = hist_pp[i]->Integral(351,650); // fabs(0.003)
        if( i == 6 ) denomi = hist_pp[i]->Integral(326,675); // fabs(0.0035)
        if( i == 8 || i == 11 ) denomi = hist_pp[i]->Integral(251,650); // fabs(0.005)
        if( i == 2 || i == 9 ) denomi = hist_pp[i]->Integral(201,800); // fabs(0.006)
        if( i == 0 ) denomi = hist_pp[i]->Integral(151,850); // fabs(0.007)
        if( i == 3 || i == 7 ) denomi = hist_pp[i]->Integral(101,900); // fabs(0.008)
        if( i == 4 || i == 10 ) denomi = hist_pp[i]->Integral(); // fabs(0.01)
        if( i == 1 ) continue; // D23 - D34 < 0.01
        
        Int_t maxbin = hist_pp[i]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_pp[i]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout.width(8);
                cout << "Maxbin: " << hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_pp[i]->GetBinCenter(maxbin+k) - hist_pp[i]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_pp[i]->GetBinCenter(maxbin) - hist_pp[i]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
   
    cout << endl << endl;
    cout << "Delta eta(PV, pixel)" << endl;
    // Loop for dEta(PV, pixel)
    for(Int_t j = 0; j < 5; j++)
    {
        cout.width(4);
        cout << j+1 << "th histogram" << endl;

        Float_t denomi;
        if( j == 4 ) denomi = hist_PV[j]->Integral(436,565); // fabs(0.0013)
        if( j == 0 || j == 2 ) denomi = hist_PV[j]->Integral(416,585); // fabs(0.0017)
        if( j == 1 || j == 3 ) denomi = hist_PV[j]->Integral(401,600); // fabs(0.002)

        Int_t maxbin = hist_PV[j]->GetMaximumBin();
        for(Int_t k = 0; k < 1000; k++)
        {
            Float_t nomi = hist_PV[j]->Integral(maxbin-k, maxbin+k);
            Float_t ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "Maxbin: " << hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "High 3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin+k) << endl;
                cout << "Low  3sigma bin: " << hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "Distance1: " << hist_PV[j]->GetBinCenter(maxbin+k) - hist_PV[j]->GetBinCenter(maxbin) << endl;
                cout << "Distance2: " << hist_PV[j]->GetBinCenter(maxbin) - hist_PV[j]->GetBinCenter(maxbin-k) << endl;
                cout << "================================" << endl;
                break;
            }
        }
    }
    cout << endl;
}

