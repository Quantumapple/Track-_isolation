#include <stdio.h>
#include <TH1.h>
#include <TFile.h>

void fit_Region2()
{
    TFile *file = new TFile("dZEtaRegion2.root","open");
    TH1F *hist[15];
    hist[0] =  (TH1F*)file->Get("FiCombi_case1");
    hist[1] =  (TH1F*)file->Get("FiCombi_case2");
    hist[2] =  (TH1F*)file->Get("FiCombi_case3");
    hist[3] =  (TH1F*)file->Get("SeCombi_case1");
    hist[4] =  (TH1F*)file->Get("SeCombi_case2");
    hist[5] =  (TH1F*)file->Get("SeCombi_case3");
    hist[6] =  (TH1F*)file->Get("ThCombi_case1");
    hist[7] =  (TH1F*)file->Get("ThCombi_case2");
    hist[8] =  (TH1F*)file->Get("ThCombi_case3");
    hist[9] =  (TH1F*)file->Get("FoCombi_case1");
    hist[10] = (TH1F*)file->Get("FoCombi_case2");
    
    for(Int_t j = 0; j < 11; j++)
    {
        Int_t maxbin = hist[j]->GetMaximumBin();
        Float_t ratio = 0.;
        Float_t denomi = 0.;
        // Loop for find 3sigma boundary 
        for(Int_t k = 1; k < 1000; k++)
        {
            if( j <= 5 ) denomi = hist[j]->Integral(181,220); // 0.02
            else denomi = hist[j]->Integral(171,230); // 0.03

            Float_t nomi = hist[j]->Integral(maxbin-k, maxbin+k);
            ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "3sigma boundary" << endl;
                cout << "  " << j << "th histogram" << endl;
                cout.precision(6);
                cout << "     Maxbin: " << hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     High 3sigma bin: " << hist[j]->GetBinCenter(maxbin+k) << endl;
                cout << "     Low  3sigma bin: " << hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << "     Distance1: " << hist[j]->GetBinCenter(maxbin+k) - hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     Distance2: " << hist[j]->GetBinCenter(maxbin) - hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << endl << endl;
                break;
            }
        }
    }
}

void fit_Region3()
{
    TFile *file = new TFile("dZEtaRegion3.root","open");
    TH1F *hist[15];
    hist[0] =  (TH1F*)file->Get("FiCombi_case1");
    hist[1] =  (TH1F*)file->Get("FiCombi_case2");
    hist[2] =  (TH1F*)file->Get("FiCombi_case3");
    hist[3] =  (TH1F*)file->Get("SeCombi_case1");
    hist[4] =  (TH1F*)file->Get("SeCombi_case2");
    hist[5] =  (TH1F*)file->Get("SeCombi_case3");
    hist[6] =  (TH1F*)file->Get("ThCombi_case1");
    hist[7] =  (TH1F*)file->Get("ThCombi_case2");
    hist[8] =  (TH1F*)file->Get("ThCombi_case3");
    hist[9] =  (TH1F*)file->Get("FoCombi_case1");
    hist[10] = (TH1F*)file->Get("FoCombi_case2");
    
    for(Int_t j = 0; j < 11; j++)
    {
        Int_t maxbin = hist[j]->GetMaximumBin();
        Float_t ratio = 0.;
        Float_t denomi = 0.;
      
        if( j == 0 || j == 1 || j == 6 || j == 7 || j == 8 ) denomi = hist[j]->Integral(181,220); // 0.02
        if( j == 2 || j == 3 || j == 4 || j == 9 ) denomi = hist[j]->Integral(171,230); // 0.03
        if( j == 5 || j == 10 ) denomi = hist[j]->Integral(141,260); // 0.06
        
        // Loop for find 3sigma boundary 
        for(Int_t k = 1; k < 1000; k++)
        {
            Float_t nomi = hist[j]->Integral(maxbin-k, maxbin+k);
            ratio = nomi/denomi;
            if( ratio > 0.997 )
            {   
                cout << "3sigma boundary" << endl;
                cout << "nominator: " << nomi << endl;
                cout << "  " << j << "th histogram" << endl;
                cout.precision(6);
                cout << "     Maxbin: " << hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     High 3sigma bin: " << hist[j]->GetBinCenter(maxbin+k) << endl;
                cout << "     Low  3sigma bin: " << hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << "     Distance1: " << hist[j]->GetBinCenter(maxbin+k) - hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     Distance2: " << hist[j]->GetBinCenter(maxbin) - hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << endl << endl;
                break;
            }
        }
    }
}

void fit_Region4()
{
    TFile *file = new TFile("dZEtaRegion4.root","open");
    TH1F *hist[15];
    hist[0] =  (TH1F*)file->Get("FiCombi_case1");
    hist[1] =  (TH1F*)file->Get("FiCombi_case2");
    hist[2] =  (TH1F*)file->Get("FiCombi_case3");
    hist[3] =  (TH1F*)file->Get("SeCombi_case1");
    hist[4] =  (TH1F*)file->Get("SeCombi_case2");
    hist[5] =  (TH1F*)file->Get("SeCombi_case3");
    hist[6] =  (TH1F*)file->Get("ThCombi_case1");
    hist[7] =  (TH1F*)file->Get("ThCombi_case2");
    hist[8] =  (TH1F*)file->Get("ThCombi_case3");
    hist[9] =  (TH1F*)file->Get("FoCombi_case1");
    hist[10] = (TH1F*)file->Get("FoCombi_case2");
    
    for(Int_t j = 0; j < 11; j++)
    {
        Int_t maxbin = hist[j]->GetMaximumBin();
        Float_t ratio = 0.;
        Float_t denomi = 0.;
       
        if( j == 7 || j == 9 ) denomi = hist[j]->Integral(161,240); // 0.04
        if( j == 0 || j == 1 || j == 6 ) denomi = hist[j]->Integral(151,250); // 0.05
        if( j == 8 ) denomi = hist[j]->Integral(121,280); // 0.08
        if( j == 4 ) denomi = hist[j]->Integral(101,300); // 0.1
        if( j == 2 || j == 3 || j == 5 || j == 10 ) denomi = hist[j]->Integral(51,350); // 0.15
        
        // Loop for find 3sigma boundary 
        for(Int_t k = 1; k < 1000; k++)
        {
            Float_t nomi = hist[j]->Integral(maxbin-k, maxbin+k);
            ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "3sigma boundary" << endl;
                cout << "  " << j << "th histogram" << endl;
                cout.precision(6);
                cout << "     Maxbin: " << hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     High 3sigma bin: " << hist[j]->GetBinCenter(maxbin+k) << endl;
                cout << "     Low  3sigma bin: " << hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << "     Distance1: " << hist[j]->GetBinCenter(maxbin+k) - hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     Distance2: " << hist[j]->GetBinCenter(maxbin) - hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << endl << endl;
                break;
            }
        }
    }
}

void fit_Region5()
{
    TFile *file = new TFile("dZEtaRegion5.root","open");
    TH1F *hist[15];
    hist[0] =  (TH1F*)file->Get("FiCombi_case1");
    hist[1] =  (TH1F*)file->Get("FiCombi_case2");
    hist[2] =  (TH1F*)file->Get("FiCombi_case3");
    hist[3] =  (TH1F*)file->Get("SeCombi_case1");
    hist[4] =  (TH1F*)file->Get("SeCombi_case2");
    hist[5] =  (TH1F*)file->Get("SeCombi_case3");
    hist[6] =  (TH1F*)file->Get("ThCombi_case1");
    hist[7] =  (TH1F*)file->Get("ThCombi_case2");
    hist[8] =  (TH1F*)file->Get("ThCombi_case3");
    hist[9] =  (TH1F*)file->Get("FoCombi_case1");
    hist[10] = (TH1F*)file->Get("FoCombi_case2");
    
    for(Int_t j = 0; j < 11; j++)
    {
        Int_t maxbin = hist[j]->GetMaximumBin();
        Float_t ratio = 0.;
        Float_t denomi = 0.;
        // Loop for find 3sigma boundary 
        for(Int_t k = 1; k < 1000; k++)
        {
            if( j == 0 || j == 2 || j == 3 || j == 5 || j == 6 || j == 10 ) continue; // deltaZ < 0.2 (cm) cut
            if( j == 7 ) denomi = hist[j]->Integral(51,350); // 0.15
            if( j == 1 || j == 4 || j == 8 || j == 9 ) denomi = hist[j]->Integral(); // 0.2

            Float_t nomi = hist[j]->Integral(maxbin-k, maxbin+k);
            ratio = nomi/denomi;
            if( ratio > 0.997 )
            {
                cout << "3sigma boundary" << endl;
                cout << "  " << j << "th histogram" << endl;
                cout.precision(6);
                cout << "     Maxbin: " << hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     High 3sigma bin: " << hist[j]->GetBinCenter(maxbin+k) << endl;
                cout << "     Low  3sigma bin: " << hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << "     Distance1: " << hist[j]->GetBinCenter(maxbin+k) - hist[j]->GetBinCenter(maxbin) << endl;
                cout << "     Distance2: " << hist[j]->GetBinCenter(maxbin) - hist[j]->GetBinCenter(maxbin-k) << endl;
                cout << endl << endl;
                break;
            }
        }
    }
}
