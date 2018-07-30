//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 22 19:38:38 2016 by ROOT version 6.04/08
// from TTree t/t
// found on file: R1_dPhi.root
//////////////////////////////////////////////////////////

#ifndef roofit_input_h
#define roofit_input_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class roofit_input {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         abs_pterr;
   Float_t         closest_egEt;
   Float_t         closest_egEta;
   Float_t         closest_egPhi;
   Float_t         closest_eg_dr;
   vector<float>   *EM_L1_dPhi;
   vector<float>   *EM_L2_dPhi;
   vector<float>   *EM_L3_dPhi;
   vector<float>   *EM_L4_dPhi;
   vector<float>   *EM_pixel_dPhi;

   vector<float>   *EM_D1_dPhi;
   vector<float>   *EM_D2_dPhi;
   vector<float>   *EM_D3_dPhi;
   vector<float>   *EM_D4_dPhi;

   // List of branches
   TBranch        *b_abs_pterr;   //!
   TBranch        *b_closest_egEt;   //!
   TBranch        *b_closest_egEta;   //!
   TBranch        *b_closest_egPhi;   //!
   TBranch        *b_closest_eg_dr;   //!
   TBranch        *b_EM_L1_dPhi;   //!
   TBranch        *b_EM_L2_dPhi;   //!
   TBranch        *b_EM_L3_dPhi;   //!
   TBranch        *b_EM_L4_dPhi;   //!
   TBranch        *b_EM_pixel_dPhi;   //!

   TBranch        *b_EM_D1_dPhi;   //!
   TBranch        *b_EM_D2_dPhi;   //!
   TBranch        *b_EM_D3_dPhi;   //!
   TBranch        *b_EM_D4_dPhi;   //!

   double getMedian(const std::vector<float> &vec);
   double getMedianErr(const std::vector<float> &vec);

   roofit_input(TTree *tree=0);
   virtual ~roofit_input();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef roofit_input_cxx
roofit_input::roofit_input(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./R1_ROI.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./R1_ROI.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

roofit_input::~roofit_input()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t roofit_input::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t roofit_input::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void roofit_input::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   EM_L1_dPhi = 0;
   EM_L2_dPhi = 0;
   EM_L3_dPhi = 0;
   EM_L4_dPhi = 0;
   EM_pixel_dPhi = 0;

   EM_D1_dPhi = 0;
   EM_D2_dPhi = 0;
   EM_D3_dPhi = 0;
   EM_D4_dPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("abs_pterr", &abs_pterr, &b_abs_pterr);
   fChain->SetBranchAddress("closest_egEt", &closest_egEt, &b_closest_egEt);
   fChain->SetBranchAddress("closest_egEta", &closest_egEta, &b_closest_egEta);
   fChain->SetBranchAddress("closest_egPhi", &closest_egPhi, &b_closest_egPhi);
   fChain->SetBranchAddress("closest_eg_dr", &closest_eg_dr, &b_closest_egEt);

   fChain->SetBranchAddress("EM_L1_dPhi", &EM_L1_dPhi, &b_EM_L1_dPhi);
   fChain->SetBranchAddress("EM_L2_dPhi", &EM_L2_dPhi, &b_EM_L2_dPhi);
   fChain->SetBranchAddress("EM_L3_dPhi", &EM_L3_dPhi, &b_EM_L3_dPhi);
   fChain->SetBranchAddress("EM_L4_dPhi", &EM_L4_dPhi, &b_EM_L4_dPhi);

   fChain->SetBranchAddress("EM_D1_dPhi", &EM_D1_dPhi, &b_EM_D1_dPhi);
   fChain->SetBranchAddress("EM_D2_dPhi", &EM_D2_dPhi, &b_EM_D2_dPhi);
   fChain->SetBranchAddress("EM_D3_dPhi", &EM_D3_dPhi, &b_EM_D3_dPhi);
   fChain->SetBranchAddress("EM_D4_dPhi", &EM_D4_dPhi, &b_EM_D4_dPhi);


   fChain->SetBranchAddress("EM_pixel_dPhi", &EM_pixel_dPhi, &b_EM_pixel_dPhi);
   Notify();
}

Bool_t roofit_input::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void roofit_input::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t roofit_input::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

double roofit_input::getMedian(const std::vector<float> &vec)
{
  // /SHarper/SHNtupliser/src/MathFuncs.cc

  double median=0.;

  //odd number, definate median
  if(vec.size() % 2 !=0) {
    int middleNr = (vec.size())/2; // index number
    median = vec[middleNr];
  }else{ //even number, take median as halfway between the two middle values
    int middleNr = (vec.size())/2;
    median= vec[middleNr];
    median+= vec[middleNr-1];
    median/=2.;
  }
  return median;

}

double roofit_input::getMedianErr(const std::vector<float> &vec)
{
  double err = 0.;

  int middleNr_down = (int)((vec.size()) * 0.45);
  int middleNr_up = (int)((vec.size()) * 0.55);

  err = (vec[middleNr_up] - vec[middleNr_down])/2.;

  return err;

}


#endif // #ifdef roofit_input_cxx
