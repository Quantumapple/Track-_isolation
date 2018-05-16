//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 25 23:38:21 2016 by ROOT version 6.04/08
// from TTree t/t
// found on file: ROI.root
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
   vector<float>   *bpix_l1_dphi_roi;
   vector<float>   *bpix_l2_dphi_roi;
   vector<float>   *bpix_l3_dphi_roi;
   vector<float>   *bpix_l4_dphi_roi;

   // List of branches
   TBranch        *b_abs_pterr;   //!
   TBranch        *b_closest_egEt;   //!
   TBranch        *b_closest_egEta;   //!
   TBranch        *b_closest_egPhi;   //!
   TBranch        *b_closest_eg_dr;   //!
   TBranch        *b_bpix_l1_dphi_roi;   //!
   TBranch        *b_bpix_l2_dphi_roi;   //!
   TBranch        *b_bpix_l3_dphi_roi;   //!
   TBranch        *b_bpix_l4_dphi_roi;   //!

   roofit_input(TTree *tree=0);
   virtual ~roofit_input();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(float low_et, float high_et);
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ROI.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ROI.root");
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
   bpix_l1_dphi_roi = 0;
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
   fChain->SetBranchAddress("bpix_l1_dphi_roi", &bpix_l1_dphi_roi, &b_bpix_l1_dphi_roi);
   fChain->SetBranchAddress("bpix_l2_dphi_roi", &bpix_l2_dphi_roi, &b_bpix_l2_dphi_roi);
   fChain->SetBranchAddress("bpix_l3_dphi_roi", &bpix_l3_dphi_roi, &b_bpix_l3_dphi_roi);
   fChain->SetBranchAddress("bpix_l4_dphi_roi", &bpix_l4_dphi_roi, &b_bpix_l4_dphi_roi);
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
#endif // #ifdef roofit_input_cxx
