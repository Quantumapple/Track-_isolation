//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 17 11:29:42 2018 by ROOT version 6.06/01
// from TTree t/t
// found on file: data.root
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
   vector<float>   *bsl1_l1l2_dphi;
   vector<float>   *genpT;

   // List of branches
   TBranch        *b_bsl1_l1l2_dphi;   //!
   TBranch        *b_genpT;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data.root");
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
   bsl1_l1l2_dphi = 0;
   genpT = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bsl1_l1l2_dphi", &bsl1_l1l2_dphi, &b_bsl1_l1l2_dphi);
   fChain->SetBranchAddress("genpT", &genpT, &b_genpT);
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
