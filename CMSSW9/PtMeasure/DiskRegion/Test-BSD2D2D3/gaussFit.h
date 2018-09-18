//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 19 00:10:39 2018 by ROOT version 6.10/05
// from TTree t/t
// found on file: data_bsd2_d2d3.root
//////////////////////////////////////////////////////////

#ifndef gaussFit_h
#define gaussFit_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class gaussFit {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *bsd2_d2d3_dphi;
   vector<float>   *genpT_bsd2_d2d3;

   // List of branches
   TBranch        *b_bsd2_d2d3_dphi;   //!
   TBranch        *b_genpT_bsd2_d2d3;   //!

   gaussFit(TTree *tree=0);
   virtual ~gaussFit();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gaussFit_cxx
gaussFit::gaussFit(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_bsd2_d2d3.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_bsd2_d2d3.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

gaussFit::~gaussFit()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gaussFit::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gaussFit::LoadTree(Long64_t entry)
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

void gaussFit::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   bsd2_d2d3_dphi = 0;
   genpT_bsd2_d2d3 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bsd2_d2d3_dphi", &bsd2_d2d3_dphi, &b_bsd2_d2d3_dphi);
   fChain->SetBranchAddress("genpT_bsd2_d2d3", &genpT_bsd2_d2d3, &b_genpT_bsd2_d2d3);
   Notify();
}

Bool_t gaussFit::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gaussFit::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gaussFit::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gaussFit_cxx
