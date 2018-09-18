//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep 16 16:39:09 2018 by ROOT version 6.10/05
// from TTree t/t
// found on file: data_bsl3_l3d1.root
//////////////////////////////////////////////////////////

#ifndef gaussFit_Special_h
#define gaussFit_Special_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class gaussFit_Special {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *bsl3_l3d1_dphi;
   vector<float>   *genpT_bsl3_l3d1;

   // List of branches
   TBranch        *b_bsl3_l3d1_dphi;   //!
   TBranch        *b_genpT_bsl3_l3d1;   //!

   gaussFit_Special(TTree *tree=0);
   virtual ~gaussFit_Special();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gaussFit_Special_cxx
gaussFit_Special::gaussFit_Special(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_bsl3_l3d1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data_bsl3_l3d1.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

gaussFit_Special::~gaussFit_Special()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gaussFit_Special::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gaussFit_Special::LoadTree(Long64_t entry)
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

void gaussFit_Special::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   bsl3_l3d1_dphi = 0;
   genpT_bsl3_l3d1 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bsl3_l3d1_dphi", &bsl3_l3d1_dphi, &b_bsl3_l3d1_dphi);
   fChain->SetBranchAddress("genpT_bsl3_l3d1", &genpT_bsl3_l3d1, &b_genpT_bsl3_l3d1);
   Notify();
}

Bool_t gaussFit_Special::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gaussFit_Special::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gaussFit_Special::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gaussFit_Special_cxx
