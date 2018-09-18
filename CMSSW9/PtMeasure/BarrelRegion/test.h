//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 24 16:29:26 2018 by ROOT version 6.06/01
// from TTree t/t
// found on file: data.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *bsl1_l1l2_dphi;
   vector<float>   *bsl1_l1l3_dphi;
   vector<float>   *bsl1_l1l4_dphi;
   vector<float>   *bsl2_l2l3_dphi;
   vector<float>   *bsl2_l2l4_dphi;
   vector<float>   *bsl3_l3l4_dphi;
   vector<float>   *genpT_case1;
   vector<float>   *genpT_case2;
   vector<float>   *genpT_case3;
   vector<float>   *genpT_case4;
   vector<float>   *genpT_case5;
   vector<float>   *genpT_case6;

   // List of branches
   TBranch        *b_bsl1_l1l2_dphi;   //!
   TBranch        *b_bsl1_l1l3_dphi;   //!
   TBranch        *b_bsl1_l1l4_dphi;   //!
   TBranch        *b_bsl2_l2l3_dphi;   //!
   TBranch        *b_bsl2_l2l4_dphi;   //!
   TBranch        *b_bsl3_l3l4_dphi;   //!
   TBranch        *b_genpT_case1;   //!
   TBranch        *b_genpT_case2;   //!
   TBranch        *b_genpT_case3;   //!
   TBranch        *b_genpT_case4;   //!
   TBranch        *b_genpT_case5;   //!
   TBranch        *b_genpT_case6;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int nth_case);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../data.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
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

void test::Init(TTree *tree)
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
   bsl1_l1l3_dphi = 0;
   bsl1_l1l4_dphi = 0;
   bsl2_l2l3_dphi = 0;
   bsl2_l2l4_dphi = 0;
   bsl3_l3l4_dphi = 0;
   genpT_case1 = 0;
   genpT_case2 = 0;
   genpT_case3 = 0;
   genpT_case4 = 0;
   genpT_case5 = 0;
   genpT_case6 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bsl1_l1l2_dphi", &bsl1_l1l2_dphi, &b_bsl1_l1l2_dphi);
   fChain->SetBranchAddress("bsl1_l1l3_dphi", &bsl1_l1l3_dphi, &b_bsl1_l1l3_dphi);
   fChain->SetBranchAddress("bsl1_l1l4_dphi", &bsl1_l1l4_dphi, &b_bsl1_l1l4_dphi);
   fChain->SetBranchAddress("bsl2_l2l3_dphi", &bsl2_l2l3_dphi, &b_bsl2_l2l3_dphi);
   fChain->SetBranchAddress("bsl2_l2l4_dphi", &bsl2_l2l4_dphi, &b_bsl2_l2l4_dphi);
   fChain->SetBranchAddress("bsl3_l3l4_dphi", &bsl3_l3l4_dphi, &b_bsl3_l3l4_dphi);
   fChain->SetBranchAddress("genpT_case1", &genpT_case1, &b_genpT_case1);
   fChain->SetBranchAddress("genpT_case2", &genpT_case2, &b_genpT_case2);
   fChain->SetBranchAddress("genpT_case3", &genpT_case3, &b_genpT_case3);
   fChain->SetBranchAddress("genpT_case4", &genpT_case4, &b_genpT_case4);
   fChain->SetBranchAddress("genpT_case5", &genpT_case5, &b_genpT_case5);
   fChain->SetBranchAddress("genpT_case6", &genpT_case6, &b_genpT_case6);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
