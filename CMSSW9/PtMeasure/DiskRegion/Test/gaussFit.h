//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep 16 16:38:32 2018 by ROOT version 6.10/05
// from TTree t/t
// found on file: data.root
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
   vector<float>   *bsl1_l1d1_dphi;
   vector<float>   *bsl1_l1d2_dphi;
   vector<float>   *bsl2_l2d1_dphi;
   vector<float>   *bsl2_l2d2_dphi;
   vector<float>   *bsl3_l3d1_dphi;
   vector<float>   *bsd1_d1d2_dphi;
   vector<float>   *bsd1_d1d3_dphi;
   vector<float>   *bsd1_d1d4_dphi;
   vector<float>   *bsd2_d2d4_dphi;
   vector<float>   *bsd2_d2d5_dphi;
   vector<float>   *bsd3_d3d4_dphi;
   vector<float>   *bsd3_d3d5_dphi;
   vector<float>   *bsd3_d3d6_dphi;
   vector<float>   *bsd4_d4d6_dphi;
   vector<float>   *genpT_bsl1_l1d1;
   vector<float>   *genpT_bsl1_l1d2;
   vector<float>   *genpT_bsl2_l2d1;
   vector<float>   *genpT_bsl2_l2d2;
   vector<float>   *genpT_bsl3_l3d1;
   vector<float>   *genpT_bsd1_d1d2;
   vector<float>   *genpT_bsd1_d1d3;
   vector<float>   *genpT_bsd1_d1d4;
   vector<float>   *genpT_bsd2_d2d4;
   vector<float>   *genpT_bsd2_d2d5;
   vector<float>   *genpT_bsd3_d3d4;
   vector<float>   *genpT_bsd3_d3d5;
   vector<float>   *genpT_bsd3_d3d6;
   vector<float>   *genpT_bsd4_d4d6;

   // List of branches
   TBranch        *b_bsl1_l1d1_dphi;   //!
   TBranch        *b_bsl1_l1d2_dphi;   //!
   TBranch        *b_bsl2_l2d1_dphi;   //!
   TBranch        *b_bsl2_l2d2_dphi;   //!
   TBranch        *b_bsl3_l3d1_dphi;   //!
   TBranch        *b_bsd1_d1d2_dphi;   //!
   TBranch        *b_bsd1_d1d3_dphi;   //!
   TBranch        *b_bsd1_d1d4_dphi;   //!
   TBranch        *b_bsd2_d2d4_dphi;   //!
   TBranch        *b_bsd2_d2d5_dphi;   //!
   TBranch        *b_bsd3_d3d4_dphi;   //!
   TBranch        *b_bsd3_d3d5_dphi;   //!
   TBranch        *b_bsd3_d3d6_dphi;   //!
   TBranch        *b_bsd4_d4d6_dphi;   //!
   TBranch        *b_genpT_bsl1_l1d1;   //!
   TBranch        *b_genpT_bsl1_l1d2;   //!
   TBranch        *b_genpT_bsl2_l2d1;   //!
   TBranch        *b_genpT_bsl2_l2d2;   //!
   TBranch        *b_genpT_bsl3_l3d1;   //!
   TBranch        *b_genpT_bsd1_d1d2;   //!
   TBranch        *b_genpT_bsd1_d1d3;   //!
   TBranch        *b_genpT_bsd1_d1d4;   //!
   TBranch        *b_genpT_bsd2_d2d4;   //!
   TBranch        *b_genpT_bsd2_d2d5;   //!
   TBranch        *b_genpT_bsd3_d3d4;   //!
   TBranch        *b_genpT_bsd3_d3d5;   //!
   TBranch        *b_genpT_bsd3_d3d6;   //!
   TBranch        *b_genpT_bsd4_d4d6;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data.root");
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
   bsl1_l1d1_dphi = 0;
   bsl1_l1d2_dphi = 0;
   bsl2_l2d1_dphi = 0;
   bsl2_l2d2_dphi = 0;
   bsl3_l3d1_dphi = 0;
   bsd1_d1d2_dphi = 0;
   bsd1_d1d3_dphi = 0;
   bsd1_d1d4_dphi = 0;
   bsd2_d2d4_dphi = 0;
   bsd2_d2d5_dphi = 0;
   bsd3_d3d4_dphi = 0;
   bsd3_d3d5_dphi = 0;
   bsd3_d3d6_dphi = 0;
   bsd4_d4d6_dphi = 0;
   genpT_bsl1_l1d1 = 0;
   genpT_bsl1_l1d2 = 0;
   genpT_bsl2_l2d1 = 0;
   genpT_bsl2_l2d2 = 0;
   genpT_bsl3_l3d1 = 0;
   genpT_bsd1_d1d2 = 0;
   genpT_bsd1_d1d3 = 0;
   genpT_bsd1_d1d4 = 0;
   genpT_bsd2_d2d4 = 0;
   genpT_bsd2_d2d5 = 0;
   genpT_bsd3_d3d4 = 0;
   genpT_bsd3_d3d5 = 0;
   genpT_bsd3_d3d6 = 0;
   genpT_bsd4_d4d6 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bsl1_l1d1_dphi", &bsl1_l1d1_dphi, &b_bsl1_l1d1_dphi);
   fChain->SetBranchAddress("bsl1_l1d2_dphi", &bsl1_l1d2_dphi, &b_bsl1_l1d2_dphi);
   fChain->SetBranchAddress("bsl2_l2d1_dphi", &bsl2_l2d1_dphi, &b_bsl2_l2d1_dphi);
   fChain->SetBranchAddress("bsl2_l2d2_dphi", &bsl2_l2d2_dphi, &b_bsl2_l2d2_dphi);
   fChain->SetBranchAddress("bsl3_l3d1_dphi", &bsl3_l3d1_dphi, &b_bsl3_l3d1_dphi);
   fChain->SetBranchAddress("bsd1_d1d2_dphi", &bsd1_d1d2_dphi, &b_bsd1_d1d2_dphi);
   fChain->SetBranchAddress("bsd1_d1d3_dphi", &bsd1_d1d3_dphi, &b_bsd1_d1d3_dphi);
   fChain->SetBranchAddress("bsd1_d1d4_dphi", &bsd1_d1d4_dphi, &b_bsd1_d1d4_dphi);
   fChain->SetBranchAddress("bsd2_d2d4_dphi", &bsd2_d2d4_dphi, &b_bsd2_d2d4_dphi);
   fChain->SetBranchAddress("bsd2_d2d5_dphi", &bsd2_d2d5_dphi, &b_bsd2_d2d5_dphi);
   fChain->SetBranchAddress("bsd3_d3d4_dphi", &bsd3_d3d4_dphi, &b_bsd3_d3d4_dphi);
   fChain->SetBranchAddress("bsd3_d3d5_dphi", &bsd3_d3d5_dphi, &b_bsd3_d3d5_dphi);
   fChain->SetBranchAddress("bsd3_d3d6_dphi", &bsd3_d3d6_dphi, &b_bsd3_d3d6_dphi);
   fChain->SetBranchAddress("bsd4_d4d6_dphi", &bsd4_d4d6_dphi, &b_bsd4_d4d6_dphi);
   fChain->SetBranchAddress("genpT_bsl1_l1d1", &genpT_bsl1_l1d1, &b_genpT_bsl1_l1d1);
   fChain->SetBranchAddress("genpT_bsl1_l1d2", &genpT_bsl1_l1d2, &b_genpT_bsl1_l1d2);
   fChain->SetBranchAddress("genpT_bsl2_l2d1", &genpT_bsl2_l2d1, &b_genpT_bsl2_l2d1);
   fChain->SetBranchAddress("genpT_bsl2_l2d2", &genpT_bsl2_l2d2, &b_genpT_bsl2_l2d2);
   fChain->SetBranchAddress("genpT_bsl3_l3d1", &genpT_bsl3_l3d1, &b_genpT_bsl3_l3d1);
   fChain->SetBranchAddress("genpT_bsd1_d1d2", &genpT_bsd1_d1d2, &b_genpT_bsd1_d1d2);
   fChain->SetBranchAddress("genpT_bsd1_d1d3", &genpT_bsd1_d1d3, &b_genpT_bsd1_d1d3);
   fChain->SetBranchAddress("genpT_bsd1_d1d4", &genpT_bsd1_d1d4, &b_genpT_bsd1_d1d4);
   fChain->SetBranchAddress("genpT_bsd2_d2d4", &genpT_bsd2_d2d4, &b_genpT_bsd2_d2d4);
   fChain->SetBranchAddress("genpT_bsd2_d2d5", &genpT_bsd2_d2d5, &b_genpT_bsd2_d2d5);
   fChain->SetBranchAddress("genpT_bsd3_d3d4", &genpT_bsd3_d3d4, &b_genpT_bsd3_d3d4);
   fChain->SetBranchAddress("genpT_bsd3_d3d5", &genpT_bsd3_d3d5, &b_genpT_bsd3_d3d5);
   fChain->SetBranchAddress("genpT_bsd3_d3d6", &genpT_bsd3_d3d6, &b_genpT_bsd3_d3d6);
   fChain->SetBranchAddress("genpT_bsd4_d4d6", &genpT_bsd4_d4d6, &b_genpT_bsd4_d4d6);
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
