//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 30 17:58:18 2018 by ROOT version 6.12/06
// from TTree t/t
// found on file: Tree_Minbias_PU140.root
//////////////////////////////////////////////////////////

#ifndef bkg_h
#define bkg_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class bkg {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           totalEvent;
   Float_t         totalEgN;
   Int_t           ntnEg2;
   Int_t           ntnEg3;
   Int_t           event_with_selection_passed;
   Int_t           event_with_trigger_passed;
   Int_t           L1TkEleN;
   Int_t           L1TkEleIsoN;
   vector<float>   *ntEgEt;
   vector<float>   *ntL1TkEleEt;
   vector<float>   *ntL1TkEleIsoEt;
   vector<float>   *ntEgEta;
   vector<bool>    *ntclusterIsEG;
   vector<float>   *ntL1TkEleEta;
   vector<float>   *ntL1TkEleIsoEta;
   vector<float>   *ntL1TkElePhi;
   vector<float>   *ntL1TkEleIsoPhi;
   vector<bool>    *ntCl_match;
   vector<bool>    *withoutEM_match;
   vector<bool>    *withEM_match;
   vector<float>   *ntEgEt_iso;
   vector<float>   *ntEgEta_iso;
   vector<float>   *ntEgPhi_iso;
   vector<bool>    *ntCl_iso_match;
   vector<bool>    *only_iso_match;
   vector<bool>    *ntCl_match_wo4thPix;
   vector<bool>    *ntCl_match_wo3thPix;
   vector<bool>    *ntCl_match_wo2thPix;
   vector<bool>    *ntCl_match_wo1thPix;
   vector<int>     *Npass_woEM_wo4thPix;
   vector<int>     *Npass_woEM_wo3thPix;
   vector<int>     *Npass_woEM_wo2thPix;
   vector<int>     *Npass_woEM_wo1thPix;
   vector<int>     *Npass_wEM_wo4thPix;
   vector<int>     *Npass_wEM_wo3thPix;
   vector<int>     *Npass_wEM_wo2thPix;
   vector<int>     *Npass_wEM_wo1thPix;
   vector<float>   *deta_L12EM_wo4thPix;
   vector<float>   *deta_L12EM_wo3thPix;
   vector<float>   *deta_L12EM_wo2thPix;
   vector<float>   *deta_L12EM_wo1thPix;
   vector<float>   *dphi_L12EM_wo4thPix;
   vector<float>   *dphi_L12EM_wo3thPix;
   vector<float>   *dphi_L12EM_wo2thPix;
   vector<float>   *dphi_L12EM_wo1thPix;
   vector<int>     *pass_Ele;
   vector<int>     *pass_Pos;
   vector<int>     *pass_ElePos;
   vector<int>     *ntfirstPix;
   vector<int>     *ntsecondPix;
   vector<int>     *ntthirdPix;
   vector<int>     *ntfourthPix;
   vector<float>   *recon_pT;

   // List of branches
   TBranch        *b_count_Entry;   //!
   TBranch        *b_EgN;   //!
   TBranch        *b_ntnEg2;   //!
   TBranch        *b_ntnEg3;   //!
   TBranch        *b_event_denominator;   //!
   TBranch        *b_event_nominator;   //!
   TBranch        *b_L1TkEleN;   //!
   TBranch        *b_L1TkEleIsoN;   //!
   TBranch        *b_ntEgEt;   //!
   TBranch        *b_ntL1TkEleEt;   //!
   TBranch        *b_ntL1TkEleIsoEt;   //!
   TBranch        *b_ntEgEta;   //!
   TBranch        *b_ntclusterIsEG;   //!
   TBranch        *b_ntL1TkEleEta;   //!
   TBranch        *b_ntL1TkEleIsoEta;   //!
   TBranch        *b_ntL1TkElePhi;   //!
   TBranch        *b_ntL1TkEleIsoPhi;   //!
   TBranch        *b_ntCl_match;   //!
   TBranch        *b_withoutEM_match;   //!
   TBranch        *b_withEM_match;   //!
   TBranch        *b_ntEgEt_iso;   //!
   TBranch        *b_ntEgEta_iso;   //!
   TBranch        *b_ntEgPhi_iso;   //!
   TBranch        *b_ntCl_iso_match;   //!
   TBranch        *b_only_iso_match;   //!
   TBranch        *b_ntCl_match_wo4thPix;   //!
   TBranch        *b_ntCl_match_wo3thPix;   //!
   TBranch        *b_ntCl_match_wo2thPix;   //!
   TBranch        *b_ntCl_match_wo1thPix;   //!
   TBranch        *b_Npass_woEM_wo4thPix;   //!
   TBranch        *b_Npass_woEM_wo3thPix;   //!
   TBranch        *b_Npass_woEM_wo2thPix;   //!
   TBranch        *b_Npass_woEM_wo1thPix;   //!
   TBranch        *b_Npass_wEM_wo4thPix;   //!
   TBranch        *b_Npass_wEM_wo3thPix;   //!
   TBranch        *b_Npass_wEM_wo2thPix;   //!
   TBranch        *b_Npass_wEM_wo1thPix;   //!
   TBranch        *b_deta_L12EM_wo4thPix;   //!
   TBranch        *b_deta_L12EM_wo3thPix;   //!
   TBranch        *b_deta_L12EM_wo2thPix;   //!
   TBranch        *b_deta_L12EM_wo1thPix;   //!
   TBranch        *b_dphi_L12EM_wo4thPix;   //!
   TBranch        *b_dphi_L12EM_wo3thPix;   //!
   TBranch        *b_dphi_L12EM_wo2thPix;   //!
   TBranch        *b_dphi_L12EM_wo1thPix;   //!
   TBranch        *b_pass_Ele;   //!
   TBranch        *b_pass_Pos;   //!
   TBranch        *b_pass_ElePos;   //!
   TBranch        *b_ntfirstPix;   //!
   TBranch        *b_ntsecondPix;   //!
   TBranch        *b_ntthirdPix;   //!
   TBranch        *b_ntfourthPix;   //!
   TBranch        *b_recon_pT;   //!

   bkg(TTree *tree=0);
   virtual ~bkg();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef bkg_cxx
bkg::bkg(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree_Minbias_PU140.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Tree_Minbias_PU140.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

bkg::~bkg()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t bkg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t bkg::LoadTree(Long64_t entry)
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

void bkg::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ntEgEt = 0;
   ntL1TkEleEt = 0;
   ntL1TkEleIsoEt = 0;
   ntEgEta = 0;
   ntclusterIsEG = 0;
   ntL1TkEleEta = 0;
   ntL1TkEleIsoEta = 0;
   ntL1TkElePhi = 0;
   ntL1TkEleIsoPhi = 0;
   ntCl_match = 0;
   withoutEM_match = 0;
   withEM_match = 0;
   ntEgEt_iso = 0;
   ntEgEta_iso = 0;
   ntEgPhi_iso = 0;
   ntCl_iso_match = 0;
   only_iso_match = 0;
   ntCl_match_wo4thPix = 0;
   ntCl_match_wo3thPix = 0;
   ntCl_match_wo2thPix = 0;
   ntCl_match_wo1thPix = 0;
   Npass_woEM_wo4thPix = 0;
   Npass_woEM_wo3thPix = 0;
   Npass_woEM_wo2thPix = 0;
   Npass_woEM_wo1thPix = 0;
   Npass_wEM_wo4thPix = 0;
   Npass_wEM_wo3thPix = 0;
   Npass_wEM_wo2thPix = 0;
   Npass_wEM_wo1thPix = 0;
   deta_L12EM_wo4thPix = 0;
   deta_L12EM_wo3thPix = 0;
   deta_L12EM_wo2thPix = 0;
   deta_L12EM_wo1thPix = 0;
   dphi_L12EM_wo4thPix = 0;
   dphi_L12EM_wo3thPix = 0;
   dphi_L12EM_wo2thPix = 0;
   dphi_L12EM_wo1thPix = 0;
   pass_Ele = 0;
   pass_Pos = 0;
   pass_ElePos = 0;
   ntfirstPix = 0;
   ntsecondPix = 0;
   ntthirdPix = 0;
   ntfourthPix = 0;
   recon_pT = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("totalEvent", &totalEvent, &b_count_Entry);
   fChain->SetBranchAddress("totalEgN", &totalEgN, &b_EgN);
   fChain->SetBranchAddress("ntnEg2", &ntnEg2, &b_ntnEg2);
   fChain->SetBranchAddress("ntnEg3", &ntnEg3, &b_ntnEg3);
   fChain->SetBranchAddress("event_with_selection_passed", &event_with_selection_passed, &b_event_denominator);
   fChain->SetBranchAddress("event_with_trigger_passed", &event_with_trigger_passed, &b_event_nominator);
   fChain->SetBranchAddress("L1TkEleN", &L1TkEleN, &b_L1TkEleN);
   fChain->SetBranchAddress("L1TkEleIsoN", &L1TkEleIsoN, &b_L1TkEleIsoN);
   fChain->SetBranchAddress("ntEgEt", &ntEgEt, &b_ntEgEt);
   fChain->SetBranchAddress("ntL1TkEleEt", &ntL1TkEleEt, &b_ntL1TkEleEt);
   fChain->SetBranchAddress("ntL1TkEleIsoEt", &ntL1TkEleIsoEt, &b_ntL1TkEleIsoEt);
   fChain->SetBranchAddress("ntEgEta", &ntEgEta, &b_ntEgEta);
   fChain->SetBranchAddress("ntclusterIsEG", &ntclusterIsEG, &b_ntclusterIsEG);
   fChain->SetBranchAddress("ntL1TkEleEta", &ntL1TkEleEta, &b_ntL1TkEleEta);
   fChain->SetBranchAddress("ntL1TkEleIsoEta", &ntL1TkEleIsoEta, &b_ntL1TkEleIsoEta);
   fChain->SetBranchAddress("ntL1TkElePhi", &ntL1TkElePhi, &b_ntL1TkElePhi);
   fChain->SetBranchAddress("ntL1TkEleIsoPhi", &ntL1TkEleIsoPhi, &b_ntL1TkEleIsoPhi);
   fChain->SetBranchAddress("ntCl_match", &ntCl_match, &b_ntCl_match);
   fChain->SetBranchAddress("withoutEM_match", &withoutEM_match, &b_withoutEM_match);
   fChain->SetBranchAddress("withEM_match", &withEM_match, &b_withEM_match);
   fChain->SetBranchAddress("ntEgEt_iso", &ntEgEt_iso, &b_ntEgEt_iso);
   fChain->SetBranchAddress("ntEgEta_iso", &ntEgEta_iso, &b_ntEgEta_iso);
   fChain->SetBranchAddress("ntEgPhi_iso", &ntEgPhi_iso, &b_ntEgPhi_iso);
   fChain->SetBranchAddress("ntCl_iso_match", &ntCl_iso_match, &b_ntCl_iso_match);
   fChain->SetBranchAddress("only_iso_match", &only_iso_match, &b_only_iso_match);
   fChain->SetBranchAddress("ntCl_match_wo4thPix", &ntCl_match_wo4thPix, &b_ntCl_match_wo4thPix);
   fChain->SetBranchAddress("ntCl_match_wo3thPix", &ntCl_match_wo3thPix, &b_ntCl_match_wo3thPix);
   fChain->SetBranchAddress("ntCl_match_wo2thPix", &ntCl_match_wo2thPix, &b_ntCl_match_wo2thPix);
   fChain->SetBranchAddress("ntCl_match_wo1thPix", &ntCl_match_wo1thPix, &b_ntCl_match_wo1thPix);
   fChain->SetBranchAddress("Npass_woEM_wo4thPix", &Npass_woEM_wo4thPix, &b_Npass_woEM_wo4thPix);
   fChain->SetBranchAddress("Npass_woEM_wo3thPix", &Npass_woEM_wo3thPix, &b_Npass_woEM_wo3thPix);
   fChain->SetBranchAddress("Npass_woEM_wo2thPix", &Npass_woEM_wo2thPix, &b_Npass_woEM_wo2thPix);
   fChain->SetBranchAddress("Npass_woEM_wo1thPix", &Npass_woEM_wo1thPix, &b_Npass_woEM_wo1thPix);
   fChain->SetBranchAddress("Npass_wEM_wo4thPix", &Npass_wEM_wo4thPix, &b_Npass_wEM_wo4thPix);
   fChain->SetBranchAddress("Npass_wEM_wo3thPix", &Npass_wEM_wo3thPix, &b_Npass_wEM_wo3thPix);
   fChain->SetBranchAddress("Npass_wEM_wo2thPix", &Npass_wEM_wo2thPix, &b_Npass_wEM_wo2thPix);
   fChain->SetBranchAddress("Npass_wEM_wo1thPix", &Npass_wEM_wo1thPix, &b_Npass_wEM_wo1thPix);
   fChain->SetBranchAddress("deta_L12EM_wo4thPix", &deta_L12EM_wo4thPix, &b_deta_L12EM_wo4thPix);
   fChain->SetBranchAddress("deta_L12EM_wo3thPix", &deta_L12EM_wo3thPix, &b_deta_L12EM_wo3thPix);
   fChain->SetBranchAddress("deta_L12EM_wo2thPix", &deta_L12EM_wo2thPix, &b_deta_L12EM_wo2thPix);
   fChain->SetBranchAddress("deta_L12EM_wo1thPix", &deta_L12EM_wo1thPix, &b_deta_L12EM_wo1thPix);
   fChain->SetBranchAddress("dphi_L12EM_wo4thPix", &dphi_L12EM_wo4thPix, &b_dphi_L12EM_wo4thPix);
   fChain->SetBranchAddress("dphi_L12EM_wo3thPix", &dphi_L12EM_wo3thPix, &b_dphi_L12EM_wo3thPix);
   fChain->SetBranchAddress("dphi_L12EM_wo2thPix", &dphi_L12EM_wo2thPix, &b_dphi_L12EM_wo2thPix);
   fChain->SetBranchAddress("dphi_L12EM_wo1thPix", &dphi_L12EM_wo1thPix, &b_dphi_L12EM_wo1thPix);
   fChain->SetBranchAddress("pass_Ele", &pass_Ele, &b_pass_Ele);
   fChain->SetBranchAddress("pass_Pos", &pass_Pos, &b_pass_Pos);
   fChain->SetBranchAddress("pass_ElePos", &pass_ElePos, &b_pass_ElePos);
   fChain->SetBranchAddress("ntfirstPix", &ntfirstPix, &b_ntfirstPix);
   fChain->SetBranchAddress("ntsecondPix", &ntsecondPix, &b_ntsecondPix);
   fChain->SetBranchAddress("ntthirdPix", &ntthirdPix, &b_ntthirdPix);
   fChain->SetBranchAddress("ntfourthPix", &ntfourthPix, &b_ntfourthPix);
   fChain->SetBranchAddress("recon_pT", &recon_pT, &b_recon_pT);
   Notify();
}

Bool_t bkg::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void bkg::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t bkg::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef bkg_cxx
