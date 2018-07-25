//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 10 17:07:06 2018 by ROOT version 6.06/01
// from TTree l1pixtrktree/l1pixtrktree
// found on file: Resultall.root
//////////////////////////////////////////////////////////

#ifndef sw_treeMaker_h
#define sw_treeMaker_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class sw_treeMaker {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nVtx;
   Int_t           nMeanPU;
   Int_t           genPartN;
   vector<float>   *genPartE;
   vector<float>   *genPartPt;
   vector<float>   *genPartEta;
   vector<float>   *genPartPhi;
   vector<int>     *genPartCharge;
   vector<int>     *genPartId;
   vector<float>   *propgenPartE;
   vector<float>   *propgenPartPt;
   vector<float>   *propgenPartEta;
   vector<float>   *propgenPartPhi;
   vector<int>     *propgenPartCharge;
   vector<int>     *propgenPartId;
   vector<float>   *propgenPartx;
   vector<float>   *propgenParty;
   vector<float>   *propgenPartz;
   Int_t           simTrkN;
   vector<float>   *simTrkPt;
   vector<float>   *simTrkEta;
   vector<float>   *simTrkPhi;
   vector<int>     *simTrkId;
   vector<int>     *simTrkType;
   vector<float>   *simTrkVx;
   vector<float>   *simTrkVy;
   vector<float>   *simTrkVz;
   vector<float>   *simVx;
   vector<float>   *simVy;
   vector<float>   *simVz;
   Float_t         lastSimtkpt;
   Float_t         initialSimtkpt;
   Int_t           bremflag;
   vector<float>   *Brempos_radius;
   vector<float>   *Brem_eLoss;
   vector<float>   *Brem_ptLoss;
   vector<float>   *Brempos_x;
   vector<float>   *Brempos_y;
   vector<float>   *Brempos_z;
   vector<int>     *bHitLayer;
   vector<int>     *bHitLadder;
   vector<int>     *bHitModule;
   vector<int>     *fHitDisk;
   vector<int>     *fHitBlade;
   vector<int>     *fHitSide;
   vector<int>     *fHitPanel;
   vector<int>     *fHitModule;
   Int_t           bHitN;
   Int_t           fHitN;
   vector<float>   *fHitGx;
   vector<float>   *fHitGy;
   vector<float>   *fHitGz;
   vector<float>   *fHitPhi;
   vector<float>   *fHitEta;
   vector<float>   *fClSize;
   vector<float>   *fClSizeX;
   vector<float>   *fClSizeY;
   vector<float>   *bHitGx;
   vector<float>   *bHitGy;
   vector<float>   *bHitGz;
   vector<float>   *bHitPhi;
   vector<float>   *bHitEta;
   vector<float>   *bClSize;
   vector<float>   *bClSizeX;
   vector<float>   *bClSizeY;
   Int_t           egCrysN;
   vector<float>   *egCrysE;
   vector<float>   *egCrysEt;
   vector<float>   *egCrysEta;
   vector<float>   *egCrysPhi;
   vector<float>   *egCrysGx;
   vector<float>   *egCrysGy;
   vector<float>   *egCrysGz;
   Int_t           egCrysClusterN;
   vector<float>   *egCrysClusterE;
   vector<float>   *egCrysClusterEt;
   vector<float>   *egCrysClusterEta;
   vector<float>   *egCrysClusterPhi;
   vector<float>   *egCrysClusterGx;
   vector<float>   *egCrysClusterGy;
   vector<float>   *egCrysClusterGz;
   vector<float>   *egCrysClusterPGx;
   vector<float>   *egCrysClusterPGy;
   vector<float>   *egCrysClusterPGz;
   vector<bool>    *isTrackMatched;
   vector<float>   *isoConeNTrack;
   vector<float>   *isoConePtTrack;
   vector<float>   *trackHighestPt;
   vector<float>   *trackHighestPtEta;
   vector<float>   *trackHighestPtPhi;
   vector<float>   *trackHighestPtChi2;
   vector<float>   *trackHighestPtCutChi2;
   vector<float>   *trackHighestPtCutChi2Eta;
   vector<float>   *trackHighestPtCutChi2Phi;
   vector<float>   *trackHighestPtCutChi2Chi2;
   vector<float>   *trackmatchingdR;
   vector<bool>    *hgcal_isTrackMatched;
   vector<float>   *hgcal_isoConeNTrack;
   vector<float>   *hgcal_isoConePtTrack;
   vector<float>   *hgcal_trackHighestPt;
   vector<float>   *hgcal_trackHighestPtEta;
   vector<float>   *hgcal_trackHighestPtPhi;
   vector<float>   *hgcal_trackHighestPtChi2;
   vector<float>   *hgcal_trackHighestPtCutChi2;
   vector<float>   *hgcal_trackHighestPtCutChi2Eta;
   vector<float>   *hgcal_trackHighestPtCutChi2Phi;
   vector<float>   *hgcal_trackHighestPtCutChi2Chi2;
   vector<float>   *hgcal_trackmatchingdR;
   Int_t           cl3d_n;
   vector<float>   *cl3d_pt;
   vector<float>   *cl3d_energy;
   vector<float>   *cl3d_eta;
   vector<float>   *cl3d_phi;
   vector<int>     *cl3d_nclu;
   vector<float>   *cl3d_x;
   vector<float>   *cl3d_y;
   vector<int>     *cl3d_z;
   UShort_t        egN;
   vector<float>   *egE;
   vector<float>   *egEt;
   vector<float>   *egEta;
   vector<float>   *egPhi;
   vector<float>   *egGx;
   vector<float>   *egGy;
   vector<float>   *egGz;
   vector<short>   *egIEt;
   vector<short>   *egIEta;
   vector<short>   *egIPhi;
   vector<short>   *egIso;
   vector<short>   *egBx;
   vector<short>   *egTowerIPhi;
   vector<short>   *egTowerIEta;
   vector<short>   *egRawEt;
   vector<short>   *egIsoEt;
   vector<short>   *egFootprintEt;
   vector<short>   *egNTT;
   vector<short>   *egShape;
   vector<short>   *egTowerHoE;

   // List of branches
   TBranch        *b_nVtx;   //!
   TBranch        *b_nMeanPU;   //!
   TBranch        *b_genPartN;   //!
   TBranch        *b_genPartE;   //!
   TBranch        *b_genPartPt;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartCharge;   //!
   TBranch        *b_genPartId;   //!
   TBranch        *b_propgenPartE;   //!
   TBranch        *b_propgenPartPt;   //!
   TBranch        *b_propgenPartEta;   //!
   TBranch        *b_propgenPartPhi;   //!
   TBranch        *b_propgenPartCharge;   //!
   TBranch        *b_propgenPartId;   //!
   TBranch        *b_propgenPartx;   //!
   TBranch        *b_propgenParty;   //!
   TBranch        *b_propgenPartz;   //!
   TBranch        *b_simTrkN;   //!
   TBranch        *b_simTrkPt;   //!
   TBranch        *b_simTrkEta;   //!
   TBranch        *b_simTrkPhi;   //!
   TBranch        *b_simTrkId;   //!
   TBranch        *b_simTrkType;   //!
   TBranch        *b_simTrkVx;   //!
   TBranch        *b_simTrkVy;   //!
   TBranch        *b_simTrkVz;   //!
   TBranch        *b_simVx;   //!
   TBranch        *b_simVy;   //!
   TBranch        *b_simVz;   //!
   TBranch        *b_lastSimtkpt;   //!
   TBranch        *b_initialSimtkpt;   //!
   TBranch        *b_bremflag;   //!
   TBranch        *b_Brempos_radius;   //!
   TBranch        *b_Brem_eLoss;   //!
   TBranch        *b_Brem_ptLoss;   //!
   TBranch        *b_Brempos_x;   //!
   TBranch        *b_Brempos_y;   //!
   TBranch        *b_Brempos_z;   //!
   TBranch        *b_bHitLayer;   //!
   TBranch        *b_bHitLadder;   //!
   TBranch        *b_bHitModule;   //!
   TBranch        *b_fHitDisk;   //!
   TBranch        *b_fHitBlade;   //!
   TBranch        *b_fHitSide;   //!
   TBranch        *b_fHitPanel;   //!
   TBranch        *b_fHitModule;   //!
   TBranch        *b_bHitN;   //!
   TBranch        *b_fHitN;   //!
   TBranch        *b_fHitGx;   //!
   TBranch        *b_fHitGy;   //!
   TBranch        *b_fHitGz;   //!
   TBranch        *b_fHitPhi;   //!
   TBranch        *b_fHitEta;   //!
   TBranch        *b_fClSize;   //!
   TBranch        *b_fClSizeX;   //!
   TBranch        *b_fClSizeY;   //!
   TBranch        *b_bHitGx;   //!
   TBranch        *b_bHitGy;   //!
   TBranch        *b_bHitGz;   //!
   TBranch        *b_bHitPhi;   //!
   TBranch        *b_bHitEta;   //!
   TBranch        *b_bClSize;   //!
   TBranch        *b_bClSizeX;   //!
   TBranch        *b_bClSizeY;   //!
   TBranch        *b_egCrysN;   //!
   TBranch        *b_egCrysE;   //!
   TBranch        *b_egCrysEt;   //!
   TBranch        *b_egCrysEta;   //!
   TBranch        *b_egCrysPhi;   //!
   TBranch        *b_egCrysGx;   //!
   TBranch        *b_egCrysGy;   //!
   TBranch        *b_egCrysGz;   //!
   TBranch        *b_egCrysClusterN;   //!
   TBranch        *b_egCrysClusterE;   //!
   TBranch        *b_egCrysClusterEt;   //!
   TBranch        *b_egCrysClusterEta;   //!
   TBranch        *b_egCrysClusterPhi;   //!
   TBranch        *b_egCrysClusterGx;   //!
   TBranch        *b_egCrysClusterGy;   //!
   TBranch        *b_egCrysClusterGz;   //!
   TBranch        *b_egCrysClusterPGx;   //!
   TBranch        *b_egCrysClusterPGy;   //!
   TBranch        *b_egCrysClusterPGz;   //!
   TBranch        *b_isTrackMatched;   //!
   TBranch        *b_isoConeNTrack;   //!
   TBranch        *b_isoConePtTrack;   //!
   TBranch        *b_trackHighestPt;   //!
   TBranch        *b_trackHighestPtEta;   //!
   TBranch        *b_trackHighestPtPhi;   //!
   TBranch        *b_trackHighestPtChi2;   //!
   TBranch        *b_trackHighestPtCutChi2;   //!
   TBranch        *b_trackHighestPtCutChi2Eta;   //!
   TBranch        *b_trackHighestPtCutChi2Phi;   //!
   TBranch        *b_trackHighestPtCutChi2Chi2;   //!
   TBranch        *b_trackmatchingdR;   //!
   TBranch        *b_hgcal_isTrackMatched;   //!
   TBranch        *b_hgcal_isoConeNTrack;   //!
   TBranch        *b_hgcal_isoConePtTrack;   //!
   TBranch        *b_hgcal_trackHighestPt;   //!
   TBranch        *b_hgcal_trackHighestPtEta;   //!
   TBranch        *b_hgcal_trackHighestPtPhi;   //!
   TBranch        *b_hgcal_trackHighestPtChi2;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Eta;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Phi;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Chi2;   //!
   TBranch        *b_hgcal_trackmatchingdR;   //!
   TBranch        *b_cl3d_n;   //!
   TBranch        *b_cl3d_pt;   //!
   TBranch        *b_cl3d_energy;   //!
   TBranch        *b_cl3d_eta;   //!
   TBranch        *b_cl3d_phi;   //!
   TBranch        *b_cl3d_nclu;   //!
   TBranch        *b_cl3d_x;   //!
   TBranch        *b_cl3d_y;   //!
   TBranch        *b_cl3d_z;   //!
   TBranch        *b_egN;   //!
   TBranch        *b_egE;   //!
   TBranch        *b_egEt;   //!
   TBranch        *b_egEta;   //!
   TBranch        *b_egPhi;   //!
   TBranch        *b_egGx;   //!
   TBranch        *b_egGy;   //!
   TBranch        *b_egGz;   //!
   TBranch        *b_egIEt;   //!
   TBranch        *b_egIEta;   //!
   TBranch        *b_egIPhi;   //!
   TBranch        *b_egIso;   //!
   TBranch        *b_egBx;   //!
   TBranch        *b_egTowerIPhi;   //!
   TBranch        *b_egTowerIEta;   //!
   TBranch        *b_egRawEt;   //!
   TBranch        *b_egIsoEt;   //!
   TBranch        *b_egFootprintEt;   //!
   TBranch        *b_egNTT;   //!
   TBranch        *b_egShape;   //!
   TBranch        *b_egTowerHoE;   //!

   sw_treeMaker(TTree *tree=0);
   virtual ~sw_treeMaker();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef sw_treeMaker_cxx
sw_treeMaker::sw_treeMaker(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Resultall.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Resultall.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Resultall.root:/l1PiXTRKTree");
      dir->GetObject("l1pixtrktree",tree);

   }
   Init(tree);
}

sw_treeMaker::~sw_treeMaker()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sw_treeMaker::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sw_treeMaker::LoadTree(Long64_t entry)
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

void sw_treeMaker::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genPartE = 0;
   genPartPt = 0;
   genPartEta = 0;
   genPartPhi = 0;
   genPartCharge = 0;
   genPartId = 0;
   propgenPartE = 0;
   propgenPartPt = 0;
   propgenPartEta = 0;
   propgenPartPhi = 0;
   propgenPartCharge = 0;
   propgenPartId = 0;
   propgenPartx = 0;
   propgenParty = 0;
   propgenPartz = 0;
   simTrkPt = 0;
   simTrkEta = 0;
   simTrkPhi = 0;
   simTrkId = 0;
   simTrkType = 0;
   simTrkVx = 0;
   simTrkVy = 0;
   simTrkVz = 0;
   simVx = 0;
   simVy = 0;
   simVz = 0;
   Brempos_radius = 0;
   Brem_eLoss = 0;
   Brem_ptLoss = 0;
   Brempos_x = 0;
   Brempos_y = 0;
   Brempos_z = 0;
   bHitLayer = 0;
   bHitLadder = 0;
   bHitModule = 0;
   fHitDisk = 0;
   fHitBlade = 0;
   fHitSide = 0;
   fHitPanel = 0;
   fHitModule = 0;
   fHitGx = 0;
   fHitGy = 0;
   fHitGz = 0;
   fHitPhi = 0;
   fHitEta = 0;
   fClSize = 0;
   fClSizeX = 0;
   fClSizeY = 0;
   bHitGx = 0;
   bHitGy = 0;
   bHitGz = 0;
   bHitPhi = 0;
   bHitEta = 0;
   bClSize = 0;
   bClSizeX = 0;
   bClSizeY = 0;
   egCrysE = 0;
   egCrysEt = 0;
   egCrysEta = 0;
   egCrysPhi = 0;
   egCrysGx = 0;
   egCrysGy = 0;
   egCrysGz = 0;
   egCrysClusterE = 0;
   egCrysClusterEt = 0;
   egCrysClusterEta = 0;
   egCrysClusterPhi = 0;
   egCrysClusterGx = 0;
   egCrysClusterGy = 0;
   egCrysClusterGz = 0;
   egCrysClusterPGx = 0;
   egCrysClusterPGy = 0;
   egCrysClusterPGz = 0;
   isTrackMatched = 0;
   isoConeNTrack = 0;
   isoConePtTrack = 0;
   trackHighestPt = 0;
   trackHighestPtEta = 0;
   trackHighestPtPhi = 0;
   trackHighestPtChi2 = 0;
   trackHighestPtCutChi2 = 0;
   trackHighestPtCutChi2Eta = 0;
   trackHighestPtCutChi2Phi = 0;
   trackHighestPtCutChi2Chi2 = 0;
   trackmatchingdR = 0;
   hgcal_isTrackMatched = 0;
   hgcal_isoConeNTrack = 0;
   hgcal_isoConePtTrack = 0;
   hgcal_trackHighestPt = 0;
   hgcal_trackHighestPtEta = 0;
   hgcal_trackHighestPtPhi = 0;
   hgcal_trackHighestPtChi2 = 0;
   hgcal_trackHighestPtCutChi2 = 0;
   hgcal_trackHighestPtCutChi2Eta = 0;
   hgcal_trackHighestPtCutChi2Phi = 0;
   hgcal_trackHighestPtCutChi2Chi2 = 0;
   hgcal_trackmatchingdR = 0;
   cl3d_pt = 0;
   cl3d_energy = 0;
   cl3d_eta = 0;
   cl3d_phi = 0;
   cl3d_nclu = 0;
   cl3d_x = 0;
   cl3d_y = 0;
   cl3d_z = 0;
   egE = 0;
   egEt = 0;
   egEta = 0;
   egPhi = 0;
   egGx = 0;
   egGy = 0;
   egGz = 0;
   egIEt = 0;
   egIEta = 0;
   egIPhi = 0;
   egIso = 0;
   egBx = 0;
   egTowerIPhi = 0;
   egTowerIEta = 0;
   egRawEt = 0;
   egIsoEt = 0;
   egFootprintEt = 0;
   egNTT = 0;
   egShape = 0;
   egTowerHoE = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nMeanPU", &nMeanPU, &b_nMeanPU);
   fChain->SetBranchAddress("genPartN", &genPartN, &b_genPartN);
   fChain->SetBranchAddress("genPartE", &genPartE, &b_genPartE);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartCharge", &genPartCharge, &b_genPartCharge);
   fChain->SetBranchAddress("genPartId", &genPartId, &b_genPartId);
   fChain->SetBranchAddress("propgenPartE", &propgenPartE, &b_propgenPartE);
   fChain->SetBranchAddress("propgenPartPt", &propgenPartPt, &b_propgenPartPt);
   fChain->SetBranchAddress("propgenPartEta", &propgenPartEta, &b_propgenPartEta);
   fChain->SetBranchAddress("propgenPartPhi", &propgenPartPhi, &b_propgenPartPhi);
   fChain->SetBranchAddress("propgenPartCharge", &propgenPartCharge, &b_propgenPartCharge);
   fChain->SetBranchAddress("propgenPartId", &propgenPartId, &b_propgenPartId);
   fChain->SetBranchAddress("propgenPartx", &propgenPartx, &b_propgenPartx);
   fChain->SetBranchAddress("propgenParty", &propgenParty, &b_propgenParty);
   fChain->SetBranchAddress("propgenPartz", &propgenPartz, &b_propgenPartz);
   fChain->SetBranchAddress("simTrkN", &simTrkN, &b_simTrkN);
   fChain->SetBranchAddress("simTrkPt", &simTrkPt, &b_simTrkPt);
   fChain->SetBranchAddress("simTrkEta", &simTrkEta, &b_simTrkEta);
   fChain->SetBranchAddress("simTrkPhi", &simTrkPhi, &b_simTrkPhi);
   fChain->SetBranchAddress("simTrkId", &simTrkId, &b_simTrkId);
   fChain->SetBranchAddress("simTrkType", &simTrkType, &b_simTrkType);
   fChain->SetBranchAddress("simTrkVx", &simTrkVx, &b_simTrkVx);
   fChain->SetBranchAddress("simTrkVy", &simTrkVy, &b_simTrkVy);
   fChain->SetBranchAddress("simTrkVz", &simTrkVz, &b_simTrkVz);
   fChain->SetBranchAddress("simVx", &simVx, &b_simVx);
   fChain->SetBranchAddress("simVy", &simVy, &b_simVy);
   fChain->SetBranchAddress("simVz", &simVz, &b_simVz);
   fChain->SetBranchAddress("lastSimtkpt", &lastSimtkpt, &b_lastSimtkpt);
   fChain->SetBranchAddress("initialSimtkpt", &initialSimtkpt, &b_initialSimtkpt);
   fChain->SetBranchAddress("bremflag", &bremflag, &b_bremflag);
   fChain->SetBranchAddress("Brempos_radius", &Brempos_radius, &b_Brempos_radius);
   fChain->SetBranchAddress("Brem_eLoss", &Brem_eLoss, &b_Brem_eLoss);
   fChain->SetBranchAddress("Brem_ptLoss", &Brem_ptLoss, &b_Brem_ptLoss);
   fChain->SetBranchAddress("Brempos_x", &Brempos_x, &b_Brempos_x);
   fChain->SetBranchAddress("Brempos_y", &Brempos_y, &b_Brempos_y);
   fChain->SetBranchAddress("Brempos_z", &Brempos_z, &b_Brempos_z);
   fChain->SetBranchAddress("bHitLayer", &bHitLayer, &b_bHitLayer);
   fChain->SetBranchAddress("bHitLadder", &bHitLadder, &b_bHitLadder);
   fChain->SetBranchAddress("bHitModule", &bHitModule, &b_bHitModule);
   fChain->SetBranchAddress("fHitDisk", &fHitDisk, &b_fHitDisk);
   fChain->SetBranchAddress("fHitBlade", &fHitBlade, &b_fHitBlade);
   fChain->SetBranchAddress("fHitSide", &fHitSide, &b_fHitSide);
   fChain->SetBranchAddress("fHitPanel", &fHitPanel, &b_fHitPanel);
   fChain->SetBranchAddress("fHitModule", &fHitModule, &b_fHitModule);
   fChain->SetBranchAddress("bHitN", &bHitN, &b_bHitN);
   fChain->SetBranchAddress("fHitN", &fHitN, &b_fHitN);
   fChain->SetBranchAddress("fHitGx", &fHitGx, &b_fHitGx);
   fChain->SetBranchAddress("fHitGy", &fHitGy, &b_fHitGy);
   fChain->SetBranchAddress("fHitGz", &fHitGz, &b_fHitGz);
   fChain->SetBranchAddress("fHitPhi", &fHitPhi, &b_fHitPhi);
   fChain->SetBranchAddress("fHitEta", &fHitEta, &b_fHitEta);
   fChain->SetBranchAddress("fClSize", &fClSize, &b_fClSize);
   fChain->SetBranchAddress("fClSizeX", &fClSizeX, &b_fClSizeX);
   fChain->SetBranchAddress("fClSizeY", &fClSizeY, &b_fClSizeY);
   fChain->SetBranchAddress("bHitGx", &bHitGx, &b_bHitGx);
   fChain->SetBranchAddress("bHitGy", &bHitGy, &b_bHitGy);
   fChain->SetBranchAddress("bHitGz", &bHitGz, &b_bHitGz);
   fChain->SetBranchAddress("bHitPhi", &bHitPhi, &b_bHitPhi);
   fChain->SetBranchAddress("bHitEta", &bHitEta, &b_bHitEta);
   fChain->SetBranchAddress("bClSize", &bClSize, &b_bClSize);
   fChain->SetBranchAddress("bClSizeX", &bClSizeX, &b_bClSizeX);
   fChain->SetBranchAddress("bClSizeY", &bClSizeY, &b_bClSizeY);
   fChain->SetBranchAddress("egCrysN", &egCrysN, &b_egCrysN);
   fChain->SetBranchAddress("egCrysE", &egCrysE, &b_egCrysE);
   fChain->SetBranchAddress("egCrysEt", &egCrysEt, &b_egCrysEt);
   fChain->SetBranchAddress("egCrysEta", &egCrysEta, &b_egCrysEta);
   fChain->SetBranchAddress("egCrysPhi", &egCrysPhi, &b_egCrysPhi);
   fChain->SetBranchAddress("egCrysGx", &egCrysGx, &b_egCrysGx);
   fChain->SetBranchAddress("egCrysGy", &egCrysGy, &b_egCrysGy);
   fChain->SetBranchAddress("egCrysGz", &egCrysGz, &b_egCrysGz);
   fChain->SetBranchAddress("egCrysClusterN", &egCrysClusterN, &b_egCrysClusterN);
   fChain->SetBranchAddress("egCrysClusterE", &egCrysClusterE, &b_egCrysClusterE);
   fChain->SetBranchAddress("egCrysClusterEt", &egCrysClusterEt, &b_egCrysClusterEt);
   fChain->SetBranchAddress("egCrysClusterEta", &egCrysClusterEta, &b_egCrysClusterEta);
   fChain->SetBranchAddress("egCrysClusterPhi", &egCrysClusterPhi, &b_egCrysClusterPhi);
   fChain->SetBranchAddress("egCrysClusterGx", &egCrysClusterGx, &b_egCrysClusterGx);
   fChain->SetBranchAddress("egCrysClusterGy", &egCrysClusterGy, &b_egCrysClusterGy);
   fChain->SetBranchAddress("egCrysClusterGz", &egCrysClusterGz, &b_egCrysClusterGz);
   fChain->SetBranchAddress("egCrysClusterPGx", &egCrysClusterPGx, &b_egCrysClusterPGx);
   fChain->SetBranchAddress("egCrysClusterPGy", &egCrysClusterPGy, &b_egCrysClusterPGy);
   fChain->SetBranchAddress("egCrysClusterPGz", &egCrysClusterPGz, &b_egCrysClusterPGz);
   fChain->SetBranchAddress("isTrackMatched", &isTrackMatched, &b_isTrackMatched);
   fChain->SetBranchAddress("isoConeNTrack", &isoConeNTrack, &b_isoConeNTrack);
   fChain->SetBranchAddress("isoConePtTrack", &isoConePtTrack, &b_isoConePtTrack);
   fChain->SetBranchAddress("trackHighestPt", &trackHighestPt, &b_trackHighestPt);
   fChain->SetBranchAddress("trackHighestPtEta", &trackHighestPtEta, &b_trackHighestPtEta);
   fChain->SetBranchAddress("trackHighestPtPhi", &trackHighestPtPhi, &b_trackHighestPtPhi);
   fChain->SetBranchAddress("trackHighestPtChi2", &trackHighestPtChi2, &b_trackHighestPtChi2);
   fChain->SetBranchAddress("trackHighestPtCutChi2", &trackHighestPtCutChi2, &b_trackHighestPtCutChi2);
   fChain->SetBranchAddress("trackHighestPtCutChi2Eta", &trackHighestPtCutChi2Eta, &b_trackHighestPtCutChi2Eta);
   fChain->SetBranchAddress("trackHighestPtCutChi2Phi", &trackHighestPtCutChi2Phi, &b_trackHighestPtCutChi2Phi);
   fChain->SetBranchAddress("trackHighestPtCutChi2Chi2", &trackHighestPtCutChi2Chi2, &b_trackHighestPtCutChi2Chi2);
   fChain->SetBranchAddress("trackmatchingdR", &trackmatchingdR, &b_trackmatchingdR);
   fChain->SetBranchAddress("hgcal_isTrackMatched", &hgcal_isTrackMatched, &b_hgcal_isTrackMatched);
   fChain->SetBranchAddress("hgcal_isoConeNTrack", &hgcal_isoConeNTrack, &b_hgcal_isoConeNTrack);
   fChain->SetBranchAddress("hgcal_isoConePtTrack", &hgcal_isoConePtTrack, &b_hgcal_isoConePtTrack);
   fChain->SetBranchAddress("hgcal_trackHighestPt", &hgcal_trackHighestPt, &b_hgcal_trackHighestPt);
   fChain->SetBranchAddress("hgcal_trackHighestPtEta", &hgcal_trackHighestPtEta, &b_hgcal_trackHighestPtEta);
   fChain->SetBranchAddress("hgcal_trackHighestPtPhi", &hgcal_trackHighestPtPhi, &b_hgcal_trackHighestPtPhi);
   fChain->SetBranchAddress("hgcal_trackHighestPtChi2", &hgcal_trackHighestPtChi2, &b_hgcal_trackHighestPtChi2);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2", &hgcal_trackHighestPtCutChi2, &b_hgcal_trackHighestPtCutChi2);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Eta", &hgcal_trackHighestPtCutChi2Eta, &b_hgcal_trackHighestPtCutChi2Eta);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Phi", &hgcal_trackHighestPtCutChi2Phi, &b_hgcal_trackHighestPtCutChi2Phi);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Chi2", &hgcal_trackHighestPtCutChi2Chi2, &b_hgcal_trackHighestPtCutChi2Chi2);
   fChain->SetBranchAddress("hgcal_trackmatchingdR", &hgcal_trackmatchingdR, &b_hgcal_trackmatchingdR);
   fChain->SetBranchAddress("cl3d_n", &cl3d_n, &b_cl3d_n);
   fChain->SetBranchAddress("cl3d_pt", &cl3d_pt, &b_cl3d_pt);
   fChain->SetBranchAddress("cl3d_energy", &cl3d_energy, &b_cl3d_energy);
   fChain->SetBranchAddress("cl3d_eta", &cl3d_eta, &b_cl3d_eta);
   fChain->SetBranchAddress("cl3d_phi", &cl3d_phi, &b_cl3d_phi);
   fChain->SetBranchAddress("cl3d_nclu", &cl3d_nclu, &b_cl3d_nclu);
   fChain->SetBranchAddress("cl3d_x", &cl3d_x, &b_cl3d_x);
   fChain->SetBranchAddress("cl3d_y", &cl3d_y, &b_cl3d_y);
   fChain->SetBranchAddress("cl3d_z", &cl3d_z, &b_cl3d_z);
   fChain->SetBranchAddress("egN", &egN, &b_egN);
   fChain->SetBranchAddress("egE", &egE, &b_egE);
   fChain->SetBranchAddress("egEt", &egEt, &b_egEt);
   fChain->SetBranchAddress("egEta", &egEta, &b_egEta);
   fChain->SetBranchAddress("egPhi", &egPhi, &b_egPhi);
   fChain->SetBranchAddress("egGx", &egGx, &b_egGx);
   fChain->SetBranchAddress("egGy", &egGy, &b_egGy);
   fChain->SetBranchAddress("egGz", &egGz, &b_egGz);
   fChain->SetBranchAddress("egIEt", &egIEt, &b_egIEt);
   fChain->SetBranchAddress("egIEta", &egIEta, &b_egIEta);
   fChain->SetBranchAddress("egIPhi", &egIPhi, &b_egIPhi);
   fChain->SetBranchAddress("egIso", &egIso, &b_egIso);
   fChain->SetBranchAddress("egBx", &egBx, &b_egBx);
   fChain->SetBranchAddress("egTowerIPhi", &egTowerIPhi, &b_egTowerIPhi);
   fChain->SetBranchAddress("egTowerIEta", &egTowerIEta, &b_egTowerIEta);
   fChain->SetBranchAddress("egRawEt", &egRawEt, &b_egRawEt);
   fChain->SetBranchAddress("egIsoEt", &egIsoEt, &b_egIsoEt);
   fChain->SetBranchAddress("egFootprintEt", &egFootprintEt, &b_egFootprintEt);
   fChain->SetBranchAddress("egNTT", &egNTT, &b_egNTT);
   fChain->SetBranchAddress("egShape", &egShape, &b_egShape);
   fChain->SetBranchAddress("egTowerHoE", &egTowerHoE, &b_egTowerHoE);
   Notify();
}

Bool_t sw_treeMaker::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sw_treeMaker::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sw_treeMaker::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sw_treeMaker_cxx
