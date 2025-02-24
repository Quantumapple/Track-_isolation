//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 28 14:34:22 2018 by ROOT version 6.10/05
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: SingleMuonNoPU_1.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class test {
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
   vector<float>   *propgenElPartE;
   vector<float>   *propgenElPartPt;
   vector<float>   *propgenElPartEta;
   vector<float>   *propgenElPartPhi;
   vector<int>     *propgenElPartCharge;
   vector<float>   *propgenElPartx;
   vector<float>   *propgenElParty;
   vector<float>   *propgenElPartz;
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
   vector<float>   *propgenPoPartE;
   vector<float>   *propgenPoPartPt;
   vector<float>   *propgenPoPartEta;
   vector<float>   *propgenPoPartPhi;
   vector<int>     *propgenPoPartCharge;
   vector<float>   *propgenPoPartx;
   vector<float>   *propgenPoParty;
   vector<float>   *propgenPoPartz;
   vector<int>     *bRecHitLayer;
   vector<int>     *bRecHitLadder;
   vector<int>     *bRecHitModule;
   vector<int>     *fRecHitDisk;
   vector<int>     *fRecHitBlade;
   vector<int>     *fRecHitSide;
   vector<int>     *fRecHitPanel;
   vector<int>     *fRecHitModule;
   Int_t           bRecHitN;
   Int_t           fRecHitN;
   vector<float>   *fRecHitGx;
   vector<float>   *fRecHitGy;
   vector<float>   *fRecHitGz;
   vector<float>   *fRhSize;
   vector<float>   *fRhSizeX;
   vector<float>   *fRhSizeY;
   vector<float>   *bRecHitGx;
   vector<float>   *bRecHitGy;
   vector<float>   *bRecHitGz;
   vector<float>   *bRhSize;
   vector<float>   *bRhSizeX;
   vector<float>   *bRhSizeY;
   Int_t           bfastsimHitN;
   Int_t           ffastsimHitN;
   vector<int>     *bfastsimHitLayer;
   vector<float>   *bfastsimHitGx;
   vector<float>   *bfastsimHitGy;
   vector<float>   *bfastsimHitGz;
   vector<int>     *ffastsimHitLayer;
   vector<float>   *ffastsimHitGx;
   vector<float>   *ffastsimHitGy;
   vector<float>   *ffastsimHitGz;
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
   vector<float>   *cl3d_hovere;
   vector<int>     *cl3d_showerlength;
   vector<int>     *cl3d_coreshowerlength;
   vector<int>     *cl3d_firstlayer;
   vector<int>     *cl3d_maxlayer;
   vector<float>   *cl3d_seetot;
   vector<float>   *cl3d_seemax;
   vector<float>   *cl3d_spptot;
   vector<float>   *cl3d_sppmax;
   vector<float>   *cl3d_szz;
   vector<float>   *cl3d_srrtot;
   vector<float>   *cl3d_srrmax;
   vector<float>   *cl3d_srrmean;
   vector<float>   *cl3d_emaxe;
   UShort_t        egN;
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
   UInt_t          me0SegNum;
   vector<unsigned int> *me0SegDetId;
   vector<float>   *me0SegPosX;
   vector<float>   *me0SegPosY;
   vector<float>   *me0SegPosZ;
   vector<float>   *me0SegDirX;
   vector<float>   *me0SegDirY;
   vector<float>   *me0SegDirZ;
   vector<int>     *me0SegNumRecHit;
   vector<float>   *me0SegDeltaPhi;

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
   TBranch        *b_propgenElPartE;   //!
   TBranch        *b_propgenElPartPt;   //!
   TBranch        *b_propgenElPartEta;   //!
   TBranch        *b_propgenElPartPhi;   //!
   TBranch        *b_propgenElPartCharge;   //!
   TBranch        *b_propgenElPartx;   //!
   TBranch        *b_propgenElParty;   //!
   TBranch        *b_propgenElPartz;   //!
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
   TBranch        *b_propgenPoPartE;   //!
   TBranch        *b_propgenPoPartPt;   //!
   TBranch        *b_propgenPoPartEta;   //!
   TBranch        *b_propgenPoPartPhi;   //!
   TBranch        *b_propgenPoPartCharge;   //!
   TBranch        *b_propgenPoPartx;   //!
   TBranch        *b_propgenPoParty;   //!
   TBranch        *b_propgenPoPartz;   //!
   TBranch        *b_bRecHitLayer;   //!
   TBranch        *b_bRecHitLadder;   //!
   TBranch        *b_bRecHitModule;   //!
   TBranch        *b_fRecHitDisk;   //!
   TBranch        *b_fRecHitBlade;   //!
   TBranch        *b_fRecHitSide;   //!
   TBranch        *b_fRecHitPanel;   //!
   TBranch        *b_fRecHitModule;   //!
   TBranch        *b_bRecHitN;   //!
   TBranch        *b_fRecHitN;   //!
   TBranch        *b_fRecHitGx;   //!
   TBranch        *b_fRecHitGy;   //!
   TBranch        *b_fRecHitGz;   //!
   TBranch        *b_fRhSize;   //!
   TBranch        *b_fRhSizeX;   //!
   TBranch        *b_fRhSizeY;   //!
   TBranch        *b_bRecHitGx;   //!
   TBranch        *b_bRecHitGy;   //!
   TBranch        *b_bRecHitGz;   //!
   TBranch        *b_bRhSize;   //!
   TBranch        *b_bRhSizeX;   //!
   TBranch        *b_bRhSizeY;   //!
   TBranch        *b_bfastsimHitN;   //!
   TBranch        *b_ffastsimHitN;   //!
   TBranch        *b_bfastsimHitLayer;   //!
   TBranch        *b_bfastsimHitGx;   //!
   TBranch        *b_bfastsimHitGy;   //!
   TBranch        *b_bfastsimHitGz;   //!
   TBranch        *b_ffastsimHitLayer;   //!
   TBranch        *b_ffastsimHitGx;   //!
   TBranch        *b_ffastsimHitGy;   //!
   TBranch        *b_ffastsimHitGz;   //!
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
   TBranch        *b_cl3d_hovere;   //!
   TBranch        *b_cl3d_showerlength;   //!
   TBranch        *b_cl3d_coreshowerlength;   //!
   TBranch        *b_cl3d_firstlayer;   //!
   TBranch        *b_cl3d_maxlayer;   //!
   TBranch        *b_cl3d_seetot;   //!
   TBranch        *b_cl3d_seemax;   //!
   TBranch        *b_cl3d_spptot;   //!
   TBranch        *b_cl3d_sppmax;   //!
   TBranch        *b_cl3d_szz;   //!
   TBranch        *b_cl3d_srrtot;   //!
   TBranch        *b_cl3d_srrmax;   //!
   TBranch        *b_cl3d_srrmean;   //!
   TBranch        *b_cl3d_emaxe;   //!
   TBranch        *b_egN;   //!
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
   TBranch        *b_me0SegNum;   //!
   TBranch        *b_me0SegDetId;   //!
   TBranch        *b_me0SegPosX;   //!
   TBranch        *b_me0SegPosY;   //!
   TBranch        *b_me0SegPosZ;   //!
   TBranch        *b_me0SegDirX;   //!
   TBranch        *b_me0SegDirY;   //!
   TBranch        *b_me0SegDirZ;   //!
   TBranch        *b_me0SegNumRecHit;   //!
   TBranch        *b_me0SegDeltaPhi;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
   class track
   {
       public:
       float DR;
       float pos_x, pos_y, pos_z;
       track() { DR = 0.; pos_x = 0.; pos_y = 0.; pos_z = 0.; }
       track(float a, float b, float c, float d) { DR = a; pos_x = b; pos_y = c; pos_z = d; }
       static bool comp(const track &t1, const track &t2)
       {
	   return ( t1.DR < t2.DR );
       }
   };

   vector<track> bLayer1; 
   vector<track> bLayer2; 
   vector<track> bLayer3; 
   vector<track> bLayer4;

   vector<track> fDisk1;
   vector<track> fDisk2;
   vector<track> fDisk3;
   vector<track> fDisk4;
   vector<track> fDisk5;
   vector<track> fDisk6;

   vector<Float_t> saveParaCase1;
   vector<Float_t> saveParaCase2;
   vector<Float_t> saveParaCase3;
   vector<Float_t> saveParaCase4;

   void StorePixelHits(Int_t eta_region, Float_t propGenPhi, Float_t propGenEta);
   void CalculateParameters(Int_t eta_region);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
   //#ifdef SINGLE_TREE
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root");
      }
      //f->GetObject("l1PiXTRKTree/L1PiXTRKTree","");
      TDirectory * dir = (TDirectory*)f->Get("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root:/l1PiXTRKTree");
      dir->GetObject("L1PiXTRKTree",tree);

   //#else   

      TChain *chain = new TChain("l1PiXTRKTree/L1PiXTRKTree","");
      //chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_1.root/l1PiXTRKTree/L1PiXTRKTree");
      //chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/SingleMuonNoPU_2.root/l1PiXTRKTree/L1PiXTRKTree");
      chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0000/*.root/l1PiXTRKTree/L1PiXTRKTree");
      chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0001/*.root/l1PiXTRKTree/L1PiXTRKTree");
      chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0002/*.root/l1PiXTRKTree/L1PiXTRKTree");
      chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0003/*.root/l1PiXTRKTree/L1PiXTRKTree");
      chain->Add("/xrootd/store/user/jhong/SingleMu_Pt2to200_Eta3p0_CMSSW_9_3_7_NoPU_D17test/crab_Muon0802/180802_123200/0004/*.root/l1PiXTRKTree/L1PiXTRKTree");

   tree = chain;
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
   genPartE = 0;
   genPartPt = 0;
   genPartEta = 0;
   genPartPhi = 0;
   genPartCharge = 0;
   genPartId = 0;
   propgenElPartE = 0;
   propgenElPartPt = 0;
   propgenElPartEta = 0;
   propgenElPartPhi = 0;
   propgenElPartCharge = 0;
   propgenElPartx = 0;
   propgenElParty = 0;
   propgenElPartz = 0;
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
   propgenPoPartE = 0;
   propgenPoPartPt = 0;
   propgenPoPartEta = 0;
   propgenPoPartPhi = 0;
   propgenPoPartCharge = 0;
   propgenPoPartx = 0;
   propgenPoParty = 0;
   propgenPoPartz = 0;
   bRecHitLayer = 0;
   bRecHitLadder = 0;
   bRecHitModule = 0;
   fRecHitDisk = 0;
   fRecHitBlade = 0;
   fRecHitSide = 0;
   fRecHitPanel = 0;
   fRecHitModule = 0;
   fRecHitGx = 0;
   fRecHitGy = 0;
   fRecHitGz = 0;
   fRhSize = 0;
   fRhSizeX = 0;
   fRhSizeY = 0;
   bRecHitGx = 0;
   bRecHitGy = 0;
   bRecHitGz = 0;
   bRhSize = 0;
   bRhSizeX = 0;
   bRhSizeY = 0;
   bfastsimHitLayer = 0;
   bfastsimHitGx = 0;
   bfastsimHitGy = 0;
   bfastsimHitGz = 0;
   ffastsimHitLayer = 0;
   ffastsimHitGx = 0;
   ffastsimHitGy = 0;
   ffastsimHitGz = 0;
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
   cl3d_hovere = 0;
   cl3d_showerlength = 0;
   cl3d_coreshowerlength = 0;
   cl3d_firstlayer = 0;
   cl3d_maxlayer = 0;
   cl3d_seetot = 0;
   cl3d_seemax = 0;
   cl3d_spptot = 0;
   cl3d_sppmax = 0;
   cl3d_szz = 0;
   cl3d_srrtot = 0;
   cl3d_srrmax = 0;
   cl3d_srrmean = 0;
   cl3d_emaxe = 0;
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
   me0SegDetId = 0;
   me0SegPosX = 0;
   me0SegPosY = 0;
   me0SegPosZ = 0;
   me0SegDirX = 0;
   me0SegDirY = 0;
   me0SegDirZ = 0;
   me0SegNumRecHit = 0;
   me0SegDeltaPhi = 0;
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
   fChain->SetBranchAddress("propgenElPartE", &propgenElPartE, &b_propgenElPartE);
   fChain->SetBranchAddress("propgenElPartPt", &propgenElPartPt, &b_propgenElPartPt);
   fChain->SetBranchAddress("propgenElPartEta", &propgenElPartEta, &b_propgenElPartEta);
   fChain->SetBranchAddress("propgenElPartPhi", &propgenElPartPhi, &b_propgenElPartPhi);
   fChain->SetBranchAddress("propgenElPartCharge", &propgenElPartCharge, &b_propgenElPartCharge);
   fChain->SetBranchAddress("propgenElPartx", &propgenElPartx, &b_propgenElPartx);
   fChain->SetBranchAddress("propgenElParty", &propgenElParty, &b_propgenElParty);
   fChain->SetBranchAddress("propgenElPartz", &propgenElPartz, &b_propgenElPartz);
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
   fChain->SetBranchAddress("propgenPoPartE", &propgenPoPartE, &b_propgenPoPartE);
   fChain->SetBranchAddress("propgenPoPartPt", &propgenPoPartPt, &b_propgenPoPartPt);
   fChain->SetBranchAddress("propgenPoPartEta", &propgenPoPartEta, &b_propgenPoPartEta);
   fChain->SetBranchAddress("propgenPoPartPhi", &propgenPoPartPhi, &b_propgenPoPartPhi);
   fChain->SetBranchAddress("propgenPoPartCharge", &propgenPoPartCharge, &b_propgenPoPartCharge);
   fChain->SetBranchAddress("propgenPoPartx", &propgenPoPartx, &b_propgenPoPartx);
   fChain->SetBranchAddress("propgenPoParty", &propgenPoParty, &b_propgenPoParty);
   fChain->SetBranchAddress("propgenPoPartz", &propgenPoPartz, &b_propgenPoPartz);
   fChain->SetBranchAddress("bRecHitLayer", &bRecHitLayer, &b_bRecHitLayer);
   fChain->SetBranchAddress("bRecHitLadder", &bRecHitLadder, &b_bRecHitLadder);
   fChain->SetBranchAddress("bRecHitModule", &bRecHitModule, &b_bRecHitModule);
   fChain->SetBranchAddress("fRecHitDisk", &fRecHitDisk, &b_fRecHitDisk);
   fChain->SetBranchAddress("fRecHitBlade", &fRecHitBlade, &b_fRecHitBlade);
   fChain->SetBranchAddress("fRecHitSide", &fRecHitSide, &b_fRecHitSide);
   fChain->SetBranchAddress("fRecHitPanel", &fRecHitPanel, &b_fRecHitPanel);
   fChain->SetBranchAddress("fRecHitModule", &fRecHitModule, &b_fRecHitModule);
   fChain->SetBranchAddress("bRecHitN", &bRecHitN, &b_bRecHitN);
   fChain->SetBranchAddress("fRecHitN", &fRecHitN, &b_fRecHitN);
   fChain->SetBranchAddress("fRecHitGx", &fRecHitGx, &b_fRecHitGx);
   fChain->SetBranchAddress("fRecHitGy", &fRecHitGy, &b_fRecHitGy);
   fChain->SetBranchAddress("fRecHitGz", &fRecHitGz, &b_fRecHitGz);
   fChain->SetBranchAddress("fRhSize", &fRhSize, &b_fRhSize);
   fChain->SetBranchAddress("fRhSizeX", &fRhSizeX, &b_fRhSizeX);
   fChain->SetBranchAddress("fRhSizeY", &fRhSizeY, &b_fRhSizeY);
   fChain->SetBranchAddress("bRecHitGx", &bRecHitGx, &b_bRecHitGx);
   fChain->SetBranchAddress("bRecHitGy", &bRecHitGy, &b_bRecHitGy);
   fChain->SetBranchAddress("bRecHitGz", &bRecHitGz, &b_bRecHitGz);
   fChain->SetBranchAddress("bRhSize", &bRhSize, &b_bRhSize);
   fChain->SetBranchAddress("bRhSizeX", &bRhSizeX, &b_bRhSizeX);
   fChain->SetBranchAddress("bRhSizeY", &bRhSizeY, &b_bRhSizeY);
   fChain->SetBranchAddress("bfastsimHitN", &bfastsimHitN, &b_bfastsimHitN);
   fChain->SetBranchAddress("ffastsimHitN", &ffastsimHitN, &b_ffastsimHitN);
   fChain->SetBranchAddress("bfastsimHitLayer", &bfastsimHitLayer, &b_bfastsimHitLayer);
   fChain->SetBranchAddress("bfastsimHitGx", &bfastsimHitGx, &b_bfastsimHitGx);
   fChain->SetBranchAddress("bfastsimHitGy", &bfastsimHitGy, &b_bfastsimHitGy);
   fChain->SetBranchAddress("bfastsimHitGz", &bfastsimHitGz, &b_bfastsimHitGz);
   fChain->SetBranchAddress("ffastsimHitLayer", &ffastsimHitLayer, &b_ffastsimHitLayer);
   fChain->SetBranchAddress("ffastsimHitGx", &ffastsimHitGx, &b_ffastsimHitGx);
   fChain->SetBranchAddress("ffastsimHitGy", &ffastsimHitGy, &b_ffastsimHitGy);
   fChain->SetBranchAddress("ffastsimHitGz", &ffastsimHitGz, &b_ffastsimHitGz);
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
   fChain->SetBranchAddress("cl3d_hovere", &cl3d_hovere, &b_cl3d_hovere);
   fChain->SetBranchAddress("cl3d_showerlength", &cl3d_showerlength, &b_cl3d_showerlength);
   fChain->SetBranchAddress("cl3d_coreshowerlength", &cl3d_coreshowerlength, &b_cl3d_coreshowerlength);
   fChain->SetBranchAddress("cl3d_firstlayer", &cl3d_firstlayer, &b_cl3d_firstlayer);
   fChain->SetBranchAddress("cl3d_maxlayer", &cl3d_maxlayer, &b_cl3d_maxlayer);
   fChain->SetBranchAddress("cl3d_seetot", &cl3d_seetot, &b_cl3d_seetot);
   fChain->SetBranchAddress("cl3d_seemax", &cl3d_seemax, &b_cl3d_seemax);
   fChain->SetBranchAddress("cl3d_spptot", &cl3d_spptot, &b_cl3d_spptot);
   fChain->SetBranchAddress("cl3d_sppmax", &cl3d_sppmax, &b_cl3d_sppmax);
   fChain->SetBranchAddress("cl3d_szz", &cl3d_szz, &b_cl3d_szz);
   fChain->SetBranchAddress("cl3d_srrtot", &cl3d_srrtot, &b_cl3d_srrtot);
   fChain->SetBranchAddress("cl3d_srrmax", &cl3d_srrmax, &b_cl3d_srrmax);
   fChain->SetBranchAddress("cl3d_srrmean", &cl3d_srrmean, &b_cl3d_srrmean);
   fChain->SetBranchAddress("cl3d_emaxe", &cl3d_emaxe, &b_cl3d_emaxe);
   fChain->SetBranchAddress("egN", &egN, &b_egN);
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
   fChain->SetBranchAddress("me0SegNum", &me0SegNum, &b_me0SegNum);
   fChain->SetBranchAddress("me0SegDetId", &me0SegDetId, &b_me0SegDetId);
   fChain->SetBranchAddress("me0SegPosX", &me0SegPosX, &b_me0SegPosX);
   fChain->SetBranchAddress("me0SegPosY", &me0SegPosY, &b_me0SegPosY);
   fChain->SetBranchAddress("me0SegPosZ", &me0SegPosZ, &b_me0SegPosZ);
   fChain->SetBranchAddress("me0SegDirX", &me0SegDirX, &b_me0SegDirX);
   fChain->SetBranchAddress("me0SegDirY", &me0SegDirY, &b_me0SegDirY);
   fChain->SetBranchAddress("me0SegDirZ", &me0SegDirZ, &b_me0SegDirZ);
   fChain->SetBranchAddress("me0SegNumRecHit", &me0SegNumRecHit, &b_me0SegNumRecHit);
   fChain->SetBranchAddress("me0SegDeltaPhi", &me0SegDeltaPhi, &b_me0SegDeltaPhi);
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
void test::StorePixelHits(Int_t eta_region, Float_t propGenPhi, Float_t propGenEta)
{
      for(Int_t i = 0; i < bRecHitN; i++)
      {
          TVector3 pixel;
          pixel.SetXYZ(bRecHitGx->at(i) - simVx->at(0), bRecHitGy->at(i) - simVy->at(0), bRecHitGz->at(i) - simVz->at(0));
          Float_t pixelPhi = pixel.Phi();
          Float_t pixelEta = pixel.PseudoRapidity();
          Float_t deltaPhi = pixelPhi - propGenPhi;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
          Float_t deltaEta = pixelEta - propGenEta;
          Float_t deltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta,2));

          if( eta_region <= 3 )
          {
              if( bRecHitLayer->at(i) == 1 )  bLayer1.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
              if( bRecHitLayer->at(i) == 2 )  bLayer2.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          }
          if( eta_region <= 2 )
          {
              if( bRecHitLayer->at(i) == 3 )  bLayer3.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          }
          if( eta_region == 1 )
          {
              if( bRecHitLayer->at(i) == 4 )  bLayer4.push_back(track(deltaR, bRecHitGx->at(i),bRecHitGy->at(i),bRecHitGz->at(i)));
          }
      }
      
      for(Int_t i = 0; i < fRecHitN; i++)
      {
          TVector3 pixel;
          pixel.SetXYZ(fRecHitGx->at(i) - simVx->at(0), fRecHitGy->at(i) - simVy->at(0), fRecHitGz->at(i) - simVz->at(0));
          Float_t pixelPhi = pixel.Phi();
          Float_t pixelEta = pixel.PseudoRapidity();
          Float_t deltaPhi = pixelPhi - propGenPhi;
          if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
          if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
          Float_t deltaEta = pixelEta - propGenEta;
          Float_t deltaR = sqrt(pow(deltaPhi,2)+pow(deltaEta,2));

          if( eta_region >= 2 && eta_region <= 4 )
          {
              if( fRecHitDisk->at(i) == 1 )  fDisk1.push_back(track(deltaR, fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)));
          }
          if( eta_region >= 3 && eta_region <= 5 )
          {
              if( fRecHitDisk->at(i) == 2 )  fDisk2.push_back(track(deltaR, fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)));
          }
          if( eta_region >= 4 )
          {
              if( fRecHitDisk->at(i) == 3 )  fDisk3.push_back(track(deltaR, fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)));
              if( fRecHitDisk->at(i) == 4 )  fDisk4.push_back(track(deltaR, fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)));
          }
          if( eta_region >= 5 )
          {
              if( fRecHitDisk->at(i) == 5 )  fDisk5.push_back(track(deltaR, fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)));
          }
          if( eta_region == 6 )
          {
              if( fRecHitDisk->at(i) == 6 )  fDisk6.push_back(track(deltaR, fRecHitGx->at(i),fRecHitGy->at(i),fRecHitGz->at(i)));
          }
      }
}

void test::CalculateParameters(Int_t eta_region)
{
    if( eta_region == 1 )
    {
        //// -------- L123 --------- //// 
        if( bLayer1.size() != 0 && bLayer2.size() != 0 && bLayer3.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = bLayer2.begin(); it2 != bLayer2.end(); ++it2)
                {
                    TVector3 pp1; // L12
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = bLayer3.begin(); it3 != bLayer3.end(); ++it3 )
                    {
                        TVector3 pp2; // L23
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L13
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVL3
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL1223 = pp1_eta - pp2_eta;
                        Float_t dEtaL1213 = pp1_eta - pp3_eta;
                        Float_t dEtaL1323 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL13 = pv3_eta - pv1_eta;
                        Float_t dEtaPVL23 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2L3
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase1.push_back(dEtaL1213);
                        saveParaCase1.push_back(dEtaL1223);
                        saveParaCase1.push_back(dEtaL1323);
                        saveParaCase1.push_back(dEtaPVL23);
                        saveParaCase1.push_back(dEtaPVL13);
                        saveParaCase1.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L124 --------- //// 
        if( bLayer1.size() != 0 && bLayer2.size() != 0 && bLayer4.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = bLayer2.begin(); it2 != bLayer2.end(); ++it2)
                {
                    TVector3 pp1; // L12
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = bLayer4.begin(); it3 != bLayer4.end(); ++it3 )
                    {
                        TVector3 pp2; // L24
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L14
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVL4
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL1224 = pp1_eta - pp2_eta;
                        Float_t dEtaL1214 = pp1_eta - pp3_eta;
                        Float_t dEtaL1424 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL14 = pv3_eta - pv1_eta;
                        Float_t dEtaPVL24 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2L4
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase2.push_back(dEtaL1214);
                        saveParaCase2.push_back(dEtaL1224);
                        saveParaCase2.push_back(dEtaL1424);
                        saveParaCase2.push_back(dEtaPVL24);
                        saveParaCase2.push_back(dEtaPVL14);
                        saveParaCase2.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L134 --------- //// 
        if( bLayer1.size() != 0 && bLayer3.size() != 0 && bLayer4.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = bLayer3.begin(); it2 != bLayer3.end(); ++it2)
                {
                    TVector3 pp1; // L13
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL3
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = bLayer4.begin(); it3 != bLayer4.end(); ++it3 )
                    {
                        TVector3 pp2; // L34
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L14
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVL4
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL1334 = pp1_eta - pp2_eta;
                        Float_t dEtaL1314 = pp1_eta - pp3_eta;
                        Float_t dEtaL1434 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL34 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L3 - L3L4
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase3.push_back(dEtaL1314);
                        saveParaCase3.push_back(dEtaL1334);
                        saveParaCase3.push_back(dEtaL1434);
                        saveParaCase3.push_back(dEtaPVL34);
                        saveParaCase3.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L234 --------- //// 
        if( bLayer2.size() != 0 && bLayer3.size() != 0 && bLayer4.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer2.begin(); it1 != bLayer2.end(); ++it1)
            {
                TVector3 pv1; // PVL2
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = bLayer3.begin(); it2 != bLayer3.end(); ++it2)
                {
                    TVector3 pp1; // L23
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    for(std::vector<track>::iterator it3 = bLayer4.begin(); it3 != bLayer4.end(); ++it3 )
                    {
                        TVector3 pp2; // L34
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L24
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        Float_t dEtaL2334 = pp1_eta - pp2_eta;
                        Float_t dEtaL2324 = pp1_eta - pp3_eta;
                        Float_t dEtaL2434 = pp3_eta - pp2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL2 - L2L3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L2L3 - L3L4
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase4.push_back(dEtaL2324);
                        saveParaCase4.push_back(dEtaL2334);
                        saveParaCase4.push_back(dEtaL2434);
                        saveParaCase4.push_back(ddPhi);
                    }
                }
            }
        }
    }

    if( eta_region == 2 )
    {
        //// -------- L123 --------- //// 
        if( bLayer1.size() != 0 && bLayer2.size() != 0 && bLayer3.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = bLayer2.begin(); it2 != bLayer2.end(); ++it2)
                {
                    TVector3 pp1; // L12
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = bLayer3.begin(); it3 != bLayer3.end(); ++it3 )
                    {
                        TVector3 pp2; // L23
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L13
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVL3
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL1223 = pp1_eta - pp2_eta;
                        Float_t dEtaL1213 = pp1_eta - pp3_eta;
                        Float_t dEtaL1323 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL13 = pv3_eta - pv1_eta;
                        Float_t dEtaPVL23 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2L3
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase1.push_back(dEtaL1213);
                        saveParaCase1.push_back(dEtaL1223);
                        saveParaCase1.push_back(dEtaL1323);
                        saveParaCase1.push_back(dEtaPVL23);
                        saveParaCase1.push_back(dEtaPVL13);
                        saveParaCase1.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L12D1 --------- //// 
        if( bLayer1.size() != 0 && bLayer2.size() != 0 && fDisk1.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = bLayer2.begin(); it2 != bLayer2.end(); ++it2)
                {
                    TVector3 pp1; // L12
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk1.begin(); it3 != fDisk1.end(); ++it3 )
                    {
                        TVector3 pp2; // L2D1
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L1D1
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD1
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL12L2D1 = pp1_eta - pp2_eta;
                        Float_t dEtaL12L1D1 = pp1_eta - pp3_eta;
                        Float_t dEtaL1D1L2D1 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL1D1 = pv3_eta - pv1_eta;
                        Float_t dEtaPVL2D1 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2D1
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase2.push_back(dEtaL12L1D1);
                        saveParaCase2.push_back(dEtaL12L2D1);
                        saveParaCase2.push_back(dEtaL1D1L2D1);
                        saveParaCase2.push_back(dEtaPVL2D1);
                        saveParaCase2.push_back(dEtaPVL1D1);
                        saveParaCase2.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L13D1 --------- //// 
        if( bLayer1.size() != 0 && bLayer3.size() != 0 && fDisk1.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = bLayer3.begin(); it2 != bLayer3.end(); ++it2)
                {
                    TVector3 pp1; // L13
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL3
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk1.begin(); it3 != fDisk1.end(); ++it3 )
                    {
                        TVector3 pp2; // L3D1
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L1D1
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD1
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL13L3D1 = pp1_eta - pp2_eta;
                        Float_t dEtaL13L1D1 = pp1_eta - pp3_eta;
                        Float_t dEtaL1D1L3D1 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL3D1 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L3 - L3D1
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase3.push_back(dEtaL13L1D1);
                        saveParaCase3.push_back(dEtaL13L3D1);
                        saveParaCase3.push_back(dEtaL1D1L3D1);
                        saveParaCase3.push_back(dEtaPVL3D1);
                        saveParaCase3.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L23D1 --------- //// 
        if( bLayer2.size() != 0 && bLayer3.size() != 0 && fDisk1.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer2.begin(); it1 != bLayer2.end(); ++it1)
            {
                TVector3 pv1; // PVL2
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = bLayer3.begin(); it2 != bLayer3.end(); ++it2)
                {
                    TVector3 pp1; // L23
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    for(std::vector<track>::iterator it3 = fDisk1.begin(); it3 != fDisk1.end(); ++it3 )
                    {
                        TVector3 pp2; // L3D1
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L2D1
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        Float_t dEtaL23L3D1 = pp1_eta - pp2_eta;
                        Float_t dEtaL23L2D1 = pp1_eta - pp3_eta;
                        Float_t dEtaL2D1L3D1 = pp3_eta - pp2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL2 - L2L3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L2L3 - L3D1
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase4.push_back(dEtaL23L2D1);
                        saveParaCase4.push_back(dEtaL23L3D1);
                        saveParaCase4.push_back(dEtaL2D1L3D1);
                        saveParaCase4.push_back(ddPhi);
                    }
                }
            }
        }
    }

    if( eta_region == 3 )
    {
        //// -------- L12D1 --------- //// 
        if( bLayer1.size() != 0 && bLayer2.size() != 0 && fDisk1.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = bLayer2.begin(); it2 != bLayer2.end(); ++it2)
                {
                    TVector3 pp1; // L12
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk1.begin(); it3 != fDisk1.end(); ++it3 )
                    {
                        TVector3 pp2; // L2D1
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L1D1
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD1
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL12L2D1 = pp1_eta - pp2_eta;
                        Float_t dEtaL12L1D1 = pp1_eta - pp3_eta;
                        Float_t dEtaL1D1L2D1 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL1D1 = pv3_eta - pv1_eta;
                        Float_t dEtaPVL2D1 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2D1
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase1.push_back(dEtaL12L1D1);
                        saveParaCase1.push_back(dEtaL12L2D1);
                        saveParaCase1.push_back(dEtaL1D1L2D1);
                        saveParaCase1.push_back(dEtaPVL2D1);
                        saveParaCase1.push_back(dEtaPVL1D1);
                        saveParaCase1.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L12D2 --------- //// 
        if( bLayer1.size() != 0 && bLayer2.size() != 0 && fDisk2.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = bLayer2.begin(); it2 != bLayer2.end(); ++it2)
                {
                    TVector3 pp1; // L12
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVL2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk2.begin(); it3 != fDisk2.end(); ++it3 )
                    {
                        TVector3 pp2; // L2D2
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L1D2
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD2
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL12L2D2 = pp1_eta - pp2_eta;
                        Float_t dEtaL12L1D2 = pp1_eta - pp3_eta;
                        Float_t dEtaL1D1L2D2 = pp3_eta - pp2_eta;
                        Float_t dEtaPVL1D2 = pv3_eta - pv1_eta;
                        Float_t dEtaPVL2D2 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1L2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1L2 - L2D2
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase2.push_back(dEtaL12L1D2);
                        saveParaCase2.push_back(dEtaL12L2D2);
                        saveParaCase2.push_back(dEtaL1D1L2D2);
                        saveParaCase2.push_back(dEtaPVL2D2);
                        saveParaCase2.push_back(dEtaPVL1D2);
                        saveParaCase2.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L1D12 --------- //// 
        if( bLayer1.size() != 0 && fDisk1.size() != 0 && fDisk2.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer1.begin(); it1 != bLayer1.end(); ++it1)
            {
                TVector3 pv1; // PVL1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk1.begin(); it2 != fDisk1.end(); ++it2)
                {
                    TVector3 pp1; // L1D1
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD1
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk2.begin(); it3 != fDisk2.end(); ++it3 )
                    {
                        TVector3 pp2; // L3D2
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L1D2
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD2
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaL1D1D12 = pp1_eta - pp2_eta;
                        Float_t dEtaL1D1L1D2 = pp1_eta - pp3_eta;
                        Float_t dEtaL1D2D12 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD12 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL1 - L1D1
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L1D1 - D1D2
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase3.push_back(dEtaL1D1L1D2);
                        saveParaCase3.push_back(dEtaL1D1D12);
                        saveParaCase3.push_back(dEtaL1D2D12);
                        saveParaCase3.push_back(dEtaPVD12);
                        saveParaCase3.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- L2D12 --------- //// 
        if( bLayer2.size() != 0 && fDisk1.size() != 0 && fDisk2.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = bLayer2.begin(); it1 != bLayer2.end(); ++it1)
            {
                TVector3 pv1; // PVL2
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk1.begin(); it2 != fDisk1.end(); ++it2)
                {
                    TVector3 pp1; // L2D1
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    for(std::vector<track>::iterator it3 = fDisk2.begin(); it3 != fDisk2.end(); ++it3 )
                    {
                        TVector3 pp2; // D1D2
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // L2D2
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        Float_t dEtaL2D1L2D2 = pp1_eta - pp2_eta;
                        Float_t dEtaL2D1D12 = pp1_eta - pp3_eta;
                        Float_t dEtaL2D2D12 = pp3_eta - pp2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVL2 - L2D1
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = L2D1 - D1D2
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase4.push_back(dEtaL2D1L2D2);
                        saveParaCase4.push_back(dEtaL2D1D12);
                        saveParaCase4.push_back(dEtaL2D2D12);
                        saveParaCase4.push_back(ddPhi);
                    }
                }
            }
        }
    }

    if( eta_region == 4 )
    {
        //// -------- D123 --------- //// 
        if( fDisk1.size() != 0 && fDisk2.size() != 0 && fDisk3.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk1.begin(); it1 != fDisk1.end(); ++it1)
            {
                TVector3 pv1; // PVD1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = fDisk2.begin(); it2 != fDisk2.end(); ++it2)
                {
                    TVector3 pp1; // D1D2
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk3.begin(); it3 != fDisk3.end(); ++it3 )
                    {
                        TVector3 pp2; // D2D3
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D1D3
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD3
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD1223 = pp1_eta - pp2_eta;
                        Float_t dEtaD1213 = pp1_eta - pp3_eta;
                        Float_t dEtaD1323 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD31 = pv3_eta - pv1_eta;
                        Float_t dEtaPVD32 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD1 - D1D2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D1D2 - D2D3
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase1.push_back(dEtaD1213);
                        saveParaCase1.push_back(dEtaD1223);
                        saveParaCase1.push_back(dEtaD1323);
                        saveParaCase1.push_back(dEtaPVD32);
                        saveParaCase1.push_back(dEtaPVD31);
                        saveParaCase1.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D124 --------- //// 
        if( fDisk1.size() != 0 && fDisk2.size() != 0 && fDisk4.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk1.begin(); it1 != fDisk1.end(); ++it1)
            {
                TVector3 pv1; // PVD1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = fDisk2.begin(); it2 != fDisk2.end(); ++it2)
                {
                    TVector3 pp1; // D1D2
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD2
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk4.begin(); it3 != fDisk4.end(); ++it3 )
                    {
                        TVector3 pp2; // D2D4
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D1D4
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD4
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD1224 = pp1_eta - pp2_eta;
                        Float_t dEtaD1214 = pp1_eta - pp3_eta;
                        Float_t dEtaD1424 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD41 = pv3_eta - pv1_eta;
                        Float_t dEtaPVD42 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD1 - D1D2
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D1D2 - D2D4
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase2.push_back(dEtaD1214);
                        saveParaCase2.push_back(dEtaD1224);
                        saveParaCase2.push_back(dEtaD1424);
                        saveParaCase2.push_back(dEtaPVD42);
                        saveParaCase2.push_back(dEtaPVD41);
                        saveParaCase2.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D134 --------- //// 
        if( fDisk1.size() != 0 && fDisk3.size() != 0 && fDisk4.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk1.begin(); it1 != fDisk1.end(); ++it1)
            {
                TVector3 pv1; // PVD1
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk3.begin(); it2 != fDisk3.end(); ++it2)
                {
                    TVector3 pp1; // D1D3
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD3
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk4.begin(); it3 != fDisk4.end(); ++it3 )
                    {
                        TVector3 pp2; // D3D4
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D1D4
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD4
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD1334 = pp1_eta - pp2_eta;
                        Float_t dEtaD1314 = pp1_eta - pp3_eta;
                        Float_t dEtaD1434 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD43 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD1 - D1D3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D1D3 - D3D4
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase3.push_back(dEtaD1314);
                        saveParaCase3.push_back(dEtaD1334);
                        saveParaCase3.push_back(dEtaD1434);
                        saveParaCase3.push_back(dEtaPVD43);
                        saveParaCase3.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D234 --------- //// 
        if( fDisk2.size() != 0 && fDisk3.size() != 0 && fDisk4.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk2.begin(); it1 != fDisk2.end(); ++it1)
            {
                TVector3 pv1; // PVD2
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk3.begin(); it2 != fDisk3.end(); ++it2)
                {
                    TVector3 pp1; // D2D3
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    for(std::vector<track>::iterator it3 = fDisk4.begin(); it3 != fDisk4.end(); ++it3 )
                    {
                        TVector3 pp2; // D3D4
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D2D4
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        Float_t dEtaD2334 = pp1_eta - pp2_eta;
                        Float_t dEtaD2324 = pp1_eta - pp3_eta;
                        Float_t dEtaD2434 = pp3_eta - pp2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD2 - D2D3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D2D3 - D3D4
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase4.push_back(dEtaD2324);
                        saveParaCase4.push_back(dEtaD2334);
                        saveParaCase4.push_back(dEtaD2434);
                        saveParaCase4.push_back(ddPhi);
                    }
                }
            }
        }
    }

    if( eta_region == 5 )
    {
        //// -------- D234 --------- //// 
        if( fDisk2.size() != 0 && fDisk3.size() != 0 && fDisk4.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk2.begin(); it1 != fDisk2.end(); ++it1)
            {
                TVector3 pv1; // PVD2
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = fDisk3.begin(); it2 != fDisk3.end(); ++it2)
                {
                    TVector3 pp1; // D2D3
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD3
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk4.begin(); it3 != fDisk4.end(); ++it3 )
                    {
                        TVector3 pp2; // D3D4
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D2D4
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD4
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD2334 = pp1_eta - pp2_eta;
                        Float_t dEtaD2324 = pp1_eta - pp3_eta;
                        Float_t dEtaD2434 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD42 = pv3_eta - pv1_eta;
                        Float_t dEtaPVD43 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD2 - D2D3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D2D3 - D3D4
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase1.push_back(dEtaD2324);
                        saveParaCase1.push_back(dEtaD2334);
                        saveParaCase1.push_back(dEtaD2434);
                        saveParaCase1.push_back(dEtaPVD43);
                        saveParaCase1.push_back(dEtaPVD42);
                        saveParaCase1.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D235 --------- //// 
        if( fDisk2.size() != 0 && fDisk3.size() != 0 && fDisk5.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk2.begin(); it1 != fDisk2.end(); ++it1)
            {
                TVector3 pv1; // PVD2
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = fDisk3.begin(); it2 != fDisk3.end(); ++it2)
                {
                    TVector3 pp1; // D2D3
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD3
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk5.begin(); it3 != fDisk5.end(); ++it3 )
                    {
                        TVector3 pp2; // D3D5
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D2D5
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD5
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD2335 = pp1_eta - pp2_eta;
                        Float_t dEtaD2325 = pp1_eta - pp3_eta;
                        Float_t dEtaD2535 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD52 = pv3_eta - pv1_eta;
                        Float_t dEtaPVD53 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD2 - D2D3
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D2D3 - D3D5
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase2.push_back(dEtaD2325);
                        saveParaCase2.push_back(dEtaD2335);
                        saveParaCase2.push_back(dEtaD2535);
                        saveParaCase2.push_back(dEtaPVD53);
                        saveParaCase2.push_back(dEtaPVD52);
                        saveParaCase2.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D245 --------- //// 
        if( fDisk2.size() != 0 && fDisk4.size() != 0 && fDisk5.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk2.begin(); it1 != fDisk2.end(); ++it1)
            {
                TVector3 pv1; // PVD2
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk4.begin(); it2 != fDisk4.end(); ++it2)
                {
                    TVector3 pp1; // D2D4
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD4
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk5.begin(); it3 != fDisk5.end(); ++it3 )
                    {
                        TVector3 pp2; // D4D5
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D2D5
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD5
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD2445 = pp1_eta - pp2_eta;
                        Float_t dEtaD2425 = pp1_eta - pp3_eta;
                        Float_t dEtaD2545 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD54 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD2 - D2D4
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D2D4 - D4D5
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase3.push_back(dEtaD2425);
                        saveParaCase3.push_back(dEtaD2445);
                        saveParaCase3.push_back(dEtaD2545);
                        saveParaCase3.push_back(dEtaPVD54);
                        saveParaCase3.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D345 --------- //// 
        if( fDisk3.size() != 0 && fDisk4.size() != 0 && fDisk5.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk3.begin(); it1 != fDisk3.end(); ++it1)
            {
                TVector3 pv1; // PVD3
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk4.begin(); it2 != fDisk4.end(); ++it2)
                {
                    TVector3 pp1; // D3D4
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    for(std::vector<track>::iterator it3 = fDisk5.begin(); it3 != fDisk5.end(); ++it3 )
                    {
                        TVector3 pp2; // D4D5
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D3D5
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        Float_t dEtaD3445 = pp1_eta - pp2_eta;
                        Float_t dEtaD3435 = pp1_eta - pp3_eta;
                        Float_t dEtaD3545 = pp3_eta - pp2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD3 - D3D4
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D3D4 - D4D5
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase4.push_back(dEtaD3435);
                        saveParaCase4.push_back(dEtaD3445);
                        saveParaCase4.push_back(dEtaD3545);
                        saveParaCase4.push_back(ddPhi);
                    }
                }
            }
        }
    }

    if( eta_region == 6 )
    {
        //// -------- D345 --------- //// 
        if( fDisk3.size() != 0 && fDisk4.size() != 0 && fDisk5.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk3.begin(); it1 != fDisk3.end(); ++it1)
            {
                TVector3 pv1; // PVD3
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = fDisk4.begin(); it2 != fDisk4.end(); ++it2)
                {
                    TVector3 pp1; // D3D4
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD4
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk5.begin(); it3 != fDisk5.end(); ++it3 )
                    {
                        TVector3 pp2; // D4D5
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D3D5
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD5
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD3445 = pp1_eta - pp2_eta;
                        Float_t dEtaD3435 = pp1_eta - pp3_eta;
                        Float_t dEtaD3545 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD53 = pv3_eta - pv1_eta;
                        Float_t dEtaPVD54 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD3 - D3D4
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D3D4 - D4D5
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase1.push_back(dEtaD3435);
                        saveParaCase1.push_back(dEtaD3445);
                        saveParaCase1.push_back(dEtaD3545);
                        saveParaCase1.push_back(dEtaPVD54);
                        saveParaCase1.push_back(dEtaPVD53);
                        saveParaCase1.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D346 --------- //// 
        if( fDisk3.size() != 0 && fDisk4.size() != 0 && fDisk6.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk3.begin(); it1 != fDisk3.end(); ++it1)
            {
                TVector3 pv1; // PVD3
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_eta = pv1.Eta();
                Float_t pv1_phi = pv1.Phi();

                for(std::vector<track>::iterator it2 = fDisk4.begin(); it2 != fDisk4.end(); ++it2)
                {
                    TVector3 pp1; // D3D4
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD4
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk6.begin(); it3 != fDisk6.end(); ++it3 )
                    {
                        TVector3 pp2; // D4D6
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D3D6
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD6
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD3446 = pp1_eta - pp2_eta;
                        Float_t dEtaD3436 = pp1_eta - pp3_eta;
                        Float_t dEtaD3646 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD63 = pv3_eta - pv1_eta;
                        Float_t dEtaPVD64 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD3 - D3D4
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D3D4 - D4D6
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase2.push_back(dEtaD3436);
                        saveParaCase2.push_back(dEtaD3446);
                        saveParaCase2.push_back(dEtaD3646);
                        saveParaCase2.push_back(dEtaPVD64);
                        saveParaCase2.push_back(dEtaPVD63);
                        saveParaCase2.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D356 --------- //// 
        if( fDisk3.size() != 0 && fDisk5.size() != 0 && fDisk6.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk3.begin(); it1 != fDisk3.end(); ++it1)
            {
                TVector3 pv1; // PVD3
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk5.begin(); it2 != fDisk5.end(); ++it2)
                {
                    TVector3 pp1; // D3D5
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    TVector3 pv2; // PVD5
                    pv2.SetXYZ( (*it2).pos_x - simVx->at(0), (*it2).pos_y - simVy->at(0), (*it2).pos_z - simVz->at(0) );
                    Float_t pv2_eta = pv2.Eta();

                    for(std::vector<track>::iterator it3 = fDisk6.begin(); it3 != fDisk6.end(); ++it3 )
                    {
                        TVector3 pp2; // D5D6
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D3D6
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        TVector3 pv3; // PVD6
                        pv3.SetXYZ( (*it3).pos_x - simVx->at(0), (*it3).pos_y - simVy->at(0), (*it3).pos_z - simVz->at(0) );
                        Float_t pv3_eta = pv3.Eta();

                        Float_t dEtaD3556 = pp1_eta - pp2_eta;
                        Float_t dEtaD3536 = pp1_eta - pp3_eta;
                        Float_t dEtaD3656 = pp3_eta - pp2_eta;
                        Float_t dEtaPVD65 = pv3_eta - pv2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD3 - D3D5
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D3D5 - D5D6
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase3.push_back(dEtaD3536);
                        saveParaCase3.push_back(dEtaD3556);
                        saveParaCase3.push_back(dEtaD3656);
                        saveParaCase3.push_back(dEtaPVD65);
                        saveParaCase3.push_back(ddPhi);
                    }
                }
            }
        }

        //// -------- D456 --------- //// 
        if( fDisk4.size() != 0 && fDisk5.size() != 0 && fDisk6.size() != 0 )
        {
            for(std::vector<track>::iterator it1 = fDisk4.begin(); it1 != fDisk4.end(); ++it1)
            {
                TVector3 pv1; // PVD4
                pv1.SetXYZ( (*it1).pos_x - simVx->at(0), (*it1).pos_y - simVy->at(0), (*it1).pos_z - simVz->at(0) );
                Float_t pv1_phi = pv1.Phi();
                for(std::vector<track>::iterator it2 = fDisk5.begin(); it2 != fDisk5.end(); ++it2)
                {
                    TVector3 pp1; // D4D5
                    pp1.SetXYZ( (*it2).pos_x - (*it1).pos_x, (*it2).pos_y - (*it1).pos_y, (*it2).pos_z - (*it1).pos_z );
                    Float_t pp1_phi = pp1.Phi();
                    Float_t pp1_eta = pp1.Eta();

                    for(std::vector<track>::iterator it3 = fDisk6.begin(); it3 != fDisk6.end(); ++it3 )
                    {
                        TVector3 pp2; // D5D6
                        pp2.SetXYZ( (*it3).pos_x - (*it2).pos_x, (*it3).pos_y - (*it2).pos_y, (*it3).pos_z - (*it2).pos_z );
                        Float_t pp2_phi = pp2.Phi();
                        Float_t pp2_eta = pp2.Eta();

                        TVector3 pp3; // D4D6
                        pp3.SetXYZ( (*it3).pos_x - (*it1).pos_x, (*it3).pos_y - (*it1).pos_y, (*it3).pos_z - (*it1).pos_z );
                        Float_t pp3_eta = pp3.Eta();

                        Float_t dEtaD4556 = pp1_eta - pp2_eta;
                        Float_t dEtaD4546 = pp1_eta - pp3_eta;
                        Float_t dEtaD4656 = pp3_eta - pp2_eta;

                        Float_t deltaPhi1 = pv1_phi - pp1_phi; // dPhi = PVD4 - D4D5
                        if( deltaPhi1 >= TMath::Pi() ) deltaPhi1 -= 2.*TMath::Pi();
                        if( deltaPhi1 < -TMath::Pi() ) deltaPhi1 += 2.*TMath::Pi();

                        Float_t deltaPhi2 = pp1_phi - pp2_phi; // dPhi = D4D5 - D5D6
                        if( deltaPhi2 >= TMath::Pi() ) deltaPhi2 -= 2.*TMath::Pi();
                        if( deltaPhi2 < -TMath::Pi() ) deltaPhi2 += 2.*TMath::Pi();

                        Float_t ddPhi = fabs(deltaPhi1) - fabs(deltaPhi2);

                        saveParaCase4.push_back(dEtaD4546);
                        saveParaCase4.push_back(dEtaD4556);
                        saveParaCase4.push_back(dEtaD4656);
                        saveParaCase4.push_back(ddPhi);
                    }
                }
            }
        }
    }
}

#endif // #ifdef test_cxx
