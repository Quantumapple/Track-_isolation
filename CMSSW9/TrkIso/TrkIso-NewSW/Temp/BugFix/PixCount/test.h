//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 17 20:40:59 2015 by ROOT version 5.34/19
// from TTree t/t
// found on file: SingleEle_ntuple_1.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
#include <iostream>
#include <fstream>
#include <string>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include "./roi_v2.h"
#include "./withEM_v2.h"
#include "./withoutEM_v2.h"

using namespace std;
class test {
private:

map<TString, TH1*> maphist;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nVtx;
   Int_t           nMeanPU;
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
   Float_t         maxBrempos;
   Float_t         lastSimtkpt;
   Float_t         initialSimtkpt;
   Int_t           bremflag;
   vector<float>   *Brempos_radius;
   vector<float>   *Brempos_x;
   vector<float>   *Brempos_y;
   vector<float>   *Brempos_z;
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
   vector<float>   *egCrysE;
   vector<float>   *egCrysEt;
   vector<float>   *egCrysEta;
   vector<float>   *egCrysPhi;
   vector<float>   *egCrysGx;
   vector<float>   *egCrysGy;
   vector<float>   *egCrysGz;
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
   vector<float>   *cl3d_egid;
   vector<float>   *cl3d_energy;
   vector<float>   *cl3d_eta;
   vector<float>   *cl3d_phi;
   vector<int>     *cl3d_nclu;
   vector<float>   *cl3d_x;
   vector<float>   *cl3d_y;
   vector<int>     *cl3d_z;
   vector<int>     *cl3d_coreshowerlength;
   vector<float>   *cl3d_srrtot;
   vector<int>     *cl3d_maxlayer;
   vector<int>     *cl3d_firstlayer;

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

   UInt_t          ntkEG;
   vector<double>  *tkEGEt;
   vector<double>  *tkEGEta;
   vector<double>  *tkEGPhi;
   vector<int>     *tkEGBx;
   vector<double>  *tkEGTrkIso;
   vector<double>  *tkEGzVtx;
   vector<double>  *tkEGHwQual;
   vector<double>  *tkEGEGRefPt;
   vector<double>  *tkEGEGRefEta;
   vector<double>  *tkEGEGRefPhi;
   UInt_t          ntkEGLoose;
   vector<double>  *tkEGLooseEt;
   vector<double>  *tkEGLooseEta;
   vector<double>  *tkEGLoosePhi;
   vector<int>     *tkEGLooseBx;
   vector<double>  *tkEGLooseTrkIso;
   vector<double>  *tkEGLoosezVtx;
   vector<double>  *tkEGLooseHwQual;
   vector<double>  *tkEGLooseEGRefPt;
   vector<double>  *tkEGLooseEGRefEta;
   vector<double>  *tkEGLooseEGRefPhi;

   Int_t           tower_n;
   vector<float>   *tower_pt;
   vector<float>   *tower_energy;
   vector<float>   *tower_eta;
   vector<float>   *tower_phi;
   vector<float>   *tower_etEm;
   vector<float>   *tower_etHad;
   vector<float>   *tower_gx;
   vector<float>   *tower_gy;
   vector<float>   *tower_gz;
   vector<int>     *tower_iEta;
   vector<int>     *tower_iPhi;


   // List of branches
   TBranch        *b_nVtx;   //!
   TBranch        *b_nMeanPU;   //!
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
   TBranch        *b_maxBrempos;   //!
   TBranch        *b_lastSimtkpt;   //!
   TBranch        *b_initialSimtkpt;   //!
   TBranch        *b_bremflag;   //!
   TBranch        *b_Brempos_radius;   //!
   TBranch        *b_Brempos_x;   //!
   TBranch        *b_Brempos_y;   //!
   TBranch        *b_Brempos_z;   //!
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
   TBranch        *b_egCrysE;   //!
   TBranch        *b_egCrysEt;   //!
   TBranch        *b_egCrysEta;   //!
   TBranch        *b_egCrysPhi;   //!
   TBranch        *b_egCrysGx;   //!
   TBranch        *b_egCrysGy;   //!
   TBranch        *b_egCrysGz;   //!
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
   TBranch        *b_trackmatchingdR;
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
   TBranch        *b_hgcal_trackmatchingdR;
   TBranch        *b_cl3d_n;   //!
   TBranch        *b_cl3d_pt;   //!
   TBranch        *b_cl3d_egid;   //!
   TBranch        *b_cl3d_energy;   //!
   TBranch        *b_cl3d_eta;   //!
   TBranch        *b_cl3d_phi;   //!
   TBranch        *b_cl3d_nclu;   //!
   TBranch        *b_cl3d_x;   //!
   TBranch        *b_cl3d_y;   //!
   TBranch        *b_cl3d_z;   //!
   TBranch        *b_cl3d_coreshowerlength;
   TBranch        *b_cl3d_srrtot;
   TBranch        *b_cl3d_maxlayer;
   TBranch        *b_cl3d_firstlayer;
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
   TBranch        *b_tower_n;   //!
   TBranch        *b_tower_pt;   //!
   TBranch        *b_tower_energy;   //!
   TBranch        *b_tower_eta;   //!
   TBranch        *b_tower_phi;   //!
   TBranch        *b_tower_etEm;   //!
   TBranch        *b_tower_etHad;   //!
   TBranch        *b_tower_gx;   //!
   TBranch        *b_tower_gy;   //!
   TBranch        *b_tower_gz;   //!
   TBranch        *b_tower_iEta;   //!
   TBranch        *b_tower_iPhi;   //!
   TBranch        *b_ntkEG;   //!
   TBranch        *b_tkEGEt;   //!
   TBranch        *b_tkEGEta;   //!
   TBranch        *b_tkEGPhi;   //!
   TBranch        *b_tkEGBx;   //!
   TBranch        *b_tkEGTrkIso;   //!
   TBranch        *b_tkEGzVtx;   //!
   TBranch        *b_tkEGHwQual;   //!
   TBranch        *b_tkEGEGRefPt;   //!
   TBranch        *b_tkEGEGRefEta;   //!
   TBranch        *b_tkEGEGRefPhi;   //!
   TBranch        *b_ntkEGLoose;   //!
   TBranch        *b_tkEGLooseEt;   //!
   TBranch        *b_tkEGLooseEta;   //!
   TBranch        *b_tkEGLoosePhi;   //!
   TBranch        *b_tkEGLooseBx;   //!
   TBranch        *b_tkEGLooseTrkIso;   //!
   TBranch        *b_tkEGLoosezVtx;   //!
   TBranch        *b_tkEGLooseHwQual;   //!
   TBranch        *b_tkEGLooseEGRefPt;   //!
   TBranch        *b_tkEGLooseEGRefEta;   //!
   TBranch        *b_tkEGLooseEGRefPhi;   //!


   int Ele, Pos;
   int skip;
   float EgN;
   int eta_region;

   TVector3 emvector;
   float EgEt;
   float EgEta;
   float EgPhi;
   float EgGx;
   float EgGy;
   float EgGz;

   void StorePixelHit( int region);
   
   std::vector<TVector3> first_hits;
   std::vector<TVector3> second_hits;
   std::vector<TVector3> third_hits;
   std::vector<TVector3> fourth_hits;
  
   inline float deltaPhi(float phi1, float phi2) { 
     float result = phi1 - phi2;
     //while (result > float(M_PI)) result -= float(2*M_PI);
     //while (result <= -float(M_PI)) result += float(2*M_PI);
     if( result >= TMath::Pi() ) result -= 2.*TMath::Pi();
     if( result < -TMath::Pi() ) result += 2.*TMath::Pi();
     return result;
   }

   inline float radius(float x1, float x2) {
       float result = sqrt(pow(x1,2)+pow(x2,2));
       return result; 
   }

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TFile *outfile;
   TTree* pixtrk_tree;

   float nt_genPhi;
   float nt_genEta;
   float nt_genPt;
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

   if (tree == 0) {
      TChain * chain = new TChain("l1PiXTRKTree/L1PiXTRKTree","");
      string line;
      ifstream myfile("SE_PU200_1.txt"); // SE_PU200_1.txt will be replaced by the name of txt file that contains the location of input files
      if (myfile.is_open())
        {
          while ( getline (myfile,line) )
          {
            if( line.length() != 0 ){
              cout << line << '\n';
              char file_path[300];
              strcpy(file_path, line.c_str());
              chain->Add(file_path);
            }
          }
          myfile.close();
        }
      tree = chain;
   }
   Init(tree);
 
   //outfile = new TFile("../output_tmp/Tree_SE_PU200_1.root","recreate");
   outfile = new TFile("test.root","RECREATE");
   pixtrk_tree = new TTree("t","t");

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
   Brempos_x = 0;
   Brempos_y = 0;
   Brempos_z = 0;
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
   cl3d_egid = 0;
   cl3d_energy = 0;
   cl3d_eta = 0;
   cl3d_phi = 0;
   cl3d_nclu = 0;
   cl3d_x = 0;
   cl3d_y = 0;
   cl3d_z = 0;
   cl3d_coreshowerlength = 0;
   cl3d_srrtot = 0;
   cl3d_maxlayer = 0;
   cl3d_firstlayer = 0;
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
   tower_pt = 0;
   tower_energy = 0;
   tower_eta = 0;
   tower_phi = 0;
   tower_etEm = 0;
   tower_etHad = 0;
   tower_gx = 0;
   tower_gy = 0;
   tower_gz = 0;
   tower_iEta = 0;
   tower_iPhi = 0;
   tkEGEt = 0;
   tkEGEta = 0;
   tkEGPhi = 0;
   tkEGBx = 0;
   tkEGTrkIso = 0;
   tkEGzVtx = 0;
   tkEGHwQual = 0;
   tkEGEGRefPt = 0;
   tkEGEGRefEta = 0;
   tkEGEGRefPhi = 0;
   tkEGLooseEt = 0;
   tkEGLooseEta = 0;
   tkEGLoosePhi = 0;
   tkEGLooseBx = 0;
   tkEGLooseTrkIso = 0;
   tkEGLoosezVtx = 0;
   tkEGLooseHwQual = 0;
   tkEGLooseEGRefPt = 0;
   tkEGLooseEGRefEta = 0;
   tkEGLooseEGRefPhi = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nMeanPU", &nMeanPU, &b_nMeanPU);
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
   fChain->SetBranchAddress("maxBrempos", &maxBrempos, &b_maxBrempos);
   fChain->SetBranchAddress("lastSimtkpt", &lastSimtkpt, &b_lastSimtkpt);
   fChain->SetBranchAddress("initialSimtkpt", &initialSimtkpt, &b_initialSimtkpt);
   fChain->SetBranchAddress("bremflag", &bremflag, &b_bremflag);
   fChain->SetBranchAddress("Brempos_radius", &Brempos_radius, &b_Brempos_radius);
   fChain->SetBranchAddress("Brempos_x", &Brempos_x, &b_Brempos_x);
   fChain->SetBranchAddress("Brempos_y", &Brempos_y, &b_Brempos_y);
   fChain->SetBranchAddress("Brempos_z", &Brempos_z, &b_Brempos_z);
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
   fChain->SetBranchAddress("egCrysE", &egCrysE, &b_egCrysE);
   fChain->SetBranchAddress("egCrysEt", &egCrysEt, &b_egCrysEt);
   fChain->SetBranchAddress("egCrysEta", &egCrysEta, &b_egCrysEta);
   fChain->SetBranchAddress("egCrysPhi", &egCrysPhi, &b_egCrysPhi);
   fChain->SetBranchAddress("egCrysGx", &egCrysGx, &b_egCrysGx);
   fChain->SetBranchAddress("egCrysGy", &egCrysGy, &b_egCrysGy);
   fChain->SetBranchAddress("egCrysGz", &egCrysGz, &b_egCrysGz);
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
   fChain->SetBranchAddress("cl3d_egid", &cl3d_egid, &b_cl3d_egid);
   fChain->SetBranchAddress("cl3d_energy", &cl3d_energy, &b_cl3d_energy);
   fChain->SetBranchAddress("cl3d_eta", &cl3d_eta, &b_cl3d_eta);
   fChain->SetBranchAddress("cl3d_phi", &cl3d_phi, &b_cl3d_phi);
   fChain->SetBranchAddress("cl3d_nclu", &cl3d_nclu, &b_cl3d_nclu);
   fChain->SetBranchAddress("cl3d_x", &cl3d_x, &b_cl3d_x);
   fChain->SetBranchAddress("cl3d_y", &cl3d_y, &b_cl3d_y);
   fChain->SetBranchAddress("cl3d_z", &cl3d_z, &b_cl3d_z);
   fChain->SetBranchAddress("cl3d_coreshowerlength", &cl3d_coreshowerlength, &b_cl3d_coreshowerlength);
   fChain->SetBranchAddress("cl3d_srrtot", &cl3d_srrtot, &b_cl3d_srrtot);
   fChain->SetBranchAddress("cl3d_maxlayer", &cl3d_maxlayer, &b_cl3d_maxlayer);
   fChain->SetBranchAddress("cl3d_firstlayer", &cl3d_firstlayer, &b_cl3d_firstlayer);
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
   fChain->SetBranchAddress("tower_n", &tower_n, &b_tower_n);
   fChain->SetBranchAddress("tower_pt", &tower_pt, &b_tower_pt);
   fChain->SetBranchAddress("tower_energy", &tower_energy, &b_tower_energy);
   fChain->SetBranchAddress("tower_eta", &tower_eta, &b_tower_eta);
   fChain->SetBranchAddress("tower_phi", &tower_phi, &b_tower_phi);
   fChain->SetBranchAddress("tower_etEm", &tower_etEm, &b_tower_etEm);
   fChain->SetBranchAddress("tower_etHad", &tower_etHad, &b_tower_etHad);
   fChain->SetBranchAddress("tower_gx", &tower_gx, &b_tower_gx);
   fChain->SetBranchAddress("tower_gy", &tower_gy, &b_tower_gy);
   fChain->SetBranchAddress("tower_gz", &tower_gz, &b_tower_gz);
   fChain->SetBranchAddress("tower_iEta", &tower_iEta, &b_tower_iEta);
   fChain->SetBranchAddress("tower_iPhi", &tower_iPhi, &b_tower_iPhi);
   fChain->SetBranchAddress("ntkEG", &ntkEG, &b_ntkEG);
   fChain->SetBranchAddress("tkEGEt", &tkEGEt, &b_tkEGEt);
   fChain->SetBranchAddress("tkEGEta", &tkEGEta, &b_tkEGEta);
   fChain->SetBranchAddress("tkEGPhi", &tkEGPhi, &b_tkEGPhi);
   fChain->SetBranchAddress("tkEGBx", &tkEGBx, &b_tkEGBx);
   fChain->SetBranchAddress("tkEGTrkIso", &tkEGTrkIso, &b_tkEGTrkIso);
   fChain->SetBranchAddress("tkEGzVtx", &tkEGzVtx, &b_tkEGzVtx);
   fChain->SetBranchAddress("tkEGHwQual", &tkEGHwQual, &b_tkEGHwQual);
   fChain->SetBranchAddress("tkEGEGRefPt", &tkEGEGRefPt, &b_tkEGEGRefPt);
   fChain->SetBranchAddress("tkEGEGRefEta", &tkEGEGRefEta, &b_tkEGEGRefEta);
   fChain->SetBranchAddress("tkEGEGRefPhi", &tkEGEGRefPhi, &b_tkEGEGRefPhi);
   fChain->SetBranchAddress("ntkEGLoose", &ntkEGLoose, &b_ntkEGLoose);
   fChain->SetBranchAddress("tkEGLooseEt", &tkEGLooseEt, &b_tkEGLooseEt);
   fChain->SetBranchAddress("tkEGLooseEta", &tkEGLooseEta, &b_tkEGLooseEta);
   fChain->SetBranchAddress("tkEGLoosePhi", &tkEGLoosePhi, &b_tkEGLoosePhi);
   fChain->SetBranchAddress("tkEGLooseBx", &tkEGLooseBx, &b_tkEGLooseBx);
   fChain->SetBranchAddress("tkEGLooseTrkIso", &tkEGLooseTrkIso, &b_tkEGLooseTrkIso);
   fChain->SetBranchAddress("tkEGLoosezVtx", &tkEGLoosezVtx, &b_tkEGLoosezVtx);
   fChain->SetBranchAddress("tkEGLooseHwQual", &tkEGLooseHwQual, &b_tkEGLooseHwQual);
   fChain->SetBranchAddress("tkEGLooseEGRefPt", &tkEGLooseEGRefPt, &b_tkEGLooseEGRefPt);
   fChain->SetBranchAddress("tkEGLooseEGRefEta", &tkEGLooseEGRefEta, &b_tkEGLooseEGRefEta);
   fChain->SetBranchAddress("tkEGLooseEGRefPhi", &tkEGLooseEGRefPhi, &b_tkEGLooseEGRefPhi);

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

void test::StorePixelHit(int region){

    for(int a=0; a<bRecHitN; a++){

        if( region < 5 ){
            if( bRecHitLayer->at(a) == 1 ){ // First layer
                first_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
        }
        if( region == 1 || region == 2 || region == 3 ){
            if( bRecHitLayer->at(a) == 2 ){ // Second layer
                second_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
        }
        if( region == 1 || region == 2 ){
            if( bRecHitLayer->at(a) == 3 ){ // Third layer 
                third_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
        }
        if( region == 1 ){
            if( bRecHitLayer->at(a) == 4 ){ // Fourth layer
                fourth_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
            }
        }
    }

    for(int a=0; a<fRecHitN; a++){

        if( region == 5 ){
            
            if( fRecHitDisk->at(a) == 1 ){ // first disk
                first_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            
            if( fRecHitDisk->at(a) == 2 ){ // second disk
                second_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            
            if( fRecHitDisk->at(a) == 3 ){ // third disk
                third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }

            if( fRecHitDisk->at(a) == 4 ){ // fourth disk
                fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
        }

        if( region == 6 ){
            
            if( fRecHitDisk->at(a) == 2 ){ // second disk
                first_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }

            if( fRecHitDisk->at(a) == 3 ){ // third disk
                second_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            
            if( fRecHitDisk->at(a) == 4 ){ // fourth disk
                third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }

            if( fRecHitDisk->at(a) == 5 ){ // fifth disk
                fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }

        }

        if( region == 4 ){
            
            if( fRecHitDisk->at(a) == 1 ){ // first disk
                second_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            
            if( fRecHitDisk->at(a) == 2 ){ // second disk
                third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            
            if( fRecHitDisk->at(a) == 3 ){ // third disk
                fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
        }

        if( region == 3 ){
            if( fRecHitDisk->at(a) == 1 ){ // first disk
                third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
            
            if( fRecHitDisk->at(a) == 2 ){ // second disk
                fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
        }

        if( region == 2 ){
            if( fRecHitDisk->at(a) == 1 ){ // first disk
                fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
            }
        }

    }
}

#endif // #ifdef test_cxx

