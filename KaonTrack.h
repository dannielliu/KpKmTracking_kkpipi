//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep  8 21:36:37 2015 by ROOT version 5.34/28
// from TTree KaonTrack/ks N-Tuple example
// found on file: KK_30800.root
//////////////////////////////////////////////////////////

#ifndef KaonTrack_h
#define KaonTrack_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TrackingAlg;
class KaonTrack {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           rec;
   Int_t           evttag;
   Int_t           indexmc;
   Int_t           pdgid[100];   //[indexmc]
   Int_t           motheridx[100];   //[indexmc]
   Int_t           ISRtag;
   Int_t           Ktag;
   Int_t           ngch;
   Int_t           ncharg;
   Int_t           nneu;
   Int_t           pi1chrg;
   Double_t        pi1px;
   Double_t        pi1py;
   Double_t        pi1pz;
   Double_t        pi1e;
   Int_t           pi2chrg;
   Double_t        pi2px;
   Double_t        pi2py;
   Double_t        pi2pz;
   Double_t        pi2e;
   Int_t           ka1chrg;
   Double_t        ka1px;
   Double_t        ka1py;
   Double_t        ka1pz;
   Double_t        ka1e;
   Int_t           ka_chrg;
   Double_t        ka_px;
   Double_t        ka_py;
   Double_t        ka_pz;
   Double_t        ka_e;
   Double_t        costheta;
   Double_t        vtxchi2;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_rec;   //!
   TBranch        *b_evttag;   //!
   TBranch        *b_indexmc;   //!
   TBranch        *b_pdgid;   //!
   TBranch        *b_motheridx;   //!
   TBranch        *b_ISRtag;   //!
   TBranch        *b_Ktag;   //!
   TBranch        *b_ngch;   //!
   TBranch        *b_ncharg;   //!
   TBranch        *b_nneu;   //!
   TBranch        *b_pi1chrg;   //!
   TBranch        *b_pi1px;   //!
   TBranch        *b_pi1py;   //!
   TBranch        *b_pi1pz;   //!
   TBranch        *b_pi1e;   //!
   TBranch        *b_pi2chrg;   //!
   TBranch        *b_pi2px;   //!
   TBranch        *b_pi2py;   //!
   TBranch        *b_pi2pz;   //!
   TBranch        *b_pi2e;   //!
   TBranch        *b_ka1chrg;   //!
   TBranch        *b_ka1px;   //!
   TBranch        *b_ka1py;   //!
   TBranch        *b_ka1pz;   //!
   TBranch        *b_ka1e;   //!
   TBranch        *b_ka_chrg;   //!
   TBranch        *b_ka_px;   //!
   TBranch        *b_ka_py;   //!
   TBranch        *b_ka_pz;   //!
   TBranch        *b_ka_e;   //!
   TBranch        *b_costheta;   //!
   TBranch        *b_vtxchi2;   //!

   KaonTrack(TTree *tree=0);
   virtual ~KaonTrack();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TrackingAlg *trkalg);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef KaonTrack_cxx
KaonTrack::KaonTrack(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("KK_30800.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("KK_30800.root");
      }
      f->GetObject("KaonTrack",tree);

   }
   Init(tree);
}

KaonTrack::~KaonTrack()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t KaonTrack::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t KaonTrack::LoadTree(Long64_t entry)
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

void KaonTrack::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("rec", &rec, &b_rec);
   fChain->SetBranchAddress("evttag", &evttag, &b_evttag);
   fChain->SetBranchAddress("indexmc", &indexmc, &b_indexmc);
   fChain->SetBranchAddress("pdgid", pdgid, &b_pdgid);
   fChain->SetBranchAddress("motheridx", motheridx, &b_motheridx);
   fChain->SetBranchAddress("ISRtag", &ISRtag, &b_ISRtag);
   fChain->SetBranchAddress("Ktag", &Ktag, &b_Ktag);
   fChain->SetBranchAddress("ngch", &ngch, &b_ngch);
   fChain->SetBranchAddress("ncharg", &ncharg, &b_ncharg);
   fChain->SetBranchAddress("nneu", &nneu, &b_nneu);
   fChain->SetBranchAddress("pi1chrg", &pi1chrg, &b_pi1chrg);
   fChain->SetBranchAddress("pi1px", &pi1px, &b_pi1px);
   fChain->SetBranchAddress("pi1py", &pi1py, &b_pi1py);
   fChain->SetBranchAddress("pi1pz", &pi1pz, &b_pi1pz);
   fChain->SetBranchAddress("pi1e", &pi1e, &b_pi1e);
   fChain->SetBranchAddress("pi2chrg", &pi2chrg, &b_pi2chrg);
   fChain->SetBranchAddress("pi2px", &pi2px, &b_pi2px);
   fChain->SetBranchAddress("pi2py", &pi2py, &b_pi2py);
   fChain->SetBranchAddress("pi2pz", &pi2pz, &b_pi2pz);
   fChain->SetBranchAddress("pi2e", &pi2e, &b_pi2e);
   fChain->SetBranchAddress("ka1chrg", &ka1chrg, &b_ka1chrg);
   fChain->SetBranchAddress("ka1px", &ka1px, &b_ka1px);
   fChain->SetBranchAddress("ka1py", &ka1py, &b_ka1py);
   fChain->SetBranchAddress("ka1pz", &ka1pz, &b_ka1pz);
   fChain->SetBranchAddress("ka1e", &ka1e, &b_ka1e);
   fChain->SetBranchAddress("ka_chrg", &ka_chrg, &b_ka_chrg);
   fChain->SetBranchAddress("ka_px", &ka_px, &b_ka_px);
   fChain->SetBranchAddress("ka_py", &ka_py, &b_ka_py);
   fChain->SetBranchAddress("ka_pz", &ka_pz, &b_ka_pz);
   fChain->SetBranchAddress("ka_e", &ka_e, &b_ka_e);
   fChain->SetBranchAddress("costheta", &costheta, &b_costheta);
   fChain->SetBranchAddress("vtxchi2", &vtxchi2, &b_vtxchi2);
   Notify();
}

Bool_t KaonTrack::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void KaonTrack::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t KaonTrack::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef KaonTrack_cxx
