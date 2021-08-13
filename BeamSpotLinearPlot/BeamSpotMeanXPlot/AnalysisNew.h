//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 16 08:18:52 2020 by ROOT version 6.10/04
// from TTree ntuple1/Data for file
// found on file: basicOriginal_Pre.root
//////////////////////////////////////////////////////////

#ifndef AnalysisNew_h
#define AnalysisNew_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class AnalysisNew {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         nFill;
   Float_t         nTime;
   Float_t         meanX;
   Float_t         sigmaX;
   Float_t         meanY;
   Float_t         sigmaY;
   Float_t         meanErrX;
   Float_t         meanErrY;

   // List of branches
   TBranch        *b_nFill;   //!
   TBranch        *b_nTime;   //!
   TBranch        *b_meanX;   //!
   TBranch        *b_sigmaX;   //!
   TBranch        *b_meanY;   //!
   TBranch        *b_sigmaY;   //!
   TBranch        *b_meanErrX;   //!
   TBranch        *b_meanErrY;   //!

   AnalysisNew(TTree *tree=0);
   virtual ~AnalysisNew();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AnalysisNew_cxx
AnalysisNew::AnalysisNew(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("basicOriginal_FinAll.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("basicOriginal_FinAll.root");
      }
      f->GetObject("ntuple1",tree);

   }
   Init(tree);
}

AnalysisNew::~AnalysisNew()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisNew::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalysisNew::LoadTree(Long64_t entry)
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

void AnalysisNew::Init(TTree *tree)
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

   fChain->SetBranchAddress("nFill", &nFill, &b_nFill);
   fChain->SetBranchAddress("nTime", &nTime, &b_nTime);
   fChain->SetBranchAddress("meanX", &meanX, &b_meanX);
   fChain->SetBranchAddress("sigmaX", &sigmaX, &b_sigmaX);
   fChain->SetBranchAddress("meanY", &meanY, &b_meanY);
   fChain->SetBranchAddress("sigmaY", &sigmaY, &b_sigmaY);
   fChain->SetBranchAddress("meanErrX", &meanErrX, &b_meanErrX);
   fChain->SetBranchAddress("meanErrY", &meanErrY, &b_meanErrY);
   Notify();
}

Bool_t AnalysisNew::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalysisNew::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalysisNew::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AnalysisNew_cxx
