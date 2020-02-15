#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TTreeFormula.h>
#include <TStyle.h>
#include <TH2.h>
#include <TH1.h>
#include <TColor.h>
#include <TPaveStats.h>
#include <PndPidCandidate.h>

#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>

#define _NDF_MAX_ 1000

void analysis_2D()
{
  // Load basic libraries;
  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");  

  // Input simulated & reconstructed files;
  TFile *ff = new TFile("simulation.root");
  TTree *cbmsim = (TTree *)ff->Get("cbmsim"); 
  cbmsim->AddFriend("cbmsim", "reconstruction.root");

  // Figure out max and most popular ndf value;
  TClonesArray *rctrk = new TClonesArray("PndPidCandidate");
  cbmsim->SetBranchAddress("PidChargedCand",&rctrk);
  unsigned nEvt = cbmsim->GetEntries(), ndfMax = 0;
  unsigned arr[1000]; memset(arr, 0x00, sizeof(arr));
  for (unsigned evt = 0; evt<nEvt; evt++) {
    cbmsim->GetEntry(evt);
    
    if (rctrk->GetEntriesFast()) {
      PndPidCandidate *rctrack = (PndPidCandidate *)rctrk->At(0);
      
      int ndf = rctrack->GetDegreesOfFreedom();
      if (ndf > ndfMax) ndfMax = ndf;
      if (ndf < _NDF_MAX_) arr[ndf]++;
    }
  }
  unsigned ndfMostPopular = 0, ndfMostPopularStat = 0;
  for(unsigned iq=0; iq<_NDF_MAX_; iq++)
    if (arr[iq] > ndfMostPopularStat) {
      ndfMostPopular = iq;
      ndfMostPopularStat = arr[iq];
    }

  const char *expression = "(PidChargedCand.GetMomentum().Mag()-MCTrack.GetMomentum().Mag())/MCTrack.GetMomentum().Mag()";
  char cut[1024];
  // Allow to lose 3 degrees of freedom compared to the max (assume standard) case;
  sprintf(cut, "EicIdealGenTrack.fNDF>=%d&&EicIdealGenTrack.fChiSquareCCDF>.001", ndfMostPopular-2);
  sprintf(cut, "");

  double par[100]; memset(par, 0x00, sizeof(par));

  const double eta_min = 1.0;
  const double eta_max = 1.2;

#if 1
  
  TTreeFormula fEtaTru("fEtaTru", "MCTrack.GetMomentum().Eta()", cbmsim);
  TTreeFormula fHit("fHit", "FstMoCaPoint@.GetEntries()+BstMoCaPoint@.GetEntries()+VstMoCaPoint@.GetEntries()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);

  TH2* hist = new TH2D("hist", Form("%0.1f < eta < %0.1f", eta_min, eta_max), 6, 3, 9, 20, -0.03, 0.03);
  for (Long64_t j = 0; j < cbmsim->GetEntries(); j++) {
 	  cbmsim->GetEntry(j);
  	  const double eta_tru = fEtaTru.EvalInstance();
	  const int nhit = fHit.EvalInstance();
	  const double res = fx.EvalInstance();
	  if (eta_tru >= eta_min && eta_tru < eta_max) {
  		  hist->Fill(nhit, res);
	  }
  } 

  TFile f("output/default.root","recreate");

  gStyle->SetOptStat(10);
  gStyle->SetStatW(0.155);
  gStyle->SetStatH(0.255);
  
  hist->GetXaxis()->SetTitle("number of hits");
  hist->GetYaxis()->SetTitle("dp/p");
  hist->Draw("colz");

#endif
}
