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

void analysis_2D_norm()
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
  const int nhit_min = 3;
  const int nhit_max = 9;
  const int nhit_bin = nhit_max - nhit_min;
  const double res_min = -0.03;
  const double res_max = 0.03;
  const int res_bin = (res_max - res_min) / 0.005;
  int nhit_arr[nhit_bin] = {0};

#if 1
  
  TTreeFormula fEtaTru("fEtaTru", "MCTrack.GetMomentum().Eta()", cbmsim);
  TTreeFormula fHit("fHit", "FstMoCaPoint@.GetEntries()+BstMoCaPoint@.GetEntries()+VstMoCaPoint@.GetEntries()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);

  TH2* hist = new TH2D("hist", Form("%0.1f < eta < %0.1f", eta_min, eta_max), nhit_bin, nhit_min, nhit_max, res_bin, res_min, res_max);
  TH2* norm = new TH2D("norm", Form("%0.1f < eta < %0.1f", eta_min, eta_max), nhit_bin, nhit_min, nhit_max, res_bin, res_min, res_max);
  for (Long64_t j = 0; j < cbmsim->GetEntries(); j++) {
 	  cbmsim->GetEntry(j);
  	  const double eta_tru = fEtaTru.EvalInstance();
	  const int nhit = fHit.EvalInstance();
	  const double res = fx.EvalInstance();
	  if (eta_tru >= eta_min && eta_tru < eta_max) {
  		  hist->Fill(nhit, res);
	  }
  }

  for (int i = 0; i < nhit_bin; i++){
	  nhit_arr[i] = hist->Integral(i+1, i+1);
  }

  for (int i = 0; i < nhit_bin; i++){
  	  for (int k = 0; k < res_bin; k++){
		  double content = hist->GetBinContent(i+1, k+1) / nhit_arr[i];
  	  	  norm->SetBinContent(i+1, k+1, content);
  	  }
  }

  TFile f("output/default.root","recreate");

  gStyle->SetOptStat(10);
  gStyle->SetStatW(0.155);
  gStyle->SetStatH(0.255);
  
  norm->GetXaxis()->SetTitle("number of hits");
  norm->GetYaxis()->SetTitle("dp/p");
  norm->Draw("colz");

#endif
}
