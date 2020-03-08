#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TTreeFormula.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TPaveStats.h>
#include <PndPidCandidate.h>

#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>

#define _NDF_MAX_ 1000

void analysis_dp_p_dist_eta_cone()
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

  const int nhit_min = 0;
  const int nhit_max = 9;
  const int nhit_bin = nhit_max - nhit_min;
  const double res_min = -0.03;
  const double res_max = 0.03;
  const int res_bin = (res_max - res_min) / 0.005;
  const double eta_min = 1.0;
  const double eta_max = 1.2;

#if 1
  
  TTreeFormula fEtaTru("fEtaTru", "MCTrack.GetMomentum().Eta()", cbmsim);
  TTreeFormula fHit("fHit", "FstMoCaPoint@.GetEntries()+BstMoCaPoint@.GetEntries()+VstMoCaPoint@.GetEntries()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);

  TObjArray proj;
  TObjArray fitArr;
  TH1D count_nhit("count_nhit", "", nhit_bin, nhit_min, nhit_max);
  count_nhit.Sumw2(); // If TH1::Sumw2() has been called before filling, the sum of squares is also stored.
  for (int i = 0; i < nhit_bin; i++) {
  	  proj.Add(new TH1D(Form("nhit=%d", i),"", res_bin,res_min,res_max));
	  fitArr.Add(new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]^2)*exp(-0.5*((x-[1])/[4])**2)", -0.03, 0.03));
  }
 
  for (Long64_t j = 0; j < cbmsim->GetEntries(); j++) {
    	  cbmsim->GetEntry(j);
  	  const double eta_tru = fEtaTru.EvalInstance();
	  const double x = fx.EvalInstance();
	  const int nhit = fHit.EvalInstance();
	  if (eta_tru > eta_min && eta_tru < eta_max) {
	  	((TH1D *)proj.At(nhit))->Fill(x);
  	  }
  } 

  TFile f("output/default.root","recreate");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.155);
  gStyle->SetStatH(0.255);
  TCanvas *canv = new TCanvas("canv", "canv", 1600, 1200);
  canv->Divide(3,2);
  
  for (int i = 3; i < nhit_bin; i++){
 	canv->cd(i-2);
 	dist = ((TH1D *)proj.At(i));
	fit = ((TF1 *)fitArr.At(i));
	fit->SetParameter(4, -0.015);
	fit->SetParameter(2, -0.015);
	fit->SetParLimits(4, -0.03, 0.03);
	fit->SetParLimits(2, -0.03, 0.03);

 	dist->Fit(fit, "rll");
 	dist->GetXaxis()->SetTitle("dp/p");
  	dist->GetYaxis()->SetTitle("counts");
  	dist->SetMarkerStyle(20);
        dist->SetLineColor(kBlack);
	dist->Write();
	fit->Write();
        dist->Draw();
  }
  canv->Write();

#endif
}
