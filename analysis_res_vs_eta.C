#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TTreeFormula.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TObjArray.h>

#include <PndPidCandidate.h>

#define _NDF_MAX_ 1000

void analysis_res_vs_eta()
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

  const int nbin = 50;
  const double eta_min = -5.0;
  const double eta_max = 5.0;

#if 1
  gStyle->SetOptStat(0);
  
  TTreeFormula fEtaTru("fEtaTru","MCTrack.GetMomentum().Eta()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);
  TObjArray proj;
  TH1D count_eta("count_eta", "", nbin, eta_min, eta_max);
  count_eta.Sumw2(); // If TH1::Sumw2() has been called before filling, the sum of squares is also stored.
  for (int i = 0; i < nbin; i++) {
  	  proj.Add(new TH1D(Form("proj%lu", i), "", 200, -0.2, 0.2));
  }
  // for (Long64_t j = 0; j < cbmsim->GetEntries(); j++) {
  for (Long64_t j = 0; j < 200000; j++) {
  	  cbmsim->GetEntry(j);
  	  const double eta_tru = fEtaTru.EvalInstance();
	  const double x = fx.EvalInstance();
  	  const Int_t root_bin = count_eta.GetXaxis()->FindFixBin(eta_tru); // TAxis * axis is the axis object, which will be returned when calling the TH1::GetAxis() method. FindFixBin() finds bin number corresponding to abscissa x.
  	  if (root_bin >= 1 && root_bin < nbin + 1) {
  		  ((TH1D *)proj.At(root_bin - 1))->Fill(x);
  	  }
  }
 
  double y_arr[nbin] = {0}; 
  double yerr_arr[nbin] = {0};

  TObjArray fit;
  TH1D *s_x = new TH1D("s_x", "", nbin, eta_min, eta_max);
  s_x->Sumw2();
  for (int i = 0; i < nbin; i++) {
	if (proj.At(i) != NULL && ((TH1D *)proj.At(i))->GetEntries() >= 20) {
	  	TF1 *f = new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2)", -0.05, 0.05);
		f->SetParameter(2, 0.01);
		fit.Add(f);
  		((TH1D *)proj.At(i))->Fit((TF1 *)fit.At(i), "rll");
		y_arr[i] = ((TF1 *)fit.At(i))->GetParameter(2);
		yerr_arr[i] = ((TF1 *)fit.At(i))->GetParError(2);
	  	s_x->SetBinContent(i + 1, ((TF1 *)fit.At(i))->GetParameter(2));
	  	s_x->SetBinError(i + 1, ((TF1 *)fit.At(i))->GetParError(2));
	  }
  }
 
  s_x->SetMarkerStyle(20);
  s_x->SetMarkerColor(kBlack);
  s_x->SetLineColor(kBlack);
  s_x->Draw("e1x0");
 
  for (int i = 0; i < nbin; i++){
	cout << y_arr[i] << endl;
  }

  cout << "++++++++++++++++++++++++++++++++" << endl;

  for (int i = 0; i < nbin; i++){
  	cout << yerr_arr[i] << endl;
  }

#endif
}
