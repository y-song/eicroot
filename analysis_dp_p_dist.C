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

void analysis_res_dist()
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
  const double interval = 0.2;
  double eta_bin[nbin + 1];
  for (int i = 0; i < nbin + 1; i++){
  	eta_bin[i] = eta_min + i * interval;
  }
  int count[nbin] = {0};
 
#if 1
  
  TTreeFormula fEtaTru("fEtaTru","MCTrack.GetMomentum().Eta()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);
  TObjArray proj;
  TObjArray fitArr;
  TH1D count_eta("count_eta", "", nbin, eta_min, eta_max);
  count_eta.Sumw2(); // If TH1::Sumw2() has been called before filling, the sum of squares is also stored.
  for (int i = 0; i < nbin; i++) {
  	  proj.Add(new TH1D(Form("%0.1f_eta_%0.1f", eta_bin[i], eta_bin[i+1]),"", 50, -0.03, 0.03));
  	  //fitArr.Add(new TF1(Form("fit%lu", i), "gaus(0)+gaus(3)", -0.03, 0.03));
 	  fitArr.Add(new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]^2)*exp(-0.5*((x-[1])/[4])**2)", -0.03, 0.03));
  }
  for (Long64_t j = 0; j < cbmsim->GetEntries(); j++) {
  //for (Long64_t j = 0; j < 20000; j++) {
  	  cbmsim->GetEntry(j);
  	  const double eta_tru = fEtaTru.EvalInstance();
	  const double x = fx.EvalInstance();
  	  const Int_t root_bin = count_eta.GetXaxis()->FindFixBin(eta_tru); // TAxis * axis is the axis object, which will be returned when calling the TH1::GetAxis() method. FindFixBin() finds bin number corresponding to abscissa x.
  	  if (root_bin >= 1 && root_bin < nbin + 1) {
  		  ((TH1D *)proj.At(root_bin - 1))->Fill(x);
	  }
	  for (int i = 0; i < nbin; i++){
	  	  if (eta_bin[i] < eta_tru && eta_tru < eta_bin[i+1]) {
	  	  	count[i] += 1;
  	  	  }
  	  }
  } 

  for (int i = 0; i < nbin; i++){
	cout << i << " - " << count[i] << endl;
  }
  TFile f("output/5_w_2gaus.root","recreate");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  gStyle->SetStatW(0.155);
  gStyle->SetStatH(0.255);
  //gStyle->SetTitleFontSize(.8);
  //gStyle->SetLabelSize(.5, "XY");
  TCanvas *canv = new TCanvas("canv", "canv", 1600, 1200);
  canv->Divide(4,4);
  
  for (int i = 17; i < 33; i++){
 	canv->cd(i-16);
 	dist = ((TH1D *)proj.At(i));
	fit = ((TF1 *)fitArr.At(i));
	fit->SetParameter(4, -0.015);
	fit->SetParameter(2, -0.015);
	if (i == 17 || i == 32){
		fit->SetParLimits(4, -0.05, 0.05);
		fit->SetParLimits(2, -0.05, 0.05);
	}
	else{
		fit->SetParLimits(4, -0.03, 0.03);
		fit->SetParLimits(2, -0.03, 0.03);
	}
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
