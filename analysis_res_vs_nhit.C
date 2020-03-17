#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TTreeFormula.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TObjArray.h>

#include <cmath>
#include <math.h>

#include <PndPidCandidate.h>

#define _NDF_MAX_ 1000

void analysis_res_vs_nhit()
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
  const double eta_min = 1.0;
  const double eta_max = 1.2;

#if 1
  gStyle->SetOptStat(0);
  
  TTreeFormula fx("fx", expression, cbmsim);
  TTreeFormula fHit("fHit", "FstMoCaPoint@.GetEntries()+BstMoCaPoint@.GetEntries()+VstMoCaPoint@.GetEntries()", cbmsim);
  TTreeFormula fRec("fRec","PidChargedCand@.GetEntries()",cbmsim);
 
  TObjArray proj;
  TObjArray fitArr;
  for (int i = 0; i < nhit_bin; i++) {
	  if (i == 4 || i == 5 || i == 7 || i == 8){
	  	proj.Add(new TH1D(Form("proj%lu", i), "", 50, -0.05, 0.05));
		fitArr.Add(new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]^2)*exp(-0.5*((x-[1])/[4])**2)", -0.05, 0.05));
  	  }
	  else if (i == 6) {
	  	proj.Add(new TH1D(Form("proj%lu", i), "", 50, -0.03, 0.03));
		fitArr.Add(new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]^2)*exp(-0.5*((x-[1])/[4])**2)", -0.03, 0.03));
  	  }
	  else {
	  	proj.Add(new TH1D(Form("proj%lu", i), "", 50, -0.5, 0.5));
		fitArr.Add(new TF1(Form("fit%lu", i), "[0]/sqrt(2*TMath::Pi()*[2]^2)*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]^2)*exp(-0.5*((x-[1])/[4])**2)", -0.5, 0.5));
  	  }
  }

  // Attempts on figuring out if x (the dp/p variable) is saved if the track is not reconstructed

  const int entry = cbmsim->GetEntries();
  cout << entry << endl;

  double x_arr[100000] = {0};

  for (Long64_t j = 0; j < cbmsim->GetEntries(); j++) {
  	  cbmsim->GetEntry(j);
  	  const double nhit = fHit.EvalInstance();
	  const double x = fx.EvalInstance();
          const int rec = fRec.EvalInstance();
	  x_arr[j] = x;
	  if (rec != 1) {
	  	 cout << x << endl;
	  }
	  ((TH1D *)proj.At(nhit))->Fill(x);
  }

  int count = 0;
  for (Long64_t k = 0; k < cbmsim->GetEntries(); k++) {
	  if (x_arr[k] != 0) {
	  	count += 1;
	  }
  }

  cout << count << endl;
 
  double y_arr[nhit_bin] = {0}; 
  double yerr_arr[nhit_bin] = {0};

  for (int i = 0; i < nhit_bin; i++) {
  	cout << "++++++++++++++++++++++++" << i << endl;
	f = (TF1 *)fitArr.At(i);
	f->SetParameter(4, -0.015);
	f->SetParameter(2, -0.015);

	if (i == 4 || i == 5 || i == 7 || i == 8){
		f->SetParLimits(4, -0.05, 0.05);
		f->SetParLimits(2, -0.05, 0.05);
	}
	else if (i == 6){
		f->SetParLimits(4, -0.03, 0.03);
		f->SetParLimits(2, -0.03, 0.03);
	}
	else {
		f->SetParLimits(4, -0.5, 0.5);
		f->SetParLimits(2, -0.5, 0.5);
	}

	((TH1D *)proj.At(i))->Fit(f, "rll");

  	p0 = f->GetParameter(0);
	p2 = fabs(f->GetParameter(2));
	p3 = f->GetParameter(3);
	p4 = fabs(f->GetParameter(4));

	if (p0 > p3){
		y_arr[i] = p2;
		yerr_arr[i] = f->GetParError(2);
	}
	else{
		y_arr[i] = p4;
		yerr_arr[i] = f->GetParError(4);
	}
  }

  for (int i = 0; i < nhit_bin; i++){
	cout << y_arr[i] << endl;
  }

  cout << "++++++++++++++++++++++++++++++++" << endl;

  for (int i = 0; i < nhit_bin; i++){
  	cout << yerr_arr[i] << endl;
  }

#endif
}
