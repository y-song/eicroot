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

void analysis_2D_dist()
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

  const int nhit_min = 3;
  const int nhit_max = 9;
  const int nhit_bin = nhit_max - nhit_min;
  const double res_min = -0.03;
  const double res_max = 0.03;
  const int res_bin = (res_max - res_min) / 0.005;
  const double eta_min = -5.0;
  const double eta_max = 5.0;
  const int eta_bin = (eta_max - eta_min) / 0.2;
  double eta_edge[51];
  for (int i = 0; i < eta_bin + 1; i++){
  	eta_edge[i] = eta_min + i * 0.2;
  }

#if 1
  
  TTreeFormula fEtaTru("fEtaTru", "MCTrack.GetMomentum().Eta()", cbmsim);
  TTreeFormula fHit("fHit", "FstMoCaPoint@.GetEntries()+BstMoCaPoint@.GetEntries()+VstMoCaPoint@.GetEntries()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);
  
  TObjArray hist_arr;
  TObjArray norm_arr;
  TH1D count_eta("count_eta", "", eta_bin, eta_min, eta_max);
  count_eta.Sumw2();

  for (int i = 0; i < eta_bin; i++) {
	hist_arr.Add(new TH2D("hist", Form("%0.1f < eta < %0.1f", eta_edge[i], eta_edge[i+1]),nhit_bin,nhit_min,nhit_max,res_bin,res_min,res_max));
  	norm_arr.Add(new TH2D("norm", Form("%0.1f < eta < %0.1f", eta_edge[i], eta_edge[i+1]),nhit_bin,nhit_min,nhit_max,res_bin,res_min,res_max));
  }

  for (Long64_t j = 0; j < cbmsim->GetEntries(); j++) {
	
	cbmsim->GetEntry(j);
       	const double eta_tru = fEtaTru.EvalInstance();
       	const int nhit = fHit.EvalInstance();
       	const double res = fx.EvalInstance();
        const Int_t root_bin = count_eta.GetXaxis()->FindFixBin(eta_tru);
  	
	if (root_bin >= 1 && root_bin < eta_bin + 1) {
  		((TH2D *)hist_arr.At(root_bin - 1))->Fill(nhit, res);
        }
  }

  TFile f("output/default.root","recreate");
  gStyle->SetOptStat(0);
  TCanvas *canv = new TCanvas("canv", "canv", 1600, 1200);
  canv->Divide(4,4);
  
  for (int i = 17; i < 33; i++){
 	canv->cd(i-16);
 	hist = (TH2D *)hist_arr.At(i);
	norm = (TH2D *)norm_arr.At(i);
	int nhit_arr[nhit_bin] = {0};
  	for (int k = 0; k < nhit_bin; k++){
        	nhit_arr[k] = hist->Integral(k+1, k+1);
	}
  	 
	for (int l = 0; l < nhit_bin; l++){
          	for (int k = 0; k < res_bin; k++){
                  	if (nhit_arr[l] != 0){
				double content = hist->GetBinContent(l+1, k+1) / nhit_arr[l];
                  	 	norm->SetBinContent(l+1, k+1, content);
          		}
		}
 	}
	norm->GetXaxis()->SetTitle("number of hits");
 	norm->GetYaxis()->SetTitle("dp/p");
  	norm->Draw("colz");
  }

  canv->Write();

#endif
}
