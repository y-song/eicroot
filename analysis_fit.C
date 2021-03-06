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

void analysis_fit()
{
  // Load basic libraries;
  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");  

  // Input simulated & reconstructed files;
  TFile *ff = new TFile("reconstruction.root");
  TTree *cbmsim = (TTree *)ff->Get("cbmsim"); 
  cbmsim->AddFriend("cbmsim", "simulation.root");

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

  char cut[1024];
  // Allow to lose 3 degrees of freedom compared to the max (assume standard) case;
  sprintf(cut, "EicIdealGenTrack.fNDF>=%d&&EicIdealGenTrack.fChiSquareCCDF>.001", ndfMostPopular-2);
  sprintf(cut, "");

  double par[100]; memset(par, 0x00, sizeof(par));

/////////////////////////////////////////////////////////////////////
// I have pretty much no idea what the code above this line means //
///////////////////////////////////////////////////////////////////

  const size_t nbin = 100;
  const double eta_min = -5.0;
  const double eta_max = 5.0;

#if 1
 
  gStyle->SetOptStat(0);
  
  TTreeFormula fRec("fRec","PidChargedCand@.GetEntries()",cbmsim);
  TTreeFormula fEtaTru("fEtaTru","MCTrack.GetMomentum().Eta()", cbmsim);
 
  double edge[nbin + 2]; // edge = [-5.0, -4.9, ..., 4.9, 5.0]
  const double interval = 0.1;

  // cast values for edge
  for (int i = 0; i < nbin + 1; i++){
	edge[i] = eta_min + i * interval;
  }

  int eta_rec_arr[nbin] = {0}; // eta_rec_arr[0] is the # of tracks whose eta_tru is in -5.0~-4.9 and rec == 1
  int eta_tru_arr[nbin] = {0}; // eta_tru_arr[0] is the # of tracks whose eta_tru is in -5.0~-4.9
  
  // cast values for eta_rec_arr & eta_tru_arr
  for (int j = 0; j < nEvt; j++) { 
	  
	  cbmsim->GetEntry(j);
	  const double eta_tru = fEtaTru.EvalInstance();
	  const int rec = fRec.EvalInstance();

          for (int i = 0; i < nbin; i++){
	     const double low = edge[i];
	     const double high = edge[i+1];
     
	     if (eta_tru >= low && eta_tru < high){
		     if (rec == 1){ // yes, eta_rec does not have to fall in the same eta interval to be reconstructed for efficiency plots
			 eta_rec_arr[i] += 1;
		     }

		     eta_tru_arr[i] += 1; 
 	     }
          }
   }
  
  // cout eta_tru_arr
  for (int i = 0; i < nbin; i++){
	cout << eta_tru_arr[i] << endl;  
  }
  
  cout << "+++++++++++++++++++++++++++++++" << endl;
  
  // cout eta_rec_arr
  for (int i = 0; i < nbin; i++){
	cout << eta_rec_arr[i] << endl;  
  }

#endif
}
