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
  const double res_min = -0.3;
  const double res_max = 0.3;
  int nhit_arr[nhit_bin] = {0};

#if 1
  
  TTreeFormula fEtaTru("fEtaTru", "MCTrack.GetMomentum().Eta()", cbmsim);
  TTreeFormula fHit("fHit", "FstMoCaPoint@.GetEntries()+BstMoCaPoint@.GetEntries()+VstMoCaPoint@.GetEntries()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);

  TH2* hist = new TH2D("hist", Form("%0.1f < eta < %0.1f", eta_min, eta_max), nhit_bin, nhit_min, nhit_max, 12, res_min, res_max);
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

  /*
 * int tot = 0; // tot < Entries that shows up in the histogram
  for (int i = 0; i < nhit_bin; i++){
  	  cout << i << "-" << nhit_arr[i] << endl;
  	  tot += nhit_arr[i];
  }
  
  cout << "tot - " << tot << endl;

  int tot1 = hist->Integral(1,1);
  cout << "Integral1 - " << tot1 << endl;

  int tot2 = hist->Integral(2,2);
  cout << "Integral2 - " << tot2 << endl;

  int tot3 = hist->Integral(3,3);
  cout << "Integral3 - " << tot3 << endl;

  int tot4 = hist->Integral(4,4);
  cout << "Integral4 - " << tot4 << endl;

  int tot5 = hist->Integral(5,5);
  cout << "Integral5 - " << tot5 << endl;

  int tot6 = hist->Integral(6,6);
  cout << "Integral6 - " << tot6 << endl;

  int tot = hist->Integral(1, 6);
  cout << "Integral - " << tot << endl;

  int test1 = hist->GetBin(3);
  cout << "GetBin (3)- " << test1 << endl;

  int test1 = hist->GetBin(4);
  cout << "GetBin (4)- " << test1 << endl;

  int test1 = hist->GetBin(4, -0.03);
  cout << "GetBin (4, -0.03)- " << test1 << endl;

  int test2 = hist->GetBin(0,0);
  cout << "GetBin (0,0)- " << test2 << endl;

  int test3 = hist->GetBinContent(3, 5);
  cout << "GetBinContent (3, 5)- " << test3 << endl;

  int test3 = hist->GetBinContent(4, 5);
  cout << "GetBinContent (4, 5)- " << test3 << endl;

  int test3 = hist->GetBinContent(5, 5);
  cout << "GetBinContent (5, 5)- " << test3 << endl;

  int test3 = hist->GetBinContent(6, 5);
  cout << "GetBinContent (6, 5)- " << test3 << endl;

  int test4 = hist->GetBinContent(0,0);
  cout << "GetBinContent (0,0)" << test4 << endl;

  int test4 = hist->GetBinContent(1,1);
  cout << "GetBinContent (1,1)" << test4 << endl;

  int test4 = hist->GetBinContent(2,2);
  cout << "GetBinContent (2,2)" << test4 << endl;

  int test4 = hist->GetBinContent(3,3);
  cout << "GetBinContent (3,3)" << test4 << endl;
*/

  TFile f("output/default.root","recreate");

  gStyle->SetOptStat(10);
  gStyle->SetStatW(0.155);
  gStyle->SetStatH(0.255);
  
  hist->GetXaxis()->SetTitle("number of hits");
  hist->GetYaxis()->SetTitle("dp/p");
  hist->Draw("colz");

#endif
}
