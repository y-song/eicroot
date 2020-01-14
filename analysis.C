
#define _NDF_MAX_ 1000

void analysis()
{
  // Load basic libraries;
  gROOT->Macro("$VMCWORKDIR/gconfig/rootlogon.C");  

  // Input simulated & reconstructed files;
  TFile *ff = new TFile("simulation.root");
  TTree *cbmsim = ff->Get("cbmsim"); 
  cbmsim->AddFriend("cbmsim", "reconstruction.root");

  // Figure out max and most popular ndf value;
  TClonesArray *rctrk = new TClonesArray("PndPidCandidate");
  cbmsim->SetBranchAddress("PidChargedCand",&rctrk);
  unsigned nEvt = cbmsim->GetEntries(), ndfMax = 0;
  unsigned arr[1000]; memset(arr, 0x00, sizeof(arr));
  for (unsigned evt = 0; evt<nEvt; evt++) {
    cbmsim->GetEntry(evt);
    
    if (rctrk->GetEntriesFast()) {
      PndPidCandidate *rctrack = rctrk->At(0);
      
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

  char *expression = "(PidChargedCand.GetMomentum().Mag()-MCTrack.GetMomentum().Mag())/MCTrack.GetMomentum().Mag()";
  //float expression1 = (PidChargedCand.GetMomentum().Mag()-MCTrack.GetMomentum().Mag())/MCTrack.GetMomentum().Mag();
  char cut[1024];
  // Allow to lose 3 degrees of freedom compared to the max (assume standard) case;
  sprintf(cut, "EicIdealGenTrack.fNDF>=%d&&EicIdealGenTrack.fChiSquareCCDF>.001", ndfMostPopular-2);
  sprintf(cut, "");

  double par[100]; memset(par, 0x00, sizeof(par));

  const size_t nbin = 50;
  const double p_min = 0.5;
  const double p_max = 50.5;

#if 0
  TTreeFormula fp("fp", "MCTrack.GetMomentum().Mag()", cbmsim);
  TTreeFormula fx("fx", expression, cbmsim);
  TH1D mean_x("mean_x", "", nbin, p_min, p_max);
  TH1D count_p("count_p", "", nbin, p_min, p_max);
  mean_x.Sumw2();
  count_p.Sumw2();
  for (Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
	  cbmsim->GetEntry(i);
	  const double p = fp.EvalInstance();
	  const double x = fx.EvalInstance();
	  if (fabs(x) < 0.2) {
		  mean_x.Fill(p, x);
		  count_p.Fill(p);
	  }
  }
  mean_x.Divide(&count_p);
  TH1D ux2("ux2", "", nbin, p_min, p_max);
  TH1D ux4("ux4", "", nbin, p_min, p_max);
  ux2.Sumw2();
  ux4.Sumw2();
  for (Long64_t i = 0; i < cbmsim->GetEntries(); i++) {
	  cbmsim->GetEntry(i);
	  const double p = fp.EvalInstance();
	  const double x = fx.EvalInstance();
	  const double mean =
		  mean_x.GetBinContent(mean_x.GetXaxis()->FindFixBin(p));
	  const double d = x - mean;
	  if (fabs(x) < 0.2) {
		  ux2.Fill(p, d * d);
		  ux4.Fill(p, d * d * d * d);
	  }
  }
  ux2.Divide(&count_p);
  ux4.Divide(&count_p);
  TH1D *h1 = new TH1D("h1", "", nbin, p_min, p_max);
  h1->Sumw2();
  for (Int_t j = 1; j < nbin + 1; j++) {
	  const double n = count_p.GetBinContent(j);
	  if (n >= 3) {
		  const double u2 = ux2.GetBinContent(j);
		  const double u4 = ux4.GetBinContent(j);

		  h1->SetBinContent(j, sqrt((n - 1) * u2 / n));
		  h1->SetBinError(j,
						  sqrt((n - 1) *
							   ((n - 1) * u4 - (n - 3) * u2 * u2) /
							   (n * n * n)));
	  }
  }
  h1->SetMarkerStyle(20);
  h1->Draw("e1x0");
#else
  TH2D *pdp = new TH2D("pdp", "", nbin, p_min, p_max, 1000, -0.2, 0.2);
  pdp->Sumw2();
  cbmsim->Project("pdp", TString(expression) + ":MCTrack.GetMomentum().Mag()", cut);
  //cbmsim->Project("pdp", "MCTrack.GetMomentum().Eta():MCTrack.GetMomentum().Mag()", cut);
  pdp->Draw("colz");
#endif
}
