void makeShowerFractionAsFnOfEnergyPlots(){

  int chosenPDG = 211;
  double massPDG = 0.1395;
  std::string namePDG = "Pion";
  int nbins = 20;
  float binlow = 0;
  float binhigh = 1.0;

  bool pfpIsPandoraTrack_nom;
  bool pfpIsPandoraShower_nom;
  int truthMatchPDGCode_nom;
  float truthMatchEnergy_nom;

  bool pfpIsPandoraTrack_dic;
  bool pfpIsPandoraShower_dic;
  int truthMatchPDGCode_dic;
  float truthMatchEnergy_dic;

  TEfficiency* trackFrac_nom = new TEfficiency("trackFrac_nom", (namePDG+std::string(";True Energy;")).c_str(), nbins, binlow, binhigh);
  TEfficiency* trackFrac_dic = new TEfficiency("trackFrac_dic", (namePDG+std::string(";True Energy;")).c_str(), nbins, binlow, binhigh);

  TTree* t0 = (TTree*)_file0->Get("ana/analysis_tree");
  TTree* t1 = (TTree*)_file1->Get("ana/analysis_tree");

  t0->SetBranchAddress("pfpIsPandoraTrack", &pfpIsPandoraTrack_nom);
  t0->SetBranchAddress("pfpIsPandoraShower", &pfpIsPandoraShower_nom);
  t0->SetBranchAddress("truthMatchPDGCode", &truthMatchPDGCode_nom);
  t0->SetBranchAddress("truthMatchEnergy", &truthMatchEnergy_nom);

  t1->SetBranchAddress("pfpIsPandoraTrack", &pfpIsPandoraTrack_dic);
  t1->SetBranchAddress("pfpIsPandoraShower", &pfpIsPandoraShower_dic);
  t1->SetBranchAddress("truthMatchPDGCode", &truthMatchPDGCode_dic);
  t1->SetBranchAddress("truthMatchEnergy", &truthMatchEnergy_dic);

  for (int i_en = 0; i_en < t0->GetEntries(); i_en++){

    t0->GetEntry(i_en);

    if (std::abs(truthMatchPDGCode_nom) == chosenPDG){
      
      trackFrac_nom->Fill(pfpIsPandoraTrack_nom, truthMatchEnergy_nom-massPDG);

    }

  }

  for (int i_en = 0; i_en < t1->GetEntries(); i_en++){

    t1->GetEntry(i_en);

    if (std::abs(truthMatchPDGCode_dic) == chosenPDG){
      
      trackFrac_dic->Fill(pfpIsPandoraTrack_dic, truthMatchEnergy_dic-massPDG);

    }

  }

  TH1* nevents_nom = (TH1*)trackFrac_nom->GetCopyTotalHisto();
 
  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  TPad *topPad = new TPad("topPad", "", 0.005, 0.3, 0.995, 0.995);
  TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
  topPad->SetBottomMargin(0.02);
  bottomPad->SetTopMargin(0.0);
  bottomPad->SetBottomMargin(0.25);
  bottomPad->SetGridy();
  topPad->Draw();
  bottomPad->Draw();
  topPad->cd();

  nevents_nom->SetFillColor(kGray);
  nevents_nom->SetLineWidth(0);
  nevents_nom->Scale(1./nevents_nom->Integral());
  nevents_nom->GetYaxis()->SetRangeUser(0,1);
  nevents_nom->GetXaxis()->SetLabelOffset(1000);
  nevents_nom->SetTitle((namePDG+";;Fraction Reconstructed as Tracks").c_str());
  nevents_nom->Draw();
  trackFrac_nom->SetMarkerColor(kAzure+1);
  trackFrac_nom->SetLineColor(kAzure+1);
  trackFrac_nom->Draw("same");
  trackFrac_dic->SetMarkerColor(kGreen+1);
  trackFrac_dic->SetLineColor(kGreen+1);
  trackFrac_dic->Draw("same");

  gPad->RedrawAxis();

  bottomPad->cd();
  TH1* nevents_total_nom  = (TH1*)trackFrac_nom->GetCopyTotalHisto();
  TH1* nevents_passed_nom = (TH1*)trackFrac_nom->GetCopyPassedHisto();

  nevents_total_nom->Sumw2();
  nevents_passed_nom->Sumw2();

  nevents_passed_nom->Divide(nevents_total_nom);

  TH1* nevents_total_dic  = (TH1*)trackFrac_dic->GetCopyTotalHisto();
  TH1* nevents_passed_dic = (TH1*)trackFrac_dic->GetCopyPassedHisto();

  nevents_total_dic->Sumw2();
  nevents_passed_dic->Sumw2();

  nevents_passed_dic->Divide(nevents_total_dic);

  nevents_passed_dic->Divide(nevents_passed_nom);
  
  nevents_passed_dic->GetYaxis()->SetRangeUser(0.8, 1.2);
  nevents_passed_dic->SetTitle(";True E_{k} (GeV);DIC/Nominal");
  nevents_passed_dic->GetXaxis()->SetTitleSize(0.12);
  nevents_passed_dic->GetYaxis()->SetTitleSize(0.12);
  nevents_passed_dic->GetXaxis()->SetLabelSize(0.12);
  nevents_passed_dic->GetYaxis()->SetLabelSize(0.12);
  nevents_passed_dic->GetYaxis()->SetTitleOffset(0.38);
  nevents_passed_dic->GetYaxis()->SetNdivisions(305);

  nevents_passed_dic->Draw("sameE1");

  c1->SaveAs((namePDG+std::string("trackFrac.png")).c_str());

}
