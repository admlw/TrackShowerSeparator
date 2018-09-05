#include "plotStackedHistograms.h"

// ROOT script to plot comparisons of cluster and spacepoint level variables from
// the TrackShowerSeparator module output

std::vector< std::vector<float> > getVarsToPlot(vars* varstoset){

  std::vector<std::vector<float>> returno;

  // convert single valued branches to vectors
  std::vector<float> spNSpacePoints_v = {(float)varstoset->spNSpacePoints};
  std::vector<float> spAverageDistance_v = {varstoset->spAverageDistance};
  std::vector<float> spMedianDistance_v = {varstoset->spMedianDistance};
  std::vector<float> spRMSDistance_v = {varstoset->spRMSDistance};
  std::vector<float> pfpIsPandoraTrack_v = {(float)((int)varstoset->pfpIsPandoraTrack)};
  std::vector<float> pfpIsPandoraShower_v = {(float)((int)varstoset->pfpIsPandoraShower)};

  // convert branches of vectors of ints to vectors of floats
  std::vector<float> clusterPlane_v((*(varstoset->clusterPlane)).begin(), ((*varstoset->clusterPlane).end()));
  std::vector<float> clusterNHits_v((*(varstoset->clusterNHits)).begin(), ((*varstoset->clusterNHits).end()));

  returno.push_back(pfpIsPandoraTrack_v);
  returno.push_back(pfpIsPandoraShower_v);
  returno.push_back(spNSpacePoints_v);
  returno.push_back(spAverageDistance_v);
  returno.push_back(spMedianDistance_v);
  returno.push_back(spRMSDistance_v);
  returno.push_back(*(varstoset->clusterStartCharge));
  returno.push_back(*(varstoset->clusterStartAngle));
  returno.push_back(*(varstoset->clusterStartOpeningAngle));
  returno.push_back(*(varstoset->clusterEndCharge));
  returno.push_back(*(varstoset->clusterEndAngle));
  returno.push_back(*(varstoset->clusterEndOpeningAngle));
  returno.push_back(*(varstoset->clusterIntegral));
  returno.push_back(*(varstoset->clusterIntegralStdDev));
  returno.push_back(*(varstoset->clusterIntegralAverage));
  returno.push_back(*(varstoset->clusterSummedADC));
  returno.push_back(*(varstoset->clusterSummedADCStdDev));
  returno.push_back(*(varstoset->clusterSummedADCAverage));
  returno.push_back(*(varstoset->clusterMultipleHitDensity));
  returno.push_back(*(varstoset->clusterWidth));
  returno.push_back(clusterPlane_v);
  returno.push_back(clusterNHits_v);

  return returno;

}

void initialiseBranches(TTree* tree, vars* varstoset){

  tree->SetBranchAddress("pfpIsPandoraTrack", &(varstoset->pfpIsPandoraTrack));
  tree->SetBranchAddress("pfpIsPandoraShower", &(varstoset->pfpIsPandoraShower));
  tree->SetBranchAddress("clusterStartCharge", &(varstoset->clusterStartCharge));
  tree->SetBranchAddress("clusterStartAngle", &(varstoset->clusterStartAngle));
  tree->SetBranchAddress("clusterStartOpeningAngle", &(varstoset->clusterStartOpeningAngle));
  tree->SetBranchAddress("clusterEndCharge", &(varstoset->clusterEndCharge));
  tree->SetBranchAddress("clusterEndAngle", &(varstoset->clusterEndAngle));
  tree->SetBranchAddress("clusterEndOpeningAngle", &(varstoset->clusterEndOpeningAngle));
  tree->SetBranchAddress("clusterIntegral", &(varstoset->clusterIntegral));
  tree->SetBranchAddress("clusterIntegralStdDev", &(varstoset->clusterIntegralStdDev));
  tree->SetBranchAddress("clusterIntegralAverage", &(varstoset->clusterIntegralAverage));
  tree->SetBranchAddress("clusterSummedADC", &(varstoset->clusterSummedADC));
  tree->SetBranchAddress("clusterSummedADCStdDev", &(varstoset->clusterSummedADCStdDev));
  tree->SetBranchAddress("clusterSummedADCAverage", &(varstoset->clusterSummedADCAverage));
  tree->SetBranchAddress("clusterMultipleHitDensity", &(varstoset->clusterMultipleHitDensity));
  tree->SetBranchAddress("clusterWidth", &(varstoset->clusterWidth));
  tree->SetBranchAddress("clusterPlane", &(varstoset->clusterPlane));
  tree->SetBranchAddress("clusterNHits", &(varstoset->clusterNHits));
  tree->SetBranchAddress("spNSpacePoints", &(varstoset->spNSpacePoints));
  tree->SetBranchAddress("spAverageDistance", &(varstoset->spAverageDistance));
  tree->SetBranchAddress("spMedianDistance", &(varstoset->spMedianDistance));
  tree->SetBranchAddress("spRMSDistance", &(varstoset->spRMSDistance));

  tree->SetBranchAddress("truthMatchPDGCode", &(varstoset->truthMatchPDGCode));
}


void plotStackedHistograms(){

  std::vector<std::string> plotNames = {
    "pfpIsPandoraTrack",
    "pfpIsPandoraShower",
    "spNSpacePoints",
    "spAverageDistance",
    "spMedianDistance",
    "spRMSDistance",
    "clusterStartCharge",
    "clusterStartAngle",
    "clusterStartOpeningAngle",
    "clusterEndCharge",
    "clusterEndAngle",
    "clusterEndOpeningAngle",
    "clusterIntegral",
    "clusterIntegralStdDev",
    "clusterIntegralAverage",
    "clusterSummedADC",
    "clusterSummedADCStdDev",
    "clusterSummedADCAverage",
    "clusterMultipleHitDensity",
    "clusterWidth",
    "clusterPlane",
    "clusterNHits"
  };

  std::vector<std::string> plotTitles = {
    ";PFP is PandoraTrack;",
    ";PFP is PandoraShower;",
    ";PFP N Spacepoints;",
    ";PFP Spacepoint Average Distance;",
    ";PFP Spacepoint Median Distance;",
    ";PFP Spacepoint RMS Distance;",
    ";Cluster Start Charge;",
    ";Cluster Start Angle;",
    ";Cluster Start Opening Angle;",
    ";Cluster End Charge;",
    ";Cluster End Angle;",
    ";Cluster End Opening Angle;",
    ";Cluster Integral;",
    ";Cluster Integral Std Dev;",
    ";Cluster Integral Average;",
    ";Cluster Summed ADC;",
    ";Cluster Summed ADC Std Dev;",
    ";Cluster Summed ADC Average;",
    ";Cluster Multiple Hit Density;",
    ";Cluster Width;",
    ";Cluster Plane;",
    ";Cluster NHits;"
  };

  std::vector<std::vector<float>> plotBinning = {
    {2, 0, 2},       // pfpIsPandoraTrack
    {2, 0, 2},       // pfpIsPandoraShower
    {100, 0, 300},   // spNSpacePoints
    {50, 0, 0.6},    // spAverageDistance
    {50, 0, 0.6},    // spMedianDistance
    {50, 0, 0.75},   // spRMSDistance
    {50, 0, 20000},  // clusterStartCharge
    {50, -1.5, 1.5}, // clusterStartAngle
    {50, 0, 3.15},   // clusterStartOpeningCharge
    {50, 0, 20000},  // clusterEndCharge
    {50, -1.5, 1.5}, // clusterEndAngle
    {50, 0, 3.15},   // clusterEndOpeningAngle
    {50, 0, 1500},   // clusterIntegral
    {50, 0, 1000},   // clusterIntegralStdDev
    {50, 0, 1500},   // clusterIntegralAverage
    {50, 0, 50000},  // clusterSummedADC
    {50, 0, 2000},   // clusterSummedADCStdDev
    {50, 0, 3000},   // clusterSummedADCAverage
    {50, 0, 1},      // clusterMultipleHitDensity
    {20, 0, 20},     // clusterWidth
    {3, 0, 3},       // clusterPlane
    {100, 0, 300}    // clusterNHits             
  };

  std::vector<bool> plotLegend = {
   true,   // pfpIsPandoraTrack          
   false,    // pfpIsPandoraShower
   false,   // spNSpacePoints
   false,   // spAverageDistance
   false,   // spMedianDistance
   false,   // spRMSDistance
   false,   // clusterStartCharge
   false,   // clusterStartAngle
   false,   // clusterStartOpeningCharge
   false,   // clusterEndCharge
   false,   // clusterEndAngle
   false,   // clusterEndOpeningAngle
   false,   // clusterIntegral
   false,   // clusterIntegralStdDev
   false,   // clusterIntegralAverage
   false,   // clusterSummedADC
   false,   // clusterSummedADCStdDev
   false,   // clusterSummedADCAverage
   false,   // clusterMultipleHitDensity
   false,   // clusterWidth
   false,   // clusterPlane
   false    // clusterNHits             
  };

  if ( (plotNames.size() != plotTitles.size()) || (plotNames.size() != plotBinning.size()) ){

    std::cout << "Something done. Not good." << std::endl;
    std::cout << "plotNames.size(): " << plotNames.size() << std::endl;
    std::cout << "plotTitles.size(): " << plotTitles.size() << std::endl;
    std::cout << "plotBinning.size(): " << plotBinning.size() << std::endl;

    throw std::exception();

  }

  // initialise histograms

  std::vector<TH1D*> plotCollection_proton;
  std::vector<TH1D*> plotCollection_muon;
  std::vector<TH1D*> plotCollection_pion;
  std::vector<TH1D*> plotCollection_kaon;
  std::vector<TH1D*> plotCollection_electron;
  std::vector<TH1D*> plotCollection_other;

  for (size_t i = 0; i < plotNames.size(); i++){

    plotCollection_proton.push_back(new TH1D(std::string(plotNames.at(i)+"_proton").c_str(), plotTitles.at(i).c_str(), (int)plotBinning.at(i).at(0), plotBinning.at(i).at(1), plotBinning.at(i).at(2)));
    plotCollection_muon.push_back(new TH1D(std::string(plotNames.at(i)+"_muon").c_str(), plotTitles.at(i).c_str(), (int)plotBinning.at(i).at(0), plotBinning.at(i).at(1), plotBinning.at(i).at(2)));
    plotCollection_pion.push_back(new TH1D(std::string(plotNames.at(i)+"_pion").c_str(), plotTitles.at(i).c_str(), (int)plotBinning.at(i).at(0), plotBinning.at(i).at(1), plotBinning.at(i).at(2)));
    plotCollection_kaon.push_back(new TH1D(std::string(plotNames.at(i)+"_kaon").c_str(), plotTitles.at(i).c_str(), (int)plotBinning.at(i).at(0), plotBinning.at(i).at(1), plotBinning.at(i).at(2)));
    plotCollection_electron.push_back(new TH1D(std::string(plotNames.at(i)+"_electron").c_str(), plotTitles.at(i).c_str(), (int)plotBinning.at(i).at(0), plotBinning.at(i).at(1), plotBinning.at(i).at(2)));
    plotCollection_other.push_back(new TH1D(std::string(plotNames.at(i)+"_other").c_str(), plotTitles.at(i).c_str(), (int)plotBinning.at(i).at(0), plotBinning.at(i).at(1), plotBinning.at(i).at(2)));

  }

  // get trees to plot from
  TTree* t1 = (TTree*) _file0->Get("ana/analysis_tree");

  vars variables_1;
  initialiseBranches(t1, &variables_1);

  for (int i = 0; i < t1->GetEntries(); i++){
    t1->GetEntry(i);
    std::vector<std::vector<float>> varsToPlot = getVarsToPlot(&variables_1);

    for (int i_var = 0; i_var < varsToPlot.size(); i_var++){

      if (std::abs(variables_1.truthMatchPDGCode) == 2212){
        for (int i_var_it = 0; i_var_it < varsToPlot.at(i_var).size(); i_var_it++){
          plotCollection_proton.at(i_var)->Fill(varsToPlot.at(i_var).at(i_var_it));
        }
      }

      else if (std::abs(variables_1.truthMatchPDGCode) == 13){
        for (int i_var_it = 0; i_var_it < varsToPlot.at(i_var).size(); i_var_it++){
          plotCollection_muon.at(i_var)->Fill(varsToPlot.at(i_var).at(i_var_it));
        }
      }

      else if (std::abs(variables_1.truthMatchPDGCode) == 211){
        for (int i_var_it = 0; i_var_it < varsToPlot.at(i_var).size(); i_var_it++){
          plotCollection_pion.at(i_var)->Fill(varsToPlot.at(i_var).at(i_var_it));
        }
      }

      else if (std::abs(variables_1.truthMatchPDGCode) == 321){
        for (int i_var_it = 0; i_var_it < varsToPlot.at(i_var).size(); i_var_it++){
          plotCollection_kaon.at(i_var)->Fill(varsToPlot.at(i_var).at(i_var_it));
        }
      }

      else if (std::abs(variables_1.truthMatchPDGCode) == 11){
        for (int i_var_it = 0; i_var_it < varsToPlot.at(i_var).size(); i_var_it++){
          plotCollection_electron.at(i_var)->Fill(varsToPlot.at(i_var).at(i_var_it));
        }
      }

      else {
        for (int i_var_it = 0; i_var_it < varsToPlot.at(i_var).size(); i_var_it++){
          plotCollection_other.at(i_var)->Fill(varsToPlot.at(i_var).at(i_var_it));
        }
      }

    }

  }

  for (int i = 0; i < plotCollection_proton.size(); i++){

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
    c1->cd();

    plotCollection_muon.at(i)->SetFillColor(229);
    plotCollection_proton.at(i)->SetFillColor(230);
    plotCollection_pion.at(i)->SetFillColor(231);
    plotCollection_kaon.at(i)->SetFillColor(232);
    plotCollection_electron.at(i)->SetFillColor(233);
    plotCollection_other.at(i)->SetFillColor(234);
    plotCollection_muon.at(i)->SetLineWidth(0);
    plotCollection_proton.at(i)->SetLineWidth(0);
    plotCollection_pion.at(i)->SetLineWidth(0);
    plotCollection_kaon.at(i)->SetLineWidth(0);
    plotCollection_electron.at(i)->SetLineWidth(0);
    plotCollection_other.at(i)->SetLineWidth(0);
    plotCollection_muon.at(i)->SetMarkerStyle(0);
    plotCollection_proton.at(i)->SetMarkerStyle(0);
    plotCollection_pion.at(i)->SetMarkerStyle(0);
    plotCollection_kaon.at(i)->SetMarkerStyle(0);
    plotCollection_electron.at(i)->SetMarkerStyle(0);
    plotCollection_other.at(i)->SetMarkerStyle(0);

    THStack *hs = new THStack("hs", "hs");
    hs->Add(plotCollection_proton.at(i));
    hs->Add(plotCollection_muon.at(i));
    hs->Add(plotCollection_pion.at(i));
    hs->Add(plotCollection_kaon.at(i));
    hs->Add(plotCollection_electron.at(i));
    hs->Add(plotCollection_other.at(i));

    TH1D* hTotal = new TH1D("hTotal", ";;", plotBinning.at(i).at(0), plotBinning.at(i).at(1), plotBinning.at(i).at(2));
    hTotal->Add(plotCollection_proton.at(i));
    hTotal->Add(plotCollection_muon.at(i));
    hTotal->Add(plotCollection_pion.at(i));
    hTotal->Add(plotCollection_kaon.at(i));
    hTotal->Add(plotCollection_electron.at(i));
    hTotal->Add(plotCollection_other.at(i));

    hs->SetTitle(plotTitles.at(i).c_str());
    
    if (plotNames.at(i) == "clusterStartAngle" || plotNames.at(i) == "clusterEndAngle" || plotNames.at(i) == "clusterPlane"){
      TAxis* yaxis = hTotal->GetYaxis();
      yaxis->SetRangeUser(0, hTotal->GetMaximum()*2);
      hTotal->Draw("E2");
      hs->Draw("same");
    }
    else hs->Draw();


    hTotal->SetFillStyle(3345);
    hTotal->SetFillColor(kBlack);
    hTotal->Draw("E2same");

    TLegend *leg = new TLegend();
    if (plotLegend.at(i) == true){
      leg->SetX1NDC(0.15);
      leg->SetX2NDC(0.5);
      leg->SetY1NDC(0.59);
      leg->SetY2NDC(0.89);
    }
    else {
      leg->SetX1NDC(0.5);
      leg->SetX2NDC(0.85);
      leg->SetY1NDC(0.59);
      leg->SetY2NDC(0.89);
    }

    leg->AddEntry(plotCollection_proton.at(i), "Proton");
    leg->AddEntry(plotCollection_muon.at(i), "Muon");
    leg->AddEntry(plotCollection_pion.at(i), "Pion");
    leg->AddEntry(plotCollection_kaon.at(i), "Kaon");
    leg->AddEntry(plotCollection_electron.at(i), "Electron");
    leg->AddEntry(plotCollection_other.at(i), "Other");
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->Draw("same");

    c1->SaveAs((std::string(plotNames.at(i)+".png")).c_str());

    c1->Delete();
  }

}


