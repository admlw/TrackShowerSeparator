////////////////////////////////////////////////////////////////////////
// Class:       TrackShowerSeparatorAna
// Plugin Type: analyzer (art v2_05_01)
// File:        TrackShowerSeparatorAna_module.cc
//
// Generated at Tue Aug 28 10:30:13 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
//
// This analyzer module uses the CC-inclusive selection to estimate 
// issues of track/shower separation for particle species of interest
// (proton, muon, pion, electrons and photons), and attempts to 
// find some variable to separate tracks from showers.
//
////////////////////////////////////////////////////////////////////////

// base 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// LArSoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"

// root
#include "TTree.h"

// ubxsec
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

// local
#include "uboone/TrackShowerSeparator/Algorithms/SPUtility.h"

class TrackShowerSeparatorAna;


class TrackShowerSeparatorAna : public art::EDAnalyzer {
  public:
    explicit TrackShowerSeparatorAna(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    TrackShowerSeparatorAna(TrackShowerSeparatorAna const &) = delete;
    TrackShowerSeparatorAna(TrackShowerSeparatorAna &&) = delete;
    TrackShowerSeparatorAna & operator = (TrackShowerSeparatorAna const &) = delete;
    TrackShowerSeparatorAna & operator = (TrackShowerSeparatorAna &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    void beginJob() override;

    std::pair< const simb::MCParticle*, double > GetAssociatedMCParticle(std::vector< art::Ptr<recob::Hit> > hits, art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit);

    void initialiseAnalysisTree(TTree* t, bool isSimulation);

    void prepareVariables();

  private:

    art::ServiceHandle< art::TFileService > tfs;
    art::ServiceHandle< geo::Geometry > geo;

    // local class instances
    trackshowerseparator::SPUtility _sputility;

    // tree
    TTree* ana_tree;

    int run;
    int subrun;
    int event;
    bool isData;
    bool isSimulation;

    // tree is filled on a per-pfp basis

    bool pfpIsPandoraTrack;
    bool pfpIsPandoraShower;

    float clusterStartCharge;
    float clusterStartAngle;
    float clusterStartOpeningAngle;
    float clusterEndCharge;
    float clusterEndAngle;
    float clusterEndOpeningAngle;
    float clusterIntegral;
    float clusterIntegralStdDev;
    float clusterIntegralAverage;
    float clusterSummedADC;
    float clusterSummedADCStdDev;
    float clusterSummedADCAverage;
    float clusterMultipleHitDensity;
    float clusterWidth;
    float clusterPlane;
    int clusterNHits;

    float averageSPDistance;

    std::vector<float> clusterStartCharge_v;
    std::vector<float> clusterStartAngle_v;
    std::vector<float> clusterStartOpeningAngle_v;
    std::vector<float> clusterEndCharge_v;
    std::vector<float> clusterEndAngle_v;
    std::vector<float> clusterEndOpeningAngle_v;
    std::vector<float> clusterIntegral_v;
    std::vector<float> clusterIntegralStdDev_v;
    std::vector<float> clusterIntegralAverage_v;
    std::vector<float> clusterSummedADC_v;
    std::vector<float> clusterSummedADCStdDev_v;
    std::vector<float> clusterSummedADCAverage_v;
    std::vector<float> clusterMultipleHitDensity_v;
    std::vector<float> clusterWidth_v;
    std::vector<float> clusterPlane_v;
    std::vector<int> clusterNHits_v;

    std::vector<std::vector<float>> hit_peakTimes_incm;
    std::vector<std::vector<float>> hit_wire_incm;

    int truthMatchPDGCode;
    float truthMatchEnergy;

    // fhicl
    std::string fSelectionLabel;
    std::string fPfpLabel;
    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fHitLabel;
    std::string fClusterLabel;

    std::string fSelectionTPCObjAssn;
    std::string fPfpTrackAssn;
    std::string fPfpShowerAssn;
    std::string fPfpClusterAssn;
    std::string fPfpSpacePointAssn;
    std::string fTrackHitAssn;
    std::string fShowerHitAssn;
    std::string fHitTruthAssn;
    std::string fClusterHitAssn;
    
    bool fIsData;

};


TrackShowerSeparatorAna::TrackShowerSeparatorAna(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");

  fSelectionLabel = p_labels.get<std::string>("SelectionLabel", "UBXSec");
  fPfpLabel       = p_labels.get<std::string>("PFPLabel", "pandoraNu::UBXSec");
  fTrackLabel     = p_labels.get<std::string>("TrackLabel", "pandoraNu::UBXSec");
  fShowerLabel    = p_labels.get<std::string>("ShowerLabel", "pandoraNu::UBXSec");
  fClusterLabel   = p_labels.get<std::string>("ClusterLabel", "pandoraNu::UBXSec");
  fHitLabel       = p_labels.get<std::string>("HitLabel", "pancoraCosmicHitRemoval::UBXSec");

  fSelectionTPCObjAssn = p_labels.get<std::string>("SelectionTPCObjAssn", "UBXSec");
  fPfpTrackAssn        = p_labels.get<std::string>("PFPTrackAssn", "pandoraNu::UBXSec");
  fPfpShowerAssn       = p_labels.get<std::string>("PFPShowerAssn", "pandoraNu::UBXSec");
  fPfpClusterAssn      = p_labels.get<std::string>("PFPClusterAssn", "pandoraNu::UBXSec");
  fPfpSpacePointAssn   = p_labels.get<std::string>("PFPSpacePointAssn", "pandoraNu::UBXSec");
  fTrackHitAssn        = p_labels.get<std::string>("TrackHitAssn", "pandoraNu::UBXSec"); 
  fShowerHitAssn       = p_labels.get<std::string>("ShowerHitAssn", "pandoraNu::UBXSec"); 
  fHitTruthAssn        = p_labels.get<std::string>("HitTruthAssn", "pandoraCosmicHitRemoval::UBXSec");
  fClusterHitAssn      = p_labels.get<std::string>("ClusterHitAssn", "pandoraNu::UBXSec");

  fIsData = p.get<bool>("IsData");

}

void TrackShowerSeparatorAna::analyze(art::Event const & e)
{
  run = e.run();
  subrun = e.subRun();
  event = e.event();
  isData = e.isRealData();
  isSimulation = !isData;

  std::cout << "[TrackShowerSepAna] ----- Processing Event " << run << "." 
    << subrun << "." << event << std::endl;
  std::cout << "isData: " << isData << " isSimulation: " << isSimulation << std::endl;

  // get handles to all of the information we're going to want to 
  // access via associations

  art::Handle< std::vector<ubana::SelectionResult> > selectionHandle;
  e.getByLabel(fSelectionLabel, selectionHandle);
  if (!selectionHandle.isValid()){
    std::cout << "[TrackShowerSepAna] ubana::SelectionResult product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] ubana::SelectionResult product not found." 
      << std::endl;
  }

  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  e.getByLabel(fPfpLabel, pfpHandle);
  if (!pfpHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::PFParticle product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::PFParticle product not found." 
      << std::endl;
  }

  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  if (!trackHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Track product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Track product not found." 
      << std::endl;
  }

  art::Handle< std::vector<recob::Shower> > showerHandle;
  e.getByLabel(fShowerLabel, showerHandle);
  if (!showerHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Shower product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Shower product not found." 
      << std::endl;
  }

  art::Handle< std::vector<recob::Hit> > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);
  if (!hitHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Hit product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Hit product not found." 
      << std::endl;
  }

  art::Handle< std::vector<recob::Cluster> > clusterHandle;
  e.getByLabel(fClusterLabel, clusterHandle);
  if (!clusterHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Cluster product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Cluster product not found." 
      << std::endl;
  }



  // get selected TPC objects from selection handle
  art::FindManyP<ubana::TPCObject> selectedTPCObjects(selectionHandle, e, fSelectionTPCObjAssn);
  art::Ptr<ubana::TPCObject> selectedTPCObject;

  art::FindManyP<recob::Track> pfpToTrack(pfpHandle, e, fPfpTrackAssn);
  art::FindManyP<recob::Shower> pfpToShower(pfpHandle, e, fPfpShowerAssn);
  art::FindManyP<recob::Cluster> pfpToCluster(pfpHandle, e, fPfpClusterAssn);
  art::FindManyP<recob::SpacePoint> pfpToSpacePoint(pfpHandle, e, fPfpSpacePointAssn);
  art::FindManyP<recob::Hit> trackToHit(trackHandle, e, fTrackHitAssn);
  art::FindManyP<recob::Hit> showerToHit(showerHandle, e, fShowerHitAssn);
  art::FindManyP<recob::Hit> clusterToHit(clusterHandle, e, fClusterHitAssn);

  int nSelectedTPCObjects = selectedTPCObjects.at(0).size();

  std::cout << "[TrackShowerSepAna] Number of selected TPC Objects: " << nSelectedTPCObjects << std::endl;

  // nSelectedTPCObjects expected to be 1
  if (nSelectedTPCObjects == 1){

    selectedTPCObject = selectedTPCObjects.at(0).at(0);

    // at time of writing, there's no handy function in the TPC object to
    // get the reconstructed showers, only tracks and PFParticles
    // to get around this, get PFParticles and then use assns to go from 
    // these to tracks and showers

    const std::vector<recob::PFParticle>& selectedPfps = selectedTPCObject->GetPFPs(); 
    const recob::Vertex vertex = selectedTPCObject->GetVertex();

    for (size_t i = 0; i < selectedPfps.size(); i++){

      prepareVariables();

      std::cout << "getting pfp" << i << std::endl;

      const recob::PFParticle& pfp = selectedPfps.at(i);

      std::cout << "Found PFParticle with id " << pfp.Self() << std::endl;

      // check if PFP is track or shower like

      // track like
      if (pfp.PdgCode() == 13){

        pfpIsPandoraTrack = true;
        pfpIsPandoraShower = false;

        std::cout << "track like, picking up track for " << pfp.Self() << std::endl;
        std::vector<art::Ptr<recob::Track>> trackFromPFP = pfpToTrack.at(pfp.Self());
        std::cout << "got vector of ptrs to tracks of size " << trackFromPFP.size() << std::endl;
        art::Ptr<recob::Track> thisTrack = trackFromPFP.at(0);
        std::cout << "got ptr to track with ID " << thisTrack->ID() << " and key " << thisTrack.key() << std::endl;
        std::vector<art::Ptr<recob::Hit>> hitsFromTrack = trackToHit.at(thisTrack.key()); 

        if (isSimulation){
          art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > 
            particlesPerHit(hitHandle, e, fHitTruthAssn);

          std::pair< const simb::MCParticle*, double > mcpInformation = 
            GetAssociatedMCParticle(hitsFromTrack, particlesPerHit);

          const simb::MCParticle* mcpMatchedParticle = mcpInformation.first;

          std::cout << "matched PDG code: " << mcpMatchedParticle->PdgCode() << std::endl;

          truthMatchPDGCode = mcpMatchedParticle->PdgCode();
          truthMatchEnergy = mcpMatchedParticle->E();

        }

      }

      // shower like
      if (pfp.PdgCode() == 11){

        pfpIsPandoraTrack = false;
        pfpIsPandoraShower = true;

        std::cout << "shower like, picking up shower for " << pfp.Self() << std::endl;
        std::vector<art::Ptr<recob::Shower>> showerFromPFP = pfpToShower.at(pfp.Self());
        std::cout << "got vector of ptrs to showers of size " << showerFromPFP.size() << std::endl;
        art::Ptr<recob::Shower> thisShower = showerFromPFP.at(0);
        std::cout << "got ptr to shower with ID " << thisShower->ID() << " and key " << thisShower.key() << std::endl;
        std::vector<art::Ptr<recob::Hit>> hitsFromShower = showerToHit.at(thisShower.key()); 

        if (isSimulation){
          art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > 
            particlesPerHit(hitHandle, e, fHitTruthAssn);

          std::cout << "getting associated mcparticle..." << std::endl;
          std::pair< const simb::MCParticle*, double > mcpInformation = 
            GetAssociatedMCParticle(hitsFromShower, particlesPerHit);
          std::cout << "done..." << std::endl;

          const simb::MCParticle* mcpMatchedParticle = mcpInformation.first;

          std::cout << "matched pdg code: " <<  mcpMatchedParticle->PdgCode() << std::endl;

          truthMatchPDGCode = mcpMatchedParticle->PdgCode();
          truthMatchEnergy = mcpMatchedParticle->E();

        }

      }

      // get clusters from PFParticle
      std::vector<art::Ptr<recob::Cluster>> clustersFromPfp = pfpToCluster.at(pfp.Self());
      std::cout << "there are " << clustersFromPfp.size() << " clusters associated to this Pfparticle" << std::endl;


      for (size_t i = 0; i < clustersFromPfp.size(); i++){

        art::Ptr<recob::Cluster> thisCluster  = clustersFromPfp.at(i);

        clusterStartCharge = thisCluster->StartCharge();
        clusterStartAngle = thisCluster->StartAngle();
        clusterStartOpeningAngle = thisCluster->StartOpeningAngle();
        clusterEndCharge = thisCluster->EndCharge();
        clusterEndAngle = thisCluster->EndAngle();
        clusterEndOpeningAngle = thisCluster->EndOpeningAngle();
        clusterIntegral = thisCluster->Integral();
        clusterIntegralStdDev = thisCluster->IntegralStdDev();
        clusterIntegralAverage = thisCluster->IntegralAverage();
        clusterSummedADC = thisCluster->SummedADC();
        clusterSummedADCStdDev = thisCluster->SummedADCstdDev();
        clusterSummedADCAverage = thisCluster->SummedADCaverage();
        clusterMultipleHitDensity = thisCluster->MultipleHitDensity();
        clusterWidth = thisCluster->Width();
        clusterPlane = thisCluster->Plane().Plane;
        clusterNHits = thisCluster->NHits();

        clusterStartCharge_v.at(clusterPlane) = clusterStartCharge;
        clusterStartAngle_v.at(clusterPlane) = clusterStartAngle;
        clusterStartOpeningAngle_v.at(clusterPlane) = clusterStartOpeningAngle;
        clusterEndCharge_v.at(clusterPlane) = clusterEndCharge;
        clusterEndAngle_v.at(clusterPlane) = clusterEndAngle;
        clusterEndOpeningAngle_v.at(clusterPlane) = clusterEndOpeningAngle;
        clusterIntegral_v.at(clusterPlane) = clusterIntegralAverage;
        clusterIntegralStdDev_v.at(clusterPlane) = clusterIntegralStdDev;
        clusterIntegralAverage_v.at(clusterPlane) = clusterIntegralAverage;
        clusterSummedADC_v.at(clusterPlane) = clusterSummedADC;
        clusterSummedADCStdDev_v.at(clusterPlane) = clusterSummedADCStdDev;
        clusterSummedADCAverage_v.at(clusterPlane) = clusterSummedADCAverage;
        clusterMultipleHitDensity_v.at(clusterPlane) = clusterMultipleHitDensity;
        clusterWidth_v.at(clusterPlane) = clusterWidth;
        clusterPlane_v.at(clusterPlane) = clusterPlane;
        clusterNHits_v.at(clusterPlane) = clusterNHits;

        std::cout << "----- Printing cluster information -----" << std::endl;
        std::cout << "clusterStartCharge " << clusterStartCharge << std::endl; 
        std::cout << "clusterStartAngle "  << clusterStartAngle << std::endl; 
        std::cout << "clusterStartOpeningAngle " << clusterStartOpeningAngle << std::endl; 
        std::cout << "clusterEndCharge "   << clusterEndCharge << std::endl; 
        std::cout << "clusterEndAngle "    << clusterEndAngle << std::endl; 
        std::cout << "clusterEndOpeningAngle "     << clusterEndOpeningAngle << std::endl; 
        std::cout << "clusterIntegral "    << clusterIntegral << std::endl;
        std::cout << "clusterIntegralStdDev " << clusterIntegralStdDev << std::endl;
        std::cout << "clusterIntegralAverage " << clusterIntegralAverage << std::endl;
        std::cout << "clusterSummedADC "    << clusterSummedADC << std::endl;
        std::cout << "clusterSummedADCStdDev " << clusterSummedADCStdDev << std::endl;
        std::cout << "clusterSummedADCAverage " << clusterSummedADCAverage << std::endl;
        std::cout << "clusterMultipleHitDensity " << clusterMultipleHitDensity << std::endl;
        std::cout << "clusterWidth " << clusterWidth << std::endl;
        std::cout << "clusterPlane " << clusterPlane << std::endl;
        std::cout << "clusterNHits " << clusterNHits << std::endl;

        std::vector<art::Ptr<recob::Hit>> hitsFromCluster = clusterToHit.at(thisCluster.key()); 

        std::vector<float> hit_peakTimes_incm_subv = {-1.};
        std::vector<float> hit_wire_incm_subv = {-1.};

        for (size_t i_hit = 0; i_hit < hitsFromCluster.size(); i_hit++){

          art::Ptr<recob::Hit> thisHit = hitsFromCluster.at(i);

          float wire_incm = -1;
          if (clusterPlane == 0)
            wire_incm = thisHit->Channel() * 0.3;
          else if (clusterPlane == 1)
            wire_incm = (thisHit->Channel()-2400) * 0.3;
          else if (clusterPlane == 2)
            wire_incm = (thisHit->Channel()-4800) * 0.3;

          hit_peakTimes_incm_subv.push_back(((thisHit->PeakTime()-(3125./2.))/2300.) *(2*geo->DetHalfWidth()));

          hit_wire_incm_subv.push_back(wire_incm);

        }

        hit_wire_incm.at(clusterPlane) = hit_wire_incm_subv;
        hit_peakTimes_incm.at(clusterPlane) = hit_peakTimes_incm_subv;


      }

      // get PFP-associated space points
      std::vector<art::Ptr<recob::SpacePoint>> spacepointFromPfparticle = pfpToSpacePoint.at(pfp.Self());
      std::vector<art::Ptr<recob::SpacePoint>> sortedSpacePointCollection = _sputility.getSortedSPList(spacepointFromPfparticle, vertex);

      averageSPDistance = 0;
      if (sortedSpacePointCollection.size() > 1){
        for (size_t i_sp = 1; i_sp < sortedSpacePointCollection.size(); i_sp++){

          art::Ptr<recob::SpacePoint> thisSP = sortedSpacePointCollection.at(i_sp);
          art::Ptr<recob::SpacePoint> prevSP = sortedSpacePointCollection.at(i_sp-1);

          double thisSPxyz[3] = {thisSP->XYZ()[0], thisSP->XYZ()[1], thisSP->XYZ()[2]};
          double prevSPxyz[3] = {prevSP->XYZ()[0], prevSP->XYZ()[1], prevSP->XYZ()[2]};

          averageSPDistance += _sputility.get3DLength(thisSPxyz, prevSPxyz);

        }

        averageSPDistance = averageSPDistance/(float)sortedSpacePointCollection.size();

      }
      else averageSPDistance = -1;

      // each entry is for one PFP
      ana_tree->Fill();
    }

  }
}

void TrackShowerSeparatorAna::beginJob()
{

  ana_tree = tfs->make<TTree>("analysis_tree", "analysis tree");
  initialiseAnalysisTree(ana_tree, !fIsData);

}

std::pair< const simb::MCParticle*, double > TrackShowerSeparatorAna::GetAssociatedMCParticle(std::vector< art::Ptr< recob::Hit > > hits, art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit)
{

  std::unordered_map<int,double> trkide;
  double maxe=-1, tote=0;
  simb::MCParticle const* maxp_me = NULL; //pointer for the particle match we will calculate

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

  //loop only over our hits
  for(size_t i_h=0; i_h<hits.size(); ++i_h){

    particle_vec.clear(); match_vec.clear();
    particlesPerHit.get(hits[i_h].key(),particle_vec,match_vec);
    //the .key() gives us the index in the original collection

    //loop over particles
    for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
      trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
      tote += match_vec[i_p]->energy; //calculate total energy deposited
      if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
        maxe = trkide[ particle_vec[i_p]->TrackId() ];
        maxp_me = particle_vec[i_p];
      }
    }//end loop over particles per hit

  }

  std::pair< const simb::MCParticle* , double > returner;
  returner.first  = maxp_me;
  returner.second = maxe/tote;

  return returner;

}

void TrackShowerSeparatorAna::initialiseAnalysisTree(TTree* tree, bool isSimulation){

  tree->Branch("run", &run);
  tree->Branch("subrun", &subrun);
  tree->Branch("event", &event);
  tree->Branch("isData", &isData);
  tree->Branch("pfpIsPandoraTrack", &pfpIsPandoraTrack);
  tree->Branch("pfpIsPandoraShower", &pfpIsPandoraShower);
  tree->Branch("clusterStartCharge", "std::vector<float>", &clusterStartCharge_v);
  tree->Branch("clusterStartAngle", "std::vector<float>", &clusterStartAngle_v);
  tree->Branch("clusterStartOpeningAngle", "std::vector<float>", &clusterStartOpeningAngle_v);
  tree->Branch("clusterEndCharge", "std::vector<float>", &clusterEndCharge_v);
  tree->Branch("clusterEndAngle", "std::vector<float>", &clusterEndAngle_v);
  tree->Branch("clusterEndOpeningAngle", "std::vector<float>", &clusterEndOpeningAngle_v);
  tree->Branch("clusterIntegral", "std::vector<float>", &clusterIntegral_v);
  tree->Branch("clusterIntegralStdDev", "std::vector<float>", &clusterIntegralStdDev_v);
  tree->Branch("clusterIntegralAverage", "std::vector<float>", &clusterIntegralAverage_v);
  tree->Branch("clusterSummedADC", "std::vector<float>", &clusterSummedADC_v);
  tree->Branch("clusterSummedADCStdDev", "std::vector<float>", &clusterSummedADCStdDev_v);
  tree->Branch("clusterSummedADCAverage", "std::vector<float>", &clusterSummedADCAverage_v);
  tree->Branch("clusterMultipleHitDensity", "std::vector<float>", &clusterMultipleHitDensity_v);
  tree->Branch("clusterWidth", "std::vector<float>", &clusterWidth_v);
  tree->Branch("clusterPlane", "std::vector<float>", &clusterPlane_v);
  tree->Branch("clusterNHits", "std::vector<int>", &clusterNHits_v);

  tree->Branch("averageSPDistance", &averageSPDistance);

  tree->Branch("hit_peakTimes_incm", "std::vector<std::vector<float>>", &hit_peakTimes_incm);
  tree->Branch("hit_wire_incm", "std::vector<std::vector<float>>", &hit_wire_incm);

  if (isSimulation){

    tree->Branch("truthMatchPDGCode", &truthMatchPDGCode);
    tree->Branch("truthMatchEnergy", &truthMatchEnergy);

  }

}

void TrackShowerSeparatorAna::prepareVariables(){

  pfpIsPandoraTrack = false;
  pfpIsPandoraShower = false;

  clusterStartCharge_v.resize(3,-1);
  clusterStartAngle_v.resize(3,-1);
  clusterStartOpeningAngle_v.resize(3,-1);
  clusterEndCharge_v.resize(3,-1);
  clusterEndAngle_v.resize(3,-1);
  clusterEndOpeningAngle_v.resize(3,-1);
  clusterIntegral_v.resize(3,-1);
  clusterIntegralStdDev_v.resize(3,-1);
  clusterIntegralAverage_v.resize(3,-1);
  clusterSummedADC_v.resize(3,-1);
  clusterSummedADCStdDev_v.resize(3,-1);
  clusterSummedADCAverage_v.resize(3,-1);
  clusterMultipleHitDensity_v.resize(3,-1);
  clusterWidth_v.resize(3,-1);
  clusterPlane_v.resize(3,-1);
  clusterNHits_v.resize(3,-1);

  hit_peakTimes_incm.resize(0);
  hit_wire_incm.resize(0);

  averageSPDistance = -1;

  std::vector<float> hit_peakTimes_incm_subv = {-1.};
  std::vector<float> hit_wire_incm_subv = {-1.};

  // push back three, one for each plane
  hit_peakTimes_incm.push_back(hit_peakTimes_incm_subv);
  hit_peakTimes_incm.push_back(hit_peakTimes_incm_subv);
  hit_peakTimes_incm.push_back(hit_peakTimes_incm_subv);
  hit_wire_incm.push_back(hit_wire_incm_subv);
  hit_wire_incm.push_back(hit_wire_incm_subv);
  hit_wire_incm.push_back(hit_wire_incm_subv);

  truthMatchPDGCode = 0;
  truthMatchEnergy = -1;

}

DEFINE_ART_MODULE(TrackShowerSeparatorAna)
