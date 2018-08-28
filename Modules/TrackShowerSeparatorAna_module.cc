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

// LArSoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ubxsec
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

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

    std::pair< const simb::MCParticle*, double > GetAssociatedMCParticle(std::vector< art::Ptr<recob::Hit> > hits, art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > particlesPerHit);

  private:

    // Declare member data here.
    int run;
    int subrun;
    int event;
    bool isData;
    bool isSimulation;

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
    std::string fTrackHitAssn;
    std::string fShowerHitAssn;
    std::string fHitTruthAssn;

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
  fPfpClusterAssn      = p_labels.get<std::string>("TrackClusterAssn", "pandoraNu::UBXSec");
  fTrackHitAssn        = p_labels.get<std::string>("TrackHitAssn", "pandoraNu::UBXSec"); 
  fShowerHitAssn       = p_labels.get<std::string>("ShowerHitAssn", "pandoraNu::UBXSec"); 
  fHitTruthAssn        = p_labels.get<std::string>("HitTruthAssn", "pandoraCosmicHitRemoval::UBXSec");

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
    throw std::exception();
  }

  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  e.getByLabel(fPfpLabel, pfpHandle);
  if (!pfpHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::PFParticle product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::PFParticle product not found." 
      << std::endl;
    throw std::exception();
  }

  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  if (!trackHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Track product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Track product not found." 
      << std::endl;
    throw std::exception();
  }

  art::Handle< std::vector<recob::Shower> > showerHandle;
  e.getByLabel(fShowerLabel, showerHandle);
  if (!showerHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Shower product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Shower product not found." 
      << std::endl;
    throw std::exception();
  }

  art::Handle< std::vector<recob::Hit> > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);
  if (!hitHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Hit product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Hit product not found." 
      << std::endl;
    throw std::exception();
  }

  art::Handle< std::vector<recob::Cluster> > clusterHandle;
  e.getByLabel(fClusterLabel, clusterHandle);
  if (!clusterHandle.isValid()){
    std::cout << "[TrackShowerSepAna] recob::Cluster product not found."
      << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) 
      << "[TrackShowerSepAna] recob::Cluster product not found." 
      << std::endl;
    throw std::exception();
  }



  // get selected TPC objects from selection handle
  art::FindManyP<ubana::TPCObject> selectedTPCObjects(selectionHandle, e, fSelectionTPCObjAssn);
  art::Ptr<ubana::TPCObject> selectedTPCObject;

  art::FindManyP<recob::Track> pfpToTrack(pfpHandle, e, fPfpTrackAssn);
  art::FindManyP<recob::Shower> pfpToShower(pfpHandle, e, fPfpShowerAssn);
  art::FindManyP<recob::Cluster> pfpToCluster(pfpHandle, e, fPfpClusterAssn);
  art::FindManyP<recob::Hit> trackToHit(trackHandle, e, fTrackHitAssn);
  art::FindManyP<recob::Hit> showerToHit(showerHandle, e, fShowerHitAssn);

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

    for (size_t i = 0; i < selectedPfps.size(); i++){

      std::cout << "getting pfp" << i << std::endl;

      const recob::PFParticle& pfp = selectedPfps.at(i);

      std::cout << "Found PFParticle with id " << pfp.Self() << std::endl;

      // check if PFP is track or shower like

      std::vector<art::Ptr<recob::Cluster>> clustersFromPfp = pfpToCluster.at(pfp.Self());
      
      // track like
      if (pfp.PdgCode() == 13){
        std::cout << "track like, picking up track for " << pfp.Self() << std::endl;
        std::vector<art::Ptr<recob::Track>> trackFromPFP = pfpToTrack.at(pfp.Self());
        std::cout << "got vector of ptrs to tracks" << std::endl;
        art::Ptr<recob::Track> thisTrack = trackFromPFP.at(0);
        std::cout << "got ptr to track with ID " << thisTrack->ID() << " and key " << thisTrack.key() << std::endl;
        std::vector<art::Ptr<recob::Hit>> hitsFromTrack = trackToHit.at(thisTrack.key()); 
        std::cout << "got vector of hits..." << std::endl;

        if (isSimulation){
          art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > 
            particlesPerHit(hitHandle, e, fHitTruthAssn);
                
          std::pair< const simb::MCParticle*, double > mcpInformation = 
            GetAssociatedMCParticle(hitsFromTrack, particlesPerHit);

          const simb::MCParticle* mcpMatchedParticle = mcpInformation.first;

          std::cout << mcpMatchedParticle->PdgCode() << std::endl;

        }

      }

      // shower like
      if (pfp.PdgCode() == 11){
        std::cout << "shower like, picking up shower for " << pfp.Self() << std::endl;
        std::vector<art::Ptr<recob::Shower>> showerFromPFP = pfpToShower.at(pfp.Self());
        std::cout << "got vector of ptrs to showers" << std::endl;
        art::Ptr<recob::Shower> thisShower = showerFromPFP.at(0);
        std::cout << "got ptr to shower with ID " << thisShower->ID() << " and key " << thisShower.key() << std::endl;
        std::vector<art::Ptr<recob::Hit>> hitsFromShower = showerToHit.at(thisShower.key()); 
        std::cout << "got vector of hits..." << std::endl;

        if (isSimulation){
          art::FindMany< simb::MCParticle, anab::BackTrackerHitMatchingData > 
            particlesPerHit(hitHandle, e, fHitTruthAssn);
               
          std::cout << "getting associated mcparticle..." << std::endl;
          std::pair< const simb::MCParticle*, double > mcpInformation = 
            GetAssociatedMCParticle(hitsFromShower, particlesPerHit);
          std::cout << "done..." << std::endl;

          const simb::MCParticle* mcpMatchedParticle = mcpInformation.first;

          std::cout << mcpMatchedParticle->PdgCode() << std::endl;

        }

      }

      std::cout << "there are " << clustersFromPfp.size() << "cluster associated to this Pfparticle" << std::endl;

    }

  }

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

DEFINE_ART_MODULE(TrackShowerSeparatorAna)
