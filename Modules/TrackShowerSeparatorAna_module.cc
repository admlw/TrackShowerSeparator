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

  private:

    // Declare member data here.
    int run;
    int subrun;
    int event;
    bool isData;
    bool isSimulation;

    std::string fSelectionLabel;

    std::string fSelectionTPCObjAssn;

};


TrackShowerSeparatorAna::TrackShowerSeparatorAna(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");

  fSelectionLabel = p_labels.get<std::string>("SelectionLabel", "UBXSec");
  fSelectionTPCObjAssn = p_labels.get<std::string>("SelectionTPCObjAssn", "UBXSec");

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

  // first use the CC-inclusive selection to get the selected PFParticles

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

  // get selected TPC objects from selection handle
  art::FindManyP<ubana::TPCObject> selectedTPCObjects(selectionHandle, e, fSelectionTPCObjAssn);
  art::Ptr<ubana::TPCObject> selectedTPCObject;

  int nSelectedTPCObjects = selectedTPCObjects.at(0).size();

  // nSelectedTPCObjects expected to be 1
  if (nSelectedTPCObjects == 1){

    selectedTPCObject = selectedTPCObjects.at(0).at(0);

    // at time of writing, there's no handy function in the TPC object to
    // get the reconstructed showers, only tracks and PFParticles
    // to get around this, get PFParticles and then use assns to go from 
    // these to tracks and showers

    const std::vector<recob::PFParticle>& selectedPfps = selectedTPCObject->GetPFPs(); 

    for (size_t i = 0; i < selectedPfps.size(); i++){

      const recob::PFParticle& pfp = selectedPfps.at(i);

      std::cout << "Found PFParticle with id " << pfp.Self() << std::endl;

    }

  }

}

DEFINE_ART_MODULE(TrackShowerSeparatorAna)
