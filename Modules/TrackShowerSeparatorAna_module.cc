////////////////////////////////////////////////////////////////////////
// Class:       TrackShowerSeparatorAna
// Plugin Type: analyzer (art v2_05_01)
// File:        TrackShowerSeparatorAna_module.cc
//
// Generated at Tue Aug 28 10:30:13 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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

};


TrackShowerSeparatorAna::TrackShowerSeparatorAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void TrackShowerSeparatorAna::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

DEFINE_ART_MODULE(TrackShowerSeparatorAna)
