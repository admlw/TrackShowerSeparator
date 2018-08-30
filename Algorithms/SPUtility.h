#ifndef SPUTILITY_H
#define SPUTILITY_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TMath.h"
#include <vector>
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

namespace trackshowerseparator{

  class SPUtility{

    public:

      /**
       * Takes in a vector of art::Ptr<recob::SpacePoint>, and returns a sorted vector
       * of the same type, where the sorting criteria is distance from the reconstructed
       * vertex
       */
      std::vector<art::Ptr<recob::SpacePoint>> getSortedSPList(std::vector<art::Ptr<recob::SpacePoint>> spColletion, const recob::Vertex vertex);

      /**
       * Return 3d length between two points
       */
      float get3DLength(double* vertex_xyz, double* sp_xyz);

  };

}

#endif
