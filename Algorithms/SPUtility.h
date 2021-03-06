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

#include "TVector3.h"
#include "TMath.h"
#include <vector>
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

namespace trackshowerseparator{

  class SPUtility{

    public:

      /**
       * checks a vector to find duplicate spacepoints and prints to screen
       */
      void checkVectorForDuplicates(std::vector<art::Ptr<recob::SpacePoint>> spCollection);

      /**
       * Takes in a vector of art::Ptr<recob::SpacePoint>, and returns a sorted vector
       * of the same type, where the sorting criteria is distance from the reconstructed
       * vertex
       */
      std::vector<art::Ptr<recob::SpacePoint>> getSortedSPList(
          std::vector<art::Ptr<recob::SpacePoint>> spCollection, 
          const recob::Vertex vertex, 
          std::string sortType);

      /**
       * Return 3d length between two points
       */
      float get3DLength(double* vertex_xyz, double* sp_xyz);

      /**
       * Returns the number of spacepoints in a chosen shell
       * of defined width "shellWidth"
       */
      int FindNSPInShell(
          float shellWidth, 
          int shellNumber, 
          std::vector<art::Ptr<recob::SpacePoint>> sortedSPCollection);

      /**
       * Returns a vector connecting two spacepoints
       */
      TVector3 GetSegment(
          art::Ptr<recob::SpacePoint> a,
          art::Ptr<recob::SpacePoint> b);

      /**
       * Gets anle between two segments
       */
      float GetAngleBetweenSegments(
          TVector3 a,
          TVector3 b);
  };

}

#endif
