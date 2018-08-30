#include "SPUtility.h"

namespace trackshowerseparator{

  std::vector<art::Ptr<recob::SpacePoint>> SPUtility::getSortedSPList(std::vector<art::Ptr<recob::SpacePoint>> spCollection, const recob::Vertex vertex){

    trackshowerseparator::SPUtility _sputility;

    // get vertex reconstructed xyz
    double vertex_xyz[3];
    vertex.XYZ(vertex_xyz);

    // for each spacepoint, calculate distance from the vertex
    // and store information with SpacePoint in a vector of pairs
    std::vector<std::pair<art::Ptr<recob::SpacePoint>, float>> sp_distance_pairs;
    for (size_t i_sp = 0; i_sp < spCollection.size(); i_sp++){

      art::Ptr<recob::SpacePoint> thisSP = spCollection.at(i_sp);

      double sp_xyz[3] = {thisSP->XYZ()[0], thisSP->XYZ()[1], thisSP->XYZ()[2]};
      float spvLength = _sputility.get3DLength(vertex_xyz, sp_xyz);

      std::pair<art::Ptr<recob::SpacePoint>, float> sp_distance_pair(thisSP, spvLength);

      sp_distance_pairs.push_back(sp_distance_pair);
    }

    // sort vector of pairs by the distance from the vertex 
    std::sort(sp_distance_pairs.begin(), sp_distance_pairs.end(), [](auto &left, auto &right) {
      return left.second < right.second;
    });
    
    // and push into new vector to return
    std::vector<art::Ptr<recob::SpacePoint>> sortedSPCollection;
    for (size_t i_pair = 0; i_pair < sp_distance_pairs.size(); i_pair++){

      sortedSPCollection.push_back(sp_distance_pairs.at(i_pair).first);

    }

    return sortedSPCollection;

  }

  float SPUtility::get3DLength(double* xyz_1, double* xyz_2){

    float dist = std::sqrt(std::pow(xyz_1[0]-xyz_2[0],2) + 
    std::pow(xyz_1[1]-xyz_2[1],2) + 
    std::pow(xyz_1[2]-xyz_2[2],2));

    return dist;

  }

}
