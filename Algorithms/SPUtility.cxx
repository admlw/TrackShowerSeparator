#include "SPUtility.h"

namespace trackshowerseparator{

  std::vector<art::Ptr<recob::SpacePoint>> SPUtility::getSortedSPList(std::vector<art::Ptr<recob::SpacePoint>> spCollection, const recob::Vertex vertex, std::string sortType){

    std::cout << "space point collection size is " << spCollection.size() << std::endl;

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

    std::cout << sortedSPCollection.size() << std::endl;

    if (sortType == "VertexDistance" || sortedSPCollection.size() == 0)
      return sortedSPCollection;
    else if (sortType == "NearestNeighbour"){

      std::vector<art::Ptr<recob::SpacePoint>> resortedSPCollection;

      // start from the first entry of the sorted vector: nearest the vertex
      art::Ptr<recob::SpacePoint> thisSP = sortedSPCollection.at(0);

      // loop the spacepoints, finding the next closest to the current spacepoint
      // and using that as the next spacepoint
      
      // vector of bools defines which spacepoints have already been used
      std::vector<bool> isUsed(sortedSPCollection.size(), false);
      isUsed.at(0) = true;
      
      int nearestSP = 0;
      float shortestDistance = std::numeric_limits<float>::max();

      for(size_t i_sp = 0; i_sp < sortedSPCollection.size(); i_sp++){
        
        // get thisSP position
        double thissp_xyz[3] = {thisSP->XYZ()[0], thisSP->XYZ()[1], thisSP->XYZ()[2]};

        for (size_t j_sp = 0; j_sp < sortedSPCollection.size(); j_sp++){

         art::Ptr<recob::SpacePoint> testSP = sortedSPCollection.at(j_sp); 

          // get testSP position
          double testsp_xyz[3] = {testSP->XYZ()[0], testSP->XYZ()[1], testSP->XYZ()[2]};

          // check if the distance between thisSP and the testSP is smaller than
          // any other differences calculated so far
          double testLength = _sputility.get3DLength(thissp_xyz, testsp_xyz);
          if (testLength < shortestDistance && isUsed.at(j_sp) == false){
            shortestDistance = testLength;
            nearestSP = j_sp;
          }

        }

        thisSP = sortedSPCollection.at(nearestSP);
        isUsed.at(nearestSP) = true;

        resortedSPCollection.push_back(thisSP);
 
      }

      return resortedSPCollection;

    }
    else{
      std::cout << "sortType can ONLY be \"VertexDistance\" or \"NearestNeighbour\"" << std::endl;
      throw std::exception();
    }

  }

  float SPUtility::get3DLength(double* xyz_1, double* xyz_2){

    float dist = std::sqrt(std::pow(xyz_1[0]-xyz_2[0],2) + 
    std::pow(xyz_1[1]-xyz_2[1],2) + 
    std::pow(xyz_1[2]-xyz_2[2],2));

    return dist;

  }

}
