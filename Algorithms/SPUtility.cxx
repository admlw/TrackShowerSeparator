#include "SPUtility.h"

namespace trackshowerseparator{

  void SPUtility::checkVectorForDuplicates(std::vector<art::Ptr<recob::SpacePoint>> spCollection){

    int duplicateCounter = 0;
    for (size_t i = 0; i < spCollection.size(); i++){

      art::Ptr<recob::SpacePoint> sp_1 = spCollection.at(i);
      double sp1_xyz[3] = {sp_1->XYZ()[0], sp_1->XYZ()[1], sp_1->XYZ()[2]};

      for (size_t j = i+1; j < spCollection.size(); j++){

        art::Ptr<recob::SpacePoint> sp_2 = spCollection.at(j);
        double sp2_xyz[3] = {sp_2->XYZ()[0], sp_2->XYZ()[1], sp_2->XYZ()[2]};

        if (sp1_xyz[0] == sp2_xyz[0] && sp1_xyz[1] == sp2_xyz[1] && sp1_xyz[2] == sp2_xyz[2])
          duplicateCounter++;

      }

    }

    std::cout << "Number of duplicates: " << duplicateCounter << std::endl;

  }

  std::vector<art::Ptr<recob::SpacePoint>> SPUtility::getSortedSPList(std::vector<art::Ptr<recob::SpacePoint>> spCollection, const recob::Vertex vertex, std::string sortType){

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

    if (sortType == "VertexDistance" || sortedSPCollection.size() == 0)
      return sortedSPCollection;
    else if (sortType == "NearestNeighbour"){

      std::vector<art::Ptr<recob::SpacePoint>> resortedSPCollection;

      // start from the first entry of the sorted vector: nearest the vertex
      art::Ptr<recob::SpacePoint> thisSP = sortedSPCollection.at(0);

      // loop the spacepoints, finding the next closest to the current spacepoint
      // and using that as the next spacepoint

      // vector of bools defines which spacepoints have already been used
      std::vector<bool> isUsed;
      isUsed.resize(sortedSPCollection.size(), false);
      isUsed.at(0) = true;

      int nearestSP = 0;

      for(size_t i_sp = 0; i_sp < sortedSPCollection.size(); i_sp++){

        float shortestDistance = std::numeric_limits<float>::max();
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
            isUsed.at(j_sp) = true;
          }

        }

        thisSP = sortedSPCollection.at(nearestSP);

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

  int SPUtility::FindNSPInShell(float shellWidth, int shellNumber, std::vector<art::Ptr<recob::SpacePoint>> sortedSPCollection){

    trackshowerseparator::SPUtility _sputility;

    // construct shell 
    float lowerShell = shellNumber*shellWidth;
    float upperShell = (shellNumber+1)*shellWidth;

    art::Ptr<recob::SpacePoint> vertexSP = sortedSPCollection.at(0);
    double vertexSP_xyz[3] = {vertexSP->XYZ()[0], vertexSP->XYZ()[1], vertexSP->XYZ()[2]};

    int nSPInShell = 1;
    for (size_t i = 1; i < sortedSPCollection.size(); i++){

      art::Ptr<recob::SpacePoint> thisSP = sortedSPCollection.at(i);
      double thisSP_xyz[3] = {thisSP->XYZ()[0], thisSP->XYZ()[1], thisSP->XYZ()[2]};

      // get 3D length
      float thisLength = _sputility.get3DLength(vertexSP_xyz, thisSP_xyz);

      if (thisLength > lowerShell && thisLength < upperShell){

        nSPInShell++;

      }

    }

    std::cout << "Number of SP in Shell " << shellNumber << " is " << nSPInShell << std::endl;

    return nSPInShell;

  }

  TVector3 SPUtility::GetSegment(art::Ptr<recob::SpacePoint> sp1, art::Ptr<recob::SpacePoint> sp2){

    double sp1_xyz[3] = {sp1->XYZ()[0], sp1->XYZ()[1], sp1->XYZ()[2]};
    double sp2_xyz[3] = {sp2->XYZ()[0], sp2->XYZ()[1], sp2->XYZ()[2]};

    TVector3 segmentVector(sp2_xyz[0]-sp1_xyz[0], sp2_xyz[1]-sp1_xyz[1], sp2_xyz[2]-sp1_xyz[1]);

    return segmentVector;

  }

  float SPUtility::GetAngleBetweenSegments(TVector3 a, TVector3 b){

    float angle = std::acos(a.Dot(b)/(a.Mag()*b.Mag()));

    return angle;

  }
}


