// file1 information
bool file1_isData = false;
std::string file1_label = "Simulated Nominal";

bool file2_isData = false;
std::string file2_label = "Simulated DIC";

struct vars{
  bool pfpIsPandoraTrack = false;
  bool pfpIsPandoraShower =false;
  int spNSpacePoints = -1;
  float spAverageDistance = -1;
  float spMedianDistance = -1;
  float spRMSDistance = -1;
  std::vector<float> *clusterStartCharge = nullptr;
  std::vector<float> *clusterStartAngle = nullptr;
  std::vector<float> *clusterStartOpeningAngle = nullptr;
  std::vector<float> *clusterEndCharge = nullptr;
  std::vector<float> *clusterEndAngle = nullptr;
  std::vector<float> *clusterEndOpeningAngle = nullptr;
  std::vector<float> *clusterIntegral = nullptr;
  std::vector<float> *clusterIntegralStdDev = nullptr;
  std::vector<float> *clusterIntegralAverage = nullptr;
  std::vector<float> *clusterSummedADC = nullptr;
  std::vector<float> *clusterSummedADCStdDev = nullptr;
  std::vector<float> *clusterSummedADCAverage = nullptr;
  std::vector<float> *clusterMultipleHitDensity = nullptr;
  std::vector<float> *clusterWidth = nullptr;
  std::vector<float> *clusterPlane = nullptr;
  std::vector<int> *clusterNHits = nullptr;
};

