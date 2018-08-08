#ifndef GALLIUM_h
#define GALLIUM_h

#include "datasets.h"

class Gallium: public dataset{
  public:
    Gallium(){};

    int Init(std::string dataLoc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nPoints = 4;
    const double radiusGallex = 1.9;
    const double heightGallex = 5.0;
    const double sourceHeightGallex[2] = {2.7,2.38};
    const double radiusSage = 0.7;
    const double heightSage = 1.47;
    const double sourceHeightSage = 0.72;
    double Radius[3] = {radiusGallex, radiusGallex, radiusSage};
    double Height[3] = {heightGallex, heightGallex, heightSage};
    double SourceHeight[3] = {sourceHeightGallex[0], sourceHeightGallex[1], sourceHeightSage};

    double obsRatioGal[4], errorGal[4], arLinesE[2], arLinesBr[2], arLinesXSec[2], crLinesE[4], crLinesBr[4], crLinesXSec[4];
    std::vector < std::vector < double > > volInt;

    std::vector < double > FullData, Prediction;
};

#endif
