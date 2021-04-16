#ifndef MICROBOONE_DIS_H
#define MICROBOONE_DIS_H

#include "datasets.h"

class MicroBooNE_dis: public dataset{
  public:
    MicroBooNE_dis(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    bool shapeonly, signalInject;
    const int nBins = 20;
    const int nMC = 6039;
    int dm2_precalc_density;

    std::vector < float > FullData, Background, Libdis_noosc,MCStatSquared_noosc;
    std::vector < std::vector < float > > Full_fractCovMatrix, Libdis_sinsq,MCStatSquared_lib;

    TMatrixT <float> covMatrix, cov;
};

class MicroBooNE_dis_2d: public dataset{
  public:
    MicroBooNE_dis_2d(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    bool shapeonly, signalInject;
    int nBins;
    const int nMC = 6039;
    int dm2_precalc_density;

    std::vector < float > FullData, Background, Libdis_noosc,MCStatSquared_noosc;
    std::vector < std::vector < float > > Full_fractCovMatrix, Libdis_sinsq,MCStatSquared_lib;

    TMatrixT <float> covMatrix, cov;
};

#endif
