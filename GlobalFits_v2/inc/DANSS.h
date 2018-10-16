#ifndef DANSS_H
#define DANSS_H

#include "datasets.h"

class DANSS: public dataset{
  public:
    DANSS(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    static const int nBins = 24;
    const double L0_down = 12.85;
    const double L0_up = 10.9;
    const int binAvg = 50;
    const int lenAvg = 10;
    std::array < double, nBins > Observed, StatsError;
    std::array < std::array < double,48 >, 48 > EnergySmearing;
};

#endif
