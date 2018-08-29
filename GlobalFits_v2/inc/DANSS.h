#ifndef DANSS_H
#define DANSS_H

#include "datasets.h"

class DANSS: public dataset{
  public:
    DANSS(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 24;
    const double L0_down = 12.85;
    const double L0_up = 10.9;
    std::vector < double > Observed, Energy, StatsError;
    std::vector < std::vector < double > > EnergySmearing;
};

#endif
