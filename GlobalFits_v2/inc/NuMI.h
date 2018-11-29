#ifndef NUMI_H
#define NUMI_H

#include "datasets.h"

class NuMI: public dataset{
  public:
    NuMI(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    static const int nBins = 10;
    const int nFOscEvts = 3323;

    std::array < int, nBins > NueData;
    std::array < float, nBins > NueBgr, NueBgr_error;
    std::vector < float > Signal, TotalError, FOsc_fracError;

    std::vector < std::vector < float > >  Lib_sin, Lib_sinsq;
};
#endif
