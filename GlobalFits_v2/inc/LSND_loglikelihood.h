#ifndef LSND_LOGLIKELIHOOD_H
#define LSND_LOGLIKELIHOOD_H

#include "datasets.h"

class LSND_loglikelihood: public dataset{
  public:
    LSND_loglikelihood(){};

    int Init(std::string dataLoc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 5;
    double norm;
    std::vector < double > dm2Vec, Bkg, Observed;
    std::vector < std::vector < double > > sinSqDeltaGrid, sinSqDeltaGrid2;
};

#endif
