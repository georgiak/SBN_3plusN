#ifndef CCFR_H
#define CCFR_H

#include "datasets.h"

class CCFR: public dataset{
  public:
    CCFR(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 18;
    const double m_front = 105.;
    const double m_back = 444.;
    double dm2Vec[dm2VecMaxDim];
    std::vector < double > Observed;
    std::vector < std::vector < double > > sinSqDeltaGrid_front, sinSqDeltaGrid_back, noOscGrid;
    TMatrix SigmaRatio;
};

#endif
