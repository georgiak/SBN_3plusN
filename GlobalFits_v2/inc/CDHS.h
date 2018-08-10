#ifndef CDHS_H
#define CDHS_H

#include "datasets.h"

class CDHS: public dataset{
  public:
    CDHS(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 15;
    std::vector < std::vector < double > > SinSqDeltaGrid_front, SinSqDeltaGrid_back, NoOscGrid;
    std::vector < double > Dm2Vec, Observed, M_front, M_back;
    TMatrix SigmaRatio;
};

#endif
