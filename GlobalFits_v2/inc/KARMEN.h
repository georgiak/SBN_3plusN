#ifndef KARMEN_H
#define KARMEN_H

#include "datasets.h"

class KARMEN: public dataset{
  public:
    KARMEN(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 9;
    double norm;
    std::vector < double > dm2Vec, Bkg, Observed;
    std::vector < std::vector < double > > sinSqDeltaGrid, sinSqDeltaGrid2;
};

struct integralFuncsKarmen{
    double _dm2, _norm;

    double sinSqFunction(const double x);
    double sinSqFunctionCPV(const double x);
    double normFunc(const double x);
};
#endif
