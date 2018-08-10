#ifndef NOMAD_H
#define NOMAD_H

#include "datasets.h"

class NOMAD: public dataset{
  public:
    NOMAD(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int maxEnergyBins = 10; const int maxRadialBins = 3;
    const double l1 = .422;       const double l2 = .837;
    const double energyMin[10] = {3.,12.,16.,20.,25.,30.,35.,40.,50.,100.};
    const double energyMax[10] = {12.,16.,20.,25.,30.,35.,40.,50.,100.,170.};
    std::vector < std::vector < double > > sinSqDeltaGrid, sinSqDeltaGrid2;
    std::vector < double > Norm,  Signal, Observed, Bkg;
    TMatrix SigmaRemu;
};

struct integralFuncsNomad{
    double _dm2, _EnuAvg;

    double sinSqFunction(const double x);
    double sinSqFunctionCPV(const double x);
};

#endif
