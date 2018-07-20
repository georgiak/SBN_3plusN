#ifndef KARMEN_H
#define KARMEN_H

#include "datasets.h"

class KARMEN: public dataset{
  public:
    KARMEN(){};
    int Init(std::string dataLoc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 9;

    float _INTdm2, _INTnorm;          // variables for integrals
    double sinSqFunction(const double x);
    double sinSqFunctionCPV(const double x);
    double normFunc(const double x);
    };
}
