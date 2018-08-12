#ifndef BUGEY_H
#define BUGEY_H

#include "datasets.h"

class Bugey: public dataset{
  public:
    Bugey(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int maxEnergyBins = 25;
    const int nBaselines = 3;
};

struct integralFuncsBugey{
  int _jB;
  double _dm2, _energy;

  double sinSqFunction(const double x);
};

void fcnBugey(int &npar, double *gin, double &fval, double  *xval, int iflag);

#endif
