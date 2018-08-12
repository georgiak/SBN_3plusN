#ifndef XSEC_H
#define XSEC_H

#include "datasets.h"

class XSec: public dataset{
  public:
    XSec(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    TMinuit *gMinuit;
};

double getSinSqTerm(double dm2, double E1, double E2, double l1, double l2);
double osc_int(double E1, double E2, double l1, double l2);
void fcnXSec(int &npar, double *gin, double &fval, double *xval, int iflag);


#endif
