#ifndef MINIBOONE_H
#define MINIBOONE_H

#include "datasets.h"

class MiniBooNE: public dataset{
  public:
    MiniBooNE(bool _nubar){
      nubar = _nubar;
    }
    using dataset::Init;
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    bool nubar;
    const short nBins_e = 11;
    const short nBins_mu = 8;
    const int nFOscEvts = 17204;

    std::vector < float > Background, Signal, FullData, Signal_BestFit;
    std::vector < std::vector < float > > Full_fractCovMatrix, Lib_sin, Lib_sinsq;

    TMatrixT <float> covMatrix, full_covMatrix, cov;
};

#endif
