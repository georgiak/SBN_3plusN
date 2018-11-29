#ifndef MINIBOONE_COMBINED_H
#define MINIBOONE_COMBINED_H

#include "datasets.h"

class MiniBooNE_combined: public dataset{
  public:
    MiniBooNE_combined(){};
    using dataset::Init;
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    static const short nBins_e = 11;
    static const short nBins_mu = 8;
    const int nFOscEvts_nu = 17204;
    const int nFOscEvts_nubar = 117949;

    std::array < double, 2*nBins_e + 2*nBins_mu > FullData;
    std::array < double, nBins_e + nBins_mu > Background_nu, Background_nubar;
    std::array < std::array < double, 2*nBins_e + 2*nBins_e + 2*nBins_mu >, 2*nBins_e + 2*nBins_e + 2*nBins_mu > Full_fractCovMatrix;
    std::array < double, 2*nBins_e > Signal_nu, Signal_nubar, Signal_BestFit_nu, Signal_BestFit_nubar;
    std::array < std::array < double, nBins_e >, 100 > Lib_sin_nu, Lib_sinsq_nu, Lib_sin_nubar, Lib_sinsq_nubar;

    TMatrixT <float> covMatrix, full_covMatrix, cov;
};

#endif
