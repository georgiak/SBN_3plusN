#ifndef NEOS_H
#define NEOS_H

#include "datasets.h"

class NEOS: public dataset{
  public:
    NEOS(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    static const int nBins = 60;
    const int binAvg = 10;
    const int lenAvg = 20;

    const double DB_norm_EH1 = 9.32293259423e+51;
    const double DB_norm_EH2 = 8.70923956974e+51;
    TMatrix Cov;

    std::array < double, nBins > Energy, Observed, StatsError, DB_Predicted_3n, NEOS_Predicted_3n;
    std::array < double, 200 > XSec;
    std::array < double, 4 > DB_eff, NEOS_f_iso, DB_mass;
    std::array < std::array < double, 200 >, 4 > Flux_iso;
    std::array < std::array < double, 4 >, 4 > DB_f_iso;
    std::array < std::array < double, 6 >, 4 > DB_l_d;
    std::array < std::array < double, 240 >, 240 > DB_smearingmatrix;
    std::array < std::array < double, 200 >, 200 > NEOS_smearingmatrix;
};

#endif
