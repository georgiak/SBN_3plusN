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

//18.68235, 23.41176, 27.36471, 30.54118, 33.78824, 37.31765, 40.42353, 43.45882, 47.05882, 48.18824, 50.94118, 52.91765, 54.82353, 55.31765, 57.08235, 57.29412, 55.74118, 56.30588, 55.88235, 55.31765, 54.96471, 54.68235, 52.56471, 50.09412, 48.32941, 46.21176, 44.23529, 43.24706, 41.27059, 38.65882, 37.60000, 36.32941, 35.05882, 33.78824, 32.30588, 30.47059, 29.34118, 27.50588, 25.45882, 24.32941, 21.78824, 20.37647, 18.40000, 16.91765, 15.01176, 13.88235, 12.40000, 11.05882, 10.70588, 9.36471, 7.88235, 7.17647, 6.18824, 5.55294, 4.77647, 4.00000, 3.36471, 2.87059, 2.23529, 1.88235
