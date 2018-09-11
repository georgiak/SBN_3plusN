#ifndef NEOS_H
#define NEOS_H

#include "datasets.h"

class NEOS: public dataset{
  public:
    NEOS(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 60;
    const int lenAvg = 10;
    const double DB_norm_EH1 = 9.32293259423e+51;
    const double DB_norm_EH2 = 8.70923956974e+51;

    std::vector < double > DB_Predicted_3n, NEOS_Predicted_3n;
    std::vector < double > Energy, DB_eff, XSec, NEOS_f_iso;
    std::vector < std::vector < double > > Flux_iso, DB_f_iso, DB_l_d, DB_smearingmatrix, NEOS_smearingmatrix;
};

#endif
