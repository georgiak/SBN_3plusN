#ifndef MINOS_H
#define MINOS_H

#include "datasets.h"

class MINOS: public dataset{
  public:
    MINOS(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    const int nBins = 12;
    const int nBins_ws = 13;
    std::vector < double > EnuQE, NumubarData, NumubarBkg, fracError, dataErr, EnuQE_ws, NumubarData_ws, NumubarBkg_ws, fracError_ws, dataErr_ws;
    std::vector < double > Signal, Prediction, TotalError;

};

#endif
