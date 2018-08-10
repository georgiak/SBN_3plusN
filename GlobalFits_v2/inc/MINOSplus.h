#ifndef MINOSPLUS_H
#define MINOSPLUS_H

#include "datasets.h"

class MINOSplus: public dataset{
  public:
    MINOSplus(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    TH1D* GetTwoDetSpectrum(TH1D* hND, TH1D* hFD);
};

#endif
