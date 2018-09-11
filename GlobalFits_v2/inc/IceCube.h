#ifndef ICECUBE_H
#define ICECUBE_H

#include "datasets.h"

class IceCube: public dataset{
  public:
    IceCube(){};
    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
};

#endif
