#ifndef FROMICSURF_H
#define FROMICSURF_H

#include "datasets.h"
#include "TH3D.h"

class IceCube_fromSurfs: public dataset{
  public:
    IceCube_fromSurfs(){};

    int Init(std::string dataLoc, Oscillator osc, bool debug);
    float Chi2(Oscillator osc, neutrinoModel nu, bool debug);

  private:
    std::vector<TH3D*> vecLogICHisto;
};

#endif
