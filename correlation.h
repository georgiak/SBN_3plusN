#ifndef CORRELATION_H_
#define CORRELATION_H_

#include <cmath>
#include <vector>
#include "TMatrixT.h"
#include <iostream>



void fake_fill(TMatrixT <double>&  M);

std::vector<double >  calc_signal_events(struct neutrinoModel &nuModel);

void contract_signal(TMatrixT <double> & M, TMatrixT <double> &Mc);


#endif
