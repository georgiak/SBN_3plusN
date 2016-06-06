#ifndef CORRELATION_H_
#define CORRELATION_H_

#include <cmath>
#include <vector>
#include "TMatrixT.h"
#include <iostream>
#include "model.h"

void sys_fill(TMatrixT <double> & M, bool detsys);
void fake_fill(TMatrixT <double>&  M);
void stats_fill(TMatrixT <double>&  M, std::vector<double> diag);

std::vector<double >  calc_signal_events(struct neutrinoModel &nuModel);


void contract_signal(TMatrixT <double> & M, TMatrixT <double> &Mc);
void contract_signal2(TMatrixT <double> & M, TMatrixT <double> &Mc);

#endif
