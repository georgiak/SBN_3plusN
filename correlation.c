#include "correlation.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#define N_m_bins 19
#define N_e_bins 11
#define N_dets 3



void fake_fill(TMatrixT <double> &M){

	//Fills a square matrix of dim matrix_size with random numbers for now.
	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,std::time(0));

	int matrix_size=M.GetNrows();

	if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}
   	 
	for(int i=0; i<matrix_size; i++){
		for (int j = i;j<matrix_size;j++){
				M(i,j)=gsl_rng_uniform(r);
				M(j,i)=M(i,j);
		}
	
	}

gsl_rng_free(r);
 return ;
}




void stats_fill(TMatrixT <double> &M, std::vector<double> diag){
	int matrix_size = M.GetNrows();

	if(matrix_size != diag.size()){std::cout<<"#ERROR: stats_fill, matrix not equal to diagonal"<<std::endl;}
	if(M.GetNrows()!=M.GetNcols()){std::cout<<"#ERROR: not a square matrix!"<<std::endl;}


	M.Zero();
	

	for(int i=0; i<matrix_size; i++)
	{
		M(i,i) = diag[i];	

	}



 return ;
}

void contract_signal(TMatrixT <double> & M, TMatrixT <double> & Mc){

	// take the lower N_e_bins x N_m_bins matrix as a start
	//

	Mc=M.GetSub(N_dets*N_e_bins,M.GetRowUpb(),N_dets*N_e_bins,M.GetColUpb());	
//	Mc=M.GetSub(M.GetRowLwb(),N_dets*N_e_bins,M.GetColLwb(),N_dets*N_e_bins);	
	//Add top left down
	for(int i=0; i < N_dets*N_e_bins;i++){
		for(int j=0; j < N_dets*N_e_bins;j++){
			Mc(i,j)+= M(i,j);
		}
	}

	//add remaining top down
	for(int i=0; i < N_dets*(N_e_bins+N_m_bins) ;i++){
		for(int j=0; j < N_dets*N_e_bins;j++){
			Mc(i,j)+= M(i+N_dets*N_e_bins,j);
		}
	}
	//add remaining left across
	for(int i=0; i < N_dets*N_e_bins ;i++){
		for(int j=0; j < N_dets*(N_e_bins+N_m_bins);j++){
			Mc(i,j)+= M(i,j+N_dets*N_e_bins);
		}
	}



return;
}







std::vector<double > calc_signal_events(struct neutrinoModel &nuModel){
	
	std::vector<double > ans;

	

	
/*
for(int i=0; i<N_e_bins;i++){
	signal[i] =
	signal[i+N_e_bins]=
	signal[i+2*N_e_bins]=	
}*/

return ans;
}

