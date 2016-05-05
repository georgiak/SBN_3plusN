#include "correlation.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <TFile.h>
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

//	std::cout<<"te "<<M.GetRowLwb()<<std::endl;

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


void sys_fill(TMatrixT <double> & Min)
{
	Min.Zero();

	TFile *fm= new TFile("rootfiles/covariance_matrices_351x351.root");
	Min = *(TMatrixT <float> *)fm->Get("TMatrixT<float>;7");

	fm->Close();
return;
}


void contract_signal2(TMatrixT <double> & M, TMatrixT <double> & Mc){

	// take the lower N_e_bins x N_m_bins matrix as a start
	//

	//these shouldnt be hardcoded
	int bblock = 117; 
	int eblock = 11;
	int mblock = 20;
	int ebnum = 7;	
	int mbnum = 2;

	int cblock = eblock+mblock; //size of each contracted matrix


	std::vector<std::vector< TMatrixT<double>  >> mtot(3, std::vector<TMatrixT<double> >(3));	
		
	for(int i = 0; i < 3; i++)
	{
	for(int j = 0; j < 3; j++)
	{
		mtot[i][j].ResizeTo(eblock+mblock,eblock+mblock);

		TMatrixT <double > Mtemp(bblock,bblock);
		Mtemp = M.GetSub(i*bblock,i*bblock+bblock-1,j*bblock,j*bblock+bblock-1);

		//std::cout<<"Gotten first subblock"<<std::endl;

		std::vector< TMatrixT<double> > ve(ebnum);
		std::vector< TMatrixT<double> > vm(mbnum);

		for(int k = 0; k<ebnum; k++)
		{
		//	std::cout<<k<<" starting!"<<std::endl;
			ve[k].ResizeTo(0,eblock-1,0,bblock-1);	
		//	std::cout<<k<<" resized!"<<std::endl;
			ve[k] = Mtemp.GetSub(k*eblock,k*eblock+eblock-1,0,bblock-1);
		//	std::cout<<k<<" done!"<<std::endl;
		}

		
		for(int k =1; k<ebnum; k++)
		{
			ve[0]=ve[0]+ve[k];		

		}
		//now ve[0] is the horiz collapsed erow



	//	std::cout<<"Gotten ve[0]"<<std::endl;

		for(int k = 0; k<mbnum; k++)
		{
			vm[k].ResizeTo(0,mblock-1,0,bblock-1);

			vm[k]=Mtemp.GetSub(ebnum*eblock+k*mblock, ebnum*eblock + k*mblock+mblock-1,  0,  bblock-1); 

		}
		for(int k =1; k<mbnum; k++)
		{
			vm[0] =vm[0]+ vm[k];		

		}
		//now vm[0] is the horiz collapsed murow


	//	std::cout<<"Gotten vm[0]"<<std::endl;


		std::vector< TMatrixT<double> > vev(ebnum);
		std::vector< TMatrixT<double> > vmv(ebnum);



		for(int k = 0; k<ebnum; k++)
		{
			vev[k].ResizeTo(0,eblock-1,0,eblock-1);
			vmv[k].ResizeTo(0,mblock-1,0,eblock-1);

			vev[k]=ve[0].GetSub(0,eblock-1  , k*eblock , k*eblock+eblock-1);
			vmv[k]=vm[0].GetSub(0,mblock-1,   k*eblock , k*eblock+eblock-1);

		}
		for(int k =1; k<ebnum; k++)
		{
			vev[0] = vev[0]+vev[k];		
			vmv[0] = vmv[0]+ vmv[k];		
		}



	//	std::cout<<"Gotten vev and vmv[0]"<<std::endl;

		std::vector< TMatrixT<double> > vev2(mbnum);
		std::vector< TMatrixT<double> > vmv2(mbnum);


		for(int k = 0; k<mbnum; k++)
		{
			vev2[k].ResizeTo(0,eblock-1,0,mblock-1);
			vmv2[k].ResizeTo(0,mblock-1,0,mblock-1);
				
			vev2[k]=ve[0].GetSub(0,eblock-1  , ebnum*eblock+k*mblock  ,  ebnum*eblock+k*mblock+mblock-1);
			vmv2[k]=vm[0].GetSub(0,mblock-1  , ebnum*eblock+k*mblock  ,  ebnum*eblock+k*mblock+mblock-1);

		}
		for(int k =1; k<mbnum; k++)
		{
			vev2[0]=vev2[0]+vev2[k];		
			vmv2[0]=vmv2[0]+vmv2[k];		
		}

	//	std::cout<<"gottena all"<<std::endl;

	/*	std::cout<<"  rows : cols"<<std::endl;
		std::cout<<"tl "<<vev[0].GetNrows()<<" "<<vev[0].GetNcols()<<std::endl;
		std::cout<<"bl "<<vmv[0].GetNrows()<<" "<<vmv[0].GetNcols()<<std::endl;
		std::cout<<"tr "<<vev2[0].GetNrows()<<" "<<vev2[0].GetNcols()<<std::endl;
		std::cout<<"br "<<vmv2[0].GetNrows()<<" "<<vmv2[0].GetNcols()<<std::endl;
	*/
		TMatrixT<double> Mdone(eblock+mblock,eblock+mblock);
		Mdone.Zero();


		Mdone.SetSub(0,0,vev[0]);
		Mdone.SetSub(eblock,eblock,vmv2[0]);
		Mdone.SetSub(0,eblock,vev2[0]);
		Mdone.SetSub(eblock,0,vmv[0]);

		mtot[i][j]=Mdone;
	}//j for
	}//i for


	for(int i =0; i<3; i++){
		for(int k =0; k<3; k++){
			Mc.SetSub(cblock*i,cblock*k,mtot[i][k]);	
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

