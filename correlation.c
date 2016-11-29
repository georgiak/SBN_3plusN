#include "correlation.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <TFile.h>


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


std::vector<std::vector<double >> to_vector(TMatrixT <double > Min)
{
	int dimension =  Min.GetNrows();

	std::vector<std::vector<double >>  ans(dimension, std::vector<double>(dimension));

	for(int i = 0; i< dimension; i++){
		for(int k = 0; k< dimension; k++){
			ans[i][k]=Min(i,k);
			if(ans[i][k]==-0){
				ans[i][k]=0;
			}
		}	
	}
return ans;


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


TMatrixT<double > sys_fill_direct(int dim, bool detsys){
	TMatrixT<float > * temp;
	TMatrixT<double > temp2(dim,dim);
	//std::cout<<"inputted dim is: "<<dim<<std::endl;
	if(dim ==   (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti  ){
		TFile *fm= new TFile("rootfiles/covariance_matrices_690x690.root");
		 temp = (TMatrixT <float>* )fm->Get("TMatrixT<float>;7");
		//std::cout<<"outputted temp dim is: "<<temp.GetNcols()<<std::endl;
		for(int i =0; i<dim; i++)
		{
			for(int j =0; j<dim; j++)
			{
		//		std::cout<<i<<" "<<j<<" "<<(*temp)(i,j)<<std::endl;
				temp2(i,j)=(*temp)(i,j);
			}
		}
		delete temp;




	//	delete fm;
		fm->Close();
		delete fm;
		return temp2;
	} else {
	if(detsys){
		TFile *fm= new TFile("rootfiles/covariance_matrices_345x345.root");
		 temp2 = *(TMatrixT <float> *)fm->Get("TMatrixT<float>;7");
	//	delete fm
//		for(int i =0; i<dim; i++)
//		{
//			for(int j =0; j<dim; j++)
//			{
//				std::cout<<i<<" "<<j<<" "<<temp(i,j)<<std::endl;
//			}
//		}


	
		fm->Close();
		return temp2;
	} else {
		std::cout<<"ERROR: This probably shouldnt run, no-detsys"<<std::endl;
		TFile *fm= new TFile("rootfiles/covariance_matrices_nodetsys_345x345.root");
		 temp2 = *(TMatrixT <float> *)fm->Get("TMatrixT<float>;7");
	//	delete fm;
		fm->Close();
		return temp2;
	}	


	}

}

void sys_fill2(int dim, TMatrixT<float> * ans)
{


	if(dim==   (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti  ){
		TFile *fm= new TFile("rootfiles/covariance_matrices_690x690.root");
		TMatrixT<float> * temp = new TMatrixT<float>(dim,dim);
		temp = (TMatrixT <float> *)fm->Get("TMatrixT<float>;7");
		ans =temp;

		fm->Close();
		delete fm;
	} else{

		TFile *fm= new TFile("rootfiles/covariance_matrices_345x345.root");
		TMatrixT <float>* temp = new TMatrixT<float>(dim,dim);
		temp= (TMatrixT <float> *)fm->Get("TMatrixT<float>;7");
		ans = temp;
		
		fm->Close();
		delete fm;
	}	

return;
}



void sys_fill(TMatrixT <double> & Min, bool detsys)
{

	Min.Zero();

	if(Min.GetNrows()==   (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti  ){
		TFile *fm= new TFile("rootfiles/covariance_matrices_690x690.root");
		Min = *(TMatrixT <float> *)fm->Get("TMatrixT<float>;7");
	//	delete fm;
		fm->Close();
		delete fm;
	} else {
	if(detsys){

		TFile *fm= new TFile("rootfiles/covariance_matrices_345x345.root");
		Min = *(TMatrixT <float> *)fm->Get("TMatrixT<float>;7");
	//	delete fm;
		fm->Close();
	} else {
		TFile *fm= new TFile("rootfiles/covariance_matrices_nodetsys_345x345.root");
		Min = *(TMatrixT <float> *)fm->Get("TMatrixT<float>;7");
	//	delete fm;
		fm->Close();
	}	


	}
return;
}




void contract_signal(TMatrixT <double> & M, TMatrixT <double> & Mc){
	/******************************************************
	 * REDUNDANT/OBSOLETE see contracr_signal2() below
	 *****************************************************/

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


void contract_signal2(TMatrixT <double> & M, TMatrixT <double> & Mc){
	// So this function takes a Matrix, which has N_detectors multiples of ( ebnum sections of eblock bins follows mbnum sections of mblock) .
	// Contacts it down to (eblock+mblock)*N_detectors.
	// Basically is a much more generic version of contract_signal() above


	// take the lower N_e_bins x N_m_bins matrix as a start
	//these shouldnt be hardcoded
	int eblock = N_e_bins;
	int mblock = N_m_bins;
	
	int ebnum = N_e_spectra;	
	int mbnum = N_m_spectra;

	int bblock = eblock*ebnum+mblock*mbnum;//big block 
	int cblock = eblock+mblock; //size of each contracted matrix


	std::vector<std::vector< TMatrixT<double>  >> mtot(N_dets, std::vector<TMatrixT<double> >(3));	
		
	for(int i = 0; i < N_dets; i++)
	{
	for(int j = 0; j < N_dets; j++)
	{
		mtot[i][j].ResizeTo(eblock+mblock,eblock+mblock);

		TMatrixT <double > Mtemp(bblock,bblock);
		//First get the one of 9 subblocks to contract
		Mtemp = M.GetSub(i*bblock,i*bblock+bblock-1,j*bblock,j*bblock+bblock-1);

		//std::cout<<"Gotten first subblock"<<std::endl;


		std::vector< TMatrixT<double> > ve(ebnum);// will keep the electronlike submatricies (whole horizontal rows)
		std::vector< TMatrixT<double> > vm(mbnum);//will keep the muonlike submatricies

		//So for each of the (7) electron like spectra
		for(int k = 0; k<ebnum; k++)
		{
			// resize it to a eblock deep and whole bblock wide thing
			ve[k].ResizeTo(0,eblock-1,0,bblock-1);	
			//And get the lots of horizontal rows 
			ve[k] = Mtemp.GetSub(k*eblock,k*eblock+eblock-1,0,bblock-1);
		//	std::cout<<k<<" done!"<<std::endl;
		}

		// Going to use the ve[0] one to store the collapsed ones
		for(int k =1; k<ebnum; k++)
		{
			ve[0]=ve[0]+ve[k];		

		}
		//now ve[0] is the horiz collapsed erow



		// Do the same for muon row
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


	for(int i =0; i<N_dets; i++){
		for(int k =0; k<N_dets; k++){
			Mc.SetSub(cblock*i,cblock*k,mtot[i][k]);	
		}
	}
	

	

return;
}


void contract_signal2_anti(TMatrixT <double> & M, TMatrixT <double> & Mc){
	//so basically takes the 4 quadrants of (N_detectors multiples of ( ebnum sections of eblock bins follows mbnum sections of mblock)) 
	// and runs each one through the contract_signal2() function above. 


	int bblock = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;		//big block 
	int cblock = (N_e_bins+N_m_bins)*N_dets; 			//size of each contracted block matrix

	int antibblock =bblock*N_anti;
	int anticblock = cblock*N_anti;
		// tr is top right, bl is bottem left

	//std::cout<<M.GetNRows()<<" "<<M.GetNCols()<<" "<<Mc.GetNRows()<<" "<<Mc.GetNCols()<<std::endl;
	//std::cout<<"usual : "<<bblock<<" "<<cblock<<" "<<" anti: "<<antibblock<<" "<<anticblock<<std::endl;

	//Ok make four sub-matricies
	TMatrixT <double > Mnu(bblock,bblock); 
	TMatrixT <double > MnuBar(bblock,bblock); 	
	TMatrixT <double > Mtr(bblock,bblock); 	
	TMatrixT <double > Mbl(bblock,bblock); 	

	//And select them	
	Mnu = 	 M.GetSub(0,bblock-1,0,bblock-1); //checked
	MnuBar = M.GetSub(bblock,antibblock-1,bblock,antibblock-1);// checked
	Mtr =    M.GetSub(0,bblock-1,bblock,antibblock-1);// checked
	Mbl =    M.GetSub(bblock,antibblock-1,0,bblock-1);//checked

	TMatrixT <double > MnuC(cblock,cblock); 
	TMatrixT <double > MnuBarC(cblock,cblock); 	
	TMatrixT <double > MtrC(cblock,cblock); 	
	TMatrixT <double > MblC(cblock,cblock); 

	MnuC.Zero();
	MnuBarC.Zero();
	MtrC.Zero();
	MblC.Zero();

	contract_signal_layer2(Mnu,MnuC);	
	contract_signal_layer2(MnuBar,MnuBarC);	
	contract_signal_layer2(Mtr,MtrC);	
	contract_signal_layer2(Mbl,MblC);

//	Mc.Zero();
	Mc.ResizeTo(2*cblock,2*cblock); 

	Mc.SetSub(0,0,MnuC); //checked
	Mc.SetSub(cblock,cblock,MnuBarC); //checked
	Mc.SetSub(cblock,0,MtrC); //chdecked  //Ok the source of the error could be the switching of these two 
	Mc.SetSub(0,cblock,MblC);//checkded


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


void contract_signal_layer1(TMatrixT <double> & M, TMatrixT <double> & Mc){


//So given just a single detector, collapses.
/*#define N_m_bins 19
#define N_e_bins 11

#define N_e_spectra 7
#define N_m_spectra 2

#define N_dets 3
#define N_anti 2
*/
bool debug = false;

	if(debug)	std::cout<<"Starting:M "<<M.GetNcols()<<" "<<M.GetNrows()<<" "<<115<<std::endl;
	if(debug)	std::cout<<"Starting:Mc "<<Mc.GetNcols()<<" "<<Mc.GetNrows()<<" "<<30<<std::endl;

		if(debug) std::cout<<"Starting elike"<<std::endl;
		std::vector< TMatrixT<double> > elike;
		TMatrixT<double> elikeSummed(N_e_bins,N_e_bins);
		elikeSummed=0.0;


		if(debug)std::cout<<"Starting elike getting"<<std::endl;
		for(int m=0; m <= N_e_spectra-1; m++){
			for(int n=0; n<= N_e_spectra-1; n++){
			elike.push_back(M.GetSub(n*N_e_bins,n*N_e_bins+N_e_bins-1,m*N_e_bins, m*N_e_bins+N_e_bins-1 ));
			}
		}
		
		if(debug)std::cout<<"Starting elike summing"<<std::endl;
		for(int i =0; i<elike.size(); i++){
			if(debug) std::cout<<elike[i].GetNcols()<<" "<<elike[i].GetNrows()<<std::endl;
			elikeSummed+=elike[i];
		}
//*******************************************************
		
		if(debug)std::cout<<"Starting mlike"<<std::endl;
		std::vector< TMatrixT<double> > mlike;
		TMatrixT<double> mlikeSummed(N_m_bins,N_m_bins);
		mlikeSummed=0.0;
		int mstart = N_e_bins*N_e_spectra;

		if(debug)std::cout<<"Starting mlike getting"<<std::endl;
		for(int m=0; m <= N_m_spectra-1; m++){
			for(int n=0; n<= N_m_spectra-1; n++){
			mlike.push_back(M.GetSub(mstart+n*N_m_bins,mstart + n*N_m_bins+N_m_bins-1,mstart + m*N_m_bins,mstart+ m*N_m_bins+N_m_bins-1 ));
			}
		}
		
		if(debug)std::cout<<"Starting mlike summing"<<std::endl;
		for(int i =0; i<mlike.size(); i++){
			mlikeSummed+=mlike[i];
		}


//*******************************************************

		if(debug) std::cout<<"Starting emlike"<<std::endl;
		std::vector< TMatrixT<double> > emlike;
		TMatrixT<double> emlikeSummed(N_m_bins,N_e_bins);
		emlikeSummed =0.0;

		if(debug) std::cout<<"Starting emlike getting"<<std::endl;
		for(int m=0; m <= N_e_spectra-1; m++){
			for(int n=0; n<= N_m_spectra-1; n++){
			emlike.push_back(M.GetSub(mstart+n*N_m_bins,mstart + n*N_m_bins+N_m_bins-1, m*N_e_bins, m*N_e_bins+N_e_bins-1 ));
			}
		}
		
		if(debug) std::cout<<"Starting emlike summing"<<std::endl;
		for(int i =0; i<emlike.size(); i++){
			emlikeSummed+=emlike[i];
		}

//*******************************************************
//
		if(debug) std::cout<<"Starting melike"<<std::endl;
		std::vector< TMatrixT<double> > melike;
		TMatrixT<double> melikeSummed(N_e_bins,N_m_bins);
		melikeSummed =0.0;

		if(debug) std::cout<<"Starting melike getting"<<std::endl;
		for(int m=0; m <= N_m_spectra-1; m++){
			for(int n=0; n<= N_e_spectra-1; n++){
			melike.push_back(M.GetSub(n*N_e_bins, n*N_e_bins+N_e_bins-1,mstart+ m*N_m_bins,mstart+ m*N_m_bins+N_m_bins-1 ));
			}
		}
		
		if(debug) std::cout<<"Starting melike summing"<<std::endl;
		for(int i =0; i<melike.size(); i++){
			melikeSummed+=melike[i];
		}


	Mc.Zero(); 
	Mc.SetSub(0,0,elikeSummed); //checked
	Mc.SetSub(N_e_bins,N_e_bins,mlikeSummed); //checked
	Mc.SetSub(0,N_e_bins,melikeSummed); //chdecked  // should be me, 100 checked of mathematica
	Mc.SetSub(N_e_bins,0,emlikeSummed);//checkded


return;
}


void contract_signal_layer2(TMatrixT <double> & M, TMatrixT <double> & Mc){
		Mc.Zero();
		int nrow = N_e_bins*N_e_spectra+N_m_bins*N_m_spectra;
		int crow=N_e_bins+N_m_bins;
//		std::vector< TMatrixT<double> > mats;
//		std::vector< TMatrixT<double> > matsC(9);
	
		for(int m =0; m<= N_dets-1; m++){
			for(int n =0; n<= N_dets-1; n++){
				//mats.push_back(M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1));
				TMatrixT<double> imat(nrow,nrow);
				TMatrixT<double> imatc(crow,crow);
				imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
				contract_signal_layer1(imat,imatc);
				Mc.SetSub(n*crow,m*crow,imatc);

			}
		}

//		for(int i=0;i<mats.size(); i++){
			//std::cout<<mats[i].GetNcols()<<" "<<mats[i].GetNrows()<<std::endl;
//			matsC[i].Zero();
//			matsC[i].ResizeTo(N_e_bins+N_m_bins,N_e_bins+N_m_bins);
//			contract_signal_layer1(mats[i],matsC[i]);
	//	}

/*
		Mc.Zero();
		
		Mc.SetSub(0,0,matsC[0]);	
		Mc.SetSub(crow,crow,matsC[4]);	
		Mc.SetSub(2*crow,2*crow,matsC[8]);	

		// 3 6 7 1 2 5

		//This will take a while.
		Mc.SetSub(0,crow,matsC[3]);	
		Mc.SetSub(0,2*crow,matsC[6]);	
		Mc.SetSub(crow,2*crow,matsC[7]);	

		Mc.SetSub(crow,0,matsC[1]);	
		Mc.SetSub(2*crow,0,matsC[2]);	
		Mc.SetSub(2*crow,crow,matsC[5]);	
*/


return;
}
	
void contract_signal_layer3(TMatrixT <double> & M, TMatrixT <double> & Mc){
		Mc.Zero();
		int nrow = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		int crow=(N_e_bins+N_m_bins)*N_dets;
//		std::vector< TMatrixT<double> > mats;
//		std::vector< TMatrixT<double> > matsC(9);
	
		for(int m =0; m<= N_anti-1; m++){
			for(int n =0; n<= N_anti-1; n++){
				//mats.push_back(M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1));
				TMatrixT<double> imat(nrow,nrow);
				TMatrixT<double> imatc(crow,crow);
				imat = M.GetSub(n*nrow,n*nrow+nrow-1, m*nrow,m*nrow+nrow-1);
				contract_signal_layer2(imat,imatc);
				Mc.SetSub(n*crow,m*crow,imatc);

			}
		}

//		for(int i=0;i<mats.size(); i++){
			//std::cout<<mats[i].GetNcols()<<" "<<mats[i].GetNrows()<<std::endl;
//			matsC[i].Zero();
//			matsC[i].ResizeTo(N_e_bins+N_m_bins,N_e_bins+N_m_bins);
//			contract_signal_layer1(mats[i],matsC[i]);
	//	}

/*
		Mc.Zero();
		
		Mc.SetSub(0,0,matsC[0]);	
		Mc.SetSub(crow,crow,matsC[4]);	
		Mc.SetSub(2*crow,2*crow,matsC[8]);	

		// 3 6 7 1 2 5

		//This will take a while.
		Mc.SetSub(0,crow,matsC[3]);	
		Mc.SetSub(0,2*crow,matsC[6]);	
		Mc.SetSub(crow,2*crow,matsC[7]);	

		Mc.SetSub(crow,0,matsC[1]);	
		Mc.SetSub(2*crow,0,matsC[2]);	
		Mc.SetSub(2*crow,crow,matsC[5]);	
*/


return;
}
	
