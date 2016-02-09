#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <cstring>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph2D.h"
#include "TSystem.h"
#include "TMatrixT.h"
//#include "background.h"
#include "model.h"
#include "correlation.h"

#define N_m_bins 19
#define N_e_bins 11
#define N_dets 3




//double SBN_chi2(

/* ########## Main function ############### */
int main(int argc, char* argv[])
{

gSystem->Load("libTree");

bool fit_flag = false;
bool hist_flag = false;
bool verbose_flag = false;

int c;
opterr = 0;

while ((c = getopt(argc, argv, "v-:")) != -1)
{
switch(c)
{
case 'v':
	verbose_flag=true;
	break;
case '-':
	if(!strcmp(optarg,"fit"))
	{ 
		fit_flag = true;
	}
	else if(!strcmp(optarg,"hist"))
	{ 	
		hist_flag =  true;
	}
	else
	{
		std::cout<<"You passed a -- flag which I didn't understand."<<std::endl;
		return -1;
	}	
	break;
case '?':
//	std::cout<<"Abandon hope all ye who enter this value. "<<std::endl;
	std::cout<<"Allowed arguments:"<<std::endl;
	std::cout<<"\t-v\t\tVerbose run, mostly debugging"<<std::endl;	
	std::cout<<"\t--fit\t\trun SBN fitting code"<<std::endl;
	std::cout<<"\t--hist\t\truns to produce root histograms"<<std::endl;
	return 0;
default:
	std::cout<<"I don't know how you got here."<<std::endl;
	return -1;
}}



if(fit_flag){


	
SBN_spectrum myspec;
//myspec.THprint();
myspec.fill_vectors();
std::vector<double > back= myspec.get_sixvector();

neutrinoModel testModel;
SBN_spectrum newspec(testModel);
newspec.fill_vectors();
std::vector<double > pred= newspec.get_sixvector();




int matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
int matrix_size_c = (N_e_bins + N_m_bins) * N_dets;

TMatrixT <double> M(matrix_size,matrix_size);
TMatrixT <double> Mc(matrix_size_c,matrix_size_c);
TMatrixT <double> McI(matrix_size_c, matrix_size_c);

fake_fill(M);
//M.Print();


//std::cout<<M.GetNrows()<<" "<<M.GetNcols()<<std::endl;

double invdet=0;
//MI=M.Invert(&invdet);

contract_signal(M,Mc);
//Mc.Print();
//std::cout<<"Det: "<<M.Determinant()<<" Det(Inv): "<<invdet<<std::endl;
  
double res=0;
McI=Mc.Invert(&invdet);

//known bug!
if(false&&matrix_size_c != pred.size()&& matrix_size_c != back.size()){
	std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
	std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<pred.size()<<" back: "<<back.size()<<std::endl;	
}

for(int i =0; i<matrix_size_c; i++){
	for(int j =0; j<matrix_size_c; j++){
		res += (back[i]-pred[i])*McI(i,j)*(back[j]-pred[j]);

	}
}

if(verbose_flag){
	std::cout<<"****************** Background prediction ************"<<std::endl;
	myspec.THprint();
	std::cout<<"****************** Oscillation prediction ************"<<std::endl;
	newspec.THprint();
}
std::cout<<res<<std::endl;

//for(int i =0;i<wonder.size(); i++){
//	std::cout<<i<<"  "<<wonder[i]<<std::endl;
//}

//TNtuple *chi2Nt = new TNtuple("chi2Nt","chi2Nt","chi2:dof:gof:m4:ue4:um4:m5:ue5:um5:m6:ue6:um6:phi45:phi46:phi56");

}

if(hist_flag){
//Nothin' here yet!


}




return 0;
}





