#include "prob.h"
#include <complex>
#include <vector>
#include <iostream>


/*********************************************
*  Constructors for neutrinoModel
* *******************************************/

/******	NULL constructor ******/
neutrinoModel::neutrinoModel(){
	zero();
	difference();
}

/******	3+3 constructor ******/
neutrinoModel::neutrinoModel(double * mn, double * ue, double * um, double * ph){
	zero();
	for(int i = 0; i < 3; i ++){
			mNu[i] = mn[i]; Ue[i] = ue[i];
			Um[i] = um[i];  phi[i] =ph [i];
	}
	difference();
}

/******	3+1 constructor ******/
neutrinoModel::neutrinoModel(double  mn, double  ue4, double  um4){
	zero();
	mNu[0] = mn;
	Ue[0]=ue4;
	Um[0]=um4;

	difference();
}




/*********************************************
*  	Other Models
* *******************************************/


void neutrinoModel::zero(){
	for(int i = 0; i < 3; i ++){
			mNu[i] = 0; Ue[i] = 0;
			Um[i] = 0;  phi[i] = 0;
	}
}

void neutrinoModel::difference(){
		dm41Sq = pow(mNu[0],2);
		dm51Sq = pow(mNu[1],2);
		dm61Sq = pow(mNu[2],2);
		dm54Sq = dm51Sq - dm41Sq;
		dm64Sq = dm61Sq - dm41Sq;
		dm65Sq = dm61Sq - dm51Sq;
	}


double neutrinoModel::oscProb(int a, int b, double Ev, double L){


	if(a == b)
	{
		return oscProb_dis(a, Ev, L);
	} 
	else
	{
		return oscProb_app(a,b, Ev, L);
	}
}

double neutrinoModel::oscProb_app(int a, int b, double Ev, double L){
	
	double nubarmod= 1.0;

	if(a < 0 && b < 0){
		nubarmod = -1.0;
	}

	double Ua4=0,Ua5=0,Ua6=0,Ub4=0,Ub5=0,Ub6=0;

	switch(a)
	{
		case 1:
			Ua4 = Ue[0];
			Ua5 = Ue[1];
			Ua6 = Ue[2];
			break;
		case 2:
			Ua4 = Um[0];
			Ua5 = Um[1];
			Ua6 = Um[2];
			break;
		case 3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}
	switch(b)
	{
		case 1:
			Ub4 = Ue[0];
			Ub5 = Ue[1];
			Ub6 = Ue[2];
			break;
		case 2:
			Ub4 = Um[0];
			Ub5 = Um[1];
			Ub6 = Um[2];
			break;
		case 3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}

	
	double phi54 = nubarmod*phi[0];
	double phi64 = nubarmod*phi[1];
	double phi65 = nubarmod*phi[2];
	
	
	double ans =0.0;

	ans =  -4.0*fabs(Ua5*Ub5*Ua4*Ub4)*cos(phi54)*pow(sin(1.27*dm54Sq*L/Ev),2.0);
	ans += -4.0*fabs(Ua6*Ub6*Ua4*Ub4)*cos(phi64)*pow(sin(1.27*dm64Sq*L/Ev),2.0);
	ans += -4.0*fabs(Ua5*Ub5*Ua6*Ub6)*cos(phi65)*pow(sin(1.27*dm65Sq*L/Ev),2.0);


	ans += 4.0*(fabs(Ua4*Ub4)+fabs(Ua5*Ub5)*cos(phi54)+fabs(Ua6*Ub6)*cos(phi64))*fabs(Ua4*Ub4)*pow(sin(1.27*dm41Sq*L/Ev),2.0);
	ans += 4.0*(fabs(Ua4*Ub4)*cos(phi54)+fabs(Ua5*Ub5)+fabs(Ua6*Ub6)*cos(phi65))*fabs(Ua5*Ub5)*pow(sin(1.27*dm51Sq*L/Ev),2.0);
	ans += 4.0*(fabs(Ua4*Ub4)*cos(phi64)+fabs(Ua5*Ub5)*cos(phi65)+fabs(Ua6*Ub6))*fabs(Ua6*Ub6)*pow(sin(1.27*dm61Sq*L/Ev),2.0);


	ans += (2.0*sin(phi54))*fabs(Ub5*Ua5*Ub4*Ua4)*sin(2.53*dm54Sq*L/Ev);
	ans += (2.0*sin(phi64))*fabs(Ub6*Ua6*Ub4*Ua4)*sin(2.53*dm64Sq*L/Ev);
	ans += (2.0*sin(phi65))*fabs(Ub6*Ua6*Ub5*Ua5)*sin(2.53*dm65Sq*L/Ev);


	ans += 2.0*(fabs(Ua5*Ub5)*sin(phi54)+fabs(Ua6*Ub6)*sin(phi64))*fabs(Ua4*Ub4)*sin(2.53*dm41Sq*L/Ev);
	ans += 2.0*(-fabs(Ua4*Ub4)*sin(phi54)+fabs(Ua6*Ub6)*sin(phi65))*fabs(Ua5*Ub5)*sin(2.53*dm51Sq*L/Ev);
	ans += 2.0*(-fabs(Ua4*Ub4)*sin(phi64)-fabs(Ua5*Ub5)*sin(phi65))*fabs(Ua6*Ub6)*sin(2.53*dm61Sq*L/Ev);

	if(ans <0 || ans >1){
		std::cout<<"#ERROR: Lets preserve probability shall we?\n#ERROR neutrinoModel::oscProb_dis @ prob.c\n#ERROR Prob: "<<ans<<" L: "<<L<<" Ev: "<<Ev<<std::endl;
	}

	return ans;
}


double neutrinoModel::oscProb_dis(int a, double Ev, double L){
	double ans =0.0;

	double Ua4 =0,Ua5=0,Ua6=0;

	switch(a)
	{
		case 1:
			Ua4 = Ue[0];
			Ua5 = Ue[1];
			Ua6 = Ue[2];
			break;
		case 2:
			Ua4 = Um[0];
			Ua5 = Um[1];
			Ua6 = Um[2];
			break;
		case 3:
			std::cout<<"#ERROR neutrinoModel::oscProb. taus not yet implemented!"<<std::endl;
			exit(EXIT_FAILURE);	
			break;
	}

	ans = 1.0 - 4.0*pow(Ua4*Ua5*sin(1.27*dm54Sq*L/Ev),2.0);
	ans += -4.0*pow(Ua4*Ua6*sin(1.27*dm64Sq*L/Ev),2.0);
        ans += -4.0*pow(Ua5*Ua5*sin(1.27*dm65Sq*L/Ev),2.0);
	ans += -4.0*(1.0-Ua4*Ua4-Ua5*Ua5-Ua6*Ua6)*( Ua4*Ua4*pow(sin(1.27*dm41Sq*L/Ev),2.0)+ Ua5*Ua5*pow(1.27*dm51Sq*L/Ev,2.0)+Ua6*Ua6*pow(sin(1.27*dm61Sq*L/Ev),2.0) );

	if(ans <0 || ans >1){
		std::cout<<"#ERROR: Lets preserve probability shall we?\n#ERROR neutrinoModel::oscProb_dis @ prob.c\n#ERROR Prob: "<<ans<<" L: "<<L<<" Ev: "<<Ev<<std::endl;
	}


	return ans;
}


/*********************************************
*  Other non- neutrinoModel functions, might be useful for arbitrary complex matrix U, non-unitarity?
* *******************************************/





double oscProb(int a, int b, double Ev, double L,  std::vector< std::vector < std::complex<double> > >    U, std::vector<std::vector<double> > dm ){

	double del = 0;
	if(a == b){ del = 1.0;}
	
	double ans = 0;

	//Allocate a 6x6 vector	
//	std::vector< std::vector < std::complex<double> > >  U(6,std::vector<std::complex<double>>(6) );
//	std::vector<std::vector < double > >  dm;


//	U[1][2] = std::complex<double>(10.0,1.0); 
//	std::cout<<real(U[1][2])<<" "<<arg(U[1][2])<<std::endl;


	for (int i=0; i<6; i++)
	{
		for(int j = 0; j < i; j++)
		{
		std::complex<double > temp = U[b][i]*conj(U[a][i])*conj(U[b][j])*U[a][j];
		ans += 4.0*real(temp)*pow(sin(1.27*dm[i][j]*L/Ev),2.0);
		ans += -2.0*imag(temp)*sin(2.53*dm[i][j]*L/Ev);		

		}
	}
	
	

	return del-ans;
}



double Pmue(double L, double E, double dm, double sin2)
{
	return sin2*pow( sin(1.27*dm*(L/1000.0)/E),2.0);
}


double Pmm(double L, double Ev, double Dm, double sinSq2thmm){

	double ans = 1.0-sinSq2thmm*pow( sin(Dm*(L/1000.0)*1.27/Ev), 2.0);

	if(ans <0 || ans >1){
		std::cout<<"frack. violatinf prob? "<<ans<<" l "<<L<<" Ev "<<Ev<<" Dm "<<Dm<< " sinSq2thmm "<<sinSq2thmm<<std::endl;
	}
	return ans;
}



