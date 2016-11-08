
#include "wrkInstance.h"
wrkInstance::~wrkInstance(){

		delete bkgspec;
		delete bkgbarspec;

		delete UBOONE;
		delete SBND;
		delete ICARUS;

}


wrkInstance::wrkInstance(int channel_mode, int fbeam_mode) : wrkInstance(channel_mode,fbeam_mode,1.0,1.0) {}

wrkInstance::wrkInstance(int channel_mode, int fbeam_mode, double pot_scale, double pot_scale_bar){
	isVerbose = true;
	which_mode = channel_mode;
	beam_mode = fbeam_mode;
	
	double chi2 = 0; //old chi to be passed in

	nullModel = neutrinoModel();
	nullModel.zero();
	workingModel= neutrinoModel();
	workingModel.zero();

	UBOONE = new SBN_detector(1);
	SBND =  new SBN_detector(0);
	ICARUS =  new SBN_detector(2);

	bkgspec = new SBN_spectrum(nullModel);
	bkgbarspec = new SBN_spectrum(nullModel);
	
	bkgbarspec->SetNuBarMode();


		bkgspec->load_bkg(ICARUS);
		bkgspec->load_bkg(SBND);
		bkgspec->load_bkg(UBOONE);
	
		bkgspec->scale_by_pot(pot_scale);

			back6 = bkgspec->get_sixvector();
			back  = bkgspec->get_vector();

		bool usedetsys = true;	
		bool stat_only = false;
	
	if (beam_mode == 0){
	

		matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
		matrix_size_c = (N_e_bins + N_m_bins) * N_dets;
		
		bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		contMsize = (N_e_bins+N_m_bins)*N_dets;

		/* create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */

		 TMatrixT <double> m(matrix_size,matrix_size);
		 TMatrixT <double>  mc(matrix_size_c,matrix_size_c);
		 TMatrixT <double>  mci(matrix_size_c, matrix_size_c);
		 TMatrixT <double>  msys(bigMsize,bigMsize);

		
		 msys = sys_fill_direct(msys.GetNcols(),usedetsys);

		for(int i =0; i<msys.GetNcols(); i++)
		{
			for(int j =0; j<msys.GetNrows(); j++)
			{
		//		std::cout<<i<<" "<<j<<" "<<msys(i,j)<<std::endl;
				msys(i,j)=msys(i,j)*back[i]*back[j];
			}
		}



		TMatrixT <double> mstat(bigMsize,bigMsize);
		stats_fill(mstat, back);

		TMatrixT <double > mtotal(bigMsize,bigMsize);

		if(stat_only){
			mtotal =  mstat;
		} else {
			mtotal = msys+mstat;
		}

		TMatrixT <double> mctotal(contMsize,contMsize);
		contract_signal2(mtotal,mctotal);


		double invdet=0; // just to hold determinant

		//	bit o inverting, root tmatrix seems perfectly fast	
		mci = mctotal.Invert(&invdet);

		vMcI = to_vector(mci);

	/*	*/
	
	//	std::cout<<i<<" input_mn: "<<m4<<" "<<m5<<" "<<m6<<" input_ue "<<ue4<<" "<<ue5<<" "<<ue6<<" input_um4: "<<um4<<" "<<um5<<" "<<um6<<" input_chi: "<<chi2<<" "<<std::endl;

	} else if(beam_mode == 1)
	{
		
		bkgbarspec->load_bkg(ICARUS);
		bkgbarspec->load_bkg(SBND);
		bkgbarspec->load_bkg(UBOONE);
		bkgbarspec->SetNuBarMode();

		//if(pot_scale_bar !=1.0){
		//just scale it as muboone isnt right
		bkgbarspec->scale_by_pot(pot_scale_bar);
		//}

		backbar6 = bkgbarspec->get_sixvector();
		backbar  = bkgbarspec->get_vector();

		back_all = back;
		back_all.insert(back_all.end(), backbar.begin(), backbar.end() );

		back_all_12 = back6;
		back_all_12.insert(back_all_12.end(), backbar6.begin(), backbar6.end() );
	

		bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti;
		contMsize = (N_e_bins+N_m_bins)*N_dets*N_anti;

		/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */

		TMatrixT <double> McI(contMsize, contMsize);
		//TMatrixT <double> McIbar(contMsize, contMsize);

		// Fill systematics from pre-computed files
		TMatrixT <double> Msys(bigMsize,bigMsize);
		Msys = sys_fill_direct(Msys.GetNcols(),usedetsys);

		/*for(int i =0; i< Msys.GetNcols(); i++){
		for(int j =0; j< Msys.GetNcols(); j++){
		std::cout<<i<<" "<<j<<" "<<Msys(i,j)<<" "<<Msys(i,j)<<std::endl;
		}}
		exit(EXIT_FAILURE);
		*/

		// systematics per scaled event
		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*back_all[i]*back_all[j];
			}
		}
			
		// Fill stats from the back ground vector
		TMatrixT <double> Mstat(bigMsize,bigMsize);
		stats_fill(Mstat, back_all);

		//And then define the total covariance matrix in all its glory
		TMatrixT <double > Mtotal(bigMsize,bigMsize);
		if(stat_only){
			Mtotal =  Mstat;
		} else {
			Mtotal = Msys+Mstat;
		}
		
		
		// Now contract back the larger antimatrix
		TMatrixT<double > Mctotal(contMsize,contMsize);

		contract_signal2_anti(Mtotal,Mctotal);

		/*for(int i =0; i< Mtotal.GetNcols(); i++){
		for(int j =0; j< Mtotal.GetNcols(); j++){
		std::cout<<i<<" "<<j<<" "<<Mtotal(i,j)<<" "<<Msys(i,j)<<std::endl;
		}}
		exit(EXIT_FAILURE);
		*/
		vMc = to_vector(Mctotal);

		// just to hold determinant
		double invdet=0; 

		// Bit o inverting, root tmatrix seems perfectly and sufficiently fast for this, even with anti_mode
		McI = Mctotal.Invert(&invdet);


		// There is currently a bug, somehow a memory leak perhaps. converting the TMatrix to a vector of vectors fixes it for now. 
		vMcI = to_vector(McI);
	} //end anti_initialiser



	delete UBOONE;
	delete ICARUS;
	delete SBND;



}//end wrkInstance constructor;


int wrkInstance::clear_all(){


	//	~SigSpec(); //make a destructor
		pred6.clear();
		predbar6.clear();
		pred_all_12.clear();
		Current_Chi = -9999;

return 1;
}

int wrkInstance::init_minim(){
	isVerbose = false;

//	min = new ROOT::Math::GSLMinimizer(ROOT::Math::kConjugateFR);	
	min = new ROOT::Math::GSLMinimizer(ROOT::Math::kVectorBFGS2);	
//	min = new ROOT::Math::GSLSimAnMinimizer();
	//min->SetMaxFunctionCalls(100); // for Minuit/Minuit2
   	min->SetMaxIterations(100);  // for GSL
	min->SetTolerance(0.001); //times 4 for normal
	min->SetPrintLevel(0);
	min->SetPrecision(0.0001);//times 4 for normal
return 1;
}

double wrkInstance::minim_calc_chi(const double * x){
	double ans = 99999;

	double Imn[3] = {x[0],x[1],x[2]};
	double Iue[3] = {pow(10,x[3]),pow(10,x[4]),x[5]};
	double Ium[3] = {pow(10,x[6]),pow(10,x[7]),x[8]};
	double Iphi[3] = {x[9],x[10],x[11]};
	neutrinoModel tmpModel(Imn,Iue,Ium,Iphi);

	ans = this->calc_chi(tmpModel,1, pot, pot_bar);

	/*if(ans<0){std::cout<<"WARNING: chi^2 is less than 0 in minimizer: "<<ans<<std::endl;
			std::cout<<"mass: "<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
			std::cout<<"ue: "<<Iue[0]<<" "<<Iue[1]<<" "<<Iue[2]<<std::endl;
			std::cout<<"um: "<<Ium[0]<<" "<<Ium[1]<<" "<<Ium[2]<<std::endl;
			std::cout<<"phi: "<<Iphi[0]<<" "<<Iphi[1]<<" "<<Iphi[2]<<std::endl;
		exit(EXIT_FAILURE);
	}*/

	return ans;

}

double wrkInstance::minimize(double inphi45, double ipot, double ipotbar){

	pot=ipot;
	pot_bar =ipotbar;

	//
        ROOT::Math::Functor f( this, &wrkInstance::minim_calc_chi,12); 
	TRandom3 *rangen    = new TRandom3(0);


	double variable[12] = {0.398107,1.0,0,log10(rangen->Uniform(0,0.2)),log10(rangen->Uniform(0,0.2)),0,log10(rangen->Uniform(0,0.2)),log10(rangen->Uniform(0,0.2)),0, inphi45,0.0,0.0};
	double step[12] = {0.01,0.01,0.01, 0.005,0.005,0.005, 0.005,0.005,0.001,  0.01,0.01,0.01};
	double lower[12] = {0,0,0,-4,-4,0,-4,-4,0,0,0,0}	;
	double myup=0.3;
	double upper[12] = {1,1,1,log10(myup),log10(myup),0.3,log10(myup),log10(myup),0.3,2*3.14159,2*3.14159,2*3.14159};	
	
	std::string name[12] ={"Dm41\0","Dm51","Dm61","Ue4\0","Ue5","Ue6","Um4","Um5","Um6","phi45","phi46","phi56"};
	int isfixed[12]={1,1,1,0,0,1,0,0,1,1,1,1};

   min->SetFunction(f);

   for(int i=0;i<12;i++){
	if(isfixed[i]){
	   	min->SetFixedVariable(i,name[i],variable[i]);
	} else {

   		min->SetLimitedVariable(i,name[i],variable[i], step[i], lower[i],upper[i]);
	}

   }
   min->Minimize(); 
   //            
   
  const double *xs = min->X();
   std::cout<<inphi45<<" Minimum: ";
   for(int i=0; i<11; i++){
	if(!isfixed[i]){
		std::cout<<" "<<xs[i];
	}
	}
	std::cout<<" : " << wrkInstance::minim_calc_chi(xs) << std::endl;
   //                           

	min->Clear();
return wrkInstance::minim_calc_chi(xs);
}



double wrkInstance::calc_chi(neutrinoModel newModel, int runnumber){ return calc_chi(newModel, runnumber, 1.0, 1.0);}

double wrkInstance::calc_chi(neutrinoModel newModel, int runnumber, double pot_scale, double pot_scale_bar){

				this->clear_all();


				double chi2 = 0;
				int i = runnumber;
		
				workingModel=newModel;	
				SigSpec = new SBN_spectrum(workingModel);
				SigSpec->which_mode = which_mode;
				
				SigBarSpec =new SBN_spectrum(workingModel);
				SigBarSpec->which_mode=which_mode;   //AHAHAHAHAHAHA!!!!
				SigBarSpec->SetNuBarMode();


				SigSpec->load_freq_3p3(ICARUS);//0 is silly app flag (get rid of this)
				SigSpec->load_freq_3p3(SBND);
				SigSpec->load_freq_3p3(UBOONE);
				if(pot_scale!=1){
					SigSpec->scale_by_pot(pot_scale);
				}

				pred6 = SigSpec->get_sixvector();

				int whatsize = vMcI[0].size();


				double mychi2=0;
		if(beam_mode == NU_MODE){

				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						mychi2 += (back6[i]-pred6[i])*vMcI[i][j]*(back6[j]-pred6[j]);
					}
				}


		}else if(beam_mode == NU_NUBAR_MODE)
		{
			
				SigBarSpec->load_freq_3p3(ICARUS);
				SigBarSpec->load_freq_3p3(SBND);
				SigBarSpec->load_freq_3p3(UBOONE);
				SigBarSpec->scale_by_pot(pot_scale_bar); // THIS IS SUPER NECESSARY, to scale muboone out.	
				
				predbar6 = SigBarSpec->get_sixvector();

				pred_all_12 = pred6;	
			     	pred_all_12.insert( pred_all_12.end(), predbar6.begin(), predbar6.end() );
				whatsize=pred_all_12.size();	
		
				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						mychi2 += (back_all_12[i]-pred_all_12[i])*vMcI[i][j]*(back_all_12[j]-pred_all_12[j]);
					//	std::cout<<i<<" "<<j<<" "<<vMcI[i][j]<<" "<<vMc[i][j]<<" "<<vMc[i][j]/(back_all_12[i]*back_all_12[j])<<std::endl;
					}
				}
			
		
		}// end anti-mode

			Current_Chi = mychi2;

		if(mychi2<0){
			std::cout<<"ERROR: wrkInstance::calc_chi Chi^2 <0: "<<mychi2<<std::endl;
	
			std::cout<<"#"<<i<<" m4: "<<SigSpec->workingModel.mNu[0]<<" m5: "<<SigSpec->workingModel.mNu[1]<<" m6: "<<SigSpec->workingModel.mNu[2]<<" ue4 "<<SigSpec->workingModel.Ue[0]<<" ue5 "<<SigSpec->workingModel.Ue[1]<<" ue6  "<<SigSpec->workingModel.Ue[2]<<" um4 "<<SigSpec->workingModel.Um[0]<<" um5 "<<SigSpec->workingModel.Um[1]<<" um6 "<<SigSpec->workingModel.Um[2]<<" phi45 "<<SigSpec->workingModel.phi[0]<<" phi46 "<<SigSpec->workingModel.phi[1]<<" phi56 "<<SigSpec->workingModel.phi[2]<<" pot "<<SigSpec->pot_scaling<<" "<<which_mode<<" "<<chi2<<" "<<Current_Chi<<" "<<std::endl;

			std::vector< double > printvecB;

			std::vector< double > printvec;

			if(beam_mode==1){ printvec = pred_all_12; printvecB = back_all_12;} else {printvec = pred6;printvecB = back6;}
			std::cout<<"Background: "<<std::endl;
			std::cout<<printvecB[0];
			for(int u=1;u< printvecB.size(); u++){
				std::cout<<" "<<printvecB[u];
			}	
			std::cout<<std::endl;


			std::cout<<"Signal: "<<std::endl;
			std::cout<<printvec[0];
			for(int u=1;u< printvec.size(); u++){
				std::cout<<" "<<printvec[u];
			}	
			std::cout<<std::endl;
/*
			std::cout<<"MATRIX: "<<std::endl;
				double tmpchi=0;
				for(int i =0; i<printvec.size(); i++){
					for(int j =0; j<printvec.size(); j++){
						double thisone = (printvecB[i]-printvec[i])*vMcI[i][j]*(printvecB[j]-printvec[j]);
						tmpchi += thisone;
						std::cout<<i<<" "<<j<<" "<<vMcI[i][j]<<" "<<vMc[i][j]<<" "<<vMc[i][j]/(printvecB[i]*printvecB[j])<<" "<<tmpchi<<" "<<thisone<<std::endl;
					}
				}
			
		*/


			exit(EXIT_FAILURE);


		}


if(isVerbose){
std::cout<<"#"<<i<<" "<<SigSpec->workingModel.mNu[0]<<" "<<SigSpec->workingModel.mNu[1]<<" "<<SigSpec->workingModel.mNu[2]<<" "<<SigSpec->workingModel.Ue[0]<<" "<<SigSpec->workingModel.Ue[1]<<" "<<SigSpec->workingModel.Ue[2]<<" "<<SigSpec->workingModel.Um[0]<<" "<<SigSpec->workingModel.Um[1]<<" "<<SigSpec->workingModel.Um[2]<<" "<<SigSpec->workingModel.phi[0]<<" "<<SigSpec->workingModel.phi[1]<<" "<<SigSpec->workingModel.phi[2]<<" "<<SigSpec->pot_scaling<<" "<<which_mode<<" "<<chi2<<" "<<Current_Chi<<" "<<std::endl;
			std::vector< double > printvec;
			if(beam_mode==1){ printvec = pred_all_12; } else {printvec = pred6;}

			std::cout<<printvec[0];
			for(int u=1;u< printvec.size(); u++){
				std::cout<<" "<<printvec[u];
			}	
			std::cout<<std::endl;
}

	delete SigBarSpec;
	delete SigSpec;

	return Current_Chi;

}

double wrkInstance::calc_chi_POT_vector(neutrinoModel newModel, std::vector<double> vecin, int runnumber , double potin, double potinbar){



				this->clear_all();

				double chi2 = 0;
				int i = runnumber;

				int whatsize = vecin.size();

				double mychi2=0;
			
				//check for previous known bug!
			/*	if(false && matrix_size_c != vecin.size() && matrix_size_c != back6.size())
				{
					std::cout<<"#ERROR, soemthing wrong lengthwise"<<std::endl;
					std::cout<<"#ERROR, matrix_size_c: "<<matrix_size_c<<" pred: "<<vecin.size()<<" back: "<<back6.size()<<std::endl;	
				}
			*/
				//Calculate the answer, ie chi square! will functionise
				// should be matrix_size_c for full app+dis

			if(beam_mode==0){
				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
					double tmp = (back6[i]-vecin[i])*vMcI[i][j]*(back6[j]-vecin[j]);
					//std::cout<<i<<" "<<j<<" "<<back6[i]<<" "<<vecin[i]<<" "<<back6[j]<<" "<<vecin[j]<<" "<<mychi2<<" "<<tmp<<" "<<vMcI[i][j]<<std::endl;						
						mychi2+= tmp;
					}
				}

			}else if(beam_mode ==1){
				for(int i =0; i<whatsize; i++){
					for(int j =0; j<whatsize; j++){
						double tmp = (back_all_12[i]-vecin[i])*vMcI[i][j]*(back_all_12[j]-vecin[j]);
				//			std::cout<<i<<" "<<j<<" "<<back_all_12[i]<<" "<<vecin[i]<<" "<<back_all_12[j]<<" "<<vecin[j]<<" "<<mychi2<<" "<<tmp<<" "<<vMcI[i][j]<<" "<<vMc[i][j]<<std::endl;						
							mychi2 +=tmp;
					}
				}
			}
			




	std::cout<<i<<" "<<newModel.mNu[0]<<" "<<newModel.mNu[1]<<" "<<newModel.mNu[2]<<" "<<newModel.Ue[0]<<" "<<newModel.Ue[1]<<" "<<newModel.Ue[2]<<" "<<newModel.Um[0]<<" "<<newModel.Um[1]<<" "<<newModel.Um[2]<<" "<<newModel.phi[0]<<" "<<newModel.phi[1]<<" "<<newModel.phi[2]<<" "<<potin<<" "<<potinbar<<" "<<0.0<<" "<<mychi2<<std::endl;
/*
			std::cout<<vecin[0];
			for(int u=1;u< vecin.size(); u++){
				std::cout<<" "<<vecin[u];
			}	
			std::cout<<std::endl;
		
			std::cout<<back_all_12[0];
			for(int u=1;u< back_all_12.size(); u++){
				std::cout<<" "<<back_all_12[u];
			}	
			std::cout<<std::endl;
	

			exit(EXIT_FAILURE);
*/
			Current_Chi = mychi2;

	return Current_Chi;

}


double wrkInstance::inject_signal(neutrinoModel signalModel, int channel_mode, int fbeam_mode, double pot_scale, double pot_scale_bar ){

	back6.clear();
	this->clear_all();
	vMcI.clear();


	which_mode = channel_mode;
	beam_mode = fbeam_mode;
	
	double chi2 = 0; //old chi to be passed in

	nullModel = neutrinoModel();
	nullModel.zero();
	workingModel= neutrinoModel();
	workingModel.zero();
	workingModel= signalModel;

	UBOONE = new SBN_detector(1);
	SBND =  new SBN_detector(0);
	ICARUS =  new SBN_detector(2);

	//background is NOW my signal, dont be confused
	bkgspec = new SBN_spectrum(signalModel);
	bkgbarspec = new SBN_spectrum(signalModel);
	
	bkgspec->which_mode=which_mode;
	bkgbarspec->which_mode=which_mode;
	bkgbarspec->SetNuBarMode();
				

	bkgspec->load_freq_3p3(ICARUS);//0 is silly app flag (get rid of this)
	bkgspec->load_freq_3p3(SBND);
	bkgspec->load_freq_3p3(UBOONE);
	
	bkgspec->scale_by_pot(pot_scale);


	back6 = bkgspec->get_sixvector();
	back  = bkgspec->get_vector();

		bool usedetsys = true;	
		bool stat_only = false;
	
	if (beam_mode == 0){
	

		matrix_size =(N_e_bins + N_e_bins + N_m_bins)*N_dets;
		matrix_size_c = (N_e_bins + N_m_bins) * N_dets;
		
		bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets;
		contMsize = (N_e_bins+N_m_bins)*N_dets;

		/* create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */

		 TMatrixT <double> m(matrix_size,matrix_size);
		 TMatrixT <double>  mc(matrix_size_c,matrix_size_c);
		 TMatrixT <double>  mci(matrix_size_c, matrix_size_c);
		 TMatrixT <double>  msys(bigMsize,bigMsize);

		
		msys=sys_fill_direct(msys.GetNcols(),usedetsys);

		for(int i =0; i<msys.GetNcols(); i++)
		{
			for(int j =0; j<msys.GetNrows(); j++)
			{
				msys(i,j)=msys(i,j)*back[i]*back[j];
			}
		}



		TMatrixT <double> mstat(bigMsize,bigMsize);
		stats_fill(mstat, back);

		TMatrixT <double > mtotal(bigMsize,bigMsize);

		if(stat_only){
			mtotal =  mstat;
		} else {
			mtotal = msys+mstat;
		}

		TMatrixT <double> mctotal(contMsize,contMsize);
		contract_signal2(mtotal,mctotal);


		double invdet=0; // just to hold determinant

		//	bit o inverting, root tmatrix seems perfectly fast	
		mci = mctotal.Invert(&invdet);

		vMcI = to_vector(mci);

	/*	*/
	
	//	std::cout<<i<<" input_mn: "<<m4<<" "<<m5<<" "<<m6<<" input_ue "<<ue4<<" "<<ue5<<" "<<ue6<<" input_um4: "<<um4<<" "<<um5<<" "<<um6<<" input_chi: "<<chi2<<" "<<std::endl;

	} else if(beam_mode == 1)
	{
		bkgbarspec->SetNuBarMode();
		bkgbarspec->load_freq_3p3(ICARUS);//0 is silly app flag (get rid of this)
		bkgbarspec->load_freq_3p3(SBND);
		bkgbarspec->load_freq_3p3(UBOONE);
		
		bkgbarspec->scale_by_pot(pot_scale_bar);

		backbar6 = bkgbarspec->get_sixvector();
		backbar  = bkgbarspec->get_vector();

		back_all = back;
		back_all.insert(back_all.end(), backbar.begin(), backbar.end() );

		back_all_12 = back6;
		back_all_12.insert(back_all_12.end(), backbar6.begin(), backbar6.end() );
	

		bigMsize = (N_e_bins*N_e_spectra+N_m_bins*N_m_spectra)*N_dets*N_anti;
		contMsize = (N_e_bins+N_m_bins)*N_dets*N_anti;

		/* Create three matricies, full 9x9 block, contracted 6x6 block, and inverted 6x6
		 * */

		TMatrixT <double> McI(contMsize, contMsize);
		//TMatrixT <double> McIbar(contMsize, contMsize);

		// Fill systematics from pre-computed files
		TMatrixT <double> Msys(bigMsize,bigMsize);
		Msys=sys_fill_direct(Msys.GetNcols(),usedetsys);

		/*for(int i =0; i< Msys.GetNcols(); i++){
		for(int j =0; j< Msys.GetNcols(); j++){
		std::cout<<i<<" "<<j<<" "<<Msys(i,j)<<" "<<Msys(i,j)<<std::endl;
		}}
		exit(EXIT_FAILURE);
		*/

		// systematics per scaled event
		for(int i =0; i<Msys.GetNcols(); i++)
		{
			for(int j =0; j<Msys.GetNrows(); j++)
			{
				Msys(i,j)=Msys(i,j)*back_all[i]*back_all[j];
			}
		}
			
		// Fill stats from the back ground vector
		TMatrixT <double> Mstat(bigMsize,bigMsize);
		stats_fill(Mstat, back_all);

		//And then define the total covariance matrix in all its glory
		TMatrixT <double > Mtotal(bigMsize,bigMsize);
		if(stat_only){
			Mtotal =  Mstat;
		} else {
			Mtotal = Msys+Mstat;
		}
		
		
		// Now contract back the larger antimatrix
		TMatrixT<double > Mctotal(contMsize,contMsize);

		contract_signal2_anti(Mtotal,Mctotal);

		/*for(int i =0; i< Mtotal.GetNcols(); i++){
		for(int j =0; j< Mtotal.GetNcols(); j++){
		std::cout<<i<<" "<<j<<" "<<Mtotal(i,j)<<" "<<Msys(i,j)<<std::endl;
		}}
		exit(EXIT_FAILURE);
		*/
		vMc = to_vector(Mctotal);

		// just to hold determinant
		double invdet=0; 

		// Bit o inverting, root tmatrix seems perfectly and sufficiently fast for this, even with anti_mode
		McI = Mctotal.Invert(&invdet);


		// There is currently a bug, somehow a memory leak perhaps. converting the TMatrix to a vector of vectors fixes it for now. 
		vMcI = to_vector(McI);
	} //end anti_initialiser

}//end the inject_signal;




