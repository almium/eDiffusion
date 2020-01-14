/*
 * =====================================================================================
 *
 *       Filename:  eDiffusion.cpp
 *
 *    Description:  contains main() for exciton diffusion 
 *
 *        Version:  1.0
 *        Created:  12/14/2010 14:06:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Oleksandr (Alex) Mikhnenko (www.mikhnenko.com), alex@mikhnennko.com
 *        Company:  University of Groningen
 *
 * =====================================================================================
 */

#include	<stdlib.h>
#include	<string>
#include	<iostream>
#include	<fstream>
#include	<sstream>
#include	<ctime>
#include	<math.h>
//#include	"randomc\randomc.h"

#include	"config.h"
#include	"log.h"
#include	"eDiffusion.h"

#include 	"PreProcessorOptions.h"

using namespace std;


ClassRandom Rnd;        /* gloabal instance of random generator */
Config *config; // instance of config Class
string global_folder; 	// where all the output is saved for the whole execution
			// ex: output/YYYYMMDD-HH.MM-experiment_name

double ClassExciton::hopsize; /* static variable has to be declared globally to allocate memory */
double ClassExciton::a1; /* static variable has to be declared globally to allocate memory */
double ClassExciton::tau1; /* static variable has to be declared globally to allocate memory */
double ClassExciton::a2; /* static variable has to be declared globally to allocate memory */
double ClassExciton::tau2; /* static variable has to be declared globally to allocate memory */
double ClassExciton::a3; /* static variable has to be declared globally to allocate memory */
double ClassExciton::tau3; /* static variable has to be declared globally to allocate memory */
double ClassExciton::radius; /* static variable has to be declared globally to allocate memory */
double ClassExciton::TRel_time; /* Time of thermal relaxation */
double ClassExciton::TRel_M; /* Hopsize of hot excitons (Gaussian median)  */
double ClassExciton::TRel_S; /* Hopsize of hot excitons (Gaussian sigma) */


double ClassQuencher::radius; /* static variable has to be declared globally to allocate memory */
double ClassQuencher::ASradius; /* static variable has to be declared globally to allocate memory */

double ClassMonteCarlo::dT; /* static variable has to be declared globally to allocate memory */

double ClassMedium::dX; /* static variable has to be declared globally to allocate memory */



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
	int
main ( int argc, char *argv[], char* envp[])
{
	ClassMonteCarlo *  mc;
	int i;
//	int j;
	int gen_num = 1;
	unsigned int radiative = 0;
	double Xbox,Ybox,Zbox;


	// file name related definitions
	char fname[256];
	char cmd[512];
	time_t rawtime;
	struct tm * timeinfo;
	char time_str[100];

	// output of several 
	double Output[50][4];



	/* load configuaration from experiment.cfg */
	config = new Config("experiment.cfg",envp);
	
	ClassMonteCarlo::dT 	= config->pDouble("dT"); // time discretization in ps
	ClassExciton::a1	= config->pDouble("a1"); /* weighting */
	ClassExciton::tau1 	= config->pDouble("tau1"); /* in ps */
	ClassExciton::a2 	= config->pDouble("a2"); /* weithging */
	ClassExciton::tau2 	= config->pDouble("tau2"); /* in ps */
	ClassExciton::a3 	= config->pDouble("a3"); /* weithging */
	ClassExciton::tau3 	= config->pDouble("tau3"); /* in ps */
	ClassExciton::radius 	= config->pDouble("Eradius"); /* in nm */
	ClassExciton::TRel_time = config->pDouble("TRel_time"); /* in ps */
	ClassExciton::TRel_M	= config->pDouble("TRel_M")*sqrt(ClassMonteCarlo::dT); /* in nm */
	ClassExciton::TRel_S	= config->pDouble("TRel_S")*sqrt(ClassMonteCarlo::dT); /* in nm */
	

	ClassQuencher::radius 	= config->pDouble("Qradius"); /* fullerene radius in nm */
	ClassQuencher::ASradius = config->pDouble("QASradius"); /* quencher action sphere radius in nm */
	ClassMedium::dX 	= config->pDouble("dX"); // time discretization in ps
	gen_num 		= config->pDouble("gen_num"); // time discretization in ps
	Xbox			= config->pDouble("X"); // box dimensions [nm]
	Ybox			= config->pDouble("Y"); // box dimensions [nm]
	Zbox			= config->pDouble("Z"); // box dimensions [nm] 

	cout.setf(ios::fixed,ios::floatfield);            // floatfield not set
	cout.precision(7);

	// create a folder for output  	
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (time_str, 100, "%Y%m%d-%H.%M.%S-", timeinfo);

	// folder in unix-like format
	global_folder = config->pString("outputFolder");
	global_folder += "/";
	global_folder += time_str;
	global_folder += config->pString("subfolder");
//_________________________________________________________
#ifdef _WIN32
	// make folder
	snprintf(cmd, 512, "mkdir %s\\%s%s",(config->pString("outputFolder")).c_str(),time_str,(config->pString("subfolder")).c_str());
	system (cmd);
	// copy config file
	snprintf(cmd, 512, "copy experiment.cfg %s\\%s%s",(config->pString("outputFolder")).c_str(),time_str,(config->pString("subfolder")).c_str());
	system (cmd);
#else // assuming unix
	// make folder
	snprintf(cmd, 512, "mkdir %s/%s%s",(config->pString("outputFolder")).c_str(),time_str,(config->pString("subfolder")).c_str());
	system (cmd);
	// copy config file
	snprintf(cmd, 512, "cp experiment.cfg %s/%s%s",(config->pString("outputFolder")).c_str(),time_str,(config->pString("subfolder")).c_str());
	system (cmd);
#endif
//_________________________________________________________

	// create a software version notice
	{
		ofstream versionFile;
		char path[512];

		snprintf(path, 512, "%s/%s%s/version.txt",(config->pString("outputFolder")).c_str(),time_str,(config->pString("subfolder")).c_str());
		versionFile.open(path);
		versionFile << VERSION;
		versionFile.close();
	}


	if (config->pBool("useConcentrations")){
		cout << "\nNOTE: CONCENTRATIONS are used instead of the volume fractions. ";
		cout << "See experiment.cfg for details.\n\n";
	}



	switch (config->pInt("experiment_type"))
	{
		case 0 : // one run without quenchers 
		{
			int num_points;	
			double Ld2; // Simulation() returns RMS displacement

			ClassExciton::hopsize = config->pDouble("hopsize") /* in nm per 1 ps */ 
						* sqrt(ClassMonteCarlo::dT); 
			num_points  = config->pInt("num_points");

			sprintf(fname, "One run no quenchers");	
			cout	<< "\n\n***** New simulation ***** \n" << fname << "\n" << endl;
			mc = new  ClassMonteCarlo;       /* create new simulation */
				
			mc->Init(fname,Xbox,Ybox,Zbox, 0); /* args: filename; box dimensions; volume fraction of quenchers */
				
			radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */
		
			cout 	<< "\n";		
			cout	<< "Number of radiatively decayd excitons: " << radiative << endl;
			cout 	<< "\n";		
			delete mc;                      /* free memory */
		}
		break;
	
#if !NO_Z_PERIODIC
		case 10 : // keep the hop size constant, vary the concentration
		{
			double vFrac_start, vFrac_end, vFrac_inc, vFrac_offset;
			int num_points;	
			double Ld2; // Simulation() returns RMS displacement

			ClassExciton::hopsize = config->pDouble("hopsize") /* in nm per 1 ps */ 
						* sqrt(ClassMonteCarlo::dT); 
			num_points  = config->pInt("num_points");
			if (config->pBool("useConcentrations")){
				vFrac_start = config->pDouble("conc_start")*( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) );
				vFrac_end   = config->pDouble("conc_end")*( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) );
				vFrac_offset   = config->pDouble("conc_offset")*( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) );
			}else{
				vFrac_start = config->pDouble("vFrac_start");
				vFrac_end   = config->pDouble("vFrac_end");
				vFrac_offset   = config->pDouble("vFrac_offset");
			}
			switch (config->pInt("spacing"))
			{
				case 1: // equally spaced on log scale
					if (num_points == 1) 
					{
						vFrac_inc = 0;
					}else{
						vFrac_inc = (log(vFrac_end) - log(vFrac_start) ) / (num_points-1);
					}
					break;

				case 2: // equally spaced on the linear scale
					if (num_points == 1) 
					{
						vFrac_inc = 0;
					}else{
						vFrac_inc =(vFrac_end - vFrac_start ) / (num_points-1);
					}
					break;	

				default:
					cout << "experiment.cfg: specing is invalid \n";
					exit(0);
			}


			for ( i = 0; i < num_points; i += 1 ) 
			{
				if (config->pInt("spacing") == 1 ) // equally spaced on log scale
				{
					vFrac = exp( i*vFrac_inc + log(vFrac_start) ) + vFrac_offset;
				}else if (config->pInt("spacing") == 2 ){ // equally spaced on the linear scale 
					vFrac = i*vFrac_inc + vFrac_start + vFrac_offset;
				}	

				sprintf(fname, "Volume fraction %2.5f", vFrac);	
				cout	<< "\n\n***** New simulation #" << i << " ***** \n" << fname << "\n" << endl;
				mc = new  ClassMonteCarlo;       /* create new simulation */
				
				mc->Init(fname,Xbox,Ybox,Zbox, vFrac); /* args: filename; box dimensions; volume fraction of quenchers */
				
				radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */
				Output[i][0]=vFrac; /* volume fraction  */
				//Output[i][1]=1 - (radiative*1.0)/(GEN*gen_num); /* relative quenchign efficiency */
				Output[i][1]=vFrac/( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) ); /* conc in cm -3 */
				
				Output[i][2]=(radiative*1.0); /* proportional to lifetime */
				Output[i][3]=(Ld2*1.0); /* <dL^2> */
				
				cout 	<< "\n";		
				cout	<< "Number of radiatively decayd excitons: " << Output[i][2] << endl;
				cout 	<< "\n";		
				delete mc;                      /* free memory */
				SaveTable(Output,i+1,"z.Radiative vs Amount of Quenchers.txt",true); 
			}
		}
		break;

		case 20 : //  vary hop size to get specific relative quenching efficiency
		{
			double max_hopsize, min_hopsize;
			double targetQ;
			double currentQ;
			double diffQ;
			double Ld2; // Simulation() returns RMS displacement

			max_hopsize = config->pDouble("max_hopsize") * sqrt(ClassMonteCarlo::dT);
			min_hopsize = config->pDouble("min_hopsize") * sqrt(ClassMonteCarlo::dT);
			if (config->pBool("useConcentrations")){
				vFrac = config->pDouble("Hopsize_qConc")*( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) );
			}else{
				vFrac= config->pDouble("Hopsize_vFrac");
			}
			targetQ= config->pDouble("Hopsize_targetQ");
			currentQ = 0;
			mc = new  ClassMonteCarlo;       /* create new simulation */
			mc->Init("Pre-Int",Xbox,Ybox,Zbox, vFrac); /* args: filename; box dimensions; volume fraction of quenchers */
			do 
			{
				ClassExciton::hopsize = (max_hopsize + min_hopsize)/2;
				sprintf(fname, "hopsize %2.5f", ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT) );	
			
				cout	<< "\n\n***** New simulation ***** \n" << fname << "\n" << endl;
				printf("Try hopsize %2.5f ..... \n", ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT) );	

				
				mc->Init2(fname); /* change the simulation ID and clear the array PL[]  */
				
				radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */

				currentQ=1 - (radiative*1.0)/(GEN*gen_num); /* relative quenchign efficiency */
				diffQ = targetQ-currentQ;
				if (diffQ>0)
				{
					// increase hopsize
					min_hopsize = ClassExciton::hopsize;
				}else {
					// decrease hopsize
					max_hopsize = ClassExciton::hopsize;
				}

				cout 	<< "\n";		
				cout	<< "Number of radiatively decayd excitons: " << radiative << endl;

				printf("targetQ-currentQ = %2.4f \n", diffQ);

			}while (fabs(diffQ) > 0.005 );
			delete mc;                      /* free memory */
	
			cout << "\n\n HopSize = " << ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT) << "\n";
		}			
		break;			

		case 21 : //  vary TRel_M to get specific relative quenching efficiency
		{
			double max_TRelM,min_TRelM;
			double prev_TRelM = -1;
			double targetQ;
			double currentQ;
			double prev_diffQ = -1;
			double diffQ;
			double Ld2; // Simulation() returns RMS displacement
			double ans; 

			if (config->pBool("thermal_relaxation") == false ){
				cout << "\n Please set thermal_relaxation to \"yes\" for this experiment type. \n";
				cout << "exiting\n";
				exit(0);
			}

			ClassExciton::hopsize = config->pDouble("TRelM_hopsize") * sqrt(ClassMonteCarlo::dT);
			max_TRelM = config->pDouble("max_TRelM") * sqrt(ClassMonteCarlo::dT);
			min_TRelM = config->pDouble("TRelM_hopsize") * sqrt(ClassMonteCarlo::dT);

			if (config->pBool("useConcentrations")){
				vFrac = config->pDouble("TRelM_qConc")*( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) );
			}else{
				vFrac= config->pDouble("TRelM_vFrac");
			}

			targetQ= config->pDouble("TRelM_targetQ");
			currentQ = 0;
			mc = new  ClassMonteCarlo;       /* create new simulation */
			mc->Init("Pre-Int",Xbox,Ybox,Zbox, vFrac); /* args: filename; box dimensions; volume fraction of quenchers */
		
			ClassExciton::TRel_M = (max_TRelM + min_TRelM)/2;

			do 
			{
				sprintf(fname, "TRel_M %2.5f", ClassExciton::TRel_M/sqrt(ClassMonteCarlo::dT) );	
				ans = ClassExciton::TRel_M;

				cout	<< "\n\n***** New simulation ***** \n" << fname << "\n" << endl;
				printf("Try TRel_M %2.5f ..... \n", ClassExciton::TRel_M/sqrt(ClassMonteCarlo::dT) );	

				// caclulate currentQ	
				mc->Init2(fname); /* change the simulation ID and clear the array PL[]  */
				
				radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */

				currentQ=1 - (radiative*1.0)/(GEN*gen_num); /* relative quenchign efficiency */
				diffQ = targetQ-currentQ;
			
				cout 	<< "\n";		
				cout	<< "Number of radiatively decayd excitons: " << radiative << endl;

				printf("targetQ-currentQ = %2.4f \n", diffQ);

				// calculate new TRel_M
				if (prev_diffQ == -1)
				{
					// less than two simulations
					if ( diffQ > 0 )
					{
						// increase hopsize
						min_TRelM = ClassExciton::TRel_M;
					}else{
						// decrease hopsize
						max_TRelM = ClassExciton::TRel_M;
					}
					prev_diffQ = diffQ;
					prev_TRelM = ClassExciton::TRel_M;	
					ClassExciton::TRel_M = (max_TRelM + min_TRelM)/2;
		
				}else{
					// make linear extrapolation of the new TRel_M
					double a1, a2, h1, h2, c;
					h1 = ClassExciton::TRel_M;
					h2 = prev_TRelM;
					a1 = diffQ;
					a2 = prev_diffQ;

					// update previous values
					prev_TRelM = ClassExciton::TRel_M;
					prev_diffQ = diffQ;
					
					c = (a1-a2)/(h1-h2);
					
					/*  cout << "a1 = " << a1 << "\n";
					cout << "a2 = " << a2 << "\n";
					cout << "h1 = " << h1 << "\n";
					cout << "h2 = " << h2 << "\n";
					cout << "c = " << c << "\n";
					*/
					ClassExciton::TRel_M = h2 - (a2/c);

					if (ClassExciton::TRel_M <= 0) {

						if ( diffQ > 0 )
						{
							// increase hopsize
							min_TRelM = h1; 
						}else{
							// decrease hopsize
							max_TRelM = h1; 
						}
						ClassExciton::TRel_M = (max_TRelM + min_TRelM)/2;
					}
					
				}

			}while (fabs(diffQ) > 0.005 );
			delete mc;                      /* free memory */
	
			cout << "\n\n TRel_M = " << ans/sqrt(ClassMonteCarlo::dT) << "\n";
		}			
		break;			

		case 30 : //  vary the quencher size to get specific relative quenching efficiency; hopsize and quencher conc are set
		{
			double Rmax, Rmin;
			double targetQ;
			double currentQ;
			double diffQ;
			double qConc;
			double Ld2; // Simulation() returns RMS displacement

			ClassExciton::hopsize = config->pDouble("qSize_hopsize") * sqrt(ClassMonteCarlo::dT);
			Rmax = config->pDouble("qSize_Rmax");
			Rmin = config->pDouble("qSize_Rmin");
			qConc= config->pDouble("qSize_qConc");
			targetQ= config->pDouble("qSize_targetQ");
			currentQ = 0;
			do 
			{
				ClassQuencher::radius = (Rmin + Rmax) /2;
				vFrac= qConc * PI * pow(ClassQuencher::radius*1E-7,3)*4.0/3.0 ;
				sprintf(fname, "qRadius %2.5f", ClassQuencher::radius);	
			
				cout	<< "\n\n***** New simulation ***** \n" << fname << "\n" << endl;
				printf("Try quencher radius %2.5f ..... \n", ClassQuencher::radius );	
				printf("New volume fraction is %2.5f \n", vFrac );	

				{ // here an instance of ClassMonteCarlo is created
					mc = new  ClassMonteCarlo;       /* create new simulation */
					mc->Init(fname,Xbox,Ybox,Zbox, vFrac); /* args: filename; box dimensions; volume fraction of quenchers */
					radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */
	
					currentQ=1 - (radiative*1.0)/(GEN*gen_num); /* relative quenchign efficiency */
					diffQ = targetQ-currentQ;
					if (diffQ>0)
					{
						// increase hopsize
						Rmin = ClassQuencher::radius;
					}else {
						// decrease hopsize
						Rmax = ClassQuencher::radius;
					}

					printf("targetQ-currentQ = %2.4f \n", diffQ);
	
					delete mc;                      /* free memory */
				}
			}while (fabs(diffQ) > 0.005 );

			cout << "\n\n qRadius = " << ClassQuencher::radius << "\n";
		}			
		break;			


		case 40 : // keep the hop size constant, vary the monoexponential decay time 
		{
			double tau_start, tau_end, tau_inc;
			double Ld2; // Simulation() returns RMS displacement
			int num_points;	

			ClassExciton::tau2=0;
			ClassExciton::a2=0;
			ClassExciton::a1=1;

			ClassExciton::hopsize = config->pDouble("tau_hopsize") /* in nm per 1 ps */ 
						* sqrt(ClassMonteCarlo::dT); 
			num_points  = config->pInt("tau_num_points");
			tau_start = config->pDouble("tau_start");
			tau_end   = config->pDouble("tau_end");
			
			

			if (config->pBool("useConcentrations")){
				vFrac   = config->pDouble("tau_qConc")*( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) );
			}else{
				vFrac   = config->pDouble("tau_vFrac");
			}
			switch (config->pInt("tau_spacing"))
			{
				case 1: // equally spaced on log scale
					tau_inc = (log(tau_end) - log(tau_start) ) / (num_points-1);
					break;

				case 2: // equally spaced on the linear scale
					tau_inc =(tau_end - tau_start ) / (num_points-1);
					break;	

				default:
					cout << "experiment.cfg: spacing is invalid \n";
					exit(0);
			}

			cout << "tau_start = " << tau_start << "\ttau_end = " << tau_end << "\t tau_inc = " << tau_inc << "\n";

			mc = new  ClassMonteCarlo;       /* create new simulation */
			mc->Init("Pre-Int",Xbox,Ybox,Zbox, vFrac); /* args: filename; box dimensions; volume fraction of quenchers */

			for ( i = 0; i < num_points; i += 1 ) 
			{
				if (config->pInt("tau_spacing") == 1 ) // equally spaced on log scale
				{
					ClassExciton::tau1 = exp( i*tau_inc + log(tau_start) );
				}else if (config->pInt("tau_spacing") == 2 ){ // equally spaced on the linear scale 
					ClassExciton::tau1 = i*tau_inc + tau_start;
				}	

				sprintf(fname, "tau=%2.5f", ClassExciton::tau1);	
				cout	<< "\n\n***** New simulation #" << i << " ***** \n" << fname << "\n" << endl;
				
				mc->Init2(fname); /* args: filename;  */
				
				radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */
				Output[i][0]=ClassExciton::tau1; /* tau  */
				
				Output[i][1]=(radiative*1.0); /* proportional to lifetime */
				Output[i][2]=Ld2; /* root mean square displacement */
				Output[i][3]=Ld2/sqrt(6); /* <dL^2> */

				cout 	<< "\n";		
				cout	<< "Number of radiatively decayd excitons: " << Output[i][1] << endl;
				cout 	<< "\n";		
				SaveTable(Output,i+1,"z.tau-radiative-dL2-Ld(1D).txt",false); 
			}
			delete mc;                      /* free memory */
		}
		break;


		case 50 : // vary the hopsize, keep tau and qConc constant 
		{
			double Hop_start, Hop_end, Hop_inc;
			double Ld2; // Simulation() returns RMS displacement
			int num_points;	


			ClassExciton::hopsize = config->pDouble("tau_hopsize") /* in nm per 1 ps */ 
						* sqrt(ClassMonteCarlo::dT); 
			num_points  = config->pInt("Hop_num_points");
			Hop_start = config->pDouble("Hop_start");
			Hop_end   = config->pDouble("Hop_end");
			
			

			if (config->pBool("useConcentrations")){
				vFrac   = config->pDouble("Hop_qConc")*( 4.0/3.0*PI*pow( (ClassQuencher::radius*1E-7), 3) );
			}else{
				vFrac   = config->pDouble("Hop_vFrac");
			}
			switch (config->pInt("Hop_spacing"))
			{
				case 1: // equally spaced on log scale
					Hop_inc = (log(Hop_end) - log(Hop_start) ) / (num_points-1);
					break;

				case 2: // equally spaced on the linear scale
					Hop_inc =(Hop_end - Hop_start ) / (num_points-1);
					break;	

				default:
					cout << "experiment.cfg: spacing is invalid \n";
					exit(0);
			}

			cout << "Hop_start = " << Hop_start << "\tHop_end = " << Hop_end << "\t Hop_inc = " << Hop_inc << "\n";

			mc = new  ClassMonteCarlo;       /* create new simulation */
			mc->Init("Pre-Int",Xbox,Ybox,Zbox, vFrac); /* args: filename; box dimensions; volume fraction of quenchers */

			for ( i = 0; i < num_points; i += 1 ) 
			{
				double Hop=0.0;
				double deltaT0, deltaT;

				if (config->pInt("Hop_spacing") == 1 ) // equally spaced on log scale
				{
					Hop = exp( i*Hop_inc + log(Hop_start) );
				}else { // equally spaced on the linear scale
					Hop = i*Hop_inc + Hop_start;
				}	
				deltaT0 = ( ClassExciton::a1 * pow (ClassExciton::tau1,2) + ClassExciton::a2 * pow (ClassExciton::tau2,2) ) 
					  /( ClassExciton::a1 * ClassExciton::tau1 + ClassExciton::a2 * ClassExciton::tau2 )
					  /100; // we would use at least tau/100 time sampling

				deltaT =  pow( ClassQuencher::radius / Hop * config->pDouble("Hop_dT_coeff") , 2);

				if (deltaT0<deltaT)  deltaT = deltaT0; 

				ClassMonteCarlo::dT=deltaT;

				ClassExciton::hopsize = Hop * sqrt(ClassMonteCarlo::dT);

				sprintf(fname, "hop=%2.5f",ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT));	
				cout	<< "\n\n***** New simulation #" << i << " ***** \n" << fname << "\n" << endl;
				cout << "hopsize = " << ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT) << "\t dT = " << ClassMonteCarlo::dT << "\n";

				mc->Init2(fname); /* args: filename;  */
				
				radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */
				Output[i][0]=ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT); /* hopsize  */
				Output[i][1]=pow (ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT),2)*1E-14/1E-12/6; /* diffusion coefficient  */
				
				Output[i][2]=(radiative*1.0); /* proportional to lifetime */
				Output[i][3]=Ld2; /* root mean square displacement */
				
				cout 	<< "\n";		
				cout	<< "Number of radiatively decayd excitons: " << Output[i][2] << endl;
				cout 	<< "\n";		
				SaveTable(Output,i+1,"z.hopsize-D-radiative-dL2.txt",false); 
			}
			delete mc;                      /* free memory */
		}
		break;

#else


		case 110 : // Interfacial quencher distribution. Keep the hopsize and qRadius constant, vary the voltage (quencher concentration)
		{
			double V_start, V_end;
			double Ld2; // Simulation() returns RMS displacement
			int num_points;
			double iConc_start_cm2, iConc_end_cm2, iConc_inc;
			ClassQuencher::radius = config->pDouble("e110_qRadius");


			ClassExciton::hopsize = config->pDouble("e110_hopsize") /* in nm per 1 ps */
						* sqrt(ClassMonteCarlo::dT);
			num_points  = config->pInt("e110_num_points");
			V_start = config->pDouble("e110_voltage_start");
			V_end   = config->pDouble("e110_voltage_end");

			double capacitance_cm2;
			capacitance_cm2 =  8.854E-14 * config->pDouble("IQ_dielectric_e")
											/ (config->pDouble("IQ_dielectric_thickness") * 1E-7);

			if (config->pBool("useConcentrations")){

				iConc_start_cm2 = config->pDouble("e110_conc_start");
				iConc_end_cm2 = config->pDouble("e110_conc_end");

			}else{
				// capacitance of simulation box in F/cm2

				iConc_start_cm2 = capacitance_cm2 * V_start / 1.6E-19;
				iConc_end_cm2 = capacitance_cm2 * V_end / 1.6E-19;
			}
			switch (config->pInt("Hop_spacing"))
			{
				case 1: // equally spaced on log scale
					iConc_inc = (log(iConc_end_cm2) - log(iConc_start_cm2) ) / (num_points-1);
					break;

				case 2: // equally spaced on the linear scale
					iConc_inc =(iConc_end_cm2 - iConc_start_cm2 ) / (num_points-1);
					break;

				default:
					cout << "experiment.cfg: spacing is invalid \n";
					exit(0);
			}

			cout << scientific << "iConc_start_cm2 = " << iConc_start_cm2 << "\t iConc_end_cm2 = " << iConc_end_cm2 << "\t iConc_inc = " << iConc_inc << "\n";
			cout.setf(ios::fixed,ios::floatfield);            // floatfield not set
			cout.precision(7);

			for ( i = 0; i < num_points; i += 1 )
			{
				double iConc=0.0;
				double deltaT0, deltaT;

				if (config->pInt("e110_spacing") == 1 ) // equally spaced on log scale
				{
					iConc = exp( i*iConc_inc + log(iConc_start_cm2) );
				}else { // equally spaced on the linear scale
					iConc = i*iConc_inc + iConc_start_cm2;
				}
				deltaT0 = ( ClassExciton::a1 * pow (ClassExciton::tau1,2) + ClassExciton::a2 * pow (ClassExciton::tau2,2) )
					  /( ClassExciton::a1 * ClassExciton::tau1 + ClassExciton::a2 * ClassExciton::tau2 )
					  /100; // we would use at least tau/100 time sampling

				deltaT =  pow( ClassQuencher::radius / config->pDouble("e110_hopsize") * config->pDouble("Hop_dT_coeff") , 2);

				if (deltaT0<deltaT)  deltaT = deltaT0;

				ClassMonteCarlo::dT=deltaT;

				ClassExciton::hopsize = config->pDouble("e110_hopsize") * sqrt(ClassMonteCarlo::dT);

				sprintf(fname, "iConc=%2.5e",iConc);
				cout	<< "\n\n***** New simulation #" << i << " ***** \n" << fname << "\n" << endl;
				cout << "hopsize = " << ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT) << "\t dT = " << ClassMonteCarlo::dT << "\n";


				mc = new  ClassMonteCarlo;       /* create new simulation */
				mc->Init_IQ(fname,Xbox,Ybox, config->pDouble("IQ_film_thickness"),iConc); /* args: filename; box dimensions; density of quenchers */

				cout << "Starting Monte Carlo....\n";
				radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */

				cout 	<< "\n";
				cout	<< "Number of radiatively decayed excitons: " << radiative << endl;
				cout 	<< "\n";


				//voltage = Q/C
				Output[i][0]= iConc * 1.6e-19 / capacitance_cm2;


				Output[i][1]=iConc; // Concentration
				Output[i][2]=(radiative*1.0); // proportional to lifetime
				Output[i][3]=Ld2; // root mean square displacement


				SaveTable(Output,i+1,"z.V-Conc_cm2-radiative-dL2.txt",true);

				cout << "Deleting the Monte Carlo object \n";
				delete mc;                      /* free memory */
				cout << "The Monte Carlo object has been deleted \n";
			}
			cout << "\nFinished.\n";
		}
		break;


		case 120 : // Interfacial quencher distribution.  Keep the hopsize and voltage constant, vary the quencher size

		{

			double Ld2; // Simulation() returns RMS displacement
			double qRadius_min, qRadius_max, qRadius;
			double iConc;
			double diffQ = 0.0;
			double targetQ = 0.0;
			double currentQ = 0.0;
			double deltaT0, deltaT;

			ClassExciton::hopsize = config->pDouble("e120_hopsize") /* in nm per 1 ps */
						* sqrt(ClassMonteCarlo::dT);
			qRadius_min = config->pDouble("e120_qRadius_min");
			qRadius_max = config->pDouble("e120_qRadius_max");
			targetQ = config->pDouble("e120_targetQ");

			double capacitance_cm2;
			capacitance_cm2 =  8.854E-14 * config->pDouble("IQ_dielectric_e")
											/ (config->pDouble("IQ_dielectric_thickness") * 1E-7);

			if (config->pBool("useConcentrations")){

				iConc = config->pDouble("e120_iConc");

			}else{
				// capacitance of simulation box in F/cm2
				iConc = capacitance_cm2 * config->pDouble("e120_voltage") / 1.6E-19;

			}

			cout << scientific << "iConc = " << iConc << "\n";
			cout.setf(ios::fixed,ios::floatfield);            // floatfield not set
			cout.precision(7);

			deltaT0 = ( ClassExciton::a1 * pow (ClassExciton::tau1,2) + ClassExciton::a2 * pow (ClassExciton::tau2,2) )
				  /( ClassExciton::a1 * ClassExciton::tau1 + ClassExciton::a2 * ClassExciton::tau2 )
				  /100; // we would use at least tau/100 time sampling

			deltaT =  pow( ClassQuencher::radius / config->pDouble("e120_hopsize") * config->pDouble("Hop_dT_coeff") , 2);

			if (deltaT0<deltaT)  deltaT = deltaT0;

			ClassMonteCarlo::dT=deltaT;
			ClassExciton::hopsize = config->pDouble("e120_hopsize") * sqrt(ClassMonteCarlo::dT);

			qRadius = (qRadius_min + qRadius_max) / 2.0;

			int i = 1;
			do {
				ClassQuencher::radius = qRadius = (qRadius_min + qRadius_max) / 2.0;

				sprintf(fname, "qRadius=%2.5f",qRadius);
				cout	<< "\n\n***** New simulation #" << i << " ***** \n" << fname << "\n" << endl;
				cout << "hopsize = " << ClassExciton::hopsize/sqrt(ClassMonteCarlo::dT) << "\t dT = " << ClassMonteCarlo::dT << "\n";
				cout << "qRadius = " << qRadius << "\n";

				mc = new  ClassMonteCarlo;       /* create new simulation */
				mc->Init_IQ(fname,Xbox,Ybox, config->pDouble("IQ_film_thickness"),iConc); /* args: filename; box dimensions; density of quenchers */

				cout << "Starting Monte Carlo....\n";
				radiative = mc->Simulation(gen_num,&Ld2);   /* args: number of generations */

				cout 	<< "\n";
				cout	<< "Number of radiatively decayed excitons: " << radiative << endl;
				cout 	<< "\n";

				currentQ = 1.0 - 1.0*radiative/(config->pDouble("gen_num") * 10000.0);
				diffQ = targetQ-currentQ;
				if (diffQ>0.0) {
					// increase radius
					qRadius_min = qRadius;
				}else{
					// decrease radius
					qRadius_max = qRadius;
				}

				cout << "targetQ-currentQ = " << diffQ << "\n";

				cout << "Deleting the Monte Carlo object \n";
				delete mc;                      /* free memory */
				cout << "The Monte Carlo object has been deleted \n";
				i++;
			}while (fabs(diffQ) > 0.005 );

			cout << "\n\n  quencher Radius is " << qRadius << "\n Finished. \n";


		}
		break;
#endif

		default:
			cout << "Wrong experiment type. \nIf you use experiment type with number, which is less than 100, then check NO_Z_PERIODIC in PreProcessorOptions.h \n";

	}

		
//	for ( j = 4; j < 5; j += 1 ) {
//		//ClassExciton::hopsize= 0.11 * j * sqrt(dT) ; 
//		ClassExciton::hopsize= 0.253 * sqrt(dT) ; 
////		ClassExciton::hopsize= 0.5 * sqrt(dT) ; 
//	
//		for ( i = 1; i < 2; i += 1 ) 
//		{
//			vFrac = i*0.004;
//			vFrac = 0.02;
//			sprintf(fname, "Hopsize %2.3f, Volume ratio %2.3f", ClassExciton::hopsize,vFrac);	
//			cout	<< "\n\n***** New simulation ***** \n" << fname << "\n" << endl;
//			mc = new  ClassMonteCarlo;       /* create new simulation */
//			
//			mc->Init(fname,50,50,50, vFrac); /* args: filename; box dimensions; volume fraction of quenchers */
//			
//			radiative = mc->Simulation(gen_num);   /* args: number of generations. Returns number of excitons, which decayed radiatively */
//			
//			
//			Output[i][0]=vFrac; /* volume ratio */
//			//Output[i][1]=1 - (radiative*1.0)/(GEN*gen_num); /* relative quenchign efficiency */
//			Output[i][1]=(radiative*1.0); /* proportional to lifetime */
//			
//			cout	<< Output[i][1] << endl;
//			delete mc;                      /* free memory */
//		}
//	
//		sprintf(fname, "Radiative vs VolumeRatio, hopsize  %2.3f.txt", ClassExciton::hopsize);	
//		SaveTable(Output,10,fname); // !! this functionality should be outside ClassMonteCarlo`
//
//	}

	delete config;

	return EXIT_SUCCESS;
}




void
SaveTable ( double  arr[][4], int n, string filename, bool UseDate )
{
	int i;
	ofstream datafile;
	string fname="";
	time_t rawtime;
	struct tm * timeinfo;
	char time_str[100];

	fname += global_folder; 
       	fname += "/";

	if (UseDate) {
		time (&rawtime);
		timeinfo = localtime (&rawtime);
		strftime (time_str, 100, "%H.%M.%S-", timeinfo);

		fname += time_str;	
	}
	fname +=filename;


	
	datafile.open (fname.c_str(),ios::out | ios::trunc);
	datafile.precision(7);
	for ( i = 0; i < n; i += 1 ) {
		datafile << arr[i][0] << "\t" <<  arr[i][1] << "\t" <<  arr[i][2] << "\t" << arr[i][3] << endl;	
	}
	datafile.close();


	return ;
}		/* -----  end of method ClassMonteCarlo::SavePL  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassRandom
 *      Method:  ClassRandom :: ClassRandom
 * Description:  Constructor of ClassRandom
 *--------------------------------------------------------------------------------------
 */

ClassRandom::ClassRandom ()
{
	seed = (int)time(0);            // random seed
	RandGen = new CRandomMersenne(seed);
	RandGen->RandomInit(seed);	
	sto = new StochasticLib1(seed);
	
	return ;
}		/* -----  end of method ClassRandom::ClassRandom  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassRandom
 *      Method:  ClassRandom :: ~ClassRandom
 * Description:  destructor of ClassRandom
 *--------------------------------------------------------------------------------------
 */

ClassRandom::~ClassRandom ( )
{
	delete RandGen;
	delete sto;
	return ;
}		/* -----  end of method ClassRandom::~ClassRandom  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassRandom
 *      Method:  ClassRandom :: dRnd
 * Description:  
 *--------------------------------------------------------------------------------------
 */

	double	
ClassRandom::dRnd ( )
{
	return 1 - RandGen->Random(); /* random value in the range 0<r<=1 */
}		/* -----  end of method ClassRandom::iRand  ----- */


	int
ClassRandom::iRnd (int max )
{
	return RandGen->IRandomX(0,max-1); /* random integer in the range 0 <= r <= max -1 */
}		/* -----  end of method ClassRandom::iRand  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassRandom
 *      Method:  ClassRandom :: dNormal
 * Description:  
 *--------------------------------------------------------------------------------------
 */
	double
ClassRandom::dNormal ( double center, double sigma )
{
	double r;
	do {
		r=sto->Normal(center, sigma);
	} while (r<=0);
	return  r;
}		/* -----  end of method ClassRandom::dNormal  ----- */


