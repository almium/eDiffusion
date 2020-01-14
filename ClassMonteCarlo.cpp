/*
 * =====================================================================================
 *
 *       Filename:  ClassMonteCarlo.cpp
 *
 *    Description:  Contains declarations of ClassMonteCarlo 
 *
 *        Version:  1.0
 *        Created:  1/20/2011 10:05:58
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
#include	"eDiffusion.h"
#include	"randomc.h"
#include	"config.h"
#include	"log.h"

#include 	"PreProcessorOptions.h"

extern ClassRandom Rnd;        /* gloabal instance of random generator */
extern Config *config;		/* global instance of config class */
extern string global_folder; 	// where all the output is saved for the whole execution
			// ex: output/YYYYMMDD-HH.MM-experiment_name

/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: ClassMonteCarlo
 * Description:  Constructor of ClassMonteCarlo
 *-------------------------------------------------------------------------------------:call C_HlMessage()
 -
 */


ClassMonteCarlo::ClassMonteCarlo ()
{
	SimulationID="";	


	return;
}		/* -----  end of method ClassMonteCarlo::ClassMonteCarlo  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: ~ClassMonteCarlo
 * Description:  Destructor of ClassMonteCarlo 
 *-------------------------------------------------------------------------------------:call C_HlMessage()
 -
 */


ClassMonteCarlo::~ClassMonteCarlo ()
{
	delete medium;
//	PL.resize (TN,0);

	return;
}		/* -----  end of method ClassMonteCarlo::ClassMonteCarlo  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: Init 
 * Description:  Clears up the instance of simulation 
 *-------------------------------------------------------------------------------------
 -
 */

void
ClassMonteCarlo::Init(string sID, double bx, double by, double bz, double frac) /* frac - volume fraction of fullerenes */
{
	int i; 
	SimulationID=sID;	
	for ( i = 0; i < TN; i += 1 ) {
		PL[i]=0;
	}
	medium = new ClassMedium(bx,by,bz,frac);     /* create new object medium 
						last argument is volume fraction of 
						quenchers*/


	return;
}		/* -----  end of method ClassMonteCarlo::ClassMonteCarlo  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: Init
 * Description:  Clears up the instance of simulation
 *-------------------------------------------------------------------------------------
 -
 */

void
ClassMonteCarlo::Init_IQ(string sID, double bx, double by, double bz, double iConc_cm2)
{
	int i;
	SimulationID=sID;
	for ( i = 0; i < TN; i += 1 ) {
		PL[i]=0;
	}
	medium = new ClassMedium(bx, by, bz);     /* create new object medium */
	medium->PlaceInterfaceQuenchers(iConc_cm2);

	return;
}		/* -----  end of method ClassMonteCarlo::ClassMonteCarlo  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: Init2 
 * Description:  Clears up the array PL, but keeps the medium 
 *-------------------------------------------------------------------------------------
 -
 */

void
ClassMonteCarlo::Init2(string sID) 
{
	int i; 
	SimulationID=sID;	
	for ( i = 0; i < TN; i += 1 ) {
		PL[i]=0;
	}

	return;
}		/* -----  end of method ClassMonteCarlo::ClassMonteCarlo  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: Simulation
 * Description:  Public method of the MonteCarlo class
 *
 * MaxGen is the maximum number of generations. It determines the duration of 
 * 	the calculation
 *--------------------------------------------------------------------------------------
 */

unsigned int	
ClassMonteCarlo::Simulation ( unsigned int MaxGen, double* ret_dL2)
{
	unsigned int i;
	unsigned int radiative = 0; /* number of excitons that decayed radiatively */
	double dL=0, dL2=0, dx=0, dx2=0, average_dL2=0, average_dL=0, average_dx2=0, average_dx=0;
	string filename;

	cout << "\n";
	cout	<< "#Gen.\tsqrt <dL^2>\tsqrt <dx^2>\t<dL>\t\t<dx>\n" ;
        cout << "------------------------------------------------------------------\n";
	for ( i = 0; i < MaxGen; i += 1 ) {
		cout << " " << i; // Gen. #
		radiative+=SimulationG1(dL2,dx2, dL,dx);	
		average_dL+=dL;
		average_dL2+=dL2;
		average_dx+=dx;
		average_dx2+=dx2;
	}
	average_dL=average_dL/MaxGen;
	average_dL2=average_dL2/MaxGen;
	average_dx=average_dx/MaxGen;
	average_dx2=average_dx2/MaxGen;
        cout << "------------------------------------------------------------------\n";
	cout << "av.\t" << average_dL2 << "\t" << average_dx2 << "\t" << average_dL << "\t" << average_dx << "\n";

	filename = SimulationID + ".txt";
	SavePL(PLtime, PL, "latest.txt");	
	SavePL(PLtime, PL, filename);	
	*ret_dL2 = average_dL2;
	return radiative;
}		/* -----  end of method ClassMonteCarlo::Simulation  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: SimulationG1
 * Description:  Runs simulatino for one generation of excitons
 * 		 Returns number of excitons that decayed radiatively
 * 		 and the average displacements dL and dx in 3D and 1D, respectively
 *--------------------------------------------------------------------------------------
 */


unsigned int	
ClassMonteCarlo::SimulationG1 ( double& dL2, double& dx2, double& dL, double& dx )
{
	int i;
	int t;          /* time variable goes from 0 to TN-1  */
	double time;    /* time in ps starting from simulation start used when dT is variable */
	ClassExciton e[GEN]; /* generation of excitons */
	unsigned int radiative= 0; /* number of exctions that decayed radiatively */
	unsigned int PLg[TN]; /* photoluminescence of this generation */
	double r;    /* biexponential decay */
	double delim1 = ClassExciton::a1 * ClassExciton::tau1 / (ClassExciton::a1 * ClassExciton::tau1 + ClassExciton::a2 * ClassExciton::tau2 + ClassExciton::a3 * ClassExciton::tau3) ; 
	double delim2 = (ClassExciton::a1 * ClassExciton::tau1 + ClassExciton::a2 * ClassExciton::tau2 )/ (ClassExciton::a1 * ClassExciton::tau1 + ClassExciton::a2 * ClassExciton::tau2 + ClassExciton::a3 * ClassExciton::tau3) ; 

	double microHop; // if ClassExciton::hopsize is larger than ClassQuencher::radius then do a couple of smaller hops per cycle
	int k, microHopNum; // k is iterator; microHopNum is the number of small hops per cylce
	int microKtime; // is microHopNum^2
	double hop_QuencherRatio;

//	double hopsize_hot;  /* exciton relaxation parameters */
	double hopsize_cold;
	double hopsize_curr;
	bool thermal_relaxation;

	thermal_relaxation = config->pBool("thermal_relaxation");
//	hopsize_hot = config->pDouble("TRel_hopsize_hot")*sqrt(dT);
	hopsize_cold = ClassExciton::hopsize;
//	TRel_hotFrac= config->pDouble("TRel_hotFrac")/100; 
	hop_QuencherRatio = config->pDouble("hop-QuencherRatio");


	microHopNum = ceil ( ClassExciton::hopsize / (ClassQuencher::radius * hop_QuencherRatio) );
	microHop = ClassExciton::hopsize / (microHopNum*1.0);
	microKtime = microHopNum * microHopNum; 

 	// cout << "microHopNum\t= " << microHopNum << "\n";
	// cout << "microHop\t= " << microHop << "\n";
	// cout << "microKtime\t= " << microKtime << "\n";

//	cout	<< "I am in SimulationG1()" << endl;
	time = 0;
	for ( i = 0; i < TN; i += 1 ) {
		PLg[i]=0;
	}

/* 	initialize excitons	
 */

	for ( i = 0; i < GEN; i += 1 ) 
	{
//
		e[i].nexti = i+1;               /* this is needed to exclude inactive exctions from the loop */
		e[i].previ = i-1;

		r=Rnd.dRnd(); /* biexponential decay */
		if (r<=delim1)
		{
			e[i].lifetime = -log(Rnd.dRnd())*ClassExciton::tau1;

		}else if (r<=delim2){
			e[i].lifetime = -log(Rnd.dRnd())*ClassExciton::tau2;
		}else{
			e[i].lifetime = -log(Rnd.dRnd())*ClassExciton::tau3;
		}

		if (thermal_relaxation == true) 
		{
			e[i].hopsize0 = Rnd.dNormal(ClassExciton::TRel_M, ClassExciton::TRel_S);
			if (e[i].hopsize0 <=0) 
			{
				cout << "Warning: hopsize0 = " << e[i].hopsize0 << "\n";
				e[i].hopsize0 = ClassExciton::hopsize;	
			}
		}else{
			e[i].hopsize0 = 0;
		}

	
		do {
			e[i].x1=e[i].x2=e[i].x=Rnd.dRnd()*medium->X; /* initial position */
			e[i].y1=e[i].y2=e[i].y=Rnd.dRnd()*medium->Y;
			e[i].z1=e[i].z2=e[i].z=Rnd.dRnd()*medium->Z;

//			e[i].x1=e[i].x2=e[i].x=medium->X*(1/2+Rnd.dRnd()/5); /* initial position */
//			e[i].y1=e[i].y2=e[i].y=medium->Y*(1/2+Rnd.dRnd()/5);
//			e[i].z1=e[i].z2=e[i].z=medium->Z*(1/2+Rnd.dRnd()/5);
		}while (medium->GridQuenching(medium->grid,&e[i]));
	}

	e[GEN-1].nexti = 0; /* close the circle */ 
	e[0].previ = GEN - 1;

/*	main simulation loop 	
 */
	t = 0;
	i = 0;
	time = 0;

	hopsize_curr = hopsize_cold;
	microHopNum = ceil ( hopsize_curr / (ClassQuencher::radius * hop_QuencherRatio) ); // recalculate hopsize
	microHop = hopsize_curr / (microHopNum*1.0);
	microKtime = microHopNum * microHopNum; 
	

	do { // loop though the time
		PLtime[t]=time;
		
		do { // loop through the exitons
			if (e[i].nexti==e[i].previ) break; /* neglect the last exciton */

			
			
/*			check for radiative decay 
 */


			if ( e[i].lifetime<=time + dT )
			{
				PLg[t]++; /* collect light */
				radiative++; /* count emitted photons */
				e[i].active=false; /* deactivate exciton */
				e[e[i].previ].nexti = e[i].nexti; /* exclude this exciton from the circle */
				e[e[i].nexti].previ= e[i].previ; /* exclude this exciton from the circle */
				i=e[i].nexti;
				continue;
			}
/* 			move exciton
 */

			if ( thermal_relaxation == true ) 
			{
				//hopsize_curr =  (e[i].hopsize0 - hopsize_cold)*2.7*exp(-exp(time/TRel_time ) ) + hopsize_cold;
				// hopsize_curr =  (e[i].hopsize0 - hopsize_cold)*exp(-time/TRel_time ) + hopsize_cold;
				//hopsize_curr =  sqrt ( (e[i].hopsize0*e[i].hopsize0 - hopsize_cold*hopsize_cold)*exp(-time/TRel_time )/pow( (time/TRel_time+1.0) ,2) + hopsize_cold*hopsize_cold );
				//hopsize_curr =  sqrt ( (e[i].hopsize0*e[i].hopsize0 - hopsize_cold*hopsize_cold)*exp(-time/TRel_time )/pow( (time+1) ,6) + hopsize_cold*hopsize_cold );
				hopsize_curr =  sqrt ( (e[i].hopsize0*e[i].hopsize0 - hopsize_cold*hopsize_cold)*exp(-time/ClassExciton::TRel_time ) + hopsize_cold*hopsize_cold );
				microHopNum = ceil ( hopsize_curr / (ClassQuencher::radius * hop_QuencherRatio) ); // recalculate hopsize
				microHop = hopsize_curr / (microHopNum*1.0);
				microKtime = microHopNum * microHopNum; 
			}		

			//cout << time << "\t" << hopsize_curr << "\n";	

			for (k=0;  k < microKtime; k++)
			{
				e[i].hop3(medium->X, medium->Y, medium->Z, microHop); /* args: box size to fulfill the boundary conditions  */
/* 				check for quenching
 */
//				if ((e[i].z>0) && (e[i].z<1) ) /* interface quenching close to x=0 */
//				if (medium->TwoInterfaceQuenching(e[i], medium->X, medium->Y, medium->Z ) ) 
							/* interface quenching close to x=0 */
	
				if ( medium->GridQuenching(medium->grid,&e[i]))
				{
					//cout << "Quenching!" << "\n";
					e[i].active=false; /* deactivate exciton */
					e[e[i].previ].nexti = e[i].nexti; /* exclude this exciton from the circle */
					e[e[i].nexti].previ= e[i].previ; /* exetlude this exciton from the circle */
					//i=e[i].nexti; // this will be executed before "while"
					break; // exit (for k=0) cycle
				}
	
			}



			i=e[i].nexti;
		} while ( e[i].nexti > e[i].previ );		/* -----  end do-while  ----- */
		i=e[i].nexti;
		

		t+=1;
		time+=dT;

	}while ( (t<TN) && (e[i].nexti!=e[i].previ) ); // simulation time finished or excitons finished 
				/* -----  end do-while  ----- */

	{ // calculate <dL^2> and <dx^2>
		dL = dL2 = 0;
		dx = dx2 = 0;
		double x,y,z;
		for ( i = 0; i < GEN; i += 1 ) 
		{

			x = e[i].x1 - e[i].x2;
			y = e[i].y1 - e[i].y2;
			z = e[i].z1 - e[i].z2;
			dx2+=x*x;
			dx+=fabs(x);
			dL2+=x*x + y*y + z*z;
			dL+=sqrt(x*x + y*y + z*z);
		}
		dL2 = sqrt(dL2/GEN);
		dL = dL/GEN;
		dx2= sqrt(dx2/GEN);
		dx= dx/GEN;
		cout << "\t" << dL2 <<  "\t" << dx2 << "\t" << dL << "\t" << dx << "\n"; 
	}

/*	copy PLg to PL 
 */

	for ( i = 0; i < TN; i += 1 ) {
		PL[i]+=PLg[i];
	}
	return radiative;
}		/* -----  end of method ClassMonteCarlo::SimulationG1  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMonteCarlo
 *      Method:  ClassMonteCarlo :: SavePL
 * Description:  Saves array of PL data to hard drive
 *--------------------------------------------------------------------------------------
 */



	void
ClassMonteCarlo::SavePL ( double * arr_t, unsigned * arr, string filename )
{
	int i;
	ofstream datafile;
	string fname="";
	time_t rawtime;
	struct tm * timeinfo;
	char time_str[100];
	int save_max; // maximum number of points to be saved 

	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (time_str, 100, "%H.%M.%S-", timeinfo);

	if (filename.compare("latest.txt")==0) // returns 0 if strings are equal 
	{
		fname = config->pString("outputFolder");
		fname += "/";
	}else{
		fname += global_folder; 
	       	fname += "/";
		fname += time_str;	
	}
	fname += filename;
	
	if (config->pInt("save_max")< TN)
	{
		save_max = config->pInt("save_max");
	}else{
		save_max = TN;
	}
	

	datafile.open (fname.c_str(),ios::out);
	for ( i = 0; i < save_max; i += 1 ) {

//		cout << i*dT << "\t" <<  arr[i] << endl;	
		if ( (arr_t[i]>0) || (i==0) ) {
			datafile << arr_t[i] << "\t" <<  arr[i] << endl;	

		}
	}
	datafile.close();

	return ;
}		/* -----  end of method ClassMonteCarlo::SavePL  ----- */


