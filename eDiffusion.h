/*
 * =====================================================================================
 *
 *       Filename:  eDiffusion.h
 *
 *    Description:  header files for exciton diffusion
 *
 *        Version:  1.0
 *        Created:  12/14/2010 14:07:19
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
#include	<fstream>
#include	<vector>
#include	"randomc.h"
#include	"stocc.h"


#define VERSION "eDiffusion.v.2.0.0"

//#define	dT 0.1			/* time discretization, ps */

#define	TN 100000			/* TN*dT = tolal simulation time */
#define	GEN 10002		/* number of exciton per generation. "+2": Two excitons are neglected per generation */
#define	PI 3.14159265358979323846			/*  PI */


using namespace std;

void
SaveTable ( double  arr[][4], int n, string filename, bool UseDate );
                                                /* saves a two column table to hard drive */

/*
 * =====================================================================================
 *        Class:  Bool3D
 *  Description:  container of 3D bool array 
 * =====================================================================================
 */
class Bool3D
{
	public:
		Bool3D (int bx,int by, int bz)       /* constructor */
			: Nbx(bx), Nby(by), Nbz(bz)
		{ if(bx > 0 && by > 0 && bz > 0) { vec.resize(0); vec.resize(bx * by * bz, false); } }

		~Bool3D(){vec.clear(); vector <bool>().swap(vec);}; /* destructor */

		const bool operator () (int x,int y, int z) const /* index operator */
		{ return vec[z*Nbx*Nby + y*Nbx + x]; }	

		bool get (int x,int y, int z) /* accesses the grid */
		{ return vec[z*Nbx*Nby + y*Nbx + x]; }	

		void set (int x,int y, int z, bool value); /* set value to x,y,z */
		
		int Nbx, Nby, Nbz;         /* dimensions of the array */
	private:
		vector<bool> vec;               /* vector<bool> takes 1 bit per item, 
						   while array bool[] takes 1 BYTE per item */

}; /* -----  end of class Bool3D  ----- */








/*
 * =====================================================================================
 *        Class:  ClassQuencher
 *  Description:  Container of quencher coordinates
 * =====================================================================================
 */
class ClassQuencher
{
	public:
		ClassQuencher () {x=y=z=dCenter=0; return;};                   /* constructor */
		double x,y,z;
		double dCenter;                 /* distance to the center of the box 
						quenchers are sorted by this parameter
						is used to speed up the quencher placement*/
		static double radius;
		static double ASradius;   /* quencher action sphere radius */

		bool operator< (const ClassQuencher &q2){ /* is needed for sort() */
			return dCenter<q2.dCenter;
		}


		double distance (ClassQuencher& b, double X, double Y, double Z); /* calculate 
						      distance between the center of 
						      this qencher and quencher b. X,Y,Z - box size */
		bool RandPos(double X, double Y, double Z); /* randomize position of a quencher */
		void SetCoordinates (double x2, double y2, double z2, double X, double Y, double Z);
							/* args: new coordinates, box size;
							 * sets new coordinates for a quencher; 
							 * checks for boundary conditions*/	

	protected:
		/* ====================  DATA MEMBERS  ======================================= */

	private:
		/* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class ClassQuencher  ----- */





/*
 * =====================================================================================
 *        Class:  Exciton
 *  Description:  Exciton with its properties
 * =====================================================================================
 */
class ClassExciton
{
	public:
		/* ====================  LIFECYCLE     ======================================= */
		ClassExciton ();                             /* constructor */

		/* ====================  DATA MEMBERS  ======================================= */
		double lifetime; /* exciton lifetime */
		double hopsize0; /* hopsize of a hot exciton */
		double x1,y1,z1; /* original coordinates */
		double x2,y2,z2; /* final coordinates, disregard the periodic boundary conditions. 
				    Needed to calculate <dL^2>, <dx^2>, etc. */
		double x,y,z;   /* current coordinates */
		bool active;    /* determines if exciton is alive */
		static double hopsize; /* hopsize, parameter that we adjust */
		static double TRel_time; /* Time of thermal relaxation */
		static double TRel_M; /* Hopsize of hot excitons (Gaussian median)  */
		static double TRel_S; /* Hopsize of hot excitons (Gaussian sigma)  */
		static double a1; /* weighting */
		static double tau1; /* radiative lifetime of an exciton */
		static double a2; /* weighting */
		static double tau2; /* radiative lifetime of an exciton */
		static double radius; /* exciton radius */
		static double a3; /* weighting */
		static double tau3; /* radiative lifetime of an exciton */
		int nexti,previ; /* array id of next and previous active exciton; */
		/* ====================  OPERATORS     ======================================= */

		void hop(double bx, double by, double bz);  /* hop for distance hopsize */
		void hop2(double bx, double by, double bz);  /* hop for distance hopsize */
		void hop3(double bx, double by, double bz, double hop);  /* hop for distance hop */
	private:
		/* ====================  DATA MEMBERS  ======================================= */

}; /* -----  end of class Exciton  ----- */

/*
 * =====================================================================================
 *        Class:  ClassMedium
 *  Description:  This class contains 
 *  			* 3D bool array for quenchers maps;
 *  			* box size;
 *  			* operations on the medium
 *
 *  medium.cpp
 * =====================================================================================
 */
class ClassMedium
{
	public:
		/* ====================  LIFECYCLE     ======================================= */
		ClassMedium (double bx, double by, double bz); /* does not allocate memroy for Bool3D */
		ClassMedium (double bx, double by, double bz, double frac);         /* constructor 
								bx, by, bz are box size in [nm]
								frac is volume fraction of quenchers*/
		~ClassMedium ();                /* destructor, frees memory for Bool3D */

		void PlaceInterfaceQuenchers ( double iConc_cm2);
		 /* places quenchers at plane at z=0 */

		bool TwoInterfaceQuenching (ClassExciton& e, double bx, double by, double bz);
                                                /* returns true if exciton is close to z=0 */

		bool GridQuenching (Bool3D* grid1, ClassExciton * e);
                                                /* returns true if exciton is quenched on the grid */


		double X;                   /* box size in x-direction [nm], set by Init() to value bsize*/
		double Y;                   /* box size in y-direction [nm], set by Init() to value bsize*/
		double Z;                   /* box size in z-direction [nm], by Init() it is set to bsize*/
		double boxdx;                   /* box discretization, set by Init, typically 0.05 nm 
						   _size/boxdx should be an integer !!! */
		int Nx,Ny,Nz;                          /* bsize/boxdx is used to allocate memory for media */
		int numQ;                       /* number of grid points which are set as ture */
		Bool3D * grid;                /* a pointer to a 3D boolean array. true==quencher; false==no quencher */
		Bool3D * ASgrid;                /* a pointer to a 3D boolean array. Action Sphere grid. */

		static double dX; /* space discritization, read from config */
	protected:
		/* ====================  DATA MEMBERS  ======================================= */

	private:
		/* ====================  DATA MEMBERS  ======================================= */
		ClassQuencher * q;              /* array of quenchers */
		int qn;                         /* total number of quencher molecules in a simulation*/
		vector <int> gi, gj, gk;        /* grid points that overap with a quencher 
						   that is placed at (0,0,0) */
		vector <int> ASgi, ASgj, ASgk;   /* ASgrid points that overap with an action sphere 
						   that is placed at (0,0,0) */

		vector <double> clusterX, clusterY, clusterZ; /* coordinates of 
								 quenchers in a cluster */
		vector <double> clusterD;       /* distance from the center */

		int PixelSphere ( double qR, vector<int>&  vi, vector<int>& vj, vector<int>& vk );
		/* Pixel Sphere sets the list of coordinates of a pixelated sphre of a given radius qR 
		 * to vector vi, vj, vk. Returns number of grid cells. */
		
		bool PlaceHomogeneousQ(); /* places quencher molecules to the grid in 
						   homogeneous manner */
		bool PlaceInterfaceQ(); /* places quenchers at z=0 */

		bool PlaceClustersQ(int num);   /* places clustered quencher molecules 
						   to the grid. Argument is a number of 
						   quenchers per cluster
						    */
		bool ReadMediumBMP(const char* folder); /*  reads medium from BMP files, stored 
						      in specified folder. Files should be named
						     1.bmp, 2.bmp, 3.bmp etc. */


		int CheckValidQposition(int a, ClassQuencher*  quencher); 
					/* Check if the quencher overlaps 
					with others in the array q[0..a]; 
					returns -2 if position is invalid, otherwise  
					returns index in array q[] after which 
					this quencher should be inserted*/
		int FindQposition(int a, ClassQuencher* quencher);  
					     /* finds the index in the ordered array q[0..a] 
					       after which new element quencher must be 
						inserted.  */

		void SetOnGrid(ClassQuencher *quencher); /* sets a quencher on the grid */
		int SetCluster (int a, ClassQuencher* centerQ, int num); 
						/* a is the index of the last element in q[];
						creates a cluster of size num at central position centerQ;
						inserts new elements into q[] array;
						sets the cluster on the grid;
						returns index of the last element in updated q[];
						DOES NOT CHECK overlaps. */

		int ReadCrystal ();            /* reads position of quenchers in a crystal from 
	       					an ordered file. Returns number of lines read*/
		void ViewGridStats();           /* prints statistics of the grid to stdout */
		
		void Save2DgridBMP(Bool3D* grid1, int zindex, char* filename); 
						/* saves 2D cut of the grid (0..Nx,0..Ny, zindex) 
						 * to bmp file */

		

}; /* -----  end of class ClassMedium  ----- */


/*
 * =====================================================================================
 *        Class:  ClassMonteCarlo
 *  Description:  This class will handle the entire MC simulation
 * =====================================================================================
 */
class ClassMonteCarlo
{
	public:
		/* ====================  LIFECYCLE     ======================================= */
		ClassMonteCarlo ();                             /* constructor */
		~ClassMonteCarlo ();                             /* constructor */

		/* ====================  ACCESSORS     ======================================= */

		/* ====================  MUTATORS      ======================================= */

		/* ====================  DATA MEMBERS  ======================================= */
		string SimulationID;            /* is used in file names */
		unsigned PL[TN];            /* resulting photoluminescence */
		double PLtime[TN];            /* resulting photoluminescence */
		unsigned int iDesired;          /* desired integrated PL intensity */
                                                /* determines the duration of calculation */
		unsigned int	
		Simulation(unsigned int MaxGen, double* ret_dL2); 
						/* Max-gen maximum number of exciton generations */
                                                /* determines the duration of calculation;
						 * dL2 returns mean MRS displacement of excitons */
		void Init(string sID, double bx, double by, double bz, double frac);
	       		/* calls Init(); Initializes the medium with fullerens
			 * frac - volume fraction of fullerens 0..1 */
		void Init_IQ(string sID, double bx, double by, double bz, double iConc_cm2);  /* this is for interface quenching */
		void Init2(string sID);  // clears array PL[], but does not touch the medium
		void InitMedium(string sID);  // clears array PL[], update medium 
		static double dT; // time discretization, read from config

	private:
		/* ====================  DATA MEMBERS  ======================================= */

		ClassMedium *medium;            /* medium, in wich excitons diffuse */

			unsigned int	
		SimulationG1(double& dL2, double& dx2, double& dL, double& dx); /* Function that actually performs simulation. 
	       					   Returns number of excitons that decayed radiatively.	
						   and the average displacements dL and dx in 3D and 1D, respectively */
			void 
		SavePL(double * arr_t, unsigned * arr, string filename); /* saves PL decay to a file */
		

}; /* -----  end of class ClassMonteCarlo  ----- */



/*
 * =====================================================================================
 *        Class:  ClassRandom
 *  Description:  Random number generator
 * =====================================================================================
 */
class ClassRandom
{
	public:
		/* ====================  LIFECYCLE     ======================================= */
		ClassRandom ();                             /* constructor */
		~ClassRandom ();                             /* destructor */ 
		double dRnd(); 		/* this function generates 
					   random nubmers between 0 and 1 */
		double dNormal(double center, double sigma); 		/* this function generates 
					   random nubmers between 0 and 1 */
		int iRnd(int max); /* returns random integer 0 <= r <= max */
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		int seed;
		CRandomMersenne* RandGen; // arbitrary seed for random init. will init later
		StochasticLib1* sto; // this is for normal distributions
}; /* -----  end of class ClassRandom  ----- */


