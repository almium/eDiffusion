/*
 * =====================================================================================
 *
 *       Filename:  medium.cpp
 *
 *    Description:  Contains definitions of classes related to medium
 *
 *        Version:  1.0
 *        Created:  1/20/2011 10:49:46
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
#include	<algorithm>                     /* needed for sort () */
#include	<iterator>                      /* probably needed for sort() */
#include	"eDiffusion.h"
#include	"randomc.h"
#include 	"EasyBMP.h"
#include	"config.h"
#include	"log.h"

extern ClassRandom Rnd;        /* gloabal instance of random generator */
extern Config *config;		/* global instance of config class */
extern string global_folder; 	// where all the output is saved for the whole execution
			// ex: output/YYYYMMDD-HH.MM-experiment_name






/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: ClassMedium
 * Description:  Constructor of ClassMedium
 *--------------------------------------------------------------------------------------
 */
ClassMedium::ClassMedium (double bx, double by, double bz) /* this constructor does not 
							      allocate memory for Bool3D*/
{
	X = bx;
	Y = by;
	Z = bz;
	
	qn=0;                                   /* number of quencher molecules */

	return ;
}

ClassMedium::ClassMedium (double bx, double by, double bz, double frac ) /* this constructor allocates
									memory for Bool3D*/
{
	int i,j,k;
	int qRadiusN;                            /* quencher radius expressed in grid points */
	double qR, qR2;                         /* quencher radius, quencher radius squared */
	int counter=0;	
	/* set the box size */
	X = bx;
	Y = by;
	Z = bz;

	/* Discretize the space	*/
	Nx = X/dX +0.5;                         /* +0.5 is needed because 10/0.1 returns 99! */
	Ny = Y/dX +0.5;
	Nz = Z/dX +0.5;
	cout	<< "(X, Y, Z)=(" << X << ", "<< Y << ", "<< Z << ")" << endl;
	cout << "dX = " << dX << endl;
	cout	<< "(Nx, Ny, Nz)=(" << Nx << ", "<< Ny << ", "<< Nz<< ")" << endl;
	//exit (0);

	numQ=0;
	grid = new Bool3D(Nx,Ny,Nz);             /* allocate memory for medium.
						   memory gets free at the destructor */

	/* generate quencher distribution */

	/* calculate number of quenchers and allocate memory */
	qn = (int) X*Y*Z * frac/ (PI * pow(ClassQuencher::radius,3)*4.0/3.0 );	
//	qn = 5000; // debug
	cout	<< "number of quencher molecules: " << qn  << endl;
	


	q = new ClassQuencher[qn];
	
	cout << "memory have been allocated for quenchers" << "\n";

	/* put quenchers into the grid */

	/* make a set of grid points that are overlapping with a quencher molecule
	 * in case if molecule is positioned at coordinates (0,0,0).
	 * Exciton radius is also included.
	 * negative values are okay */

	qR =ClassQuencher::radius + ClassExciton::radius; 
//	qR2=(qR+dX*3.6/5.0)*(qR+dX*3.6/5.0)/(dX*dX);    /* quencher radius squared divided by dX*dX */
        						/* extra factor +dX*a/b is to make sure that volume fraction is right */	
	qR2=qR*qR/(dX*dX);
	qRadiusN =qR/dX + 1;                    /* +1 to be sure */

	//cout << "\nqR2 " << qR2 << "\n";	

	counter = 0;
	for ( i = -qRadiusN; i <= qRadiusN; i += 1 ) {
		//cout << "\n\ni = " << i << "\t gi.size()" << gi.size() << "\n";
		for ( j = -qRadiusN; j <= qRadiusN; j += 1 ) {
			for ( k = -qRadiusN; k <= qRadiusN; k += 1 ) {
				if ( /* any corner of the cube i,j,k crosses the sphere with radius qR */
						( i*i + j*j + k*k <= qR2 )
						|| ( (i+1)*(i+1) + j*j + k*k <= qR2 )
						|| ( (i+1)*(i+1) + (j+1)*(j+1) + k*k <= qR2 )
						|| ( (i+1)*(i+1) + (j+1)*(j+1) + (k+1)*(k+1)  <= qR2 )
						|| ( i*i + (j+1)*(j+1) + (k+1)*(k+1)  <= qR2 )
						|| ( i*i + j*j + (k+1)*(k+1)  <= qR2 )
						|| ( i*i + (j+1)*(j+1) + k*k  <= qR2 )
						|| ( (i+1)*(i+1) + j*j + (k+1)*(k+1)  <= qR2 )

				   )
//				if ( /* center of cell i,j,k is within the sphere of radius qR */
//						( (i+copysign(0.5,i))*(i+copysign(0.5,i)) 
//						  + (j+copysign(0.5,j))*(j+copysign(0.5,j))
//						  + (k+copysign(0.5,k))*(k+copysign(0.5,k)) ) <=qR2
//				   )

				{
					gi.push_back(i); /* memorize coordinates i,j,k */
					gj.push_back(j);
					gk.push_back(k);
					counter++;
					//cout << 'x';
				}else{
					//cout << '-';
				}
			}
			// cout << "\n";
		}
	}

	cout << "quencher includes " << counter << " grid cells \n";

	/* Actually place quenchers to the grid */
	if (config->pInt("perCluster")>1){
		PlaceClustersQ(config->pInt("perCluster"));    /* argument is number of quenchers per cluster */
	}else{
		PlaceHomogeneousQ();
	}

//	ReadMediumBMP("medium.P3HT-ZnO");	


	cout	<< "quenchers have been placed to the grid" << endl;
	
	/* Save 2D BMP image of the grid nz=Nz/2 */

	{	
		char fname[128];
		sprintf(fname,  "Grid Image.Volume frac %2.5f.bmp",frac);
		Save2DgridBMP(Nz/2,fname);
	}

	return ;
}		/* -----  end of method ClassMedium::ClassMedium  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: ~ClassMedium
 * Description:  Destructor of Class medium
 *--------------------------------------------------------------------------------------
 */


ClassMedium::~ClassMedium ( )
{
	delete [] q;                            /* free memory */
	delete grid;
	return ;
}		/* -----  end of method ClassMedium::~ClassMedium  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: PlaceHomogeneousQ
 * Description:  Places the quencher molecules to the grid in homogeneous manner
 *--------------------------------------------------------------------------------------
 */

bool
ClassMedium::PlaceHomogeneousQ ()
{
//	unsigned i; 
	int i, a;                                /* iterators */
	int cnt;                                /* counter */
	bool goodPos;                           /* good position flag */
	int pos;                     /* indexes for sorting */
	ClassQuencher quencher;                 /* is needed to resort the array */

	
	if (qn>0){                              /* position the first quencher */
		q[0].RandPos(X,Y,Z);
		SetOnGrid(&q[0]);
	}
	
	for ( a = 1; a < qn; a += 1 ) 
	{
		if (a%1000 == 0)
			cout << "total quenchers: " << qn << "\t\tcurrent quencher: " << a << "\n";

		/* randomly set a quencher and check the position validity 
		 * quenchers should not overlap in Homogeneous Placement*/

		cnt=0;
		do { /* repeat until find a good position */
			cnt++;	
			//cout << "try position for quencher " << a << "\n";
			goodPos = q[a].RandPos(X,Y,Z); /* random position */
			
			pos = CheckValidQposition (a-1,&q[a]);

			if (pos == -2)
			{
				goodPos=false;
				continue;
			}

		}while (goodPos==false);


		/* if position is valid, resort the array */	
		if (goodPos){
			/* mark this quencher on the grid */
			SetOnGrid (&q[a]);
			quencher = q[a];

			for ( i = a; i > pos+1; i += -1 ) { /*  resort the array */
				q[i]=q[i-1];
			}
			q[pos+1]=quencher;

		}

	
		if (cnt > 100 )
			cout << cnt << " false positions for quencher " << a << "\n";

	}


//	ViewGridStats();                        /* shows statistics on the grid */

	return true;
}		/* -----  end of method ClassMedium::PlaceHomogeneousQ  ----- */








/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: PlaceClustersQ
 * Description:  Places clustered quenchers into the Grid
 *--------------------------------------------------------------------------------------
 */


bool
ClassMedium::PlaceClustersQ (int num)     /* num is number of quenchers per cluster */
{
	int cn;  /* number of clusters */
	int i, k, a ;   /* iterators */
	int pos; /* index position which is returned by CheckValidQposition() */
	//int cnt; /* counter */
	int anum;  /* actual number of quenchers per cluster */
	bool cPosFlag;   /* flag of a good cluster position */
	ClassQuencher centerQ,quencher; /* variables to test random quencher position */



	anum = ReadCrystal (); // reads the coordinates of num quenchers in a cluster 
			// from external file and stores them in vecotrs
			// clusterX, clusterY, clusterZ.

	if (anum < num) // num is larger than number of records in the file
	{
		
		cout	<< "\n* WARNING: the number of quenchers per cluster exceeds the number of lines in the file with crystal description.\n" ;
	
	}

	cn = qn/num; // number of clusters


	// put the first cluster on the grid
	a=0;
	
	centerQ.RandPos(X,Y,Z); /* randomly select coordinates of the center of first cluster */
	//centerQ.SetCoordinates (5,5,5,X,Y,Z); // debug

	cout << "I am going to place the first cluster ....  ";
	if ( qn > cn * num)  /* check if all the clusters should be of the size num */
	{
		/* place a smaller cluster*/
		cout << "of " <<  ( qn- cn * num )<< " elements ... ";
		a = SetCluster(-1, &centerQ, qn - cn * num );

	}else{
		/* place a full size cluster */
		a = SetCluster(-1, &centerQ, num );
		cn--;	
	}
	cout << "OK\n";

	// loop through the rest of the clusters
		
	for ( k = 1; k <= cn; k += 1 ) {
		//cout << "place cluster number " << k << "\n";
		do // find a good position for cluster k 
		{
			cPosFlag=centerQ.RandPos(X,Y,Z); /* randomly choose the center of the cluster */
			
			for ( i = num-1; i >= 0; i += -1 ) /* check if the randomly chosen position is good */
			{
				quencher.SetCoordinates(centerQ.x + clusterX[i], 
							centerQ.y + clusterY[i], 
							centerQ.z + clusterZ[i], 
							X,Y,Z);
				pos=CheckValidQposition(a,&quencher);  	
				if (pos ==-2)
				{
					//cout << "false position cluster " << k << "\n";
					cPosFlag=false;
					break;
				}
			}
		}while (cPosFlag==false);

		/* at this moment centerQ should have a valid coordinates */

		/* insert quenchers into the q[] array */
		/* mark this quencher on the grid */
		a = SetCluster(a,&centerQ, num);

	}

	



	return true;
}		/* -----  end of method ClassMedium::PlaceClustersQ  ----- */




/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: FindQposition
 * Description:  finds the index in the ordered array q[]
 * 		after whcih new element q[a] must be inserted
 *--------------------------------------------------------------------------------------
 */

int 
ClassMedium::FindQposition ( int a , ClassQuencher* quencher)
{
	int left, right, mid;

	left=0;
	right=a;
	// cout << "\n*  a = " << a << "; new value = " << quencher->dCenter << "\n";
	
	if (a == -1 ) return -1;  /* empty array, place new element at a=-1+1=0 */

	if ( quencher->dCenter < q[left].dCenter ) {
		// check if this point  with very
		// right element of the sorted array q[i]

		mid = -1;
	}else if ( quencher->dCenter >= q[right].dCenter ){
		// then check this point with very
		// left element of the sorted array q[i]
		mid = a; 
	}else{
		// quencher quencher should be in between of q[left] 
		// and q[right]
		/* find position of quencher quencher in sorted array */

		while ( right-left>1)
	       	{
			mid = (right+left)/2;
			// cout << "\tleft = " << left << "; right = " << right << "; mid=" << mid<< "\n";
			if (quencher->dCenter >= q[mid].dCenter)
			{
				left = mid;
			}else{
				right = mid;
			}			
		} // quencher should be at positon mid+1 
		mid = left;
	}
	return mid;
}		/* -----  end of method ClassMedium::FindQposition  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: CheckValidQposition
 * Description: Check if the quencher overlaps 
					with others in the array q[0..a]; 
					returns -2 if position is invalid, otherwise  
					returns index in array q[] after which 
					this quencher should be inserted
 *--------------------------------------------------------------------------------------
 */

int
ClassMedium::CheckValidQposition (int a, ClassQuencher* quencher)
{
	int i;                                /* iterators */
	bool goodPos;                           /* good position flag */
	int pos;                     /* index in q[] array after which the element quencher 
					should be inserted */

	goodPos=true;
	pos = FindQposition(a, quencher);

	/* check validity of the quencher */

	/* first on the right side */ 			
	i = pos+1;
	while ( i<=a )
	{
		if (q[i].dCenter-quencher->dCenter > ClassQuencher::radius*2){
			break;
		}

		if (quencher->distance(q[i], X, Y, Z)<=ClassQuencher::radius*2)
		{
			goodPos=false; /* ovelapping quenchers */
			break;
		}
		i++;
		
	}


	/* then on the left side */
	i = pos;
	if (goodPos)
		while (i>=0)
		{
			if (quencher->dCenter-q[i].dCenter > ClassQuencher::radius*2){
				break;
			}

			if (quencher->distance(q[i], X, Y, Z)<=ClassQuencher::radius*2)
			{
				goodPos=false; /* overlapping quenchers */
				break;
			}
			i--;
		}	

	if (goodPos){
		return pos;
	}else{
		return -2;
	}



}		/* -----  end of method ClassMedium::GetValidQposion  ----- */






/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: SetOnGrid
 * Description:  sets a quencher on the grid
 *--------------------------------------------------------------------------------------
 */

void
ClassMedium::SetOnGrid ( ClassQuencher  *quencher  )
{
	int i;
	int nx,ny,nz;                           /* quencher positon in grid numbers */

	nx=quencher->x/dX;
	ny=quencher->y/dX;
	nz=quencher->z/dX;

//	{
//		int i,j;
//	       
//		if (grid->get(nx,ny,nz)==true ){
//			cout << "\n* WARNING: wrong quencher center on the grid!" << " "<< "(nx, ny, nz)=(" << nx <<"," << ny << "," << nz <<")\n" ;
//			{
//				for ( i = nx-10; i < nx+10; i += 1 ) 
//				{
//					
//					for ( j = ny-10; j < ny+10; j += 1 ) 
//					{
//						if (grid->get(i,j,nz)==true){
//							cout << "x";
//						}else{
//							cout << "-";
//						}	
//					}
//					cout << "\n";
//				}
//			}	
//		}	
//
//	}


	for ( i = 0; i < (int) gi.size() ; i += 1 )
       	{
		numQ++;
		grid->set(nx+gi[i], ny+gj[i], nz+gk[i], true); 
								/* Class Bool3D takes care 
			       					about periodic boundary conditions */
		//cout	<< "(nx+gi[i], ny+gj[i], nz+gk[i])=(" << nx+gi[i] <<"," << ny+gj[i] << "," << nz+gk[i] <<")\n";  
		// cout	<< "(gi[i], gj[i], gk[i])=(" << gi[i] <<"," << gj[i] << "," << gk[i] <<")\n";  
	}

	return ;
}		/* -----  end of method ClassMedium::SetOnGrid  ----- */




/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: SetCluster
 * Description:  sets a cluster on the grid and updates the q[] array;
 * 		returns total number of quenchers that are set on the grid
 *--------------------------------------------------------------------------------------
 */


int
ClassMedium::SetCluster ( int a, ClassQuencher* centerQ, int num ) 
{

	int i,j;
	int pos;    /* index of a new quencher to be inserted into q[] array */ 
//	int a;          /* index of the last element in the q[] array */
	ClassQuencher quencher;

	/* check if cluster size is not too big */
	if (num > (int) clusterX.size())
	{
		cout << "* WARNING: attemptint to create cluster, which contains more quenchers than number of lines in the original crystal file \n";
		cout << "\tnumber of lines read from the file = " << (int) clusterX.size() << "\t";
		cout << "attempting cluster size = " << num << "\n";
		cout << "\twill use cluster size ";
	}

	


	for ( i = 0; i < num; i += 1 ) 
	{
		a++;
		q[a].SetCoordinates(	centerQ->x + clusterX[i], 
					centerQ->y + clusterY[i], 
					centerQ->z + clusterZ[i], 
					X,Y,Z);
		
		//cout << "setting quecher" << " "<< "(clusterX, clusterY, clusterZ)=(" << clusterX[i] <<"," << clusterY[i] << "," << clusterZ[i] <<")\n" ;

	
		SetOnGrid (&q[a]); /* mark this quencher on the grid */
		
		pos = FindQposition(a-1,&q[a]); /*  place it in the ordered array */

		// resort the array
		quencher = q[a];

		for ( j = a; j > pos+1; j += -1 ) {
			q[j]=q[j-1];
		}
		q[pos+1]=quencher;
	
	}

	return a; // returns the index of last quencher in the array
}		/* -----  end of method ClassMedium::SetClusterOnGrid  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: ReadMediumBMP
 * Description: reads medium from BMP files, stored 
		      in specified folder. Files should be named
		     1.bmp, 2.bmp, 3.bmp etc. 
 *--------------------------------------------------------------------------------------
 */


bool
ClassMedium::ReadMediumBMP (const char* folder )
{
	int i,j,k,m;
	char buf[256];
	BMP slice;

	

	
	for ( k = 1; k <= 28 ; k += 1 ) {

		sprintf(buf, "%s/%i.bmp", folder, k);
		slice.ReadFromFile(buf);
	
		
		for ( i = 0; i < slice.TellWidth(); i += 1 ) {
			for ( j = 0; j < slice.TellHeight(); j += 1 ) {
				if (slice(i,j)->Red==0 ){ // all three channels are 0 to make it black. 
							  // no need to check them all though
					for ( m = 0; m < (int) gi.size() ; m += 1 )
				       	{
						numQ++;
						grid->set(i+gi[m], j+gj[m], k+gk[m], true); 
					}
	
					grid->set(i, j, k, true); 
				}
			}
		}

	}




	return true;
}		/* -----  end of method ClassMedium::ReadMediumBMP  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: GridQuenching
 * Description:  Checks if exciton on the grid; returns true if exciton has to be quenched
 *--------------------------------------------------------------------------------------
 */


bool
ClassMedium::GridQuenching ( ClassExciton* e)
{
	int nx, ny, nz;

	nx=e->x/dX;
	ny=e->y/dX;
	nz=e->z/dX;


//	cout	<< "(nx, ny, nz)=(" << nx <<"," << ny << "," << nz <<")\t"  
//		<< "(Nbx, Nby, Nbz)=(" << grid->Nbx <<"," << grid->Nby << "," << grid->Nbz <<")" <<endl;
	return grid->get(nx,ny,nz);

}		/* -----  end of method ClassMedium::GridQuenching  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: TwoInterfaceQuenching
 * Description:  Returns true if exciton is close to z=0
 *--------------------------------------------------------------------------------------
 */
bool
ClassMedium::TwoInterfaceQuenching (ClassExciton& e, double bx, double by, double bz)
{
	if (e.z>0 && e.z<1)
	{
		return true;
	}else{
		return false;
	}

}		/* -----  end of method ClassMedium::TwoInterfaceQuenching  ----- */


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: ReadCrystal
 * Description:  Reads position of quenchers in crystal from an ordered file
 * 		 The quenchers are ordered in the file for their distance from the center
 *--------------------------------------------------------------------------------------
 */

int
ClassMedium::ReadCrystal ()
{
	int i,j,cnt;
	ifstream inpF("crystals/PCBM-triclinic.txt");
	char line[256];
	double buf;


	i=j=cnt=0;
	while ( !(inpF.eof())  )
	{
		cnt++;

		inpF.getline(line, sizeof(line), '\t');
		buf = atof(line);
		clusterX.push_back(buf);
		//cout	<< buf << "\t";

		inpF.getline(line, sizeof(line), '\t');
		buf = atof(line);
		clusterY.push_back(buf);
		//cout	<< buf << "\t";

		inpF.getline(line, sizeof(line), '\t');
		buf = atof(line);
		clusterZ.push_back(buf);
		//cout	<< buf << "\t" ;
	
		inpF.getline(line, sizeof(line), '\n');
		buf = atof(line);
		clusterD.push_back(buf);
		//cout	<< buf << "\n";
	}


	return cnt;

}		/* -----  end of method ClassMedium::ReadCrystal  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: ViewGridStats
 * Description:  Prints statistics of the grid to stdout
 *--------------------------------------------------------------------------------------
 */


void
ClassMedium::ViewGridStats ( )
{// statistics about he grid
	int i, j, k, a;
	int total,marked, total2;
	double ratio;

	total = grid->Nbx * grid->Nby * grid->Nbz;
	marked = total2= 0;

	for ( i = 0; i < grid->Nbx; i += 1 ) {
		for ( j = 0; j < grid->Nby; j += 1 ) {
			for ( k=0;  k < grid->Nbz;  k += 1 ) {
				if (grid->get(i,j,k)==true)
				{
					marked++;
				}	
				total2++;
			}
		}
	}
	//ratio = marked/total;
	ratio = marked*1.0/total2;

	cout << "gi.size()*qn = " << gi.size()*qn << "\n";	
	cout << "numQ = " << numQ << "\n";	
	cout << "marked/totoal = " << marked << "/" << total << "\n"; 
	cout << "marked/totoal2 = " << marked << "/" << total2 << "\n"; 
	cout << "ratio = " << ratio << "\n"; 

	for ( a = 0; a < qn; a += 1 ) {
		cout << a << ". \t"<< q[a].dCenter << "\n";
	}
	return ;
}		/* -----  end of method ClassMedium::ViewGridStats  ----- */




/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassMedium
 *      Method:  ClassMedium :: Save3DgridBMP
 * Description:  saves 2D cut of the grid (0..Nx,0..Ny, zindex) 
 *		 to bmp file
 *--------------------------------------------------------------------------------------
 */


void
ClassMedium::Save2DgridBMP ( int zindex, char* filename )
{
	int i,j;
	string fname="";
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[100];
	
	BMP image;
	RGBApixel pixel_dark,pixel_white;

	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (buffer, 100, "%H.%M.%S-", timeinfo);

	
	fname += global_folder; 
       	fname += "/";
       	fname += buffer;
       	fname += filename;

	
	/* setup the BMP file */	
	image.SetSize(Nx,Ny);	
	image.SetBitDepth(1);

	/* set up the colors */
	pixel_dark.Red=0;
	pixel_dark.Green=0;
	pixel_dark.Blue=0;

	pixel_white.Red=127;
	pixel_white.Green=127;
	pixel_white.Blue=127;
	

	/* raster image */
	for ( i = 0; i < Nx; i += 1 ) {
		for ( j = 0; j < Ny; j += 1 ) {
			if ( grid->get(i,j,zindex) == true ) 
			{
				image.SetPixel(i,j,pixel_dark);
			}else{
				// white is default
				//	image.SetPixel(i,j,pixel_white);
			}

		}
	}
	
	
	/* write file */	
	image.WriteToFile(fname.c_str());

	return ;
}		/* -----  end of method ClassMedium::Save3DgridBMP  ----- */

