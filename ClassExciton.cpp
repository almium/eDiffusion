/*
 * =====================================================================================
 *
 *       Filename:  ClassExciton.cpp
 *
 *    Description:  Contains declarations of ClassExciton
 *
 *        Version:  1.0
 *        Created:  1/20/2011 10:47:14
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

#include 	"PreProcessorOptions.h"

extern ClassRandom Rnd;        /* gloabal instance of random generator */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassExciton
 *      Method:  ClassExciton :: ClassExciton
 * Description:  Constructor of exctions
 *--------------------------------------------------------------------------------------
 */
ClassExciton::ClassExciton ()
{
	lifetime = tau1; 
	x1=y1=z1=x2=y2=z2=x=z=y=0;
	active = true;
	return ;
}		/* -----  end of method ClassExciton::ClassExciton  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassExciton
 *      Method:  ClassExciton :: hop
 * Description:  Moves exction for hopsize in random direction (slow!!!!)
 *--------------------------------------------------------------------------------------
 */

	void
ClassExciton::hop (double bx, double by, double bz)
{
	double theta, phi;

	theta = Rnd.dRnd()*PI-PI/2;
	phi = Rnd.dRnd()*2*PI;
	x = x + hopsize * cos(theta) * cos(phi);
	y = y + hopsize * cos(theta) * sin(phi);
	z = z + hopsize * sin(phi);

	if (x>bx) x-=bx;
	else if (x<0) x+=bx;
	if (y>by) y-=by;
	else if (y<0) y+=by;
	if (z>bz) z-=bz;
	else if (z<0) z+=bz;

	return ;
}		/* -----  end of method ClassExciton::hop  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassExciton
 *      Method:  ClassExciton :: hop2
 * Description:  Moves exction for hopsize in random direction
 *--------------------------------------------------------------------------------------
 */

	void
ClassExciton::hop2 (double X, double Y, double Z)
{
	double dx, dy, dz, metric; 

	dx = Rnd.dRnd()-0.5;
	dy = Rnd.dRnd()-0.5;
	dz = Rnd.dRnd()-0.5;
	//dy=dz=0; // debug - 1D diffusion
	metric = sqrt (dx*dx + dy*dy + dz*dz);

	dx=hopsize*dx/metric;
	dy=hopsize*dy/metric;
	dz=hopsize*dz/metric;

	x+=dx;
	y+=dy;
	z+=dz;

	x2+=dx;
	y2+=dy;
	z2+=dz;

	if (x>=X) x-=X;     /* periodic boundary conditions */
	else if (x<0) x+=X;
	if (y>=Y) y-=Y;
	else if (y<0) y+=Y;
	if (z>=Z) z-=Z;
	else if (z<0) z+=Z;

	return ;
}		/* -----  end of method ClassExciton::hop  ----- */

/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassExciton
 *      Method:  ClassExciton :: hop3
 * Description:  Moves exction for hopsize in random direction
 *--------------------------------------------------------------------------------------
 */

	void
ClassExciton::hop3 (double X, double Y, double Z, double hop)
{
	double dx, dy, dz, metric; 

	dx = Rnd.dRnd()-0.5;
	dy = Rnd.dRnd()-0.5;
	dz = Rnd.dRnd()-0.5;
	//dy=dz=0; // debug - 1D diffusion
	metric = hop*1.0/sqrt (dx*dx + dy*dy + dz*dz);



	dx=dx*metric;
	dy=dy*metric;
	dz=dz*metric;

	x+=dx;
	y+=dy;
	z+=dz;

	x2+=dx;
	y2+=dy;

	if (x>=X) x-=X;     /* periodic boundary conditions */
	else if (x<0.0) x+=X;
	if (y>=Y) y-=Y;
	else if (y<0.0) y+=Y;

#if !NO_Z_PERIODIC
	if (z>=Z) z-=Z;
	else if (z<0.0) z+=Z;
	z2+=dz;
#else
	if (z>=Z) z=Z-0.1;
	else if (z<0.0) z=0.0;
	else z2+=dz;
#endif

	return ;
}		/* -----  end of method ClassExciton::hop3  ----- */



