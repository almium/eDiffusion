/*
 * =====================================================================================
 *
 *       Filename:  ClassQuencher.cpp
 *
 *    Description:  Contains implementations of ClassQuencher methods
 *
 *        Version:  1.0
 *        Created:  10/12/2011 09:57:59
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
#include	<ctime>
#include	<math.h>
#include	<algorithm>                     /* needed for sort () */
#include	<iterator>                      /* probably needed for sort() */
#include	"eDiffusion.h"
#include	"randomc.h"
#include	"config.h"
#include	"log.h"

extern ClassRandom Rnd;        /* gloabal instance of random generator */
extern Config *config;		/* global instance of config class */
extern string global_folder; 	// where all the output is saved for the whole execution
			// ex: output/YYYYMMDD-HH.MM-experiment_name


/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassQuencher
 *      Method:  ClassQuencher :: RandPos
 * Description:  Ranomize position of a quencher. 
 *--------------------------------------------------------------------------------------
 */


bool
ClassQuencher::RandPos ( double X, double Y, double Z )
{
	x=X*Rnd.dRnd();
	y=X*Rnd.dRnd();
	z=Z*Rnd.dRnd();
	dCenter= sqrt (  pow(X/2-x,2) + pow (Y/2-y,2) + pow (Z/2-z,2) );

	return true;
}		/* -----  end of method ClassMedium::RandPos  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassQuencher
 *      Method:  ClassQuencher :: dist
 * Description: calculates distance between the center of this qencher and quencher b 
 *--------------------------------------------------------------------------------------
 */

	double
ClassQuencher::distance ( ClassQuencher&  b, double X, double Y, double Z)
{
	double x2, y2, z2;      /* coordinates of b
	       			when this quencher would be shifted to the center*/

	/* shift quencher accordingly  */
	x2 =b.x + X/2 - x;
	y2 =b.y + Y/2 - y;
	z2 =b.z + Z/2 - z;


	if (x2>X) x2-=X;     /* periodic boundary conditions */
	else if (x2<0) x2+=X;
	if (y2>Y) y2-=Y;
	else if (y2<0) y2+=Y;
	if (z2>Z) z2-=Z;
	else if (z2<0) z2+=Z;

	return sqrt ( pow(X/2-x2,2) + pow(Y/2-y2,2) + pow(Z/2-z2,2) ) ;
}		/* -----  end of method ClassQuencher::dist  ----- */



/*
 *--------------------------------------------------------------------------------------
 *       Class:  ClassQuencher
 *      Method:  ClassQuencher :: SetCoordinates
 * Description:  sets new coordinates for a quencher; 
 *		  checks for boundary conditions;
 *		  recaclulates distance to the center dDistance
 *--------------------------------------------------------------------------------------
 */

	void
ClassQuencher::SetCoordinates ( double x2, double y2, double z2, double X, double Y, double Z )
{

	if (x2>X) x2-=X;     /* periodic boundary conditions */
	else if (x2<0) x2+=X;
	if (y2>Y) y2-=Y;
	else if (y2<0) y2+=Y;
	if (z2>Z) z2-=Z;
	else if (z2<0) z2+=Z;
	
	x = x2;
	y = y2;
	z = z2;
	
	dCenter= sqrt (  pow(X/2-x,2) + pow (Y/2-y,2) + pow (Z/2-z,2) );

	return ;
}		/* -----  end of method ClassQuencher::SetCoordinates  ----- */

