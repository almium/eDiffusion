/*
 * =====================================================================================
 *
 *       Filename:  ClassBool3D.cpp
 *
 *    Description:  Contains methods of ClassBool3D
 *
 *        Version:  1.0
 *        Created:  10/12/2011 10:00:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Oleksandr (Alex) Mikhnenko (www.mikhnenko.com), alex@mikhnennko.com
 *        Company:  University of Groningen
 *
 * =====================================================================================
 */

#include	<stdlib.h>
#include	<math.h>
#include	"eDiffusion.h"



/*
 *--------------------------------------------------------------------------------------
 *       Class:  Bool3D
 *      Method:  Bool3D :: set
 * Description:  Sets of a specific cell on the grid;
 * 		 accounts periodic boundary conditions
 *--------------------------------------------------------------------------------------
 */

void
Bool3D::set (int x,int y, int z, bool value)
{
	int nx, ny, nz; 
	nx=x;
	ny=y;
	nz=z;

	if (nx>=Nbx) nx-=Nbx;     /* periodic boundary conditions */
	else if (nx<0) nx+=Nbx;
	if (ny>=Nby) ny-=Nby;
	else if (ny<0) ny+=Nby;
	if (nz>=Nbz) nz-=Nbz;
	else if (nz<0) nz+=Nbz;


	if (vec[nz*Nbx*Nby + ny*Nbx + nx]==true)
	{
		//cout << "reasigning the spot" << " "<< "(nx, ny, nz)=(" << nx <<"," << ny << "," << nz <<")\n" ; 
	}

	vec[nz*Nbx*Nby + ny*Nbx + nx]=value;
	return ;
}		/* -----  end of method Bool3D::set  ----- */



