/*
 * PreProcessorOptions.h
 *
 *  Created on: Sep 21, 2016
 *      Author: almium
 */

#ifndef PREPROCESSOROPTIONS_H_
#define PREPROCESSOROPTIONS_H_

// set NO_Z_PERIODIC to 1 if periodic conditions along z-axis are not needed.
// This is useful for interface quenching. In this way the two interfaces of the box, which are normal
// to the z-axis will have exciton reflecting conditions at z=0, and z=film_thickness
// If you need periodic conditions along z-axis, set NO_Z_PERIODIC to 0

#define NO_Z_PERIODIC 0



#endif /* PREPROCESSOROPTIONS_H_ */
