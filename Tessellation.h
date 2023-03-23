#ifndef _Tessellation_H_
#define _Tessellation_H_

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>
#include <string>

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<algorithm>
#include<array>
#include<vector>
#include<Eigen>
#include<math.h>
#include<stdlib.h>
#include<queue>
#include<map>
#include<utility>
#include<time.h>
#include<ctime>
#include <iomanip>



using namespace std;
using namespace Eigen;


namespace Tessellation {

	// Change double to float for single precision
	typedef double Real;
	//Real PI = 3.14159265;


	

	// Equilibrium distance R_e (If not acquired from an MC simulation, predefine some value here)

	//1002//const Real R_e = 0.5;
	//396//const Real R_e = 0.16;
	//962//const Real R_e = 26.0;
	//920//const Real R_e = 0.31;
	//300//const Real R_e = 0.8;
	//162//const Real R_e = 0.5;
	//49//const Real R_e = 2.0;
	//36//const Real R_e = 0.3;
	//MC//
	const Real R_e = 0.2;
	const Real R_c = 2.23 * R_e;        // Cut-off distance for particle search
	const Real delr = 0.25 * R_e;	    // Window width
	const Real dr = delr / 40;		    // Radius change
	const Real initialr = delr / 2; 	// Starting radius

	extern vector<vector<int>>pos_neighbor;		// neighbour nodes of target particles
	extern vector<Vector3d>pos;					// position of particles
}; // namespace Tessellation

#endif // _Tessellation_H_

