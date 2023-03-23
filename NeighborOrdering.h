#ifndef _NEIGHBORORDERING_H_
#define _NEIGHBORORDERING_H_
#include "HelperFunctions.h"
#include "Tessellation.h"


namespace Tessellation {
	//-------- Function to sort the list of neighbors (using the global list of neighbors- "pos_neighbor") according to its angle with the target particle \n //
	// and an arbitrary direction (vector between target particle and the first neighbor) ...................//
	vector<int> neighborordering(int i, vector<vector<int>>pos_neighbor, vector<Vector3d>pos);

	//-------- Same function as above but the list of neigbhors is provided through "givencontainer" instead of using the global list .........//	
	vector<int> neighborordering(int i, vector<int>givencontainer, vector<vector<int>>pos_neighbor, vector<Vector3d>pos);
};

#endif // _NEIGHBORORDERING_H_
