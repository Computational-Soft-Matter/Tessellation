#ifndef _NEIGHBORORDERING_H_
#define _NEIGHBORORDERING_H_
#include "HelperFunctions.h"
#include "Tessellation.h"


namespace Tessellation {
	//-------- Function to sort the list of neighbors (using the global list of neighbors- "pos_neighbor") according to their angle with the target particle \n //
	// and an arbitrary direction (vector between target particle and the first neighbor) (Algorithm 4) ...................//
	vector<int> neighborordering(int i, vector<vector<int>>&pos_neighbor, vector<Vector3d>&pos);

	//-------- Same function as above but the list of neigbhors is provided through "givencontainer" instead of using the global list \n //	
	//-------- This allows the use of a modified list instead of the initial input -------- //
	vector<int> neighborordering(int i, vector<int>&givencontainer, vector<vector<int>>&pos_neighbor, vector<Vector3d>&pos);
	
	//-------- Same function as above but the reference vector direction is provided through t_ind -------- //
	vector<int> neighborordering(int i, vector<vector<int>>& pos_neighbor, vector<Vector3d>& pos, int t_ind);

	//............... Function to update the deletable connections list (See Sec. 2.1.2).............//
	Real maxdiffneighbor(int i, int node_j,int del, vector<vector<int>>& pos_neighbor, vector<Vector3d>& pos);
};

#endif // _NEIGHBORORDERING_H_
