#include "NeighborOrdering.h"

namespace Tessellation {
	//-------- Function to sort the list of neighbors (using the global list of neighbors- "pos_neighbor") according to its angle with the target particle \n //
	// and an arbitrary direction (vector between target particle and the first neighbor) ...................//
	vector<int> neighborordering(int i, vector<vector<int>>pos_neighbor, vector<Vector3d>pos) {
		if (pos_neighbor[i].size() == 1) {
			vector<int>finall;
			finall.push_back(pos_neighbor[i][0]);
			finall.push_back(pos_neighbor[i][0]);
			return finall;
		}
		vector<int>finallist;
		int i_size = pos_neighbor[i].size();
		map<Real, Real>p_angles;					// Container for particles on positive side of connection line between i and reference
		map<Real, Real>n_angles;					// Container for particles on negative side of connection line between i and reference
		//Building vector from i-th node to 1st neighbour node
		int t_ind = pos_neighbor[i][0];
		Vector3d ot(pos.at(t_ind)(0) - pos.at(i)(0), pos.at(t_ind)(1) - pos.at(i)(1), pos.at(t_ind)(2) - pos.at(i)(2));
		Real testminimum = 999;						// container for closest angle to vector 'ot'
		int test_n = 1;								// safeguard if minimum angle not found
		int change = 0;								// safeguard if minimum angle not found
		Real testangle;
		Real nochangetestangle;
		//find angle with the initial vector with the rest of the nodes and find the minimum angle constructing vector on the positive side
		for (int j = 1; j < i_size; j++) {
			int n_ind = pos_neighbor[i][j];
			Vector3d on(pos.at(n_ind)(0) - pos.at(i)(0), pos.at(n_ind)(1) - pos.at(i)(1), pos.at(n_ind)(2) - pos.at(i)(2));
			pair<Real, Real>resultangles;
			resultangles = angle(ot, on);
			//------------Angle between two vectors - two outputs - first-> smallest angle(regardless of direction) ;  second-> positive angle(0,PI)-----------//
			if (j == 1)nochangetestangle = resultangles.second;
			if (resultangles.first > 0) {
				if (resultangles.first < testminimum) {
					testminimum = resultangles.first;
					testangle = resultangles.second;
					test_n = n_ind;
					change++;
				}
			}
		}
		// if no angle found on the positive side, then put the first neighbor as reference
		if (change == 0) {
			test_n = pos_neighbor[i][1];
			p_angles[nochangetestangle] = test_n;
		}
		Vector3d on(pos.at(test_n)(0) - pos.at(i)(0), pos.at(test_n)(1) - pos.at(i)(1), pos.at(test_n)(2) - pos.at(i)(2));				// reference neighbor with which other neighbors will be tested to find which side they are on
		// creating a cross vector to test sides 
		Vector3d t_cross = ot.cross(on);
		if (change != 0) {
			p_angles[testangle] = test_n;
		}
		// iterating through neighbors to determine their side with respect to the reference neighbor
		for (int j = 1; j < i_size; j++) {
			int n_ind = pos_neighbor[i][j];
			if (n_ind != test_n) {
				Vector3d on_f(pos.at(n_ind)(0) - pos.at(i)(0), pos.at(n_ind)(1) - pos.at(i)(1), pos.at(n_ind)(2) - pos.at(i)(2));
				pair<Real, Real>resultangles;
				resultangles = angle(ot, on_f);
				Vector3d n_cross = ot.cross(on_f);
				if (t_cross.dot(n_cross) >= 0) {
					if (isnan(resultangles.second)) {

						p_angles[180] = n_ind;
					}
					else {
						p_angles[resultangles.second] = n_ind;

					}
				}
				else {
					n_angles[resultangles.second] = n_ind;
				}
			}
		}
		vector<Real>vector_p;
		vector<Real>vector_n;
		vector_p = sort(p_angles);					// Sorted list of particles on positive side of connection line between i and reference
		vector_n = sort(n_angles);					// Sorted list of particles on negative side of connection line between i and reference
		// Building up of the final ordered list
		if (vector_p.empty()) {
			finallist.push_back(t_ind);
			int nsize = vector_n.size();
			for (int t = 0; t < nsize; t++) {
				finallist.push_back(vector_n.at(t));
			}
			finallist.push_back(t_ind);
		}
		else if (vector_n.empty()) {
			finallist.push_back(t_ind);
			int psize = vector_p.size();
			for (int t = 0; t < psize; t++) {
				finallist.push_back(vector_p.at(t));
			}
			finallist.push_back(t_ind);
		}
		else {
			int nsize = vector_n.size();
			int psize = vector_p.size();
			for (int t = psize - 1; t >= 0; t--) {
				finallist.push_back(vector_p.at(t));
			}
			finallist.push_back(t_ind);
			for (int t = 0; t < nsize; t++) {
				finallist.push_back(vector_n.at(t));
			}
			finallist.push_back(vector_p.at(psize - 1));
		}
		return finallist;
	}


	//-------- Same function as above but the list of neigbhors is provided through "givencontainer" instead of using the global list .........//
	vector<int> neighborordering(int i, vector<int>givencontainer, vector<vector<int>>pos_neighbor, vector<Vector3d>pos) {
		vector<int>finallist;
		int currentnode = i;
		Real isize = givencontainer.size();
		map<Real, Real>p_angles;
		map<Real, Real>n_angles;
		int t_ind = givencontainer[0];
		Vector3d ot(pos.at(t_ind)(0) - pos.at(currentnode)(0), pos.at(t_ind)(1) - pos.at(currentnode)(1), pos.at(t_ind)(2) - pos.at(currentnode)(2));
		Real testminimum = 999;						// container for closest angle to vector 'ot'
		int test_n = 1;								// safeguard if minimum angle not found
		int change = 0;								// safeguard if minimum angle not found
		Real testangle;
		Real nochangetestangle;
		//find angle with the initial vector with the rest of the nodes and find the minimum angle constructing vector on the positive side
		for (int jj = 1; jj < isize; jj++) {
			int n_ind = givencontainer[jj];
			Vector3d on(pos.at(n_ind)(0) - pos.at(currentnode)(0), pos.at(n_ind)(1) - pos.at(currentnode)(1), pos.at(n_ind)(2) - pos.at(currentnode)(2));
			pair<Real, Real>resultangles;
			resultangles = angle(ot, on);
			if (jj == 1)nochangetestangle = resultangles.second;
			//------------Angle between two vectors - two outputs - first-> smallest angle(regardless of direction) ;  second-> positive angle(0,PI)-----------//
			if (resultangles.first > 0) {
				if (resultangles.first < testminimum) {
					testminimum = resultangles.first;
					testangle = resultangles.second;
					test_n = n_ind;
					change++;
				}
			}
		}
		// if no angle found on the positive side, then put the first neighbor as reference
		if (change == 0) {
			test_n = givencontainer[1];
			p_angles[nochangetestangle] = test_n;
		}
		Vector3d on(pos.at(test_n)(0) - pos.at(currentnode)(0), pos.at(test_n)(1) - pos.at(currentnode)(1), pos.at(test_n)(2) - pos.at(currentnode)(2));
		Vector3d t_cross = ot.cross(on);
		if (change != 0) {
			p_angles[testangle] = test_n;
		}
		// iterating through neighbors to determine their side with respect to the reference neighbor
		for (int jj = 1; jj < isize; jj++) {
			int n_ind = givencontainer[jj];
			if (n_ind != test_n) {
				Vector3d on_f(pos.at(n_ind)(0) - pos.at(currentnode)(0), pos.at(n_ind)(1) - pos.at(currentnode)(1), pos.at(n_ind)(2) - pos.at(currentnode)(2));
				pair<Real, Real>resultangles;
				resultangles = angle(ot, on_f);
				Vector3d n_cross = ot.cross(on_f);
				if (t_cross.dot(n_cross) >= 0) {
					if (isnan(resultangles.second)) {

						p_angles[180] = n_ind;
					}
					else {
						p_angles[resultangles.second] = n_ind;

					}
				}
				else {
					n_angles[resultangles.second] = n_ind;
				}
			}
		}
		vector<Real>vector_p;
		vector<Real>vector_n;
		vector_p = sort(p_angles);					// Sorted list of particles on positive side of connection line between i and reference
		vector_n = sort(n_angles);					// Sorted list of particles on negative side of connection line between i and reference
		// Building up of the final ordered list
		if (vector_p.empty()) {
			finallist.push_back(t_ind);
			int nsize = vector_n.size();
			for (int t = 0; t < nsize; t++) {
				finallist.push_back(vector_n.at(t));
			}
			finallist.push_back(t_ind);
		}
		else if (vector_n.empty()) {
			finallist.push_back(t_ind);
			int psize = vector_p.size();
			for (int t = 0; t < psize; t++) {
				finallist.push_back(vector_p.at(t));
			}
			finallist.push_back(t_ind);
		}
		else {
			int nsize = vector_n.size();
			int psize = vector_p.size();
			for (int t = psize - 1; t >= 0; t--) {
				finallist.push_back(vector_p.at(t));
			}
			finallist.push_back(t_ind);
			for (int t = 0; t < nsize; t++) {
				finallist.push_back(vector_n.at(t));
			}
			finallist.push_back(vector_p.at(psize - 1));
		}
		return finallist;
	}

}