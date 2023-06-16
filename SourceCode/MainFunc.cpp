#include "Tessellation.h"
#include "HelperFunctions.h"
#include "NeighborOrdering.h"
using namespace Tessellation;
//---------Struct to contain information about hole formation-----------//

struct store
{
	int node;
	int neigh;
	int totfix;
	Real avgdist;
};

//---------Struct to contain information about deleting process-----------//

struct storedelete
{
	Real dist;
	int deletefrom;
	int deletenode;
};

//---------Struct to contain information about average point deduction-----------//

struct founds
{
	vector<int>alreadyfound;
	Vector3d Point;
};

//----------Compare function to sort based on distance from struct-"storedelete".........//

bool CompareData(store a, store b)
{
	return a.avgdist < b.avgdist;
}

bool CompareDatadlt(storedelete a, storedelete b)
{
	return a.dist > b.dist;
}
//--------------Starting Main Function-------------------//

int main() {
	string input_file = "und5exp";								// File contains the 3D coordinates of the colloidal particles
	string input_extension = ".vtk";
	string neighbor_input = "und5expdata.txt";								// If neighbor list has been calculated previously, then instead of NULL, the file will contain the index of a particle, number of neighbors and the radius of search to acquire those neighbors
	string output_capsomer = input_file + "caps.vtk";			// Output file name that has the information of the capsomer cells
	string output_outline = input_file + "outline.vtk";			// Output file name that has the information of capsomer outline
	clock_t starttime = clock();
	int totalpoints = 0;					// Total number of points in the final VTK file containing all the capsomers (this includes the initial set of points and the vertices of the capsomers)
	founds Emp[1000];						// Data structure for averaging points in empty region - (container to store the vertices of quadrilaterals that uses average capsomer vertices position described in section 2.1.4
	int empcount = 0;						// Data structure for averaging points in empty region
	int totalcaps = 0;						// total number of capsomers in the final vtk file
	vector<Real>n_neighbor;		            // neighbour information of target particles - first -> number of n_neighbor, second-> radius range of search
	vector<int>vtknodes;					// particle IDs for final VTK file (including capsomer vertices)
	vector<Real>avg_neighbor;				// average distances for n_neighbor
	vector<pair<int, int>>single_edges;		// pairs of single edge connections
	map<int, int>marked;                    // marked node which has single edge connections and changed capsomers
	vector<vector<int>>pos_neighbor;		// neighbour nodes of target particles
	vector<Vector3d>pos;				// position of particles
	//n_part : number of particles -> can be changed according to the system

	const int n_part = 1134;			// Total Number of Particles in the inputfile
	const int nonrepeat_part = 1000;

	//---------------------Reading particle positions-------------------//
	string myText;
	ifstream MyReadFile(input_file + input_extension);
	while (getline(MyReadFile, myText)) {
		string word;
		istringstream ss(myText);
		Real temp_pos[3] = {};
		int gp = 0;
		while (ss >> word)
		{
			temp_pos[gp] = stof(word);
			gp++;
			if (gp == 3) {
				pos.push_back({ temp_pos[0], temp_pos[1], temp_pos[2] });
				gp = 0;
			}
		}
	}
	MyReadFile.close();

	//-----------------Reading neighbour information if given-------------//
	if (neighbor_input != "NULL") {
		string myText1;
		ifstream MyReadFile1(neighbor_input);
		while (getline(MyReadFile1, myText1)) {
			string word;
			int count = 0;
			istringstream ss(myText1);
			Real temp_pos[3] = {};
			int gp = 0;
			while (ss >> word)
			{
				temp_pos[gp] = stof(word);
				gp++;
				if (gp == 3) {
					n_neighbor.push_back(temp_pos[2]);
					gp = 0;
					count++;
				}
			}
		}
		MyReadFile1.close();
	}
	//ofstream m;									
	//m.open("PositionsFile2.txt");						//............... File to keep the record of initial neighboring particles' container for debugging purposes....//
	//-------------Neighbour search algorithm using moving window-----------------//
	if (neighbor_input == "NULL") {
		//m.open(input_file + "pos.txt");
		for (int i = 0; i < n_part; i++) {
			Real isstop = 0;
			Real init_r = initialr;
			vector<int>list_p;
			vector<Real>list_n;
			Real it = 0;
			Real isconst = 0;
			map<Real, Real>constcount;
			map<Real, Real>trackcount;
			Real maxn = -10;
			Real maxr;
			Real lastcount_n = 0;
			while (isstop != 2) {
				Real count_n = 0;
				for (int j = 0; j < n_part; j++) {
					Real r = dist(pos.at(i), pos.at(j));
					if (r <= init_r + delr / 2 && r >= init_r - delr / 2 && i != j) {
						count_n++;
					}
				}
				Real countdist = count_n - maxn;
				lastcount_n = count_n;
				if (maxn < count_n) {
					maxn = count_n;
					maxr = init_r;
				}
				list_p.push_back(count_n);
				list_n.push_back(init_r);
				if (init_r > R_c || count_n > 15 || countdist < -2) {
					isstop = 2;
				}
				it++;
				init_r += dr;
			}
			n_neighbor.push_back(maxr);
			//m << i << " " << maxn << " " << maxr << endl;
			//cout << i << " " << maxn << " " << maxr << endl;				//..................Printing the radius of search for each particle and the number of particles found as initial neighbors for debugging...............// 
		}
		//m.close();
	}
	//-------------Build-up of particle neighbour containers------------//
	for (int i = 0; i < n_part; i++) {
		pos_neighbor.push_back(vector<int>());
	}
	for (int i = 0; i < n_part; i++) {
		Real range = n_neighbor.at(i);
		for (int j = 0; j < n_part; j++) {
			Real r = dist(pos.at(i), pos.at(j));
			if (r <= range + delr / 2 && i != j) {
				int isize = pos_neighbor[i].size();
				int jsize = pos_neighbor[j].size();

				// If 'i' contains 'j' then the vice versa should be true, made sure in the next part of the code

				int ci = 0, cj = 0;
				for (int k = 0; k < isize; k++) {
					if (pos_neighbor[i].at(k) == j) {
						ci++;
					}
				}
				for (int k = 0; k < jsize; k++) {
					if (pos_neighbor[j].at(k) == i) {
						cj++;
					}
				}
				if (ci == 0) {
					pos_neighbor[i].push_back(j);
				}
				if (cj == 0) {
					pos_neighbor[j].push_back(i);
				}
			}
		}
	}

	//-------------------Calculating average neighbour distances-------------//
	for (int i = 0; i < n_part; i++) {
		int nsize = pos_neighbor[i].size();
		Real total = 0;
		for (int j = 0; j < nsize; j++) {
			Real r = dist(pos.at(i), pos.at(pos_neighbor[i].at(j)));
			total += r;
		}
		avg_neighbor.push_back(total / nsize);
	}


	//-----------------Deleting overlapping particle edge-connections--------------//
	for (int i = 0; i < n_part; i++) {
		int deleted = 0;											// Number of particles already deleted from neighbor container
		int i_size = pos_neighbor[i].size();						// Size of neighbor container for particle 'i'
		for (int j = 0; j < i_size - deleted; j++) {
			int node_j = pos_neighbor[i].at(j);						// Iterating through particle i's neighborhood and establishing multiple shared particle connections
			vector<int>deletelist;
			for (int k = 0; k < i_size - deleted; k++) {
				int nodenow = pos_neighbor[i].at(k);
				if (nodenow != i || nodenow != node_j) {
					if (vector_search(pos_neighbor[i], nodenow) == 1 && vector_search(pos_neighbor[node_j], nodenow) == 1) {
						deletelist.push_back(nodenow);				// For particle j in particle i's neighborhood, store particles that are commonly shared by i and j
					}
				}
			}
			map<Real, Real>p_angles;								// Container for particles on positive side of connection line between i and j
			map<Real, Real>n_angles;								// Container for particles on negative side of connection line between i and j
			//Building vector from i-th node to 1st neighbour node
			int t_ind = node_j;
			Vector3d ot(pos.at(t_ind)(0) - pos.at(i)(0), pos.at(t_ind)(1) - pos.at(i)(1), pos.at(t_ind)(2) - pos.at(i)(2));
			Real testminimum = 999;									// container for closest angle to vector 'ot'
			int test_n = 1;											// safeguard if minimum angle not found
			int change = 0;											// Check value to see if minimum angle found		
			Real testangle;
			Real testangletest;
			//find angle with the initial vector with the rest of the nodes and find the minimum angle constructing vector
			for (int j = 0; j < deletelist.size(); j++) {
				int n_ind = deletelist[j];
				Vector3d on(pos.at(n_ind)(0) - pos.at(i)(0), pos.at(n_ind)(1) - pos.at(i)(1), pos.at(n_ind)(2) - pos.at(i)(2));
				pair<Real, Real>resultangles;
				resultangles = angle(ot, on);     //----Angle between two vectors - two outputs - first-> smallest angle(regardless of direction) ;  second-> positive angle(0,PI)----//
				if (j == 1)testangletest = resultangles.second;
				if (resultangles.first > 0) {
					if (resultangles.first < testminimum) {
						testminimum = resultangles.first;
						testangle = resultangles.second;
						test_n = n_ind;
						change++;
					}
				}
			}
			//Building vector from i-th node to 1st neighbour node
			Vector3d on(pos.at(test_n)(0) - pos.at(i)(0), pos.at(test_n)(1) - pos.at(i)(1), pos.at(test_n)(2) - pos.at(i)(2));
			Vector3d t_cross = ot.cross(on);
			if (change != 0) {
				p_angles[testangle] = test_n;
			}
			// Finding on which side each particle belong to
			for (int j = 0; j < deletelist.size(); j++) {
				int n_ind = deletelist[j];
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
			vector_p = sort(p_angles);							// sorted list of particles on the positive side of i and j according to angle
			vector_n = sort(n_angles);							// sorted list of particles on the negative side of i and j according to angle
			//----------------Deleting overlapping particles from positive side-----------------//
			int psize = vector_p.size();
			storedelete P[10];
			int co = 0;
			for (int ite = 0; ite < vector_p.size(); ite++) {
				int pnow = vector_p[ite];
				if (dist(pos.at(i), pos.at(pnow)) - avg_neighbor[i] < dist(pos.at(node_j), pos.at(pnow)) - avg_neighbor[node_j]) {
					P[co].dist = (dist(pos.at(node_j), pos.at(pnow)) - avg_neighbor[node_j]) + (dist(pos.at(i), pos.at(pnow)) - avg_neighbor[i]);
					P[co].deletefrom = node_j;
					P[co].deletenode = pnow;
					co++;
				}
				else {
					P[co].dist = (dist(pos.at(i), pos.at(pnow)) - avg_neighbor[i]) + (dist(pos.at(node_j), pos.at(pnow)) - avg_neighbor[node_j]);
					P[co].deletefrom = i;
					P[co].deletenode = pnow;
					co++;
				}
			}
			for (int ite = co; ite < 10; ite++) {
				P[ite].dist = -999;
				P[ite].deletefrom = -999;
				P[ite].deletenode = -999;
			}
			sort(P, P + 10, CompareDatadlt);
			for (int ite = 0; ite < psize - 1; ite++) {
				int deletnow = P[ite].deletenode;
				int dltfr = P[ite].deletefrom;
				int deletat;
				for (int y = 0; y < pos_neighbor[dltfr].size(); y++) {
					if (deletnow == pos_neighbor[dltfr].at(y)) {
						deletat = y;
					}
				}
				deleted++;
				pos_neighbor[dltfr].erase(pos_neighbor[dltfr].begin() + deletat);
				int dltfr2 = P[ite].deletenode;
				int deletnow2 = P[ite].deletefrom;
				int deletat2;
				for (int y = 0; y < pos_neighbor[dltfr2].size(); y++) {
					if (deletnow2 == pos_neighbor[dltfr2].at(y)) {
						deletat2 = y;
					}
				}
				pos_neighbor[dltfr2].erase(pos_neighbor[dltfr2].begin() + deletat2);
			}
			//----------------Deleting overlapping particles from negative side-----------------//
			int nsize = vector_n.size();
			storedelete N[10];
			int cor = 0;
			for (int ite = 0; ite < vector_n.size(); ite++) {
				int nnow = vector_n[ite];
				if ((dist(pos.at(i), pos.at(nnow)) - avg_neighbor[i]) < (dist(pos.at(node_j), pos.at(nnow)) - avg_neighbor[node_j])) {
					N[cor].dist = (dist(pos.at(node_j), pos.at(nnow)) - avg_neighbor[node_j]) + (dist(pos.at(i), pos.at(nnow)) - avg_neighbor[i]);
					N[cor].deletefrom = node_j;
					N[cor].deletenode = nnow;
					cor++;
				}
				else {
					N[cor].dist = (dist(pos.at(i), pos.at(nnow)) - avg_neighbor[i]) + (dist(pos.at(node_j), pos.at(nnow)) - avg_neighbor[node_j]);
					N[cor].deletefrom = i;
					N[cor].deletenode = nnow;
					cor++;
				}
			}
			for (int ite = cor; ite < 10; ite++) {
				N[ite].dist = -999;
				N[ite].deletefrom = -999;
				N[ite].deletenode = -999;
			}
			sort(N, N + 10, CompareDatadlt);
			for (int ite = 0; ite < nsize - 1; ite++) {
				int deletnow = N[ite].deletenode;
				int dltfr = N[ite].deletefrom;
				int deletat;
				for (int y = 0; y < pos_neighbor[dltfr].size(); y++) {
					if (deletnow == pos_neighbor[dltfr].at(y)) {
						deletat = y;
					}
				}
				deleted++;
				pos_neighbor[dltfr].erase(pos_neighbor[dltfr].begin() + deletat);
				int dltfr1 = N[ite].deletenode;
				int deletnow1 = N[ite].deletefrom;
				int deletat1;
				for (int y = 0; y < pos_neighbor[dltfr1].size(); y++) {
					if (deletnow1 == pos_neighbor[dltfr1].at(y)) {
						deletat1 = y;
					}
				}
				pos_neighbor[dltfr1].erase(pos_neighbor[dltfr1].begin() + deletat1);
			}
		}
	}

	// For each particle going through its ordered neighbors and deleting overlaps due to distantly placed particles (if two adjacent neighbors of a neighbor are connected, then that is deleted)
	for (int i = 0; i < n_part; i++) {
		if (pos_neighbor[i].size() < 3)continue;
		vector<int>finallist;
		finallist = neighborordering(i, pos_neighbor, pos);
		finallist.pop_back();
		if (finallist.size() == 3) {
			int nextcheck = finallist[2];
			int pastcheck = finallist[0];
			if (vector_search(pos_neighbor[nextcheck], pastcheck) == 1) {
				if (dist(pos.at(finallist.at(1)), pos.at(i)) < dist(pos.at(nextcheck), pos.at(i)) && dist(pos.at(finallist.at(1)), pos.at(i)) < dist(pos.at(pastcheck), pos.at(i)))continue;
				int deletenow = finallist[1];
				int deletfrom = i;
				if (vector_search(pos_neighbor[deletfrom], deletenow) == 0)continue;
				int deletat = find_ind(pos_neighbor[deletfrom], deletenow);
				int deletat2 = find_ind(pos_neighbor[deletenow], deletfrom);
			}
		}
		else {
			for (int j = 0; j < finallist.size(); j++) {
				int nextcheck = finallist[(j + 1) % finallist.size()];
				int pastcheck = finallist[((j - 1) + finallist.size()) % finallist.size()];
				if (vector_search(pos_neighbor[pastcheck], nextcheck) == 1) {
					int deletenow = finallist[j];
					int deletfrom = i;
					if (vector_search(pos_neighbor[deletfrom], deletenow) == 0)continue;
					int deletat = find_ind(pos_neighbor[deletfrom], deletenow);
					pos_neighbor[deletfrom].erase(pos_neighbor[deletfrom].begin() + deletat);
					int deletat2 = find_ind(pos_neighbor[deletenow], deletfrom);
					pos_neighbor[deletenow].erase(pos_neighbor[deletenow].begin() + deletat2);
				}
			}
		}
	}

	//--------------------Finding single edge-connections for basis of hole formation--------------//
	vector<vector<int>>single_conn;								// Connectivity of single edges
	vector<pair<int, int>>single_edge;							// Container for pairs of single edges
	map<int, bool>visited;										// Check if node already visited
	vector<int>listofsingles;									// list of single edged nodes which needs to be repaired
	for (int i = 0; i < n_part; i++) {
		single_conn.push_back(vector<int>());
	}
	for (int i = 0; i < n_part; i++) {
		int deleted = 0;
		int i_size = pos_neighbor[i].size();
		for (int j = 0; j < i_size + deleted; j++) {
			int eiba = 0;
			int noden = pos_neighbor[i].at(j);
			int commoncount = 0;
			for (int k = 0; k < i_size + deleted; k++) {
				int nodenow = pos_neighbor[i].at(k);
				if (nodenow != i && nodenow != noden) {
					if (vector_search(pos_neighbor[i], nodenow) == 1 && vector_search(pos_neighbor[noden], nodenow) == 1) {
						commoncount++;
					}
				}
			}
			int baddibo = 0;
			if (commoncount == 1 && i_size >= 2 && pos_neighbor[noden].size() >= 2) {
				if (visited[i] != 1) {
					listofsingles.push_back(i);
					visited[i] = 1;
				}
				if (visited[noden] != 1) {
					listofsingles.push_back(noden);
					visited[noden] = 1;
				}
				single_edge.push_back(make_pair(i, noden));
				single_conn.at(i).push_back(noden);
			}
		}
	}

	//---------------------Hole connectivity formation-----------------//
	map<int, bool>worked;									// container to indicate how many times a single edged node has been worked with
	vector<vector<int>>holes;								// container for all the chain of holes
	for (int i = 0; i < listofsingles.size(); i++) {
		int workingnow = listofsingles.at(i);				// single edged node that is being used to create chain of hole 
		if (single_conn.at(workingnow).size() >= 2) {
			if ((worked[workingnow] == 0 && single_conn.at(workingnow).size() >= 2) || (worked[workingnow] == 1 && single_conn.at(workingnow).size() > 2)) {
				vector<int>hole;							// container for the chain of hole
				hole.push_back(workingnow);
				// Finding the next node in chain of hole depending on past hole constructions
				int next = -5;
				int past;
				if (worked[workingnow] == 0) {
					next = single_conn.at(workingnow).at(0);
					past = workingnow;
					if (worked[next] == 1) {
						if (single_conn.at(next).size() > 2) {
							worked[next] == 2;
						}
					}
					else (worked[next] == 1);
				}
				else {
					for (int ll = 0; ll < single_conn.at(workingnow).size(); ll++) {
						if (worked[single_conn.at(workingnow).at(ll)] == 0) {
							next = single_conn.at(workingnow).at(ll);
							worked[next] = 1;
							break;
						}
						else if (worked[single_conn.at(workingnow).at(ll)] == 1 && single_conn.at(single_conn.at(workingnow).at(ll)).size() > 2) {
							next = single_conn.at(workingnow).at(ll);
							worked[next] = 2;
							break;
						}
					}
					if (next == -5)continue;
					past = workingnow;
				}
				hole.push_back(next);

				// Start creating the chain of holes through this while loop
				int stop = 0;					// Stopping crtieria for the while loop
				while (next != workingnow && stop == 0) {
					if (single_conn.at(next).size() < 2) {
						stop = 1;											// chain ends if no more single edged nodes to discover
					}
					else if (single_conn.at(next).size() == 2) {			// this condition is visited if their are one single edged node to discover from a position in the chain of hole
						if (single_conn.at(next).at(0) == past) {
							past = next;
							next = single_conn.at(next).at(1);
							if (worked[next] == 1) {
								if (single_conn.at(next).size() > 2) {
									worked[next] == 2;
								}
								else stop = 1;
							}
							else if (worked[next] == 2) {
								stop = 1;
							}
							else {
								worked[next] = 1;
							}
							hole.push_back(next);
						}
						else {
							past = next;
							next = single_conn.at(next).at(0);
							if (worked[next] == 1) {
								if (single_conn.at(next).size() > 2) {
									worked[next] == 2;
								}
								else stop = 1;
							}
							else if (worked[next] == 2) {
								stop = 1;
							}
							else {
								worked[next] = 1;
							}
							worked[next] = 1;
							hole.push_back(next);
						}
					}
					else {													// this condition is visited if their are multiple single edged nodes to discover from a position in the chain of hole
						int nextdirect = 0;
						vector<int>finallist;
						finallist = neighborordering(next, pos_neighbor, pos);
						int option1;
						int option2;
						int pastpos = find_ind(finallist, past);
						int ns = finallist.size();
						if (pastpos == 0) {
							int previous = finallist.size() - 2;
							int poroborti = 1;
							if (vector_search(pos_neighbor[past], finallist[previous]) == 1) {
								if (vector_search(single_conn.at(next), finallist[poroborti]) == 1)nextdirect = finallist[poroborti];
								else {
									int poroporoborti = (poroborti + 1) % ns;
									nextdirect = finallist[poroporoborti];
								}
							}
							else {
								if (vector_search(single_conn.at(next), finallist[previous]) == 1)nextdirect = finallist[previous];
								else {
									int preprevious = ((previous - 1) % ns + ns) % ns;
									nextdirect = finallist[preprevious];
								}
							}
						}
						else {
							int previous = ((pastpos - 1) % ns + ns) % ns;
							int poroborti = (pastpos + 1) % ns;
							if (vector_search(pos_neighbor[past], finallist[previous]) == 1) {
								if (vector_search(single_conn.at(next), finallist[poroborti]) == 1)nextdirect = finallist[poroborti];
								else {
									int poroporoborti;
									if (poroborti == finallist.size() - 1) poroporoborti = 1;
									else poroporoborti = (poroborti + 1) % ns;
									nextdirect = finallist[poroporoborti];
								}
							}
							else {
								if (vector_search(single_conn.at(next), finallist[previous]) == 1)nextdirect = finallist[previous];
								else {
									int preprevious;
									if (previous == 0) preprevious = finallist.size() - 2;
									else { preprevious = ((previous - 1) % ns + ns) % ns; }
									nextdirect = finallist[preprevious];
								}
							}
						}
						past = next;
						next = nextdirect;
						if (worked[next] == 1) {
							if (single_conn.at(next).size() > 2) {
								worked[next] == 2;
							}
							else stop = 1;
						}
						else if (worked[next] == 2) {
							stop = 1;
						}
						else {
							worked[next] = 1;
						}
						worked[next] = 1;
						hole.push_back(next);
					}
				}
				holes.push_back(hole);
			}
		}
	}


	for (int i = 0; i < n_part; i++) {
		int nsize = pos_neighbor[i].size();
		Real total = 0;
		for (int j = 0; j < nsize; j++) {
			Real r = dist(pos.at(i), pos.at(pos_neighbor[i].at(j)));
			total += r;
		}
		avg_neighbor[i] = (total / nsize);
	}

	//-----------------------Hole Repair-------------------------//
	for (int i = 0; i < holes.size(); i++) {
		int holesize = holes.at(i).size() - 1;				// size of current chain of hole
		store nodes[1000];
		map<int, bool>done_fixing;
		vector<vector<int>>fixes;							// container of list of fixes for each node in hole
		for (int j = 0; j < holes.at(i).size() - 1; j++) {
			int nodenow = holes.at(i).at(j);
			vector<int>possiblefixes;						// container of fixes for nodenow - j-th node in chain of hole
			possiblefixes.push_back(nodenow);				// for j-th node in chain of hole, find which nodes can be patched up with j
			storedelete NN[300];
			int NNco = 0;
			int total = 0;
			Real disttot = 0;
			for (int k = 0; k < holes.at(i).size() - 1; k++) {
				Real r = dist(pos.at(nodenow), pos[holes[i][k]]);
				Real diff = r - avg_neighbor[nodenow];
				vector<int>cont;
				cont = pos_neighbor[nodenow];
				cont.push_back(holes[i][k]);
				vector<int>finalfix;
				finalfix = neighborordering(nodenow, cont, pos_neighbor, pos);
				int opt1;
				int opt2;
				int pastpos = find_ind(finalfix, holes[i][k]);
				int ns = finalfix.size();
				if (pastpos == 0) {
					opt1 = finalfix[1];
					opt2 = finalfix[finalfix.size() - 2];
				}
				else {
					opt1 = finalfix[pastpos + 1];
					opt2 = finalfix[pastpos - 1];
				}
				int overlapholefix = 0;
				if (vector_search(pos_neighbor[opt1], opt2) == 1)overlapholefix = 1;
				if (nodenow != holes[i][k] && vector_search(pos_neighbor.at(nodenow), holes[i][k]) == 0 && r < R_c && overlapholefix == 0) {
					NN[NNco].dist = dist(pos.at(nodenow), pos.at(k));
					NN[NNco].deletefrom = nodenow;
					NN[NNco].deletenode = holes[i][k];
					NNco++;
					disttot += r;
					total++;
				}
			}
			int actualnnco = NNco;
			for (int ite = NNco; ite < 300; ite++) {
				NN[ite].dist = 0;
				NN[ite].deletefrom = 0;
				NN[ite].deletenode = 0;
			}
			sort(NN, NN + 300, CompareDatadlt);
			for (int mm = actualnnco - 1; mm >= 0; mm--) {
				possiblefixes.push_back(NN[mm].deletenode);
			}
			fixes.push_back(possiblefixes);
			nodes[j].node = nodenow;
			nodes[j].neigh = pos_neighbor.at(nodenow).size();
			nodes[j].totfix = total;
			nodes[j].avgdist = disttot / total;
		}
		//sorting the list of nodes to be fixed according to which ones has the closest fixes
		sort(nodes, nodes + holes.at(i).size() - 1, CompareData);

		for (int j = 0; j < holesize; j++) {
			int fixing_now = nodes[j].node;
			done_fixing[fixing_now] = 5;
			int fixing_now_ind = find_ind(holes.at(i), nodes[j].node);
			// Finding which two nodes are adjacent to fixing_now node to check if there is still a single edge connection there
			vector<pair<int, int>>singleedge_fixingnow;
			if (fixing_now_ind != 0) {
				singleedge_fixingnow.push_back(make_pair(holes.at(i).at(fixing_now_ind), holes.at(i).at(fixing_now_ind + 1)));
				singleedge_fixingnow.push_back(make_pair(holes.at(i).at(fixing_now_ind), holes.at(i).at(fixing_now_ind - 1)));
			}
			else
			{
				singleedge_fixingnow.push_back(make_pair(holes.at(i).at(fixing_now_ind), holes.at(i).at(fixing_now_ind + 1)));
				singleedge_fixingnow.push_back(make_pair(holes.at(i).at(fixing_now_ind), holes.at(i).at(holes.at(i).size() - 2)));
			}
			for (int l = 0; l < 2; l++) {
				int currentnode_fix = singleedge_fixingnow.at(l).first;    // the node that is being foxed
				int adjacentnode = singleedge_fixingnow.at(l).second;	   // adjacent node in the hole chain
				int deleted = 0;
				int currentsize = pos_neighbor[currentnode_fix].size();
				int commoncount = 0;
				for (int k = 0; k < currentsize + deleted; k++) {
					int nodenow = pos_neighbor[currentnode_fix].at(k);		// checking if this node is common
					if (nodenow != currentnode_fix && nodenow != adjacentnode) {
						if (vector_search(pos_neighbor[currentnode_fix], nodenow) == 1 && vector_search(pos_neighbor[adjacentnode], nodenow) == 1) {
							commoncount++;
						}
					}
				}
				if (commoncount == 1 && currentsize > 2 && pos_neighbor[adjacentnode].size() > 2) {			// if fixing node is still a single edged connection then keep going for fix	
					storedelete PP[300];
					int PPco = 0;
					int indnow;
					for (int mm = 0; mm < fixes.size(); mm++) {
						if (fixes[mm][0] == currentnode_fix)indnow = mm;
					}
					for (int mm = 1; mm < fixes[indnow].size(); mm++) {
						int currentfix = fixes[indnow][mm];
						int checkvalueforfix = 0;
						if (currentfix != currentnode_fix && vector_search(pos_neighbor.at(currentnode_fix), currentfix) == 0)
						{
							if (dist(pos.at(currentnode_fix), pos.at(currentfix)) < 1.0 * R_c) {
								PP[PPco].dist = dist(pos.at(currentnode_fix), pos.at(currentfix));
								PP[PPco].deletefrom = currentnode_fix;
								PP[PPco].deletenode = currentfix;
								PPco++;
							}
						}
					}
					int actualppco = PPco;
					for (int ite = PPco; ite < 300; ite++) {
						PP[ite].dist = 0;
						PP[ite].deletefrom = 0;
						PP[ite].deletenode = 0;
					}
					sort(PP, PP + 300, CompareDatadlt);
					int confirmfixednode = -5;							// If this is not -5 after checks, then this will contain the node that will be connected to the current node that is being fixed
					for (int mm = actualppco - 1; mm >= 0; mm--) {
						int currentfix = PP[mm].deletenode;
						int checkvalueforfix = 0;
						if (currentfix != currentnode_fix && vector_search(pos_neighbor.at(currentnode_fix), currentfix) == 0)
						{
							if (dist(pos.at(currentnode_fix), pos.at(currentfix)) < 1.0 * R_c) {
								confirmfixednode = currentfix;					// passed 1st check
								// checking problems for currentnode_fix for connecting with confirmfixednode
								vector<int>newcontainer_currentnode_fix;
								newcontainer_currentnode_fix = pos_neighbor[currentnode_fix];
								newcontainer_currentnode_fix.push_back(confirmfixednode);
								vector<int>finallist_currentnode_fix;
								finallist_currentnode_fix = neighborordering(currentnode_fix, newcontainer_currentnode_fix, pos_neighbor, pos);
								int option1;
								int option2;
								int pastpos = find_ind(finallist_currentnode_fix, confirmfixednode);
								int ns = finallist_currentnode_fix.size();
								if (pastpos == 0) {
									option1 = finallist_currentnode_fix[1];
									option2 = finallist_currentnode_fix[finallist_currentnode_fix.size() - 2];
								}
								else {
									option1 = finallist_currentnode_fix[pastpos + 1];
									option2 = finallist_currentnode_fix[pastpos - 1];
								}


								// checking problems for confirmfixednode for connecting with currentnode_fix
								vector<int>newcontainer_confirmfixednode;
								newcontainer_confirmfixednode = pos_neighbor[confirmfixednode];
								newcontainer_confirmfixednode.push_back(currentnode_fix);
								vector<int>finallist_confirmfixednode;
								finallist_confirmfixednode = neighborordering(confirmfixednode, newcontainer_confirmfixednode, pos_neighbor, pos);
								int option3;
								int option4;
								int pastpos1 = find_ind(finallist_confirmfixednode, currentnode_fix);
								int ns1 = finallist_confirmfixednode.size();
								if (pastpos1 == 0) {
									option3 = finallist_confirmfixednode[1];
									option4 = finallist_confirmfixednode[finallist_confirmfixednode.size() - 2];
								}
								else {
									option3 = finallist_confirmfixednode[pastpos1 + 1];
									option4 = finallist_confirmfixednode[pastpos1 - 1];
								}
								//Finding if this fix is creating any problem with any other common nodes of the two nodes being fixed
								vector<int>common;
								for (int iter = 0; iter < pos_neighbor[currentnode_fix].size(); iter++) {
									if (vector_search(pos_neighbor[confirmfixednode], pos_neighbor[currentnode_fix][iter]) == 1) {
										common.push_back(pos_neighbor[currentnode_fix][iter]);
									}
								}
								int stopfix = 0;
								if (common.size() == 0)stopfix = 1;
								for (int il = 0; il < common.size(); il++) {
									vector<int>newcontainer_checkind;
									newcontainer_checkind = pos_neighbor[common[il]];
									vector<int>finallist_checkind;
									int checkind = common[il];
									finallist_checkind = neighborordering(checkind, newcontainer_checkind, pos_neighbor, pos);
									int option5;
									int option6;
									int pastpos2 = find_ind(finallist_checkind, currentnode_fix);
									int ns2 = finallist_checkind.size();
									if (pastpos2 == 0) {
										option5 = finallist_checkind[1];
										option6 = finallist_checkind[finallist_checkind.size() - 2];
									}
									else {
										option5 = finallist_checkind[pastpos2 + 1];
										option6 = finallist_checkind[pastpos2 - 1];
									}
									if (option5 == confirmfixednode || option6 == confirmfixednode) {
										stopfix = 1;
									}
								}
								int finalcheck = 0;
								if (dist(pos.at(currentnode_fix), pos.at(confirmfixednode)) < dist(pos.at(option1), pos.at(currentnode_fix)) && dist(pos.at(currentnode_fix), pos.at(confirmfixednode)) < dist(pos.at(option2), pos.at(currentnode_fix)))finalcheck = 1;
								if (vector_search(pos_neighbor[option1], option2) == 0 && vector_search(pos_neighbor[option3], option4) == 0 && stopfix == 1 || finalcheck == 1) {
									pos_neighbor[currentnode_fix].push_back(confirmfixednode);
									if (vector_search(pos_neighbor[confirmfixednode], currentnode_fix) == 0) {
										pos_neighbor[confirmfixednode].push_back(currentnode_fix);
									}
									deleted++;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i = 0; i < n_part; i++) {
		if (pos_neighbor[i].size() < 3)continue;
		vector<int>finallist;
		finallist = neighborordering(i, pos_neighbor, pos);
		finallist.pop_back();
		if (finallist.size() == 3) {
			int nextcheck = finallist[2];
			int pastcheck = finallist[0];
			if (vector_search(pos_neighbor[nextcheck], pastcheck) == 1) {
				if (dist(pos.at(finallist.at(1)), pos.at(i)) < dist(pos.at(nextcheck), pos.at(i)) && dist(pos.at(finallist.at(1)), pos.at(i)) < dist(pos.at(pastcheck), pos.at(i)))continue;
				int deletenow = finallist[1];
				int deletfrom = i;
				if (vector_search(pos_neighbor[deletfrom], deletenow) == 0)continue;
				int deletat = find_ind(pos_neighbor[deletfrom], deletenow);
				int deletat2 = find_ind(pos_neighbor[deletenow], deletfrom);
			}
		}
		else {
			for (int j = 0; j < finallist.size(); j++) {
				int nextcheck = finallist[(j + 1) % finallist.size()];
				int pastcheck = finallist[((j - 1) + finallist.size()) % finallist.size()];
				if (vector_search(pos_neighbor[pastcheck], nextcheck) == 1) {
					int deletenow = finallist[j];
					int deletfrom = i;
					if (vector_search(pos_neighbor[deletfrom], deletenow) == 0)continue;
					int deletat = find_ind(pos_neighbor[deletfrom], deletenow);
					pos_neighbor[deletfrom].erase(pos_neighbor[deletfrom].begin() + deletat);
					int deletat2 = find_ind(pos_neighbor[deletenow], deletfrom);
					pos_neighbor[deletenow].erase(pos_neighbor[deletenow].begin() + deletat2);
				}
			}
		}
	}

	//--------------------Finding empty connections--------------------------//
	vector<pair<int, int>>empty_connection;								// pairs with no common neighbor between them
	for (int i = 0; i < n_part; i++) {
		int isize = pos_neighbor[i].size();
		for (int j = 0; j < isize; j++) {
			int testnode = pos_neighbor[i].at(j);
			int commoncount = 0;
			for (int k = 0; k < isize; k++) {
				int nodenow = pos_neighbor[i].at(k);
				if (nodenow != i && nodenow != testnode) {
					if (vector_search(pos_neighbor[i], nodenow) == 1 && vector_search(pos_neighbor[testnode], nodenow) == 1) {
						commoncount++;
					}
				}
			}
			int baddibo = 0;
			if (commoncount == 0 && isize > 0 && pos_neighbor[testnode].size() > 0) {
				empty_connection.push_back(make_pair(i, testnode));
			}
		}
	}


	//--------- Finding capsomer positions for locations of empty nodes ----------------//

	for (int q = 0; q < empty_connection.size(); q++) {
		int emptynode1 = empty_connection[q].first;
		int emptynode2 = empty_connection[q].second;
		int i = emptynode1;
		Real isize = pos_neighbor[i].size();
		vector<int>finallist;
		//positive angles and negative angles data structure
		map<Real, Real>p_angles;
		map<Real, Real>n_angles;
		//Building vector from i-th node to 1st neighbour node
		int t_ind = emptynode2;
		Vector3d ot(pos.at(t_ind)(0) - pos.at(i)(0), pos.at(t_ind)(1) - pos.at(i)(1), pos.at(t_ind)(2) - pos.at(i)(2));
		Real testminimum = 999;
		int test_n = 1;
		int change = 0;
		Real testangle;
		Real testangletest;
		//find angle with the initial vector with the rest of the nodes and find the minimum angle constructing vector
		for (int j = 0; j < isize; j++) {
			int n_ind = pos_neighbor[i][j];
			if (n_ind == emptynode2)continue;
			Vector3d on(pos.at(n_ind)(0) - pos.at(i)(0), pos.at(n_ind)(1) - pos.at(i)(1), pos.at(n_ind)(2) - pos.at(i)(2));
			pair<Real, Real>resultangles;
			resultangles = angle(ot, on);
			if (change != 0)continue;
			if (j == 1)testangletest = resultangles.second;
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
		Vector3d on(pos.at(test_n)(0) - pos.at(i)(0), pos.at(test_n)(1) - pos.at(i)(1), pos.at(test_n)(2) - pos.at(i)(2));
		Vector3d t_cross = ot.cross(on);
		if (change != 0) {
			p_angles[testangle] = test_n;
		}
		for (int j = 0; j < isize; j++) {
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
		vector_p = sort(p_angles);
		vector_n = sort(n_angles);
		int a1, a2;
		if (vector_p.size() + vector_n.size() < 2)continue;
		if (vector_n.size() == 0) {
			a1 = vector_p[0];
			a2 = vector_p[1];
		}
		else if (vector_p.size() == 0) {
			a1 = vector_n[0];
			a2 = vector_n[1];
		}
		else {
			a1 = vector_p[0];
			a2 = vector_n[0];
		}
		int i1 = emptynode2;
		Real isize1 = pos_neighbor[i1].size();
		//positive angles and negative angles data structure
		map<Real, Real>p_angles1;
		map<Real, Real>n_angles1;
		//Building vector from i-th node to 1st neighbour node
		int t_ind1 = emptynode1;
		Vector3d ot1(pos.at(t_ind1)(0) - pos.at(i1)(0), pos.at(t_ind1)(1) - pos.at(i1)(1), pos.at(t_ind1)(2) - pos.at(i1)(2));
		Real testminimum1 = 999;
		int test_n1 = 1;
		int change1 = 0;
		Real testangle1;
		Real testangletest1;
		//find angle with the initial vector with the rest of the nodes and find the minimum angle constructing vector
		for (int j = 1; j < isize1; j++) {
			int n_ind1 = pos_neighbor[i1][j];
			if (n_ind1 == emptynode2)continue;
			Vector3d on1(pos.at(n_ind1)(0) - pos.at(i1)(0), pos.at(n_ind1)(1) - pos.at(i1)(1), pos.at(n_ind1)(2) - pos.at(i1)(2));
			pair<Real, Real>resultangles1;
			resultangles1 = angle(ot1, on1);
			if (change != 0)continue;
			if (j == 1)testangletest1 = resultangles1.second;
			//------------Angle between two vectors - two outputs - first-> smallest angle(regardless of direction) ;  second-> positive angle(0,PI)-----------//
			if (resultangles1.first > 0) {
				if (resultangles1.first < testminimum1) {
					testminimum1 = resultangles1.first;
					testangle1 = resultangles1.second;
					test_n1 = n_ind1;
					change1++;
				}
			}
		}
		Vector3d on1(pos.at(test_n1)(0) - pos.at(i1)(0), pos.at(test_n1)(1) - pos.at(i1)(1), pos.at(test_n1)(2) - pos.at(i1)(2));
		Vector3d t_cross1 = ot1.cross(on1);
		if (change1 != 0) {
			p_angles1[testangle1] = test_n1;
		}
		for (int j = 0; j < isize1; j++) {
			int n_ind1 = pos_neighbor[i1][j];
			if (n_ind1 != test_n1) {

				Vector3d on_f1(pos.at(n_ind1)(0) - pos.at(i1)(0), pos.at(n_ind1)(1) - pos.at(i1)(1), pos.at(n_ind1)(2) - pos.at(i1)(2));
				pair<Real, Real>resultangles1;
				resultangles1 = angle(ot1, on_f1);
				Vector3d n_cross1 = ot1.cross(on_f1);
				if (t_cross1.dot(n_cross1) >= 0) {
					if (isnan(resultangles1.second)) {

						p_angles1[180] = n_ind1;
					}
					else {
						p_angles1[resultangles1.second] = n_ind1;

					}
				}
				else {
					n_angles1[resultangles1.second] = n_ind1;

				}
			}
		}
		vector<Real>vector_p1;
		vector<Real>vector_n1;
		vector_p1 = sort(p_angles1);
		vector_n1 = sort(n_angles1);
		if (vector_p1.size() + vector_n1.size() < 2)continue;
		int b1, b2;
		if (vector_n1.size() == 0) {
			b1 = vector_p1[0];
			b2 = vector_p1[1];
		}
		else if (vector_p1.size() == 0) {
			b1 = vector_n1[0];
			b2 = vector_n1[1];
		}
		else {
			b1 = vector_p1[0];
			b2 = vector_n1[0];
		}
		if (dist(pos[a1], pos[b1]) < dist(pos[a1], pos[b2])) {
			vector<int>sortfound;
			sortfound.push_back(a1);
			sortfound.push_back(b1);
			sortfound.push_back(emptynode2);
			sortfound.push_back(emptynode1);
			sort(sortfound.begin(), sortfound.begin() + 4);
			int mile = 0;
			for (int find = 0; find < empcount; find++) {
				if (Emp[find].alreadyfound[0] == sortfound[0] && Emp[find].alreadyfound[1] == sortfound[1] && Emp[find].alreadyfound[2] == sortfound[2] && Emp[find].alreadyfound[3] == sortfound[3]) {
					mile++;
				}
			}
			Vector3d Aa((pos[a1][0] + pos[b1][0] + pos[emptynode1][0] + pos[emptynode2][0]) / 4, (pos[a1][1] + pos[b1][1] + pos[emptynode1][1] + pos[emptynode2][1]) / 4, (pos[a1][2] + pos[b1][2] + pos[emptynode1][2] + pos[emptynode2][2]) / 4);
			if (mile == 0) {
				Emp[empcount].alreadyfound = sortfound;
				Emp[empcount].Point = Aa;
				empcount++;
			}
			vector<int>sortfound1;
			sortfound1.push_back(a2);
			sortfound1.push_back(b2);
			sortfound1.push_back(emptynode2);
			sortfound1.push_back(emptynode1);
			sort(sortfound1.begin(), sortfound1.begin() + 4);
			mile = 0;
			for (int find = 0; find < empcount; find++) {
				if (Emp[find].alreadyfound[0] == sortfound1[0] && Emp[find].alreadyfound[1] == sortfound1[1] && Emp[find].alreadyfound[2] == sortfound1[2] && Emp[find].alreadyfound[3] == sortfound1[3]) {
					mile++;
				}

			}
			Vector3d AAa((pos[a2][0] + pos[b2][0] + pos[emptynode1][0] + pos[emptynode2][0]) / 4, (pos[a2][1] + pos[b2][1] + pos[emptynode1][1] + pos[emptynode2][1]) / 4, (pos[a2][2] + pos[b2][2] + pos[emptynode1][2] + pos[emptynode2][2]) / 4);
			if (mile == 0) {
				Emp[empcount].alreadyfound = sortfound1;
				Emp[empcount].Point = AAa;
				empcount++;
			}
		}
		else {
			vector<int>sortfound;
			sortfound.push_back(a1);
			sortfound.push_back(b2);
			sortfound.push_back(emptynode2);
			sortfound.push_back(emptynode1);
			sort(sortfound.begin(), sortfound.begin() + 4);
			int mile = 0;
			for (int find = 0; find < empcount; find++) {
				if (Emp[find].alreadyfound[0] == sortfound[0] && Emp[find].alreadyfound[1] == sortfound[1] && Emp[find].alreadyfound[2] == sortfound[2] && Emp[find].alreadyfound[3] == sortfound[3]) {
					mile++;
				}
			}
			Vector3d Aa((pos[a1][0] + pos[b2][0] + pos[emptynode1][0] + pos[emptynode2][0]) / 4, (pos[a1][1] + pos[b2][1] + pos[emptynode1][1] + pos[emptynode2][1]) / 4, (pos[a1][2] + pos[b2][2] + pos[emptynode1][2] + pos[emptynode2][2]) / 4);
			if (mile == 0) {
				Emp[empcount].alreadyfound = sortfound;
				Emp[empcount].Point = Aa;
				empcount++;
			}
			vector<int>sortfound1;
			sortfound1.push_back(a2);
			sortfound1.push_back(b1);
			sortfound1.push_back(emptynode2);
			sortfound1.push_back(emptynode1);
			sort(sortfound1.begin(), sortfound1.begin() + 4);
			mile = 0;
			for (int find = 0; find < empcount; find++) {
				if (Emp[find].alreadyfound[0] == sortfound1[0] && Emp[find].alreadyfound[1] == sortfound1[1] && Emp[find].alreadyfound[2] == sortfound1[2] && Emp[find].alreadyfound[3] == sortfound1[3]) {
					mile++;
				}
			}
			Vector3d AAa((pos[a2][0] + pos[b1][0] + pos[emptynode1][0] + pos[emptynode2][0]) / 4, (pos[a2][1] + pos[b1][1] + pos[emptynode1][1] + pos[emptynode2][1]) / 4, (pos[a2][2] + pos[b1][2] + pos[emptynode1][2] + pos[emptynode2][2]) / 4);
			if (mile == 0) {
				Emp[empcount].alreadyfound = sortfound1;
				Emp[empcount].Point = AAa;
				empcount++;
			}
		}
	}

	for (int i = 0; i < n_part; i++) {
		int nowsize = pos_neighbor[i].size();
		for (int j = 0; j < nowsize; j++) {
			int nownode = pos_neighbor[i][j];
			if (vector_search(pos_neighbor[nownode], i) == 0) {
				pos_neighbor[nownode].push_back(i);
			}
		}
	}

	//-------------------------------DataStructure for edge length---------------//
	vector<vector<Vector3d>>edgedone;
	vector<Real>edges;


	//-------------------------------Build-up of capsomers-----------------------------//
	for (int i = 0; i < nonrepeat_part; i++) {
		if (pos_neighbor[i].size() > 3) {       //.....................Only making capsomer cells for cell valance 4 and greater (can be changed if needed) ..................//
			totalpoints += pos_neighbor[i].size();
			totalcaps++;
			vtknodes.push_back(i);
			//}
		}
	}

	ofstream f;
	ofstream fo;
	ofstream ft;
	ofstream fp;
	int trianglecount = 0;
	f.open(output_capsomer);															//.......................This file creates the capsomers in vtk..........................//
	f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	f << nonrepeat_part + totalpoints;
	f << " float\n";

	fo.open(output_outline);															//.......................This file creates the capsomers' outline in vtk..........................//
	fo << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	fo << nonrepeat_part + totalpoints;
	fo << " float\n";

	ft.open(input_file + "_triangle.vtk");												//.......................This file creates the triangulation among the particles in vtk..........................//
	ft << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	ft << nonrepeat_part + totalpoints;
	ft << " float\n";

	fp.open(input_file + "data.txt");													//.......................This file creates the final neighboring particles information for each particle in vtk..........................//

	for (int i = 0; i < nonrepeat_part; i++) {
		f << pos.at(i)(0) << " " << pos.at(i)(1) << " " << pos.at(i)(2) << endl;
		fo << pos.at(i)(0) << " " << pos.at(i)(1) << " " << pos.at(i)(2) << endl;
	}

	for (int i = 0; i < nonrepeat_part; i++) {
		Real isize = pos_neighbor[i].size();
		if (isize > 3) {
			vector<int>finalcaps;
			vector<int>finallist;
			finallist = neighborordering(i, pos_neighbor, pos);
			for (int fsize = 0; fsize < finallist.size() - 1; fsize++) {
				fp << finallist[fsize] << " ";
			}
			fp << endl;
			vector<Vector3d>bcs;
			for (int j = 0; j < isize; j++) {
				int alreadypushed = 0;
				for (int empiter = 0; empiter < empcount; empiter++) {
					if (vector_search(Emp[empiter].alreadyfound, i) == 1) {
						if (vector_search(Emp[empiter].alreadyfound, finallist.at(j)) == 1) {
							if (vector_search(Emp[empiter].alreadyfound, finallist.at(j + 1)) == 1) {
								f << Emp[empiter].Point(0) << " " << Emp[empiter].Point(1) << " " << Emp[empiter].Point(2) << " " << endl;
								fo << Emp[empiter].Point(0) << " " << Emp[empiter].Point(1) << " " << Emp[empiter].Point(2) << " " << endl;
								alreadypushed++;
								break;
							}
						}
					}

				}
				Vector3d test = findbc(pos.at(i), pos.at(finallist.at(j)), pos.at(finallist.at(j + 1)));
				bcs.push_back(test);
				ft << pos.at(i)(0) << " " << pos.at(i)(1) << " " << pos.at(i)(2) << endl;
				ft << pos.at(finallist.at(j))(0) << " " << pos.at(finallist.at(j))(1) << " " << pos.at(finallist.at(j))(2) << endl;
				ft << pos.at(finallist.at(j + 1))(0) << " " << pos.at(finallist.at(j + 1))(1) << " " << pos.at(finallist.at(j + 1))(2) << endl;
				trianglecount++;
				if (alreadypushed == 0) {
					f << test(0) << " " << test(1) << " " << test(2) << " " << endl;
					fo << test(0) << " " << test(1) << " " << test(2) << " " << endl;
				}
			}
			int bcsize = bcs.size();
			for (int bi = 0; bi < bcsize; bi++) {
				Vector3d A = bcs[bi];
				Vector3d B = bcs[(bi + 1) % bcsize];
				bool found = 0;
				for (int eds = 0; eds < edgedone.size(); eds++) {
					if ((A == edgedone[eds][0] && B == edgedone[eds][1]) || (A == edgedone[eds][1] && B == edgedone[eds][0])) {
						found = 1;
					}
				}
				if (found == 0) {
					edges.push_back(dist(A, B));
					edgedone.push_back({ A,B });
				}
			}
		}
	}
	Real distances = 0;
	int disttotal = 0;
	for (int i = 0; i < nonrepeat_part; i++) {
		int isize = pos_neighbor[i].size();
		for (int j = 0; j < isize; j++) {
			distances += dist(pos[i], pos[pos_neighbor[i][j]]);
			disttotal++;
		}
	}
	Real avgdistpart = distances / disttotal;
	Real avgedgelen = 0.5 * avgdistpart;

	int totalbad = 0;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i] < 0.5 * avgedgelen)totalbad++;
	}

	//-------------------------Post processing VTK files for visualisation--------------------------//

	f << "\n";
	f << "CELLS ";
	f << totalcaps << " " << totalpoints + 2 * totalcaps << endl;
	fo << "\n";
	fo << "CELLS ";

	fo << totalcaps << " " << totalpoints + 2 * totalcaps << endl;

	ft << "\n";
	ft << "CELLS ";
	ft << trianglecount << " " << 4 * trianglecount << endl;

	int triangle = 0;
	for (int i = 0; i < trianglecount; i++) {
		ft << "3 " << triangle << " " << triangle + 1 << " " << triangle + 2 << endl;
		triangle = triangle + 3;
	}
	ft << "\n";
	ft << "CELL_TYPES " << trianglecount << endl;
	for (int i = 0; i < trianglecount; i++) {
		ft << "5" << endl;
	}

	int cellnow = nonrepeat_part;
	for (int i = 0; i < totalcaps; i++) {
		int nodenow = vtknodes.at(i);
		int csize = pos_neighbor[nodenow].size();
		f << csize + 1 << " ";
		//" " << nodenow <<
		fo << csize + 1 << " ";
		for (int j = 0; j < csize; j++) {
			f << cellnow << " ";
			fo << cellnow << " ";
			cellnow++;
		}
		f << cellnow - csize << endl;
		fo << cellnow - csize << endl;
	}
	f << "\n";
	f << "CELL_TYPES " << totalcaps << endl;
	fo << "\n";
	fo << "CELL_TYPES " << totalcaps << endl;
	for (int i = 0; i < totalcaps; i++) {
		f << "7" << endl;
		fo << "7" << endl;
	}
	f << "\nCELL_DATA " << totalcaps << "\n";
	f << "SCALARS neighbour float\n";
	f << "LOOKUP_TABLE default\n";
	for (int i = 0; i < totalcaps; i++) {
		int nodenow = vtknodes.at(i);
		int csize = pos_neighbor[nodenow].size();
		if (csize == 5) {
			f << "-3\n";
		}
		else if (csize == 7) {
			f << "3\n";
		}
		else if (csize == 6) {
			f << "0\n";
		}
		else if (csize == 4) {
			f << "-2\n";
		}
		else if (csize == 8) {
			f << "2\n";
		}
		else if (csize == 9) {
			f << "1\n";
		}
		else if (csize == 10) {
			f << "0.5\n";
		}
		else if (csize == 11) {
			f << "0.2\n";
		}
		else if (csize == 12) {
			f << "0.15\n";
		}
		else {
			f << "0.3\n";
		}
	}
	cout << triangle << " is the number of vertices in triangulation information file" << endl;			//................................Use this as the number of particles in the triangulation file............................//
	f.close();
	fo.close();
	ft.close();
	fp.close();
	clock_t stoptime = clock();
	Real sum = ((double)(stoptime - starttime) / CLOCKS_PER_SEC);
	cout << "RunningTime: " << sum << endl;
}
