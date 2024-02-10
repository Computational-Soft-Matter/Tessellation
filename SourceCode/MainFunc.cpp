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

	//.............................User Inputs.......................................//
	string input_file = "Figure18br";		// File that contains the 3D coordinates of the colloidal particles
	string input_extension = ".vtk";		// Extension of the file that contains the 3D coordinates of the colloidal particles
	string neighbor_input = "NULL";			// If neighbor list has been calculated previously, then instead of NULL, the file will contain the index of a particle, number of neighbors and the radius of search to acquire those neighbors
	const int n_part = 1000;			// Total Number of Particles in the inputfile
	const int nonrepeat_part = 1000;		// Non repeating number of Particles in the inputfile (This is equal to n_part for non periodic boundary cases) 	
	bool swdel1 = 1;				// 1 if deleting process-1 is activated, 0 otherwise
	bool swdel2 = 1;				// 1 if deleting process-2 is activated, 0 otherwise
	bool swhole = 1;				// 1 if holes patch up is activated, 0 otherwise
	bool sweempty = 1;				// 1 if empty regions fix is activated, 0 otherwise
	

	//.............................Initializations.......................................//
	string output_capsomer = input_file + "caps.vtk";		// Output file name that has the information of the capsomer cells
	string output_outline = input_file + "otln.vtk";		// Output file name that has the information of capsomer outline
	clock_t starttime = clock();					// Function to calculate algorithm's processor run time
	int totalpoints = 0;						// Total number of points in the final VTK file containing all the capsomers (this includes the initial set of points and the vertices of the capsomers)
	founds Emp[3000];						// Data structure for averaging points in empty region - (container to store the vertices of quadrilaterals that uses average capsomer vertices position described in section 2.1.4
	int empcount = 0;						// Data structure for averaging points in empty region
	int totalcaps = 0;						// Total number of capsomers in the final vtk file
	vector<Real>n_neighbor;						// Neighbour information of target particles - first -> number of n_neighbor, second-> radius range of search
	vector<int>vtknodes;						// Particle IDs for final VTK file (including capsomer vertices)
	vector<Real>avg_neighbor;					// Average distances for n_neighbor
	vector<pair<int, int>>single_edges;				// Pairs of single edge connections
	map<int, int>marked;						// Marked node which has single edge connections and changed capsomers
	vector<vector<int>>pos_neighbor;				// Neighbour nodes of target particles
	vector<Vector3d>pos;						// Position of particles

	
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
	vector<vector<int>>RcPos;
	if (neighbor_input == "NULL")RcPos=rccontainer(pos, n_part);

	//-------------Neighbour search algorithm using moving window (Sec. 2.1.1 in manuscript - Algorithm 1)-----------------//
	if (neighbor_input == "NULL") {
		for (int i = 0; i <n_part; i++) {
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
			vector<Real>distances;
			for (int j = 0; j < RcPos[i].size(); j++) {
				Real r = dist(pos.at(i), pos.at(RcPos[i][j]));
				distances.push_back(r);
			}
			sort(distances.begin(), distances.end());
			int leftind, rightind;
			for (int j = 0; j < distances.size(); j++) {
				if (distances[j] >= init_r - delr / 2.0) {
					leftind = j;
				}
				if (distances[j] > init_r + delr / 2.0) {
					rightind = j;
					break;
				}
			}
			int cont = rightind - leftind;
			bool start = 1;
			while (isstop != 2) {
				int count_n = 0;
				if (start == 1) {
					count_n = cont;
					start = 0;
				}
				else {
					bool stopwhile = 0;
					while (distances[leftind] < (init_r - delr / 2.0) && stopwhile==0) {
						leftind++;
						if (leftind == distances.size()) {
							stopwhile = 1;
							leftind = 0;
						}
					}
					stopwhile = 0;
					while (stopwhile == 0 && distances[rightind] < (init_r + delr / 2.0)) {
						rightind++;
						if (rightind == distances.size()) {
							stopwhile = 1;
							rightind = 0;
						}
					}
					count_n = rightind - leftind;
				}
				int countdist = count_n - maxn;
				lastcount_n = count_n;
				if (maxn < count_n) {
					maxn = count_n;
					maxr = init_r;
				}
				list_p.push_back(count_n);
				list_n.push_back(init_r);
				if (init_r >= (R_c-delr/2.0) || countdist < NC) {
					isstop = 2;
				}
				it++;
				init_r += dr;
			}
			n_neighbor.push_back(maxr);
			//cout << i << " " << maxn << " " << maxr << endl;	//..................Printing the radius of search for each particle and the number of particles found as initial neighbors for debugging...............// 
		}
	}
	
	//-------------Build-up of particle neighbour containers------------//
	for (int i = 0; i < n_part; i++) {
		pos_neighbor.push_back(vector<int>());
	}
	for (int i = 0; i < n_part; i++) {
		Real range = n_neighbor.at(i);
		if (NC == -100) {    // if NC==-100, then no NC will be used, rather only R_c will be used as cut-off
			range = R_c - delr / 2.0;;
		}
		for (int j = 0; j < RcPos[i].size(); j++) {
			Real r = dist(pos.at(i), pos.at(RcPos[i][j]));
			if (r <= range + delr / 2.0 && i != RcPos[i][j]) {
				int isize = pos_neighbor[i].size();
				int jsize = pos_neighbor[RcPos[i][j]].size();
				int ci = 0, cj = 0;
				for (int k = 0; k < isize; k++) {
					if (pos_neighbor[i].at(k) == RcPos[i][j]) {
						ci++;
					}
				}
				for (int k = 0; k < jsize; k++) {
					if (pos_neighbor[RcPos[i][j]].at(k) == i) {
						cj++;
					}
				}
				if (ci == 0) {
					pos_neighbor[i].push_back(RcPos[i][j]);
				}
				if (cj == 0) {
					pos_neighbor[RcPos[i][j]].push_back(i);
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


	//-----------------Deleting overlapping particle edge-connections- Method 1 (Sec. 2.1.2 in manuscript - Algorithm 2)--------------//
	vector<pair<Real,vector<int>>>corrbef;								
	map<pair<int, int>, bool>donebef;									
	int Corrcount = 0;
	if (swdel1 == 1) {
		for (int i = 0; i < n_part; i++) {
			int i_size = pos_neighbor[i].size();						
			for (int j = 0; j < i_size; j++) {
				int node_j = pos_neighbor[i].at(j);						
				if (donebef[make_pair(node_j, i)] == 1)continue;
				vector<int>deletelist;
				for (int k = 0; k < i_size; k++) {
					int nodenow = pos_neighbor[i].at(k);
					if (nodenow != i || nodenow != node_j) {
						if (vector_search(pos_neighbor[i], nodenow) == 1 && vector_search(pos_neighbor[node_j], nodenow) == 1) {
							deletelist.push_back(nodenow);				
						}
					}
				}
				map<Real, Real>p_angles;								
				map<Real, Real>n_angles;								
				int t_ind = node_j;
				Vector3d ot(pos.at(t_ind)(0) - pos.at(i)(0), pos.at(t_ind)(1) - pos.at(i)(1), pos.at(t_ind)(2) - pos.at(i)(2));
				Real testminimum = 999;									
				int test_n = 1;											
				int change = 0;													
				Real testangle;
				Real testangletest;
				for (int j = 0; j < deletelist.size(); j++) {
					int n_ind = deletelist[j];
					Vector3d on(pos.at(n_ind)(0) - pos.at(i)(0), pos.at(n_ind)(1) - pos.at(i)(1), pos.at(n_ind)(2) - pos.at(i)(2));
					pair<Real, Real>resultangles;
					resultangles = angle(ot, on);     
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
				Vector3d on(pos.at(test_n)(0) - pos.at(i)(0), pos.at(test_n)(1) - pos.at(i)(1), pos.at(test_n)(2) - pos.at(i)(2));
				Vector3d t_cross = ot.cross(on);
				if (change != 0) {
					p_angles[testangle] = test_n;
				}
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
				vector_p = sort(p_angles);							
				vector_n = sort(n_angles);							
				int psize = vector_p.size();
				if (psize > 1) {
					for (int ite = 0; ite < vector_p.size(); ite++) {
						int pnow = vector_p[ite];
						if (dist(pos.at(i), pos.at(pnow)) - avg_neighbor[i] < dist(pos.at(node_j), pos.at(pnow)) - avg_neighbor[node_j]) {
							vector<int>tem;
							tem = { node_j,i,pnow };
							corrbef.push_back(make_pair((dist(pos.at(node_j), pos.at(pnow)) - avg_neighbor[node_j]) + (dist(pos.at(i), pos.at(pnow)) - avg_neighbor[i]), tem));
						}
						else {
							vector<int>tem;
							tem = { i,node_j,pnow };
							corrbef.push_back(make_pair((dist(pos.at(node_j), pos.at(pnow)) - avg_neighbor[node_j]) + (dist(pos.at(i), pos.at(pnow)) - avg_neighbor[i]), tem));
						}
					}
				}
				int nsize = vector_n.size();
				if (nsize > 1) {
					for (int ite = 0; ite < vector_n.size(); ite++) {
						int nnow = vector_n[ite];
						if ((dist(pos.at(i), pos.at(nnow)) - avg_neighbor[i]) < (dist(pos.at(node_j), pos.at(nnow)) - avg_neighbor[node_j])) {
							vector<int>tem;
							tem = { node_j,i,nnow };
							corrbef.push_back(make_pair((dist(pos.at(node_j), pos.at(nnow)) - avg_neighbor[node_j]) + (dist(pos.at(i), pos.at(nnow)) - avg_neighbor[i]), tem));
						}
						else {
							vector<int>tem;
							tem = { i,node_j,nnow };
							corrbef.push_back(make_pair((dist(pos.at(node_j), pos.at(nnow)) - avg_neighbor[node_j]) + (dist(pos.at(i), pos.at(nnow)) - avg_neighbor[i]), tem));
						}
					}
				}
			}
		}
	}
	sort(corrbef.begin(), corrbef.end());
	if (swdel1 == 1) {
		while(corrbef.size()){
			int cornow = corrbef.size() - 1;
			int i = corrbef[cornow].second[0];
			int i_size = pos_neighbor[i].size();
			int node_j = corrbef[cornow].second[1];						
			int del = corrbef[cornow].second[2];
			vector<int>deletelist;
			for (int k = 0; k < i_size; k++) {
				int nodenow = pos_neighbor[i].at(k);
				if (nodenow != i || nodenow != node_j) {
					if (vector_search(pos_neighbor[i], nodenow) == 1 && vector_search(pos_neighbor[node_j], nodenow) == 1) {
						deletelist.push_back(nodenow);				
					}
				}
			}
			
			map<Real, Real>p_angles;								
			map<Real, Real>n_angles;								
			int t_ind = node_j;
			Vector3d ot(pos.at(t_ind)(0) - pos.at(i)(0), pos.at(t_ind)(1) - pos.at(i)(1), pos.at(t_ind)(2) - pos.at(i)(2));
			Real testminimum = 999;									
			int test_n = 1;											
			int change = 0;													
			Real testangle;
			Real testangletest;
			for (int j = 0; j < deletelist.size(); j++) {
				int n_ind = deletelist[j];
				Vector3d on(pos.at(n_ind)(0) - pos.at(i)(0), pos.at(n_ind)(1) - pos.at(i)(1), pos.at(n_ind)(2) - pos.at(i)(2));
				pair<Real, Real>resultangles;
				resultangles = angle(ot, on);     
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
			Vector3d on(pos.at(test_n)(0) - pos.at(i)(0), pos.at(test_n)(1) - pos.at(i)(1), pos.at(test_n)(2) - pos.at(i)(2));
			Vector3d t_cross = ot.cross(on);
			if (change != 0) {
				p_angles[testangle] = test_n;
			}
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
			bool deleted = 0;
			vector<Real>vector_p;
			vector<Real>vector_n;
			vector_p = sort(p_angles);							
			vector_n = sort(n_angles);							
			int psize = vector_p.size();
			int nsize = vector_n.size();
			bool pyes = 0;
			for (int ite = 0; ite < psize; ite++) {
				int pnow = vector_p[ite];
				if (pnow == del) {
					pyes = 1;
					break;
				}
			}
			if (pyes == 1 && psize > 1) {
				deleted = 1;
				int deletnow = del;
				int dltfr = i;
				int deletat=-5;
				for (int y = 0; y < pos_neighbor[dltfr].size(); y++) {
					if (deletnow == pos_neighbor[dltfr].at(y)) {
						deletat = y;
					}
				}
				pos_neighbor[dltfr].erase(pos_neighbor[dltfr].begin() + deletat);
				int dltfr2 = del;
				int deletnow2 = i;
				int deletat2=-5;
				for (int y = 0; y < pos_neighbor[dltfr2].size(); y++) {
					if (deletnow2 == pos_neighbor[dltfr2].at(y)) {
						deletat2 = y;
					}
				}
				pos_neighbor[dltfr2].erase(pos_neighbor[dltfr2].begin() + deletat2);
			}
			
			bool nyes = 0;
			for (int ite = 0; ite < nsize; ite++) {
				int nnow = vector_n[ite];
				if (nnow == del) {
					pyes = 1;
					break;
				}
			}
			if (nyes == 1 && nsize > 1) {
				deleted = 1;
				int deletnow = del;
				int dltfr = i;
				int deletat=-5;
				for (int y = 0; y < pos_neighbor[dltfr].size(); y++) {
					if (deletnow == pos_neighbor[dltfr].at(y)) {
						deletat = y;
					}
				}
				pos_neighbor[dltfr].erase(pos_neighbor[dltfr].begin() + deletat);
				int dltfr2 = del;
				int deletnow2 = i;
				int deletat2=-5;
				for (int y = 0; y < pos_neighbor[dltfr2].size(); y++) {
					if (deletnow2 == pos_neighbor[dltfr2].at(y)) {
						deletat2 = y;
					}
				}
				pos_neighbor[dltfr2].erase(pos_neighbor[dltfr2].begin() + deletat2);
			}
			corrbef.pop_back();
			if (deleted == 1) {
				for (int delse = 0; delse < corrbef.size(); delse++) {
					if (vector_search(corrbef[delse].second, i) == 1) {
						Real dfavg = maxdiffneighbor(corrbef[delse].second[0], corrbef[delse].second[1], corrbef[delse].second[2], pos_neighbor, pos);
						corrbef[delse].first = dfavg;
					}
					else if ((vector_search(corrbef[delse].second, del) == 1)) {
						Real dfavg = maxdiffneighbor(corrbef[delse].second[0], corrbef[delse].second[1], corrbef[delse].second[2], pos_neighbor, pos);
						corrbef[delse].first = dfavg;
					}
				}
				sort(corrbef.begin(), corrbef.end());
			}
		}	
	}	
	//-----------------Deleting overlapping particle edge-connections- Method 2 (Sec. 2.1.2 in manuscript - Algorithm 3) --------------//
	if (swdel2 == 1) {
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
						int deletenow;
						int deletfrom;
						if (dist(pos[pastcheck], pos[nextcheck]) < dist(pos[i], pos[finallist[j]])) {
							deletenow = finallist[j];
							deletfrom = i;
						}
						else {
							deletenow = pastcheck;
							deletfrom = nextcheck;
						}
						if (vector_search(pos_neighbor[deletfrom], deletenow) == 0)continue;
						int deletat = find_ind(pos_neighbor[deletfrom], deletenow);
						pos_neighbor[deletfrom].erase(pos_neighbor[deletfrom].begin() + deletat);
						int deletat2 = find_ind(pos_neighbor[deletenow], deletfrom);
						pos_neighbor[deletenow].erase(pos_neighbor[deletenow].begin() + deletat2);
					}
				}
			}
		}
	}
	
	//-----------------Search for pairs creating gaps (Sec. 2.1.3 in manuscript - Algorithm 5) --------------//
	vector<vector<int>>single_conn;								
	vector<pair<int, int>>single_edge;							
	map<int, bool>visited;										
	vector<int>listofsingles;									
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
	//-----------------Forming particles chains surrounding gaps in the tessellation (Sec. 2.1.3 in manuscript - Algorithm 6) --------------//
	map<int, bool>worked;									
	vector<vector<int>>holes;								
	for (int i = 0; i < listofsingles.size(); i++) {
		int workingnow = listofsingles.at(i);				 
		if (single_conn.at(workingnow).size() >= 2) {
			if ((worked[workingnow] == 0 && single_conn.at(workingnow).size() >= 2) || (worked[workingnow] == 1 && single_conn.at(workingnow).size() > 2)) {
				vector<int>hole;							
				hole.push_back(workingnow);
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

				int stop = 0;					
				while (next != workingnow && stop == 0) {
					if (single_conn.at(next).size() < 2) {
						stop = 1;											
					}
					else if (single_conn.at(next).size() == 2) {			
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
					else {													
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
	//-----------------Repairing gaps in the tessellation (Sec. 2.1.3 in manuscript - Algorithm 7) --------------//
	if (swhole == 1) {
		for (int i = 0; i < holes.size(); i++) {
			int holesize = holes.at(i).size() - 1;				
			vector<vector<int>>fixes;							
			storedelete NN[1000];
			int NNco = 0;
			for (int j = 0; j < holes.at(i).size() - 1; j++) {
				int nodenow = holes.at(i).at(j);
				vector<int>possiblefixes;						
				possiblefixes.push_back(nodenow);				
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
						NN[NNco].dist = dist(pos.at(nodenow), pos.at(holes[i][k]));
						NN[NNco].deletefrom = nodenow;
						NN[NNco].deletenode = holes[i][k];
						NNco++;
						disttot += r;
						total++;
					}
				}
			}
			int actualnnco = NNco;
			for (int ite = NNco; ite < 1000; ite++) {
				NN[ite].dist = 0;
				NN[ite].deletefrom = 0;
				NN[ite].deletenode = 0;
			}
			sort(NN, NN + 1000, CompareDatadlt);
			for (int j = NNco - 1; j >= 0; j -= 2) {
				int fixing_now = NN[j].deletefrom;
				int fixing_now_ind = find_ind(holes.at(i), fixing_now);
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
					int currentnode_fix = singleedge_fixingnow.at(l).first;    
					int adjacentnode = singleedge_fixingnow.at(l).second;	   
					int deleted = 0;
					int currentsize = pos_neighbor[currentnode_fix].size();
					int commoncount = 0;
					for (int k = 0; k < currentsize + deleted; k++) {
						int nodenow = pos_neighbor[currentnode_fix].at(k);		
						if (nodenow != currentnode_fix && nodenow != adjacentnode) {
							if (vector_search(pos_neighbor[currentnode_fix], nodenow) == 1 && vector_search(pos_neighbor[adjacentnode], nodenow) == 1) {
								commoncount++;
							}
						}
					}
					if (currentsize > 2 && pos_neighbor[adjacentnode].size() > 2) {				
						int confirmfixednode = -5;							
						int currentfix = NN[j].deletenode;
						int checkvalueforfix = 0;
						if (currentfix != currentnode_fix && vector_search(pos_neighbor.at(currentnode_fix), currentfix) == 0)
						{
							if (dist(pos.at(currentnode_fix), pos.at(currentfix)) < 1.0 * R_c) {
								confirmfixednode = currentfix;					
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

	//-----------------Averaging Cells in Empty Regions (Sec. 2.1.4 in manuscript - Algorithm 8) --------------//
	vector<pair<int, int>>empty_connection;								
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
	if (sweempty == 1) {
		for (int q = 0; q < empty_connection.size(); q++) {
			int emptynode1 = empty_connection[q].first;
			int emptynode2 = empty_connection[q].second;
			if (emptynode1 > nonrepeat_part && emptynode2 > nonrepeat_part)continue;     
			int i = emptynode1;
			Real isize = pos_neighbor[i].size();
			vector<int>finallisti;
			int a1, a2;
			finallisti = neighborordering(i, pos_neighbor, pos);
			for (int flo = 0; flo < finallisti.size(); flo++) {
				if (finallisti[flo] == emptynode2) {
					if (flo - 1 == -1)a1 = finallisti[finallisti.size() - 2];
					else a1 = finallisti[flo - 1];
					a2 = finallisti[flo + 1];
					break;
				}
			}
			int i1 = emptynode2;
			Real isize1 = pos_neighbor[i1].size();			
			vector<int>finallisti1;
			int b1, b2;
			finallisti1 = neighborordering(i1, pos_neighbor, pos);
			for (int fl = 0; fl < finallisti1.size(); fl++) {
				if (finallisti1[fl] == i) {
					if (fl - 1 == -1)b1 = finallisti1[finallisti1.size() - 2];
					else b1 = finallisti1[fl - 1];
					b2 = finallisti1[fl + 1];
					break;
				}
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

	//-----------------Construction of tessellation cells (Sec. 2.1.5 in manuscript - Algorithm 9) --------------//
	for (int i = 0; i < nonrepeat_part; i++) {
		if (pos_neighbor[i].size() > 3) {       //.....................Only making capsomer cells for cell valance 4 and greater (can be changed if needed) ..................//
			totalpoints += pos_neighbor[i].size();
			totalcaps++;
			vtknodes.push_back(i);
		}
	}
	clock_t starttimePRINT = clock();
	ofstream f;
	ofstream fo;
	int trianglecount = 0;
	f.open(output_capsomer);															//.......................This file creates the capsomers in vtk..........................//
	f << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	f << nonrepeat_part + totalpoints;
	f << " float\n";

	fo.open(output_outline);															//.......................This file creates the capsomers' outline in vtk..........................//
	fo << "# vtk DataFile Version 1.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ";
	fo << nonrepeat_part + totalpoints;
	fo << " float\n";

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
			vector<Vector3d>bcs;
			Real TotTri = 0;
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
				Vector3d A = pos[finallist.at(j)] - pos[i];
				Vector3d B = pos[finallist.at(j+1)] - pos[i];
				Real tri = A.cross(B).norm()/2.0;
				TotTri += tri;
				pair<int, int>pAB;
				pAB = angle(A, B);
				bcs.push_back(test);
				trianglecount++;
				if (alreadypushed == 0) {
					f << test(0) << " " << test(1) << " " << test(2) << " " << endl;
					fo << test(0) << " " << test(1) << " " << test(2) << " " << endl;
				}
			}
		}
	}
	//-------------------------Post processing VTK files for visualisation--------------------------//
	f << "\n";
	f << "CELLS ";
	f << totalcaps << " " << totalpoints + 2 * totalcaps << endl;
	fo << "\n";
	fo << "CELLS ";

	fo << totalcaps << " " << totalpoints + 2 * totalcaps << endl;

	int cellnow = nonrepeat_part;
	for (int i = 0; i < totalcaps; i++) {
		int nodenow = vtknodes.at(i);
		int csize = pos_neighbor[nodenow].size();
		f << csize + 1 << " ";
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
	f.close();
	fo.close();
	clock_t stoptimePRINT = clock();
	Real sumPRINT = ((double)(stoptimePRINT - starttimePRINT) / CLOCKS_PER_SEC);
	clock_t stoptime = clock();
	Real sum = ((double)(stoptime - starttime) / CLOCKS_PER_SEC);
	cout << "ProcessorRunningTime: " << sum << endl;
}
