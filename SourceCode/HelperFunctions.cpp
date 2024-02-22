#include "HelperFunctions.h"

namespace Tessellation {
	//---------Distance between two particles with given 3D co-ordinates-----------//
	Real dist(Vector3d A, Vector3d B) {
		return sqrt(pow(A(0) - B(0), 2.0) + pow(A(1) - B(1), 2.0) + pow(A(2) - B(2), 2.0));
	}

	//----------------Vector search for 'a' inside Vector A-----------------//
	int vector_search(vector<int>&A, int a) {
		int Asize = A.size();
		int counta = 0;
		for (int i = 0; i < Asize; i++) {
			if (A.at(i) == a) {
				counta++;
			}
		}
		if (counta != 0) {
			return 1;
		}
		else return 0;
	}

	//----------------Vector search to find index of key--------//
	int find_ind(vector<int>&A, int a) {
		int Asize = A.size();
		int counta = 0;
		for (int i = 0; i < Asize; i++) {
			if (A.at(i) == a) {
				return i;
			}
		}
		return 0;
	}

	//------------Circumcenter for a triangle with vertices given-------//
	Vector3d findcc(Vector3d A, Vector3d B, Vector3d C) {
		Real den = 2.0 * pow(((B - A).cross(C - A)).norm(), 2);
		Vector3d num = pow((C - A).norm(), 2.0) * ((B - A).cross(C - A)).cross(B - A) + pow((B - A).norm(), 2.0) * ((C - A).cross(B - A)).cross(C - A);
		return A + num / den;
	}

	//------------Barycenter for a triangle with vertices given-------//
	Vector3d findbc(Vector3d A, Vector3d B, Vector3d C) {
		return (A + B + C) / 3.0;
	}

	//------------Angle between two vectors - two outputs - first-> smallest angle(regardless of direction) ;  second-> positive angle(0,PI)-----------//
	pair<Real, Real> angle(Vector3d A, Vector3d B) {
		pair<Real, Real>Result;
		Real DP = A.dot(B);
		//cout << DP << "DP" << endl;
		Vector3d CP = A.cross(B);
		Real NormCP = CP.norm();
		Result.first = atan(NormCP / DP) * 180 / 3.14159265;
		Real sign = NormCP * CP.dot(CP) / abs(CP.dot(CP));
		Result.second = atan2(sign, DP) * 180 / 3.14159265;
		return Result;
	}

	//----------------Map sorting with overloading------------//

	bool cmp(pair<Real, Real>& a,
		pair<Real, Real>& b)
	{
		return a.first < b.first;
	}

	bool cmp2(pair<int, int>& a,
		pair<int, int>& b)
	{
		if (a.second == b.second)
			return 0;
		return a.second < b.second;
	}

	//Real,Real maps sorted by first key

	vector<Real> sort(map<Real, Real>& M)
	{

		vector<pair<Real, Real> > A;

		for (auto& it : M) {
			A.push_back(it);
		}


		sort(A.begin(), A.end(), cmp);

		vector<Real>Result;

		for (auto& it : A) {

			Result.push_back(it.second);
		}

		return Result;
	}

	//int,int maps sorted by second key

	vector<int> sort(map<int, int>& M)
	{

		vector<pair<int, int> > A;

		for (auto& it : M) {
			A.push_back(it);
		}


		sort(A.begin(), A.end(), cmp2);

		vector<int>Result;

		for (auto& it : A) {

			Result.push_back(it.first);
		}

		return Result;
	}

	//-------------Debugging Functions---------------//
	void print_vector(vector<int>&A) {
		int Asize = A.size();
		for (int i = 0; i < Asize; i++) {
			cout << A.at(i) << ", ";
		}
		cout << "printed vector of ints" << endl;
	}

	void print_vector(vector<Real>&A) {
		int Asize = A.size();
		for (int i = 0; i < Asize; i++) {
			cout << A.at(i) << ", ";
		}
		cout << "printed vector of Reals" << endl;
	}

	void print_vector(vector<pair<int, int>>&A) {
		int Asize = A.size();
		for (int i = 0; i < Asize; i++) {
			cout << A.at(i).first << ", " << A.at(i).second << endl;
		}
		cout << "printed vector of pairs" << endl;
	}

	void test() {
		cout << "test" << endl;
	}

	//----------Function for searching particles that are inside a cube of 2*R_c length keeping the target particle in the center-------------//
	vector<vector<int>> rccontainer(vector<Vector3d>&pos, int n_part) {
		vector<vector<int>>finalcontainerrc;
		finalcontainerrc.resize(n_part, vector<int>());
		vector<pair<Real,int>>xdir;
		for (int i = 0; i < n_part; i++) {
			xdir.push_back(make_pair(pos[i][0], i));
		}
		sort(xdir.begin(), xdir.end());
		deque<int>dqx;
		for (int i = 0; i < n_part; i++) {
			int xnow = xdir[i].second;
			Real yRr = pos[xnow][1] + R_c;
			Real yRl = pos[xnow][1] - R_c;
			Real zRr = pos[xnow][2] + R_c;
			Real zRl = pos[xnow][2] - R_c;
			Real xRr = xdir[i].first + R_c;
			Real xRl = xdir[i].first - R_c;
			int xrind = i+1;
			int xlind = i - 1;
			while (xrind<n_part) {
				if (xdir[xrind].first < xRr) {
					if (pos[xdir[xrind].second][1]<yRr && pos[xdir[xrind].second][1] > yRl && pos[xdir[xrind].second][2]<zRr && pos[xdir[xrind].second][2] > zRl) {
						finalcontainerrc[xnow].push_back(xdir[xrind].second);
					}
					xrind++;
				}
				else break;
			}
			while (xlind > -1) {
				if (xdir[xlind].first > xRl) {
					if (pos[xdir[xlind].second][1]<yRr && pos[xdir[xlind].second][1] > yRl && pos[xdir[xlind].second][2]<zRr && pos[xdir[xlind].second][2] > zRl) {
						finalcontainerrc[xnow].push_back(xdir[xlind].second);
					}
					xlind--;
				}
				else break;
			}
		}	
		return finalcontainerrc;
	}


}
