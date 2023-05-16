#ifndef _HELPERFUNCTIONS_H_
#define _HELPERFUNCTIONS_H_
#include "Tessellation.h"


namespace Tessellation {
	//---------Distance between two particles with given 3D co-ordinates-----------//
	Real dist(Vector3d A, Vector3d B);

	//----------------Vector search for 'a' inside Vector A-----------------//
	int vector_search(vector<int>A, int a);

	//----------------Vector search to find index of key--------//
	int find_ind(vector<int>A, int a);

	//------------Circumcenter for a triangle with vertices given-------//
	Vector3d findcc(Vector3d A, Vector3d B, Vector3d C);

	//------------Barycenter for a triangle with vertices given-------//
	Vector3d findbc(Vector3d A, Vector3d B, Vector3d C);

	//------------Angle between two vectors - two outputs - first-> smallest angle(regardless of direction) ;  second-> positive angle(0,PI)-----------//
	pair<Real, Real> angle(Vector3d A, Vector3d B);

	//----------------Map sorting with overloading------------//
	vector<int> sort(map<int, int>& M);

	vector<Real> sort(map<Real, Real>& M);

	//-------------Debugging Functions---------------//
	void print_vector(vector<int>A);

	void print_vector(vector<Real>A);

	void print_vector(vector<pair<int, int>>A);

	void test();

	
};

#endif // _HELPERFUNCTIONS_H_