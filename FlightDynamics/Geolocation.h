#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits.h>
#include <iostream>
#include <vector>

using namespace std;

class Geolocation
{
public:
	Geolocation(double radius, double flattening);
	~Geolocation();

//	void forwardCalc(double lat, double lon, double h, double &x, double &y, double &z) const;
	void Forward(double lat, double lon, double h, double& X, double& Y, double& Z, std::vector<double>& M) const {
		if (!Init())
			return;
		if (M.end() == M.begin() + dim2) {
			double t[dim2];
			IntForward(lat, lon, h, X, Y, Z, t);
			std::copy(t, t + dim2, M.begin());
		}
		else
			IntForward(lat, lon, h, X, Y, Z, NULL);
	}

	void Forward(double lat, double lon, double h, double& X, double& Y, double& Z)
		const {
		if (Init())
		IntForward(lat, lon, h, X, Y, Z, NULL);
	}

	void Reverse(double X, double Y, double Z, double& lat, double& lon, double& h)  const {
	if (Init())
		IntReverse(X, Y, Z, lat, lon, h, NULL);
	}

	void Reverse(double x, double y, double z, double &lat, double &lon, double &h, std::vector<double>& M) const {
		if (!Init())
			return;
		if (M.end() == M.begin() + dim2) {
			double t[dim2];
			IntReverse(x, y, z, lat, lon, h, t);
			std::copy(t, t + dim2, M.begin());
		}
		else
			IntReverse(x, y, z, lat, lon, h, NULL);
	}

	bool Init() const { return major_radius > 0; }

	static const Geolocation& WGS84();

private:
	friend class LCartesian;
	double lat;
	double lon;
	double x;
	double y;
	double z;
	double h;
	double flattening;
	double major_radius;
	double minor_radius;
	double latOrigin;
	double lonOrigin;
	double heightOrigin;
	double firstEccentric;
	double e2m;
	double e2a;
	double e4a;
	double maxRad;
	static const size_t dim = 3;
	static const size_t dim2 = dim * dim;
	static void Rotation(double sphi, double cphi, double slam, double clam, double M[dim2]);
	static void Rotate(double M[dim2], double x, double y, double z, double& X, double& Y, double& Z) {
	// Perform [X,Y,Z]^t = M.[x,y,z]^t
	// (typically local cartesian to geocentric)
	X = M[0] * x + M[1] * y + M[2] * z;
	Y = M[3] * x + M[4] * y + M[5] * z;
	Z = M[6] * x + M[7] * y + M[8] * z;
	}
	static void Unrotate(double M[dim2], double X, double Y, double Z,
		double& x, double& y, double& z) {
		// Perform [x,y,z]^t = M^t.[X,Y,Z]^t
		// (typically geocentric to local cartesian)
		x = M[0] * X + M[3] * Y + M[6] * Z;
		y = M[1] * X + M[4] * Y + M[7] * Z;
		z = M[2] * X + M[5] * Y + M[8] * Z;
	}
	void IntForward(double lat, double lon, double h, double& X, double& Y, double& Z,
		double M[dim2]) const;
	void IntReverse(double X, double Y, double Z, double& lat, double& lon, double& h,
		double M[dim2]) const;
};

