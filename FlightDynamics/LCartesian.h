#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include "Geolocation.h"
#include <iostream>
using namespace std;
class LCartesian
{
public:
	LCartesian(double lat0, double lon0, double h0 = 0, const Geolocation& earth = Geolocation::WGS84())
		: _earth(earth)
		{ Reset(lat0, lon0, h0); }

	explicit LCartesian(const Geolocation& earth = Geolocation::WGS84()) : _earth(earth)
		{ Reset(double(0), double(0), double(0)); }

	void Reset(double lat0, double lon0, double h0 = 0);

	void Forward(double lat, double lon, double h, double& x, double& y, double& z) const {
		IntForward(lat, lon, h, x, y, z, NULL);
	}

	void Forward(double lat, double lon, double h, double& x, double& y, double& z, std::vector<double>& M) const  {
		if (M.end() == M.begin() + dim2_) {
			double t[dim2_];
			IntForward(lat, lon, h, x, y, z, t);
			std::copy(t, t + dim2_, M.begin());
		}
		else
			IntForward(lat, lon, h, x, y, z, NULL);
	}

	void Reverse(double x, double y, double z, double& lat, double& lon, double& h) const {
		IntReverse(x, y, z, lat, lon, h, NULL);
	}

	void Reverse(double x, double y, double z, double& lat, double& lon, double& h, std::vector<double>& M) const {
		if (M.end() == M.begin() + dim2_) {
			double t[dim2_];
			IntReverse(x, y, z, lat, lon, h, t);
			std::copy(t, t + dim2_, M.begin());
		}
		else
			IntReverse(x, y, z, lat, lon, h, NULL);
	}

private:
	static const size_t dim_ = 3;
	static const size_t dim2_ = dim_ * dim_;
	Geolocation _earth;
	double _lat0, _lon0, _h0;
	double _x0, _y0, _z0, _r[dim2_];
	void IntForward(double lat, double lon, double h, double& x, double& y, double& z, double M[dim2_]) const;
	void IntReverse(double x, double y, double z, double& lat, double& lon, double& h, double M[dim2_]) const;
	void MatrixMultiply(double M[dim2_]) const;
};

