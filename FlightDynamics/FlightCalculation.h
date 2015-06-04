#pragma once

#define _USE_MATH_DEFINES
//#include <math.h>
#include <cmath>
#include <iostream>
#include <GeographicLib/TransverseMercator.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#include <boost\geometry.hpp>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <boost\format.hpp>
#include <iomanip>
#include <GeographicLib/LocalCartesian.hpp>
#include <fstream>
using std::ofstream;

using namespace GeographicLib;
using namespace std;
namespace bg = boost::geometry;

class FlightCalculation
{
public:
	FlightCalculation();
	~FlightCalculation();

	double getvx(double x2, double x1, double dt);
	double getvy(double y2, double y1, double dt);
	double getvz(double z2, double z1, double dt);
	double getLat(double y);
	void updateLatLong(double centerLat, double centerLong, double x, double y, double z, double &newLat, double &newLong);
	double getHeading(double x1, double x2, double z1, double z2);
	double getPitch(double x1, double x2, double y1, double y2, double z1, double z2);
	double getRoll(double y1, double y2, double z1, double z2);
	double getDistanceToWaypoint(double lat1, double lon1, double lat2, double lon2);
	void getVelocity(double x1, double x2, double y1, double y2, double z1, double z2, double &vx1, double &vy1, double &vz1);
};
