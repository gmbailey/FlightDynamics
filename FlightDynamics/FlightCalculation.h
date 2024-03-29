#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <GeographicLib\TransverseMercator.hpp>
#include <GeographicLib\Geocentric.hpp>
#include <GeographicLib\LocalCartesian.hpp>
#include <iomanip>
#include <GeographicLib\LocalCartesian.hpp>
#include <fstream>
#include "Geolocation.h"
#include "LCartesian.h"

using std::ofstream;

using namespace GeographicLib;
using namespace std;

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
	double getHeading(double x1, double z1, double x2, double z2);
	double getPitch(double x1, double y1, double z1, double x2, double y2, double z2);
	double getRoll(double y1, double z1, double y2, double z2);
	double getDistanceToWaypoint(double lat1, double lon1, double lat2, double lon2, double x1, double y1, double z1, double x2, double y2, double z2);
	void getVelocity(double x1, double y1, double z1, double x2, double y2, double z2, double &vx1, double &vy1, double &vz1, double dist, double &heading, double &pitch, double &roll, double prevVx, double prevVy);
	void getNextCoords(double centerLat, double centerLong, double lat, double lon, double &x1, double &y1, double &z1); //, double &x2, double &y2, double &z2);
};
