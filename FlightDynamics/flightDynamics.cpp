//#include <boost/format.hpp>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include "FlightCalculation.h"

using namespace std;

struct cartesianCoords {
	double x = 0;
	double y = 0;
	double z = 0;
};

struct flightData {
	double time = 0;
	double x = 0;
	double y = 0;
	double z = 0;
	double vx = 0;
	double vy = 0;
	double vz = 0;
	double lat = 0;
	double lon = 0;
	double alt = 0;
	double spd = 0;
	double heading = 0;
	double pitch = 0;
	double roll = 0;
};

void loadWaypoints(cartesianCoords&, vector<cartesianCoords>&);
void calculateFlightData(vector<cartesianCoords>, vector<flightData>&);

int main(){
	cartesianCoords coord;
	vector<cartesianCoords> waypoints;
	vector<flightData> newData;
	loadWaypoints(coord, waypoints);
/*
	double w, x, y, z, x2, y2, z2;
	w = 1;
	x = -99882.50;
	y = 9980;
	z = 1200;

	x2 = -50000;
	y2 = 9900;
	z2 = 1200;

	double dX = x - -50000;
	double dY = y - 9900;
	double dZ = z - 1200;

	double yaw = 450 - atan2(dZ, dX);
	double heading = fmod(yaw, 360.0) * (M_PI / 180.0);
	heading = atan(sqrt(pow(y, 2) + pow(x, 2)) / z);
	cout << "HEADING?! : " << heading << endl;

	//double yaw = atan2(2.0*(x*y + w*z), w*w * 2 + x*x * 2 - y*y * 2 - z*z * 2);
//	double pitch = -asin(2 * w*y - 2 * x*z) * 180.0 / M_PI;
//	double pitch = atan2(sqrt(dZ * dZ + dX * dX), dY) + M_PI;
//	pitch = fmod(pitch, 360.0) * (M_PI / 180.0);
	double roll = atan2(dY, dZ);

	double arg2 = std::sqrt(dX*dX + dY*dY);
	double pitch = std::atan2(dZ, arg2);

//	Quaterniond q(x, y, z, w);
//	Quaterniond q2(-50000, 9900, 1200, w);
	double test = q.x() * q.y() + q.z() * q.w();
//	q.norm();
//	q2.norm();
	if (test > 0.499) { // singularity at north pole
		heading = 2 * atan2(q.x(), q.w());
		double attitude = M_PI / 2;
		double bank = 0;
		cout << bank << "   " << attitude << endl;
	}
	if (test < -0.499) { // singularity at south pole
		heading = -2 * atan2(q.x(), q.w());
		double attitude = -M_PI / 2;
		double bank = 0;
		cout << bank << "   " << attitude << endl;
	}
	cout << heading << endl;

	double sqx = q.x() * q.x();
	double sqy = q.y() * q.y();
	double sqz = q.z() * q.z();
	double bank = atan2(2 * q.x() * q.w() - 2 * q.y()* q.z() , 1 - 2 * sqx - 2 * sqz);
	cout << bank << " bank" << endl;
	bank = atan2(2 * q.x() * q.w() - 2 * q.y() * q.z(), 1 - 2 * q2.x() - 2 * q2.z());

	roll = atan2(2 * q.x() * q.w() - 2 * q.y() * q.z(), 1 - 2 * pow(q.x(), 2) - 2 * pow(q.z(), 2));			// heading : atan2(2 * q.y() * q.w() - 2 * q.x() * q.z(), 1 - 2 * q2.y() - 2 * q2.z());
	cout << roll;

	Vector3d p1(x, y, z);
	Vector3d p2(-50000, 9900, 1200);
	Vector3d angles(0,0,0);
	FlightCalculation calc;
	angles = calc.AngleTo(p1, p2);
	cout << angles << endl;
	
	Vector3d m1(x, y, z);
	Vector3d m2(-50000, 9900, 1200);
	Vector3d m3;
	m3.UnitZ();
	
	Quaterniond q3(w, x, y, z);
	Quaterniond q4(w, -50000, 9900, 1200);
//	Quaterniond q5;
//	q5 = q3.slerp(1, q4);
//	cout << "slerp: " << q5.y() << endl;
//	Quaterniond qa, qb, qres;
//	qa = m1;
//	qb = m2;
	double t = 1;
	Quaterniond qa;
	qa = q3.slerp(t, q4);
	cout << qa.w() << "  " << qa.x() << "  " << qa.y() << "   " << qa.z() << endl;
	qa.dot(q3);
	cout << qa.w() << "  " << qa.x() << "  " << qa.y() << "   " << qa.z() << endl;
	double roll2 = atan2(2 * q3.x() * q3.w() - 2 * q3.y() * q3.z(), 1 - 2 * q3.x() - 2 * q3.z());
	cout << roll2 << endl;
	
	Vector3d relPos = p1 - p2;
	relPos.normalize();
	roll2 = asin(relPos.z());
	double yaw2 = 450 -(acos(relPos.x() / cos(roll2)));
	if (relPos.y() < 0)
		yaw2 = M_PI_2 - yaw2;

	yaw2 = fmod(yaw2, 360.0) * (M_PI / 180.0);
	double pitch2 = atan2(relPos.z(), sqrt(relPos.x() * relPos.x() + relPos.y() * relPos.y()));

	cout << "roll 2 = " << roll2 << "    yaw2 = " << yaw2 << "   pitch: " << pitch2 << endl;
	yaw2 = atan2(relPos.y(), relPos.x());
	cout << "yaw2 again: " << yaw2 << endl;

	double local_x = -99882.50;
	double local_y = 9980;
	double local_z = 1200;
	double transx = -50000 - x;
	double transy = 9900 - y;
	double transz = 1200 - z;

	double heading1 = tan(34.92187 / 33.58434);
	cout << "new hdg " << heading1 << endl;


//	Vector3d r1;
//	r1.eulerAngles(yaw2, pitch2, roll2);
	double sina = sin(yaw2); double sinb = sin(pitch2); double siny = sin(roll2);
	double cosa = cos(yaw2); double cosb = cos(pitch2); double cosy = cos(roll2);
	local_x = transx*(cosa*cosy - sina*cosb*siny) + transy*(-sina*cosy - cosa*cosb*siny) + transz*(sinb*siny);
	local_y = transx*(cosa*siny + sina*cosb*cosy) + transy*(-sina*siny + cosa*cosb*cosy) + transz*(-sinb*cosy);
	local_z = transx*(sina*sinb) + transy*(cosa*sinb) + transz*(cosb);
	
	cout << "localx : " << local_x << "localy : " << local_y << "localz : " << local_z << endl;

//	qres = qa.slerp(t, qb);
//	cout << qres.x() << endl;
	
//	q5.FromTwoVectors(m1, m2);

//	gte::Quaternion<double> qq(x, y, z, w);
//	gte::EulerAngles<double> eularAng;
//	gte::Rotation<4, double> rotate1(qq);
//	EulerAngles<double> eTest;
//	Rotation<4, gte::Quaternion<double>> rTest(gte::Quaternion<double> qq);
	
	
//	cout << "eular angle   "  << eularAng.angle << endl;
/*	Matrix3d mT;
	
	double cp, sp, cr, sr, cy, sy;

	cp = cos(pitch2); sp = sin(pitch2);
	cr = cos(roll2);  sr = sin(roll2);
	cy = cos(yaw2);   sy = sin(yaw2);

	mT(1, 1) = cp*cy;
	mT(1, 2) = cp*sy;
	mT(1, 3) = -sp;

	mT(2, 1) = sr*sp*cy - cr*sy;
	mT(2, 2) = sr*sp*sy + cr*cy;
	mT(2, 3) = sr*cp;

	mT(3, 1) = cr*sp*cy + sr*sy;
	mT(3, 2) = cr*sp*sy - sr*cy;
	mT(3, 3) = cr*cp;

//	cout << mT.rows() << endl;

	Quaterniond qHead(0, 1, 1, yaw2);
	Quaterniond qPitch(1, 0, 0, pitch2);
	Quaterniond qRoll(0, 0, 0, roll2);
	Quaterniond qH, qP, qR;

	qH = q3 * qHead; 
	qP = q3 * qPitch;
	qR = q3 * qRoll;

	Vector3d v1;
	
	v1 = p2 - p1;
	Vector3d vp(v1.x(), v1.y(), 0);
	//double newYaw = acos((v1*vp)  v1.cwiseAbs() *vp.cwiseAbs())// / v1.cwiseAbs() *vp.cwiseAbs());
	Vector3d v2, v3, v4, v5;
//	v2 = v1*vp;
	double roll3 = atan2((x2, y), (x2, x));
	double argP = atan2(x2, z);
	double argP2 = atan2(z2, z);
	double pitch3 = atan2(argP, argP2);
	cout << "roll 3  " << roll3 << "    " << pitch3 << endl;
//	cout << v4.array().acos() << endl;
*/

	calculateFlightData(waypoints, newData);

	  	return 0;
}

void loadWaypoints(cartesianCoords &coord, vector<cartesianCoords> &waypoints) {
	ifstream indata("waypoints.txt");
	string line;
	if(!indata){
		cerr << "File could not be opened!" << endl;
		exit(1);
	}

	while (!indata.eof()) {
		getline(indata, line, ',');
		coord.x = stod(line);
		getline(indata, line, ',');
		coord.y = stod(line);
		getline(indata, line);
		coord.z = stoi(line);
		waypoints.push_back(coord);
	}
	indata.close();

	cout << waypoints.size() << endl;
	for (int i = 0; i < waypoints.size(); i++) {
		cout << "x: " << waypoints[i].x << ", y: " << waypoints[i].y << ", z: " << waypoints[i].z << endl;
	}
}

void calculateFlightData(vector<cartesianCoords> waypoints, vector<flightData> &newData) {
	FlightCalculation calc;
	flightData curData;

	double centerLat = 33.5;
	double centerLong = 36.0;	
	double dt = 1.0;
	double distance = 0;
	double nextLat, nextLon;
	double xPt, yPt, zPt;
	double speed, altitude, time = 0;

	/*
	//Quaterniond qPos1, qPos2, qVel, qHdg, qPitch, qRoll;
	double rotationRate = .047;
	double angleUp = rotationRate * -1;;   //angles for the rotations
	double angleDown = rotationRate;
	double angle2up = rotationRate * -1;
	double angle2down = rotationRate;
	*/

	double prevVx;
	double prevVy = 0;
	prevVx = 235;

	speed = 235;
	altitude = 3937.00781;
	time = 1;
	curData.alt = altitude;
	curData.spd = speed;
	curData.time = time;
	curData.x = waypoints[0].x;
	curData.y = waypoints[0].y;
	curData.z = waypoints[0].z;
	calc.updateLatLong(centerLat, centerLong, curData.x, curData.y, curData.z, curData.lat, curData.lon);
	calc.updateLatLong(centerLat, centerLong, waypoints[1].x, waypoints[1].y, waypoints[1].z, nextLat, nextLon);
	distance = calc.getDistanceToWaypoint(curData.lat, curData.lon, nextLat, nextLon, curData.x, curData.y, curData.z, waypoints[1].x, waypoints[1].y, waypoints[1].z);
	calc.getVelocity(waypoints[0].x, waypoints[0].y, waypoints[0].z, waypoints[1].x, waypoints[1].y, waypoints[1].z, curData.vx, curData.vy, curData.vz, distance, curData.heading, curData.pitch, curData.roll, prevVx, prevVy);

	newData.push_back(curData);

	xPt = curData.x;
	yPt = curData.y;
	zPt = curData.z;


	for (int i = 0; i < waypoints.size() - 1; i++) {
		calc.updateLatLong(centerLat, centerLong, curData.x, curData.y, curData.z, curData.lat, curData.lon);
		calc.updateLatLong(centerLat, centerLong, waypoints[i+1].x, waypoints[i + 1].y, waypoints[i + 1].z, nextLat, nextLon);
		distance = calc.getDistanceToWaypoint(curData.lat, curData.lon, nextLat, nextLon, curData.x, curData.y, curData.z, waypoints[i + 1].x, waypoints[i + 1].y, waypoints[i + 1].z);
		while (distance >= 235) {
			time++;
			curData.time = time;
			curData.spd = speed;
			curData.alt = altitude;
			prevVx = curData.vx;
			prevVy = curData.vy;
			calc.getVelocity(curData.x, curData.y, curData.z, waypoints[i + 1].x, waypoints[i + 1].y, waypoints[i + 1].z, curData.vx, curData.vy, curData.vz, distance, curData.heading, curData.pitch, curData.roll, prevVx, prevVy);
			xPt += curData.vx;
			yPt += curData.vy;
			zPt += curData.vz;

			curData.x = xPt;
			curData.y = yPt;
			curData.z = zPt;
			calc.updateLatLong(centerLat, centerLong, curData.x, curData.y, curData.z, curData.lat, curData.lon);
			
			
			distance = calc.getDistanceToWaypoint(curData.lat, curData.lon, nextLat, nextLon, curData.x, curData.y, curData.z, waypoints[i + 1].x, waypoints[i + 1].y, waypoints[i + 1].z);

			newData.push_back(curData);	
		}
	}

		ofstream outdata;
		outdata.open("result.txt");

	cout << newData.size() << "total size " << endl;
	outdata << left << setw(15) << "time(s)" << setw(15) << "x(m)" << setw(15) << "y(m)" << setw(15) << "x(m)" << setw(15) << "lat(degs)" << setw(15) << "lon(degs)" << setw(15) << "vx(m/s)" << setw(15)
			<< "vy(m/s)" << setw(15) << "vz(m/s)" << setw(15) << "alt_msl(feet)" << setw(15) << "speed(m/s)" << setw(15) << "heading(rads)" << setw(15) << "pitch(rads)" << setw(15) << "roll(rads)" << endl;
	for (int i = 0; i < newData.size(); i++) {
		outdata << setprecision(5) << fixed << left << setw(15) << newData[i].time << setw(15) << newData[i].x << setw(15) << newData[i].y << setw(15) << newData[i].z << setw(15) << newData[i].lat << setw(15) << newData[i].lon
			<< setw(15) << newData[i].vx << setw(15) << newData[i].vy << setw(15) << newData[i].vz << setw(15) << newData[i].alt << setw(15) << newData[i].spd << setw(15) << newData[i].heading << setw(15) << newData[i].pitch << setw(15) << newData[i].roll  << endl;
	}
}
