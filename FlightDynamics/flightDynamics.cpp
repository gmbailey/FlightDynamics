#include <iostream>
#include <boost/format.hpp>
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

//	cout << waypoints.size() << endl;
//	for (int i = 0; i < waypoints.size(); i++) {
//		cout << "x: " << waypoints[i].x << ", y: " << waypoints[i].y << ", z: " << waypoints[i].z << endl;
//	}
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

	xPt = waypoints[2].x;
	yPt = waypoints[2].y;
	zPt = waypoints[2].z;

	calc.getVelocity(xPt, yPt, zPt, waypoints[3].x, waypoints[3].y, waypoints[3].z, curData.vx, curData.vy, curData.vz, distance, curData.heading, curData.pitch, curData.roll);
	
	/*
	curData.x = waypoints[0].x;
	curData.y = waypoints[0].y;
	curData.z = waypoints[0].z;
	calc.updateLatLong(centerLat, centerLong, curData.x, curData.y, curData.z, curData.lat, curData.lon);
	calc.updateLatLong(centerLat, centerLong, waypoints[1].x, waypoints[1].y, waypoints[1].z, nextLat, nextLon);
	distance = calc.getDistanceToWaypoint(curData.lat, curData.lon, nextLat, nextLon, curData.x, curData.y, curData.z, waypoints[1].x, waypoints[1].y, waypoints[1].z);
	calc.getVelocity(waypoints[0].x, waypoints[0].y, waypoints[0].z, waypoints[1].x, waypoints[1].y, waypoints[1].z, curData.vx, curData.vy, curData.vz, distance, curData.heading, curData.pitch, curData.roll);
	curData.heading = calc.getHeading(curData.x, curData.z, waypoints[1].x, waypoints[1].z);
	curData.pitch = calc.getPitch(curData.x, curData.y, curData.z, waypoints[1].x, waypoints[1].y, waypoints[1].z);
	curData.roll = calc.getRoll(curData.y, curData.z, waypoints[1].y, waypoints[1].z);

	xPt = curData.x;
	yPt = curData.y;
	zPt = curData.z;
	newData.push_back(curData);

	for (int i = 0; i < waypoints.size() - 1; i++) {
		calc.updateLatLong(centerLat, centerLong, curData.x, curData.y, curData.z, curData.lat, curData.lon);
		calc.updateLatLong(centerLat, centerLong, waypoints[i+1].x, waypoints[i + 1].y, waypoints[i + 1].z, nextLat, nextLon);
		distance = calc.getDistanceToWaypoint(curData.lat, curData.lon, nextLat, nextLon, curData.x, curData.y, curData.z, waypoints[i + 1].x, waypoints[i + 1].y, waypoints[i + 1].z);
		while (distance >= 235) {
			calc.getVelocity(curData.x, curData.y, curData.z, waypoints[i + 1].x, waypoints[i + 1].y, waypoints[i + 1].z, curData.vx, curData.vy, curData.vz, distance, curData.heading, curData.pitch, curData.roll);
			xPt += curData.vx;
			yPt += curData.vy;
			zPt += curData.vz;
			
			curData.x = xPt;
			curData.y = yPt;
			curData.z = zPt;
			calc.updateLatLong(centerLat, centerLong, curData.x, curData.y, curData.z, curData.lat, curData.lon);
			
//			curData.heading = calc.getHeading(xPt, zPt, waypoints[i + 1].x, waypoints[i + 1].z);
//			curData.pitch = calc.getPitch(xPt, yPt, zPt, waypoints[i + 1].x, waypoints[i + 1].y, waypoints[i + 1].z);
//			curData.roll = calc.getRoll(yPt, zPt, waypoints[i + 1].y, waypoints[i + 1].z);
			distance = calc.getDistanceToWaypoint(curData.lat, curData.lon, nextLat, nextLon, curData.x, curData.y, curData.z, waypoints[i + 1].x, waypoints[i + 1].y, waypoints[i + 1].z);
//			cout << distance << endl;
			newData.push_back(curData);	
//			distance -= 235;
		}
	}

		ofstream outdata;
		outdata.open("sample.txt");

	cout << newData.size() << "total size " << endl;
	outdata << left << setw(15) << "x(m)" << setw(15) << "y(m)" << setw(15) << "x(m)" << setw(15) << "lat(degs)" << setw(15) << "lon(degs)" << setw(15) << "heading(rads)" << setw(15)
		<< "pitch(rads)" << setw(15) << "roll(rads)" << endl;
	for (int i = 0; i < newData.size(); i++) {
		outdata << setprecision(5) << fixed << left << setw(15) << newData[i].x << setw(15) << newData[i].y << setw(15) << newData[i].z << setw(15) << newData[i].lat << setw(15) << newData[i].lon << setw(15) << newData[i].heading <<
			setw(15) << newData[i].pitch << setw(15) << newData[i].roll << endl;
//		outdata  << setprecision(5) << fixed << i+1 << " : " << newData[i].x << "  " << newData[i].y << "  " << newData[i].z << "  " << newData[i].lat << "  " << newData[i].lon << "  " 
//			<< "  " << newData[i].vx << "  " << newData[i].vy << "  " << newData[i].vz << "   "  << newData[i].heading << "  "  << newData[i].pitch << "  " << newData[i].roll << endl;
	}


*/


//	curData.x = waypoints[0].x;
//	curData.y = waypoints[0].y;
//	curData.z = waypoints[0].z;
	//curData.lat = calc.getLat(9900);
	//calc.updateLatLong(33.5, 36.0, -99882.50, 9900, 1200, newLong, newLat);
//	double hdng = calc.getHeading(-99882.50, -99647.50, 1200, 1200);
//	double pitch = calc.getPitch(-99882.50, -99647.50, 9900.00, 9900.00, 1200, 1200);
//	double roll = calc.getRoll(9900.00, 9900.00, 1200, 1200);
//	calc.getDistanceToWaypoint(33.58434, 34.92187, 33.58784, 35.45961);
//	calc.getVelocity(-99882.50, 9900, 1200, -50000, 9900, 1200, curData.vx, curData.vy, curData.vz);
//	cout << "curData : " << curData.vx << "  " << curData.vy << "  " << curData.vz << endl;
}
