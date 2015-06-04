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

	cout << waypoints.size() << endl;
	for (int i = 0; i < waypoints.size(); i++) {
		cout << "x: " << waypoints[i].x << ", y: " << waypoints[i].y << ", z: " << waypoints[i].z << endl;
	}
}

void calculateFlightData(vector<cartesianCoords> waypoints, vector<flightData> &newData) {
	FlightCalculation calc;
	flightData curData;

	double newLat = 0;
	double newLong = 0;
	double centerLat = 33.5;
	double centerLong = 36.0;	
	double dt = 1.0;
	double distance = 0;

//	for (int i = 0; i < waypoints.size(); i++) {
//		distance = sqrt(waypoints[i].x - )
//		while(wayp)
//	}
	curData.x = waypoints[0].x;
	curData.y = waypoints[0].y;
	curData.z = waypoints[0].z;
	//curData.lat = calc.getLat(9900);
	//calc.updateLatLong(33.5, 36.0, -99882.50, 9900, 1200, newLong, newLat);
	double hdng = calc.getHeading(-99882.50, -99647.50, 1200, 1200);
	double pitch = calc.getPitch(-99882.50, -99647.50, 9900.00, 9900.00, 1200, 1200);
	double roll = calc.getRoll(9900.00, 9900.00, 1200, 1200);
	calc.getDistanceToWaypoint(33.58434, 34.92187, 33.58784, 35.45961);
	calc.getVelocity(-99882.50, 9900, 1200, -50000, 9900, 1200, curData.vx, curData.vy, curData.vz);
	cout << "curData : " << curData.vx << "  " << curData.vy << "  " << curData.vz << endl;
}
