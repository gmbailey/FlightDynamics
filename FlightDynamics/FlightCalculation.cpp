#include "FlightCalculation.h"

FlightCalculation::FlightCalculation(){
}

FlightCalculation::~FlightCalculation(){
}

double FlightCalculation::getvx(double x2, double x1, double dt){
	double vx = 0;
	vx = (x2 - x1) / dt;
	return vx;
}

double FlightCalculation::getvy(double y2, double y1, double dt){
	double vy = 0;
	vy = (y2 - y1) / dt;
	return vy;
}

double FlightCalculation::getvz(double z2, double z1, double dt){
	double vz = 0;
	vz = (z2 - z1) / dt;
	return vz;
}

double FlightCalculation::getLat(double y){
	double deg2rad = M_PI / 180.0;
	double RAD2Deg = 180.0 / M_PI;
	double lat = 0;
//	lat = 180.0 / M_PI * (2.0 * atan(exp(y * (M_PI / 180))) - M_PI / 2);
	
	
	cout << setprecision(5) << fixed;
	double lat0 = 33.5 - .000235; //33.4997;
	double lonCent = -60.88;		//higher decreases lower increases
	double long0 = 36 - .00246;//35.99755;
	double x = -99882.50;
	double z = 1200.00;
	double lat1, lon, h, lat2;

	TransverseMercator proj(Constants::WGS84_a(), Constants::WGS84_f(),
		Constants::UTM_k0());

	proj.Reverse(lonCent, x, y, lat1, lon);
//	cout << "real Lat1 - " << lat1 + 33.5 << endl;
	lon = lon * (M_PI / 180.0) + lat0;
//	cout << "real long - " << lon << endl;
	double _f = 30000.001;

	TransverseMercator proj1(100000.0, Constants::WGS84_f(),						//6776767.0		//33.58434			33.58892
		35);
	lonCent = 36;//lonCent = -60.9270;
	for (int i = 0; i < 4; i++) {
		proj1.Reverse(lonCent, x, y, lat1, lon);
	//	cout << "Lat - " << lat1 + 33.5 << endl;
	//	lon = lon * (M_PI / 180.0) + 36;
		cout << i+1 << " : lat/long - " << lat1 + 33.50 << ", " << lon << endl;
		x += 235;
	}

	ofstream outdata;
	outdata.open("sample.txt");

	Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
	LocalCartesian proj5(lat0, long0, 0, earth);
	x = -99882.50;
	for (int i = 1; i < 25; i++) {
		proj5.Reverse(x, 9900.00, 1200, lat1, lon, h);
		//	lon = lon * (M_PI / 180.0) + 36;
		cout << i  << " : " << lat1 << " " << lon << " " << x << "\n";
		outdata << setprecision(5) << fixed << i << " | " <<  "x: " << x << " y: " << y << " z: " << z << " lat: " << lat1 << " long : " << lon << endl;
		x += 235;
	}
	outdata.close();

//	TransverseMercator proj3(Constants::WGS84_a(), Constants::WGS84_f(),							//33.58434			33.58892
//		.8996);
//	lonCent = 35.8;
//	proj3.Reverse(lonCent, x, y, lat1, lon);
//	lon = lon * (M_PI / 180.0) + 36;
//	cout << lat1 + 33.5 << " " << lon << "\n";

/*	Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
	LocalCartesian proj1(lat0, lon0, 0, earth);
	proj1.Reverse(x, y, z, lat2, lon, h);
	cout << "new Lat2 - " << lat2 << endl;
	cout << "new lon - " << lon << endl;

	double R_MAJOR = 6378137.0;
	double R_MINOR = 6356752.3142;
	double RATIO = R_MINOR / R_MAJOR;
	double ECCENT = 0.0818;
	double COM = 0.5 * ECCENT;
	double DEG2RAD = M_PI / 180.0;
	double RAD2Deg = 180.0 / M_PI;
	double PI_2 = M_PI / 2.0;


	double ts = exp(-y / R_MAJOR);
	double phi = PI_2 - 2 * atan(ts);
	double dphi = 1.0;
	int i = 0;
	while ((abs(dphi) > 0.000000001) && (i < 15)) {
		double con = ECCENT * sin(phi);
		dphi = PI_2 - 2 * atan(ts * pow((1.0 - con) / (1.0 + con), COM)) - phi;
		phi += dphi;
		i++;
	}
	cout << "new lat0: " << phi << endl;
	phi = phi * RAD2Deg;
	cout << "new lat: " << phi << endl;
*/

	x = -50000;
	double nextLat, nextLon;
	proj5.Reverse(x, 9900.00, 1200, nextLat, nextLon, h);
	cout << nextLat << " " << nextLon << " " << x << "\n";

	double latOrig = deg2rad * lat1;
	double longOrig = deg2rad * lon;
	//	cout << "orig lat2? :" << latOrig << endl;
	//	cout << "orig lon2? :" << longOrig << endl;

	nextLat = deg2rad * nextLat;
	nextLon = deg2rad * nextLon;

	//	cout << "next lat2? :" << nextLat << endl;
	//	cout << "next lon2? :" << nextLong << endl;

	double R = 6378137;
	double delLat = (latOrig - nextLat);
	double delLon = (longOrig - nextLon);
	//	cout << "del lat2? :" << delLat << endl;
	//	cout << "del lon2? :" << delLon << endl;

	double a = sin(delLat / 2) * sin(delLat / 2) + cos(latOrig) * cos(nextLat) * sin(delLon / 2) * sin(delLon / 2);
	//	cout << "a : " << a << endl;
	double c = 2 * atan2(sqrt(a), sqrt(1 - a));
	//	cout << "c : " << c << endl;
	double distance = R * c;

	return lat;
}

void FlightCalculation::updateLatLong(double centerLat, double centerLong, double x, double y, double z, double &newLat, double &newLong){
	double deg2rad = M_PI / 180.0;
	double RAD2Deg = 180.0 / M_PI;

	cout << setprecision(5) << fixed;
//	double lat0 = 33.5 - .000235; //33.4997;
//	double lonCent = -60.88;		//higher decreases lower increases
//	double long0 = 36 - .00246;//35.99755;
//	double x = -99882.50;
//	double z = 1200.00;
	double h;

	Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
	LocalCartesian projection(centerLat - .000235, centerLong - .00246, 0, earth);
	projection.Reverse(x, y, z, newLat, newLong, h);

/*	bg::model::point<long double, 3, bg::cs::cartesian> p1(x, y, z);
//	cout << "initial P1: " << bg::dsv(p1) << endl;
	bg::model::point<long double, 3, bg::cs::spherical<bg::radian>> p2;
//	cout << "initial P2: " << bg::dsv(p2) << endl;
	bg::transform(p1, p2);
//	cout << "new P1: " << bg::dsv(p1) << endl;
//	cout << "new P2: " << bg::dsv(p2) << endl;

	newLat = acos(p2.get<1>() / 100) + centerLat;
	newLong = atan(p2.get<0>() / z) + centerLong;
//	cout << "new lat1? :" << newLat << endl;
//	cout << "new lon1? :" << newLong << endl;
//	cout << "new lat? :" << lat1 << endl;
//	cout << "new lon? :" << lon1 << endl;

	bg::model::point<long double, 3, bg::cs::cartesian> p3(-50000, 9900, 1200);
	bg::model::point<long double, 3, bg::cs::spherical<bg::degree>> p4;
	bg::transform(p3, p4);

	double lat2 = acos(p4.get<1>() / 6371000) + centerLat;
	double long2 = atan(p4.get<0>() / z) + centerLong;

//	cout << "new lat2? :" << lat2 << endl;
//	cout << "new lon2? :" << long2 << endl;

	double latOrig = deg2rad * newLat;
	double longOrig = deg2rad * newLong;
//	cout << "orig lat2? :" << latOrig << endl;
//	cout << "orig lon2? :" << longOrig << endl;

	double nextLat = deg2rad * lat2;
	double nextLong = deg2rad * long2;

//	cout << "next lat2? :" << nextLat << endl;
//	cout << "next lon2? :" << nextLong << endl;
	
	double R = 6371000;
	double delLat = (latOrig - nextLat);
	double delLon = (longOrig - nextLong);
//	cout << "del lat2? :" << delLat << endl;
//	cout << "del lon2? :" << delLon << endl;

	double a = sin(delLat / 2) * sin(delLat / 2) + cos(latOrig) * cos(nextLat) * sin(delLon / 2) * sin(delLon / 2);
//	cout << "a : " << a << endl;
	double c = 2 * atan2(sqrt(a), sqrt(1-a));
//	cout << "c : " << c << endl;
	double distance = R * c;
//	cout << distance << endl;

*/

}

double FlightCalculation::getHeading(double x1, double z1, double x2, double z2){
	double dx = x2 - x1;
//	double dy = y2 - y1;
	double deg2rad = M_PI / 180.0;
	double heading;

//	double arg1 = 450 - atan2(dy, dx);
//	heading = fmod(arg1, 360.0) * (M_PI / 180.0);
//	double dx = x2 - x1;
	double dz = z2 - z1;

	heading = M_PI / 180 * (atan2(dz, dx) + 90);

	return heading;
}

double FlightCalculation::getPitch(double x1, double y1, double z1, double x2, double y2, double z2){
	double dx = x2 - x1;
	double dy = y2 - y1;
	double dz = z2 - z1;
//	double pitch = M_PI / 180 * (-atan2(dy, sqrt(dx * dx + dz * dz)));

//	double arg2 = sqrt(dx*dx + dy*dy);
	double pitch = atan2(dz, sqrt(dx*dx + dy*dy)) * (M_PI / 180.0);

	return pitch;
}

double FlightCalculation::getRoll(double y1, double z1, double y2, double z2){
	double dy = y2 - y1;
	double dz = z2 - z1;
	double roll = M_PI / 180 * (atan2(dy, dz));
	return roll;
}

double FlightCalculation::getDistanceToWaypoint(double lat1, double lon1, double lat2, double lon2, double x1, double y1, double z1, double x2, double y2, double z2){
	double distance;
	double RAD2Deg = 180.0 / M_PI;
	double deg2rad = M_PI / 180.0;
//	distance = acos(sin(lat1)*sin(lat2) + cos(lat1) + cos(lat2) + cos(lon1 - lon2));
//	cout << "distance: " << distance << endl;

	double dLat, dLong, dLat1, dLong1, dLat2, dLong2, a, c;
	dLat1 = lat1*(M_PI / 180);
	dLong1 = lon1 * (M_PI / 180);	
	dLat2 = lat2*(M_PI / 180);
	dLong2 = lon2 * (M_PI / 180);
	
	dLong = dLong1 - dLong2;
	dLat = dLat1 - dLat2;

	a = pow(sin(dLat / 2.0), 2.0) + cos(dLat1) * cos(dLat2) * pow(sin(dLong / 2), 2);
	c = 2 * atan2(sqrt(a), sqrt(1.0 - a));
	distance = 6378137 * c;
//	cout << "distance : " << distance << endl;

//	bearing = tc1 * 180 / M_PI;
//	cout << "tc1 : " << tc1 << " bearing : " << bearing << endl;

	double dt = 1;
	double dx = x2 - x1;
	double dy = y2 - y1;
	double dz = z2 - z1;
	double dist = sqrt(dx*dx + dy*dy + dz*dz);
	if (dist <= 0)
		cout << dist << " less than 0 " << endl;

	double tc1;
	if (sin(lon2 - lon1) < 0)
		tc1 = acos((sin(lat2) - sin(lat1) * cos(dist)) / (sin(dist) * cos(lat1)));
	else
		tc1 = 2 * M_PI - acos((sin(lat2) - sin(lat1) * cos(dist)) / (sin(dist) * cos(lat1)));

//	cout << tc1 << endl;

	return dist;

	/*
	IF sin(lon2-lon1)<0
	tc1=acos((sin(lat2)-sin(lat1)*cos(d))/(sin(d)*cos(lat1)))
	ELSE
	tc1=2*pi-acos((sin(lat2)-sin(lat1)*cos(d))/(sin(d)*cos(lat1)))
	ENDIF
	*/
}

void FlightCalculation::getVelocity(double x1, double y1, double z1, double x2, double y2, double z2, double &vx1, double &vy1, double &vz1, double dist, double &heading, double &pitch, double &roll){
	double vx, vy, vz, dt, dx, dy, dz, ux, uy, uz;
	double RAD2Deg = 180.0 / M_PI;
	double deg2rad = M_PI / 180.0;
	dt = 1;
	dx = (x2 - x1) / dt;
	dy = (y2 - y1) / dt;
	dz = (z2 - z1) / dt;

	ux = dx / dist;
	uy = dy / dist;
	uz = dz / dist;

	vx1 = 235 * ux;
	vy1 = 235 * uy;
	vz1 = 235 * uz;
	

	double arg1 = 450 - atan2(dy, dx);
//	heading = fmod(arg1, 360.0) * deg2rad;
//	cout << "heading : " << hdg << endl;
	heading = M_PI / 180 * (atan2(dy, dx) + 90);
	double arg2 = sqrt(dx*dx + dy*dy);
	pitch = atan2(dz, arg2);
	roll = M_PI / 180 * (atan2(dy, dz));

	cout << heading << "  " << pitch << "  " << roll << endl;


/*	while (dist > 235) {
		dist -= 235;
		dx = (x2 - x1) / dt;
		dy = (y2 - y1) / dt;
		dz = (z2 - z1) / dt;

		ux = dx / dist;
		uy = dy / dist;
		uz = dz / dist;

		vx = 235 * ux;
		vy = 235 * uy;
		vz = 235 * uz;

		arg1 = 450 - atan2(dy, dx);
		hdg = fmod(arg1, 360.0) * deg2rad;
		arg2 = sqrt(dx*dx + dy*dy);
		pitch = atan2(dz, arg2);
		roll = M_PI / 180 * (atan2(dy, dz));

		x1 += vx;
		y1 += vy;
		z1 += vz;

		cout << x1 << "  " << y1 << "  " << z1 << "      " << hdg << "      " << pitch << "   " << roll << endl;
	}
	*/
}

void FlightCalculation::qRotation(){
	double w, x, y, z;
	gte::Quaternion<double> q(x, y, z, w);
	
	double yaw = asin(-2 * (q[x] * q[z] - q[w] * q[y]));


	

}
