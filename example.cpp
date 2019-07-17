#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include "Libs/polar_coordinate.h"
#include "Libs/coordinate.h"

using namespace std;

#define INACCURACY_DISTANCE_RANGE 5.

// Triangulate
ifstream fcoordinate;
ofstream fdistance;
ofstream ftruedistance;

vector<Coordinate> beacon;
vector<double> dist_to_beacon;
vector<double> dist_to_beacon_true;
Coordinate receiver;

// Test coordinate conversion
ifstream fcoordinateb;
ifstream fpolarb;
ifstream fpairb;

string buf;

int main(){
    /*  Test triangulate */
    fcoordinate.open("../input_coordinate.txt");
    fdistance.open("../output_distance_rand.txt");
    ftruedistance.open("../output_distance_true.txt");
    
    // Read coordinate from file
    Coordinate receiver;
    Coordinate temp;
    char eater[255];
    while (getline(fcoordinate, buf)){
        istringstream iss(buf);
        if(iss >> eater >> temp){
            beacon.push_back(temp);
        }
    }
    receiver = beacon.back();
    beacon.pop_back();

    // Randomize distance to emulate inaccuracy in measurement and output to file
    for (int i = 0; i < beacon.size(); i++){
        double ran = (double) rand() / RAND_MAX * INACCURACY_DISTANCE_RANGE * 2 - INACCURACY_DISTANCE_RANGE;
        double dist = Coordinate::euclideanDistance(receiver, beacon[i]);
        ftruedistance << dist << endl;
        dist_to_beacon_true.push_back(dist);
        dist += ran;
        if (dist < 0) {
            cout << "." << dist << endl;
            dist -= 2*ran;
        }
        fdistance << dist << endl;
        dist_to_beacon.push_back(dist);
    }

    // Triangulate 
    Coordinate ret_triangulate = Coordinate::triangulate(dist_to_beacon, beacon);

    // Result
    cout << "true coordinate: \n" << receiver << endl;  
    cout << "distance between true coordinate and triangulate coordinate: \n" << receiver.euclideanDistance(ret_triangulate) << endl;
    double sum_dist_diff = 0, sum_dist_diff_true = 0;
    for (int i = 0; i < beacon.size(); i++){
        sum_dist_diff += pow(ret_triangulate.euclideanDistance(beacon[i]) - dist_to_beacon[i], 2);
        sum_dist_diff_true += pow(receiver.euclideanDistance(beacon[i]) - dist_to_beacon[i], 2);
    }
    cout << "Fitness of true coordinate with given distance:\n" << sum_dist_diff_true << endl;
    cout << "Fitness of triangulate coordinate with given distance:\n" << sum_dist_diff << endl;
    cout << "**Lower is better" << endl;
    

    /* Test coordinate conversion */
    fcoordinateb.open("../convert_coordinate.txt");
    fpolarb.open("../convert_polar.txt");
    fpairb.open("../convert_ENU_to_WGS84.txt");

    vector<Coordinate> test_coordinate;
    vector<polarCoordinate> test_polar;


    Coordinate tempC;
    polarCoordinate tempP;
    Coordinate refPoint;
    Coordinate localENU;
    // Read coordinate from file
    while (getline(fcoordinateb, buf)){
        istringstream iss(buf);
        if(iss >> tempC){
            test_coordinate.push_back(tempC);
        }
    }
    while (getline(fpolarb, buf)){
        istringstream iss(buf);
        if(iss >> tempP){
            test_polar.push_back(tempP);
        }
    }
    getline(fpairb, buf);
    istringstream iss(buf);
    iss >> eater >> refPoint >> refPoint.getPolar();
    getline(fpairb, buf);
    iss.str(buf);
    iss >> eater >> localENU;

    // Convert coordinate
    cout << "\n***************" << endl;
    cout << "ECEF to WGS84 to ECEF:" << endl;
    for (int i = 0; i < test_coordinate.size(); i++){
        polarCoordinate ret = polarCoordinate::ECEFtoWGS84(test_coordinate[i]);
        cout    << fixed << setprecision(3) 
                << test_coordinate[i] 
                << " >> "
                << ret
                << " >> "
                << Coordinate::WGS84toECEF(ret)
                << endl;
    }
    cout << "\nWGS84 to ECEF to WGS84:" << endl;
    for (int i = 0; i < test_polar.size(); i++){
        Coordinate ret = Coordinate::WGS84toECEF(test_polar[i]);
        cout    << fixed << setprecision(3) 
                << test_polar[i] 
                << " >> "
                << ret 
                << " >> "
                << polarCoordinate::ECEFtoWGS84(ret)
                << endl;
    }
    cout << "\nENU to WGS84 to ENU" << endl;
    polarCoordinate ret = polarCoordinate::ENUtoWGS84(refPoint, localENU);
    cout    << fixed << setprecision(8) 
            << refPoint
            << " >> "
            << ret
            << " >> "
            << Coordinate::WGS84toENU(refPoint, ret)
            << endl; 
    
    fdistance.close();
    fcoordinate.close();
    ftruedistance.close();
    fpolarb.close();
    fcoordinateb.close();
    return 0;
}