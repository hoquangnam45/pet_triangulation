#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include "Libs/coordinate.h"
using namespace std;
ifstream fcoordinate;
ifstream fdistance;

Coordinate gateway[5];
double dist_to_gateway[5];
Coordinate receiver;

char buffer[255];
int main(){
    fcoordinate.open("../input coordinate.txt");
    fdistance.open("../output distance.txt");
    if (!fcoordinate.is_open()){
        cout << "Không mở được file input" << endl;
        return 1;
    }
    for (int i = 0; i < 6; i++){
        char eater[255];
        if (i < 5) {
            fcoordinate >> eater >> gateway[i];
            // cout << gateway[i] << endl;
        }
        else {
            fcoordinate >> eater >> receiver;
            // cout << receiver << endl;
        }
    }
    for (int i = 0; i < 5; i++){
        fdistance >> dist_to_gateway[i];
        //fdistance << fixed << setprecision(12) << dist_to_gateway[i] << endl;
    }
    vector<double> dist_to_beacon(dist_to_gateway, dist_to_gateway + 5);
    vector<Coordinate> beacon(gateway, gateway + 5);
    Coordinate::triangulate(dist_to_beacon, beacon);
    //cout << dist_to_beacon << endl;
    // for(int i = 0; i < 2; i++){
    //     cin >> gateway[i];
    // }
    // cout << gateway[0].euclideanDistance(gateway[1]) << endl;
    fdistance.close();
    fcoordinate.close();
    return 0;
}