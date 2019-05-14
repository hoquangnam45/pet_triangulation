#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include "Libs/coordinate.h"
using namespace std;
ifstream fcoordinate;
ofstream fdistance;

Coordinate gateway[4];
Coordinate receiver;

char buffer[255];
int main(){
    fcoordinate.open("../input coordinate.txt");
    fdistance.open("../output distance.txt");
    if (!fcoordinate.is_open()){
        cout << "Không mở được file input" << endl;
        return 1;
    }
    for (int i = 0; i < 5; i++){
        char eater[255];
        if (i < 4) {
            fcoordinate >> eater >> gateway[i];
            cout << gateway[i] << endl;
        }
        else {
            fcoordinate >> eater >> receiver;
            cout << receiver << endl;
        }
    }
    for (int i = 0; i < 4; i++){
        fdistance << fixed << setprecision(12) << gateway[i].euclideanDistance(receiver) << endl;
    }
    // for(int i = 0; i < 2; i++){
    //     cin >> gateway[i];
    // }
    // cout << gateway[0].euclideanDistance(gateway[1]) << endl;
    fdistance.close();
    fcoordinate.close();
    return 0;
}