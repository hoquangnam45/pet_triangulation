#include <iostream>
#include <fstream>
using namespace std;
ifstream fcoordinate;
int main(){
    cout << "Chuẩn bị file coordinate.txt" << endl;
    fcoordinate.open("coordinate.txt");
    if (!fcoordinate){
        cerr << "Unable to open file datafile.txt";
        exit(1);   // call system to stop
    }
    for (int i = 0; i < 4; i++){
        
    }
    cout << "1:"
}