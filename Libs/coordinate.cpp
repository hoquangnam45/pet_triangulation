#include "coordinate.h"
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <string.h>
Coordinate::Coordinate() : x(0), y(0), z(0) {};
Coordinate::~Coordinate(){};
Coordinate::Coordinate(int x, int y, int z) : x(x), y(y), z(z) {};

void Coordinate::setCoordinate(int x, int y, int z){
    this->x = x;
    this->y = y;
    this->z = z;
}
double Coordinate::getX(){
    return x;
}
double Coordinate::getY(){
    return y;
};
double Coordinate::getZ(){
    return z;
};
double Coordinate::euclideanDistance(Coordinate &a, Coordinate& b){
    double  diffx = a.x - b.x,
            diffy = a.y - b.y,
            diffz = a.z - b.z;
    return sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
}
double Coordinate::euclideanDistance(Coordinate& b){
    return euclideanDistance(*this, b);
}
void Coordinate::printOut(){
    std::cout << *this << std::endl;
}
Coordinate& Coordinate::operator=(const Coordinate& b){
    this->x = b.x;
    this->y = b.y;
    this->z = b.z;
};
void Coordinate::parseBuf(const char* buffer){
    // const char* temp = nullptr;
    // while (*buffer){
    //     if (*buffer != ' '){
    //         temp = buffer;
    //         break;
    //     }
    //     buffer++;
    // }
    // if (temp == nullptr) return; 
    int count = sscanf(buffer, " (%lf;%lf;%lf)", &this->x, &this->y, &this->z);
    if (count != 3) std::cout << "Input bị lỗi, vui lòng kiểm tra lại" << std::endl;
}

std::ostream & operator << (std::ostream &out, const Coordinate &c){
    // (x; y; z)
    out << std::setprecision(12) << "(" << c.x << "; " << c.y << "; " << c.z << ")";
    return out;
}
std::istream & operator >> (std::istream &in,  Coordinate &c){
    // (x; y; z)
    char buffer[255];
    char eater;
    int i = 0;
    while(in.peek() == ' ') in.get();
    while (in.peek() && i < sizeof(buffer) - 1){
        if ((buffer[i] = in.get()) == ')') {
            i++;
            break;
        }
        i++;
    }
    buffer[i] = '\0';
    c.parseBuf(buffer);
    return in;
}

Coordinate Coordinate::triangulate(Coordinate& ref_pointd1, double d1, double d2, double d3, double d4){

}