#include "base_coordinate.h"

baseCoordinate::baseCoordinate() : x(0), y(0), z(0) {};
baseCoordinate::~baseCoordinate(){};
baseCoordinate::baseCoordinate(double x, double y, double z) : x(x), y(y), z(z) {};

void baseCoordinate::setCoordinate(int x, int y, int z){
    this->x = x;
    this->y = y;
    this->z = z;
}

baseCoordinate& baseCoordinate::operator=(const baseCoordinate& b){
    this->x = b.x;
    this->y = b.y;
    this->z = b.z;
};

void baseCoordinate::parseBuf(const char* buffer){
    int count = sscanf(buffer, " (%lf;%lf;%lf)", &this->x, &this->y, &this->z);
    if (count != 3) std::cout << "Input bị lỗi, vui lòng kiểm tra lại" << std::endl;
}
std::istream & operator >> (std::istream &in, baseCoordinate &c){
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

std::ostream & operator << (std::ostream &out, const baseCoordinate &c){
    // (x; y; z)
    out << std::setprecision(12) << "(" << c.x << "; " << c.y << "; " << c.z << ")";
    return out;
}

void baseCoordinate::printOut(){
    std::cout << *this << std::endl;
}

double baseCoordinate::getX(){
    return x;
}
double baseCoordinate::getY(){
    return y;
};
double baseCoordinate::getZ(){
    return z;
};


