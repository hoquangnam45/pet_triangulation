#include "base_coordinate.h"

baseCoordinate::baseCoordinate(){
    setCoordinate(0,0,0);
}

baseCoordinate::baseCoordinate(const baseCoordinate& b){
    *this = b;
}

baseCoordinate::~baseCoordinate(){};

baseCoordinate::baseCoordinate(double x, double y, double z){
    setCoordinate(x,y,z);
}

void baseCoordinate::setCoordinate(double x, double y, double z){
    this->x = x;
    this->y = y;
    this->z = z;
}

baseCoordinate& baseCoordinate::operator=(const baseCoordinate& b){
    setCoordinate(b.x, b.y, b.z);
};

bool baseCoordinate::parseBuf(const char* buffer){
    double x, y, z;
    int count = sscanf(buffer, " (%lf;%lf;%lf)", &x, &y, &z);
    if (count != 3){ 
        std::cout << "Input error, please check coordinate input" << std::endl;
        return false;
    }
    else{
        setCoordinate(x, y, z);
        return true;
    }
}

std::istream & operator >> (std::istream &in, baseCoordinate &c){
    // (x; y; z)
    while(in.peek() == ' ') in.get();
    if (in.peek() != '('){
        in.setstate(std::ios_base::failbit);
        return in;
    }
    int i = 0;
    char buffer[255];
    while (i < sizeof(buffer) - 1){
        if ((buffer[i] = in.get()) == ')') {
            i++;
            break;
        }
        i++;
    }
    buffer[i] = '\0';
    if(!c.parseBuf(buffer))
        in.setstate(std::ios_base::failbit);
    return in;
}

std::ostream & operator << (std::ostream &out, const baseCoordinate &c){
    // (x; y; z)
    out << "(" << c.x << "; " << c.y << "; " << c.z << ")";
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

double baseCoordinate::sind(double degree){
    return sin(degree / 180. * PI);
}

double baseCoordinate::cosd(double degree){
    return cos(degree / 180. * PI);
}