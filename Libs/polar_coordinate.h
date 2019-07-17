#ifndef POLAR_COORDINATE_H
#define POLAR_COORDINATE_H

#include "base_coordinate.h"

class Coordinate;
class polarCoordinate: public baseCoordinate{
    public:
        polarCoordinate(double lat, double lon, double height): baseCoordinate(lat, lon, height){};
        polarCoordinate(): baseCoordinate(){};
        polarCoordinate(const polarCoordinate&);
        ~polarCoordinate(){};

        polarCoordinate& operator=(const polarCoordinate& b){
            baseCoordinate::operator=(b);
        }

        double getLat();
        double getLong();
        double getHeight();

        static polarCoordinate ENUtoWGS84(Coordinate& refPoint, Coordinate& localENU);
        static polarCoordinate ECEFtoWGS84(Coordinate&);
    private:
        double& lat = x;
        double& lon = y;
        double& height = z;
};
#endif