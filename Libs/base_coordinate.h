#ifndef BASE_COORDINATE_H
#define BASE_COORDINATE_H

#include <iostream>
#include <iomanip>
#include <cmath>    
#include <stdio.h>
#include <string.h>

// Parameter reference: GPS.G1-X-00006.pdf
// X, Y, Z measure in meters
// lat, long, height measure in degree
class baseCoordinate{
    public:
        // constructor
        baseCoordinate(double x, double y, double z);
        baseCoordinate(const baseCoordinate&);
        baseCoordinate();
        virtual ~baseCoordinate();
        // set
        void setCoordinate(double x, double y, double z);
        // get
        double getX();
        double getY();
        double getZ();
        // copy
        baseCoordinate& operator=(const baseCoordinate& b);
        // read;
        friend std::istream & operator >> (std::istream &in,  baseCoordinate &c);
        // output
        friend std::ostream & operator << (std::ostream &out,  const baseCoordinate &c);
    protected:
        double x, y, z;
        static double sind(double degree);
        static double cosd(double degree);
        static constexpr double EQUATORIAL_RADIUS_METER = 6378137;
        static constexpr double POLAR_RADIUS_METER = 6356752.31424518;
        static constexpr double PI = 3.14159265358979323846264;
        static constexpr double NEWTON_DIFF_THRESH = 1e-6;
        static const int MAX_NEWTON_ITER = 20;
    private:
        bool parseBuf(const char*);
        void printOut();
};
#endif