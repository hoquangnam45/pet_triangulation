#ifndef COORDINATE_H
#define COORDINATE_H

#include <vector>
#include "armadillo_include/armadillo"
#include "base_coordinate.h"
#include "polar_coordinate.h"
#define PRINT_NEWTON_DIFF 0

class polarCoordinate;

class Coordinate: public baseCoordinate{
    public:
        Coordinate(double x, double y, double z);
        Coordinate();
        Coordinate(const Coordinate&);
        ~Coordinate();
        
        using baseCoordinate::setCoordinate;
        void setCoordinate(Coordinate&);

        void setPolar(polarCoordinate& b);
        polarCoordinate& getPolar();

        Coordinate& operator=(const Coordinate& b);
        bool operator==(Coordinate& b);
        bool operator!=(Coordinate& b);

        static double euclideanDistance(Coordinate &a, Coordinate& b);
        double euclideanDistance(Coordinate& b);
        static double euclideanDistanceSquared(Coordinate &a, Coordinate& b);
        double euclideanDistanceSquared(Coordinate& b);

        static Coordinate triangulate(
            std::vector<double>& dist_to_beacon,
            std::vector<Coordinate>& beacon 
        );

        arma::mat toMat();
        static Coordinate toCoordinate(arma::mat&);

        static Coordinate WGS84toENU(Coordinate& refPoint, polarCoordinate& WGS84Point);
        static Coordinate WGS84toECEF(polarCoordinate&);
        static Coordinate findOriginENU(Coordinate& refPoint);
    private:
        polarCoordinate self_polar;
        static constexpr double RADIUS_SAME_POINT = 1e-1;
};
#endif