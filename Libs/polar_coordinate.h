#include "coordinate.h"

class polarCoordinate: public baseCoordinate{
    public:
        polarCoordinate(double lat, double lon, double height): baseCoordinate(lat, lon, height){};
        polarCoordinate(): baseCoordinate(){};
        ~polarCoordinate(){};

        double getLat();
        double getLong();
        double getHeight();

        polarCoordinate CartersiantoPolar(Coordinate& originPoint, polarCoordinate& originPointPolar, Coordinate& convertPoint);
        Coordinate PolartoCartersian(Coordinate& originPoint, polarCoordinate& originPointPolar, Coordinate& convertPoint);
    private:
        double angleLon(double baseLon, double lon);
        double angleLat(double baseLat, double lat);
        double& lat = x;
        double& lon = y;
        double& height = z;
        const double EARTH_RADIUS = 0;
};