#include "coordinate.h"

// Parameter reference: 14gravity1_2.pdf
class polarCoordinate: public baseCoordinate{
    public:
        polarCoordinate(double lat, double lon, double height): baseCoordinate(lat, lon, height){};
        polarCoordinate(): baseCoordinate(){};
        ~polarCoordinate(){};

        double getLat();
        double getLong();
        double getHeight();

        static polarCoordinate ENUtoWSG84(polarCoordinate& refPolarWSG84, Coordinate& localCartPoint);
        static polarCoordinate ENUtoWSG84();
        static polarCoordinate ECEFtoWSG84(Coordinate);
    private:
        double& lat = x;
        double& lon = y;
        double& height = z;
        static double sind(double degree);
        static double cosd(double degree);
        static constexpr double EQUATORIAL_RADIUS_METER = 6378137;
        static constexpr double POLAR_RADIUS_METER = 6356752.3;
        static constexpr double PI = 3.14159265358979323846264;
};