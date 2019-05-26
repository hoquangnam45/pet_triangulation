#include "coordinate.h"

class polarCoordinate{
    public:
        polarCoordinate(double, double, double);
        void setCoordinate(double, double, double);
        double getLat();
        double getLong();
        double getHeight();
        void setLat(double);
        void setLong(double);
        void setHeight(double);
        polarCoordinate CartersiantoPolar()
    private:
        double lat;
        double lon;
        double height;
}