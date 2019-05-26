#include <iostream>
#include <vector>
#include "armadillo_include/armadillo"
class Coordinate{
    public:
        Coordinate(int x, int y, int z);
        Coordinate();
        ~Coordinate();
        void setCoordinate(int x, int y, int z);
        double getX();
        double getY();
        double getZ();
        static double euclideanDistance(Coordinate &a, Coordinate& b);
        double euclideanDistance(Coordinate& b);
        static double euclideanDistanceSquared(Coordinate &a, Coordinate& b);
        double euclideanDistanceSquared(Coordinate& b);
        void printOut();
        void parseBuf(const char*);
        static Coordinate triangulate(Coordinate& ref_pointd1, Coordinate& ref_point2, Coordinate& ref_point3, Coordinate& ref_pointd4, double d1, double d2, double d3, double d4);
        Coordinate& operator=(const Coordinate& b);
        friend std::istream & operator >> (std::istream &in,  Coordinate &c);
        friend std::ostream & operator << (std::ostream &out,  const Coordinate &c);
        //Coordinate& convertMode(Coordinate& ref_point_GPS,)
        //Coordinate& convertMode()
        //coordinateMode getMode()
        //Coordinate& setMode(coordinateMode)
        //triangulateWithTime();
        static Coordinate triangulate(
            std::vector<double>& dist_to_beacon,
            std::vector<Coordinate>& beacon 
        );
        arma::mat toMat();
    private:
        double x, y, z;
        const double EARTH_RADIUS_IN_KM = 6371;
        const double NEWTON_ERR_THRESH;
        // coordinateMode = FLAT;
        // enum coordinateMode{
        //     FLAT,
        //     EARTH_GPS
        // }
};