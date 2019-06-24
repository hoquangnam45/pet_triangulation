#include <vector>
#include "armadillo_include/armadillo"
#include "base_coordinate.h"
class Coordinate: public baseCoordinate{
    public:
        Coordinate(double x, double y, double z): baseCoordinate(x, y, z){};
        Coordinate(): baseCoordinate(){};
        ~Coordinate(){};
        
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
    private:
        static constexpr double NEWTON_DIFF_THRESH = 1e-15;
        static const int MAX_NEWTON_ITER = 1000;
};