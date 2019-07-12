#include "polar_coordinate.h"

double polarCoordinate::getLat(){
    return lat;
};
double polarCoordinate::getLong(){
    return lon;
}
double polarCoordinate::getHeight(){
    return height;
}

double polarCoordinate::sind(double degree){
    return sin(degree / 180. * PI);
}
double polarCoordinate::cosd(double degree){
    return cos(degree / 180. * PI);
}

// WSG84 model
// Reference: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
// local Cartersian is ENU
// *Extra read: the gravity vector WILL point directly towards the centre of mass 
//              of the oblate spheroid (not normal to the surface). It is the resultant force 
//              of the gravity vector and centrifugal force vector that causes the force that 
//              we experience to be normal to the surface of the Earth. This is called "Effective" 
//              or "Apparent" gravity.
// A reference: http://www-gpsg.mit.edu/12.201_12.501/BOOK/chapter2.pdf

polarCoordinate polarCoordinate::ENUtoWSG84(polarCoordinate& refPolarWSG84, 
                                            Coordinate& localCartPoint){
    double  lat = refPolarWSG84.getLat(),
            lon = refPolarWSG84.getLong(),
            height = refPolarWSG84.getHeight();
    // double  lat = 45.,
    //         lon = 0.,
    //         height = 0.;

    // Calculate ECEF coordinate of reference point
    double N =  pow(EQUATORIAL_RADIUS_METER,2) / 
                sqrt(
                    pow(EQUATORIAL_RADIUS_METER,2) * pow(cosd(lat), 2) +
                    pow(POLAR_RADIUS_METER,2) * pow(sind(lat), 2)
                );
    double refOx_ECEF = (N + height) * cosd(lat) * cosd(lon);
    double refOy_ECEF = (N + height) * cosd(lat) * sind(lon);
    double refOz_ECEF = (N * pow(POLAR_RADIUS_METER / EQUATORIAL_RADIUS_METER,2) + height) * sind(lat);

    arma::mat refECEF = {refOx_ECEF, refOy_ECEF, refOz_ECEF};
    std::cout << refECEF << std::endl;

    // Calculate rotation matrix axis
    arma::mat rot = {
        {-sind(lon), -cosd(lon) * sind(lat), cosd(lon) * cosd(lat)},
        {cosd(lon), -sind(lon) * sind(lat), sind(lon) * cosd(lat)},
        {0, cosd(lat), sind(lat)}
    };

    // Convert local Cart to ECEF
    Coordinate localCartPoint(15.,45.,1000.);
    arma::mat localECEF = (rot * localCartPoint.toMat()).t() + refECEF;
    std::cout << localECEF << std::endl;

    return ECEFtoWSG84(Coordinate::toCoordinate(localECEF));
};
polarCoordinate polarCoordinate::ECEFtoWSG84(Coordinate ECEFcoor){
    
}
// Coordinate polarCoordinate::PolartoCartersian(Coordinate& refPoint, 
//                                               polarCoordinate& refPolar, 
//                                               polarCoordinate& convertPoint){
//     double Ox;
//     double Oy;
//     double Oz;

//     return Coordinate(Ox, Oy, Oz);
// };
