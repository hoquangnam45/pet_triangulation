#include "polar_coordinate.h"
#include "coordinate.h"

double polarCoordinate::getLat(){
    return lat;
};

double polarCoordinate::getLong(){
    return lon;
}

double polarCoordinate::getHeight(){
    return height;
}

polarCoordinate::polarCoordinate(const polarCoordinate& b){
    *this = b;
}

polarCoordinate polarCoordinate::ECEFtoWGS84(Coordinate& ECEFCoor){
    double  x = ECEFCoor.getX(),
            y = ECEFCoor.getY(),
            z = ECEFCoor.getZ();
    double a = EQUATORIAL_RADIUS_METER;
    double b = POLAR_RADIUS_METER;
    double e = sqrt(1 - pow(b/a, 2));
    double ep = sqrt(pow(a/b, 2) - 1);
    double p = sqrt(x*x + y*y);
    double longtitude = atan2(ECEFCoor.getY(), ECEFCoor.getX());

    double k0 = 1. / (1 - e*e);
    double k = k0;
    int idx = 0;
    double diff = 1000;
    double min_diff = diff;
    while(diff > NEWTON_DIFF_THRESH && idx < MAX_NEWTON_ITER){
        double c = pow(p*p + (1-e*e) * z*z * k*k, 3./2.) / (a * e*e);
        double k_temp = 1 + (p*p + (1-e*e) * z*z * k*k*k) / (c - p*p);
        diff = std::abs(k_temp-k);
        if (diff < min_diff){
            k = k_temp;
            min_diff = diff;
        }
        else 
            break;
        idx++;
    }
    double height = 1./(e*e) * (1./k - 1./k0) * sqrt(p*p + z*z * k*k);
    double latitude = atan2(k, p/z);

    // convert rad to degree
    longtitude = longtitude / PI * 180.;
    latitude = latitude / PI * 180.; 
    return polarCoordinate(latitude, longtitude, height);
}

polarCoordinate polarCoordinate::ENUtoWGS84(Coordinate& refPoint, 
                                            Coordinate& localENU){
    Coordinate originPoint = Coordinate::findOriginENU(refPoint);
    
    double  lat = originPoint.getPolar().getLat(),
            lon = originPoint.getPolar().getLong(),
            height = originPoint.getPolar().getHeight();

    // Calculate ECEF coordinate of origin point
    polarCoordinate originPolar = originPoint.getPolar();
    Coordinate originECEF = Coordinate::WGS84toECEF(originPolar);
    // std::cout << originECEF << std::endl;

    // Calculate rotation matrix axis
    arma::mat rot = {
        {-sind(lon), -cosd(lon) * sind(lat), cosd(lon) * cosd(lat)},
        {cosd(lon), -sind(lon) * sind(lat), sind(lon) * cosd(lat)},
        {0, cosd(lat), sind(lat)}
    };

    // Convert local ENU to ECEF
    arma::mat mlocalECEF= (rot * localENU.toMat()) + originECEF.toMat(); 
    Coordinate localECEF = Coordinate::toCoordinate(mlocalECEF);
    // std::cout << localECEF << std::endl;

    return ECEFtoWGS84(localECEF);
};


