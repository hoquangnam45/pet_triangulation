#include "coordinate.h"

Coordinate::Coordinate(double x, double y, double z){
    setCoordinate(x, y, z);
}

Coordinate::Coordinate(){
    setCoordinate(0,0,0);
}

Coordinate::Coordinate(const Coordinate& b){
    *this = b;
}

Coordinate::~Coordinate(){}

void Coordinate::setPolar(polarCoordinate& b){
    self_polar = b;
}

void Coordinate::setCoordinate(Coordinate& b){
    setCoordinate(b.x, b.y, b.z);
}

polarCoordinate& Coordinate::getPolar(){
    return self_polar;
}

Coordinate& Coordinate::operator=(const Coordinate& b){
    baseCoordinate::operator=(b);
    this->self_polar = b.self_polar;
}

bool Coordinate::operator==(Coordinate& b){
    if (this->euclideanDistance(b) <= RADIUS_SAME_POINT) return true;
    return false;
}

bool Coordinate::operator!=(Coordinate& b){
    return !(*this == b);
}

double Coordinate::euclideanDistance(Coordinate &a, Coordinate& b){
    double  diffx = a.x - b.x,
            diffy = a.y - b.y,
            diffz = a.z - b.z;
    return sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
}

double Coordinate::euclideanDistance(Coordinate& b){
    return euclideanDistance(*this, b);
}

double Coordinate::euclideanDistanceSquared(Coordinate &a, Coordinate& b){
    double  diffx = a.x - b.x,
            diffy = a.y - b.y,
            diffz = a.z - b.z;
    return diffx * diffx + diffy * diffy + diffz * diffz;
}

double Coordinate::euclideanDistanceSquared(Coordinate& b){
    return euclideanDistanceSquared(*this, b);
}

arma::mat Coordinate::toMat(){
    arma::mat ret(3,1, arma::fill::zeros);
    ret(0,0) = x;
    ret(1,0) = y;
    ret(2,0) = z;
    return ret;
}

Coordinate Coordinate::toCoordinate(arma::mat& Matrix){
    if (arma::size(Matrix) == arma::size(3,1)){
        double  x = Matrix(0,0),
                y = Matrix(1,0),
                z = Matrix(2,0);
        Coordinate ret(x, y, z);
        //std::cout << "*" << ret << std::endl;
        return ret;  
    }
    else if(arma::size(Matrix) == arma::size(1,3)){
        double  x = Matrix(0,0),
                y = Matrix(0,1),
                z = Matrix(0,2);
        Coordinate ret(x, y, z);
        //std::cout << "*" << ret << std::endl;
        return ret;  
    }
    else{
        throw std::runtime_error("Can't convert matrix to coordinate");
    }

}
// WGS84 model
// Reference: https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
Coordinate Coordinate::WGS84toECEF(polarCoordinate& coor){
    double  lat = coor.getLat(),
            lon = coor.getLong(),
            height = coor.getHeight();
    double a = EQUATORIAL_RADIUS_METER;
    double b = POLAR_RADIUS_METER;
    double e = sqrt(1 - pow(b/a, 2));
    double N =  a / sqrt(1-pow(e * sind(lat), 2));
    double x = (N + height) * cosd(lat) * cosd(lon);
    double y = (N + height) * cosd(lat) * sind(lon);
    double z = (N * pow(b / a, 2) + height) * sind(lat);
    return Coordinate(x, y, z);
}

Coordinate Coordinate::WGS84toENU(Coordinate& refPoint, polarCoordinate& WGS84Point){
    Coordinate originPoint = findOriginENU(refPoint);

    double  lat = originPoint.getPolar().getLat(),
            lon = originPoint.getPolar().getLong(),
            height = originPoint.getPolar().getHeight();

    // Calculate ECEF coordinate of origin point
    polarCoordinate originPolar = originPoint.getPolar();
    Coordinate originECEF = Coordinate::WGS84toECEF(originPolar);
    //std::cout << originECEF << std::endl;

    // Calculate rotation matrix axis
    arma::mat rot = {
        {-sind(lon), -cosd(lon) * sind(lat), cosd(lon) * cosd(lat)},
        {cosd(lon), -sind(lon) * sind(lat), sind(lon) * cosd(lat)},
        {0, cosd(lat), sind(lat)}
    };

    // Convert local WGS84 to ECEF
    Coordinate localECEF = Coordinate::WGS84toECEF(WGS84Point);

    arma::mat localENU = rot.t() * (localECEF.toMat() - originECEF.toMat());
    return Coordinate::toCoordinate(localENU);
}

Coordinate Coordinate::findOriginENU(Coordinate& refPoint){
    Coordinate originPoint;
    if (originPoint != refPoint){
        Coordinate newRefPoint;
        
        newRefPoint.setCoordinate(0,0,0);
        polarCoordinate temp = refPoint.getPolar();
        newRefPoint.setPolar(temp);

        originPoint.setCoordinate(-refPoint.getX(), -refPoint.getY(), -refPoint.getZ());

        polarCoordinate originPolar = polarCoordinate::ENUtoWGS84(newRefPoint, originPoint);
        originPoint.setPolar(originPolar);
        originPoint.setCoordinate(0,0,0);
    }
    else{
        originPoint = refPoint;
    }
    return originPoint;
}

void reduceMatrix(arma::mat& Matrix);
void backSubstitution(arma::mat& equation, arma::mat& ret);
void evaluateJacobian(
    //std::vector<double>& dist_to_beacon,
    std::vector<Coordinate>& beacon,
    arma::mat& theta,
    arma::mat& ret
);
void evaluateF(
    std::vector<double>& dist_to_beacon,
    std::vector<Coordinate>& beacon,
    arma::mat& theta,
    arma::mat& ret
);

// Implementation detail: TurgutOzal-11-Trilateration.pdf
Coordinate Coordinate::triangulate(
    std::vector<double>& dist_to_beacon,
    std::vector<Coordinate>& beacon
){
    if (dist_to_beacon.size() != beacon.size() | dist_to_beacon.size() < 4){
        if (dist_to_beacon.size() != beacon.size()) 
            throw std::runtime_error("Inconsistent number of beacon and distance");
        if (dist_to_beacon.size() < 4)
            throw std::runtime_error("Not enough beacon to triangulate");
    }
    /* Calculate initial guess */
    // Set ref point
    std::vector<double> clone_dist_to_beacon(dist_to_beacon);
    std::vector<Coordinate> clone_beacon(beacon);
    Coordinate ref_point = clone_beacon.back();
    double d_target_ref = clone_dist_to_beacon.back();
    clone_beacon.pop_back();
    clone_dist_to_beacon.pop_back();

    // Construct matrix b
    arma::mat d_beacon_target(clone_dist_to_beacon);
    arma::mat d_beacon_ref(clone_dist_to_beacon.size(), 1, arma::fill::zeros);
    for(int i = 0; i < clone_dist_to_beacon.size(); i++){
        d_beacon_ref(i,0) = ref_point.euclideanDistance(clone_beacon[i]);
    }
    arma::mat b = 1./2 * (pow(d_target_ref, 2) - pow(d_beacon_target,2) + pow(d_beacon_ref,2));

    // Construct matrix A
    arma::mat A(clone_dist_to_beacon.size(), 3, arma::fill::zeros);
    for (int i = 0; i < clone_dist_to_beacon.size(); i++){
        A(i, 0) = clone_beacon[i].getX() - ref_point.getX();
        A(i, 1) = clone_beacon[i].getY() - ref_point.getY();
        A(i, 2) = clone_beacon[i].getZ() - ref_point.getZ();
    }
    
    // Construct matrix condition and ret
    arma::mat cond_matrix = A.t() * A;
    arma::mat b_cond_matrix = A.t() * b;
    arma::mat X;

    // Test condition of A
    if (cond(cond_matrix) < 5){
        X = inv(cond_matrix) * b_cond_matrix + ref_point.toMat();
    }
    else{
        // Construct Q and R
        arma::mat U;
        arma::vec s;
        arma::mat V;

        arma::svd(U,s,V,cond_matrix);
        V = arma::diagmat(s) * V.t();
        b_cond_matrix = U.t() * b_cond_matrix;
        arma::mat concat = arma::join_rows(V, b_cond_matrix);
        reduceMatrix(concat);
        backSubstitution(concat, X);
        X = X + ref_point.toMat();
    }
    std::cout << "initial guessed coordinate: \n" << toCoordinate(X) << std::endl;
    
    arma::mat f;
    arma::mat jacobian;
    double diff = 1000;
    double min_diff = diff;
    int iter = 0;
    // Solve by newton method
    while (diff > NEWTON_DIFF_THRESH && iter < MAX_NEWTON_ITER){
        evaluateJacobian(beacon, X, jacobian);
        evaluateF(dist_to_beacon, beacon, X, f);
        //std::cout << jacobian << f << beacon.size() << dist_to_beacon.size()<< std::endl;
        arma::mat temp = X - inv(jacobian.t() * jacobian) * jacobian.t() * f;
        diff = norm(X - temp);
        #if PRINT_NEWTON_DIFF != 0
        std::cout << "iter " << iter + 1 << " diff: "<< diff << std::endl;//"\n>>>>\n" << X << "\n***\n" << temp << std::endl;
        #endif
        if (diff < min_diff) {
            min_diff = diff;
            X = temp;
        }
        else
            break;
        iter++;
        //std::cout << X << std::endl;
    }
    //std::cout << iter << std::endl;
    std::cout << "refined coordinate: \n" << toCoordinate(X) << std::endl;
    Coordinate ret = toCoordinate(X); 
    return ret;
}

void reduceMatrix(arma::mat& Matrix){
    int nrow = Matrix.n_rows;
    int ncol = Matrix.n_cols;
    int min_dimension = nrow < ncol ? nrow : ncol;
    int lead = 0;
    //std::cout << Matrix << std::endl;
    while (lead < min_dimension){
        for (int i = lead; i < nrow; i++){ 
            double temp_lead = Matrix(i, lead),
                   temp_ratio = temp_lead / Matrix(lead, lead);
            for (int j = lead; j < ncol; j++){
                if(i == lead) Matrix(i, j) /= temp_lead;
                else Matrix(i, j) -=  temp_ratio * Matrix(lead, j);
            }
        }
        //std::cout << Matrix << std::endl;
        lead++;
    }
}

void backSubstitution(arma::mat& equation, arma::mat& ret){
    // Check full rank
    int nrow = equation.n_rows;
    int ncol = equation.n_cols;
    //int min_dimension = nrow < ncol ? nrow : ncol;
    if (rank(equation) < ncol - 1) throw std::runtime_error("Can't back substitute since matrix is so singular");
    
    if (ret.empty() || arma::size(ret) != arma::size(ncol - 1, 1)) {
        ret.clear();
        ret = arma::mat(ncol - 1, 1);
    }
    for (int i = ncol - 2; i >= 0; i--){
        if (i == ncol - 2) ret(i, 0) = equation(i, ncol - 1);
        else {
            ret(i, 0) = equation(i, ncol - 1);
            for (int j = ncol - 2; j >= i + 1; j--) 
                ret(i, 0) -= equation(i, j) * ret(j, 0);
            //ret(i, 0) /= equation(i, i);
        }
    }
}

void evaluateJacobian(
    //std::vector<double>& dist_to_beacon,
    std::vector<Coordinate>& beacon,
    arma::mat& theta,
    arma::mat& ret
){
    if (beacon.size() < 4){
        throw std::runtime_error("not enough beacon");
    }

    if (arma::size(theta) != arma::size(3,1)) throw std::runtime_error("Theta dimension is illegal");
         
    if (ret.empty() || arma::size(ret) != arma::size(beacon.size(), 3)){
        ret.clear();
        ret = arma::mat(beacon.size(), 3);
    }

    for (int i = 0; i < beacon.size(); i++){
        //std::cout << theta <<std::endl;
        double denominator =    sqrt(
                                    pow(theta(0, 0) - beacon[i].getX(), 2) + 
                                    pow(theta(1, 0) - beacon[i].getY(), 2) +
                                    pow(theta(2, 0) - beacon[i].getZ(), 2) 
                                );
        ret(i, 0) = (theta(0, 0) - beacon[i].getX()) / denominator;
        ret(i, 1) = (theta(1, 0) - beacon[i].getY()) / denominator;
        ret(i, 2) = (theta(2, 0) - beacon[i].getZ()) / denominator;
    }
}

void evaluateF(
    std::vector<double>& dist_to_beacon,
    std::vector<Coordinate>& beacon,
    arma::mat& theta,
    arma::mat& ret
){
    if (dist_to_beacon.size() != beacon.size() | dist_to_beacon.size() < 4){
        throw std::runtime_error("Inconsistent amount of distance and beacon or not enough beacon");
    }

    if (arma::size(theta) != arma::size(3,1)) throw std::runtime_error("Theta dimension is illegal");
         
    if (ret.empty() || arma::size(ret) != arma::size(beacon.size(), 1)){
        ret.clear();
        ret = arma::mat(beacon.size(), 1);
    }

    for (int i = 0; i < beacon.size(); i++){
        ret(i, 0) = sqrt(
                            pow(theta(0, 0) - beacon[i].getX(), 2) + 
                            pow(theta(1, 0) - beacon[i].getY(), 2) +
                            pow(theta(2, 0) - beacon[i].getZ(), 2) 
                        ) - dist_to_beacon[i];
    }
}

