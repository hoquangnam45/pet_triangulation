#include "coordinate.h"
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <string.h>

Coordinate::Coordinate() : x(0), y(0), z(0) {};
Coordinate::~Coordinate(){};
Coordinate::Coordinate(int x, int y, int z) : x(x), y(y), z(z) {};

void Coordinate::setCoordinate(int x, int y, int z){
    this->x = x;
    this->y = y;
    this->z = z;
}
double Coordinate::getX(){
    return x;
}
double Coordinate::getY(){
    return y;
};
double Coordinate::getZ(){
    return z;
};
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
void Coordinate::printOut(){
    std::cout << *this << std::endl;
}
Coordinate& Coordinate::operator=(const Coordinate& b){
    this->x = b.x;
    this->y = b.y;
    this->z = b.z;
};
void Coordinate::parseBuf(const char* buffer){
    int count = sscanf(buffer, " (%lf;%lf;%lf)", &this->x, &this->y, &this->z);
    if (count != 3) std::cout << "Input bị lỗi, vui lòng kiểm tra lại" << std::endl;
}

std::ostream & operator << (std::ostream &out, const Coordinate &c){
    // (x; y; z)
    out << std::setprecision(12) << "(" << c.x << "; " << c.y << "; " << c.z << ")";
    return out;
}
std::istream & operator >> (std::istream &in,  Coordinate &c){
    // (x; y; z)
    char buffer[255];
    char eater;
    int i = 0;
    while(in.peek() == ' ') in.get();
    while (in.peek() && i < sizeof(buffer) - 1){
        if ((buffer[i] = in.get()) == ')') {
            i++;
            break;
        }
        i++;
    }
    buffer[i] = '\0';
    c.parseBuf(buffer);
    return in;
}
arma::mat Coordinate::toMat(){
    arma::mat ret(3,1, arma::fill::zeros);
    ret(0,0) = x;
    ret(1,0) = y;
    ret(2,0) = z;
    return ret;
}

void reduceMatrix(arma::mat& Matrix);
void backSubstitution(arma::mat& equation, arma::mat& ret);

// Implementation detail: TurgutOzal-11-Trilateration.pdf
Coordinate Coordinate::triangulate(
    std::vector<double>& dist_to_beacon,
    std::vector<Coordinate>& beacon
){
    if (dist_to_beacon.size() != beacon.size() | dist_to_beacon.size() < 4){
        Coordinate dummy(-1, -1, -1);
        return dummy;
    }
    /* Calculate initial guess */
    // Set ref point
    Coordinate ref_point = beacon.back();
    double d_target_ref = dist_to_beacon.back();
    beacon.pop_back();
    dist_to_beacon.pop_back();

    // Construct matrix b
    arma::mat d_beacon_target(dist_to_beacon);
    arma::mat d_beacon_ref(dist_to_beacon.size(), 1, arma::fill::zeros);
    for(int i = 0; i < dist_to_beacon.size(); i++){
        d_beacon_ref(i,0) = ref_point.euclideanDistance(beacon[i]);
    }
    arma::mat b = 1./2 * (pow(d_target_ref, 2) - pow(d_beacon_target,2) + pow(d_beacon_ref,2));

    // Construct matrix A
    arma::mat A(dist_to_beacon.size(), 3, arma::fill::zeros);
    for (int i = 0; i < dist_to_beacon.size(); i++){
        A(i, 0) = beacon[i].getX() - ref_point.getX();
        A(i, 1) = beacon[i].getY() - ref_point.getY();
        A(i, 2) = beacon[i].getZ() - ref_point.getZ();
    }
    
    // Construct matrix condition and ret
    arma::mat cond_matrix = A.t() * A;
    arma::mat b_cond_matrix = A.t() * b;
    arma::mat X;

    // Check full rank
    int nrow = cond_matrix.n_rows;
    int ncol = cond_matrix.n_cols;
    int min_dimension = nrow < ncol ? nrow : ncol;
    if (rank(cond_matrix) < min_dimension) throw std::runtime_error("Singular matrix, not full rank");

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

    // Solve by newton method
    // while (err > ERR_THRESH){
    //     temp = X = X - jacobian(X).inv() * f(X);
    //     err = temp - X;
    // }

    Coordinate dummy(-1, -1, -1);
    return dummy;
}

void reduceMatrix(arma::mat& Matrix){
    int nrow = Matrix.n_rows;
    int ncol = Matrix.n_cols;
    int min_dimension = nrow < ncol ? nrow : ncol;
    int lead = 0;
    std::cout << Matrix << std::endl;
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
    ret = arma::mat(3,1);
    ret(2, 0) = equation(2, 3);
    ret(1, 0) = equation(1, 3) - ret(2, 0) * equation(1, 2);
    ret(0, 0) = equation(0, 3) - ret(2, 0) * equation(0, 2) - ret(1, 0) * equation(0, 1);
}
