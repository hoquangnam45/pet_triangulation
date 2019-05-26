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
    // const char* temp = nullptr;
    // while (*buffer){
    //     if (*buffer != ' '){
    //         temp = buffer;
    //         break;
    //     }
    //     buffer++;
    // }
    // if (temp == nullptr) return; 
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
arma::Mat<double> Coordinate::toMat(){
    arma::Mat<double> ret(3,1, arma::fill::zeros);
    ret(0,0) = x;
    ret(1,0) = y;
    ret(2,0) = z;
    return ret;
}


// template<class Matrix>
// void print_matrix(Matrix matrix) {
//     matrix.print(std::cout);
// }

//provide explicit instantiations of the template function for
//every matrix type you use somewhere in your program.
//template void print_matrix<arma::mat>(arma::mat matrix);
//template void print_matrix<arma::Mat<double>>(arma::Mat<double> matrix);
//template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);
void reduceMatrix(arma::mat& Matrix);
arma::mat backSubstitution(arma::mat& A, arma::mat& b);
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
    arma::Mat<double> d_beacon_target(dist_to_beacon);
    arma::Mat<double> d_beacon_ref(dist_to_beacon.size(), 1, arma::fill::zeros);
    for(int i = 0; i < dist_to_beacon.size(); i++){
        d_beacon_ref(i,0) = ref_point.euclideanDistance(beacon[i]);
    }
    arma::Mat<double> b = 1./2 * (pow(d_target_ref, 2) - pow(d_beacon_target,2) + pow(d_beacon_ref,2));
    // std::cout << std::fixed << pow(d_target_ref, 2) << std::endl;
    // std::cout << std::fixed << d_beacon_target << std::fixed << pow(d_beacon_target,2) << std::endl;
    // std::cout << std::fixed << d_beacon_ref << std::fixed << pow(d_beacon_ref,2) << std::endl;

    // Construct matrix A
    arma::Mat<double> A(dist_to_beacon.size(), 3, arma::fill::zeros);
    for (int i = 0; i < dist_to_beacon.size(); i++){
        //std::cout << beacon[i] << " " << refprint_matrix<arma::Mat<double> >_point << std::endl;
        //std::cout << beacon[i].getX() << " " << beacon[i].getY() << " " << beacon[i].getZ() << std::endl;
        A(i, 0) = beacon[i].getX() - ref_point.getX();
        A(i, 1) = beacon[i].getY() - ref_point.getY();
        A(i, 2) = beacon[i].getZ() - ref_point.getZ();
    }
    // std::cout << ref_point << std::endl;

    // Check full rank
    //std::cout << rank(Matrix) << std::endl;
    //std::cout << inv(Matrix) << std::endl;
    // int nrow = A.t() * A().n_rows;
    // int ncol = (A.t() * A).n_cols;
    // int min_dimension = nrow < ncol ? nrow : ncol;
    // throw std::runtime_error("TEST");

    // Test condition of A
    arma::Mat<double> X;
    // std::cout << cond(A.t() * A) << std::endl;
    if (cond(A.t() * A) < 10){
        // Least squared
        X = inv(A.t() * A) * A.t() * b + ref_point.toMat();
    }
    else{
        // Construct Q and R
        arma::mat U;
        arma::vec s;
        arma::mat V;

        arma::svd(U,s,V,A);
        V = arma::diagmat(s) * V;
        b = U.t() * b;
        V = {{1, 2, 3}, {1, 3, 5}, {1, 4, 7}};
        reduceMatrix(V);
        //std::cout << V << std::endl;
        // X = backSubstitution(V, b);
        // X += ref_point.toMat();
        // std::cout << U << std::endl;
        // std::cout << arma::diagmat(s) << std::endl;
        // std::cout << V << std::endl;
        // std::cout << V.t() * V << std::endl;
    }

    // Solve by newton method

    Coordinate dummy(-1, -1, -1);
    return dummy;
    // 
    // arma::Mat<double> X;
    // arma::Mat<double> A;
    // arma::Mat<double> b;
    // // check condition number
    // double condition_number = arma::cond(A);
    // // well condition solution
    // if(well_condition(condition_number))
    // // illed condition solution
    // else{
    //     arma::vec temp = arma::svd(A);

    // /*
    // // Iterate until error is small enough
    // while(error > ERR_THRESH){
    //     arma::Mat<double> Ftheta;
    //     arma::Mat<double> ftheta;
    //     arma::Mat<double> jacobian;
    //     arma::Mat<double> X;
    //     X_temp = X - jacobian(X).inv*ftheta(X); 
    //     error = abs(X_temp - X);
    // }
    // */
    // Coordinate new
}

void reduceMatrix(arma::mat& Matrix){
    int nrow = Matrix.n_rows;
    int ncol = Matrix.n_cols;
    int min_dimension = nrow < ncol ? nrow : ncol;
    int lead = 0;
    // arma::mat U;
    // arma::vec s;
    // arma::mat V;
    // arma::svd(U,s,V,Matrix);
    // std::cout << U << std::endl;
    // std::cout << arma::diagmat(s) << std::endl;
    // std::cout << V << std::endl;
    // std::cout << cond(Matrix) << std::endl;
    while (lead < min_dimension){
        for (int i = lead; i < nrow; i++){ 
            double temp_lead = Matrix(i, lead),
                   temp_ratio = temp_lead / Matrix(lead, lead);
            for (int j = lead; j < ncol; j++){
                //if (!temp_lead || std::isnan(temp_ratio)) continue; 
                if(i == lead) Matrix(i, j) /= temp_lead;
                else Matrix(i, j) -=  temp_ratio * Matrix(lead, j);
            }
        }
        std::cout << Matrix << std::endl;
        lead++;
    }
}

arma::mat backSubstitution(arma::mat& A, arma::mat& b){

}
