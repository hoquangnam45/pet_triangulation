#include "coordinate.h"
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
    double  x = Matrix(0,0),
            y = Matrix(1,0),
            z = Matrix(2,0);
    Coordinate ret(x, y, z);
    //std::cout << "*" << ret << std::endl;
    return ret;
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
        Coordinate dummy(-1, -1, -1);
        return dummy;
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
    std::cout << X << std::endl;
    
    arma::mat f;
    arma::mat jacobian;
    double diff = 1000;
    int iter = 0;
    // Solve by newton method
    while (diff > NEWTON_DIFF_THRESH && iter < MAX_NEWTON_ITER){
        evaluateJacobian(beacon, X, jacobian);
        evaluateF(dist_to_beacon, beacon, X, f);
        //std::cout << jacobian << f << beacon.size() << dist_to_beacon.size()<< std::endl;
        arma::mat temp = X;
        X -= inv(jacobian.t() * jacobian) * jacobian.t() * f;
        diff = norm(X - temp); 
        std::cout << "iter " << iter + 1 << ": "<< diff << std::endl;//"\n>>>>\n" << X << "\n***\n" << temp << std::endl;
        iter++;
        //std::cout << X << std::endl;
    }
    std::cout << X << std::endl;
    Coordinate ret = toCoordinate(X); 
    std::cout << ret << std::endl;
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
    if (ret.empty()) ret = arma::mat(3,1);
    ret(2, 0) = equation(2, 3);
    ret(1, 0) = equation(1, 3) - ret(2, 0) * equation(1, 2);
    ret(0, 0) = equation(0, 3) - ret(2, 0) * equation(0, 2) - ret(1, 0) * equation(0, 1);
}

void evaluateJacobian(
    //std::vector<double>& dist_to_beacon,
    std::vector<Coordinate>& beacon,
    arma::mat& theta,
    arma::mat& ret
){
    //std::cout << ret.empty() << std::endl;
    if (ret.empty()) ret = arma::mat(beacon.size(), 3);
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
    if (ret.empty()) ret = arma::mat(beacon.size(), 1);
    for (int i = 0; i < beacon.size(); i++){
        ret(i, 0) = sqrt(
                            pow(theta(0, 0) - beacon[i].getX(), 2) + 
                            pow(theta(1, 0) - beacon[i].getY(), 2) +
                            pow(theta(2, 0) - beacon[i].getZ(), 2) 
                        ) - dist_to_beacon[i];
    }
}