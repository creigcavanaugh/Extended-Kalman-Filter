#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include <iostream>
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  VectorXd h = VectorXd(3);

  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  //Prevent div by zero
  if((fabs(px) < 0.0001) && (fabs(py) < 0.0001)){
    px = 0.0001;
    py = 0.0001;
  }
  else if (fabs(px) < 0.0001) {
    px = 0.0001;
  }

  float rho = sqrt((px * px) + (py * py));

  //Prevent div by zero
  if(fabs(rho) < 0.0001){
    rho = 0.0001;
  }

  float phi = atan2(py, px);

  float rho_dot = ((px * vx) + (py * vy))/rho;

  VectorXd hx(3);
  hx << rho, phi, rho_dot;

  VectorXd y = z - hx;

  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;


  //Sub-optimal loop
  // while (y[1] > M_PI){
  //   y[1] -= (2 * M_PI);
  // }

  // while (phi < (M_PI * -1)){
  //   y[1] -= (2 * M_PI);
  // }

  //Optimized angle normalization
  //https://stackoverflow.com/questions/24234609/standard-way-to-normalize-an-angle-to-%CF%80-radians-in-java
  y[1] -= (2 * M_PI) * floor((y[1] + M_PI) / (2 * M_PI));


  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;


}
