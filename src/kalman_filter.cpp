#include "kalman_filter.h"

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

  /**
 * No external motion
 */
  //VectorXd u_;	// external motion
  //u_ = VectorXd(2);
  //u_ << 0, 0;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
    * Udacity sample code
  */

  //x_ = F_ * x_ + u_; //with external motion
  x_ = F_ * x_; //without external motion

  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
    * Udacity sample code
  */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
    * using the state transition function and measurement function
  */

  // get x vector from state transition function
  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];

  //Measurement function f(x')
  double rho = sqrt(px * px + py *py); // radial distance from origin
  double phi = atan(py / px);// angle between rho and x
  double rho_dot = (px * vx + py * vy)/(sqrt(px * px + py *py)); //change of rho (range rate)

  // Create z_pred vector
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
