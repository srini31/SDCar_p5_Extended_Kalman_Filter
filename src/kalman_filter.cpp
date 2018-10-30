#include "kalman_filter.h"
#include <iostream>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  // TODO: predict the state
	x_ = F_ * x_; // + u;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}
////x_, P_, F_, Q_, H_, R_
//this is for Lidar -- H_laser, R_laser
void KalmanFilter::Update(const VectorXd &z) {
  // TODO:update the state by using Kalman Filter equations
  //VectorXd y = z - H * x;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_; //H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  //MatrixXd K =  P * Ht * Si;
  //MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

	//new state 
	x_ = x_ + (K * y);
  	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
// This is for Radar  - Hj, R_radar
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // TODO: update the state by using Extended Kalman Filter equations
  //z containes ro, theta, ro_dot
  
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float val = px*px + py*py;
  //check if px*px + py*py is = 0, div_by_zero err
  if(fabs(val) < 0.0001){
		cout << "UpdateEKF () - Error - Division by Zero" << endl;
		return ;
	}
  float rho = sqrt(val);
  float phi = atan2(py,px);
  float rho_dot = (px*vx+py*vy)/rho; //val;
  
  VectorXd ztemp = VectorXd(3);
  ztemp << rho, phi, rho_dot;
  
  VectorXd y = z - ztemp;
  
  //make sure phi is between -PI and PI
  if(y[1] < -M_PI) {
  	y[1] = y[1] + 2*M_PI;
  } else if(y[1] > M_PI) {
  	y[1] = y[1] - 2*M_PI;
  }
  
  //return y;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_; //H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  //MatrixXd K =  P * Ht * Si;
  //MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new state 
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


