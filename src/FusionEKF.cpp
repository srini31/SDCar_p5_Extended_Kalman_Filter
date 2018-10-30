#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/* Constructor */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  //x, P , R_laser_, R_radar_, H_laser_, Hj_, F
  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 0, 0, 0, 0;
  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
			0, 0, 1000, 0,
			0, 0, 0, 1000;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  //measurement covariance matrix - radar -- 3x3 covariance matrix, sigma2 of rho, P, rho_dot on the diagonal
  //measurement vector has 3 values
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  //measurement matrix laser
  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;
  //measurement matrix
  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
	     1, 1, 1, 1;
  ////x_, P_, F_, Q_, H_, R_ are the variables in kalman_filter.h
  //to initialize, x, P, F,  -- done H, Q, R
  /** TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises */
  // remaining - L 23 laser meas part 3 - tracking.cpp
  //create a 4D state vector, we don't know yet the values of the x state
  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
		  0, 1, 0, 1,
		  0, 0, 1, 0,
		  0, 0, 0, 1;
  
  //set the acceleration noise components
  noise_ax = 9; //5;
  noise_ay = 9; //5;
}

/** Destructor */
FusionEKF::~FusionEKF() {}

/*
first time, initialize state (F), x and covariance  (Q) matrices?
predict - dt, and then new EKF, F, Q - predict x, P
update L or R
L - setup L matrices
R - linearize  meas function h(x_dash), setup R matrices
state update with new  measurement z
*/ 
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) { //only first time
    /** TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
	cout << "Kalman Filter Initialization " << endl;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      //set the state with the initial location and zero velocity -- ok for Laser; what abour Radar?
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) { 
      
      float ro = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];

      ekf_.x_(0) = ro * cos(theta);
      ekf_.x_(1) = ro * sin(theta);
      ekf_.x_(2) = ro_dot * cos(theta);
      ekf_.x_(3)  = ro_dot * sin(theta);
    }
	
	previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  /********************** Prediction ****************************************/
  /** TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  if(dt > 0.001) {
    //Modify the F matrix so that the time is integrated
    ekf_.F_(0, 2) = dt; //kf_
    ekf_.F_(1, 3) = dt; //kf_
  
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    //set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4); //kf_
    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0, //kf_
			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
  
    //F and Q are required for predict
    ekf_.Predict(); //should be same for Radar and Lidar
  }
  /****************** Update *****************************************/
  /**   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.    */
  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) { 
    // Laser updates
    // -- comment this section to test only for Radar
     ekf_.H_ = H_laser_;
     ekf_.R_ = R_laser_;
     ekf_.Update(measurement_pack.raw_measurements_); //(z)
    
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    /** Convert radar from polar to cartesian coordinates and initialize state  */
    //* -- comment this section to test only for Laser
    Hj_  = tools.CalculateJacobian(ekf_.x_); 
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_); //(z)
    //*/
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
} //process measurement
