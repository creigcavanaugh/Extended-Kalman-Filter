#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  MatrixXd P_ = MatrixXd(4,4);
  MatrixXd F_ = MatrixXd(4,4);
  MatrixXd Q_ = MatrixXd(4,4);
  VectorXd x_ = VectorXd(4);

  //measurement covariance matrix - laser
  R_laser_ <<   0.0225, 0,
                0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ <<   0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

  //state covariance matrix P
  P_ <<   1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  //measurement matrix
  H_laser_ <<   1, 0, 0, 0,
                0, 1, 0, 0;

  //the initial transition matrix F_
  F_ <<   1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  //initialize ekf
  ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);

}



/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      ro,theta, ro_dot
      */

      //set the state with the initial location and zero velocity
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];

      float x = rho * cos(phi);
      float y = rho * sin(phi);

      //Initialize state
      ekf_.x_ << x, y, 0, 0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      px, py
      */

      float x = measurement_pack.raw_measurements_[0];
      float y = measurement_pack.raw_measurements_[1];

      //Initialize state
      ekf_.x_ << x, y, 0, 0;

    }

    //Update previous timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;


  //Set the process covariance matrix Q
  float t4 = (dt*dt*dt*dt)/4;
  float t3 = (dt*dt*dt)/2;
  float t2 = (dt*dt);
  
  //set the acceleration noise components
  int noise_ax = 9;
  int noise_ay = 9;

  ekf_.Q_ << (t4*noise_ax), 0, (t3*noise_ax), 0,
  0, (t4*noise_ay), 0, (t3*noise_ay),
  (t3*noise_ax), 0, (t2*noise_ax), 0,
  0, (t3*noise_ay), 0, (t2*noise_ay);


  if (dt < 0.001) {
    cout << "dt < 0.001" << endl;
    return;
  }

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    if ((fabs(ekf_.x_[0]) != 0 ) && (fabs(ekf_.x_[1]) != 0)) {


      ekf_.R_ = MatrixXd(3, 3);
      ekf_.R_ <<  R_radar_;

      ekf_.H_ = MatrixXd(3, 4);
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ << Hj_;

      //cout << endl << "Hj_ = " << ekf_.H_ << endl << endl;

      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }

  } else {
    // Laser updates
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ <<  R_laser_;
    
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ <<  H_laser_;

    //cout << endl << "H_ = " << ekf_.H_ << endl << endl;

    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
