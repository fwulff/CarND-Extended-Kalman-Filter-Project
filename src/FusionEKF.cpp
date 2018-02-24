#include "FusionEKF.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

  /**
  * Finish initializing the FusionEKF.
  * Udacity sample code
  */
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);

  //no correlation
  R_laser_ << 0.0225, 0.,
              0., 0.0225;

  R_radar_ = MatrixXd(3, 3);

  //no correlation
  R_radar_ << 0.09, 0., 0.,
              0., 0.0009, 0.,
              0., 0., 0.09;

  H_laser_ = MatrixXd(2, 4);

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ = MatrixXd(3, 4);
  Hj_ <<  1, 1, 0, 0,
          1, 1, 0, 0,
          1, 1, 1, 1;

  //state covariance matrix P
  P_ = MatrixXd(4, 4);
  P_ <<  1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  //the initial transition matrix F_
  F_ = MatrixXd(4, 4);
  F_ <<  1, 0, 1, 0,
         0, 1, 0, 1,
         0, 0, 1, 0,
         0, 0, 0, 1;

  //set the process covariance matrix Q
  Q_ = MatrixXd(4, 4);

  //create a 4D state vector, we don't know yet the values of the x state
  x_ = VectorXd(4);
  x_ << 1, 1, 1, 1;

  //set the acceleration noise components
  noise_ax_ = 10;
  noise_ay_ = 10;

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

    //first measurement
    previous_timestamp_ = 0;

    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * convert radar from polar to cartesian coordinates.
    */

    // first measurement
    cout << "EKF: " << endl;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "Kalman Filter Initialization with RADAR" << endl;

      //Init EKF
      ekf_.Init(x_, P_, F_, Hj_, R_radar_, Q_);

      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      double rho_ = measurement_pack.raw_measurements_[0];
      double phi_ = measurement_pack.raw_measurements_[1];
      double rho_dot_ = measurement_pack.raw_measurements_[2];

      //set the state with the initial location and zero velocity
      ekf_.x_ << rho_ * cos(phi_), rho_ * sin(phi_), rho_dot_ * cos(phi_), rho_dot_ * sin(phi_);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Kalman Filter Initialization with LIDAR" << endl;

      /**
      Initialize state.
      */

      //Init EKF
      ekf_.Init(x_, P_, F_, H_laser_, R_laser_, Q_);

      //set the state with the initial location and zero velocity directly
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    //set timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    //check if initialisation was successful
    // done initializing, no need to predict or update
    if ((ekf_.x_[0] != 0.) and (ekf_.x_[1] != 0.)) {
        is_initialized_ = true;
    }

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
     * Time is measured in seconds.
     * Update the process noise covariance matrix.
   */


  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  cout << "Timestamp " <<  measurement_pack.timestamp_ << endl;
  previous_timestamp_ = measurement_pack.timestamp_;
  cout << "Elapsed time " <<  dt << endl;

  // precompute values
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax_, 0, dt_3/2*noise_ax_, 0,
          0, dt_4/4*noise_ay_, 0, dt_3/2*noise_ay_,
          dt_3/2*noise_ax_, 0, dt_2*noise_ax_, 0,
          0, dt_3/2*noise_ay_, 0, dt_2*noise_ay_;

  //predict
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    /**
     * Radar update
     * measurement update
     * Use EKF - Jacobian for Hj
     */

    VectorXd z = ekf_.x_;
    Hj_ = tools.CalculateJacobian(z);

    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else {
    /**
     * Radar update
     * measurement update
     * Use EKF - Jacobian for Hj
     */

    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
