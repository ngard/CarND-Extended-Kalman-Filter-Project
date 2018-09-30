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
FusionEKF::FusionEKF(int argc, char** argv) {
  is_initialized_ = false;

  use_laser_ = true;
  use_radar_ = true;

  for (int argi=0; argi<argc; ++argi) {
    string arg(argv[argi]);
    if (arg=="-no-laser")
      use_laser_ = false;
    if (arg=="-no-radar")
      use_radar_ = false;
  }

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  P_ = MatrixXd(4,4);
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  P_ <<
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1000, 0,
    0, 0, 0, 1000;
  
  H_laser_ <<
    1, 0, 0, 0,
    0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
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

    if (use_radar_ && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];

      cerr << "ro:" << ro << endl;
      cerr << "phi:" << phi << endl;
      cerr << "ro_dot:" << ro_dot << endl;
      cerr << "px:" << ro * cos(phi) << endl;     // px
      cerr << "py:" << ro * sin(phi) << endl;     // py
      cerr << "vx:" << ro_dot * cos(phi) << endl; // vx
      cerr << "vy:" << ro_dot * sin(phi) << endl; // vy
      
      ekf_.x_ <<
	ro * cos(phi),     // px
	ro * sin(phi),     // py
	ro_dot * cos(phi), // vx
	ro_dot * sin(phi); // vy

      is_initialized_ = true;
    } else if (use_laser_ && measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ <<
	measurement_pack.raw_measurements_[0], // px
	measurement_pack.raw_measurements_[1], // py
	0,                                     // vx
	0;                                     // vy

      is_initialized_ = true;
    }

    // done initializing, no need to predict or update
    ekf_.P_ = P_;
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
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

  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1.0e6;
  previous_timestamp_ = measurement_pack.timestamp_;

  cerr << "dt:" << dt << endl;

  float dt2 = dt*dt, dt3 = dt2*dt, dt4 = dt2*dt2;

  // Laser updates
  Eigen::MatrixXd F = MatrixXd(4,4);;
  F <<
    1, 0,dt, 0,
    0, 1, 0,dt,
    0, 0, 1, 0,
    0, 0, 0, 1;
  
  float noise_ax = 9., noise_ay = 9.;
  Eigen::MatrixXd Q = MatrixXd(4,4);
  Q <<
    dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
    0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
    dt3/2*noise_ax, 0, dt2*noise_ax, 0,
    0, dt3/2*noise_ay, 0, dt2*noise_ay;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    MatrixXd Hj = tools.CalculateJacobian(ekf_.x_);
    ekf_.Init(ekf_.x_, ekf_.P_, F, Hj, R_radar_, Q);

    if (use_radar_) {
      cerr << "radar:" << endl;
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];
      cerr << "ro:" << ro << endl;
      cerr << "phi:" << phi << endl;
      cerr << "ro_dot:" << ro_dot << endl;
      cerr << "px:" << ro * cos(phi) << endl;     // px
      cerr << "py:" << ro * sin(phi) << endl;     // py
      cerr << "vx:" << ro_dot * cos(phi) << " vx*dt:" << ro_dot*cos(phi)*dt << endl; // vx
      cerr << "vy:" << ro_dot * sin(phi) << " vy*dt:" << ro_dot*sin(phi)*dt << endl; // vy
    }
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    ekf_.Init(ekf_.x_, ekf_.P_, F, H_laser_, R_laser_, Q);

    if (use_laser_) {
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      cerr << "laser:" << endl;
      cerr << "px:" << px << endl;     // px
      cerr << "py:" << py << endl;     // py
    }
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
  if (use_radar_ && measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else if (use_laser_ && measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  cerr << "update:" << endl;
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
