#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include<math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
* Initializes Unscented Kalman filter
*/
UKF::UKF() {
	is_initialized_ = false;

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	n_x_ = 5;

	n_aug_ = 7;

	n_sigma_ = 2 * n_aug_ + 1;

	lambda_ = 3 - n_aug_;


	// initial state vector
	x_ = VectorXd(n_x_);
	x_.fill(0.0);

	// initial covariance matrix
	P_ = MatrixXd(n_x_, n_x_);
	P_ << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0,
		0, 0, 1, 0, 0,
		0, 0, 0, 1, 0,
		0, 0, 0, 0, 1;

	Xsig_pred_ = MatrixXd(n_x_, n_sigma_);

	time_us_ = 0;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 0.2;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.55;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	weights_ = VectorXd(n_sigma_);
	// set weights
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i<n_sigma_; i++) {  //2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda_);
		weights_(i) = weight;
	}

	///* the current NIS for radar
	NIS_radar_ = 0.0;

	///* the current NIS for laser
	NIS_laser_ = 0.0;

}

UKF::~UKF() {}

/**
* @param {MeasurementPackage} meas_package The latest measurement data of
* either radar or laser.
*/
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	


	/*****************************************************************************
	*  Initialization
	****************************************************************************/
	if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == false) {
		return;
	}
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == false) {
		return;
	}

	if (!is_initialized_) {

		time_us_ = meas_package.timestamp_;
		

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/

			double rho = meas_package.raw_measurements_[0];
			double phi = meas_package.raw_measurements_[1];
			double rho_dot = meas_package.raw_measurements_[2];

			double px = rho * cos(phi);
			double py = rho * sin(phi);

			x_ << px, py, 0, 0, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			double px = meas_package.raw_measurements_[0];
			double py = meas_package.raw_measurements_[1];

			x_ << px, py, 0, 0, 0;


		}

		// done initializing
		if (x_[0] != 0 || x_[1] != 0) {
			is_initialized_ = true;
		}
		
		return;
	}

	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		Prediction(delta_t);
		UpdateRadar(meas_package);
	}
	else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
		Prediction(delta_t);
		UpdateLidar(meas_package);
	}

}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Prediction(double delta_t) {

	/*****************************************************************************
	*  Generate Augmented Sigma Points
	****************************************************************************/
	
	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);
	

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

	//create augmented mean state
	x_aug.fill(0);
	x_aug.head(5) = x_;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	/*****************************************************************************
	*  Predict Sigma Points
	****************************************************************************/
	//predict sigma points
	for (int i = 0; i< n_sigma_; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p = 0.0;
		double py_p = 0.0;

		//avoid division by zero
		if (fabs(yawd) > 0.00001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}
	/*****************************************************************************
	*  Predict Mean and Covariance
	****************************************************************************/


	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P_.fill(0.0);

	for (int i = 0; i < n_sigma_; i++) {  //iterate over sigma points

										  // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)));
		
		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateLidar(MeasurementPackage meas_package) {

	int n_z_ = 2;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_, n_sigma_);

	Zsig.row(0) = Xsig_pred_.row(0);
	Zsig.row(1) = Xsig_pred_.row(1);




	//transform sigma points into measurement space
	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	z_pred.fill(0.0);
	z_pred = Zsig * weights_;

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);
	S.fill(0.0);
	for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points
										  //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}
	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_, n_z_);
	R << std_laspx_*std_laspx_, 0,
		0, std_radphi_*std_radphi_;
	S = S + R;

	MatrixXd Tc = MatrixXd(n_x_, n_z_);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points

										  //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {MeasurementPackage} meas_package
*/
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	int n_z_ = 3;
	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z_, n_sigma_);

	//transform sigma points into measurement space
	for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points

										  // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);
		double yawd = Xsig_pred_(4, i);

		double rho = sqrt(p_x*p_x + p_y*p_y);
		
		if (p_x < 0.000001)
			p_x = 0.000001;
		double phi = atan2(p_y, p_x);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;
		double r_dot = 0.0;

		if (rho > 0.01)
			r_dot = (p_x*v1 + p_y*v2) / rho;
		
		// measurement model
		Zsig(0, i) = rho;       //rho
		Zsig(1, i) = phi;     //phi
		Zsig(2, i) = r_dot;   //r_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z_);
	z_pred.fill(0.0);
	for (int i = 0; i < n_sigma_; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z_, n_z_);
	S.fill(0.0);
	for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points
										  //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z_, n_z_);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	S = S + R;


	MatrixXd Tc = MatrixXd(n_x_, n_z_);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < n_sigma_; i++) {  //2n+1 simga points

										  //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z = VectorXd(n_z_);
	z <<meas_package.raw_measurements_[0],   //rho in m
		meas_package.raw_measurements_[1],   //phi in rad
		meas_package.raw_measurements_[2];   //rho_dot in m/s
	VectorXd z_diff = z - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
