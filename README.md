# Unscented Kalman Filter Project
Now it is time to implement an unscented Kalman filter using the CTRV motion model. We will be using the same bicycle simulation data set from the extended Kalman filter project. That way we can compare our results with the EKF project.

Remember that all Kalman filters have the same three steps:

- Initialization
- Prediction
- Update
A standard Kalman filter can only handle linear equations. Both the extended Kalman filter and the unscented Kalman filter allow us to use non-linear equations; the difference between EKF and UKF is how they handle non-linear equations. But the basics are the same: initialize, predict, update.

## Algorithm
- Your algorithm uses the first measurements to initialize the state vectors and covariance matrices.
- Upon receiving a measurement after the first, the algorithm predicts object position to the current timestep and then update the prediction using the new measurement.
- Algorithm sets up the appropriate matrices given the type of measurement and calls the correct measurement function for a given sensor type.

## Results
- Our algorithm is run against "obj_pose-laser-radar-synthetic-input.txt". We'll collect the positions that our algorithm outputs and compare them to ground truth data. Our px, py, vx, and vy, RMSE are less than or equal to the values [.09, .10, .40, .30].

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/_rHwX3TeDA4/0.jpg)](https://www.youtube.com/watch?v=_rHwX3TeDA4)

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`
    
