# Satellite State Estimation Project
This project simulates a satellite in low Earth orbit and compares the performance of the Extended and Unscented Kalman Filter using optical and range measurements. A Monte Carlo analyses is also conducted and estimation error plots with 3σ bounds to assess filter consistency. A high fidelity and low fidelity is also developed to analyze the effects of model mismatch with atmospheric drag and solar radiation pressure perterbations. A report is included that details the analysis and results over 100 Monte Carlo simulations. Code developed for various homework assignments in AAE 590FE at Purdue are also included.

## Known Issues:
There are still issues with the measurement limits in the meas_func function. An outline for the logic can be found, but there were issues in implementation where the covariance matrices would yield all "NaN" values and is currently not yet resolved.
