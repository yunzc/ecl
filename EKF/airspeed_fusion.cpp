/****************************************************************************
 *
 *   Copyright (c) 2015 Estimation and Control Library (ECL). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name ECL nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file airspeed_fusion.cpp
 * airspeed fusion methods.
 *
 * @author Carl Olsson <carlolsson.co@gmail.com>
 * @author Roman Bast <bapstroman@gmail.com>
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */
#include "../ecl.h"
#include "ekf.h"
#include "mathlib.h"

void Ekf::fuseAirspeed()
{

	// Initialize variables
	float vn; // Velocity in north direction
	float ve; // Velocity in east direction
	float vd; // Velocity in downwards direction
	float vwn; // Wind speed in north direction
	float vwe; // Wind speed in east direction
	float v_tas_pred; // Predicted measurement
	float R_TAS = sq(math::constrain(_params.eas_noise, 0.5f, 5.0f) * math::constrain(_airspeed_sample_delayed.eas2tas, 0.9f,
			 10.0f)); // Variance for true airspeed measurement - (m/sec)^2
	float SH_TAS[3] = {}; // Varialbe used to optimise calculations of measurement jacobian
	float H_TAS[24] = {}; // Observation Jacobian
	float SK_TAS[2] = {}; // Varialbe used to optimise calculations of the Kalman gain vector
	float Kfusion[24] = {}; // Kalman gain vector

	// Copy required states to local variable names
	vn = _state.vel(0);
	ve = _state.vel(1);
	vd = _state.vel(2);
	vwn = _state.wind_vel(0);
	vwe = _state.wind_vel(1);

	// Calculate the predicted airspeed
	v_tas_pred = sqrtf((ve - vwe) * (ve - vwe) + (vn - vwn) * (vn - vwn) + vd * vd);

	// Perform fusion of True Airspeed measurement
	if (v_tas_pred > 1.0f) {
		// determine if we need the airspeed fusion to correct states other than wind
		bool update_wind_only = !_is_dead_reckoning;

		// Calculate the observation jacobian
		// intermediate variable from algebraic optimisation
		SH_TAS[0] = 1.0f/v_tas_pred;
		SH_TAS[1] = (SH_TAS[0]*(2.0f*ve - 2.0f*vwe))*0.5f;
		SH_TAS[2] = (SH_TAS[0]*(2.0f*vn - 2.0f*vwn))*0.5f;

		for (uint8_t i = 0; i < _k_num_states; i++) { H_TAS[i] = 0.0f; }

		H_TAS[4] = SH_TAS[2];
		H_TAS[5] = SH_TAS[1];
		H_TAS[6] = vd*SH_TAS[0];
		H_TAS[22] = -SH_TAS[2];
		H_TAS[23] = -SH_TAS[1];

		// We don't want to update the innovation variance if the calculation is ill conditioned
		float _airspeed_innov_var_temp = (R_TAS + SH_TAS[2]*(P[4][4]*SH_TAS[2] + P[5][4]*SH_TAS[1] - P[22][4]*SH_TAS[2] - P[23][4]*SH_TAS[1] + P[6][4]*vd*SH_TAS[0]) + SH_TAS[1]*(P[4][5]*SH_TAS[2] + P[5][5]*SH_TAS[1] - P[22][5]*SH_TAS[2] - P[23][5]*SH_TAS[1] + P[6][5]*vd*SH_TAS[0]) - SH_TAS[2]*(P[4][22]*SH_TAS[2] + P[5][22]*SH_TAS[1] - P[22][22]*SH_TAS[2] - P[23][22]*SH_TAS[1] + P[6][22]*vd*SH_TAS[0]) - SH_TAS[1]*(P[4][23]*SH_TAS[2] + P[5][23]*SH_TAS[1] - P[22][23]*SH_TAS[2] - P[23][23]*SH_TAS[1] + P[6][23]*vd*SH_TAS[0]) + vd*SH_TAS[0]*(P[4][6]*SH_TAS[2] + P[5][6]*SH_TAS[1] - P[22][6]*SH_TAS[2] - P[23][6]*SH_TAS[1] + P[6][6]*vd*SH_TAS[0]));

		if (_airspeed_innov_var_temp >= R_TAS) { // Check for badly conditioned calculation
			SK_TAS[0] = 1.0f / _airspeed_innov_var_temp;
			_fault_status.flags.bad_airspeed = false;

		} else { // Reset the estimator covarinace matrix
			_fault_status.flags.bad_airspeed = true;

			// if we are getting aiding from other sources, warn and reset the wind states and covariances only
			if (update_wind_only) {
				resetWindStates();
				resetWindCovariance();
				ECL_ERR("EKF airspeed fusion badly conditioned - wind covariance reset");

			} else {
				initialiseCovariance();
				_state.wind_vel.setZero();
				ECL_ERR("EKF airspeed fusion badly conditioned - full covariance reset");
			}

			return;
		}

		SK_TAS[1] = SH_TAS[1];

		if (update_wind_only) {
			// If we are getting aiding from other sources, then don't allow the airspeed measurements to affect the non-windspeed states
			for (unsigned row = 0; row <= 21; row++) {
				Kfusion[row] = 0.0f;
			}
		} else {
			// we have no other source of aiding, so use airspeed measurements to correct states
			Kfusion[0] = SK_TAS[0]*(P[0][4]*SH_TAS[2] - P[0][22]*SH_TAS[2] + P[0][5]*SK_TAS[1] - P[0][23]*SK_TAS[1] + P[0][6]*vd*SH_TAS[0]);
			Kfusion[1] = SK_TAS[0]*(P[1][4]*SH_TAS[2] - P[1][22]*SH_TAS[2] + P[1][5]*SK_TAS[1] - P[1][23]*SK_TAS[1] + P[1][6]*vd*SH_TAS[0]);
			Kfusion[2] = SK_TAS[0]*(P[2][4]*SH_TAS[2] - P[2][22]*SH_TAS[2] + P[2][5]*SK_TAS[1] - P[2][23]*SK_TAS[1] + P[2][6]*vd*SH_TAS[0]);
			Kfusion[3] = SK_TAS[0]*(P[3][4]*SH_TAS[2] - P[3][22]*SH_TAS[2] + P[3][5]*SK_TAS[1] - P[3][23]*SK_TAS[1] + P[3][6]*vd*SH_TAS[0]);
			Kfusion[4] = SK_TAS[0]*(P[4][4]*SH_TAS[2] - P[4][22]*SH_TAS[2] + P[4][5]*SK_TAS[1] - P[4][23]*SK_TAS[1] + P[4][6]*vd*SH_TAS[0]);
			Kfusion[5] = SK_TAS[0]*(P[5][4]*SH_TAS[2] - P[5][22]*SH_TAS[2] + P[5][5]*SK_TAS[1] - P[5][23]*SK_TAS[1] + P[5][6]*vd*SH_TAS[0]);
			Kfusion[6] = SK_TAS[0]*(P[6][4]*SH_TAS[2] - P[6][22]*SH_TAS[2] + P[6][5]*SK_TAS[1] - P[6][23]*SK_TAS[1] + P[6][6]*vd*SH_TAS[0]);
			Kfusion[7] = SK_TAS[0]*(P[7][4]*SH_TAS[2] - P[7][22]*SH_TAS[2] + P[7][5]*SK_TAS[1] - P[7][23]*SK_TAS[1] + P[7][6]*vd*SH_TAS[0]);
			Kfusion[8] = SK_TAS[0]*(P[8][4]*SH_TAS[2] - P[8][22]*SH_TAS[2] + P[8][5]*SK_TAS[1] - P[8][23]*SK_TAS[1] + P[8][6]*vd*SH_TAS[0]);
			Kfusion[9] = SK_TAS[0]*(P[9][4]*SH_TAS[2] - P[9][22]*SH_TAS[2] + P[9][5]*SK_TAS[1] - P[9][23]*SK_TAS[1] + P[9][6]*vd*SH_TAS[0]);
			Kfusion[10] = SK_TAS[0]*(P[10][4]*SH_TAS[2] - P[10][22]*SH_TAS[2] + P[10][5]*SK_TAS[1] - P[10][23]*SK_TAS[1] + P[10][6]*vd*SH_TAS[0]);
			Kfusion[11] = SK_TAS[0]*(P[11][4]*SH_TAS[2] - P[11][22]*SH_TAS[2] + P[11][5]*SK_TAS[1] - P[11][23]*SK_TAS[1] + P[11][6]*vd*SH_TAS[0]);
			Kfusion[12] = SK_TAS[0]*(P[12][4]*SH_TAS[2] - P[12][22]*SH_TAS[2] + P[12][5]*SK_TAS[1] - P[12][23]*SK_TAS[1] + P[12][6]*vd*SH_TAS[0]);
			Kfusion[13] = SK_TAS[0]*(P[13][4]*SH_TAS[2] - P[13][22]*SH_TAS[2] + P[13][5]*SK_TAS[1] - P[13][23]*SK_TAS[1] + P[13][6]*vd*SH_TAS[0]);
			Kfusion[14] = SK_TAS[0]*(P[14][4]*SH_TAS[2] - P[14][22]*SH_TAS[2] + P[14][5]*SK_TAS[1] - P[14][23]*SK_TAS[1] + P[14][6]*vd*SH_TAS[0]);
			Kfusion[15] = SK_TAS[0]*(P[15][4]*SH_TAS[2] - P[15][22]*SH_TAS[2] + P[15][5]*SK_TAS[1] - P[15][23]*SK_TAS[1] + P[15][6]*vd*SH_TAS[0]);
			Kfusion[16] = SK_TAS[0]*(P[16][4]*SH_TAS[2] - P[16][22]*SH_TAS[2] + P[16][5]*SK_TAS[1] - P[16][23]*SK_TAS[1] + P[16][6]*vd*SH_TAS[0]);
			Kfusion[17] = SK_TAS[0]*(P[17][4]*SH_TAS[2] - P[17][22]*SH_TAS[2] + P[17][5]*SK_TAS[1] - P[17][23]*SK_TAS[1] + P[17][6]*vd*SH_TAS[0]);
			Kfusion[18] = SK_TAS[0]*(P[18][4]*SH_TAS[2] - P[18][22]*SH_TAS[2] + P[18][5]*SK_TAS[1] - P[18][23]*SK_TAS[1] + P[18][6]*vd*SH_TAS[0]);
			Kfusion[19] = SK_TAS[0]*(P[19][4]*SH_TAS[2] - P[19][22]*SH_TAS[2] + P[19][5]*SK_TAS[1] - P[19][23]*SK_TAS[1] + P[19][6]*vd*SH_TAS[0]);
			Kfusion[20] = SK_TAS[0]*(P[20][4]*SH_TAS[2] - P[20][22]*SH_TAS[2] + P[20][5]*SK_TAS[1] - P[20][23]*SK_TAS[1] + P[20][6]*vd*SH_TAS[0]);
			Kfusion[21] = SK_TAS[0]*(P[21][4]*SH_TAS[2] - P[21][22]*SH_TAS[2] + P[21][5]*SK_TAS[1] - P[21][23]*SK_TAS[1] + P[21][6]*vd*SH_TAS[0]);
		}
		Kfusion[22] = SK_TAS[0]*(P[22][4]*SH_TAS[2] - P[22][22]*SH_TAS[2] + P[22][5]*SK_TAS[1] - P[22][23]*SK_TAS[1] + P[22][6]*vd*SH_TAS[0]);
		Kfusion[23] = SK_TAS[0]*(P[23][4]*SH_TAS[2] - P[23][22]*SH_TAS[2] + P[23][5]*SK_TAS[1] - P[23][23]*SK_TAS[1] + P[23][6]*vd*SH_TAS[0]);


		// Calculate measurement innovation
		_airspeed_innov = v_tas_pred -
				  _airspeed_sample_delayed.true_airspeed;

		// Calculate the innovation variance
		_airspeed_innov_var = 1.0f / SK_TAS[0];

		// Compute the ratio of innovation to gate size
		_tas_test_ratio = sq(_airspeed_innov) / (sq(fmaxf(_params.tas_innov_gate, 1.0f)) * _airspeed_innov_var);

		// If the innovation consistency check fails then don't fuse the sample and indicate bad airspeed health
		if (_tas_test_ratio > 1.0f) {
			_innov_check_fail_status.flags.reject_airspeed = true;
			return;
		} else {
			_innov_check_fail_status.flags.reject_airspeed = false;
		}

		// Airspeed measurement sample has passed check so record it
		_time_last_arsp_fuse = _time_last_imu;

		// apply covariance correction via P_new = (I -K*H)*P
		// first calculate expression for KHP
		// then calculate P - KHP
		float KHP[_k_num_states][_k_num_states];
		float KH[5];
		for (unsigned row = 0; row < _k_num_states; row++) {

			KH[0] = Kfusion[row] * H_TAS[4];
			KH[1] = Kfusion[row] * H_TAS[5];
			KH[2] = Kfusion[row] * H_TAS[6];
			KH[3] = Kfusion[row] * H_TAS[22];
			KH[4] = Kfusion[row] * H_TAS[23];

			for (unsigned column = 0; column < _k_num_states; column++) {
				float tmp = KH[0] * P[4][column];
				tmp += KH[1] * P[5][column];
				tmp += KH[2] * P[6][column];
				tmp += KH[3] * P[22][column];
				tmp += KH[4] * P[23][column];
				KHP[row][column] = tmp;
			}
		}

		// if the covariance correction will result in a negative variance, then
		// the covariance marix is unhealthy and must be corrected
		bool healthy = true;
		_fault_status.flags.bad_airspeed = false;
		for (int i = 0; i < _k_num_states; i++) {
			if (P[i][i] < KHP[i][i]) {
				// zero rows and columns
				zeroRows(P,i,i);
				zeroCols(P,i,i);

				//flag as unhealthy
				healthy = false;

				// update individual measurement health status
				_fault_status.flags.bad_airspeed = true;

			}
		}

		// only apply covariance and state corrrections if healthy
		if (healthy) {
			// apply the covariance corrections
			for (unsigned row = 0; row < _k_num_states; row++) {
				for (unsigned column = 0; column < _k_num_states; column++) {
					P[row][column] = P[row][column] - KHP[row][column];
				}
			}

			// correct the covariance marix for gross errors
			fixCovarianceErrors();

			// apply the state corrections
			fuse(Kfusion, _airspeed_innov);

		}
	}
}

void Ekf::fuseAirspeedBodyX()
{
	
	// Initialize variables
	float vn; // Velocity in north direction
	float ve; // Velocity in east direction
	float vd; // Velocity in downwards direction
	float vwn; // Wind speed in north direction
	float vwe; // Wind speed in east direction
	float v_tas_pred; // Predicted measurement
	float R_TAS = sq(math::constrain(_params.eas_noise, 0.5f, 5.0f) * math::constrain(_airspeed_sample_delayed.eas2tas, 0.9f,
			 10.0f)); // Variance for true airspeed measurement - (m/sec)^2
	float SH_TASX[2] = {}; // Varialbe used to optimise calculations of measurement jacobian
	float H_TASX[24] = {}; // Observation Jacobian
	float SK_TASX[8] = {}; // Varialbe used to optimise calculations of the Kalman gain vector
	float Kfusion[24] = {}; // Kalman gain vector

	// get latest estimated orientation
	float q0 = _state.quat_nominal(0);
	float q1 = _state.quat_nominal(1);
	float q2 = _state.quat_nominal(2);
	float q3 = _state.quat_nominal(3);

	// Copy required states to local variable names
	vn = _state.vel(0);
	ve = _state.vel(1);
	vd = _state.vel(2);
	vwn = _state.wind_vel(0);
	vwe = _state.wind_vel(1);

	// relative wind velocity in earth frame
	Vector3f rel_wind;
	rel_wind(0) = vn - vwn;
	rel_wind(1) = ve - vwe;
	rel_wind(2) = vd;

	Dcmf earth_to_body(_state.quat_nominal);
	earth_to_body = earth_to_body.transpose();

	// rotate into body axes
	rel_wind = earth_to_body * rel_wind;

	// The predicted airspeed measurement is the relative wind along body x
	v_tas_pred = rel_wind(0);

	// Perform fusion of True Airspeed measurement along body x
	if (v_tas_pred > 1.0f) {
		// determine if we need the airspeed fusion to correct states other than wind
		bool update_wind_only = !_is_dead_reckoning;

		// Calculate the observation jacobian
		// intermediate variable from algebraic optimisation
		SH_TASX[0] = vn - vwn;
		SH_TASX[1] = ve - vwe;

		for (uint8_t i = 0; i < _k_num_states; i++) { H_TASX[i] = 0.0f; }

		H_TASX[0] = 2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd;
		H_TASX[1] = 2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd;
		H_TASX[2] = 2*q1*SH_TASX[1] - 2*q2*SH_TASX[0] - 2*q0*vd;
		H_TASX[3] = 2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd;
		H_TASX[4] = sq(q0) + sq(q1) - sq(q2) - sq(q3);
		H_TASX[5] = 2*q0*q3 + 2*q1*q2;
		H_TASX[6] = 2*q1*q3 - 2*q0*q2;
		H_TASX[22] = sq(q2) - sq(q1) - sq(q0) + sq(q3);
		H_TASX[23] = - 2*q0*q3 - 2*q1*q2;

		// We don't want to update the innovation variance if the calculation is ill conditioned
		float _airspeed_innov_var_temp = (R_TAS + (sq(q0) + sq(q1) - sq(q2) - sq(q3))*(P[4][4]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][4]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][4]*(2*q0*q3 + 2*q1*q2) - P[6][4]*(2*q0*q2 - 2*q1*q3) - P[23][4]*(2*q0*q3 + 2*q1*q2) + P[0][4]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][4]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][4]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][4]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) - (sq(q0) + sq(q1) - sq(q2) - sq(q3))*(P[4][22]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][22]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][22]*(2*q0*q3 + 2*q1*q2) - P[6][22]*(2*q0*q2 - 2*q1*q3) - P[23][22]*(2*q0*q3 + 2*q1*q2) + P[0][22]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][22]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][22]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][22]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) + (2*q0*q3 + 2*q1*q2)*(P[4][5]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][5]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][5]*(2*q0*q3 + 2*q1*q2) - P[6][5]*(2*q0*q2 - 2*q1*q3) - P[23][5]*(2*q0*q3 + 2*q1*q2) + P[0][5]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][5]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][5]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][5]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) - (2*q0*q2 - 2*q1*q3)*(P[4][6]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][6]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][6]*(2*q0*q3 + 2*q1*q2) - P[6][6]*(2*q0*q2 - 2*q1*q3) - P[23][6]*(2*q0*q3 + 2*q1*q2) + P[0][6]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][6]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][6]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][6]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) - (2*q0*q3 + 2*q1*q2)*(P[4][23]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][23]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][23]*(2*q0*q3 + 2*q1*q2) - P[6][23]*(2*q0*q2 - 2*q1*q3) - P[23][23]*(2*q0*q3 + 2*q1*q2) + P[0][23]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][23]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][23]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][23]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) + (2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd)*(P[4][0]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][0]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][0]*(2*q0*q3 + 2*q1*q2) - P[6][0]*(2*q0*q2 - 2*q1*q3) - P[23][0]*(2*q0*q3 + 2*q1*q2) + P[0][0]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][0]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][0]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][0]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) + (2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd)*(P[4][1]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][1]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][1]*(2*q0*q3 + 2*q1*q2) - P[6][1]*(2*q0*q2 - 2*q1*q3) - P[23][1]*(2*q0*q3 + 2*q1*q2) + P[0][1]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][1]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][1]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][1]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) - (2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd)*(P[4][2]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][2]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][2]*(2*q0*q3 + 2*q1*q2) - P[6][2]*(2*q0*q2 - 2*q1*q3) - P[23][2]*(2*q0*q3 + 2*q1*q2) + P[0][2]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][2]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][2]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][2]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)) + (2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)*(P[4][3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - P[22][3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) + P[5][3]*(2*q0*q3 + 2*q1*q2) - P[6][3]*(2*q0*q2 - 2*q1*q3) - P[23][3]*(2*q0*q3 + 2*q1*q2) + P[0][3]*(2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd) + P[1][3]*(2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd) - P[2][3]*(2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd) + P[3][3]*(2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd)));

		if (_airspeed_innov_var_temp >= R_TAS) { // Check for badly conditioned calculation
			SK_TASX[0] = 1.0f / _airspeed_innov_var_temp;
			_fault_status.flags.bad_airspeed = false;

		} else { // Reset the estimator covarinace matrix
			_fault_status.flags.bad_airspeed = true;

			// if we are getting aiding from other sources, warn and reset the wind states and covariances only
			if (update_wind_only) {
				resetWindStates();
				resetWindCovariance();
				ECL_ERR("EKF airspeed fusion badly conditioned - wind covariance reset");

			} else {
				initialiseCovariance();
				_state.wind_vel.setZero();
				ECL_ERR("EKF airspeed fusion badly conditioned - full covariance reset");
			}

			return;
		}

		SK_TASX[1] = sq(q0) + sq(q1) - sq(q2) - sq(q3);
		SK_TASX[2] = 2*q0*q3 + 2*q1*q2;
		SK_TASX[3] = 2*q0*SH_TASX[1] - 2*q3*SH_TASX[0] + 2*q1*vd;
		SK_TASX[4] = 2*q2*SH_TASX[0] - 2*q1*SH_TASX[1] + 2*q0*vd;
		SK_TASX[5] = 2*q0*SH_TASX[0] + 2*q3*SH_TASX[1] - 2*q2*vd;
		SK_TASX[6] = 2*q1*SH_TASX[0] + 2*q2*SH_TASX[1] + 2*q3*vd;
		SK_TASX[7] = 2*q0*q2 - 2*q1*q3;

		if (update_wind_only) {
			// If we are getting aiding from other sources, then don't allow the airspeed measurements to affect the non-windspeed states
			for (unsigned row = 0; row <= 21; row++) {
				Kfusion[row] = 0.0f;
			}
		} else {
			// we have no other source of aiding, so use airspeed measurements to correct states
			Kfusion[0] = SK_TASX[0]*(P[0][0]*SK_TASX[5] + P[0][4]*SK_TASX[1] - P[0][2]*SK_TASX[4] + P[0][3]*SK_TASX[3] + P[0][1]*SK_TASX[6] + P[0][5]*SK_TASX[2] - P[0][6]*SK_TASX[7] - P[0][22]*SK_TASX[1] - P[0][23]*SK_TASX[2]);
			Kfusion[1] = SK_TASX[0]*(P[1][0]*SK_TASX[5] + P[1][4]*SK_TASX[1] - P[1][2]*SK_TASX[4] + P[1][3]*SK_TASX[3] + P[1][1]*SK_TASX[6] + P[1][5]*SK_TASX[2] - P[1][6]*SK_TASX[7] - P[1][22]*SK_TASX[1] - P[1][23]*SK_TASX[2]);
			Kfusion[2] = SK_TASX[0]*(P[2][0]*SK_TASX[5] + P[2][4]*SK_TASX[1] - P[2][2]*SK_TASX[4] + P[2][3]*SK_TASX[3] + P[2][1]*SK_TASX[6] + P[2][5]*SK_TASX[2] - P[2][6]*SK_TASX[7] - P[2][22]*SK_TASX[1] - P[2][23]*SK_TASX[2]);
			Kfusion[3] = SK_TASX[0]*(P[3][0]*SK_TASX[5] + P[3][4]*SK_TASX[1] - P[3][2]*SK_TASX[4] + P[3][3]*SK_TASX[3] + P[3][1]*SK_TASX[6] + P[3][5]*SK_TASX[2] - P[3][6]*SK_TASX[7] - P[3][22]*SK_TASX[1] - P[3][23]*SK_TASX[2]);
			Kfusion[4] = SK_TASX[0]*(P[4][0]*SK_TASX[5] + P[4][4]*SK_TASX[1] - P[4][2]*SK_TASX[4] + P[4][3]*SK_TASX[3] + P[4][1]*SK_TASX[6] + P[4][5]*SK_TASX[2] - P[4][6]*SK_TASX[7] - P[4][22]*SK_TASX[1] - P[4][23]*SK_TASX[2]);
			Kfusion[5] = SK_TASX[0]*(P[5][0]*SK_TASX[5] + P[5][4]*SK_TASX[1] - P[5][2]*SK_TASX[4] + P[5][3]*SK_TASX[3] + P[5][1]*SK_TASX[6] + P[5][5]*SK_TASX[2] - P[5][6]*SK_TASX[7] - P[5][22]*SK_TASX[1] - P[5][23]*SK_TASX[2]);
			Kfusion[6] = SK_TASX[0]*(P[6][0]*SK_TASX[5] + P[6][4]*SK_TASX[1] - P[6][2]*SK_TASX[4] + P[6][3]*SK_TASX[3] + P[6][1]*SK_TASX[6] + P[6][5]*SK_TASX[2] - P[6][6]*SK_TASX[7] - P[6][22]*SK_TASX[1] - P[6][23]*SK_TASX[2]);
			Kfusion[7] = SK_TASX[0]*(P[7][0]*SK_TASX[5] + P[7][4]*SK_TASX[1] - P[7][2]*SK_TASX[4] + P[7][3]*SK_TASX[3] + P[7][1]*SK_TASX[6] + P[7][5]*SK_TASX[2] - P[7][6]*SK_TASX[7] - P[7][22]*SK_TASX[1] - P[7][23]*SK_TASX[2]);
			Kfusion[8] = SK_TASX[0]*(P[8][0]*SK_TASX[5] + P[8][4]*SK_TASX[1] - P[8][2]*SK_TASX[4] + P[8][3]*SK_TASX[3] + P[8][1]*SK_TASX[6] + P[8][5]*SK_TASX[2] - P[8][6]*SK_TASX[7] - P[8][22]*SK_TASX[1] - P[8][23]*SK_TASX[2]);
			Kfusion[9] = SK_TASX[0]*(P[9][0]*SK_TASX[5] + P[9][4]*SK_TASX[1] - P[9][2]*SK_TASX[4] + P[9][3]*SK_TASX[3] + P[9][1]*SK_TASX[6] + P[9][5]*SK_TASX[2] - P[9][6]*SK_TASX[7] - P[9][22]*SK_TASX[1] - P[9][23]*SK_TASX[2]);
			Kfusion[10] = SK_TASX[0]*(P[10][0]*SK_TASX[5] + P[10][4]*SK_TASX[1] - P[10][2]*SK_TASX[4] + P[10][3]*SK_TASX[3] + P[10][1]*SK_TASX[6] + P[10][5]*SK_TASX[2] - P[10][6]*SK_TASX[7] - P[10][22]*SK_TASX[1] - P[10][23]*SK_TASX[2]);
			Kfusion[11] = SK_TASX[0]*(P[11][0]*SK_TASX[5] + P[11][4]*SK_TASX[1] - P[11][2]*SK_TASX[4] + P[11][3]*SK_TASX[3] + P[11][1]*SK_TASX[6] + P[11][5]*SK_TASX[2] - P[11][6]*SK_TASX[7] - P[11][22]*SK_TASX[1] - P[11][23]*SK_TASX[2]);
			Kfusion[12] = SK_TASX[0]*(P[12][0]*SK_TASX[5] + P[12][4]*SK_TASX[1] - P[12][2]*SK_TASX[4] + P[12][3]*SK_TASX[3] + P[12][1]*SK_TASX[6] + P[12][5]*SK_TASX[2] - P[12][6]*SK_TASX[7] - P[12][22]*SK_TASX[1] - P[12][23]*SK_TASX[2]);
			Kfusion[13] = SK_TASX[0]*(P[13][0]*SK_TASX[5] + P[13][4]*SK_TASX[1] - P[13][2]*SK_TASX[4] + P[13][3]*SK_TASX[3] + P[13][1]*SK_TASX[6] + P[13][5]*SK_TASX[2] - P[13][6]*SK_TASX[7] - P[13][22]*SK_TASX[1] - P[13][23]*SK_TASX[2]);
			Kfusion[14] = SK_TASX[0]*(P[14][0]*SK_TASX[5] + P[14][4]*SK_TASX[1] - P[14][2]*SK_TASX[4] + P[14][3]*SK_TASX[3] + P[14][1]*SK_TASX[6] + P[14][5]*SK_TASX[2] - P[14][6]*SK_TASX[7] - P[14][22]*SK_TASX[1] - P[14][23]*SK_TASX[2]);
			Kfusion[15] = SK_TASX[0]*(P[15][0]*SK_TASX[5] + P[15][4]*SK_TASX[1] - P[15][2]*SK_TASX[4] + P[15][3]*SK_TASX[3] + P[15][1]*SK_TASX[6] + P[15][5]*SK_TASX[2] - P[15][6]*SK_TASX[7] - P[15][22]*SK_TASX[1] - P[15][23]*SK_TASX[2]);
			Kfusion[16] = SK_TASX[0]*(P[16][0]*SK_TASX[5] + P[16][4]*SK_TASX[1] - P[16][2]*SK_TASX[4] + P[16][3]*SK_TASX[3] + P[16][1]*SK_TASX[6] + P[16][5]*SK_TASX[2] - P[16][6]*SK_TASX[7] - P[16][22]*SK_TASX[1] - P[16][23]*SK_TASX[2]);
			Kfusion[17] = SK_TASX[0]*(P[17][0]*SK_TASX[5] + P[17][4]*SK_TASX[1] - P[17][2]*SK_TASX[4] + P[17][3]*SK_TASX[3] + P[17][1]*SK_TASX[6] + P[17][5]*SK_TASX[2] - P[17][6]*SK_TASX[7] - P[17][22]*SK_TASX[1] - P[17][23]*SK_TASX[2]);
			Kfusion[18] = SK_TASX[0]*(P[18][0]*SK_TASX[5] + P[18][4]*SK_TASX[1] - P[18][2]*SK_TASX[4] + P[18][3]*SK_TASX[3] + P[18][1]*SK_TASX[6] + P[18][5]*SK_TASX[2] - P[18][6]*SK_TASX[7] - P[18][22]*SK_TASX[1] - P[18][23]*SK_TASX[2]);
			Kfusion[19] = SK_TASX[0]*(P[19][0]*SK_TASX[5] + P[19][4]*SK_TASX[1] - P[19][2]*SK_TASX[4] + P[19][3]*SK_TASX[3] + P[19][1]*SK_TASX[6] + P[19][5]*SK_TASX[2] - P[19][6]*SK_TASX[7] - P[19][22]*SK_TASX[1] - P[19][23]*SK_TASX[2]);
			Kfusion[20] = SK_TASX[0]*(P[20][0]*SK_TASX[5] + P[20][4]*SK_TASX[1] - P[20][2]*SK_TASX[4] + P[20][3]*SK_TASX[3] + P[20][1]*SK_TASX[6] + P[20][5]*SK_TASX[2] - P[20][6]*SK_TASX[7] - P[20][22]*SK_TASX[1] - P[20][23]*SK_TASX[2]);
			Kfusion[21] = SK_TASX[0]*(P[21][0]*SK_TASX[5] + P[21][4]*SK_TASX[1] - P[21][2]*SK_TASX[4] + P[21][3]*SK_TASX[3] + P[21][1]*SK_TASX[6] + P[21][5]*SK_TASX[2] - P[21][6]*SK_TASX[7] - P[21][22]*SK_TASX[1] - P[21][23]*SK_TASX[2]);

		}
		Kfusion[22] = SK_TASX[0]*(P[22][0]*SK_TASX[5] + P[22][4]*SK_TASX[1] - P[22][2]*SK_TASX[4] + P[22][3]*SK_TASX[3] + P[22][1]*SK_TASX[6] + P[22][5]*SK_TASX[2] - P[22][6]*SK_TASX[7] - P[22][22]*SK_TASX[1] - P[22][23]*SK_TASX[2]);
		Kfusion[23] = SK_TASX[0]*(P[23][0]*SK_TASX[5] + P[23][4]*SK_TASX[1] - P[23][2]*SK_TASX[4] + P[23][3]*SK_TASX[3] + P[23][1]*SK_TASX[6] + P[23][5]*SK_TASX[2] - P[23][6]*SK_TASX[7] - P[23][22]*SK_TASX[1] - P[23][23]*SK_TASX[2]);


		// Calculate measurement innovation
		_airspeed_innov = v_tas_pred -
				  _airspeed_sample_delayed.true_airspeed;

		// Calculate the innovation variance
		_airspeed_innov_var = 1.0f / SK_TASX[0];

		// Compute the ratio of innovation to gate size
		_tas_test_ratio = sq(_airspeed_innov) / (sq(fmaxf(_params.tas_innov_gate, 1.0f)) * _airspeed_innov_var);

		// If the innovation consistency check fails then don't fuse the sample and indicate bad airspeed health
		if (_tas_test_ratio > 1.0f) {
			_innov_check_fail_status.flags.reject_airspeed = true;
			return;
		} else {
			_innov_check_fail_status.flags.reject_airspeed = false;
		}

		// Airspeed measurement sample has passed check so record it
		_time_last_arsp_fuse = _time_last_imu;

		// apply covariance correction via P_new = (I -K*H)*P
		// first calculate expression for KHP
		// then calculate P - KHP
		float KHP[_k_num_states][_k_num_states];
		float KH[9];
		for (unsigned row = 0; row < _k_num_states; row++) {
			KH[0] = Kfusion[row] * H_TASX[0];
			KH[1] = Kfusion[row] * H_TASX[1];
			KH[2] = Kfusion[row] * H_TASX[2];
			KH[3] = Kfusion[row] * H_TASX[3];
			KH[4] = Kfusion[row] * H_TASX[4];
			KH[5] = Kfusion[row] * H_TASX[5];
			KH[6] = Kfusion[row] * H_TASX[6];
			KH[7] = Kfusion[row] * H_TASX[22];
			KH[8] = Kfusion[row] * H_TASX[23];

			for (unsigned column = 0; column < _k_num_states; column++) {
				float tmp = KH[0] * P[0][column];
				tmp += KH[1] * P[1][column];
				tmp += KH[2] * P[2][column];
				tmp += KH[3] * P[3][column];
				tmp += KH[4] * P[4][column];
				tmp += KH[5] * P[5][column];
				tmp += KH[6] * P[6][column];
				tmp += KH[7] * P[22][column];
				tmp += KH[8] * P[23][column];
				KHP[row][column] = tmp;
			}
		}

		// if the covariance correction will result in a negative variance, then
		// the covariance marix is unhealthy and must be corrected
		bool healthy = true;
		_fault_status.flags.bad_airspeed = false;
		for (int i = 0; i < _k_num_states; i++) {
			if (P[i][i] < KHP[i][i]) {
				// zero rows and columns
				zeroRows(P,i,i);
				zeroCols(P,i,i);

				//flag as unhealthy
				healthy = false;

				// update individual measurement health status
				_fault_status.flags.bad_airspeed = true;

			}
		}

		// only apply covariance and state corrrections if healthy
		if (healthy) {
			// apply the covariance corrections
			for (unsigned row = 0; row < _k_num_states; row++) {
				for (unsigned column = 0; column < _k_num_states; column++) {
					P[row][column] = P[row][column] - KHP[row][column];
				}
			}

			// correct the covariance marix for gross errors
			fixCovarianceErrors();

			// apply the state corrections
			fuse(Kfusion, _airspeed_innov);

		}
	}
}

void Ekf::get_wind_velocity(float *wind)
{
	wind[0] = _state.wind_vel(0);
	wind[1] = _state.wind_vel(1);
}

void Ekf::get_wind_velocity_var(float *wind_var)
{
	wind_var[0] = P[22][22];
	wind_var[1] = P[23][23];
}

void Ekf::get_true_airspeed(float *tas)
{
	float tempvar = sqrtf(sq(_state.vel(0) - _state.wind_vel(0)) + sq(_state.vel(1) - _state.wind_vel(1)) + sq(_state.vel(2)));
	memcpy(tas, &tempvar, sizeof(float));
}

/*
 * Reset the wind states using the current airspeed measurement, ground relative nav velocity, yaw angle and assumption of zero sideslip
*/
void Ekf::resetWindStates()
{
	// get euler yaw angle
	Eulerf euler321(_state.quat_nominal);
	float euler_yaw = euler321(2);

	if (_tas_data_ready && (_imu_sample_delayed.time_us - _airspeed_sample_delayed.time_us < (uint64_t)5e5)) {
		// estimate wind using zero sideslip assumption and airspeed measurement if airspeed available
		_state.wind_vel(0) = _state.vel(0) - _airspeed_sample_delayed.true_airspeed * cosf(euler_yaw);
		_state.wind_vel(1) = _state.vel(1) - _airspeed_sample_delayed.true_airspeed * sinf(euler_yaw);

	} else {
		// If we don't have an airspeed measurement, then assume the wind is zero
		_state.wind_vel(0) = 0.0f;
		_state.wind_vel(1) = 0.0f;
	}
}
