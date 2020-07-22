#include "rot_pim_factor.h"
#include <iostream>

using namespace std;

namespace qcalib
{
    //------------------------------------------------------------------------------
    // Inner class PreintegratedMeasurements
    //------------------------------------------------------------------------------
    void PreintegratedRot3Measurements::print(const string &s) const
    {
        PreintegratedRotation::print(s);
        cout << "biasHat [" << biasHat_.transpose() << "]" << endl;
        cout << " PreintMeasCov [ " << preintMeasCov_ << " ]" << endl;
    }

    //------------------------------------------------------------------------------
    bool PreintegratedRot3Measurements::equals(
        const PreintegratedRot3Measurements &other, double tol) const
    {
        return PreintegratedRotation::equals(other, tol) &&
               equal_with_abs_tol(biasHat_, other.biasHat_, tol);
    }

    //------------------------------------------------------------------------------
    void PreintegratedRot3Measurements::resetIntegration()
    {
        PreintegratedRotation::resetIntegration();
        preintMeasCov_.setZero();
    }

    //------------------------------------------------------------------------------
    void PreintegratedRot3Measurements::integrateMeasurement(
        const Vector3 &measuredOmega, double deltaT)
    {

        Matrix3 D_incrR_integratedOmega, Fr;
        PreintegratedRotation::integrateMeasurement(measuredOmega,
                                                    biasHat_, deltaT, &D_incrR_integratedOmega, &Fr);

        // first order uncertainty propagation
        // the deltaT allows to pass from continuous time noise to discrete time noise
        preintMeasCov_ = Fr * preintMeasCov_ * Fr.transpose() + p().gyroscopeCovariance * deltaT;
    }

    //------------------------------------------------------------------------------
    Vector3 PreintegratedRot3Measurements::predict(const Vector3 &bias,
                                                   OptionalJacobian<3, 3> H) const
    {
        const Vector3 biasOmegaIncr = bias - biasHat_;
        const Rot3 biascorrected = biascorrectedDeltaRij(biasOmegaIncr, H);
        Matrix3 D_omega_biascorrected;
        const Vector3 omega = Rot3::Logmap(biascorrected, H ? &D_omega_biascorrected : 0);
        if (H)
            (*H) = D_omega_biascorrected * (*H);
        return omega;
    }
    //------------------------------------------------------------------------------
    Vector PreintegratedRot3Measurements::DeltaAngles(
        const Vector &msr_gyro_t, const double msr_dt,
        const Vector3 &delta_angles)
    {

        // Note: all delta terms refer to an IMU\sensor system at t0

        // Calculate the corrected measurements using the Bias object
        Vector body_t_omega_body = msr_gyro_t;

        Rot3 R_t_to_t0 = Rot3::Expmap(delta_angles);

        R_t_to_t0 = R_t_to_t0 * Rot3::Expmap(body_t_omega_body * msr_dt);
        return Rot3::Logmap(R_t_to_t0);
    }

    Vector Rot3Factor::evaluateError(const Rot3 &rot_i, const Vector3 &bias,
                                     boost::optional<Matrix &> HR,
                                     boost::optional<Matrix &> HB) const
    {
        Vector error;
        Matrix33 H1, H2, H3, H4, H5, H6, H7;
        Vector3 bias_inr = bias - _PIM_.biasHat();
        Rot3 rot_delta_new = _PIM_.biascorrectedDeltaRij(bias_inr, H3);
        error = Rot3::Logmap(q_w_i_.compose(rot_i.compose(rot_delta_new, H2, H4), boost::none, H1), H5) -
                Rot3::Logmap(q_w_j_.compose(rot_i, boost::none, H6), H7);

        if (HR)
        {
            *HR = H5 * H1 * H2 + H7 * H6;
        }

        if (HB)
        {
            *HB = H5 * H1 * H4 * -H3;
        }

        return error;
    }

} // namespace qcalib