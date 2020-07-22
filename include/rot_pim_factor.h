#pragma once

/* GTSAM includes */
#include <gtsam/navigation/PreintegratedRotation.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <boost/make_shared.hpp>

namespace qcalib
{
    using namespace gtsam;

    class GTSAM_EXPORT PreintegratedRot3Measurements : public PreintegratedRotation
    {

    protected:
        gtsam::Vector3 biasHat_; ///< Angular rate bias values used during preintegration.
        Matrix3 preintMeasCov_;  ///< Covariance matrix of the preintegrated measurements (first-order propagation from *measurementCovariance*)

        friend class Rot3Factor;

    public:
        /// Default constructor, only for serialization and Cython wrapper
        PreintegratedRot3Measurements() {}

        /**
         *  Default constructor, initialize with no measurements
         *  @param bias Current estimate of acceleration and rotation rate biases
         */
        PreintegratedRot3Measurements(const boost::shared_ptr<PreintegratedRotationParams> &p,
                                      const Vector3 &biasHat) : PreintegratedRotation(p), biasHat_(biasHat)
        {
            resetIntegration();
        }

        const PreintegratedRotationParams &p() const { return *boost::static_pointer_cast<const PreintegratedRotationParams>(p_); }
        const Vector3 &biasHat() const { return biasHat_; }
        const Matrix3 &preintMeasCov() const { return preintMeasCov_; }

        /// print
        void print(const std::string &s = "Preintegrated Measurements: ") const;

        /// equals
        bool equals(const PreintegratedRot3Measurements &, double tol = 1e-9) const;

        /// Reset inetgrated quantities to zero
        void resetIntegration();

        /**
         * Add a single Gyroscope measurement to the preintegration.
         * @param measureOmedga Measured angular velocity (in body frame)
         * @param deltaT Time step
         */
        void integrateMeasurement(const Vector3 &measuredOmega, double deltaT);

        /// Predict bias-corrected incremental rotation
        /// TODO: The matrix Hbias is the derivative of predict? Unit-test?
        Vector3 predict(const Vector3 &bias, OptionalJacobian<3, 3> H = boost::none) const;

        // This function is only used for test purposes
        // (compare numerical derivatives wrt analytic ones)
        static Vector DeltaAngles(const Vector &msr_gyro_t, const double msr_dt,
                                  const Vector3 &delta_angles);

        /// @deprecated constructor
        PreintegratedRot3Measurements(const Vector3 &biasHat,
                                      const Matrix3 &measuredOmegaCovariance)
            : PreintegratedRotation(boost::make_shared<PreintegratedRotationParams>()),
              biasHat_(biasHat)
        {
            p_->gyroscopeCovariance = measuredOmegaCovariance;
            resetIntegration();
        }

    private:
        /** Serialization function */
        friend class boost::serialization::access;
        template <class ARCHIVE>
        void serialize(ARCHIVE &ar, const unsigned int /*version*/)
        {
            ar &BOOST_SERIALIZATION_BASE_OBJECT_NVP(PreintegratedRotation);
            ar &BOOST_SERIALIZATION_NVP(p_);
            ar &BOOST_SERIALIZATION_NVP(biasHat_);
        }
    };

    class GTSAM_EXPORT Rot3Factor : public NoiseModelFactor2<Rot3, Vector3>
    {

        typedef Rot3Factor This;
        typedef NoiseModelFactor2<Rot3, Vector3> Base;

        PreintegratedRot3Measurements _PIM_;

        Rot3Factor() {}

    public:
        Rot3Factor(Key rot_i, Key bias, gtsam::Rot3 q_w_i, gtsam::Rot3 q_w_j,
                   const PreintegratedRot3Measurements &preintegratedMeasurements)
            : Base(noiseModel::Gaussian::Covariance(preintegratedMeasurements.preintMeasCov_),
                   rot_i, bias),
              q_w_i_(q_w_i), q_w_j_(q_w_j), _PIM_(preintegratedMeasurements)
        {
        }

        virtual ~Rot3Factor()
        {
        }
        const PreintegratedRot3Measurements &preintegratedMeasurements() const
        {
            return _PIM_;
        }

        /// vector of errors
        Vector evaluateError(const Rot3 &rot_i, const Vector3 &bias,
                             boost::optional<Matrix &> HR = boost::none,
                             boost::optional<Matrix &> HB = boost::none) const;

    private:
        Rot3 q_w_i_;
        Rot3 q_w_j_;
    };
} // namespace qcalib
