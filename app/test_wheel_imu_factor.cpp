#include "wheel_imu_factor.h"

#include <gtsam/nonlinear/Values.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/LevenbergMarquardtParams.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/slam/BetweenFactor.h>

#include <fstream>
#include <vector>

using namespace wio;

using namespace std;
using namespace gtsam;

using symbol_shorthand::B;
using symbol_shorthand::P;
using symbol_shorthand::R;
using symbol_shorthand::V;
using symbol_shorthand::X;

int main(void)
{
    std::cout << "###################### congratulation ######################" << std::endl;

    gtsam::LevenbergMarquardtParams parameters;
    parameters.absoluteErrorTol = 1e-8;
    parameters.relativeErrorTol = 1e-8;
    parameters.maxIterations = 500;
    parameters.verbosity = gtsam::NonlinearOptimizerParams::ERROR;
    parameters.verbosityLM = gtsam::LevenbergMarquardtParams::SUMMARY;

    const double kGravity = 9.81;
    auto params = PreintegrationParams::MakeSharedU(kGravity);
    params->setAccelerometerCovariance(I_3x3 * 1e-3);
    params->setGyroscopeCovariance(I_3x3 * 1e-4);
    params->setIntegrationCovariance(I_3x3 * 1e-8);
    params->setUse2ndOrderCoriolis(false);
    params->setOmegaCoriolis(Vector3(0, 0, 0));

    imuBias::ConstantBias prior_imu_bias; // assume zero initial bias
    Vector3 wheel_speed_noise(0.01, 0, 0);
    Matrix3 wheel_cov = wheel_speed_noise.asDiagonal();
    // Vector3 rot_xyz(0, 0.2, 0.3);
    Quaternion rot_b_o_q = Quaternion(0.996613, 0, 0.043584, 0.0697343); // 0 5deg 8deg
    // gtsam::Rot3 bRo = Rot3::RzRyRx(0, 0, 0);                             // 0 0 0
    gtsam::Rot3 bRo = Rot3(rot_b_o_q.toRotationMatrix());
    Vector3 bPo(0, 1, 0); // real 0 1 0

    PreintegratedImuWheelMeasurements i_w_pim(params, prior_imu_bias, bRo, wheel_cov);
    NonlinearFactorGraph graph;

    Vector3 position, acc, gyro, imu_velocity, wheel_speed;
    vector<Vector3> positions, accs, gyros, imu_velocitys, wheel_speeds;

    Pose3 pose_bak;
    Vector3 vel_bak;
    Pose3 pose_0;
    Vector3 vel_0;

    Quaternion imu_quat;
    vector<Quaternion> imu_quats;
    double dt = 0.01f;
    Values initial_value;
    ifstream imu_file;
    imu_file.open("imu_pose_noise.txt");
    ofstream est_file;
    est_file.open("est_pose.txt");

    double t, qw, qx, qy, qz, t0, t1, t2, v0, v1, v2, g0, g1, g2, a0, a1, a2, wv;
    uint64_t j = 0, key = 0;

    auto velnoise = noiseModel::Diagonal::Sigmas(Vector3(0.001, 0.001, 0.001));
    auto wheelexttonoise = noiseModel::Diagonal::Sigmas(Vector3(0.001, 0.001, 0.001));
    auto wheelextronoise = noiseModel::Diagonal::Sigmas(Vector3(0, 100, 100));
    auto posnoise = noiseModel::Diagonal::Sigmas(Vector6::Constant(0.001));
    auto biasnoise = noiseModel::Diagonal::Sigmas(Vector6::Constant(1e-6));

    imuBias::ConstantBias init_bias(Vector3(0.019, 0.019, 0.019), Vector3(0.015, 0.015, 0.015));
    imuBias::ConstantBias zero_bias(Vector3(0, 0, 0), Vector3(0, 0, 0));

    while (imu_file >> t >> qw >> qx >> qy >> qz >> t0 >> t1 >> t2 >> v0 >> v1 >> v2 >> g0 >> g1 >> g2 >> a0 >> a1 >> a2 >> wv)
    {
        acc << a0, a1, a2;
        gyro << g0, g1, g2;
        imu_velocity << v0, v1, v2;
        position << t0, t1, t2;
        imu_quat = Quaternion(qw, qx, qy, qz);
        wheel_speed << wv, 0, 0;

        accs.push_back(acc);
        gyros.push_back(gyro);
        imu_velocitys.push_back(imu_velocity);
        positions.push_back(position);
        imu_quats.push_back(imu_quat);
        wheel_speeds.push_back(wheel_speed);
    }
    imu_file.close();

    for (uint32_t index_t = 0; index_t < accs.size(); index_t++)
    {
        acc = accs[index_t];
        gyro = gyros[index_t];
        imu_velocity = imu_velocitys[index_t];
        position = positions[index_t];
        imu_quat = imu_quats[index_t];
        wheel_speed = wheel_speeds[index_t];

        if (j == 0)
        {
            pose_0 = Pose3(Rot3(imu_quat.toRotationMatrix()), position);
            vel_0 = imu_velocity;
        }

        // std::cout << "j: " << j << " ,acc: " << acc << " ,gyro: " << gyro << " ,imu_vel: " << imu_velocity << " ,pos: " << position << " ,ws: " << wheel_speed << std::endl;
        if (j % 10 == 0)
        {
            pose_bak = Pose3(Rot3(imu_quat.toRotationMatrix()), position);
            vel_bak = imu_velocity;
        }

        i_w_pim.integrateMeasurement(acc, gyro, wheel_speed, dt);
        if (j % 10 == 9)
        {
            if (key == 0)
            {
                graph.addPrior(X(key), pose_bak, posnoise);
                graph.addPrior(V(key), vel_bak, velnoise);
                graph.addPrior(B(key), init_bias, biasnoise);
            }
            initial_value.insert(X(key), pose_bak);
            initial_value.insert(V(key), vel_bak);
            initial_value.insert(B(key), init_bias);

            WheelImuFactor wheel_imu_factor(X(key), V(key), X(key + 1), V(key + 1), B(key), R(0), P(0), i_w_pim);

            Pose3 pose_end = Pose3(Rot3(imu_quat.toRotationMatrix()), position);
            Vector3 vel_end = imu_velocity;

            Matrix jac_posei, jac_posej, jac_bias;
            Matrix jac_veli, jac_velj, jac_r, jac_p;

            graph.add(wheel_imu_factor);
            if (j != 1999) //499
            {
                graph.add(BetweenFactor<imuBias::ConstantBias>(B(key), B(key + 1), zero_bias, biasnoise));
            }
            i_w_pim.resetIntegration();
            key++;
            if (key == 200)
            {
                break;
            };
        }
        j++;
        std::cout << "------------- j: " << j << " -------------key: " << key << std::endl;
    }

    Pose3 pose_end = Pose3(Rot3(imu_quat.toRotationMatrix()), position);
    Vector3 vel_end = imu_velocity;
    initial_value.insert(X(key), pose_end);
    initial_value.insert(V(key), vel_end);

    initial_value.insert(R(0), bRo);
    initial_value.insert(P(0), bPo);
    graph.addPrior(P(0), bPo, wheelexttonoise);
    graph.addPrior(R(0), bRo, wheelextronoise);

    std::cout << "###################### init optimizer ######################" << std::endl;
    LevenbergMarquardtOptimizer optimizer(graph, initial_value, parameters);

    std::cout << "###################### begin optimize ######################" << std::endl;
    Values result = optimizer.optimize();

    Rot3 optbRo = result.at<Rot3>(R(0));
    Vector3 optbPo = result.at<Vector3>(P(0));

    std::cout << "optbRo: \n"
              << Rot3::LocalCoordinates(optbRo) << std::endl;
    std::cout << "optbPo: \n"
              << optbPo << std::endl;
    for (uint32_t ikey = 0; ikey < key; ikey++)
    {
        // std::cout << "------------- pose check --------------\n";

        if (ikey != key)
        {
            imuBias::ConstantBias i_bias = result.at<imuBias::ConstantBias>(B(ikey));
            // std::cout << "ikey: " << ikey << "\n";
            // std::cout << "estimated i bias: \n"
            //           << i_bias << " \n";
        }

        Pose3 i_pose = result.at<Pose3>(X(ikey));
        Pose3 g_pose = Pose3(Rot3(imu_quats[ikey * 10].toRotationMatrix()), positions[ikey * 10]);
        // std::cout
        //     << "groud pose: \n"
        //     << Pose3::Logmap(g_pose) << "\n estimated pose: \n"
        //     << Pose3::Logmap(i_pose) << "\n";

        // est_file << Pose3::Logmap(g_pose).transpose() << " " << Pose3::Logmap(i_pose).transpose() << "\n";
        est_file << g_pose.translation().transpose() << " " << i_pose.translation().transpose() << "\n";
        // std::cout << "---------------------------------------\n";

        // Logmap vs LocalCoordinates
    }
    est_file.close();
    return 0;
}