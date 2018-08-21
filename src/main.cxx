//=========================================================================
//
// Copyright 2018 Kitware, Inc.
// Author: Guilbert Pierre (spguilbert@gmail.com)
// Data: 03-27-2018
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//=========================================================================

// This program propose an automatic multisensor system calibration
//
// Let's suppose that you have a set of sensors (GPS / IMU, lidar, camera, ...)
// on a solid frame. The system sensors + solid frames (Solid) constitute a solid
// which means that for any couple of points X1, X2 belong to Solid dX1X2 / dt = 0
// at any time.
//
// hence, all sensors si have a constant position Tsi and a constant orientation
// Rsi in the solid frame coordinate system. We also suppose that you have the trajectory
// and the orientation over the time of all sensors relative to their own orientation / position
// coordinate system at time t0.
//
// This tool will provide you an M-estimation of the calibration of all sensors using
// their trajectory / orientation over the time. To do that, we use the solid constraint
// of all sensors trajectory and orientation over the time. these constraints will lead
// to a non linear least square cost function solved using a Levenberg-Marquardt algorithm
//
// In the end, the relative Orientation Rsi/sj and position Tsi/sj of the sensors
// will be estimated. The algorithm is proof to gaussian noise since it is a
// maximum likelihood estimator

// EIGEN
#include <Eigen/Dense>

// STD
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <random>

// VTK
#include <vtkQuaternion.h>

#define PI 3.14159265359


//-----------------------------------------------------------------------------
Eigen::Matrix<double, 3, 3> GetRotationMatrix(Eigen::Matrix<double, 3, 1> T);

//-----------------------------------------------------------------------------
void AutoCalib(std::string transformsFolder);

//-----------------------------------------------------------------------------
std::vector<Eigen::Matrix<double, 6, 1> > LoadTransforms1(std::string inputFile);

//-----------------------------------------------------------------------------
std::vector<Eigen::Matrix<double, 6, 1> > LoadTransforms2(std::string inputFile);

//-----------------------------------------------------------------------------
void GlobalMethod(std::vector<Eigen::Matrix<double, 6, 1> > X, std::vector<Eigen::Matrix<double, 6, 1> > Y);

//-----------------------------------------------------------------------------
Eigen::Matrix<double, 3, 3> LinearMethod(Eigen::MatrixXd S);

//-----------------------------------------------------------------------------
Eigen::Matrix<double, 3, 3> NonLinearMethod(Eigen::MatrixXd S);

//-----------------------------------------------------------------------------
void CreateSyntheticData(std::string outputFilename);

//-----------------------------------------------------------------------------
Eigen::Matrix<double, 3, 1> GetAngleFromRotationMatrix(Eigen::Matrix<double, 3, 3> R);

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cout << "Please indicate what you want to do:" << std::endl;
    std::cout << "- autocalib ${transforms_folder}" << std::endl;
    std::cout << "- createdata ${folder_to_store}" << std::endl;

    return EXIT_FAILURE;
  }

  std::string operation(argv[1]);

  if (operation == "autocalib")
  {
    AutoCalib(std::string(argv[2]));
  }
  else if (operation == "createdata")
  {
    CreateSyntheticData(std::string(argv[2]));
  }
  else
  {
    std::cout << "Supported operation are: convert and benchmark:" << std::endl;
    std::cout << "- autocalib ${kuka_transforms} ${velodyne_transforms}" << std::endl;
    std::cout << "- createdata ${folder_to_store}" << std::endl;

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
void AutoCalib(std::string transformsFolder)
{
  // Load the transforms
  std::stringstream kuka1, kuka2, lidar1, lidar2;
  kuka1 << transformsFolder << "/KukaTransforms1.csv";
  kuka2 << transformsFolder << "/KukaTransforms2.csv";
  lidar1 << transformsFolder << "/SlamTransforms1.csv";
  lidar2 << transformsFolder << "/SlamTransforms2.csv";

  std::vector<Eigen::Matrix<double, 6, 1> > sensor1_data1 = LoadTransforms1(kuka1.str().c_str());
  std::vector<Eigen::Matrix<double, 6, 1> > sensor2_data1 = LoadTransforms2(lidar1.str().c_str());

  std::vector<Eigen::Matrix<double, 6, 1> > sensor1_data2 = LoadTransforms1(kuka2.str().c_str());
  std::vector<Eigen::Matrix<double, 6, 1> > sensor2_data2 = LoadTransforms2(lidar2.str().c_str());

  std::vector<Eigen::Matrix<double, 6, 1> > Kuka;
  std::vector<Eigen::Matrix<double, 6, 1> > Velodyne;
  for (unsigned int k = 0; k < sensor1_data1.size(); ++k)
  {
    Kuka.push_back(sensor1_data1[k]);
    Velodyne.push_back(sensor2_data1[k]);
  }
  for (unsigned int k = 0; k < sensor1_data2.size(); ++k)
  {
    Kuka.push_back(sensor1_data2[k]);
    Velodyne.push_back(sensor2_data2[k]);
  }

  std::cout << "Kuka: " << Kuka.size() << " first: " << Kuka[0].transpose() << std::endl;
  std::cout << "Velodyne: " << Velodyne.size() << " first: " << Velodyne[1].transpose() << std::endl;

  GlobalMethod(Kuka, Velodyne);
}

//----------------------------------------------------------------------------
void GlobalMethod(std::vector<Eigen::Matrix<double, 6, 1> > X, std::vector<Eigen::Matrix<double, 6, 1> > Y)
{
  // Matrix of contraints between the orientations
  Eigen::MatrixXd S(9 * X.size(), 9);

  Eigen::Matrix<double, 3, 1> anglesX0;
  anglesX0 << X[0](0), X[0](1), X[0](2);
  Eigen::Matrix<double, 3, 3> B = GetRotationMatrix(anglesX0);
  Eigen::Matrix<double, 3, 3> Rp0 = B;
  Eigen::Matrix<double, 3, 1> Tp0;
  Tp0 << X[0](3), X[0](4), X[0](5);

  Eigen::Matrix<double, 3, 1> anglesCalib;
  anglesCalib << 0 / 180.0 * PI,
                 0 / 180.0 * PI,
                 -90.0 / 180.0 * PI;

  Eigen::Matrix<double, 3, 3> Rcalib = GetRotationMatrix(anglesCalib);

  Eigen::Matrix<double, 3, 3> I;
  I << 1, 0, 0,
       0, 1, 0,
       0, 0, 1;

  Eigen::Matrix<double, 9, 2> EquationMapping;
  EquationMapping << 0, 0,
                     0, 1,
                     0, 2,
                     1, 0,
                     1, 1,
                     1, 2,
                     2, 0,
                     2, 1,
                     2, 2;

  Eigen::Matrix<double, 3, 3> ReverseMapping;
  ReverseMapping << 0, 1, 2,
                    3, 4, 5,
                    6, 7, 8;

  Eigen::Matrix<double, 3, 9> PatternA;
  PatternA << 1, 0, 0, 1, 0, 0, 1, 0, 0,
              0, 1, 0, 0, 1, 0, 0, 1, 0,
              0, 0, 1, 0, 0, 1, 0, 0, 1;

  Eigen::Matrix<double, 3, 3> A, C;

  Eigen::Matrix<double, 3, 3> ErrorRealCalib = Eigen::Matrix<double, 3, 3>::Zero();

  // For all time sample, add the constraint to the system
  for (unsigned int k = 0; k < X.size(); ++k)
  {
    Eigen::Matrix<double, 3, 1> anglesX, anglesY;
    anglesX << X[k](0), X[k](1), X[k](2);
    anglesY << Y[k](0), Y[k](1), Y[k](2);

    A = GetRotationMatrix(anglesX);
    C = GetRotationMatrix(anglesY);

    // check:
    Eigen::Matrix<double, 3, 3> Test = A * Rcalib - B * Rcalib * C;

    ErrorRealCalib = ErrorRealCalib + A * Rcalib - B * Rcalib * C;

    for (unsigned int eq = 0; eq < 9; ++eq)
    {
      // eq correspond to the contraint of the (i, j) coordinate
      // of the Matirx. Hence, convert eq to i, j coords
      unsigned int i = EquationMapping(eq, 0);
      unsigned int j = EquationMapping(eq, 1);

      // equation corresponding to the unknow coordinates
      // i.e Rk,p
      for (unsigned int xcoord = 0; xcoord < 9; ++xcoord)
      {
        unsigned int kk = EquationMapping(xcoord, 0);
        unsigned int pp = EquationMapping(xcoord, 1);

        S(9 * k + eq, xcoord) = A(i, kk) * PatternA(j, xcoord) - B(i, kk) * C(pp, j);
      }
    }
  }

  // Now, we will solve the system using a levenberg marquardt algorithm
  Eigen::Matrix<double, 3, 3> Restimation = NonLinearMethod(S);
  
  Eigen::Matrix<double, 9, 1> V1;
  V1 << Restimation(0, 0), Restimation(0, 1), Restimation(0, 2),
        Restimation(1, 0), Restimation(1, 1), Restimation(1, 2),
        Restimation(2, 0), Restimation(2, 1), Restimation(2, 2);

  Eigen::MatrixXd residualErros = S * V1;

  Eigen::Matrix<double, 3, 3> ErrorSolution = Eigen::Matrix<double, 3, 3>::Zero();
  for (unsigned int k = 0; k < X.size(); ++k)
  {
    Eigen::Matrix<double, 3, 1> anglesX, anglesY;
    anglesX << X[k](0), X[k](1), X[k](2);
    anglesY << Y[k](0), Y[k](1), Y[k](2);

    A = GetRotationMatrix(anglesX);
    C = GetRotationMatrix(anglesY);
    // check:
    ErrorSolution = A * Restimation - B * Restimation * C;
  }

  // Now retrieve T
  Eigen::MatrixXd St(3 * X.size(), 6);
  Eigen::MatrixXd Yu(3 * X.size(), 1);
  for (unsigned int eq = 0; eq < X.size(); ++eq)
  {
    Eigen::Matrix<double, 3, 1> anglesX;
    anglesX << X[eq](0), X[eq](1), X[eq](2);
    Eigen::Matrix<double, 3, 3> Rp = GetRotationMatrix(anglesX);

    Eigen::Matrix<double, 3, 1> Tp, Tl;
    Tp << X[eq](3), X[eq](4), X[eq](5);
    Tl << Y[eq](3), Y[eq](4), Y[eq](5);

    Eigen::Matrix<double, 3, 3> K = Restimation.transpose() * Rp0.transpose() * Rp - Restimation.transpose();

    Eigen::Matrix<double, 3, 1> Yui = Tl + Restimation.transpose() * Rp0.transpose() * (Tp0 - Tp);
    Eigen::Matrix<double, 3, 3> B1 = K;
    Eigen::Matrix<double, 3, 3> A1 = K + Restimation.transpose() - Restimation.transpose() * Rp0.transpose() * Rp;

    for (unsigned int coord = 0; coord < 3; ++coord)
    {
      Yu(3 * eq + coord) = Yui(coord);

      St(3 * eq + coord, 0) = A1(coord, 0);
      St(3 * eq + coord, 1) = A1(coord, 1);
      St(3 * eq + coord, 2) = A1(coord, 2);

      St(3 * eq + coord, 3) = B1(coord, 0);
      St(3 * eq + coord, 4) = B1(coord, 1);
      St(3 * eq + coord, 5) = B1(coord, 2);
    }
  }

  Eigen::Matrix<double, 6, 1> Testimation = (St.transpose() * St).inverse() * St.transpose() * Yu;
  std::cout << "Estimated psotion params: " << Testimation.transpose() << std::endl;
  
  std::cout << "Error with real calib: " << std::endl << ErrorRealCalib / static_cast<double>(X.size()) << std::endl;
  std::cout << "Error with solution: " << std::endl << ErrorSolution / static_cast<double>(X.size()) << std::endl;
  std::cout << "mean residual error: " << residualErros.sum() / static_cast<double>(9 * X.size()) << std::endl;
}

//----------------------------------------------------------------------------
Eigen::Matrix<double, 3, 3> NonLinearMethod(Eigen::MatrixXd S)
{
  Eigen::Matrix<double, 3, 3> Restimation;
  unsigned int maxIteration = 100;

  Eigen::Matrix<double, 3, 1> angles;
  angles << 0, 0, 0;

  std::vector<double> errors;
  std::vector<double> jacobianNorm;
  std::vector<double> lambdaValues;

  double lambda = 0.1;

  for (unsigned int iteration = 0; iteration < maxIteration; ++iteration)
  {
    // rx
    double cA = std::cos(angles(0));
    double sA = std::sin(angles(0));
    // ry
    double cB = std::cos(angles(1));
    double sB = std::sin(angles(1));
    // rz
    double cC = std::cos(angles(2));
    double sC = std::sin(angles(2));

    // Get rotation matrix for this position
    Eigen::Matrix<double, 3, 3> R = GetRotationMatrix(angles);
    Eigen::Matrix<double, 9, 1> X;
    X << R(0, 0), R(0, 1), R(0, 2),
         R(1, 0), R(1, 1), R(1, 2),
         R(2, 0), R(2, 1), R(2, 2);

    // Compute jacobian for this position
    Eigen::MatrixXd Jacobian(S.rows(), 3);
    Eigen::MatrixXd JacobianR(9, 3);
    Eigen::MatrixXd JacobianM(1, 9);
    Eigen::MatrixXd residualJacobian(1, 3);
    Eigen::MatrixXd residualValues(S.rows(), 1);

    for (unsigned int eq = 0; eq < S.rows(); ++eq)
    {
      // dX1 / drot
      JacobianR(0, 0) = 0;
      JacobianR(0, 1) = -sB * cC;
      JacobianR(0, 2) = -sC * cB;
      // dX2 / drot
      JacobianR(1, 0) = sC * sA + cC * sB * cA;
      JacobianR(1, 1) = cC * cB * sA;
      JacobianR(1, 2) = -cC * cA - sC * sB * sA;
      // dX3 / drot
      JacobianR(2, 0) = sC * cA - cC * sB * sA;
      JacobianR(2, 1) = cC * cB * cA;
      JacobianR(2, 2) = cC * sA - sC * sB * cA;
      // dX4 / drot
      JacobianR(3, 0) = 0;
      JacobianR(3, 1) = -sC * sB;
      JacobianR(3, 2) = cC * cB;
      // dX5 / drot
      JacobianR(4, 0) = -cC * sA + sC * sB * cA;
      JacobianR(4, 1) = sC * cB * sA;
      JacobianR(4, 2) = -sC * sA + cC * sB * sA;
      // dX6 / drot
      JacobianR(5, 0) = -cC * cA - sC * sB * sA;
      JacobianR(5, 1) = sC * cB * cA;
      JacobianR(5, 2) = sC * sA + cC * sB * cA;
      // dX7 / drot
      JacobianR(6, 0) = 0;
      JacobianR(6, 1) = -cB;
      JacobianR(6, 2) = 0;
      // dX8 / drot
      JacobianR(7, 0) = cB * cA;
      JacobianR(7, 1) = -sB * sA;
      JacobianR(7, 2) = 0;
      // dX9 / drot
      JacobianR(8, 0) = -cB * sA;
      JacobianR(8, 1) = -sB * cA;
      JacobianR(8, 2) = 0;

      // Jacobian of the linear form
      JacobianM = S.row(eq);

      // Jacobian of the composition
      residualJacobian = JacobianM * JacobianR;

      // fill the current Jacobian of residual
      for (unsigned int k = 0; k < 3; ++k)
      {
        Jacobian(eq, k) = residualJacobian(0, k);
      }

      residualValues(eq, 0) = S.row(eq) * X;
    } // System ready

    // Compute current RMSE
    double RMSE = std::sqrt(1.0 / static_cast<double>(residualValues.rows()) * (residualValues.transpose() * residualValues)(0));
    errors.push_back(RMSE);

    // Compute the next step direction
    Eigen::MatrixXd Jt = Jacobian.transpose();
    Eigen::MatrixXd JtJ = Jt * Jacobian;
    Eigen::MatrixXd JtY = Jt * residualValues;
    Eigen::Matrix<double, 3, 3> diagJtJ;
    diagJtJ << JtJ(0, 0), 0, 0,
               0, JtJ(1, 1), 0,
               0, 0, JtJ(2, 2);

    jacobianNorm.push_back(Jacobian.norm());
    lambdaValues.push_back(lambda);

    // The next step of the L-M algorithm is computed by solving
    // (JtJ + lambda * diagJtJ) = Jt * Y. To avoid the computation
    // of the inverse of (JtJ + lambda * diagJtJ) we use a gauss-pivot
    // algorithm to solve the linear equation for this particular point
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(JtJ + lambda * diagJtJ);
    Eigen::Matrix<double, 3, 1> dAngles = dec.solve(JtY);

    // Check if the cost function has not increase
    // in the last iteration. If it does, we are too
    // away from the solution to use the Gauss-Newton
    // algorithm. Increase lambda to drift toward gradient descent
    Eigen::Matrix<double, 3, 1> candidateAngles = angles - dAngles;
    // Get rotation matrix for this position
    R = GetRotationMatrix(angles);
    X << R(0, 0), R(0, 1), R(0, 2),
         R(1, 0), R(1, 1), R(1, 2),
         R(2, 0), R(2, 1), R(2, 2);
    for (unsigned int eq = 0; eq < S.rows(); ++eq)
    {
      residualValues(eq, 0) = S.row(eq) * X;
    }
    // Compute current RMSE
    double newRMSE = std::sqrt(1.0 / static_cast<double>(residualValues.rows()) * (residualValues.transpose() * residualValues)(0));
    
    if (newRMSE > RMSE)
    {
      lambda = 1.1 * lambda;
    }
    else
    {
      angles = candidateAngles;
      lambda = 1 / 1.1 * lambda;
    }
  }

  std::cout << "Cost function goes from: " << errors[0] << " to: " << errors[errors.size() - 1] << std::endl;
  std::cout << "Jacobian norm goes from: " << jacobianNorm[0] << " to: " << jacobianNorm[jacobianNorm.size() - 1] << std::endl;
  std::cout << "lambda value goes from: " << lambdaValues[0] << " to: " << lambdaValues[lambdaValues.size() - 1] << std::endl;
  std::cout << "final solution: " << angles.transpose() / PI * 180.0 << std::endl;

  for (unsigned int k = 0; k < maxIteration; ++k)
  {
    //std::cout << "error at iteration: " << k << " is : " << errors[k] << std::endl;
  }

  Restimation = GetRotationMatrix(angles);
  return Restimation;
}

//----------------------------------------------------------------------------
std::vector<Eigen::Matrix<double, 6, 1> > LoadVelodyneTransforms(std::string inputFile)
{
  std::vector<Eigen::Matrix<double, 6, 1> > transforms(0);

  std::ifstream file;
  file.open(inputFile);

  if (!file.is_open())
  {
    std::cout << "can't open file: " << inputFile << std::endl;
    return transforms;
  }

  std::string line;
  std::getline(file, line);
  if (line != "time, X, Y, Z, Roll, Pitch, Yaw")
  {
    std::cout << "Transform file header is not expected, must be a problem" << std::endl;
    return transforms;
  }

  // get the transforms
  while(std::getline(file, line))
  {
    std::vector<std::string> params;
    boost::split(params, line, boost::is_any_of(","));

    if (params.size() != 7)
    {
      std::cout << "Problem parsing the file, wrong number of transforms parameters" << std::endl;
    }

    Eigen::Matrix<double, 6, 1> transform;
    transform(0) = std::atof(params[4].c_str());
    transform(1) = std::atof(params[5].c_str());
    transform(2) = std::atof(params[6].c_str());
    transform(3) = std::atof(params[1].c_str());
    transform(4) = std::atof(params[2].c_str());
    transform(5) = std::atof(params[3].c_str());
    transforms.push_back(transform);
  }
  
  return transforms;
}

//----------------------------------------------------------------------------
std::vector<Eigen::Matrix<double, 6, 1> > LoadKukaTransforms(std::string inputFile)
{
  std::vector<Eigen::Matrix<double, 6, 1> > transforms(0);

  std::ifstream file;
  file.open(inputFile);

  if (!file.is_open())
  {
    std::cout << "can't open file: " << inputFile << std::endl;
    return transforms;
  }

  std::string line;
  std::getline(file, line);
  if (line != "X,Y,Z,Roll,Pitch,Yaw")
  {
    std::cout << "Transform file header is not expected, must be a problem" << std::endl;
    return transforms;
  }

  // get the transforms
  while(std::getline(file, line))
  {
    std::vector<std::string> params;
    boost::split(params, line, boost::is_any_of(","));

    if (params.size() != 6)
    {
      std::cout << "Problem parsing the file, wrong number of transforms parameters" << std::endl;
    }

    Eigen::Matrix<double, 6, 1> transform;
    transform(0) = std::atof(params[3].c_str()) / 180.0 * PI;
    transform(1) = std::atof(params[4].c_str()) / 180.0 * PI;
    transform(2) = std::atof(params[5].c_str()) / 180.0 * PI;
    transform(3) = std::atof(params[0].c_str()) / 1000.0;
    transform(4) = std::atof(params[1].c_str()) / 1000.0;
    transform(5) = std::atof(params[2].c_str()) / 1000.0;
    transforms.push_back(transform);
  }
  
  return transforms;
}

//-----------------------------------------------------------------------------
Eigen::Matrix<double, 3, 3> GetRotationMatrix(Eigen::Matrix<double, 3, 1> T)
{
  Eigen::Matrix3d R(
        Eigen::AngleAxisd(T(2), Eigen::Vector3d::UnitZ())       /* rotation around Z-axis */
        * Eigen::AngleAxisd(T(1), Eigen::Vector3d::UnitY())     /* rotation around Y-axis */
        * Eigen::AngleAxisd(T(0), Eigen::Vector3d::UnitX()));   /* rotation around X-axis */

  return R;
}

//----------------------------------------------------------------------------
Eigen::Matrix<double, 3, 1> GetAngleFromRotationMatrix(Eigen::Matrix<double, 3, 3> R)
{
  Eigen::Matrix<double, 3, 1> angles;
  angles(0) = std::atan2(R(2, 1), R(2, 2));
  angles(1) = -std::asin(R(2, 0));
  angles(2) = std::atan2(R(1, 0), R(0, 0));

  return angles;
}

//----------------------------------------------------------------------------
void CreateSyntheticData(std::string outputFilename)
{
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, 1.0);

  double sigmaAngle = 2.5 / 180.0 * PI;
  double sigmaMeters = 0.04;

  // degrees of rotations
  double angle0 = 0;
  double angle1 = 360.0 / 180.0 * PI;
  
  // Params
  Eigen::Matrix<double, 3, 3> Rparam;
  Eigen::Matrix<double, 3, 1> Tsolid, Tpince, Tparam, anglesParam;
  anglesParam << 26.7 / 180.0 * PI, -12.89 / 180.0 * PI, 90.0 / 180.0 * PI;
  Rparam = GetRotationMatrix(anglesParam);
  Tsolid << 2.5, 1.0, 4.5;
  Tpince << 0.2, 0.05, 0.1;
  Tparam << 0.05, 0.07, 0.08;

  // Create first rotaion
  Eigen::Matrix<double, 3, 1> axe1;
  axe1 << 0.26, -0.54, 0.14;
  axe1.normalize();
  
  // Create second rotation
  Eigen::Matrix<double, 3, 1> axe2;
  axe2 << -0.75, -0.64, 0.04;
  axe2.normalize();

  Eigen::Matrix<double, 3, 3> Rpince01, Rpince02;

  std::stringstream kuka1;
  std::stringstream kuka2;
  std::stringstream lidar1;
  std::stringstream lidar2;

  kuka1 << outputFilename << "/KukaTransforms1.csv";
  kuka2 << outputFilename << "/KukaTransforms2.csv";
  lidar1 << outputFilename << "/SlamTransforms1.csv";
  lidar2 << outputFilename << "/SlamTransforms2.csv";

  std::ofstream fileKuka1, fileKuka2, fileLidar1, fileLidar2;
  fileKuka1.open(kuka1.str().c_str());
  fileKuka2.open(kuka2.str().c_str());
  fileLidar1.open(lidar1.str().c_str());
  fileLidar2.open(lidar2.str().c_str());

  if (!fileKuka1.is_open() || !fileKuka2.is_open() || !fileLidar1.is_open() || !fileLidar2.is_open())
  {
    return;
  }

  fileKuka1 << "X,Y,Z,Roll,Pitch,Yaw" << std::endl;
  fileKuka2 << "X,Y,Z,Roll,Pitch,Yaw" << std::endl;
  fileLidar1 << "time, X, Y, Z, Roll, Pitch, Yaw" << std::endl;
  fileLidar2 << "time, X, Y, Z, Roll, Pitch, Yaw" << std::endl;

  double Nsample = 150;
  for (unsigned int k = 0; k < Nsample; ++k)
  {
    double time = static_cast<double>(k) / static_cast<double>(Nsample - 1);
    double angle = (1 - time) * angle0 + time * angle1;

    // error position
    Eigen::Matrix<double, 3, 1> dX1, dX2;
    dX1 << sigmaMeters * distribution(generator), sigmaMeters * distribution(generator), sigmaMeters * distribution(generator);
    dX2 << sigmaMeters * distribution(generator), sigmaMeters * distribution(generator), sigmaMeters * distribution(generator);
    // error angle
    Eigen::Matrix<double, 3, 1> dAngles1, dAngles2;
    dAngles1 << sigmaAngle * distribution(generator), sigmaAngle * distribution(generator), sigmaAngle * distribution(generator);
    dAngles2 << sigmaAngle * distribution(generator), sigmaAngle * distribution(generator), sigmaAngle * distribution(generator);


    Eigen::Matrix<double, 3, 3> Rpince, Rlidar;
    Eigen::Matrix<double, 3, 1> Tp, Tlidar;

    // Rotation 1
    double ax1[3] = {axe1(0), axe1(1), axe1(2)};
    vtkQuaterniond quat1;
    quat1.SetRotationAngleAndAxis(angle, ax1);
    double A1[3][3];
    quat1.ToMatrix3x3(A1);
    Eigen::Matrix<double, 3, 3> R1;
    R1 << A1[0][0], A1[0][1], A1[0][2],
          A1[1][0], A1[1][1], A1[1][2],
          A1[2][0], A1[2][1], A1[2][2];
    

    // rotation 2
    double ax2[3] = {axe2(0), axe2(1), axe2(2)};
    vtkQuaterniond quat2;
    quat2.SetRotationAngleAndAxis(angle, ax2);
    double A2[3][3];
    quat2.ToMatrix3x3(A2);
    Eigen::Matrix<double, 3, 3> R2;
    R2 << A2[0][0], A2[0][1], A2[0][2],
          A2[1][0], A2[1][1], A2[1][2],
          A2[2][0], A2[2][1], A2[2][2];

    if (k == 0)
    {
      Rpince01 = R1;
      Rpince02 = R2;
    }

    // Transform 1
    Rpince = R1;
    Rlidar = Rparam.transpose() * Rpince01.transpose() * R1 * Rparam;
    Tp = Tsolid + R1 * Tpince;
    Tlidar = (Rparam.transpose() * Rpince01.transpose() * R1 - Rparam.transpose()) * (Tpince + Tparam);

    Eigen::Matrix<double, 3, 1> anglePince, angleLidar;
    anglePince = GetAngleFromRotationMatrix(Rpince);
    angleLidar = GetAngleFromRotationMatrix(Rlidar);

    anglePince = anglePince / PI * 180.0;
    Tp *= 1000.0;

    angleLidar += dAngles1;
    Tlidar += dX1;

    // save it
    fileKuka1 << Tp(0) << "," << Tp(1) << "," << Tp(2) << "," << anglePince(0) << "," << anglePince(1) << "," << anglePince(2) << std::endl;
    fileLidar1 << time << "," << Tlidar(0) << "," << Tlidar(1) << "," << Tlidar(2) << "," << angleLidar(0) << "," << angleLidar(1) << "," << angleLidar(2) << std::endl;

    // Transform 2
    Rpince = R2;
    Rlidar = Rparam.transpose() * Rpince02.transpose() * R2 * Rparam;
    Tp = Tsolid + R2 * Tpince;
    Tlidar = (Rparam.transpose() * Rpince02.transpose() * R2 - Rparam.transpose()) * (Tpince + Tparam);

    anglePince = GetAngleFromRotationMatrix(Rpince);
    angleLidar = GetAngleFromRotationMatrix(Rlidar);

    anglePince = anglePince / PI * 180.0;
    Tp *= 1000.0;

    angleLidar += dAngles2;
    Tlidar += dX2;

    // save it
    fileKuka2 << Tp(0) << "," << Tp(1) << "," << Tp(2) << "," << anglePince(0) << "," << anglePince(1) << "," << anglePince(2) << std::endl;
    fileLidar2 << time << "," << Tlidar(0) << "," << Tlidar(1) << "," << Tlidar(2) << "," << angleLidar(0) << "," << angleLidar(1) << "," << angleLidar(2) << std::endl;
  }

  fileKuka1.close();
  fileKuka2.close();
  fileLidar1.close();
  fileLidar2.close();
  return;
}