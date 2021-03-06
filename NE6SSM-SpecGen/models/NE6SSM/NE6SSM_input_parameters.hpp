// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Fri 20 Jun 2014 21:37:54

#ifndef NE6SSM_INPUT_PARAMETERS_H
#define NE6SSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct NE6SSM_input_parameters {
   double TanBeta;
   double vSInput;
   double muPrimeInput;
   double BmuPrimeInput;
   double QS;
   double KappaPrInput;
   double SigmaxInput;
   Eigen::Matrix<double,3,3> KappaInput;
   double LambdaxInput;
   double MuPhiInput;
   double XiFInput;
   double g1pInput;
   double vsbInput;
   double vphiInput;
   Eigen::Matrix<double,3,3> TYdInput;
   Eigen::Matrix<double,3,3> TYeInput;
   double TKappaPrInput;
   double TSigmaxInput;
   Eigen::Matrix<double,3,3> TKappaInput;
   double TLambdaxInput;
   Eigen::Matrix<double,3,3> TYuInput;
   double BMuPhiInput;
   double LXiFInput;
   Eigen::Matrix<double,3,3> mq2Input;
   Eigen::Matrix<double,3,3> ml2Input;
   Eigen::Matrix<double,3,3> md2Input;
   Eigen::Matrix<double,3,3> mu2Input;
   Eigen::Matrix<double,3,3> me2Input;
   Eigen::Matrix<double,3,3> mDx2Input;
   Eigen::Matrix<double,3,3> mDxbar2Input;
   double mHp2Input;
   double mHpbar2Input;
   double MassBInput;
   double MassWBInput;
   double MassGInput;
   double MassBpInput;

   NE6SSM_input_parameters()
      : TanBeta(0), vSInput(0), muPrimeInput(0), BmuPrimeInput(0), QS(0),
   KappaPrInput(0), SigmaxInput(0), KappaInput(Eigen::Matrix<double,3,3>::Zero(
   )), LambdaxInput(0), MuPhiInput(0), XiFInput(0), g1pInput(0), vsbInput(0),
   vphiInput(0), TYdInput(Eigen::Matrix<double,3,3>::Zero()), TYeInput(
   Eigen::Matrix<double,3,3>::Zero()), TKappaPrInput(0), TSigmaxInput(0),
   TKappaInput(Eigen::Matrix<double,3,3>::Zero()), TLambdaxInput(0), TYuInput(
   Eigen::Matrix<double,3,3>::Zero()), BMuPhiInput(0), LXiFInput(0), mq2Input(
   Eigen::Matrix<double,3,3>::Zero()), ml2Input(Eigen::Matrix<double,3,3>::Zero
   ()), md2Input(Eigen::Matrix<double,3,3>::Zero()), mu2Input(Eigen::Matrix<
   double,3,3>::Zero()), me2Input(Eigen::Matrix<double,3,3>::Zero()), mDx2Input
   (Eigen::Matrix<double,3,3>::Zero()), mDxbar2Input(Eigen::Matrix<double,3,3>
   ::Zero()), mHp2Input(0), mHpbar2Input(0), MassBInput(0), MassWBInput(0),
   MassGInput(0), MassBpInput(0)

   {}
};

} // namespace flexiblesusy

#endif
