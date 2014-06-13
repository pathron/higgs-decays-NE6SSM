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

// File generated at Fri 13 Jun 2014 15:36:48

#ifndef NE6SSM_INFO_H
#define NE6SSM_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace NE6SSM_info {
   enum Particles : unsigned {Glu, Fv, ChaP, VZ, VZp, Sd, Sv, Su, Se, SDX, hh,
      Ah, Hpm, Chi, Cha, Fe, Fd, Fu, FDX, SHp0, SHpp, ChiP, VG, VP, VWm,
      NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21,
      Yd22, Ye00, Ye01, Ye02, Ye10, Ye11, Ye12, Ye20, Ye21, Ye22, KappaPr, Sigmax,
      Kappa00, Kappa01, Kappa02, Kappa10, Kappa11, Kappa12, Kappa20, Kappa21,
      Kappa22, Lambdax, Yu00, Yu01, Yu02, Yu10, Yu11, Yu12, Yu20, Yu21, Yu22, MuPr
      , MuPhi, XiF, g1, g2, g3, g1p, vd, vu, vs, vsb, vphi, TYd00, TYd01, TYd02,
      TYd10, TYd11, TYd12, TYd20, TYd21, TYd22, TYe00, TYe01, TYe02, TYe10, TYe11,
      TYe12, TYe20, TYe21, TYe22, TKappaPr, TSigmax, TKappa00, TKappa01, TKappa02
      , TKappa10, TKappa11, TKappa12, TKappa20, TKappa21, TKappa22, TLambdax,
      TYu00, TYu01, TYu02, TYu10, TYu11, TYu12, TYu20, TYu21, TYu22, BMuPr, BMuPhi
      , LXiF, mq200, mq201, mq202, mq210, mq211, mq212, mq220, mq221, mq222, ml200
      , ml201, ml202, ml210, ml211, ml212, ml220, ml221, ml222, mHd2, mHu2, md200,
      md201, md202, md210, md211, md212, md220, md221, md222, mu200, mu201, mu202
      , mu210, mu211, mu212, mu220, mu221, mu222, me200, me201, me202, me210,
      me211, me212, me220, me221, me222, ms2, msbar2, mDx200, mDx201, mDx202,
      mDx210, mDx211, mDx212, mDx220, mDx221, mDx222, mDxbar200, mDxbar201,
      mDxbar202, mDxbar210, mDxbar211, mDxbar212, mDxbar220, mDxbar221, mDxbar222,
      mHp2, mHpbar2, mphi2, MassB, MassWB, MassG, MassBp, NUMBER_OF_PARAMETERS};

   extern const unsigned particle_multiplicities[NUMBER_OF_PARTICLES];
   extern const char* particle_names[NUMBER_OF_PARTICLES];
   extern const char* particle_latex_names[NUMBER_OF_PARTICLES];
   extern const char* parameter_names[NUMBER_OF_PARAMETERS];
   extern const char* model_name;
   extern const bool is_low_energy_model;

   void print(std::ostream&);
}

} // namespace flexiblesusy

#endif
