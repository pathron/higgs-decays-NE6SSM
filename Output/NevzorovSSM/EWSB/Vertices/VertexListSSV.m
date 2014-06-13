{{{Ah[{gt1}], hh[{gt2}], VZ[{lt3}]}, 
  {(conj[ZH[gt2, 1]]*(g2*Cos[ThetaW]*Cos[ThetaWp] + g1*Cos[ThetaWp]*
        Sin[ThetaW] - 6*g1p*Sin[ThetaWp])*ZA[gt1, 1] - 
     conj[ZH[gt2, 2]]*(g2*Cos[ThetaW]*Cos[ThetaWp] + 
       g1*Cos[ThetaWp]*Sin[ThetaW] + 4*g1p*Sin[ThetaWp])*ZA[gt1, 2] + 
     2*g1p*QS*Sin[ThetaWp]*(conj[ZH[gt2, 3]]*ZA[gt1, 3] - 
       conj[ZH[gt2, 4]]*ZA[gt1, 4]))/2, Mom[Ah[{gt1}], lt3] - 
    Mom[hh[{gt2}], lt3]}}, {{Ah[{gt1}], hh[{gt2}], VZp[{lt3}]}, 
  {(-(conj[ZH[gt2, 1]]*(6*g1p*Cos[ThetaWp] + 
        (g2*Cos[ThetaW] + g1*Sin[ThetaW])*Sin[ThetaWp])*ZA[gt1, 1]) - 
     conj[ZH[gt2, 2]]*(4*g1p*Cos[ThetaWp] - (g2*Cos[ThetaW] + g1*Sin[ThetaW])*
        Sin[ThetaWp])*ZA[gt1, 2] - 2*g1p*QS*Cos[ThetaWp]*
      (-(conj[ZH[gt2, 3]]*ZA[gt1, 3]) + conj[ZH[gt2, 4]]*ZA[gt1, 4]))/2, 
   Mom[Ah[{gt1}], lt3] - Mom[hh[{gt2}], lt3]}}, 
 {{Ah[{gt1}], Hpm[{gt2}], conj[VWm[{lt3}]]}, 
  {(g2*(conj[ZP[gt2, 1]]*ZA[gt1, 1] + conj[ZP[gt2, 2]]*ZA[gt1, 2]))/2, 
   Mom[Ah[{gt1}], lt3] - Mom[Hpm[{gt2}], lt3]}}, 
 {{Ah[{gt1}], conj[Hpm[{gt2}]], VWm[{lt3}]}, 
  {(g2*(ZA[gt1, 1]*ZP[gt2, 1] + ZA[gt1, 2]*ZP[gt2, 2]))/2, 
   Mom[Ah[{gt1}], lt3] - Mom[conj[Hpm[{gt2}]], lt3]}}, 
 {{hh[{gt1}], Hpm[{gt2}], conj[VWm[{lt3}]]}, 
  {(I/2)*g2*(conj[ZH[gt1, 1]]*conj[ZP[gt2, 1]] - conj[ZH[gt1, 2]]*
      conj[ZP[gt2, 2]]), Mom[hh[{gt1}], lt3] - Mom[Hpm[{gt2}], lt3]}}, 
 {{hh[{gt1}], conj[Hpm[{gt2}]], VWm[{lt3}]}, 
  {(-I/2)*g2*(conj[ZH[gt1, 1]]*ZP[gt2, 1] - conj[ZH[gt1, 2]]*ZP[gt2, 2]), 
   -Mom[conj[Hpm[{gt2}]], lt3] + Mom[hh[{gt1}], lt3]}}, 
 {{Hpm[{gt1}], conj[Hpm[{gt2}]], VP[{lt3}]}, 
  {(I/2)*Delta[gt1, gt2]*(g1*Cos[ThetaW] + g2*Sin[ThetaW]), 
   -Mom[conj[Hpm[{gt2}]], lt3] + Mom[Hpm[{gt1}], lt3]}}, 
 {{Hpm[{gt1}], conj[Hpm[{gt2}]], VZ[{lt3}]}, 
  {(I/2)*(conj[ZP[gt1, 1]]*(g2*Cos[ThetaW]*Cos[ThetaWp] - 
       g1*Cos[ThetaWp]*Sin[ThetaW] + 6*g1p*Sin[ThetaWp])*ZP[gt2, 1] + 
     conj[ZP[gt1, 2]]*(g2*Cos[ThetaW]*Cos[ThetaWp] - 
       g1*Cos[ThetaWp]*Sin[ThetaW] - 4*g1p*Sin[ThetaWp])*ZP[gt2, 2]), 
   -Mom[conj[Hpm[{gt2}]], lt3] + Mom[Hpm[{gt1}], lt3]}}, 
 {{Hpm[{gt1}], conj[Hpm[{gt2}]], VZp[{lt3}]}, 
  {(I/2)*(conj[ZP[gt1, 1]]*(6*g1p*Cos[ThetaWp] + 
       (-(g2*Cos[ThetaW]) + g1*Sin[ThetaW])*Sin[ThetaWp])*ZP[gt2, 1] + 
     conj[ZP[gt1, 2]]*(-4*g1p*Cos[ThetaWp] + 
       (-(g2*Cos[ThetaW]) + g1*Sin[ThetaW])*Sin[ThetaWp])*ZP[gt2, 2]), 
   -Mom[conj[Hpm[{gt2}]], lt3] + Mom[Hpm[{gt1}], lt3]}}, 
 {{Sd[{gt1, ct1}], conj[Sd[{gt2, ct2}]], VG[{ct3, lt3}]}, 
  {(-I/2)*g3*Delta[gt1, gt2]*Lam[ct3, ct2, ct1], 
   -Mom[conj[Sd[{gt2, ct2}]], lt3] + Mom[Sd[{gt1, ct1}], lt3]}}, 
 {{Sd[{gt1, ct1}], conj[Sd[{gt2, ct2}]], VP[{lt3}]}, 
  {(-I/6)*Delta[ct1, ct2]*((g1*Cos[ThetaW] - 3*g2*Sin[ThetaW])*
      sum[j1, 1, 3, conj[ZD[gt1, j1]]*ZD[gt2, j1]] - 
     2*g1*Cos[ThetaW]*sum[j1, 1, 3, conj[ZD[gt1, 3 + j1]]*ZD[gt2, 3 + j1]]), 
   -Mom[conj[Sd[{gt2, ct2}]], lt3] + Mom[Sd[{gt1, ct1}], lt3]}}, 
 {{Sd[{gt1, ct1}], conj[Sd[{gt2, ct2}]], VZ[{lt3}]}, 
  {(I/6)*Delta[ct1, ct2]*((3*g2*Cos[ThetaW]*Cos[ThetaWp] + 
       g1*Cos[ThetaWp]*Sin[ThetaW] - 6*g1p*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZD[gt1, j1]]*ZD[gt2, j1]] - 
     2*(g1*Cos[ThetaWp]*Sin[ThetaW] - 6*g1p*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZD[gt1, 3 + j1]]*ZD[gt2, 3 + j1]]), 
   -Mom[conj[Sd[{gt2, ct2}]], lt3] + Mom[Sd[{gt1, ct1}], lt3]}}, 
 {{Sd[{gt1, ct1}], conj[Sd[{gt2, ct2}]], VZp[{lt3}]}, 
  {(-I/6)*Delta[ct1, ct2]*
    ((6*g1p*Cos[ThetaWp] + (3*g2*Cos[ThetaW] + g1*Sin[ThetaW])*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZD[gt1, j1]]*ZD[gt2, j1]] - 
     2*(6*g1p*Cos[ThetaWp] + g1*Sin[ThetaW]*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZD[gt1, 3 + j1]]*ZD[gt2, 3 + j1]]), 
   -Mom[conj[Sd[{gt2, ct2}]], lt3] + Mom[Sd[{gt1, ct1}], lt3]}}, 
 {{Sd[{gt1, ct1}], conj[Su[{gt2, ct2}]], conj[VWm[{lt3}]]}, 
  {((-I)*g2*Delta[ct1, ct2]*sum[j1, 1, 3, conj[ZD[gt1, j1]]*ZU[gt2, j1]])/
    Sqrt[2], -Mom[conj[Su[{gt2, ct2}]], lt3] + Mom[Sd[{gt1, ct1}], lt3]}}, 
 {{SDX[{gt1, ct1}], conj[SDX[{gt2, ct2}]], VG[{ct3, lt3}]}, 
  {(-I/2)*g3*Delta[gt1, gt2]*Lam[ct3, ct2, ct1], 
   -Mom[conj[SDX[{gt2, ct2}]], lt3] + Mom[SDX[{gt1, ct1}], lt3]}}, 
 {{SDX[{gt1, ct1}], conj[SDX[{gt2, ct2}]], VP[{lt3}]}, 
  {(I/3)*g1*Cos[ThetaW]*Delta[ct1, ct2]*Delta[gt1, gt2], 
   -Mom[conj[SDX[{gt2, ct2}]], lt3] + Mom[SDX[{gt1, ct1}], lt3]}}, 
 {{SDX[{gt1, ct1}], conj[SDX[{gt2, ct2}]], VZ[{lt3}]}, 
  {(-I/3)*Delta[ct1, ct2]*((g1*Cos[ThetaWp]*Sin[ThetaW] - 6*g1p*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZDX[gt1, j1]]*ZDX[gt2, j1]] + 
     (g1*Cos[ThetaWp]*Sin[ThetaW] + 9*g1p*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZDX[gt1, 3 + j1]]*ZDX[gt2, 3 + j1]]), 
   -Mom[conj[SDX[{gt2, ct2}]], lt3] + Mom[SDX[{gt1, ct1}], lt3]}}, 
 {{SDX[{gt1, ct1}], conj[SDX[{gt2, ct2}]], VZp[{lt3}]}, 
  {(I/3)*Delta[ct1, ct2]*((6*g1p*Cos[ThetaWp] + g1*Sin[ThetaW]*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZDX[gt1, j1]]*ZDX[gt2, j1]] + 
     (-9*g1p*Cos[ThetaWp] + g1*Sin[ThetaW]*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZDX[gt1, 3 + j1]]*ZDX[gt2, 3 + j1]]), 
   -Mom[conj[SDX[{gt2, ct2}]], lt3] + Mom[SDX[{gt1, ct1}], lt3]}}, 
 {{Se[{gt1}], conj[Se[{gt2}]], VP[{lt3}]}, 
  {(I/2)*((g1*Cos[ThetaW] + g2*Sin[ThetaW])*sum[j1, 1, 3, 
       conj[ZE[gt1, j1]]*ZE[gt2, j1]] + 2*g1*Cos[ThetaW]*
      sum[j1, 1, 3, conj[ZE[gt1, 3 + j1]]*ZE[gt2, 3 + j1]]), 
   -Mom[conj[Se[{gt2}]], lt3] + Mom[Se[{gt1}], lt3]}}, 
 {{Se[{gt1}], conj[Se[{gt2}]], VZ[{lt3}]}, 
  {(I/2)*((g2*Cos[ThetaW]*Cos[ThetaWp] - g1*Cos[ThetaWp]*Sin[ThetaW] - 
       4*g1p*Sin[ThetaWp])*sum[j1, 1, 3, conj[ZE[gt1, j1]]*ZE[gt2, j1]] + 
     2*(-(g1*Cos[ThetaWp]*Sin[ThetaW]) + g1p*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZE[gt1, 3 + j1]]*ZE[gt2, 3 + j1]]), 
   -Mom[conj[Se[{gt2}]], lt3] + Mom[Se[{gt1}], lt3]}}, 
 {{Se[{gt1}], conj[Se[{gt2}]], VZp[{lt3}]}, 
  {(-I/2)*((4*g1p*Cos[ThetaWp] + (g2*Cos[ThetaW] - g1*Sin[ThetaW])*
        Sin[ThetaWp])*sum[j1, 1, 3, conj[ZE[gt1, j1]]*ZE[gt2, j1]] - 
     2*(g1p*Cos[ThetaWp] + g1*Sin[ThetaW]*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZE[gt1, 3 + j1]]*ZE[gt2, 3 + j1]]), 
   -Mom[conj[Se[{gt2}]], lt3] + Mom[Se[{gt1}], lt3]}}, 
 {{Se[{gt1}], conj[Sv[{gt2}]], conj[VWm[{lt3}]]}, 
  {((-I)*g2*sum[j1, 1, 3, conj[ZE[gt1, j1]]*ZV[gt2, j1]])/Sqrt[2], 
   -Mom[conj[Sv[{gt2}]], lt3] + Mom[Se[{gt1}], lt3]}}, 
 {{SHp0[{gt1}], conj[SHpp[{gt2}]], VWm[{lt3}]}, 
  {((-I)*g2*(conj[UHp0[gt1, 1]]*UHpp[gt2, 1] - conj[UHp0[gt1, 2]]*
       UHpp[gt2, 2]))/Sqrt[2], -Mom[conj[SHpp[{gt2}]], lt3] + 
    Mom[SHp0[{gt1}], lt3]}}, {{SHp0[{gt1}], conj[SHp0[{gt2}]], VZ[{lt3}]}, 
  {(-I/2)*Delta[gt1, gt2]*(g2*Cos[ThetaW]*Cos[ThetaWp] + 
     g1*Cos[ThetaWp]*Sin[ThetaW] + 4*g1p*Sin[ThetaWp]), 
   -Mom[conj[SHp0[{gt2}]], lt3] + Mom[SHp0[{gt1}], lt3]}}, 
 {{SHp0[{gt1}], conj[SHp0[{gt2}]], VZp[{lt3}]}, 
  {(-I/2)*Delta[gt1, gt2]*(4*g1p*Cos[ThetaWp] - 
     (g2*Cos[ThetaW] + g1*Sin[ThetaW])*Sin[ThetaWp]), 
   -Mom[conj[SHp0[{gt2}]], lt3] + Mom[SHp0[{gt1}], lt3]}}, 
 {{SHpp[{gt1}], conj[SHpp[{gt2}]], VP[{lt3}]}, 
  {(I/2)*Delta[gt1, gt2]*(g1*Cos[ThetaW] + g2*Sin[ThetaW]), 
   -Mom[conj[SHpp[{gt2}]], lt3] + Mom[SHpp[{gt1}], lt3]}}, 
 {{SHpp[{gt1}], conj[SHpp[{gt2}]], VZ[{lt3}]}, 
  {(I/2)*Delta[gt1, gt2]*(g2*Cos[ThetaW]*Cos[ThetaWp] - 
     g1*Cos[ThetaWp]*Sin[ThetaW] - 4*g1p*Sin[ThetaWp]), 
   -Mom[conj[SHpp[{gt2}]], lt3] + Mom[SHpp[{gt1}], lt3]}}, 
 {{SHpp[{gt1}], conj[SHpp[{gt2}]], VZp[{lt3}]}, 
  {(-I/2)*Delta[gt1, gt2]*(4*g1p*Cos[ThetaWp] + 
     (g2*Cos[ThetaW] - g1*Sin[ThetaW])*Sin[ThetaWp]), 
   -Mom[conj[SHpp[{gt2}]], lt3] + Mom[SHpp[{gt1}], lt3]}}, 
 {{SHpp[{gt1}], conj[SHp0[{gt2}]], conj[VWm[{lt3}]]}, 
  {((-I)*g2*(conj[UHpp[gt1, 1]]*UHp0[gt2, 1] - conj[UHpp[gt1, 2]]*
       UHp0[gt2, 2]))/Sqrt[2], -Mom[conj[SHp0[{gt2}]], lt3] + 
    Mom[SHpp[{gt1}], lt3]}}, {{Su[{gt1, ct1}], conj[Su[{gt2, ct2}]], 
   VG[{ct3, lt3}]}, {(-I/2)*g3*Delta[gt1, gt2]*Lam[ct3, ct2, ct1], 
   -Mom[conj[Su[{gt2, ct2}]], lt3] + Mom[Su[{gt1, ct1}], lt3]}}, 
 {{Su[{gt1, ct1}], conj[Su[{gt2, ct2}]], VP[{lt3}]}, 
  {(-I/6)*Delta[ct1, ct2]*((g1*Cos[ThetaW] + 3*g2*Sin[ThetaW])*
      sum[j1, 1, 3, conj[ZU[gt1, j1]]*ZU[gt2, j1]] + 
     4*g1*Cos[ThetaW]*sum[j1, 1, 3, conj[ZU[gt1, 3 + j1]]*ZU[gt2, 3 + j1]]), 
   -Mom[conj[Su[{gt2, ct2}]], lt3] + Mom[Su[{gt1, ct1}], lt3]}}, 
 {{Su[{gt1, ct1}], conj[Sd[{gt2, ct2}]], VWm[{lt3}]}, 
  {((-I)*g2*Delta[ct1, ct2]*sum[j1, 1, 3, conj[ZU[gt1, j1]]*ZD[gt2, j1]])/
    Sqrt[2], -Mom[conj[Sd[{gt2, ct2}]], lt3] + Mom[Su[{gt1, ct1}], lt3]}}, 
 {{Su[{gt1, ct1}], conj[Su[{gt2, ct2}]], VZ[{lt3}]}, 
  {(-I/6)*Delta[ct1, ct2]*((3*g2*Cos[ThetaW]*Cos[ThetaWp] - 
       g1*Cos[ThetaWp]*Sin[ThetaW] + 6*g1p*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZU[gt1, j1]]*ZU[gt2, j1]] - 
     2*(2*g1*Cos[ThetaWp]*Sin[ThetaW] + 3*g1p*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZU[gt1, 3 + j1]]*ZU[gt2, 3 + j1]]), 
   -Mom[conj[Su[{gt2, ct2}]], lt3] + Mom[Su[{gt1, ct1}], lt3]}}, 
 {{Su[{gt1, ct1}], conj[Su[{gt2, ct2}]], VZp[{lt3}]}, 
  {(-I/6)*Delta[ct1, ct2]*
    ((6*g1p*Cos[ThetaWp] + (-3*g2*Cos[ThetaW] + g1*Sin[ThetaW])*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZU[gt1, j1]]*ZU[gt2, j1]] + 
     2*(-3*g1p*Cos[ThetaWp] + 2*g1*Sin[ThetaW]*Sin[ThetaWp])*
      sum[j1, 1, 3, conj[ZU[gt1, 3 + j1]]*ZU[gt2, 3 + j1]]), 
   -Mom[conj[Su[{gt2, ct2}]], lt3] + Mom[Su[{gt1, ct1}], lt3]}}, 
 {{Sv[{gt1}], conj[Se[{gt2}]], VWm[{lt3}]}, 
  {((-I)*g2*sum[j1, 1, 3, conj[ZV[gt1, j1]]*ZE[gt2, j1]])/Sqrt[2], 
   -Mom[conj[Se[{gt2}]], lt3] + Mom[Sv[{gt1}], lt3]}}, 
 {{Sv[{gt1}], conj[Sv[{gt2}]], VZ[{lt3}]}, 
  {(-I/2)*Delta[gt1, gt2]*(g2*Cos[ThetaW]*Cos[ThetaWp] + 
     g1*Cos[ThetaWp]*Sin[ThetaW] + 4*g1p*Sin[ThetaWp]), 
   -Mom[conj[Sv[{gt2}]], lt3] + Mom[Sv[{gt1}], lt3]}}, 
 {{Sv[{gt1}], conj[Sv[{gt2}]], VZp[{lt3}]}, 
  {(-I/2)*Delta[gt1, gt2]*(4*g1p*Cos[ThetaWp] - 
     (g2*Cos[ThetaW] + g1*Sin[ThetaW])*Sin[ThetaWp]), 
   -Mom[conj[Sv[{gt2}]], lt3] + Mom[Sv[{gt1}], lt3]}}}
