{{{aB, Ah[{gt2}], Ah[{gt3}]}, 
  {-(g1*(ZA[gt2, 1]*ZA[gt3, 1] - ZA[gt2, 2]*ZA[gt3, 2]))/2, 1}}, 
 {{aB, hh[{gt2}], hh[{gt3}]}, 
  {-(g1*(conj[ZH[gt2, 1]]*conj[ZH[gt3, 1]] - conj[ZH[gt2, 2]]*
        conj[ZH[gt3, 2]]))/2, 1}}, {{aB, Hpm[{gt2}], conj[Hpm[{gt3}]]}, 
  {-(g1*(conj[ZP[gt2, 1]]*ZP[gt3, 1] - conj[ZP[gt2, 2]]*ZP[gt3, 2]))/2, 1}}, 
 {{aB, Sd[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {(g1*Delta[ct2, ct3]*(sum[j1, 1, 3, conj[ZD[gt2, j1]]*ZD[gt3, j1]] + 
      2*sum[j1, 1, 3, conj[ZD[gt2, 3 + j1]]*ZD[gt3, 3 + j1]]))/6, 1}}, 
 {{aB, SDX[{gt2, ct2}], conj[SDX[{gt3, ct3}]]}, 
  {-(g1*Delta[ct2, ct3]*(sum[j1, 1, 3, conj[ZDX[gt2, j1]]*ZDX[gt3, j1]] - 
       sum[j1, 1, 3, conj[ZDX[gt2, 3 + j1]]*ZDX[gt3, 3 + j1]]))/3, 1}}, 
 {{aB, Se[{gt2}], conj[Se[{gt3}]]}, 
  {-(g1*(sum[j1, 1, 3, conj[ZE[gt2, j1]]*ZE[gt3, j1]] - 
       2*sum[j1, 1, 3, conj[ZE[gt2, 3 + j1]]*ZE[gt3, 3 + j1]]))/2, 1}}, 
 {{aB, SHp0[{gt2}], conj[SHp0[{gt3}]]}, 
  {-(g1*(conj[UHp0[gt2, 1]]*UHp0[gt3, 1] - conj[UHp0[gt2, 2]]*UHp0[gt3, 2]))/
    2, 1}}, {{aB, SHpp[{gt2}], conj[SHpp[{gt3}]]}, 
  {-(g1*(conj[UHpp[gt2, 1]]*UHpp[gt3, 1] - conj[UHpp[gt2, 2]]*UHpp[gt3, 2]))/
    2, 1}}, {{aB, Su[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {(g1*Delta[ct2, ct3]*(sum[j1, 1, 3, conj[ZU[gt2, j1]]*ZU[gt3, j1]] - 
      4*sum[j1, 1, 3, conj[ZU[gt2, 3 + j1]]*ZU[gt3, 3 + j1]]))/6, 1}}, 
 {{aB, Sv[{gt2}], conj[Sv[{gt3}]]}, {-(g1*Delta[gt2, gt3])/2, 1}}, 
 {{aBp, Ah[{gt2}], Ah[{gt3}]}, 
  {-(g1p*(3*ZA[gt2, 1]*ZA[gt3, 1] + 2*ZA[gt2, 2]*ZA[gt3, 2] - 
      QS*ZA[gt2, 3]*ZA[gt3, 3] + QS*ZA[gt2, 4]*ZA[gt3, 4])), 1}}, 
 {{aBp, hh[{gt2}], hh[{gt3}]}, 
  {-(g1p*(3*conj[ZH[gt2, 1]]*conj[ZH[gt3, 1]] + 2*conj[ZH[gt2, 2]]*
       conj[ZH[gt3, 2]] - QS*conj[ZH[gt2, 3]]*conj[ZH[gt3, 3]] + 
      QS*conj[ZH[gt2, 4]]*conj[ZH[gt3, 4]])), 1}}, 
 {{aBp, Hpm[{gt2}], conj[Hpm[{gt3}]]}, 
  {-(g1p*(3*conj[ZP[gt2, 1]]*ZP[gt3, 1] + 2*conj[ZP[gt2, 2]]*ZP[gt3, 2])), 
   1}}, {{aBp, Sd[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {g1p*Delta[ct2, ct3]*(sum[j1, 1, 3, conj[ZD[gt2, j1]]*ZD[gt3, j1]] + 
     2*sum[j1, 1, 3, conj[ZD[gt2, 3 + j1]]*ZD[gt3, 3 + j1]]), 1}}, 
 {{aBp, SDX[{gt2, ct2}], conj[SDX[{gt3, ct3}]]}, 
  {-(g1p*Delta[ct2, ct3]*(2*sum[j1, 1, 3, conj[ZDX[gt2, j1]]*ZDX[gt3, j1]] + 
      3*sum[j1, 1, 3, conj[ZDX[gt2, 3 + j1]]*ZDX[gt3, 3 + j1]])), 1}}, 
 {{aBp, Se[{gt2}], conj[Se[{gt3}]]}, 
  {g1p*(2*sum[j1, 1, 3, conj[ZE[gt2, j1]]*ZE[gt3, j1]] + 
     sum[j1, 1, 3, conj[ZE[gt2, 3 + j1]]*ZE[gt3, 3 + j1]]), 1}}, 
 {{aBp, SHp0[{gt2}], conj[SHp0[{gt3}]]}, 
  {2*g1p*(conj[UHp0[gt2, 1]]*UHp0[gt3, 1] - conj[UHp0[gt2, 2]]*UHp0[gt3, 2]), 
   1}}, {{aBp, SHpp[{gt2}], conj[SHpp[{gt3}]]}, 
  {2*g1p*(conj[UHpp[gt2, 1]]*UHpp[gt3, 1] - conj[UHpp[gt2, 2]]*UHpp[gt3, 2]), 
   1}}, {{aBp, Su[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {g1p*Delta[ct2, ct3]*Delta[gt2, gt3], 1}}, 
 {{aBp, Sv[{gt2}], conj[Sv[{gt3}]]}, {2*g1p*Delta[gt2, gt3], 1}}, 
 {{AdL[{gt1, ct1}], Ah[{gt2}], conj[Sd[{gt3, ct3}]]}, 
  {(I*Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yd[gt1, gt1]*ZA[gt2, 1]*
     ZD[gt3, 3 + gt1])/Sqrt[2], 1}}, 
 {{AdL[{gt1, ct1}], hh[{gt2}], conj[Sd[{gt3, ct3}]]}, 
  {(conj[ZH[gt2, 1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yd[gt1, gt1]*
     ZD[gt3, 3 + gt1])/Sqrt[2], 1}}, 
 {{AdL[{gt1, ct1}], conj[Hpm[{gt2}]], conj[Su[{gt3, ct3}]]}, 
  {-(Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yu[gt1, gt1]*ZP[gt2, 2]*
     ZU[gt3, 3 + gt1]), 1}}, {{AdR[{gt1, ct1}], Ah[{gt2}], 
   conj[Sd[{gt3, ct3}]]}, {((-I)*conj[Yd[gt1, gt1]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*ZA[gt2, 1]*ZD[gt3, gt1])/Sqrt[2], 1}}, 
 {{AdR[{gt1, ct1}], hh[{gt2}], conj[Sd[{gt3, ct3}]]}, 
  {(conj[Yd[gt1, gt1]]*conj[ZH[gt2, 1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*
     ZD[gt3, gt1])/Sqrt[2], 1}}, {{AdR[{gt1, ct1}], conj[Hpm[{gt2}]], 
   conj[Su[{gt3, ct3}]]}, {-(conj[Yd[gt1, gt1]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*ZP[gt2, 1]*ZU[gt3, gt1]), 1}}, 
 {{ADxbarR[{gt1, ct1}], Ah[{gt2}], conj[SDX[{gt3, ct3}]]}, 
  {((-I)*conj[\[Kappa][gt1, gt1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*
     ZA[gt2, 3]*ZDX[gt3, gt1])/Sqrt[2], 1}}, 
 {{ADxbarR[{gt1, ct1}], hh[{gt2}], conj[SDX[{gt3, ct3}]]}, 
  {(conj[ZH[gt2, 3]]*conj[\[Kappa][gt1, gt1]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*ZDX[gt3, gt1])/Sqrt[2], 1}}, 
 {{ADxL[{gt1, ct1}], Ah[{gt2}], conj[SDX[{gt3, ct3}]]}, 
  {(I*Delta[ct1, ct3]*ThetaStep[gt1, 3]*ZA[gt2, 3]*ZDX[gt3, 3 + gt1]*
     \[Kappa][gt1, gt1])/Sqrt[2], 1}}, 
 {{ADxL[{gt1, ct1}], hh[{gt2}], conj[SDX[{gt3, ct3}]]}, 
  {(conj[ZH[gt2, 3]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*ZDX[gt3, 3 + gt1]*
     \[Kappa][gt1, gt1])/Sqrt[2], 1}}, 
 {{AeL[{gt1}], Ah[{gt2}], conj[Se[{gt3}]]}, 
  {(I*ThetaStep[gt1, 3]*Ye[gt1, gt1]*ZA[gt2, 1]*ZE[gt3, 3 + gt1])/Sqrt[2], 
   1}}, {{AeL[{gt1}], hh[{gt2}], conj[Se[{gt3}]]}, 
  {(conj[ZH[gt2, 1]]*ThetaStep[gt1, 3]*Ye[gt1, gt1]*ZE[gt3, 3 + gt1])/
    Sqrt[2], 1}}, {{AeR[{gt1}], Ah[{gt2}], conj[Se[{gt3}]]}, 
  {((-I)*conj[Ye[gt1, gt1]]*ThetaStep[gt1, 3]*ZA[gt2, 1]*ZE[gt3, gt1])/
    Sqrt[2], 1}}, {{AeR[{gt1}], hh[{gt2}], conj[Se[{gt3}]]}, 
  {(conj[Ye[gt1, gt1]]*conj[ZH[gt2, 1]]*ThetaStep[gt1, 3]*ZE[gt3, gt1])/
    Sqrt[2], 1}}, {{AeR[{gt1}], conj[Hpm[{gt2}]], conj[Sv[{gt3}]]}, 
  {-(conj[Ye[gt1, gt1]]*ThetaStep[gt1, 3]*ZP[gt2, 1]*ZV[gt3, gt1]), 1}}, 
 {{aG[{ct1}], Sd[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {g3*(sum[j1, 1, 3, conj[ZD[gt2, j1]]*ZD[gt3, j1]] - 
     sum[j1, 1, 3, conj[ZD[gt2, 3 + j1]]*ZD[gt3, 3 + j1]]), 1}}, 
 {{aG[{ct1}], SDX[{gt2, ct2}], conj[SDX[{gt3, ct3}]]}, 
  {g3*(sum[j1, 1, 3, conj[ZDX[gt2, j1]]*ZDX[gt3, j1]] - 
     sum[j1, 1, 3, conj[ZDX[gt2, 3 + j1]]*ZDX[gt3, 3 + j1]]), 1}}, 
 {{aG[{ct1}], Su[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {g3*(sum[j1, 1, 3, conj[ZU[gt2, j1]]*ZU[gt3, j1]] - 
     sum[j1, 1, 3, conj[ZU[gt2, 3 + j1]]*ZU[gt3, 3 + j1]]), 1}}, 
 {{AHd0, Ah[{gt2}], Ah[{gt3}]}, 
  {(\[Lambda]*(ZA[gt2, 3]*ZA[gt3, 2] + ZA[gt2, 2]*ZA[gt3, 3]))/2, 1}}, 
 {{AHu0, Ah[{gt2}], Ah[{gt3}]}, 
  {(\[Lambda]*(ZA[gt2, 3]*ZA[gt3, 1] + ZA[gt2, 1]*ZA[gt3, 3]))/2, 1}}, 
 {{AphiR, Ah[{gt2}], Ah[{gt3}]}, 
  {(\[Sigma]*ZA[gt2, 4]*ZA[gt3, 3] + \[Sigma]*ZA[gt2, 3]*ZA[gt3, 4] - 
     2*\[Kappa]Pr*ZA[gt2, 5]*ZA[gt3, 5])/2, 1}}, 
 {{AsbarR, Ah[{gt2}], Ah[{gt3}]}, 
  {(\[Sigma]*(ZA[gt2, 5]*ZA[gt3, 3] + ZA[gt2, 3]*ZA[gt3, 5]))/2, 1}}, 
 {{AsR, Ah[{gt2}], Ah[{gt3}]}, 
  {(\[Lambda]*ZA[gt2, 2]*ZA[gt3, 1] + \[Lambda]*ZA[gt2, 1]*ZA[gt3, 2] + 
     \[Sigma]*ZA[gt2, 5]*ZA[gt3, 4] + \[Sigma]*ZA[gt2, 4]*ZA[gt3, 5])/2, 1}}, 
 {{aWB[{gt1}], Ah[{gt2}], Ah[{gt3}]}, 
  {(g2*(Sig[gt1, 1, 1]*ZA[gt2, 1]*ZA[gt3, 1] + Sig[gt1, 2, 2]*ZA[gt2, 2]*
       ZA[gt3, 2]))/2, 1}}, {{conj[AHd0], Ah[{gt2}], Ah[{gt3}]}, 
  {(conj[\[Lambda]]*(ZA[gt2, 3]*ZA[gt3, 2] + ZA[gt2, 2]*ZA[gt3, 3]))/2, 1}}, 
 {{conj[AHu0], Ah[{gt2}], Ah[{gt3}]}, 
  {(conj[\[Lambda]]*(ZA[gt2, 3]*ZA[gt3, 1] + ZA[gt2, 1]*ZA[gt3, 3]))/2, 1}}, 
 {{conj[AphiR], Ah[{gt2}], Ah[{gt3}]}, 
  {(conj[\[Sigma]]*(ZA[gt2, 4]*ZA[gt3, 3] + ZA[gt2, 3]*ZA[gt3, 4]) - 
     2*conj[\[Kappa]Pr]*ZA[gt2, 5]*ZA[gt3, 5])/2, 1}}, 
 {{conj[AsbarR], Ah[{gt2}], Ah[{gt3}]}, 
  {(conj[\[Sigma]]*(ZA[gt2, 5]*ZA[gt3, 3] + ZA[gt2, 3]*ZA[gt3, 5]))/2, 1}}, 
 {{conj[AsR], Ah[{gt2}], Ah[{gt3}]}, 
  {(conj[\[Lambda]]*(ZA[gt2, 2]*ZA[gt3, 1] + ZA[gt2, 1]*ZA[gt3, 2]) + 
     conj[\[Sigma]]*(ZA[gt2, 5]*ZA[gt3, 4] + ZA[gt2, 4]*ZA[gt3, 5]))/2, 1}}, 
 {{AHd0, Ah[{gt2}], hh[{gt3}]}, 
  {(-I/2)*\[Lambda]*(conj[ZH[gt3, 3]]*ZA[gt2, 2] + 
     conj[ZH[gt3, 2]]*ZA[gt2, 3]), 1}}, {{AHdm, Ah[{gt2}], conj[Hpm[{gt3}]]}, 
  {(I*\[Lambda]*ZA[gt2, 3]*ZP[gt3, 2])/Sqrt[2], 1}}, 
 {{AHu0, Ah[{gt2}], hh[{gt3}]}, 
  {(-I/2)*\[Lambda]*(conj[ZH[gt3, 3]]*ZA[gt2, 1] + 
     conj[ZH[gt3, 1]]*ZA[gt2, 3]), 1}}, {{AHup, Ah[{gt2}], Hpm[{gt3}]}, 
  {(I*\[Lambda]*conj[ZP[gt3, 1]]*ZA[gt2, 3])/Sqrt[2], 1}}, 
 {{AphiR, Ah[{gt2}], hh[{gt3}]}, 
  {I*(-(\[Sigma]*conj[ZH[gt3, 4]]*ZA[gt2, 3])/2 - 
     (\[Sigma]*conj[ZH[gt3, 3]]*ZA[gt2, 4])/2 + \[Kappa]Pr*conj[ZH[gt3, 5]]*
      ZA[gt2, 5]), 1}}, {{AsbarR, Ah[{gt2}], hh[{gt3}]}, 
  {(-I/2)*\[Sigma]*(conj[ZH[gt3, 5]]*ZA[gt2, 3] + 
     conj[ZH[gt3, 3]]*ZA[gt2, 5]), 1}}, {{AsR, Ah[{gt2}], hh[{gt3}]}, 
  {(I/2)*(-(\[Lambda]*conj[ZH[gt3, 2]]*ZA[gt2, 1]) - 
     \[Lambda]*conj[ZH[gt3, 1]]*ZA[gt2, 2] - 
     \[Sigma]*(conj[ZH[gt3, 5]]*ZA[gt2, 4] + conj[ZH[gt3, 4]]*ZA[gt2, 5])), 
   1}}, {{AuL[{gt1, ct1}], Ah[{gt2}], conj[Su[{gt3, ct3}]]}, 
  {(I*Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yu[gt1, gt1]*ZA[gt2, 2]*
     ZU[gt3, 3 + gt1])/Sqrt[2], 1}}, 
 {{AuR[{gt1, ct1}], Ah[{gt2}], conj[Su[{gt3, ct3}]]}, 
  {((-I)*conj[Yu[gt1, gt1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*ZA[gt2, 2]*
     ZU[gt3, gt1])/Sqrt[2], 1}}, {{aWB[{gt1}], Ah[{gt2}], Hpm[{gt3}]}, 
  {((I/2)*g2*Sig[gt1, 1, 2]*(-(conj[ZP[gt3, 1]]*ZA[gt2, 1]) + 
      conj[ZP[gt3, 2]]*ZA[gt2, 2]))/Sqrt[2], 1}}, 
 {{aWB[{gt1}], Ah[{gt2}], conj[Hpm[{gt3}]]}, 
  {((I/2)*g2*Sig[gt1, 2, 1]*(ZA[gt2, 1]*ZP[gt3, 1] - ZA[gt2, 2]*ZP[gt3, 2]))/
    Sqrt[2], 1}}, {{conj[AHd0], Ah[{gt2}], hh[{gt3}]}, 
  {(I/2)*conj[\[Lambda]]*(conj[ZH[gt3, 3]]*ZA[gt2, 2] + 
     conj[ZH[gt3, 2]]*ZA[gt2, 3]), 1}}, {{conj[AHu0], Ah[{gt2}], hh[{gt3}]}, 
  {(I/2)*conj[\[Lambda]]*(conj[ZH[gt3, 3]]*ZA[gt2, 1] + 
     conj[ZH[gt3, 1]]*ZA[gt2, 3]), 1}}, {{conj[AphiR], Ah[{gt2}], hh[{gt3}]}, 
  {(I/2)*(conj[\[Sigma]]*(conj[ZH[gt3, 4]]*ZA[gt2, 3] + 
       conj[ZH[gt3, 3]]*ZA[gt2, 4]) - 2*conj[\[Kappa]Pr]*conj[ZH[gt3, 5]]*
      ZA[gt2, 5]), 1}}, {{conj[AsbarR], Ah[{gt2}], hh[{gt3}]}, 
  {(I/2)*conj[\[Sigma]]*(conj[ZH[gt3, 5]]*ZA[gt2, 3] + 
     conj[ZH[gt3, 3]]*ZA[gt2, 5]), 1}}, {{conj[AsR], Ah[{gt2}], hh[{gt3}]}, 
  {(I/2)*(conj[\[Lambda]]*(conj[ZH[gt3, 2]]*ZA[gt2, 1] + 
       conj[ZH[gt3, 1]]*ZA[gt2, 2]) + conj[\[Sigma]]*
      (conj[ZH[gt3, 5]]*ZA[gt2, 4] + conj[ZH[gt3, 4]]*ZA[gt2, 5])), 1}}, 
 {{conj[AHdm], Ah[{gt2}], Hpm[{gt3}]}, 
  {((-I)*conj[\[Lambda]]*conj[ZP[gt3, 2]]*ZA[gt2, 3])/Sqrt[2], 1}}, 
 {{conj[AdL[{gt1, ct1}]], Ah[{gt2}], Sd[{gt3, ct3}]}, 
  {((-I)*conj[Yd[gt1, gt1]]*conj[ZD[gt3, 3 + gt1]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*ZA[gt2, 1])/Sqrt[2], 1}}, 
 {{conj[AdR[{gt1, ct1}]], Ah[{gt2}], Sd[{gt3, ct3}]}, 
  {(I*conj[ZD[gt3, gt1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yd[gt1, gt1]*
     ZA[gt2, 1])/Sqrt[2], 1}}, {{conj[ADxbarR[{gt1, ct1}]], Ah[{gt2}], 
   SDX[{gt3, ct3}]}, {(I*conj[ZDX[gt3, gt1]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*ZA[gt2, 3]*\[Kappa][gt1, gt1])/Sqrt[2], 1}}, 
 {{conj[ADxL[{gt1, ct1}]], Ah[{gt2}], SDX[{gt3, ct3}]}, 
  {((-I)*conj[ZDX[gt3, 3 + gt1]]*conj[\[Kappa][gt1, gt1]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*ZA[gt2, 3])/Sqrt[2], 1}}, 
 {{conj[AeL[{gt1}]], Ah[{gt2}], Se[{gt3}]}, 
  {((-I)*conj[Ye[gt1, gt1]]*conj[ZE[gt3, 3 + gt1]]*ThetaStep[gt1, 3]*
     ZA[gt2, 1])/Sqrt[2], 1}}, {{conj[AeR[{gt1}]], Ah[{gt2}], Se[{gt3}]}, 
  {(I*conj[ZE[gt3, gt1]]*ThetaStep[gt1, 3]*Ye[gt1, gt1]*ZA[gt2, 1])/Sqrt[2], 
   1}}, {{conj[AuL[{gt1, ct1}]], Ah[{gt2}], Su[{gt3, ct3}]}, 
  {((-I)*conj[Yu[gt1, gt1]]*conj[ZU[gt3, 3 + gt1]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*ZA[gt2, 2])/Sqrt[2], 1}}, 
 {{conj[AuR[{gt1, ct1}]], Ah[{gt2}], Su[{gt3, ct3}]}, 
  {(I*conj[ZU[gt3, gt1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yu[gt1, gt1]*
     ZA[gt2, 2])/Sqrt[2], 1}}, {{conj[AHup], Ah[{gt2}], conj[Hpm[{gt3}]]}, 
  {((-I)*conj[\[Lambda]]*ZA[gt2, 3]*ZP[gt3, 1])/Sqrt[2], 1}}, 
 {{AHd0, hh[{gt2}], hh[{gt3}]}, 
  {-(\[Lambda]*(conj[ZH[gt2, 3]]*conj[ZH[gt3, 2]] + conj[ZH[gt2, 2]]*
        conj[ZH[gt3, 3]]))/2, 1}}, 
 {{AHd0, Sd[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {Delta[ct2, ct3]*sum[j1, 1, 3, conj[ZD[gt2, j1]]*Yd[j1, j1]*
      ZD[gt3, 3 + j1]], 1}}, {{AHd0, Se[{gt2}], conj[Se[{gt3}]]}, 
  {sum[j1, 1, 3, conj[ZE[gt2, j1]]*Ye[j1, j1]*ZE[gt3, 3 + j1]], 1}}, 
 {{AHdm, hh[{gt2}], conj[Hpm[{gt3}]]}, 
  {(\[Lambda]*conj[ZH[gt2, 3]]*ZP[gt3, 2])/Sqrt[2], 1}}, 
 {{AHdm, Su[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {-(Delta[ct2, ct3]*sum[j1, 1, 3, conj[ZU[gt2, j1]]*Yd[j1, j1]*
       ZD[gt3, 3 + j1]]), 1}}, {{AHdm, Sv[{gt2}], conj[Se[{gt3}]]}, 
  {-sum[j1, 1, 3, conj[ZV[gt2, j1]]*Ye[j1, j1]*ZE[gt3, 3 + j1]], 1}}, 
 {{AHu0, hh[{gt2}], hh[{gt3}]}, 
  {-(\[Lambda]*(conj[ZH[gt2, 3]]*conj[ZH[gt3, 1]] + conj[ZH[gt2, 1]]*
        conj[ZH[gt3, 3]]))/2, 1}}, 
 {{AHu0, Su[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {Delta[ct2, ct3]*sum[j1, 1, 3, conj[ZU[gt2, j1]]*Yu[j1, j1]*
      ZU[gt3, 3 + j1]], 1}}, {{AHup, hh[{gt2}], Hpm[{gt3}]}, 
  {(\[Lambda]*conj[ZH[gt2, 3]]*conj[ZP[gt3, 1]])/Sqrt[2], 1}}, 
 {{AHup, Sd[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {-(Delta[ct2, ct3]*sum[j1, 1, 3, conj[ZD[gt2, j1]]*Yu[j1, j1]*
       ZU[gt3, 3 + j1]]), 1}}, {{AphiR, hh[{gt2}], hh[{gt3}]}, 
  {(-(\[Sigma]*conj[ZH[gt2, 4]]*conj[ZH[gt3, 3]]) - 
     \[Sigma]*conj[ZH[gt2, 3]]*conj[ZH[gt3, 4]] + 
     2*\[Kappa]Pr*conj[ZH[gt2, 5]]*conj[ZH[gt3, 5]])/2, 1}}, 
 {{AsbarR, hh[{gt2}], hh[{gt3}]}, 
  {-(\[Sigma]*(conj[ZH[gt2, 5]]*conj[ZH[gt3, 3]] + conj[ZH[gt2, 3]]*
        conj[ZH[gt3, 5]]))/2, 1}}, {{AsR, hh[{gt2}], hh[{gt3}]}, 
  {(-(\[Lambda]*conj[ZH[gt2, 2]]*conj[ZH[gt3, 1]]) - 
     \[Lambda]*conj[ZH[gt2, 1]]*conj[ZH[gt3, 2]] - \[Sigma]*conj[ZH[gt2, 5]]*
      conj[ZH[gt3, 4]] - \[Sigma]*conj[ZH[gt2, 4]]*conj[ZH[gt3, 5]])/2, 1}}, 
 {{AsR, Hpm[{gt2}], conj[Hpm[{gt3}]]}, 
  {\[Lambda]*conj[ZP[gt2, 1]]*ZP[gt3, 2], 1}}, 
 {{AsR, SDX[{gt2, ct2}], conj[SDX[{gt3, ct3}]]}, 
  {Delta[ct2, ct3]*sum[j1, 1, 3, conj[ZDX[gt2, j1]]*ZDX[gt3, 3 + j1]*
      \[Kappa][j1, j1]], 1}}, {{AuL[{gt1, ct1}], hh[{gt2}], 
   conj[Su[{gt3, ct3}]]}, 
  {(conj[ZH[gt2, 2]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yu[gt1, gt1]*
     ZU[gt3, 3 + gt1])/Sqrt[2], 1}}, 
 {{AuL[{gt1, ct1}], Hpm[{gt2}], conj[Sd[{gt3, ct3}]]}, 
  {-(conj[ZP[gt2, 1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*Yd[gt1, gt1]*
     ZD[gt3, 3 + gt1]), 1}}, {{AuR[{gt1, ct1}], hh[{gt2}], 
   conj[Su[{gt3, ct3}]]}, {(conj[Yu[gt1, gt1]]*conj[ZH[gt2, 2]]*
     Delta[ct1, ct3]*ThetaStep[gt1, 3]*ZU[gt3, gt1])/Sqrt[2], 1}}, 
 {{AuR[{gt1, ct1}], Hpm[{gt2}], conj[Sd[{gt3, ct3}]]}, 
  {-(conj[Yu[gt1, gt1]]*conj[ZP[gt2, 2]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*
     ZD[gt3, gt1]), 1}}, {{AvL[{gt1}], Hpm[{gt2}], conj[Se[{gt3}]]}, 
  {-(conj[ZP[gt2, 1]]*ThetaStep[gt1, 3]*Ye[gt1, gt1]*ZE[gt3, 3 + gt1]), 1}}, 
 {{aWB[{gt1}], hh[{gt2}], hh[{gt3}]}, 
  {(g2*(conj[ZH[gt2, 1]]*conj[ZH[gt3, 1]]*Sig[gt1, 1, 1] + 
      conj[ZH[gt2, 2]]*conj[ZH[gt3, 2]]*Sig[gt1, 2, 2]))/2, 1}}, 
 {{aWB[{gt1}], hh[{gt2}], Hpm[{gt3}]}, 
  {(g2*(conj[ZH[gt2, 1]]*conj[ZP[gt3, 1]] + conj[ZH[gt2, 2]]*
       conj[ZP[gt3, 2]])*Sig[gt1, 1, 2])/(2*Sqrt[2]), 1}}, 
 {{aWB[{gt1}], hh[{gt2}], conj[Hpm[{gt3}]]}, 
  {(g2*Sig[gt1, 2, 1]*(conj[ZH[gt2, 1]]*ZP[gt3, 1] + 
      conj[ZH[gt2, 2]]*ZP[gt3, 2]))/(2*Sqrt[2]), 1}}, 
 {{aWB[{gt1}], Hpm[{gt2}], conj[Hpm[{gt3}]]}, 
  {(g2*(conj[ZP[gt2, 1]]*Sig[gt1, 2, 2]*ZP[gt3, 1] + 
      conj[ZP[gt2, 2]]*Sig[gt1, 1, 1]*ZP[gt3, 2]))/2, 1}}, 
 {{aWB[{gt1}], Sd[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {(g2*Delta[ct2, ct3]*Sig[gt1, 2, 2]*sum[j1, 1, 3, 
      conj[ZD[gt2, j1]]*ZD[gt3, j1]])/2, 1}}, 
 {{aWB[{gt1}], Sd[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {(g2*Delta[ct2, ct3]*Sig[gt1, 1, 2]*sum[j1, 1, 3, 
      conj[ZD[gt2, j1]]*ZU[gt3, j1]])/2, 1}}, 
 {{aWB[{gt1}], Se[{gt2}], conj[Se[{gt3}]]}, 
  {(g2*Sig[gt1, 2, 2]*sum[j1, 1, 3, conj[ZE[gt2, j1]]*ZE[gt3, j1]])/2, 1}}, 
 {{aWB[{gt1}], Se[{gt2}], conj[Sv[{gt3}]]}, 
  {(g2*Sig[gt1, 1, 2]*sum[j1, 1, 3, conj[ZE[gt2, j1]]*ZV[gt3, j1]])/2, 1}}, 
 {{aWB[{gt1}], SHp0[{gt2}], conj[SHp0[{gt3}]]}, 
  {(g2*(conj[UHp0[gt2, 1]]*Sig[gt1, 1, 1]*UHp0[gt3, 1] + 
      conj[UHp0[gt2, 2]]*Sig[gt1, 2, 2]*UHp0[gt3, 2]))/2, 1}}, 
 {{aWB[{gt1}], SHp0[{gt2}], conj[SHpp[{gt3}]]}, 
  {(g2*Sig[gt1, 2, 1]*(conj[UHp0[gt2, 1]]*UHpp[gt3, 1] + 
      conj[UHp0[gt2, 2]]*UHpp[gt3, 2]))/2, 1}}, 
 {{aWB[{gt1}], SHpp[{gt2}], conj[SHp0[{gt3}]]}, 
  {(g2*Sig[gt1, 1, 2]*(conj[UHpp[gt2, 1]]*UHp0[gt3, 1] + 
      conj[UHpp[gt2, 2]]*UHp0[gt3, 2]))/2, 1}}, 
 {{aWB[{gt1}], SHpp[{gt2}], conj[SHpp[{gt3}]]}, 
  {(g2*(conj[UHpp[gt2, 1]]*Sig[gt1, 2, 2]*UHpp[gt3, 1] + 
      conj[UHpp[gt2, 2]]*Sig[gt1, 1, 1]*UHpp[gt3, 2]))/2, 1}}, 
 {{aWB[{gt1}], Su[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {(g2*Delta[ct2, ct3]*Sig[gt1, 2, 1]*sum[j1, 1, 3, 
      conj[ZU[gt2, j1]]*ZD[gt3, j1]])/2, 1}}, 
 {{aWB[{gt1}], Su[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {(g2*Delta[ct2, ct3]*Sig[gt1, 1, 1]*sum[j1, 1, 3, 
      conj[ZU[gt2, j1]]*ZU[gt3, j1]])/2, 1}}, 
 {{aWB[{gt1}], Sv[{gt2}], conj[Se[{gt3}]]}, 
  {(g2*Sig[gt1, 2, 1]*sum[j1, 1, 3, conj[ZV[gt2, j1]]*ZE[gt3, j1]])/2, 1}}, 
 {{aWB[{gt1}], Sv[{gt2}], conj[Sv[{gt3}]]}, 
  {(g2*Delta[gt2, gt3]*Sig[gt1, 1, 1])/2, 1}}, 
 {{conj[AHd0], hh[{gt2}], hh[{gt3}]}, 
  {-(conj[\[Lambda]]*(conj[ZH[gt2, 3]]*conj[ZH[gt3, 2]] + 
       conj[ZH[gt2, 2]]*conj[ZH[gt3, 3]]))/2, 1}}, 
 {{conj[AHu0], hh[{gt2}], hh[{gt3}]}, 
  {-(conj[\[Lambda]]*(conj[ZH[gt2, 3]]*conj[ZH[gt3, 1]] + 
       conj[ZH[gt2, 1]]*conj[ZH[gt3, 3]]))/2, 1}}, 
 {{conj[AphiR], hh[{gt2}], hh[{gt3}]}, 
  {(-(conj[\[Sigma]]*(conj[ZH[gt2, 4]]*conj[ZH[gt3, 3]] + 
        conj[ZH[gt2, 3]]*conj[ZH[gt3, 4]])) + 2*conj[\[Kappa]Pr]*
      conj[ZH[gt2, 5]]*conj[ZH[gt3, 5]])/2, 1}}, 
 {{conj[AsbarR], hh[{gt2}], hh[{gt3}]}, 
  {-(conj[\[Sigma]]*(conj[ZH[gt2, 5]]*conj[ZH[gt3, 3]] + 
       conj[ZH[gt2, 3]]*conj[ZH[gt3, 5]]))/2, 1}}, 
 {{conj[AsR], hh[{gt2}], hh[{gt3}]}, 
  {(-(conj[\[Lambda]]*(conj[ZH[gt2, 2]]*conj[ZH[gt3, 1]] + 
        conj[ZH[gt2, 1]]*conj[ZH[gt3, 2]])) - conj[\[Sigma]]*
      (conj[ZH[gt2, 5]]*conj[ZH[gt3, 4]] + conj[ZH[gt2, 4]]*
        conj[ZH[gt3, 5]]))/2, 1}}, {{conj[AHdm], hh[{gt2}], Hpm[{gt3}]}, 
  {(conj[\[Lambda]]*conj[ZH[gt2, 3]]*conj[ZP[gt3, 2]])/Sqrt[2], 1}}, 
 {{conj[AdL[{gt1, ct1}]], hh[{gt2}], Sd[{gt3, ct3}]}, 
  {(conj[Yd[gt1, gt1]]*conj[ZD[gt3, 3 + gt1]]*conj[ZH[gt2, 1]]*
     Delta[ct1, ct3]*ThetaStep[gt1, 3])/Sqrt[2], 1}}, 
 {{conj[AdR[{gt1, ct1}]], hh[{gt2}], Sd[{gt3, ct3}]}, 
  {(conj[ZD[gt3, gt1]]*conj[ZH[gt2, 1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*
     Yd[gt1, gt1])/Sqrt[2], 1}}, {{conj[ADxbarR[{gt1, ct1}]], hh[{gt2}], 
   SDX[{gt3, ct3}]}, {(conj[ZDX[gt3, gt1]]*conj[ZH[gt2, 3]]*Delta[ct1, ct3]*
     ThetaStep[gt1, 3]*\[Kappa][gt1, gt1])/Sqrt[2], 1}}, 
 {{conj[ADxL[{gt1, ct1}]], hh[{gt2}], SDX[{gt3, ct3}]}, 
  {(conj[ZDX[gt3, 3 + gt1]]*conj[ZH[gt2, 3]]*conj[\[Kappa][gt1, gt1]]*
     Delta[ct1, ct3]*ThetaStep[gt1, 3])/Sqrt[2], 1}}, 
 {{conj[AeL[{gt1}]], hh[{gt2}], Se[{gt3}]}, 
  {(conj[Ye[gt1, gt1]]*conj[ZE[gt3, 3 + gt1]]*conj[ZH[gt2, 1]]*
     ThetaStep[gt1, 3])/Sqrt[2], 1}}, 
 {{conj[AeR[{gt1}]], hh[{gt2}], Se[{gt3}]}, 
  {(conj[ZE[gt3, gt1]]*conj[ZH[gt2, 1]]*ThetaStep[gt1, 3]*Ye[gt1, gt1])/
    Sqrt[2], 1}}, {{conj[AuL[{gt1, ct1}]], hh[{gt2}], Su[{gt3, ct3}]}, 
  {(conj[Yu[gt1, gt1]]*conj[ZH[gt2, 2]]*conj[ZU[gt3, 3 + gt1]]*
     Delta[ct1, ct3]*ThetaStep[gt1, 3])/Sqrt[2], 1}}, 
 {{conj[AuR[{gt1, ct1}]], hh[{gt2}], Su[{gt3, ct3}]}, 
  {(conj[ZH[gt2, 2]]*conj[ZU[gt3, gt1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*
     Yu[gt1, gt1])/Sqrt[2], 1}}, {{conj[AHup], hh[{gt2}], conj[Hpm[{gt3}]]}, 
  {(conj[\[Lambda]]*conj[ZH[gt2, 3]]*ZP[gt3, 1])/Sqrt[2], 1}}, 
 {{conj[AdL[{gt1, ct1}]], Hpm[{gt2}], Su[{gt3, ct3}]}, 
  {-(conj[Yu[gt1, gt1]]*conj[ZP[gt2, 2]]*conj[ZU[gt3, 3 + gt1]]*
     Delta[ct1, ct3]*ThetaStep[gt1, 3]), 1}}, 
 {{conj[AdR[{gt1, ct1}]], Hpm[{gt2}], Su[{gt3, ct3}]}, 
  {-(conj[ZP[gt2, 1]]*conj[ZU[gt3, gt1]]*Delta[ct1, ct3]*ThetaStep[gt1, 3]*
     Yd[gt1, gt1]), 1}}, {{conj[AeR[{gt1}]], Hpm[{gt2}], Sv[{gt3}]}, 
  {-(conj[ZP[gt2, 1]]*conj[ZV[gt3, gt1]]*ThetaStep[gt1, 3]*Ye[gt1, gt1]), 
   1}}, {{conj[AsR], Hpm[{gt2}], conj[Hpm[{gt3}]]}, 
  {conj[\[Lambda]]*conj[ZP[gt2, 2]]*ZP[gt3, 1], 1}}, 
 {{conj[AHd0], Sd[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {Delta[ct2, ct3]*sum[j1, 1, 3, conj[Yd[j1, j1]]*conj[ZD[gt2, 3 + j1]]*
      ZD[gt3, j1]], 1}}, {{conj[AHdm], Sd[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {-(Delta[ct2, ct3]*sum[j1, 1, 3, conj[Yd[j1, j1]]*conj[ZD[gt2, 3 + j1]]*
       ZU[gt3, j1]]), 1}}, {{conj[AuL[{gt1, ct1}]], Sd[{gt2, ct2}], 
   conj[Hpm[{gt3}]]}, {-(conj[Yd[gt1, gt1]]*conj[ZD[gt2, 3 + gt1]]*
     Delta[ct1, ct2]*ThetaStep[gt1, 3]*ZP[gt3, 1]), 1}}, 
 {{conj[AuR[{gt1, ct1}]], Sd[{gt2, ct2}], conj[Hpm[{gt3}]]}, 
  {-(conj[ZD[gt2, gt1]]*Delta[ct1, ct2]*ThetaStep[gt1, 3]*Yu[gt1, gt1]*
     ZP[gt3, 2]), 1}}, {{conj[AsR], SDX[{gt2, ct2}], conj[SDX[{gt3, ct3}]]}, 
  {Delta[ct2, ct3]*sum[j1, 1, 3, conj[ZDX[gt2, 3 + j1]]*
      conj[\[Kappa][j1, j1]]*ZDX[gt3, j1]], 1}}, 
 {{conj[AHd0], Se[{gt2}], conj[Se[{gt3}]]}, 
  {sum[j1, 1, 3, conj[Ye[j1, j1]]*conj[ZE[gt2, 3 + j1]]*ZE[gt3, j1]], 1}}, 
 {{conj[AHdm], Se[{gt2}], conj[Sv[{gt3}]]}, 
  {-sum[j1, 1, 3, conj[Ye[j1, j1]]*conj[ZE[gt2, 3 + j1]]*ZV[gt3, j1]], 1}}, 
 {{conj[AvL[{gt1}]], Se[{gt2}], conj[Hpm[{gt3}]]}, 
  {-(conj[Ye[gt1, gt1]]*conj[ZE[gt2, 3 + gt1]]*ThetaStep[gt1, 3]*ZP[gt3, 1]), 
   1}}, {{conj[AHu0], Su[{gt2, ct2}], conj[Su[{gt3, ct3}]]}, 
  {Delta[ct2, ct3]*sum[j1, 1, 3, conj[Yu[j1, j1]]*conj[ZU[gt2, 3 + j1]]*
      ZU[gt3, j1]], 1}}, {{conj[AHup], Su[{gt2, ct2}], conj[Sd[{gt3, ct3}]]}, 
  {-(Delta[ct2, ct3]*sum[j1, 1, 3, conj[Yu[j1, j1]]*conj[ZU[gt2, 3 + j1]]*
       ZD[gt3, j1]]), 1}}}
