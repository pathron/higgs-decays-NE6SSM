{{mq2[i1, i2], (-2*g1^2*MassB*conj[MassB]*Kronecker[i1, i2])/15 - 
   (g1p^2*MassBp*conj[MassBp]*Kronecker[i1, i2])/5 - 
   (32*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/3 - 
   6*g2^2*MassWB*conj[MassWB]*Kronecker[i1, i2] + 
   (g1*Kronecker[i1, i2]*Tr1[1])/Sqrt[15] + (g1p*Kronecker[i1, i2]*Tr1[4])/
    Sqrt[10] + 2*mHd2*MatMul[Adj[Yd], Yd][i1, i2] + 
   2*mHu2*MatMul[Adj[Yu], Yu][i1, i2] + 2*MatMul[Adj[T[Yd]], T[Yd]][i1, i2] + 
   2*MatMul[Adj[T[Yu]], T[Yu]][i1, i2] + MatMul[mq2, Adj[Yd], Yd][i1, i2] + 
   MatMul[mq2, Adj[Yu], Yu][i1, i2] + 2*MatMul[Adj[Yd], md2, Yd][i1, i2] + 
   MatMul[Adj[Yd], Yd, mq2][i1, i2] + 2*MatMul[Adj[Yu], mu2, Yu][i1, i2] + 
   MatMul[Adj[Yu], Yu, mq2][i1, i2], 
  (16*g1^2*g3^2*MassB*conj[MassG]*Kronecker[i1, i2])/45 + 
   (8*g1p^2*g3^2*MassBp*conj[MassG]*Kronecker[i1, i2])/15 + 
   (32*g1^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/45 + 
   (16*g1p^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/15 + 
   32*g2^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2] + 
   (160*g3^4*MassG*conj[MassG]*Kronecker[i1, i2])/3 + 
   16*g2^2*g3^2*MassWB*conj[MassG]*Kronecker[i1, i2] + 
   (g1^2*g2^2*MassB*conj[MassWB]*Kronecker[i1, i2])/5 + 
   (3*g1p^2*g2^2*MassBp*conj[MassWB]*Kronecker[i1, i2])/10 + 
   16*g2^2*g3^2*MassG*conj[MassWB]*Kronecker[i1, i2] + 
   (2*g1^2*g2^2*MassWB*conj[MassWB]*Kronecker[i1, i2])/5 + 
   (3*g1p^2*g2^2*MassWB*conj[MassWB]*Kronecker[i1, i2])/5 + 
   51*g2^4*MassWB*conj[MassWB]*Kronecker[i1, i2] + 
   32*g2^2*g3^2*MassWB*conj[MassWB]*Kronecker[i1, i2] + 
   6*g2^4*Kronecker[i1, i2]*Tr2[2] + (32*g3^4*Kronecker[i1, i2]*Tr2[3])/3 + 
   (2*g1^2*Kronecker[i1, i2]*Tr2U1[1, 1])/15 + 
   (Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[1, 4])/5 + 
   (Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[4, 1])/5 + 
   (g1p^2*Kronecker[i1, i2]*Tr2U1[4, 4])/5 + (4*g1*Kronecker[i1, i2]*Tr3[1])/
    Sqrt[15] + 2*Sqrt[2/5]*g1p*Kronecker[i1, i2]*Tr3[4] + 
   (4*g1^2*mHd2*MatMul[Adj[Yd], Yd][i1, i2])/5 + 
   (6*g1p^2*mHd2*MatMul[Adj[Yd], Yd][i1, i2])/5 - 
   4*mHd2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   2*mHu2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   2*ms2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   2*conj[T[\[Lambda]]]*T[\[Lambda]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   12*mHd2*trace[Yd, Adj[Yd]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   4*mHd2*trace[Ye, Adj[Ye]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   6*trace[conj[T[Yd]], Tp[T[Yd]]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   2*trace[conj[T[Ye]], Tp[T[Ye]]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   6*trace[md2, Yd, Adj[Yd]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   2*trace[me2, Ye, Adj[Ye]]*MatMul[Adj[Yd], Yd][i1, i2] - 
   2*trace[ml2, Adj[Ye], Ye]*MatMul[Adj[Yd], Yd][i1, i2] - 
   6*trace[mq2, Adj[Yd], Yd]*MatMul[Adj[Yd], Yd][i1, i2] - 
   2*\[Lambda]*conj[T[\[Lambda]]]*MatMul[Adj[Yd], T[Yd]][i1, i2] - 
   6*trace[conj[T[Yd]], Tp[Yd]]*MatMul[Adj[Yd], T[Yd]][i1, i2] - 
   2*trace[conj[T[Ye]], Tp[Ye]]*MatMul[Adj[Yd], T[Yd]][i1, i2] + 
   (8*g1^2*mHu2*MatMul[Adj[Yu], Yu][i1, i2])/5 + 
   (2*g1p^2*mHu2*MatMul[Adj[Yu], Yu][i1, i2])/5 - 
   2*mHd2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yu], Yu][i1, i2] - 
   4*mHu2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yu], Yu][i1, i2] - 
   2*ms2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yu], Yu][i1, i2] - 
   2*conj[T[\[Lambda]]]*T[\[Lambda]]*MatMul[Adj[Yu], Yu][i1, i2] - 
   12*mHu2*trace[Yu, Adj[Yu]]*MatMul[Adj[Yu], Yu][i1, i2] - 
   6*trace[conj[T[Yu]], Tp[T[Yu]]]*MatMul[Adj[Yu], Yu][i1, i2] - 
   6*trace[mq2, Adj[Yu], Yu]*MatMul[Adj[Yu], Yu][i1, i2] - 
   6*trace[mu2, Yu, Adj[Yu]]*MatMul[Adj[Yu], Yu][i1, i2] + 
   (g1^2*conj[MassB]*((1518*g1^2*MassB - 51*g1p^2*(2*MassB + MassBp) + 
        10*(16*g3^2*(2*MassB + MassG) + 9*g2^2*(2*MassB + MassWB)))*
       Kronecker[i1, i2] + 360*(2*MassB*MatMul[Adj[Yd], Yd][i1, i2] - 
        MatMul[Adj[Yd], T[Yd]][i1, i2] + 4*MassB*MatMul[Adj[Yu], Yu][i1, 
          i2] - 2*MatMul[Adj[Yu], T[Yu]][i1, i2])))/450 + 
   (g1p^2*conj[MassBp]*((-68*g1^2*(MassB + 2*MassBp) + 
        20*(16*g3^2*(2*MassBp + MassG) + 9*g2^2*(2*MassBp + MassWB)) + 
        9*g1p^2*MassBp*(251 + 2*QS^2))*Kronecker[i1, i2] + 
      240*(6*MassBp*MatMul[Adj[Yd], Yd][i1, i2] - 
        3*MatMul[Adj[Yd], T[Yd]][i1, i2] + 2*MassBp*MatMul[Adj[Yu], Yu][i1, 
          i2] - MatMul[Adj[Yu], T[Yu]][i1, i2])))/600 - 
   2*\[Lambda]*conj[T[\[Lambda]]]*MatMul[Adj[Yu], T[Yu]][i1, i2] - 
   6*trace[conj[T[Yu]], Tp[Yu]]*MatMul[Adj[Yu], T[Yu]][i1, i2] - 
   (4*g1^2*MassB*MatMul[Adj[T[Yd]], Yd][i1, i2])/5 - 
   (6*g1p^2*MassBp*MatMul[Adj[T[Yd]], Yd][i1, i2])/5 - 
   2*conj[\[Lambda]]*T[\[Lambda]]*MatMul[Adj[T[Yd]], Yd][i1, i2] - 
   6*trace[Adj[Yd], T[Yd]]*MatMul[Adj[T[Yd]], Yd][i1, i2] - 
   2*trace[Adj[Ye], T[Ye]]*MatMul[Adj[T[Yd]], Yd][i1, i2] + 
   (4*g1^2*MatMul[Adj[T[Yd]], T[Yd]][i1, i2])/5 + 
   (6*g1p^2*MatMul[Adj[T[Yd]], T[Yd]][i1, i2])/5 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[T[Yd]], T[Yd]][i1, i2] - 
   6*trace[Yd, Adj[Yd]]*MatMul[Adj[T[Yd]], T[Yd]][i1, i2] - 
   2*trace[Ye, Adj[Ye]]*MatMul[Adj[T[Yd]], T[Yd]][i1, i2] - 
   (8*g1^2*MassB*MatMul[Adj[T[Yu]], Yu][i1, i2])/5 - 
   (2*g1p^2*MassBp*MatMul[Adj[T[Yu]], Yu][i1, i2])/5 - 
   2*conj[\[Lambda]]*T[\[Lambda]]*MatMul[Adj[T[Yu]], Yu][i1, i2] - 
   6*trace[Adj[Yu], T[Yu]]*MatMul[Adj[T[Yu]], Yu][i1, i2] + 
   (8*g1^2*MatMul[Adj[T[Yu]], T[Yu]][i1, i2])/5 + 
   (2*g1p^2*MatMul[Adj[T[Yu]], T[Yu]][i1, i2])/5 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[T[Yu]], T[Yu]][i1, i2] - 
   6*trace[Yu, Adj[Yu]]*MatMul[Adj[T[Yu]], T[Yu]][i1, i2] + 
   (2*g1^2*MatMul[mq2, Adj[Yd], Yd][i1, i2])/5 + 
   (3*g1p^2*MatMul[mq2, Adj[Yd], Yd][i1, i2])/5 - 
   \[Lambda]*conj[\[Lambda]]*MatMul[mq2, Adj[Yd], Yd][i1, i2] - 
   3*trace[Yd, Adj[Yd]]*MatMul[mq2, Adj[Yd], Yd][i1, i2] - 
   trace[Ye, Adj[Ye]]*MatMul[mq2, Adj[Yd], Yd][i1, i2] + 
   (4*g1^2*MatMul[mq2, Adj[Yu], Yu][i1, i2])/5 + 
   (g1p^2*MatMul[mq2, Adj[Yu], Yu][i1, i2])/5 - \[Lambda]*conj[\[Lambda]]*
    MatMul[mq2, Adj[Yu], Yu][i1, i2] - 3*trace[Yu, Adj[Yu]]*
    MatMul[mq2, Adj[Yu], Yu][i1, i2] + 
   (4*g1^2*MatMul[Adj[Yd], md2, Yd][i1, i2])/5 + 
   (6*g1p^2*MatMul[Adj[Yd], md2, Yd][i1, i2])/5 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yd], md2, Yd][i1, i2] - 
   6*trace[Yd, Adj[Yd]]*MatMul[Adj[Yd], md2, Yd][i1, i2] - 
   2*trace[Ye, Adj[Ye]]*MatMul[Adj[Yd], md2, Yd][i1, i2] + 
   (2*g1^2*MatMul[Adj[Yd], Yd, mq2][i1, i2])/5 + 
   (3*g1p^2*MatMul[Adj[Yd], Yd, mq2][i1, i2])/5 - 
   \[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yd], Yd, mq2][i1, i2] - 
   3*trace[Yd, Adj[Yd]]*MatMul[Adj[Yd], Yd, mq2][i1, i2] - 
   trace[Ye, Adj[Ye]]*MatMul[Adj[Yd], Yd, mq2][i1, i2] + 
   (8*g1^2*MatMul[Adj[Yu], mu2, Yu][i1, i2])/5 + 
   (2*g1p^2*MatMul[Adj[Yu], mu2, Yu][i1, i2])/5 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Yu], mu2, Yu][i1, i2] - 
   6*trace[Yu, Adj[Yu]]*MatMul[Adj[Yu], mu2, Yu][i1, i2] + 
   (4*g1^2*MatMul[Adj[Yu], Yu, mq2][i1, i2])/5 + 
   (g1p^2*MatMul[Adj[Yu], Yu, mq2][i1, i2])/5 - \[Lambda]*conj[\[Lambda]]*
    MatMul[Adj[Yu], Yu, mq2][i1, i2] - 3*trace[Yu, Adj[Yu]]*
    MatMul[Adj[Yu], Yu, mq2][i1, i2] - 
   8*mHd2*MatMul[Adj[Yd], Yd, Adj[Yd], Yd][i1, i2] - 
   4*MatMul[Adj[Yd], Yd, Adj[T[Yd]], T[Yd]][i1, i2] - 
   4*MatMul[Adj[Yd], T[Yd], Adj[T[Yd]], Yd][i1, i2] - 
   8*mHu2*MatMul[Adj[Yu], Yu, Adj[Yu], Yu][i1, i2] - 
   4*MatMul[Adj[Yu], Yu, Adj[T[Yu]], T[Yu]][i1, i2] - 
   4*MatMul[Adj[Yu], T[Yu], Adj[T[Yu]], Yu][i1, i2] - 
   4*MatMul[Adj[T[Yd]], Yd, Adj[Yd], T[Yd]][i1, i2] - 
   4*MatMul[Adj[T[Yd]], T[Yd], Adj[Yd], Yd][i1, i2] - 
   4*MatMul[Adj[T[Yu]], Yu, Adj[Yu], T[Yu]][i1, i2] - 
   4*MatMul[Adj[T[Yu]], T[Yu], Adj[Yu], Yu][i1, i2] - 
   2*MatMul[mq2, Adj[Yd], Yd, Adj[Yd], Yd][i1, i2] - 
   2*MatMul[mq2, Adj[Yu], Yu, Adj[Yu], Yu][i1, i2] - 
   4*MatMul[Adj[Yd], md2, Yd, Adj[Yd], Yd][i1, i2] - 
   4*MatMul[Adj[Yd], Yd, mq2, Adj[Yd], Yd][i1, i2] - 
   4*MatMul[Adj[Yd], Yd, Adj[Yd], md2, Yd][i1, i2] - 
   2*MatMul[Adj[Yd], Yd, Adj[Yd], Yd, mq2][i1, i2] - 
   4*MatMul[Adj[Yu], mu2, Yu, Adj[Yu], Yu][i1, i2] - 
   4*MatMul[Adj[Yu], Yu, mq2, Adj[Yu], Yu][i1, i2] - 
   4*MatMul[Adj[Yu], Yu, Adj[Yu], mu2, Yu][i1, i2] - 
   2*MatMul[Adj[Yu], Yu, Adj[Yu], Yu, mq2][i1, i2]}, 
 {ml2[i1, i2], (-6*g1^2*MassB*conj[MassB]*Kronecker[i1, i2])/5 - 
   (4*g1p^2*MassBp*conj[MassBp]*Kronecker[i1, i2])/5 - 
   6*g2^2*MassWB*conj[MassWB]*Kronecker[i1, i2] - 
   Sqrt[3/5]*g1*Kronecker[i1, i2]*Tr1[1] + Sqrt[2/5]*g1p*Kronecker[i1, i2]*
    Tr1[4] + 2*mHd2*MatMul[Adj[Ye], Ye][i1, i2] + 
   2*MatMul[Adj[T[Ye]], T[Ye]][i1, i2] + MatMul[ml2, Adj[Ye], Ye][i1, i2] + 
   2*MatMul[Adj[Ye], me2, Ye][i1, i2] + MatMul[Adj[Ye], Ye, ml2][i1, i2], 
  (9*g1^2*g2^2*MassB*conj[MassWB]*Kronecker[i1, i2])/5 + 
   (6*g1p^2*g2^2*MassBp*conj[MassWB]*Kronecker[i1, i2])/5 + 
   (18*g1^2*g2^2*MassWB*conj[MassWB]*Kronecker[i1, i2])/5 + 
   (12*g1p^2*g2^2*MassWB*conj[MassWB]*Kronecker[i1, i2])/5 + 
   51*g2^4*MassWB*conj[MassWB]*Kronecker[i1, i2] + 
   6*g2^4*Kronecker[i1, i2]*Tr2[2] + (6*g1^2*Kronecker[i1, i2]*Tr2U1[1, 1])/
    5 - (2*Sqrt[6]*g1*g1p*Kronecker[i1, i2]*Tr2U1[1, 4])/5 - 
   (2*Sqrt[6]*g1*g1p*Kronecker[i1, i2]*Tr2U1[4, 1])/5 + 
   (4*g1p^2*Kronecker[i1, i2]*Tr2U1[4, 4])/5 - 
   4*Sqrt[3/5]*g1*Kronecker[i1, i2]*Tr3[1] + 
   4*Sqrt[2/5]*g1p*Kronecker[i1, i2]*Tr3[4] + 
   (12*g1^2*mHd2*MatMul[Adj[Ye], Ye][i1, i2])/5 + 
   (3*g1p^2*mHd2*MatMul[Adj[Ye], Ye][i1, i2])/5 - 
   4*mHd2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   2*mHu2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   2*ms2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   2*conj[T[\[Lambda]]]*T[\[Lambda]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   12*mHd2*trace[Yd, Adj[Yd]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   4*mHd2*trace[Ye, Adj[Ye]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   6*trace[conj[T[Yd]], Tp[T[Yd]]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   2*trace[conj[T[Ye]], Tp[T[Ye]]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   6*trace[md2, Yd, Adj[Yd]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   2*trace[me2, Ye, Adj[Ye]]*MatMul[Adj[Ye], Ye][i1, i2] - 
   2*trace[ml2, Adj[Ye], Ye]*MatMul[Adj[Ye], Ye][i1, i2] - 
   6*trace[mq2, Adj[Yd], Yd]*MatMul[Adj[Ye], Ye][i1, i2] + 
   (3*g1^2*conj[MassB]*((261*g1^2*MassB + 8*g1p^2*(2*MassB + MassBp) + 
        15*g2^2*(2*MassB + MassWB))*Kronecker[i1, i2] + 
      40*MassB*MatMul[Adj[Ye], Ye][i1, i2] - 
      20*MatMul[Adj[Ye], T[Ye]][i1, i2]))/25 + 
   (3*g1p^2*conj[MassBp]*((16*g1^2*(MassB + 2*MassBp) + 
        20*g2^2*(2*MassBp + MassWB) + g1p^2*MassBp*(257 + 2*QS^2))*
       Kronecker[i1, i2] + 20*MassBp*MatMul[Adj[Ye], Ye][i1, i2] - 
      10*MatMul[Adj[Ye], T[Ye]][i1, i2]))/50 - 2*\[Lambda]*conj[T[\[Lambda]]]*
    MatMul[Adj[Ye], T[Ye]][i1, i2] - 6*trace[conj[T[Yd]], Tp[Yd]]*
    MatMul[Adj[Ye], T[Ye]][i1, i2] - 2*trace[conj[T[Ye]], Tp[Ye]]*
    MatMul[Adj[Ye], T[Ye]][i1, i2] - 
   (12*g1^2*MassB*MatMul[Adj[T[Ye]], Ye][i1, i2])/5 - 
   (3*g1p^2*MassBp*MatMul[Adj[T[Ye]], Ye][i1, i2])/5 - 
   2*conj[\[Lambda]]*T[\[Lambda]]*MatMul[Adj[T[Ye]], Ye][i1, i2] - 
   6*trace[Adj[Yd], T[Yd]]*MatMul[Adj[T[Ye]], Ye][i1, i2] - 
   2*trace[Adj[Ye], T[Ye]]*MatMul[Adj[T[Ye]], Ye][i1, i2] + 
   (12*g1^2*MatMul[Adj[T[Ye]], T[Ye]][i1, i2])/5 + 
   (3*g1p^2*MatMul[Adj[T[Ye]], T[Ye]][i1, i2])/5 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[T[Ye]], T[Ye]][i1, i2] - 
   6*trace[Yd, Adj[Yd]]*MatMul[Adj[T[Ye]], T[Ye]][i1, i2] - 
   2*trace[Ye, Adj[Ye]]*MatMul[Adj[T[Ye]], T[Ye]][i1, i2] + 
   (6*g1^2*MatMul[ml2, Adj[Ye], Ye][i1, i2])/5 + 
   (3*g1p^2*MatMul[ml2, Adj[Ye], Ye][i1, i2])/10 - 
   \[Lambda]*conj[\[Lambda]]*MatMul[ml2, Adj[Ye], Ye][i1, i2] - 
   3*trace[Yd, Adj[Yd]]*MatMul[ml2, Adj[Ye], Ye][i1, i2] - 
   trace[Ye, Adj[Ye]]*MatMul[ml2, Adj[Ye], Ye][i1, i2] + 
   (12*g1^2*MatMul[Adj[Ye], me2, Ye][i1, i2])/5 + 
   (3*g1p^2*MatMul[Adj[Ye], me2, Ye][i1, i2])/5 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[Adj[Ye], me2, Ye][i1, i2] - 
   6*trace[Yd, Adj[Yd]]*MatMul[Adj[Ye], me2, Ye][i1, i2] - 
   2*trace[Ye, Adj[Ye]]*MatMul[Adj[Ye], me2, Ye][i1, i2] + 
   (6*g1^2*MatMul[Adj[Ye], Ye, ml2][i1, i2])/5 + 
   (3*g1p^2*MatMul[Adj[Ye], Ye, ml2][i1, i2])/10 - 
   \[Lambda]*conj[\[Lambda]]*MatMul[Adj[Ye], Ye, ml2][i1, i2] - 
   3*trace[Yd, Adj[Yd]]*MatMul[Adj[Ye], Ye, ml2][i1, i2] - 
   trace[Ye, Adj[Ye]]*MatMul[Adj[Ye], Ye, ml2][i1, i2] - 
   8*mHd2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye][i1, i2] - 
   4*MatMul[Adj[Ye], Ye, Adj[T[Ye]], T[Ye]][i1, i2] - 
   4*MatMul[Adj[Ye], T[Ye], Adj[T[Ye]], Ye][i1, i2] - 
   4*MatMul[Adj[T[Ye]], Ye, Adj[Ye], T[Ye]][i1, i2] - 
   4*MatMul[Adj[T[Ye]], T[Ye], Adj[Ye], Ye][i1, i2] - 
   2*MatMul[ml2, Adj[Ye], Ye, Adj[Ye], Ye][i1, i2] - 
   4*MatMul[Adj[Ye], me2, Ye, Adj[Ye], Ye][i1, i2] - 
   4*MatMul[Adj[Ye], Ye, ml2, Adj[Ye], Ye][i1, i2] - 
   4*MatMul[Adj[Ye], Ye, Adj[Ye], me2, Ye][i1, i2] - 
   2*MatMul[Adj[Ye], Ye, Adj[Ye], Ye, ml2][i1, i2]}, 
 {mHd2, (-6*g1^2*MassB*conj[MassB])/5 - (9*g1p^2*MassBp*conj[MassBp])/5 - 
   6*g2^2*MassWB*conj[MassWB] + 2*mHd2*\[Lambda]*conj[\[Lambda]] + 
   2*mHu2*\[Lambda]*conj[\[Lambda]] + 2*ms2*\[Lambda]*conj[\[Lambda]] + 
   2*conj[T[\[Lambda]]]*T[\[Lambda]] - Sqrt[3/5]*g1*Tr1[1] - 
   (3*g1p*Tr1[4])/Sqrt[10] + 6*mHd2*trace[Yd, Adj[Yd]] + 
   2*mHd2*trace[Ye, Adj[Ye]] + 6*trace[conj[T[Yd]], Tp[T[Yd]]] + 
   2*trace[conj[T[Ye]], Tp[T[Ye]]] + 6*trace[md2, Yd, Adj[Yd]] + 
   2*trace[me2, Ye, Adj[Ye]] + 2*trace[ml2, Adj[Ye], Ye] + 
   6*trace[mq2, Adj[Yd], Yd], 
  (4*g1^2*conj[MassB]*(1566*g1^2*MassB - 54*g1p^2*MassB + 180*g2^2*MassB - 
      27*g1p^2*MassBp + 90*g2^2*MassWB - 80*MassB*trace[Yd, Adj[Yd]] + 
      240*MassB*trace[Ye, Adj[Ye]] + 40*trace[Adj[Yd], T[Yd]] - 
      120*trace[Adj[Ye], T[Ye]]) + g1p^2*conj[MassBp]*
     (-108*g1^2*MassB - 216*g1^2*MassBp + 7209*g1p^2*MassBp + 
      1080*g2^2*MassBp + 540*g2^2*MassWB + 54*g1p^2*MassBp*QS^2 + 
      20*(-5 + QS^2)*conj[\[Lambda]]*(2*MassBp*\[Lambda] - T[\[Lambda]]) - 
      480*MassBp*trace[Yd, Adj[Yd]] - 160*MassBp*trace[Ye, Adj[Ye]] + 
      240*trace[Adj[Yd], T[Yd]] + 80*trace[Adj[Ye], T[Ye]]) + 
    20*(3*g2^2*(170*g2^2*MassWB + 6*g1^2*(MassB + 2*MassWB) + 
        9*g1p^2*(MassBp + 2*MassWB))*conj[MassWB] - 120*(mHd2 + mHu2 + ms2)*
       \[Lambda]^2*conj[\[Lambda]]^2 + 5*g1p^2*MassBp*\[Lambda]*
       conj[T[\[Lambda]]] - g1p^2*MassBp*QS^2*\[Lambda]*conj[T[\[Lambda]]] - 
      5*g1p^2*conj[T[\[Lambda]]]*T[\[Lambda]] + g1p^2*QS^2*conj[T[\[Lambda]]]*
       T[\[Lambda]] - 20*\[Sigma]*conj[\[Sigma]]*conj[T[\[Lambda]]]*
       T[\[Lambda]] - 20*\[Lambda]*conj[\[Sigma]]*conj[T[\[Lambda]]]*
       T[\[Sigma]] + 60*g2^4*Tr2[2] + 12*g1^2*Tr2U1[1, 1] + 
      6*Sqrt[6]*g1*g1p*Tr2U1[1, 4] + 6*Sqrt[6]*g1*g1p*Tr2U1[4, 1] + 
      18*g1p^2*Tr2U1[4, 4] - 8*Sqrt[15]*g1*Tr3[1] - 12*Sqrt[10]*g1p*Tr3[4] - 
      8*g1^2*mHd2*trace[Yd, Adj[Yd]] - 12*g1p^2*mHd2*trace[Yd, Adj[Yd]] + 
      320*g3^2*mHd2*trace[Yd, Adj[Yd]] + 640*g3^2*MassG*conj[MassG]*
       trace[Yd, Adj[Yd]] + 24*g1^2*mHd2*trace[Ye, Adj[Ye]] - 
      4*g1p^2*mHd2*trace[Ye, Adj[Ye]] - 60*conj[T[\[Lambda]]]*T[\[Lambda]]*
       trace[Yu, Adj[Yu]] - 60*conj[T[\[Lambda]]]*T[\[Lambda]]*
       trace[\[Kappa], Adj[\[Kappa]]] - 320*g3^2*conj[MassG]*
       trace[Adj[Yd], T[Yd]] - 60*\[Lambda]*conj[T[\[Lambda]]]*
       trace[Adj[Yu], T[Yu]] - 60*\[Lambda]*conj[T[\[Lambda]]]*
       trace[Adj[\[Kappa]], T[\[Kappa]]] + 8*g1^2*MassB*
       trace[conj[T[Yd]], Tp[Yd]] + 12*g1p^2*MassBp*trace[conj[T[Yd]], 
        Tp[Yd]] - 320*g3^2*MassG*trace[conj[T[Yd]], Tp[Yd]] - 
      8*g1^2*trace[conj[T[Yd]], Tp[T[Yd]]] - 
      12*g1p^2*trace[conj[T[Yd]], Tp[T[Yd]]] + 
      320*g3^2*trace[conj[T[Yd]], Tp[T[Yd]]] - 24*g1^2*MassB*
       trace[conj[T[Ye]], Tp[Ye]] + 4*g1p^2*MassBp*trace[conj[T[Ye]], 
        Tp[Ye]] + 24*g1^2*trace[conj[T[Ye]], Tp[T[Ye]]] - 
      4*g1p^2*trace[conj[T[Ye]], Tp[T[Ye]]] - 
      8*g1^2*trace[md2, Yd, Adj[Yd]] - 12*g1p^2*trace[md2, Yd, Adj[Yd]] + 
      320*g3^2*trace[md2, Yd, Adj[Yd]] + 24*g1^2*trace[me2, Ye, Adj[Ye]] - 
      4*g1p^2*trace[me2, Ye, Adj[Ye]] + 24*g1^2*trace[ml2, Adj[Ye], Ye] - 
      4*g1p^2*trace[ml2, Adj[Ye], Ye] - 8*g1^2*trace[mq2, Adj[Yd], Yd] - 
      12*g1p^2*trace[mq2, Adj[Yd], Yd] + 320*g3^2*trace[mq2, Adj[Yd], Yd] + 
      conj[\[Lambda]]*(-5*g1p^2*mHd2*\[Lambda] - 5*g1p^2*mHu2*\[Lambda] - 
        5*g1p^2*ms2*\[Lambda] + g1p^2*mHd2*QS^2*\[Lambda] + 
        g1p^2*mHu2*QS^2*\[Lambda] + g1p^2*ms2*QS^2*\[Lambda] - 
        20*(mHd2 + mHu2 + mphi2 + 2*ms2 + msbar2)*\[Lambda]*\[Sigma]*
         conj[\[Sigma]] - 240*\[Lambda]*conj[T[\[Lambda]]]*T[\[Lambda]] - 
        20*\[Sigma]*conj[T[\[Sigma]]]*T[\[Lambda]] - 
        20*\[Lambda]*conj[T[\[Sigma]]]*T[\[Sigma]] - 60*mHd2*\[Lambda]*
         trace[Yu, Adj[Yu]] - 120*mHu2*\[Lambda]*trace[Yu, Adj[Yu]] - 
        60*ms2*\[Lambda]*trace[Yu, Adj[Yu]] - 60*mHd2*\[Lambda]*
         trace[\[Kappa], Adj[\[Kappa]]] - 60*mHu2*\[Lambda]*
         trace[\[Kappa], Adj[\[Kappa]]] - 120*ms2*\[Lambda]*
         trace[\[Kappa], Adj[\[Kappa]]] - 60*T[\[Lambda]]*
         trace[conj[T[Yu]], Tp[Yu]] - 60*\[Lambda]*trace[conj[T[Yu]], 
          Tp[T[Yu]]] - 60*T[\[Lambda]]*trace[conj[T[\[Kappa]]], 
          Tp[\[Kappa]]] - 60*\[Lambda]*trace[conj[T[\[Kappa]]], 
          Tp[T[\[Kappa]]]] - 60*\[Lambda]*trace[mq2, Adj[Yu], Yu] - 
        60*\[Lambda]*trace[mu2, Yu, Adj[Yu]] - 60*\[Lambda]*
         trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] - 
        60*\[Lambda]*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]]) - 
      360*mHd2*trace[Yd, Adj[Yd], Yd, Adj[Yd]] - 
      360*trace[Yd, Adj[Yd], T[Yd], Adj[T[Yd]]] - 
      60*mHd2*trace[Yd, Adj[Yu], Yu, Adj[Yd]] - 
      60*mHu2*trace[Yd, Adj[Yu], Yu, Adj[Yd]] - 
      60*trace[Yd, Adj[Yu], T[Yu], Adj[T[Yd]]] - 
      360*trace[Yd, Adj[T[Yd]], T[Yd], Adj[Yd]] - 
      60*trace[Yd, Adj[T[Yu]], T[Yu], Adj[Yd]] - 
      120*mHd2*trace[Ye, Adj[Ye], Ye, Adj[Ye]] - 
      120*trace[Ye, Adj[Ye], T[Ye], Adj[T[Ye]]] - 
      120*trace[Ye, Adj[T[Ye]], T[Ye], Adj[Ye]] - 
      60*trace[Yu, Adj[Yd], T[Yd], Adj[T[Yu]]] - 
      60*trace[Yu, Adj[T[Yd]], T[Yd], Adj[Yu]] - 
      360*trace[md2, Yd, Adj[Yd], Yd, Adj[Yd]] - 
      60*trace[md2, Yd, Adj[Yu], Yu, Adj[Yd]] - 
      120*trace[me2, Ye, Adj[Ye], Ye, Adj[Ye]] - 
      120*trace[ml2, Adj[Ye], Ye, Adj[Ye], Ye] - 
      360*trace[mq2, Adj[Yd], Yd, Adj[Yd], Yd] - 
      60*trace[mq2, Adj[Yd], Yd, Adj[Yu], Yu] - 
      60*trace[mq2, Adj[Yu], Yu, Adj[Yd], Yd] - 
      60*trace[mu2, Yu, Adj[Yd], Yd, Adj[Yu]]))/200}, 
 {mHu2, (-6*g1^2*MassB*conj[MassB])/5 - (4*g1p^2*MassBp*conj[MassBp])/5 - 
   6*g2^2*MassWB*conj[MassWB] + 2*mHd2*\[Lambda]*conj[\[Lambda]] + 
   2*mHu2*\[Lambda]*conj[\[Lambda]] + 2*ms2*\[Lambda]*conj[\[Lambda]] + 
   2*conj[T[\[Lambda]]]*T[\[Lambda]] + Sqrt[3/5]*g1*Tr1[1] - 
   Sqrt[2/5]*g1p*Tr1[4] + 6*mHu2*trace[Yu, Adj[Yu]] + 
   6*trace[conj[T[Yu]], Tp[T[Yu]]] + 6*trace[mq2, Adj[Yu], Yu] + 
   6*trace[mu2, Yu, Adj[Yu]], (9*g1^2*g2^2*MassB*conj[MassWB])/5 + 
   (6*g1p^2*g2^2*MassBp*conj[MassWB])/5 + (18*g1^2*g2^2*MassWB*conj[MassWB])/
    5 + (12*g1p^2*g2^2*MassWB*conj[MassWB])/5 + 51*g2^4*MassWB*conj[MassWB] + 
   (g1p^2*mHd2*\[Lambda]*conj[\[Lambda]])/2 + 
   (g1p^2*mHu2*\[Lambda]*conj[\[Lambda]])/2 + 
   (g1p^2*ms2*\[Lambda]*conj[\[Lambda]])/2 + 
   (g1p^2*mHd2*QS^2*\[Lambda]*conj[\[Lambda]])/10 + 
   (g1p^2*mHu2*QS^2*\[Lambda]*conj[\[Lambda]])/10 + 
   (g1p^2*ms2*QS^2*\[Lambda]*conj[\[Lambda]])/10 - 
   12*mHd2*\[Lambda]^2*conj[\[Lambda]]^2 - 12*mHu2*\[Lambda]^2*
    conj[\[Lambda]]^2 - 12*ms2*\[Lambda]^2*conj[\[Lambda]]^2 - 
   2*mHd2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   2*mHu2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   2*mphi2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   4*ms2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   2*msbar2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   (g1p^2*MassBp*\[Lambda]*conj[T[\[Lambda]]])/2 - 
   (g1p^2*MassBp*QS^2*\[Lambda]*conj[T[\[Lambda]]])/10 + 
   (g1p^2*conj[T[\[Lambda]]]*T[\[Lambda]])/2 + 
   (g1p^2*QS^2*conj[T[\[Lambda]]]*T[\[Lambda]])/10 - 
   24*\[Lambda]*conj[\[Lambda]]*conj[T[\[Lambda]]]*T[\[Lambda]] - 
   2*\[Sigma]*conj[\[Sigma]]*conj[T[\[Lambda]]]*T[\[Lambda]] - 
   2*\[Sigma]*conj[\[Lambda]]*conj[T[\[Sigma]]]*T[\[Lambda]] - 
   2*\[Lambda]*conj[\[Sigma]]*conj[T[\[Lambda]]]*T[\[Sigma]] - 
   2*\[Lambda]*conj[\[Lambda]]*conj[T[\[Sigma]]]*T[\[Sigma]] + 
   6*g2^4*Tr2[2] + (6*g1^2*Tr2U1[1, 1])/5 - (2*Sqrt[6]*g1*g1p*Tr2U1[1, 4])/
    5 - (2*Sqrt[6]*g1*g1p*Tr2U1[4, 1])/5 + (4*g1p^2*Tr2U1[4, 4])/5 + 
   4*Sqrt[3/5]*g1*Tr3[1] - 4*Sqrt[2/5]*g1p*Tr3[4] - 
   12*mHd2*\[Lambda]*conj[\[Lambda]]*trace[Yd, Adj[Yd]] - 
   6*mHu2*\[Lambda]*conj[\[Lambda]]*trace[Yd, Adj[Yd]] - 
   6*ms2*\[Lambda]*conj[\[Lambda]]*trace[Yd, Adj[Yd]] - 
   6*conj[T[\[Lambda]]]*T[\[Lambda]]*trace[Yd, Adj[Yd]] - 
   4*mHd2*\[Lambda]*conj[\[Lambda]]*trace[Ye, Adj[Ye]] - 
   2*mHu2*\[Lambda]*conj[\[Lambda]]*trace[Ye, Adj[Ye]] - 
   2*ms2*\[Lambda]*conj[\[Lambda]]*trace[Ye, Adj[Ye]] - 
   2*conj[T[\[Lambda]]]*T[\[Lambda]]*trace[Ye, Adj[Ye]] + 
   (8*g1^2*mHu2*trace[Yu, Adj[Yu]])/5 - (3*g1p^2*mHu2*trace[Yu, Adj[Yu]])/5 + 
   32*g3^2*mHu2*trace[Yu, Adj[Yu]] + 64*g3^2*MassG*conj[MassG]*
    trace[Yu, Adj[Yu]] - 6*mHd2*\[Lambda]*conj[\[Lambda]]*
    trace[\[Kappa], Adj[\[Kappa]]] - 6*mHu2*\[Lambda]*conj[\[Lambda]]*
    trace[\[Kappa], Adj[\[Kappa]]] - 12*ms2*\[Lambda]*conj[\[Lambda]]*
    trace[\[Kappa], Adj[\[Kappa]]] - 6*conj[T[\[Lambda]]]*T[\[Lambda]]*
    trace[\[Kappa], Adj[\[Kappa]]] - 6*\[Lambda]*conj[T[\[Lambda]]]*
    trace[Adj[Yd], T[Yd]] - 2*\[Lambda]*conj[T[\[Lambda]]]*
    trace[Adj[Ye], T[Ye]] + (g1^2*conj[MassB]*(783*g1^2*MassB + 
      48*g1p^2*MassB + 90*g2^2*MassB + 24*g1p^2*MassBp + 45*g2^2*MassWB + 
      80*MassB*trace[Yu, Adj[Yu]] - 40*trace[Adj[Yu], T[Yu]]))/25 - 
   32*g3^2*conj[MassG]*trace[Adj[Yu], T[Yu]] + 
   (g1p^2*conj[MassBp]*(48*g1^2*MassB + 96*g1^2*MassBp + 771*g1p^2*MassBp + 
      120*g2^2*MassBp + 60*g2^2*MassWB + 6*g1p^2*MassBp*QS^2 + 
      5*(5 + QS^2)*conj[\[Lambda]]*(2*MassBp*\[Lambda] - T[\[Lambda]]) - 
      60*MassBp*trace[Yu, Adj[Yu]] + 30*trace[Adj[Yu], T[Yu]]))/50 - 
   6*\[Lambda]*conj[T[\[Lambda]]]*trace[Adj[\[Kappa]], T[\[Kappa]]] - 
   6*conj[\[Lambda]]*T[\[Lambda]]*trace[conj[T[Yd]], Tp[Yd]] - 
   6*\[Lambda]*conj[\[Lambda]]*trace[conj[T[Yd]], Tp[T[Yd]]] - 
   2*conj[\[Lambda]]*T[\[Lambda]]*trace[conj[T[Ye]], Tp[Ye]] - 
   2*\[Lambda]*conj[\[Lambda]]*trace[conj[T[Ye]], Tp[T[Ye]]] - 
   (8*g1^2*MassB*trace[conj[T[Yu]], Tp[Yu]])/5 + 
   (3*g1p^2*MassBp*trace[conj[T[Yu]], Tp[Yu]])/5 - 
   32*g3^2*MassG*trace[conj[T[Yu]], Tp[Yu]] + 
   (8*g1^2*trace[conj[T[Yu]], Tp[T[Yu]]])/5 - 
   (3*g1p^2*trace[conj[T[Yu]], Tp[T[Yu]]])/5 + 
   32*g3^2*trace[conj[T[Yu]], Tp[T[Yu]]] - 6*conj[\[Lambda]]*T[\[Lambda]]*
    trace[conj[T[\[Kappa]]], Tp[\[Kappa]]] - 6*\[Lambda]*conj[\[Lambda]]*
    trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] - 6*\[Lambda]*conj[\[Lambda]]*
    trace[md2, Yd, Adj[Yd]] - 2*\[Lambda]*conj[\[Lambda]]*
    trace[me2, Ye, Adj[Ye]] - 2*\[Lambda]*conj[\[Lambda]]*
    trace[ml2, Adj[Ye], Ye] - 6*\[Lambda]*conj[\[Lambda]]*
    trace[mq2, Adj[Yd], Yd] + (8*g1^2*trace[mq2, Adj[Yu], Yu])/5 - 
   (3*g1p^2*trace[mq2, Adj[Yu], Yu])/5 + 32*g3^2*trace[mq2, Adj[Yu], Yu] + 
   (8*g1^2*trace[mu2, Yu, Adj[Yu]])/5 - (3*g1p^2*trace[mu2, Yu, Adj[Yu]])/5 + 
   32*g3^2*trace[mu2, Yu, Adj[Yu]] - 6*\[Lambda]*conj[\[Lambda]]*
    trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] - 6*\[Lambda]*conj[\[Lambda]]*
    trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]] - 
   6*mHd2*trace[Yd, Adj[Yu], Yu, Adj[Yd]] - 
   6*mHu2*trace[Yd, Adj[Yu], Yu, Adj[Yd]] - 
   6*trace[Yd, Adj[Yu], T[Yu], Adj[T[Yd]]] - 
   6*trace[Yd, Adj[T[Yu]], T[Yu], Adj[Yd]] - 
   6*trace[Yu, Adj[Yd], T[Yd], Adj[T[Yu]]] - 
   36*mHu2*trace[Yu, Adj[Yu], Yu, Adj[Yu]] - 
   36*trace[Yu, Adj[Yu], T[Yu], Adj[T[Yu]]] - 
   6*trace[Yu, Adj[T[Yd]], T[Yd], Adj[Yu]] - 
   36*trace[Yu, Adj[T[Yu]], T[Yu], Adj[Yu]] - 
   6*trace[md2, Yd, Adj[Yu], Yu, Adj[Yd]] - 
   6*trace[mq2, Adj[Yd], Yd, Adj[Yu], Yu] - 
   6*trace[mq2, Adj[Yu], Yu, Adj[Yd], Yd] - 
   36*trace[mq2, Adj[Yu], Yu, Adj[Yu], Yu] - 
   6*trace[mu2, Yu, Adj[Yd], Yd, Adj[Yu]] - 
   36*trace[mu2, Yu, Adj[Yu], Yu, Adj[Yu]]}, 
 {md2[i1, i2], (-8*g1^2*MassB*conj[MassB]*Kronecker[i1, i2])/15 - 
   (4*g1p^2*MassBp*conj[MassBp]*Kronecker[i1, i2])/5 - 
   (32*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/3 + 
   (2*g1*Kronecker[i1, i2]*Tr1[1])/Sqrt[15] + Sqrt[2/5]*g1p*Kronecker[i1, i2]*
    Tr1[4] + 4*mHd2*MatMul[Yd, Adj[Yd]][i1, i2] + 
   4*MatMul[T[Yd], Adj[T[Yd]]][i1, i2] + 2*MatMul[md2, Yd, Adj[Yd]][i1, i2] + 
   4*MatMul[Yd, mq2, Adj[Yd]][i1, i2] + 2*MatMul[Yd, Adj[Yd], md2][i1, i2], 
  (64*g1^2*g3^2*MassB*conj[MassG]*Kronecker[i1, i2])/45 + 
   (32*g1p^2*g3^2*MassBp*conj[MassG]*Kronecker[i1, i2])/15 + 
   (128*g1^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/45 + 
   (64*g1p^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/15 + 
   (160*g3^4*MassG*conj[MassG]*Kronecker[i1, i2])/3 + 
   (32*g3^4*Kronecker[i1, i2]*Tr2[3])/3 + 
   (8*g1^2*Kronecker[i1, i2]*Tr2U1[1, 1])/15 + 
   (4*Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[1, 4])/5 + 
   (4*Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[4, 1])/5 + 
   (4*g1p^2*Kronecker[i1, i2]*Tr2U1[4, 4])/5 + 
   (8*g1*Kronecker[i1, i2]*Tr3[1])/Sqrt[15] + 
   4*Sqrt[2/5]*g1p*Kronecker[i1, i2]*Tr3[4] + 
   (4*g1^2*mHd2*MatMul[Yd, Adj[Yd]][i1, i2])/5 + 
   (6*g1p^2*mHd2*MatMul[Yd, Adj[Yd]][i1, i2])/5 + 
   12*g2^2*mHd2*MatMul[Yd, Adj[Yd]][i1, i2] + 24*g2^2*MassWB*conj[MassWB]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 8*mHd2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 4*mHu2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 4*ms2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 4*conj[T[\[Lambda]]]*T[\[Lambda]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 24*mHd2*trace[Yd, Adj[Yd]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 8*mHd2*trace[Ye, Adj[Ye]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 12*trace[conj[T[Yd]], Tp[T[Yd]]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 4*trace[conj[T[Ye]], Tp[T[Ye]]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 12*trace[md2, Yd, Adj[Yd]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 4*trace[me2, Ye, Adj[Ye]]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 4*trace[ml2, Adj[Ye], Ye]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 12*trace[mq2, Adj[Yd], Yd]*
    MatMul[Yd, Adj[Yd]][i1, i2] - 
   (4*g1^2*MassB*MatMul[Yd, Adj[T[Yd]]][i1, i2])/5 - 
   (6*g1p^2*MassBp*MatMul[Yd, Adj[T[Yd]]][i1, i2])/5 - 
   12*g2^2*MassWB*MatMul[Yd, Adj[T[Yd]]][i1, i2] - 
   4*conj[\[Lambda]]*T[\[Lambda]]*MatMul[Yd, Adj[T[Yd]]][i1, i2] - 
   12*trace[Adj[Yd], T[Yd]]*MatMul[Yd, Adj[T[Yd]]][i1, i2] - 
   4*trace[Adj[Ye], T[Ye]]*MatMul[Yd, Adj[T[Yd]]][i1, i2] + 
   (g1p^2*conj[MassBp]*((-56*g1^2*(MassB + 2*MassBp) + 
        320*g3^2*(2*MassBp + MassG) + 9*g1p^2*MassBp*(257 + 2*QS^2))*
       Kronecker[i1, i2] + 180*(2*MassBp*MatMul[Yd, Adj[Yd]][i1, i2] - 
        MatMul[T[Yd], Adj[Yd]][i1, i2])))/150 + 
   (4*g1^2*conj[MassB]*((768*g1^2*MassB - 21*g1p^2*(2*MassB + MassBp) + 
        80*g3^2*(2*MassB + MassG))*Kronecker[i1, i2] + 
      90*MassB*MatMul[Yd, Adj[Yd]][i1, i2] - 
      45*MatMul[T[Yd], Adj[Yd]][i1, i2]))/225 - 
   12*g2^2*conj[MassWB]*MatMul[T[Yd], Adj[Yd]][i1, i2] - 
   4*\[Lambda]*conj[T[\[Lambda]]]*MatMul[T[Yd], Adj[Yd]][i1, i2] - 
   12*trace[conj[T[Yd]], Tp[Yd]]*MatMul[T[Yd], Adj[Yd]][i1, i2] - 
   4*trace[conj[T[Ye]], Tp[Ye]]*MatMul[T[Yd], Adj[Yd]][i1, i2] + 
   (4*g1^2*MatMul[T[Yd], Adj[T[Yd]]][i1, i2])/5 + 
   (6*g1p^2*MatMul[T[Yd], Adj[T[Yd]]][i1, i2])/5 + 
   12*g2^2*MatMul[T[Yd], Adj[T[Yd]]][i1, i2] - 4*\[Lambda]*conj[\[Lambda]]*
    MatMul[T[Yd], Adj[T[Yd]]][i1, i2] - 12*trace[Yd, Adj[Yd]]*
    MatMul[T[Yd], Adj[T[Yd]]][i1, i2] - 4*trace[Ye, Adj[Ye]]*
    MatMul[T[Yd], Adj[T[Yd]]][i1, i2] + 
   (2*g1^2*MatMul[md2, Yd, Adj[Yd]][i1, i2])/5 + 
   (3*g1p^2*MatMul[md2, Yd, Adj[Yd]][i1, i2])/5 + 
   6*g2^2*MatMul[md2, Yd, Adj[Yd]][i1, i2] - 2*\[Lambda]*conj[\[Lambda]]*
    MatMul[md2, Yd, Adj[Yd]][i1, i2] - 6*trace[Yd, Adj[Yd]]*
    MatMul[md2, Yd, Adj[Yd]][i1, i2] - 2*trace[Ye, Adj[Ye]]*
    MatMul[md2, Yd, Adj[Yd]][i1, i2] + 
   (4*g1^2*MatMul[Yd, mq2, Adj[Yd]][i1, i2])/5 + 
   (6*g1p^2*MatMul[Yd, mq2, Adj[Yd]][i1, i2])/5 + 
   12*g2^2*MatMul[Yd, mq2, Adj[Yd]][i1, i2] - 4*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yd, mq2, Adj[Yd]][i1, i2] - 12*trace[Yd, Adj[Yd]]*
    MatMul[Yd, mq2, Adj[Yd]][i1, i2] - 4*trace[Ye, Adj[Ye]]*
    MatMul[Yd, mq2, Adj[Yd]][i1, i2] + 
   (2*g1^2*MatMul[Yd, Adj[Yd], md2][i1, i2])/5 + 
   (3*g1p^2*MatMul[Yd, Adj[Yd], md2][i1, i2])/5 + 
   6*g2^2*MatMul[Yd, Adj[Yd], md2][i1, i2] - 2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yd, Adj[Yd], md2][i1, i2] - 6*trace[Yd, Adj[Yd]]*
    MatMul[Yd, Adj[Yd], md2][i1, i2] - 2*trace[Ye, Adj[Ye]]*
    MatMul[Yd, Adj[Yd], md2][i1, i2] - 
   8*mHd2*MatMul[Yd, Adj[Yd], Yd, Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, Adj[Yd], T[Yd], Adj[T[Yd]]][i1, i2] - 
   4*mHd2*MatMul[Yd, Adj[Yu], Yu, Adj[Yd]][i1, i2] - 
   4*mHu2*MatMul[Yd, Adj[Yu], Yu, Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, Adj[Yu], T[Yu], Adj[T[Yd]]][i1, i2] - 
   4*MatMul[Yd, Adj[T[Yd]], T[Yd], Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, Adj[T[Yu]], T[Yu], Adj[Yd]][i1, i2] - 
   4*MatMul[T[Yd], Adj[Yd], Yd, Adj[T[Yd]]][i1, i2] - 
   4*MatMul[T[Yd], Adj[Yu], Yu, Adj[T[Yd]]][i1, i2] - 
   4*MatMul[T[Yd], Adj[T[Yd]], Yd, Adj[Yd]][i1, i2] - 
   4*MatMul[T[Yd], Adj[T[Yu]], Yu, Adj[Yd]][i1, i2] - 
   2*MatMul[md2, Yd, Adj[Yd], Yd, Adj[Yd]][i1, i2] - 
   2*MatMul[md2, Yd, Adj[Yu], Yu, Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, mq2, Adj[Yd], Yd, Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, mq2, Adj[Yu], Yu, Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, Adj[Yd], md2, Yd, Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, Adj[Yd], Yd, mq2, Adj[Yd]][i1, i2] - 
   2*MatMul[Yd, Adj[Yd], Yd, Adj[Yd], md2][i1, i2] - 
   4*MatMul[Yd, Adj[Yu], mu2, Yu, Adj[Yd]][i1, i2] - 
   4*MatMul[Yd, Adj[Yu], Yu, mq2, Adj[Yd]][i1, i2] - 
   2*MatMul[Yd, Adj[Yu], Yu, Adj[Yd], md2][i1, i2]}, 
 {mu2[i1, i2], (-32*g1^2*MassB*conj[MassB]*Kronecker[i1, i2])/15 - 
   (g1p^2*MassBp*conj[MassBp]*Kronecker[i1, i2])/5 - 
   (32*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/3 - 
   (4*g1*Kronecker[i1, i2]*Tr1[1])/Sqrt[15] + (g1p*Kronecker[i1, i2]*Tr1[4])/
    Sqrt[10] + 4*mHu2*MatMul[Yu, Adj[Yu]][i1, i2] + 
   4*MatMul[T[Yu], Adj[T[Yu]]][i1, i2] + 2*MatMul[mu2, Yu, Adj[Yu]][i1, i2] + 
   4*MatMul[Yu, mq2, Adj[Yu]][i1, i2] + 2*MatMul[Yu, Adj[Yu], mu2][i1, i2], 
  (256*g1^2*g3^2*MassB*conj[MassG]*Kronecker[i1, i2])/45 + 
   (8*g1p^2*g3^2*MassBp*conj[MassG]*Kronecker[i1, i2])/15 + 
   (512*g1^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/45 + 
   (16*g1p^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/15 + 
   (160*g3^4*MassG*conj[MassG]*Kronecker[i1, i2])/3 + 
   (32*g3^4*Kronecker[i1, i2]*Tr2[3])/3 + 
   (32*g1^2*Kronecker[i1, i2]*Tr2U1[1, 1])/15 - 
   (4*Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[1, 4])/5 - 
   (4*Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[4, 1])/5 + 
   (g1p^2*Kronecker[i1, i2]*Tr2U1[4, 4])/5 - (16*g1*Kronecker[i1, i2]*Tr3[1])/
    Sqrt[15] + 2*Sqrt[2/5]*g1p*Kronecker[i1, i2]*Tr3[4] - 
   (4*g1^2*mHu2*MatMul[Yu, Adj[Yu]][i1, i2])/5 + 
   (4*g1p^2*mHu2*MatMul[Yu, Adj[Yu]][i1, i2])/5 + 
   12*g2^2*mHu2*MatMul[Yu, Adj[Yu]][i1, i2] + 24*g2^2*MassWB*conj[MassWB]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 4*mHd2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 8*mHu2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 4*ms2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 4*conj[T[\[Lambda]]]*T[\[Lambda]]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 24*mHu2*trace[Yu, Adj[Yu]]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 12*trace[conj[T[Yu]], Tp[T[Yu]]]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 12*trace[mq2, Adj[Yu], Yu]*
    MatMul[Yu, Adj[Yu]][i1, i2] - 12*trace[mu2, Yu, Adj[Yu]]*
    MatMul[Yu, Adj[Yu]][i1, i2] + 
   (4*g1^2*MassB*MatMul[Yu, Adj[T[Yu]]][i1, i2])/5 - 
   (4*g1p^2*MassBp*MatMul[Yu, Adj[T[Yu]]][i1, i2])/5 - 
   12*g2^2*MassWB*MatMul[Yu, Adj[T[Yu]]][i1, i2] - 
   4*conj[\[Lambda]]*T[\[Lambda]]*MatMul[Yu, Adj[T[Yu]]][i1, i2] - 
   12*trace[Adj[Yu], T[Yu]]*MatMul[Yu, Adj[T[Yu]]][i1, i2] + 
   (g1p^2*conj[MassBp]*((352*g1^2*(MassB + 2*MassBp) + 
        320*g3^2*(2*MassBp + MassG) + 9*g1p^2*MassBp*(251 + 2*QS^2))*
       Kronecker[i1, i2] + 480*(2*MassBp*MatMul[Yu, Adj[Yu]][i1, i2] - 
        MatMul[T[Yu], Adj[Yu]][i1, i2])))/600 - 
   12*g2^2*conj[MassWB]*MatMul[T[Yu], Adj[Yu]][i1, i2] - 
   4*\[Lambda]*conj[T[\[Lambda]]]*MatMul[T[Yu], Adj[Yu]][i1, i2] - 
   12*trace[conj[T[Yu]], Tp[Yu]]*MatMul[T[Yu], Adj[Yu]][i1, i2] + 
   (4*g1^2*conj[MassB]*((3216*g1^2*MassB + 33*g1p^2*(2*MassB + MassBp) + 
        320*g3^2*(2*MassB + MassG))*Kronecker[i1, i2] + 
      45*(-2*MassB*MatMul[Yu, Adj[Yu]][i1, i2] + MatMul[T[Yu], Adj[Yu]][i1, 
         i2])))/225 - (4*g1^2*MatMul[T[Yu], Adj[T[Yu]]][i1, i2])/5 + 
   (4*g1p^2*MatMul[T[Yu], Adj[T[Yu]]][i1, i2])/5 + 
   12*g2^2*MatMul[T[Yu], Adj[T[Yu]]][i1, i2] - 4*\[Lambda]*conj[\[Lambda]]*
    MatMul[T[Yu], Adj[T[Yu]]][i1, i2] - 12*trace[Yu, Adj[Yu]]*
    MatMul[T[Yu], Adj[T[Yu]]][i1, i2] - 
   (2*g1^2*MatMul[mu2, Yu, Adj[Yu]][i1, i2])/5 + 
   (2*g1p^2*MatMul[mu2, Yu, Adj[Yu]][i1, i2])/5 + 
   6*g2^2*MatMul[mu2, Yu, Adj[Yu]][i1, i2] - 2*\[Lambda]*conj[\[Lambda]]*
    MatMul[mu2, Yu, Adj[Yu]][i1, i2] - 6*trace[Yu, Adj[Yu]]*
    MatMul[mu2, Yu, Adj[Yu]][i1, i2] - 
   (4*g1^2*MatMul[Yu, mq2, Adj[Yu]][i1, i2])/5 + 
   (4*g1p^2*MatMul[Yu, mq2, Adj[Yu]][i1, i2])/5 + 
   12*g2^2*MatMul[Yu, mq2, Adj[Yu]][i1, i2] - 4*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yu, mq2, Adj[Yu]][i1, i2] - 12*trace[Yu, Adj[Yu]]*
    MatMul[Yu, mq2, Adj[Yu]][i1, i2] - 
   (2*g1^2*MatMul[Yu, Adj[Yu], mu2][i1, i2])/5 + 
   (2*g1p^2*MatMul[Yu, Adj[Yu], mu2][i1, i2])/5 + 
   6*g2^2*MatMul[Yu, Adj[Yu], mu2][i1, i2] - 2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Yu, Adj[Yu], mu2][i1, i2] - 6*trace[Yu, Adj[Yu]]*
    MatMul[Yu, Adj[Yu], mu2][i1, i2] - 
   4*mHd2*MatMul[Yu, Adj[Yd], Yd, Adj[Yu]][i1, i2] - 
   4*mHu2*MatMul[Yu, Adj[Yd], Yd, Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, Adj[Yd], T[Yd], Adj[T[Yu]]][i1, i2] - 
   8*mHu2*MatMul[Yu, Adj[Yu], Yu, Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, Adj[Yu], T[Yu], Adj[T[Yu]]][i1, i2] - 
   4*MatMul[Yu, Adj[T[Yd]], T[Yd], Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, Adj[T[Yu]], T[Yu], Adj[Yu]][i1, i2] - 
   4*MatMul[T[Yu], Adj[Yd], Yd, Adj[T[Yu]]][i1, i2] - 
   4*MatMul[T[Yu], Adj[Yu], Yu, Adj[T[Yu]]][i1, i2] - 
   4*MatMul[T[Yu], Adj[T[Yd]], Yd, Adj[Yu]][i1, i2] - 
   4*MatMul[T[Yu], Adj[T[Yu]], Yu, Adj[Yu]][i1, i2] - 
   2*MatMul[mu2, Yu, Adj[Yd], Yd, Adj[Yu]][i1, i2] - 
   2*MatMul[mu2, Yu, Adj[Yu], Yu, Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, mq2, Adj[Yd], Yd, Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, mq2, Adj[Yu], Yu, Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, Adj[Yd], md2, Yd, Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, Adj[Yd], Yd, mq2, Adj[Yu]][i1, i2] - 
   2*MatMul[Yu, Adj[Yd], Yd, Adj[Yu], mu2][i1, i2] - 
   4*MatMul[Yu, Adj[Yu], mu2, Yu, Adj[Yu]][i1, i2] - 
   4*MatMul[Yu, Adj[Yu], Yu, mq2, Adj[Yu]][i1, i2] - 
   2*MatMul[Yu, Adj[Yu], Yu, Adj[Yu], mu2][i1, i2]}, 
 {me2[i1, i2], (-24*g1^2*MassB*conj[MassB]*Kronecker[i1, i2])/5 - 
   (g1p^2*MassBp*conj[MassBp]*Kronecker[i1, i2])/5 + 
   2*Sqrt[3/5]*g1*Kronecker[i1, i2]*Tr1[1] + (g1p*Kronecker[i1, i2]*Tr1[4])/
    Sqrt[10] + 4*mHd2*MatMul[Ye, Adj[Ye]][i1, i2] + 
   4*MatMul[T[Ye], Adj[T[Ye]]][i1, i2] + 2*MatMul[me2, Ye, Adj[Ye]][i1, i2] + 
   4*MatMul[Ye, ml2, Adj[Ye]][i1, i2] + 2*MatMul[Ye, Adj[Ye], me2][i1, i2], 
  (24*g1^2*Kronecker[i1, i2]*Tr2U1[1, 1])/5 + 
   (2*Sqrt[6]*g1*g1p*Kronecker[i1, i2]*Tr2U1[1, 4])/5 + 
   (2*Sqrt[6]*g1*g1p*Kronecker[i1, i2]*Tr2U1[4, 1])/5 + 
   (g1p^2*Kronecker[i1, i2]*Tr2U1[4, 4])/5 + 8*Sqrt[3/5]*g1*Kronecker[i1, i2]*
    Tr3[1] + 2*Sqrt[2/5]*g1p*Kronecker[i1, i2]*Tr3[4] - 
   (12*g1^2*mHd2*MatMul[Ye, Adj[Ye]][i1, i2])/5 + 
   (12*g1p^2*mHd2*MatMul[Ye, Adj[Ye]][i1, i2])/5 + 
   12*g2^2*mHd2*MatMul[Ye, Adj[Ye]][i1, i2] + 24*g2^2*MassWB*conj[MassWB]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 8*mHd2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 4*mHu2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 4*ms2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 4*conj[T[\[Lambda]]]*T[\[Lambda]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 24*mHd2*trace[Yd, Adj[Yd]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 8*mHd2*trace[Ye, Adj[Ye]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 12*trace[conj[T[Yd]], Tp[T[Yd]]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 4*trace[conj[T[Ye]], Tp[T[Ye]]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 12*trace[md2, Yd, Adj[Yd]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 4*trace[me2, Ye, Adj[Ye]]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 4*trace[ml2, Adj[Ye], Ye]*
    MatMul[Ye, Adj[Ye]][i1, i2] - 12*trace[mq2, Adj[Yd], Yd]*
    MatMul[Ye, Adj[Ye]][i1, i2] + 
   (12*g1^2*MassB*MatMul[Ye, Adj[T[Ye]]][i1, i2])/5 - 
   (12*g1p^2*MassBp*MatMul[Ye, Adj[T[Ye]]][i1, i2])/5 - 
   12*g2^2*MassWB*MatMul[Ye, Adj[T[Ye]]][i1, i2] - 
   4*conj[\[Lambda]]*T[\[Lambda]]*MatMul[Ye, Adj[T[Ye]]][i1, i2] - 
   12*trace[Adj[Yd], T[Yd]]*MatMul[Ye, Adj[T[Ye]]][i1, i2] - 
   4*trace[Adj[Ye], T[Ye]]*MatMul[Ye, Adj[T[Ye]]][i1, i2] + 
   (3*g1p^2*conj[MassBp]*((-32*g1^2*(MassB + 2*MassBp) + 
        g1p^2*MassBp*(251 + 2*QS^2))*Kronecker[i1, i2] + 
      160*(2*MassBp*MatMul[Ye, Adj[Ye]][i1, i2] - MatMul[T[Ye], Adj[Ye]][i1, 
         i2])))/200 - 12*g2^2*conj[MassWB]*MatMul[T[Ye], Adj[Ye]][i1, i2] - 
   4*\[Lambda]*conj[T[\[Lambda]]]*MatMul[T[Ye], Adj[Ye]][i1, i2] - 
   12*trace[conj[T[Yd]], Tp[Yd]]*MatMul[T[Ye], Adj[Ye]][i1, i2] - 
   4*trace[conj[T[Ye]], Tp[Ye]]*MatMul[T[Ye], Adj[Ye]][i1, i2] + 
   (12*g1^2*conj[MassB]*((288*g1^2*MassB - g1p^2*(2*MassB + MassBp))*
       Kronecker[i1, i2] + 5*(-2*MassB*MatMul[Ye, Adj[Ye]][i1, i2] + 
        MatMul[T[Ye], Adj[Ye]][i1, i2])))/25 - 
   (12*g1^2*MatMul[T[Ye], Adj[T[Ye]]][i1, i2])/5 + 
   (12*g1p^2*MatMul[T[Ye], Adj[T[Ye]]][i1, i2])/5 + 
   12*g2^2*MatMul[T[Ye], Adj[T[Ye]]][i1, i2] - 4*\[Lambda]*conj[\[Lambda]]*
    MatMul[T[Ye], Adj[T[Ye]]][i1, i2] - 12*trace[Yd, Adj[Yd]]*
    MatMul[T[Ye], Adj[T[Ye]]][i1, i2] - 4*trace[Ye, Adj[Ye]]*
    MatMul[T[Ye], Adj[T[Ye]]][i1, i2] - 
   (6*g1^2*MatMul[me2, Ye, Adj[Ye]][i1, i2])/5 + 
   (6*g1p^2*MatMul[me2, Ye, Adj[Ye]][i1, i2])/5 + 
   6*g2^2*MatMul[me2, Ye, Adj[Ye]][i1, i2] - 2*\[Lambda]*conj[\[Lambda]]*
    MatMul[me2, Ye, Adj[Ye]][i1, i2] - 6*trace[Yd, Adj[Yd]]*
    MatMul[me2, Ye, Adj[Ye]][i1, i2] - 2*trace[Ye, Adj[Ye]]*
    MatMul[me2, Ye, Adj[Ye]][i1, i2] - 
   (12*g1^2*MatMul[Ye, ml2, Adj[Ye]][i1, i2])/5 + 
   (12*g1p^2*MatMul[Ye, ml2, Adj[Ye]][i1, i2])/5 + 
   12*g2^2*MatMul[Ye, ml2, Adj[Ye]][i1, i2] - 4*\[Lambda]*conj[\[Lambda]]*
    MatMul[Ye, ml2, Adj[Ye]][i1, i2] - 12*trace[Yd, Adj[Yd]]*
    MatMul[Ye, ml2, Adj[Ye]][i1, i2] - 4*trace[Ye, Adj[Ye]]*
    MatMul[Ye, ml2, Adj[Ye]][i1, i2] - 
   (6*g1^2*MatMul[Ye, Adj[Ye], me2][i1, i2])/5 + 
   (6*g1p^2*MatMul[Ye, Adj[Ye], me2][i1, i2])/5 + 
   6*g2^2*MatMul[Ye, Adj[Ye], me2][i1, i2] - 2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Ye, Adj[Ye], me2][i1, i2] - 6*trace[Yd, Adj[Yd]]*
    MatMul[Ye, Adj[Ye], me2][i1, i2] - 2*trace[Ye, Adj[Ye]]*
    MatMul[Ye, Adj[Ye], me2][i1, i2] - 
   8*mHd2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye]][i1, i2] - 
   4*MatMul[Ye, Adj[Ye], T[Ye], Adj[T[Ye]]][i1, i2] - 
   4*MatMul[Ye, Adj[T[Ye]], T[Ye], Adj[Ye]][i1, i2] - 
   4*MatMul[T[Ye], Adj[Ye], Ye, Adj[T[Ye]]][i1, i2] - 
   4*MatMul[T[Ye], Adj[T[Ye]], Ye, Adj[Ye]][i1, i2] - 
   2*MatMul[me2, Ye, Adj[Ye], Ye, Adj[Ye]][i1, i2] - 
   4*MatMul[Ye, ml2, Adj[Ye], Ye, Adj[Ye]][i1, i2] - 
   4*MatMul[Ye, Adj[Ye], me2, Ye, Adj[Ye]][i1, i2] - 
   4*MatMul[Ye, Adj[Ye], Ye, ml2, Adj[Ye]][i1, i2] - 
   2*MatMul[Ye, Adj[Ye], Ye, Adj[Ye], me2][i1, i2]}, 
 {ms2, -(g1p^2*MassBp*QS^2*conj[MassBp])/5 + 4*(mHd2 + mHu2 + ms2)*\[Lambda]*
    conj[\[Lambda]] + 2*mphi2*\[Sigma]*conj[\[Sigma]] + 
   2*ms2*\[Sigma]*conj[\[Sigma]] + 2*msbar2*\[Sigma]*conj[\[Sigma]] + 
   4*conj[T[\[Lambda]]]*T[\[Lambda]] + 2*conj[T[\[Sigma]]]*T[\[Sigma]] + 
   (g1p*QS*Tr1[4])/Sqrt[10] + 6*ms2*trace[\[Kappa], Adj[\[Kappa]]] + 
   6*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] + 
   6*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] + 
   6*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]], 
  (g1p^2*conj[MassBp]*(-40*(-13 + QS^2)*conj[\[Lambda]]*
       (2*MassBp*\[Lambda] - T[\[Lambda]]) + 
      3*(g1p^2*MassBp*QS^2*(249 + 4*QS^2) - 40*MassBp*(-13 + QS^2)*
         trace[\[Kappa], Adj[\[Kappa]]] + 20*(-13 + QS^2)*
         trace[Adj[\[Kappa]], T[\[Kappa]]])) + 
    20*(-160*(mHd2 + mHu2 + ms2)*\[Lambda]^2*conj[\[Lambda]]^2 - 
      80*mphi2*\[Sigma]^2*conj[\[Sigma]]^2 - 80*ms2*\[Sigma]^2*
       conj[\[Sigma]]^2 - 80*msbar2*\[Sigma]^2*conj[\[Sigma]]^2 - 
      24*g1^2*MassB*\[Lambda]*conj[T[\[Lambda]]] - 26*g1p^2*MassBp*\[Lambda]*
       conj[T[\[Lambda]]] - 120*g2^2*MassWB*\[Lambda]*conj[T[\[Lambda]]] + 
      2*g1p^2*MassBp*QS^2*\[Lambda]*conj[T[\[Lambda]]] - 
      40*\[Sigma]*conj[\[Sigma]]*conj[T[\[Kappa]Pr]]*T[\[Kappa]Pr] + 
      24*g1^2*conj[T[\[Lambda]]]*T[\[Lambda]] + 26*g1p^2*conj[T[\[Lambda]]]*
       T[\[Lambda]] + 120*g2^2*conj[T[\[Lambda]]]*T[\[Lambda]] - 
      2*g1p^2*QS^2*conj[T[\[Lambda]]]*T[\[Lambda]] - 
      40*\[Kappa]Pr*conj[\[Sigma]]*conj[T[\[Kappa]Pr]]*T[\[Sigma]] - 
      160*\[Sigma]*conj[\[Sigma]]*conj[T[\[Sigma]]]*T[\[Sigma]] - 
      40*conj[\[Kappa]Pr]*((4*mphi2 + ms2 + msbar2)*\[Kappa]Pr*\[Sigma]*
         conj[\[Sigma]] + conj[T[\[Sigma]]]*(\[Sigma]*T[\[Kappa]Pr] + 
          \[Kappa]Pr*T[\[Sigma]])) + 2*g1p^2*QS^2*Tr2U1[4, 4] + 
      4*Sqrt[10]*g1p*QS*Tr3[4] - 120*conj[T[\[Lambda]]]*T[\[Lambda]]*
       trace[Yd, Adj[Yd]] - 40*conj[T[\[Lambda]]]*T[\[Lambda]]*
       trace[Ye, Adj[Ye]] - 120*conj[T[\[Lambda]]]*T[\[Lambda]]*
       trace[Yu, Adj[Yu]] + 16*g1^2*ms2*trace[\[Kappa], Adj[\[Kappa]]] + 
      39*g1p^2*ms2*trace[\[Kappa], Adj[\[Kappa]]] + 
      320*g3^2*ms2*trace[\[Kappa], Adj[\[Kappa]]] - 
      3*g1p^2*ms2*QS^2*trace[\[Kappa], Adj[\[Kappa]]] + 
      32*g1^2*MassB*conj[MassB]*trace[\[Kappa], Adj[\[Kappa]]] + 
      640*g3^2*MassG*conj[MassG]*trace[\[Kappa], Adj[\[Kappa]]] - 
      120*\[Lambda]*conj[T[\[Lambda]]]*trace[Adj[Yd], T[Yd]] - 
      40*\[Lambda]*conj[T[\[Lambda]]]*trace[Adj[Ye], T[Ye]] - 
      120*\[Lambda]*conj[T[\[Lambda]]]*trace[Adj[Yu], T[Yu]] - 
      16*g1^2*conj[MassB]*trace[Adj[\[Kappa]], T[\[Kappa]]] - 
      320*g3^2*conj[MassG]*trace[Adj[\[Kappa]], T[\[Kappa]]] - 
      16*g1^2*MassB*trace[conj[T[\[Kappa]]], Tp[\[Kappa]]] - 
      39*g1p^2*MassBp*trace[conj[T[\[Kappa]]], Tp[\[Kappa]]] - 
      320*g3^2*MassG*trace[conj[T[\[Kappa]]], Tp[\[Kappa]]] + 
      3*g1p^2*MassBp*QS^2*trace[conj[T[\[Kappa]]], Tp[\[Kappa]]] + 
      16*g1^2*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] + 
      39*g1p^2*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] + 
      320*g3^2*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] - 
      3*g1p^2*QS^2*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] + 
      2*conj[\[Lambda]]*(12*g1^2*mHd2*\[Lambda] + 13*g1p^2*mHd2*\[Lambda] + 
        60*g2^2*mHd2*\[Lambda] + 12*g1^2*mHu2*\[Lambda] + 
        13*g1p^2*mHu2*\[Lambda] + 60*g2^2*mHu2*\[Lambda] + 
        12*g1^2*ms2*\[Lambda] + 13*g1p^2*ms2*\[Lambda] + 
        60*g2^2*ms2*\[Lambda] - g1p^2*mHd2*QS^2*\[Lambda] - 
        g1p^2*mHu2*QS^2*\[Lambda] - g1p^2*ms2*QS^2*\[Lambda] + 
        12*g1^2*conj[MassB]*(2*MassB*\[Lambda] - T[\[Lambda]]) + 
        60*g2^2*conj[MassWB]*(2*MassWB*\[Lambda] - T[\[Lambda]]) - 
        160*\[Lambda]*conj[T[\[Lambda]]]*T[\[Lambda]] - 
        120*mHd2*\[Lambda]*trace[Yd, Adj[Yd]] - 60*mHu2*\[Lambda]*
         trace[Yd, Adj[Yd]] - 60*ms2*\[Lambda]*trace[Yd, Adj[Yd]] - 
        40*mHd2*\[Lambda]*trace[Ye, Adj[Ye]] - 20*mHu2*\[Lambda]*
         trace[Ye, Adj[Ye]] - 20*ms2*\[Lambda]*trace[Ye, Adj[Ye]] - 
        60*mHd2*\[Lambda]*trace[Yu, Adj[Yu]] - 120*mHu2*\[Lambda]*
         trace[Yu, Adj[Yu]] - 60*ms2*\[Lambda]*trace[Yu, Adj[Yu]] - 
        60*T[\[Lambda]]*trace[conj[T[Yd]], Tp[Yd]] - 
        60*\[Lambda]*trace[conj[T[Yd]], Tp[T[Yd]]] - 20*T[\[Lambda]]*
         trace[conj[T[Ye]], Tp[Ye]] - 20*\[Lambda]*trace[conj[T[Ye]], 
          Tp[T[Ye]]] - 60*T[\[Lambda]]*trace[conj[T[Yu]], Tp[Yu]] - 
        60*\[Lambda]*trace[conj[T[Yu]], Tp[T[Yu]]] - 
        60*\[Lambda]*trace[md2, Yd, Adj[Yd]] - 20*\[Lambda]*
         trace[me2, Ye, Adj[Ye]] - 20*\[Lambda]*trace[ml2, Adj[Ye], Ye] - 
        60*\[Lambda]*trace[mq2, Adj[Yd], Yd] - 60*\[Lambda]*
         trace[mq2, Adj[Yu], Yu] - 60*\[Lambda]*trace[mu2, Yu, Adj[Yu]]) + 
      16*g1^2*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] + 
      39*g1p^2*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] + 
      320*g3^2*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] - 
      3*g1p^2*QS^2*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] + 
      16*g1^2*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]] + 
      39*g1p^2*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]] + 
      320*g3^2*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]] - 
      3*g1p^2*QS^2*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]] - 
      240*ms2*trace[\[Kappa], Adj[\[Kappa]], \[Kappa], Adj[\[Kappa]]] - 
      240*trace[\[Kappa], Adj[\[Kappa]], T[\[Kappa]], Adj[T[\[Kappa]]]] - 
      240*trace[\[Kappa], Adj[T[\[Kappa]]], T[\[Kappa]], Adj[\[Kappa]]] - 
      120*trace[\[Kappa], Adj[\[Kappa]], \[Kappa], Adj[\[Kappa]], 
        conj[mDx2]] - 120*trace[\[Kappa], Adj[\[Kappa]], \[Kappa], 
        conj[mDxbar2], Adj[\[Kappa]]] - 120*trace[\[Kappa], Adj[\[Kappa]], 
        conj[mDx2], \[Kappa], Adj[\[Kappa]]] - 
      120*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]], \[Kappa], 
        Adj[\[Kappa]]]))/200}, {msbar2, -(g1p^2*MassBp*QS^2*conj[MassBp])/5 + 
   2*(mphi2 + ms2 + msbar2)*\[Sigma]*conj[\[Sigma]] + 
   2*conj[T[\[Sigma]]]*T[\[Sigma]] - (g1p*QS*Tr1[4])/Sqrt[10], 
  (3*g1p^4*MassBp*QS^2*(249 + 4*QS^2)*conj[MassBp])/200 - 
   4*mHd2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   4*mHu2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   4*mphi2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   8*ms2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   4*msbar2*\[Lambda]*\[Sigma]*conj[\[Lambda]]*conj[\[Sigma]] - 
   8*mphi2*\[Sigma]^2*conj[\[Sigma]]^2 - 8*ms2*\[Sigma]^2*conj[\[Sigma]]^2 - 
   8*msbar2*\[Sigma]^2*conj[\[Sigma]]^2 - 4*\[Sigma]*conj[\[Sigma]]*
    conj[T[\[Kappa]Pr]]*T[\[Kappa]Pr] - 4*\[Sigma]*conj[\[Sigma]]*
    conj[T[\[Lambda]]]*T[\[Lambda]] - 4*\[Sigma]*conj[\[Lambda]]*
    conj[T[\[Sigma]]]*T[\[Lambda]] - 4*\[Kappa]Pr*conj[\[Sigma]]*
    conj[T[\[Kappa]Pr]]*T[\[Sigma]] - 4*\[Lambda]*conj[\[Sigma]]*
    conj[T[\[Lambda]]]*T[\[Sigma]] - 4*\[Lambda]*conj[\[Lambda]]*
    conj[T[\[Sigma]]]*T[\[Sigma]] - 16*\[Sigma]*conj[\[Sigma]]*
    conj[T[\[Sigma]]]*T[\[Sigma]] - 4*conj[\[Kappa]Pr]*
    ((4*mphi2 + ms2 + msbar2)*\[Kappa]Pr*\[Sigma]*conj[\[Sigma]] + 
     conj[T[\[Sigma]]]*(\[Sigma]*T[\[Kappa]Pr] + \[Kappa]Pr*T[\[Sigma]])) + 
   (g1p^2*QS^2*Tr2U1[4, 4])/5 - 2*Sqrt[2/5]*g1p*QS*Tr3[4] - 
   6*mphi2*\[Sigma]*conj[\[Sigma]]*trace[\[Kappa], Adj[\[Kappa]]] - 
   12*ms2*\[Sigma]*conj[\[Sigma]]*trace[\[Kappa], Adj[\[Kappa]]] - 
   6*msbar2*\[Sigma]*conj[\[Sigma]]*trace[\[Kappa], Adj[\[Kappa]]] - 
   6*conj[T[\[Sigma]]]*T[\[Sigma]]*trace[\[Kappa], Adj[\[Kappa]]] - 
   6*\[Sigma]*conj[T[\[Sigma]]]*trace[Adj[\[Kappa]], T[\[Kappa]]] - 
   6*conj[\[Sigma]]*T[\[Sigma]]*trace[conj[T[\[Kappa]]], Tp[\[Kappa]]] - 
   6*\[Sigma]*conj[\[Sigma]]*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] - 
   6*\[Sigma]*conj[\[Sigma]]*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] - 
   6*\[Sigma]*conj[\[Sigma]]*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]]}, 
 {mDx2[i1, i2], (-8*g1^2*MassB*conj[MassB]*Kronecker[i1, i2])/15 - 
   (4*g1p^2*MassBp*conj[MassBp]*Kronecker[i1, i2])/5 - 
   (32*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/3 - 
   (2*g1*Kronecker[i1, i2]*Tr1[1])/Sqrt[15] - Sqrt[2/5]*g1p*Kronecker[i1, i2]*
    Tr1[4] + 2*ms2*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] + 
   2*MatMul[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]][i1, i2] + 
   MatMul[mDx2, conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] + 
   2*MatMul[conj[\[Kappa]], mDxbar2, Tp[\[Kappa]]][i1, i2] + 
   MatMul[conj[\[Kappa]], Tp[\[Kappa]], mDx2][i1, i2], 
  (4*g1^2*(768*g1^2*MassB - 21*g1p^2*(2*MassB + MassBp) + 
      80*g3^2*(2*MassB + MassG))*conj[MassB]*Kronecker[i1, i2])/225 + 
   (64*g1^2*g3^2*MassB*conj[MassG]*Kronecker[i1, i2])/45 + 
   (32*g1p^2*g3^2*MassBp*conj[MassG]*Kronecker[i1, i2])/15 + 
   (128*g1^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/45 + 
   (64*g1p^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/15 + 
   (160*g3^4*MassG*conj[MassG]*Kronecker[i1, i2])/3 + 
   (32*g3^4*Kronecker[i1, i2]*Tr2[3])/3 + 
   (8*g1^2*Kronecker[i1, i2]*Tr2U1[1, 1])/15 + 
   (4*Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[1, 4])/5 + 
   (4*Sqrt[2/3]*g1*g1p*Kronecker[i1, i2]*Tr2U1[4, 1])/5 + 
   (4*g1p^2*Kronecker[i1, i2]*Tr2U1[4, 4])/5 - 
   (8*g1*Kronecker[i1, i2]*Tr3[1])/Sqrt[15] - 
   4*Sqrt[2/5]*g1p*Kronecker[i1, i2]*Tr3[4] + 
   (g1p^2*ms2*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2])/2 + 
   (g1p^2*ms2*QS^2*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2])/10 - 
   4*mHd2*\[Lambda]*conj[\[Lambda]]*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, 
     i2] - 4*mHu2*\[Lambda]*conj[\[Lambda]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
   8*ms2*\[Lambda]*conj[\[Lambda]]*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, 
     i2] - 2*mphi2*\[Sigma]*conj[\[Sigma]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
   4*ms2*\[Sigma]*conj[\[Sigma]]*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, 
     i2] - 2*msbar2*\[Sigma]*conj[\[Sigma]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
   4*conj[T[\[Lambda]]]*T[\[Lambda]]*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, 
     i2] - 2*conj[T[\[Sigma]]]*T[\[Sigma]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
   12*ms2*trace[\[Kappa], Adj[\[Kappa]]]*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][
     i1, i2] - 6*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
   6*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
   6*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] + 
   (g1p^2*conj[MassBp]*((-56*g1^2*(MassB + 2*MassBp) + 
        320*g3^2*(2*MassBp + MassG) + 9*g1p^2*MassBp*(257 + 2*QS^2))*
       Kronecker[i1, i2] + 15*(5 + QS^2)*
       (2*MassBp*MatMul[conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
        MatMul[conj[\[Kappa]], Tp[T[\[Kappa]]]][i1, i2])))/150 - 
   4*\[Lambda]*conj[T[\[Lambda]]]*MatMul[conj[\[Kappa]], Tp[T[\[Kappa]]]][i1, 
     i2] - 2*\[Sigma]*conj[T[\[Sigma]]]*
    MatMul[conj[\[Kappa]], Tp[T[\[Kappa]]]][i1, i2] - 
   6*trace[conj[T[\[Kappa]]], Tp[\[Kappa]]]*
    MatMul[conj[\[Kappa]], Tp[T[\[Kappa]]]][i1, i2] - 
   (g1p^2*MassBp*MatMul[conj[T[\[Kappa]]], Tp[\[Kappa]]][i1, i2])/2 - 
   (g1p^2*MassBp*QS^2*MatMul[conj[T[\[Kappa]]], Tp[\[Kappa]]][i1, i2])/10 - 
   4*conj[\[Lambda]]*T[\[Lambda]]*MatMul[conj[T[\[Kappa]]], Tp[\[Kappa]]][i1, 
     i2] - 2*conj[\[Sigma]]*T[\[Sigma]]*
    MatMul[conj[T[\[Kappa]]], Tp[\[Kappa]]][i1, i2] - 
   6*trace[Adj[\[Kappa]], T[\[Kappa]]]*
    MatMul[conj[T[\[Kappa]]], Tp[\[Kappa]]][i1, i2] + 
   (g1p^2*MatMul[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]][i1, i2])/2 + 
   (g1p^2*QS^2*MatMul[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]][i1, i2])/10 - 
   4*\[Lambda]*conj[\[Lambda]]*MatMul[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]][i1, 
     i2] - 2*\[Sigma]*conj[\[Sigma]]*MatMul[conj[T[\[Kappa]]], 
      Tp[T[\[Kappa]]]][i1, i2] - 6*trace[\[Kappa], Adj[\[Kappa]]]*
    MatMul[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]][i1, i2] + 
   (g1p^2*MatMul[mDx2, conj[\[Kappa]], Tp[\[Kappa]]][i1, i2])/4 + 
   (g1p^2*QS^2*MatMul[mDx2, conj[\[Kappa]], Tp[\[Kappa]]][i1, i2])/20 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[mDx2, conj[\[Kappa]], Tp[\[Kappa]]][i1, 
     i2] - \[Sigma]*conj[\[Sigma]]*MatMul[mDx2, conj[\[Kappa]], Tp[\[Kappa]]][
     i1, i2] - 3*trace[\[Kappa], Adj[\[Kappa]]]*
    MatMul[mDx2, conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] + 
   (g1p^2*MatMul[conj[\[Kappa]], mDxbar2, Tp[\[Kappa]]][i1, i2])/2 + 
   (g1p^2*QS^2*MatMul[conj[\[Kappa]], mDxbar2, Tp[\[Kappa]]][i1, i2])/10 - 
   4*\[Lambda]*conj[\[Lambda]]*MatMul[conj[\[Kappa]], mDxbar2, Tp[\[Kappa]]][
     i1, i2] - 2*\[Sigma]*conj[\[Sigma]]*
    MatMul[conj[\[Kappa]], mDxbar2, Tp[\[Kappa]]][i1, i2] - 
   6*trace[\[Kappa], Adj[\[Kappa]]]*MatMul[conj[\[Kappa]], mDxbar2, 
      Tp[\[Kappa]]][i1, i2] + 
   (g1p^2*MatMul[conj[\[Kappa]], Tp[\[Kappa]], mDx2][i1, i2])/4 + 
   (g1p^2*QS^2*MatMul[conj[\[Kappa]], Tp[\[Kappa]], mDx2][i1, i2])/20 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[conj[\[Kappa]], Tp[\[Kappa]], mDx2][i1, 
     i2] - \[Sigma]*conj[\[Sigma]]*MatMul[conj[\[Kappa]], Tp[\[Kappa]], mDx2][
     i1, i2] - 3*trace[\[Kappa], Adj[\[Kappa]]]*
    MatMul[conj[\[Kappa]], Tp[\[Kappa]], mDx2][i1, i2] - 
   4*ms2*MatMul[conj[\[Kappa]], Tp[\[Kappa]], conj[\[Kappa]], Tp[\[Kappa]]][
     i1, i2] - 2*MatMul[conj[\[Kappa]], Tp[\[Kappa]], conj[T[\[Kappa]]], 
      Tp[T[\[Kappa]]]][i1, i2] - 
   2*MatMul[conj[\[Kappa]], Tp[T[\[Kappa]]], conj[T[\[Kappa]]], Tp[\[Kappa]]][
     i1, i2] - 2*MatMul[conj[T[\[Kappa]]], Tp[\[Kappa]], conj[\[Kappa]], 
      Tp[T[\[Kappa]]]][i1, i2] - 
   2*MatMul[conj[T[\[Kappa]]], Tp[T[\[Kappa]]], conj[\[Kappa]], Tp[\[Kappa]]][
     i1, i2] - MatMul[mDx2, conj[\[Kappa]], Tp[\[Kappa]], conj[\[Kappa]], 
     Tp[\[Kappa]]][i1, i2] - 2*MatMul[conj[\[Kappa]], mDxbar2, Tp[\[Kappa]], 
      conj[\[Kappa]], Tp[\[Kappa]]][i1, i2] - 
   2*MatMul[conj[\[Kappa]], Tp[\[Kappa]], mDx2, conj[\[Kappa]], Tp[\[Kappa]]][
     i1, i2] - 2*MatMul[conj[\[Kappa]], Tp[\[Kappa]], conj[\[Kappa]], 
      mDxbar2, Tp[\[Kappa]]][i1, i2] - 
   MatMul[conj[\[Kappa]], Tp[\[Kappa]], conj[\[Kappa]], Tp[\[Kappa]], mDx2][
    i1, i2]}, {mDxbar2[i1, i2], (-8*g1^2*MassB*conj[MassB]*Kronecker[i1, i2])/
    15 - (9*g1p^2*MassBp*conj[MassBp]*Kronecker[i1, i2])/5 - 
   (32*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/3 + 
   (2*g1*Kronecker[i1, i2]*Tr1[1])/Sqrt[15] - 
   (3*g1p*Kronecker[i1, i2]*Tr1[4])/Sqrt[10] + 
   2*ms2*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] + 
   2*MatMul[Tp[T[\[Kappa]]], conj[T[\[Kappa]]]][i1, i2] + 
   MatMul[mDxbar2, Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] + 
   2*MatMul[Tp[\[Kappa]], mDx2, conj[\[Kappa]]][i1, i2] + 
   MatMul[Tp[\[Kappa]], conj[\[Kappa]], mDxbar2][i1, i2], 
  (8*g1^2*(384*g1^2*MassB + 27*g1p^2*(2*MassB + MassBp) + 
      40*g3^2*(2*MassB + MassG))*conj[MassB]*Kronecker[i1, i2])/225 + 
   (64*g1^2*g3^2*MassB*conj[MassG]*Kronecker[i1, i2])/45 + 
   (24*g1p^2*g3^2*MassBp*conj[MassG]*Kronecker[i1, i2])/5 + 
   (128*g1^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/45 + 
   (48*g1p^2*g3^2*MassG*conj[MassG]*Kronecker[i1, i2])/5 + 
   (160*g3^4*MassG*conj[MassG]*Kronecker[i1, i2])/3 + 
   (32*g3^4*Kronecker[i1, i2]*Tr2[3])/3 + 
   (8*g1^2*Kronecker[i1, i2]*Tr2U1[1, 1])/15 - 
   (2*Sqrt[6]*g1*g1p*Kronecker[i1, i2]*Tr2U1[1, 4])/5 - 
   (2*Sqrt[6]*g1*g1p*Kronecker[i1, i2]*Tr2U1[4, 1])/5 + 
   (9*g1p^2*Kronecker[i1, i2]*Tr2U1[4, 4])/5 + 
   (8*g1*Kronecker[i1, i2]*Tr3[1])/Sqrt[15] - 
   6*Sqrt[2/5]*g1p*Kronecker[i1, i2]*Tr3[4] - 
   (g1p^2*ms2*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2])/2 + 
   (g1p^2*ms2*QS^2*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2])/10 - 
   4*mHd2*\[Lambda]*conj[\[Lambda]]*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, 
     i2] - 4*mHu2*\[Lambda]*conj[\[Lambda]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
   8*ms2*\[Lambda]*conj[\[Lambda]]*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, 
     i2] - 2*mphi2*\[Sigma]*conj[\[Sigma]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
   4*ms2*\[Sigma]*conj[\[Sigma]]*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, 
     i2] - 2*msbar2*\[Sigma]*conj[\[Sigma]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
   4*conj[T[\[Lambda]]]*T[\[Lambda]]*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, 
     i2] - 2*conj[T[\[Sigma]]]*T[\[Sigma]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
   12*ms2*trace[\[Kappa], Adj[\[Kappa]]]*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][
     i1, i2] - 6*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
   6*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
   6*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] + 
   (g1p^2*MassBp*MatMul[Tp[\[Kappa]], conj[T[\[Kappa]]]][i1, i2])/2 - 
   (g1p^2*MassBp*QS^2*MatMul[Tp[\[Kappa]], conj[T[\[Kappa]]]][i1, i2])/10 - 
   4*conj[\[Lambda]]*T[\[Lambda]]*MatMul[Tp[\[Kappa]], conj[T[\[Kappa]]]][i1, 
     i2] - 2*conj[\[Sigma]]*T[\[Sigma]]*
    MatMul[Tp[\[Kappa]], conj[T[\[Kappa]]]][i1, i2] - 
   6*trace[Adj[\[Kappa]], T[\[Kappa]]]*
    MatMul[Tp[\[Kappa]], conj[T[\[Kappa]]]][i1, i2] + 
   (g1p^2*conj[MassBp]*(3*(64*g1^2*(MassB + 2*MassBp) + 
        320*g3^2*(2*MassBp + MassG) + 9*g1p^2*MassBp*(267 + 2*QS^2))*
       Kronecker[i1, i2] + 20*(-5 + QS^2)*
       (2*MassBp*MatMul[Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
        MatMul[Tp[T[\[Kappa]]], conj[\[Kappa]]][i1, i2])))/200 - 
   4*\[Lambda]*conj[T[\[Lambda]]]*MatMul[Tp[T[\[Kappa]]], conj[\[Kappa]]][i1, 
     i2] - 2*\[Sigma]*conj[T[\[Sigma]]]*
    MatMul[Tp[T[\[Kappa]]], conj[\[Kappa]]][i1, i2] - 
   6*trace[conj[T[\[Kappa]]], Tp[\[Kappa]]]*
    MatMul[Tp[T[\[Kappa]]], conj[\[Kappa]]][i1, i2] - 
   (g1p^2*MatMul[Tp[T[\[Kappa]]], conj[T[\[Kappa]]]][i1, i2])/2 + 
   (g1p^2*QS^2*MatMul[Tp[T[\[Kappa]]], conj[T[\[Kappa]]]][i1, i2])/10 - 
   4*\[Lambda]*conj[\[Lambda]]*MatMul[Tp[T[\[Kappa]]], conj[T[\[Kappa]]]][i1, 
     i2] - 2*\[Sigma]*conj[\[Sigma]]*
    MatMul[Tp[T[\[Kappa]]], conj[T[\[Kappa]]]][i1, i2] - 
   6*trace[\[Kappa], Adj[\[Kappa]]]*
    MatMul[Tp[T[\[Kappa]]], conj[T[\[Kappa]]]][i1, i2] - 
   (g1p^2*MatMul[mDxbar2, Tp[\[Kappa]], conj[\[Kappa]]][i1, i2])/4 + 
   (g1p^2*QS^2*MatMul[mDxbar2, Tp[\[Kappa]], conj[\[Kappa]]][i1, i2])/20 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[mDxbar2, Tp[\[Kappa]], conj[\[Kappa]]][
     i1, i2] - \[Sigma]*conj[\[Sigma]]*
    MatMul[mDxbar2, Tp[\[Kappa]], conj[\[Kappa]]][i1, i2] - 
   3*trace[\[Kappa], Adj[\[Kappa]]]*MatMul[mDxbar2, Tp[\[Kappa]], 
      conj[\[Kappa]]][i1, i2] - 
   (g1p^2*MatMul[Tp[\[Kappa]], mDx2, conj[\[Kappa]]][i1, i2])/2 + 
   (g1p^2*QS^2*MatMul[Tp[\[Kappa]], mDx2, conj[\[Kappa]]][i1, i2])/10 - 
   4*\[Lambda]*conj[\[Lambda]]*MatMul[Tp[\[Kappa]], mDx2, conj[\[Kappa]]][i1, 
     i2] - 2*\[Sigma]*conj[\[Sigma]]*MatMul[Tp[\[Kappa]], mDx2, 
      conj[\[Kappa]]][i1, i2] - 6*trace[\[Kappa], Adj[\[Kappa]]]*
    MatMul[Tp[\[Kappa]], mDx2, conj[\[Kappa]]][i1, i2] - 
   (g1p^2*MatMul[Tp[\[Kappa]], conj[\[Kappa]], mDxbar2][i1, i2])/4 + 
   (g1p^2*QS^2*MatMul[Tp[\[Kappa]], conj[\[Kappa]], mDxbar2][i1, i2])/20 - 
   2*\[Lambda]*conj[\[Lambda]]*MatMul[Tp[\[Kappa]], conj[\[Kappa]], mDxbar2][
     i1, i2] - \[Sigma]*conj[\[Sigma]]*MatMul[Tp[\[Kappa]], conj[\[Kappa]], 
      mDxbar2][i1, i2] - 3*trace[\[Kappa], Adj[\[Kappa]]]*
    MatMul[Tp[\[Kappa]], conj[\[Kappa]], mDxbar2][i1, i2] - 
   4*ms2*MatMul[Tp[\[Kappa]], conj[\[Kappa]], Tp[\[Kappa]], conj[\[Kappa]]][
     i1, i2] - 2*MatMul[Tp[\[Kappa]], conj[\[Kappa]], Tp[T[\[Kappa]]], 
      conj[T[\[Kappa]]]][i1, i2] - 
   2*MatMul[Tp[\[Kappa]], conj[T[\[Kappa]]], Tp[T[\[Kappa]]], conj[\[Kappa]]][
     i1, i2] - 2*MatMul[Tp[T[\[Kappa]]], conj[\[Kappa]], Tp[\[Kappa]], 
      conj[T[\[Kappa]]]][i1, i2] - 
   2*MatMul[Tp[T[\[Kappa]]], conj[T[\[Kappa]]], Tp[\[Kappa]], conj[\[Kappa]]][
     i1, i2] - MatMul[mDxbar2, Tp[\[Kappa]], conj[\[Kappa]], Tp[\[Kappa]], 
     conj[\[Kappa]]][i1, i2] - 
   2*MatMul[Tp[\[Kappa]], mDx2, conj[\[Kappa]], Tp[\[Kappa]], conj[\[Kappa]]][
     i1, i2] - 2*MatMul[Tp[\[Kappa]], conj[\[Kappa]], mDxbar2, Tp[\[Kappa]], 
      conj[\[Kappa]]][i1, i2] - 
   2*MatMul[Tp[\[Kappa]], conj[\[Kappa]], Tp[\[Kappa]], mDx2, conj[\[Kappa]]][
     i1, i2] - MatMul[Tp[\[Kappa]], conj[\[Kappa]], Tp[\[Kappa]], 
     conj[\[Kappa]], mDxbar2][i1, i2]}, 
 {mHp2, (-6*g1^2*MassB*conj[MassB] - 4*g1p^2*MassBp*conj[MassBp] - 
    30*g2^2*MassWB*conj[MassWB] - Sqrt[15]*g1*Tr1[1] + Sqrt[10]*g1p*Tr1[4])/
   5, (6*g1^2*(261*g1^2*MassB + 8*g1p^2*(2*MassB + MassBp) + 
      15*g2^2*(2*MassB + MassWB))*conj[MassB] + 
    3*g1p^2*(16*g1^2*(MassB + 2*MassBp) + 20*g2^2*(2*MassBp + MassWB) + 
      g1p^2*MassBp*(257 + 2*QS^2))*conj[MassBp] + 
    10*(3*g2^2*(85*g2^2*MassWB + 3*g1^2*(MassB + 2*MassWB) + 
        2*g1p^2*(MassBp + 2*MassWB))*conj[MassWB] + 30*g2^4*Tr2[2] + 
      6*g1^2*Tr2U1[1, 1] - 2*Sqrt[6]*g1*g1p*Tr2U1[1, 4] - 
      2*Sqrt[6]*g1*g1p*Tr2U1[4, 1] + 4*g1p^2*Tr2U1[4, 4] - 
      4*Sqrt[15]*g1*Tr3[1] + 4*Sqrt[10]*g1p*Tr3[4]))/50}, 
 {mHpbar2, (-6*g1^2*MassB*conj[MassB] - 4*g1p^2*MassBp*conj[MassBp] - 
    30*g2^2*MassWB*conj[MassWB] + Sqrt[15]*g1*Tr1[1] - Sqrt[10]*g1p*Tr1[4])/
   5, (6*g1^2*(261*g1^2*MassB + 8*g1p^2*(2*MassB + MassBp) + 
      15*g2^2*(2*MassB + MassWB))*conj[MassB] + 
    3*g1p^2*(16*g1^2*(MassB + 2*MassBp) + 20*g2^2*(2*MassBp + MassWB) + 
      g1p^2*MassBp*(257 + 2*QS^2))*conj[MassBp] + 
    10*(3*g2^2*(85*g2^2*MassWB + 3*g1^2*(MassB + 2*MassWB) + 
        2*g1p^2*(MassBp + 2*MassWB))*conj[MassWB] + 30*g2^4*Tr2[2] + 
      6*g1^2*Tr2U1[1, 1] - 2*Sqrt[6]*g1*g1p*Tr2U1[1, 4] - 
      2*Sqrt[6]*g1*g1p*Tr2U1[4, 1] + 4*g1p^2*Tr2U1[4, 4] + 
      4*Sqrt[15]*g1*Tr3[1] - 4*Sqrt[10]*g1p*Tr3[4]))/50}, 
 {mphi2, 2*(6*mphi2*\[Kappa]Pr*conj[\[Kappa]Pr] + (mphi2 + ms2 + msbar2)*
     \[Sigma]*conj[\[Sigma]] + 2*conj[T[\[Kappa]Pr]]*T[\[Kappa]Pr] + 
    conj[T[\[Sigma]]]*T[\[Sigma]]), 
  -96*mphi2*\[Kappa]Pr^2*conj[\[Kappa]Pr]^2 - 8*(mphi2 + ms2 + msbar2)*
    \[Sigma]^2*conj[\[Sigma]]^2 - 8*conj[\[Kappa]Pr]*
    ((4*mphi2 + ms2 + msbar2)*\[Kappa]Pr*\[Sigma]*conj[\[Sigma]] + 
     8*\[Kappa]Pr*conj[T[\[Kappa]Pr]]*T[\[Kappa]Pr] + 
     conj[T[\[Sigma]]]*(\[Sigma]*T[\[Kappa]Pr] + \[Kappa]Pr*T[\[Sigma]])) - 
   (conj[T[\[Sigma]]]*(g1p^2*MassBp*QS^2*\[Sigma] - g1p^2*QS^2*T[\[Sigma]] + 
      20*conj[\[Lambda]]*(\[Sigma]*T[\[Lambda]] + \[Lambda]*T[\[Sigma]]) + 
      30*T[\[Sigma]]*trace[\[Kappa], Adj[\[Kappa]]] + 
      30*\[Sigma]*trace[Adj[\[Kappa]], T[\[Kappa]]]))/5 + 
   (conj[\[Sigma]]*(g1p^2*mphi2*QS^2*\[Sigma] + g1p^2*ms2*QS^2*\[Sigma] + 
      g1p^2*msbar2*QS^2*\[Sigma] - 20*(mHd2 + mHu2 + mphi2 + 2*ms2 + msbar2)*
       \[Lambda]*\[Sigma]*conj[\[Lambda]] - 40*\[Sigma]*conj[T[\[Kappa]Pr]]*
       T[\[Kappa]Pr] - 20*\[Sigma]*conj[T[\[Lambda]]]*T[\[Lambda]] + 
      g1p^2*QS^2*conj[MassBp]*(2*MassBp*\[Sigma] - T[\[Sigma]]) - 
      40*\[Kappa]Pr*conj[T[\[Kappa]Pr]]*T[\[Sigma]] - 
      20*\[Lambda]*conj[T[\[Lambda]]]*T[\[Sigma]] - 
      80*\[Sigma]*conj[T[\[Sigma]]]*T[\[Sigma]] - 30*mphi2*\[Sigma]*
       trace[\[Kappa], Adj[\[Kappa]]] - 60*ms2*\[Sigma]*
       trace[\[Kappa], Adj[\[Kappa]]] - 30*msbar2*\[Sigma]*
       trace[\[Kappa], Adj[\[Kappa]]] - 30*T[\[Sigma]]*
       trace[conj[T[\[Kappa]]], Tp[\[Kappa]]] - 
      30*\[Sigma]*trace[conj[T[\[Kappa]]], Tp[T[\[Kappa]]]] - 
      30*\[Sigma]*trace[\[Kappa], Adj[\[Kappa]], conj[mDx2]] - 
      30*\[Sigma]*trace[\[Kappa], conj[mDxbar2], Adj[\[Kappa]]]))/5}}
