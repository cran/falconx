#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Print.h>
#include <Rmath.h>
#include <math.h>

// new functions in falcon based on the EM algorithm
SEXP ScanIGSGridCumSumNewC(SEXP ATS, SEXP gridCurS) {
  double *AT = REAL(ATS);
  double *gridCur = REAL(gridCurS);
  long long gridCurLen = length(gridCurS);
  SEXP ATCumSum;
  PROTECT(ATCumSum = allocVector(REALSXP, gridCurLen-1));
  double *ATCumSumPtr = REAL(ATCumSum);
  long long i, j;
  for (i=0; i<gridCurLen-1; i++){
    ATCumSumPtr[i] = 0;
    for (j=gridCur[i]-1; j<gridCur[i+1]-1; j++){
      ATCumSumPtr[i] += AT[j];
    }
  }
  UNPROTECT(1);
  return(ATCumSum);
}

// falcon2: GetP changed to GetC; Lik, LikH, ScanStatNewCompBinom2dEMC, ScanStatRefineCompBinom2dEMC are updated

SEXP GetC(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP WS, SEXP errorS, SEXP maxIterS, SEXP COriS) {
  double *AT = REAL(ATS);
  double *BT = REAL(BTS); 
  double *AN = REAL(ANS);
  double *BN = REAL(BNS);
  // double *sT = REAL(sTS);
  // double *sN = REAL(sNS);
  double *W = REAL(WS);
  // double Tdep = REAL(TdepS)[0];
  // double Ndep = REAL(NdepS)[0];
  double *COri = REAL(COriS);
  double error = REAL(errorS)[0];
  double maxIter = REAL(maxIterS)[0];
  long long len = length(ATS);
  SEXP rS, CS;
  PROTECT(CS = allocVector(REALSXP, 2));
  PROTECT(rS = allocVector(REALSXP, len));
  
  double *r = REAL(rS);
  double *C = REAL(CS);
  double nIter = 0;
  double curError = 1;
  double Ca, Cb, CaOld, CbOld, Ca1, Ca2, Cb1, Cb2;
  Ca = COri[0];
  Cb = COri[1];
  while (curError>error && nIter<maxIter){
    CaOld = Ca;
    CbOld = Cb;
    for (long long i=0; i<len; i++){
      double temp = (AT[i]-BT[i])*log(Cb/Ca) + (AT[i]+AN[i]-BT[i]-BN[i])*log((W[i]*Ca+1)/(W[i]*Cb+1));
      if (temp>100){
        r[i] = exp(-temp);
      }else{
        r[i] = 1.0/(1.0+exp(temp));
      }
        // r[i] = 1.0/(1.0+pow(pb/pa, AT[i]-BT[i])*pow((1-pb)/(1-pa), AN[i]-BN[i]));
    }
    Ca1 = Ca2 = Cb1 = Cb2 = 0.0;
    for (long long i=0; i<len; i++){
      Ca1 += (AT[i]*r[i] + BT[i]*(1-r[i]));
      Ca2 += W[i] * ((AT[i]+AN[i])*r[i] + (BT[i]+BN[i])*(1-r[i])) / (W[i]*Ca+1);
      Cb1 += (AT[i]*(1-r[i]) + BT[i]*r[i]);
      Cb2 += W[i] * ((AT[i]+AN[i])*(1-r[i]) + (BT[i]+BN[i])*r[i]) / (W[i]*Cb+1);
    }
    Ca = Ca1/Ca2;
    Cb = Cb1/Cb2;
    curError = sqrt(pow(Ca-CaOld,2) + pow(Cb-CbOld,2));
    nIter = nIter + 1;
  }
  C[0] = Ca;
  C[1] = Cb;
  
  UNPROTECT(2);
  return(CS);
}

SEXP Lik(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP WS, SEXP CS) {
  double *AT = REAL(ATS);
  double *BT = REAL(BTS);
  double *AN = REAL(ANS);
  double *BN = REAL(BNS);
  double *W = REAL(WS);
  double *C = REAL(CS);
  double Ca = C[0];
  double Cb = C[1];
  long long len = length(ATS);
  SEXP likS;
  PROTECT(likS = allocVector(REALSXP, 1));
  double *lik = REAL(likS);
  double likVal = 0;
  for (long long i=0; i<len; i++){
    // likVal += log(pow(pa,AT[i])*pow(1-pa,AN[i])*pow(pb,BT[i])*pow(1-pb,BN[i]) + pow(pa,BT[i])*pow(1-pa,BN[i])*pow(pb,AT[i])*pow(1-pb,AN[i])); // has problem wen AT too large
    if (Ca*Cb !=0){
      // temp = log(f2/f1)
      double temp = (BT[i]-AT[i])*log(Ca/Cb) + (BN[i]+BT[i]-AN[i]-AT[i])*log((W[i]*Cb+1)/(W[i]*Ca+1));
      if (temp<100){
        // likVal = log(f1) + log(1+f2/f1)
        likVal += AT[i]*log(Ca) - (AN[i]+AT[i])*log(W[i]*Ca+1) + BT[i]*log(Cb) - (BN[i]+BT[i])*log(W[i]*Cb+1) + log(1 + exp(temp));
      }else{
        likVal += AT[i]*log(Ca) - (AN[i]+AT[i])*log(W[i]*Ca+1) + BT[i]*log(Cb) - (BN[i]+BT[i])*log(W[i]*Cb+1) + temp;
      }
    }else{
      if (Ca==0 && AT[i]==0){
        likVal += BT[i]*log(Cb) - (BT[i]+BN[i])*log(W[i]*Cb+1);
      }else if (Ca==0 && BT[i]==0){
        likVal += AT[i]*log(Cb) - (AT[i]+AN[i])*log(W[i]*Cb+1);
      }else if (Cb==0 && AT[i]==0){
        likVal += BT[i]*log(Ca) - (BT[i]+BN[i])*log(W[i]*Ca+1);
      }else if (Cb==0 && BT[i]==0){
        likVal += AT[i]*log(Ca) - (AT[i]+AN[i])*log(W[i]*Ca+1);
      }else{
        likVal += log(0);
      }
    }
  }
  lik[0] = likVal;
  UNPROTECT(1);
  return(likS);
}
 
SEXP LikH(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP WS, SEXP CS) {
  double *AT = REAL(ATS);
  double *BT = REAL(BTS);
  double *AN = REAL(ANS);
  double *BN = REAL(BNS);
  double *W = REAL(WS);
  double *C = REAL(CS);
  double Ca = C[0];
  double Cb = C[1];
  long long len = length(ATS);
  SEXP likh;
  PROTECT(likh = allocVector(REALSXP, 2));
  double *likhPtr = REAL(likh);

  SEXP lik;
  PROTECT(lik = Lik(ATS, BTS, ANS, BNS, WS, CS));
  likhPtr[0] = REAL(lik)[0];

  if (Ca*Cb ==0){
    likhPtr[1] = 0;
  }else{  
    // double H = 0;
    double ta = log(Ca);
    double tb = log(Cb);
    // double h1a = pa; // exp(ta)/(1+exp(ta));
    // double h1b = pb; // exp(tb)/(1+exp(tb));
    // double h2a = pa*(1-pa); // exp(ta)/pow(1+exp(ta),2);
    // double h2b = pb*(1-pb); // exp(tb)/pow(1+exp(tb),2);
    double l2a, l2ab, l2b, f12; // f12 = f1/f2
    l2a = l2ab = l2b = 0;
    for (long long i=0; i<len; i++){
      f12 = exp((AT[i]-BT[i])*(ta-tb) - (AT[i]+AN[i]-BT[i]-BN[i])*log((Cb*W[i]+1)/(Ca*W[i]+1)));
      // f1 = exp(AT[i]*ta - AN[i]*log(1+exp(ta)) + BT[i]*tb - BN[i]*log(1+exp(tb)));
      // f2 = exp(BT[i]*ta - BN[i]*log(1+exp(ta)) + AT[i]*tb - AN[i]*log(1+exp(tb)));
      double h1a = W[i]*Ca/(W[i]*Ca+1);
      double h1b = W[i]*Cb/(W[i]*Cb+1);
      double h2a = W[i]*Ca/pow((W[i]*Ca+1),2);
      double h2b = W[i]*Cb/pow((W[i]*Cb+1),2);
      l2a += pow((AT[i]-BT[i]-(AT[i]+AN[i]-BT[i]-BN[i])*h1a),2)/(1+1/f12)/(f12+1) - h2a*((AT[i]+AN[i])/(1+1/f12)+(BT[i]+BN[i])/(1+f12));
      l2ab += (AT[i]-BT[i]-(AT[i]+AN[i]-BT[i]-BN[i])*h1a)*(BT[i]-AT[i]-(BT[i]+BN[i]-AT[i]-AN[i])*h1b)/(1+1/f12)/(1+f12);
      l2b += pow((BT[i]-AT[i]-(BT[i]+BN[i]-AT[i]-AN[i])*h1b),2)/(1+1/f12)/(1+f12) - h2b*((BT[i]+BN[i])/(1+1/f12)+(AT[i]+AN[i])/(1+f12));
    }
    // printf("l2a: %f, l2b %f, l2ab %f\n", l2a, l2b, l2ab);
    likhPtr[1] = log(l2a*l2b - pow(l2ab,2));
  }
  UNPROTECT(2);
  return(likh);
}

SEXP SubSeq(SEXP ATS, long long begin, long long end) {
  // counting from 0, return ATS[begin, end)
  double *AT = REAL(ATS);
  SEXP newATS;
  PROTECT(newATS = allocVector(REALSXP, end-begin));
  double *newAT = REAL(newATS);
  for (long long i=0; i<(end-begin); i++){
    newAT[i] = AT[begin+i];
  }
  UNPROTECT(1);
  // double *test = REAL(newATS);
  
  return(newATS);
}

SEXP SubSeq2(SEXP ATS, long long begin, long long end) {
  // return the complement of ATS[begin, end)
  double *AT = REAL(ATS);
  long long len = length(ATS);
  SEXP newATS;
  PROTECT(newATS = allocVector(REALSXP, len-end+begin));
  double *newAT = REAL(newATS);
  for (long long i=0; i<begin; i++){
    newAT[i] = AT[i];
  }
  for (long long i=begin; i<len-end+begin; i++){
    newAT[i] = AT[end+i-begin];
  }
  UNPROTECT(1);
  // double *test = REAL(newATS);  
  return(newATS);
}

SEXP ScanStatNewCompBinom2dEMC(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP WS, SEXP errorS, SEXP maxIterS, SEXP CS, SEXP gridCurS, SEXP maxWinS) {
  /* double *AT = REAL(ATS); */
  /* double *BT = REAL(BTS); */
  /* double *AN = REAL(ANS); */
  /* double *BN = REAL(BNS); */
  //  long long ATlen = length(ATS);
  SEXP C0;
  PROTECT(C0 = GetC(ATS, BTS, ANS, BNS, WS, errorS, maxIterS, CS));
  long long maxWin = floor(REAL(maxWinS)[0]);
  double *gridCur = REAL(gridCurS);
  
  long long gridCurLen = length(gridCurS);
  long long gridCurMaxInd = gridCurLen - 1;
  long long rawi, rawj, i, j, jMax, bestWinI, bestWinJ;
  double l0, bestWinR, Rij;
  SEXP Cij, COut, newRes; //, ATij, BTij, ANij, BNij, ATOut, BTOut, ANOut, BNOut;
  SEXP l0S;
  PROTECT(l0S = Lik(ATS, BTS, ANS, BNS, WS, C0));
  l0 = REAL(l0S)[0];
  UNPROTECT(2);  
  PROTECT(newRes = allocMatrix(REALSXP, gridCurMaxInd, 3));
  double *newResPtr = REAL(newRes);
  int newIter = 1;
  
  for (rawi=0; rawi<gridCurMaxInd; rawi++){
    jMax = rawi + maxWin;
    if (jMax > gridCurMaxInd)  jMax = gridCurMaxInd;
    bestWinI = gridCur[rawi];
    bestWinJ = gridCur[jMax];
    bestWinR = 0.0;
    newIter = 1;
    for (rawj=rawi+1; rawj<=jMax; rawj++){
      if (rawj-rawi == gridCurMaxInd) break;
      i = gridCur[rawi];
      j = gridCur[rawj]-1;
      SEXP ATij; PROTECT(ATij = SubSeq(ATS,i,j));
      SEXP BTij; PROTECT(BTij = SubSeq(BTS,i,j));
      SEXP ANij; PROTECT(ANij = SubSeq(ANS,i,j));
      SEXP BNij; PROTECT(BNij = SubSeq(BNS,i,j));
      SEXP Wij; PROTECT(Wij = SubSeq(WS,i,j));
      PROTECT(Cij = GetC(ATij, BTij, ANij, BNij, Wij, errorS, maxIterS, CS));
      SEXP lijS; 
      PROTECT(lijS = Lik(ATij, BTij, ANij, BNij, Wij, Cij));
      Rij = REAL(lijS)[0];
      UNPROTECT(7); // ATij, BTij, ANij, BNij, pij, lij(Orig)
      SEXP ATOut; PROTECT(ATOut = SubSeq2(ATS,i,j));
      SEXP BTOut; PROTECT(BTOut = SubSeq2(BTS,i,j));
      SEXP ANOut; PROTECT(ANOut = SubSeq2(ANS,i,j));
      SEXP BNOut; PROTECT(BNOut = SubSeq2(BNS,i,j));
      SEXP WOut; PROTECT(WOut = SubSeq2(WS,i,j));
      PROTECT(COut = GetC(ATOut, BTOut, ANOut, BNOut, WOut, errorS, maxIterS, CS));
      SEXP lOutS;
      PROTECT(lOutS = Lik(ATOut, BTOut, ANOut, BNOut, WOut, COut));
      Rij += REAL(lOutS)[0];
      UNPROTECT(7); // ATOut, BTOut, ANOut, BNOut, pOut, lOut(Orig)
      if (Rij > bestWinR || newIter == 1) {
        bestWinI = i;
        bestWinJ = j+1;
        bestWinR = Rij;
        // printf("i: %lld, j: %lld, lij:%f, lOut:%f, Rij%f\n", i,j, lij, lOut, Rij);
      }
      newIter = 0;
    }
    bestWinR = bestWinR - l0;
    // if (bestWinR < 0) bestWinR = 0;
    newResPtr[rawi] = bestWinI;
    newResPtr[rawi + gridCurMaxInd] = bestWinJ;
    newResPtr[rawi + 2*gridCurMaxInd] = bestWinR;
  }
  UNPROTECT(1);
  return(newRes);
}
      
SEXP ScanStatRefineCompBinom2dEMC(SEXP ATS, SEXP BTS, SEXP ANS, SEXP BNS, SEXP WS, SEXP errorS, SEXP maxIterS, SEXP CS, SEXP gridCurS, SEXP idLS, SEXP idRS) {
  /* double *AT = REAL(ATS); */
  /* double *BT = REAL(BTS); */
  /* double *AN = REAL(ANS); */
  /* double *BN = REAL(BNS); */
  SEXP C0;
  PROTECT(C0 = GetC(ATS, BTS, ANS, BNS, WS, errorS, maxIterS, CS));
  long long gridCurLen = length(gridCurS);
  long long gridCurMaxInd = gridCurLen-1;
  double *gridCur = REAL(gridCurS);
  double *idL = REAL(idLS);
  double *idR = REAL(idRS);
  long long rawi, rawj, i, j, nRows, bestWinI, bestWinJ, rCt;
  double l0, lij, lOut, bestWinR, Rij;
  SEXP Cij, COut; //, ATij, BTij, ANij, BNij, ATOut, BTOut, ANOut, BNOut;
  int newIter = 1;
  SEXP l0S;
  bestWinI = 0;
  bestWinJ = 0;
  bestWinR = 0;
  PROTECT(l0S = Lik(ATS, BTS, ANS, BNS, WS, C0));
  l0 = REAL(l0S)[0];
  UNPROTECT(2);
  nRows = length(idLS);
  if (idL[nRows-1] == gridCurMaxInd) nRows = nRows-1;
  // printf("nRows: %lld, length(idRS): %lld\n", nRows, (long long) length(idRS));
  
  SEXP newRes;
  PROTECT(newRes = allocMatrix(REALSXP, nRows, 3));
  double * newResPtr = REAL(newRes);
  rCt = 0;
  for (rawi=idL[0]; rawi<=idL[nRows-1]; rawi++){
    newIter = 1;
    for (rawj=idR[0]; rawj<=idR[length(idRS)-1]; rawj++){
      while(rawj < rawi+1){
        rawj++;
      }
      // printf("rawi:%lld, rawj:%lld\n", rawi, rawj);
      if (rawj-rawi == (length(gridCurS)-1)) break;
      i = gridCur[rawi];
      j = gridCur[rawj]-1;
      SEXP ATij; PROTECT(ATij = SubSeq(ATS,i,j));
      SEXP BTij; PROTECT(BTij = SubSeq(BTS,i,j));
      SEXP ANij; PROTECT(ANij = SubSeq(ANS,i,j));
      SEXP BNij; PROTECT(BNij = SubSeq(BNS,i,j));
      // SEXP sTij; PROTECT(sTij = SubSeq(sTS,i,j));
      // SEXP sNij; PROTECT(sNij = SubSeq(sNS,i,j));
      SEXP Wij; PROTECT(Wij = SubSeq(WS,i,j));
      PROTECT(Cij = GetC(ATij, BTij, ANij, BNij, Wij, errorS, maxIterS, CS));
      SEXP lijS; PROTECT(lijS = Lik(ATij, BTij, ANij, BNij, Wij, Cij));
      lij = REAL(lijS)[0];
      UNPROTECT(7);
      Rij = lij;
      SEXP ATOut; PROTECT(ATOut = SubSeq2(ATS,i,j));
      SEXP BTOut; PROTECT(BTOut = SubSeq2(BTS,i,j));
      SEXP ANOut; PROTECT(ANOut = SubSeq2(ANS,i,j));
      SEXP BNOut; PROTECT(BNOut = SubSeq2(BNS,i,j));
      // SEXP sTOut; PROTECT( sTOut = SubSeq2(sTS,i,j));
      // SEXP sNOut; PROTECT( sNOut = SubSeq2(sNS,i,j));
      SEXP WOut; PROTECT( WOut = SubSeq2(WS,i,j));
      PROTECT(COut = GetC(ATOut, BTOut, ANOut, BNOut, WOut, errorS, maxIterS, CS));
      SEXP lOutS; PROTECT(lOutS = Lik(ATOut, BTOut, ANOut, BNOut, WOut, COut));
      lOut = REAL(lOutS)[0];
      UNPROTECT(7);
      Rij += lOut;
      // printf("Rij:%f\n", Rij);
      if (Rij > bestWinR || newIter == 1) {
        bestWinI = i;
        bestWinJ = j+1;
        bestWinR = Rij;
        // printf("i:%lld, j:%lld, Rij:%f\n", i, j, Rij);
     }
      newIter = 0;
    }
    bestWinR = bestWinR - l0;
    // if (bestWinR < 0) bestWinR = 0;
    newResPtr[rCt] = bestWinI;
    newResPtr[rCt + nRows] = bestWinJ;
    newResPtr[rCt + 2*nRows] = bestWinR;
    rCt++;
  }
  UNPROTECT(1);
  return(newRes);
}
