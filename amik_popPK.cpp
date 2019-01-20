$PARAM

// typical parameters

TVCL = 0.0493, TVVC = 0.833, TVVP1 = 0.833, TVQ1 = 0.415

// covariate relationship parameters
EBW = 1.134, SPNA = 0.213, ECW = 0.919, IBU = 0.838

// initial covariates 
BW = 1000, CW = 1000, PNA = 1, NSAID = 0

$CMT CENT P1


$MAIN

double CL = TVCL * pow(BW/1750, EBW) * (1 + SPNA*(PNA/2)) * exp(ETA(1));
double VC = TVVC * pow(CW/1760,ECW) * exp(ETA(2));
double VP1 = TVVP1* pow(CW/1760,ECW) * exp(ETA(3));
double Q1 = TVQ1 * TVCL * pow(BW/1750, EBW) * (1 + SPNA*(PNA/2))*exp(ETA(4));

if(NSAID == 1) { 
    CL = CL * IBU;
    Q1 = Q1 * IBU;
  }

double k10 = CL/VC;

double k12 = Q1/VC;
double k21 = Q1/VP1;


$OMEGA 0 0 0 0

$SIGMA @labels PROP ADD
0 0

$ODE
dxdt_CENT = -k12*CENT+k21*P1 - k10*CENT;
dxdt_P1 = k12*CENT-k21*P1;


$TABLE

double IPRED = CENT/VC;
double DV = IPRED*(1+PROP)+ADD;

while(DV < 0) {
  simeps();
  DV = IPRED*(1+PROP)+ADD;
}


$CAPTURE IPRED DV
