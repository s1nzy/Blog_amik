$PARAM
// typical parameters
//TVKA = 0, 
TVCL = 0.062, TVVC = 0.918, TVVP1 = 0.918, TVQ1 = 0.062

// covariates 
BW = 1000, CW = 1000, PNA = 1, NSAID = 0

//$CMT GUT CENT P1
$CMT CENT P1


$MAIN
//double KA = TVKA * exp(ETA(1));
double CL = TVCL * pow(BW/1750, 1.1) * (1 + 0.128*(PNA/2)) * exp(ETA(1));
double VC = TVVC * pow(CW/1760,0.929) * exp(ETA(2));
double VP1 = TVVP1* pow(CW/1760,0.929) * exp(ETA(3));
double Q1 = 0.785 * TVQ1 * pow(BW/1750, 1.1) * (1 + 0.128*(PNA/2))*exp(ETA(4));

if(NSAID == 1) { 
    CL = CL * 0.838;
    Q1 = Q1 * 0.838;
  }

double k10 = CL/VC;

double k12 = Q1/VC;
double k21 = Q1/VP1;


$OMEGA 0 0 0 0

$SIGMA @labels PROP ADD
0 0

$ODE
//dxdt_GUT = 0;
//dxdt_CENT = KA*GUT -k12*CENT+k21*P1 - k10*CENT;
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
