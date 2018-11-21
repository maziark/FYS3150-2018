double analyticEavg(double Z, double beta){
    return (16*exp(-8*beta) - 16*exp(8*beta))/Z;
}

double analyticMavg(double Z, double beta){
    return (8*exp(8*beta) + 16)/Z;
}

double analyticE2avg(double Z, double beta){
    return 256*cosh(8*beta)/Z;
}

double analyticM2avg(double Z, double beta){
    return (32 + 32*exp(8*beta))/Z;
}

double analyticSusceptibility(double analytic_Mavg, double analytic_M2avg, double temp){
    return (analytic_M2avg - analytic_Mavg*analytic_Mavg)/temp;
}

double analyticHeatCapacity(double analytic_Eavg, double analytic_E2avg, double temp){
    return (analytic_E2avg - analytic_Eavg*analytic_Eavg)/(temp*temp);
}

