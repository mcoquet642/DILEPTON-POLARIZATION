#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#define M_HBARC 0.197

double alphaEM=1.0/137.0;
double mllSqr=0.0;
double qFSqrSum=1.0/9.0+4.0/9.0+1.0/9.0;

double g=1.94163;

int QUARK_SUPPRESSION=1;

double Nc=3.0;
double Nf=3.0;

double nuG=2.0*(Nc*Nc-1.0);
double nuQ=4.0*Nc*Nf;

double dA=Nc*Nc-1.0;
double dF=Nc;

double CA=Nc;
double CF=(Nc*Nc-1.0)/(2.0*Nc);

// ENERGY DENSITIES OF QUARKS AND GLUONS //
double eGlEq(double T){
    return (M_PI*M_PI/30.0)*nuG*T*T*T*T;
}

double eQuEq(double T){
    return (7.0*M_PI*M_PI/240.0)*nuQ*T*T*T*T;
}

double eEq(double T){
    return eGlEq(T)+eQuEq(T);
}

double eNEq(double T,double qSupp){
    return eGlEq(T)+qSupp*eQuEq(T);
}

// STANDARD RANDOM NUMBER GENERATOR //
double rng(){
    return drand48();
}

#include "HydroAttractor.cpp"
#include "PhaseSpaceDistribution.cpp"
#include "DileptonRates.cpp"


// COMMANDLINE OPTIONS //
#include "IO/cfile.c"

int main(int argc, char **argv) {
    
    // SET COMMANDLINE ARGUMENTS //
    Konfig CommandlineArguments(argc,argv);
    
    // COLLISION PARAMETERS //
    double EtaOverS=0.16; double dNchdEta=1900; double Area=110;
    
    CommandlineArguments.Getval("etas",EtaOverS);
    CommandlineArguments.Getval("Nch",dNchdEta);
    CommandlineArguments.Getval("area",Area);
    CommandlineArguments.Getval("Q",QUARK_SUPPRESSION);
    
    std::cerr << "#CALCULATING DILEPTON PRODUCTION FOR dNchdEta=" << dNchdEta << " Area=" << Area << " fm^2 AND Eta/s=" << EtaOverS << " QUARK SUPPRESION " << QUARK_SUPPRESSION << std::endl;
    
    
    // DILEPTON PARAMTERS //
    double QMin=0.0; double QMax=10.0; int NQBins=100;

    CommandlineArguments.Getval("QMin",QMin);
    CommandlineArguments.Getval("QMax",QMax);
    CommandlineArguments.Getval("NQBins",NQBins);
    
    double qTMin=0.0; double qTMax=10.0; double TauMin=0.0; double TauMax=2.0;
    
    CommandlineArguments.Getval("qTMin",qTMin);
    CommandlineArguments.Getval("qTMax",qTMax);
    CommandlineArguments.Getval("TauMin",TauMin);
    CommandlineArguments.Getval("TauMax",TauMax);
    
    std::cerr << "#KINEMATIC CUTS ARE qT=" << qTMin << " - " << qTMax << " AND  tau=" << TauMin << "-" << TauMax << " fm" << std::endl;
    std::cerr << "##### WARNING NO qT CUT IMPLEMENTED HERE ##########" << std::endl;
    
    // MONTE CARLO SAMPLING //
    int NSamples=512000;
    CommandlineArguments.Getval("NSamples",NSamples);
    
    // SEED RANDOM NUMBER GENERATOR //
    srand48(time(0));
    
    // SETUP INTERPOLATORS FOR HYDRO ATTRACTORS //
    HydroAttractor::Setup();
    
    
//    // CALCULATE dN/d^4xdQ FOR FIXED TIMES //
//    for(int iTau=1;iTau<20;iTau++){
//        
//        // SET TIME //
//        double Tau=0.1*iTau;
//        
//        // GET HYDRO-ATTRACTOR VALUES //
//        double T,wTilde,e,pL,eQOvereG;
//        HydroAttractor::GetValues(dNchdEta,Area,EtaOverS,Tau,T,wTilde,e,pL,eQOvereG);
//        
//        // GET PARAMETERS OF PHASE-SPACE DISTRIBUTION //
//        double Xi,TEff,qSupp;
//        PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,TEff,qSupp);
//        
//        // CALCULATE LO AND NLO CONTRIBUTIONS TO DILEPTON PRODUCTION //
//        std::stringstream ss;
//        ss << "OUTPUT/DileptonRateTAU" << Tau << ".txt";
//        std::string fname=ss.str();
//        
//        DileptonRates::CreatedNdQOutput(QMin,QMax,NQBins,NSamples,Xi,TEff,qSupp,wTilde,fname);
//    }
    
    // CALCULATE dN/dQ INTEGRATED OVER TIME EVOLUTION OF THE SYSTEM //
    int NTau=200;
    
    DileptonRates::CalculatedNdQ(QMin,QMax,NQBins,NSamples,dNchdEta,Area,EtaOverS,TauMin,TauMax,NTau,"OUTPUT/dNdQ.txt");

    
    // EXIT //
    exit(1);
    
    
}
