namespace DileptonRates{
    
    
    // SCREENING PARAMETER //
    static const double xi0g=std::exp(5.0/6.0)/(2.0*M_SQRT2);
    static const double xi0q=std::exp(1.0)/M_SQRT2;
    
    
    double GetAnalyticLORate(double Q,double TEff){
        
        double QSqr=Q*Q;
        
        double Nq=8192;
        
        double qMin=0.0; double qMax=12.0*TEff; double dq=(qMax-qMin)/double(Nq);
        
        double Value=0.0;
        
        for(int iq=0;iq<Nq;iq++){
            
            double q=qMin+(iq+0.5)*dq;
            
            double q0=sqrt(QSqr+q*q);
            
            double PreFactor=alphaEM*alphaEM/(12.0*M_PI*M_PI*M_PI*M_PI)*(1.0+mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum*Nc;
            double Jacobian=(4.0*M_PI*q*q)*dq;
            
            Value+=Jacobian*PreFactor*(Q/q0)*(2.0*TEff/q)*log(cosh((q0+q)/(4.0*TEff))/cosh((q0-q)/(4.0*TEff)))/(exp(q0/TEff)-1.0);
            
        }
        
        return Value;
        
        
        
    }
    
    
    void CalculateRate_dNdQ(double QMin,double QMax,int NQBins,int NSamples,double Weight,double Xi,double TEff,double qSupp,double wTilde,double *dNlld4xd4Q_LO,double *dNlld4xd4Q_NLO_A,double *dNlld4xd4Q_NLO_B){
        
        // GET SCREENING MASSES //
        double mDSqr=PhaseSpaceDistribution::mDSqr(Xi,TEff,qSupp);
        double mQSqr=PhaseSpaceDistribution::mQSqr(Xi,TEff,qSupp);
        
        double dQ=(QMax-QMin)/double(NQBins);
        
        // SAMPLE DILEPTON PRODUCTION //
        for(int i=0;i<NSamples;i++){
            
            double p1Min=0.0; double p1Max=12.0*TEff;
            double p2Min=0.0; double p2Max=12.0*TEff;
            double p3Min=0.0; double p3Max=12.0*TEff;
            
            // SAMPLE MOMENTA p1,p2,p3 //
            double p1=(p1Max-p1Min)*drand48()+p1Min;
            double p2=(p2Max-p2Min)*drand48()+p2Min;
            double p3=(p3Max-p3Min)*drand48()+p3Min;
            
            double Cos1=2.0*drand48()-1.0;
            double Cos2=2.0*drand48()-1.0;
            double Cos3=2.0*drand48()-1.0;
            
            double Phi1=2.0*M_PI*drand48();
            double Phi2=2.0*M_PI*drand48();
            double Phi3=2.0*M_PI*drand48();
            
            double p1Vec[3]={p1*sqrt(1.0-Cos1*Cos1)*cos(Phi1),p1*sqrt(1.0-Cos1*Cos1)*sin(Phi1),p1*Cos1};
            double p2Vec[3]={p2*sqrt(1.0-Cos2*Cos2)*cos(Phi2),p2*sqrt(1.0-Cos2*Cos2)*sin(Phi2),p2*Cos2};
            double p3Vec[3]={p3*sqrt(1.0-Cos3*Cos3)*cos(Phi3),p3*sqrt(1.0-Cos3*Cos3)*sin(Phi3),p3*Cos3};
            
            double pT1=sqrt(p1Vec[0]*p1Vec[0]+p1Vec[1]*p1Vec[1]);
            double yP1=atanh(p1Vec[2]/p1);
            
            double pT2=sqrt(p2Vec[0]*p2Vec[0]+p2Vec[1]*p2Vec[1]);
            double yP2=atanh(p2Vec[2]/p2);
            
            double pT3=sqrt(p3Vec[0]*p3Vec[0]+p3Vec[1]*p3Vec[1]);
            double yP3=atanh(p3Vec[2]/p3);
            
            // GET MANDELSTAM VARIABLES //
            double s=+2.0*(p1*p2-(p1Vec[0]*p2Vec[0]+p1Vec[1]*p2Vec[1]+p1Vec[2]*p2Vec[2]));
            double u=-2.0*(p1*p3-(p1Vec[0]*p3Vec[0]+p1Vec[1]*p3Vec[1]+p1Vec[2]*p3Vec[2]));
            double t=-2.0*(p2*p3-(p2Vec[0]*p3Vec[0]+p2Vec[1]*p3Vec[1]+p2Vec[2]*p3Vec[2]));
            
            double wu=(p1-p3);
            double wt=(p2-p3);
            
            double qu=std::sqrt((p1Vec[0]-p3Vec[0])*(p1Vec[0]-p3Vec[0])+(p1Vec[1]-p3Vec[1])*(p1Vec[1]-p3Vec[1])+(p1Vec[2]-p3Vec[2])*(p1Vec[2]-p3Vec[2]));
            double qt=std::sqrt((p2Vec[0]-p3Vec[0])*(p2Vec[0]-p3Vec[0])+(p2Vec[1]-p3Vec[1])*(p2Vec[1]-p3Vec[1])+(p2Vec[2]-p3Vec[2])*(p2Vec[2]-p3Vec[2]));
            
            
            // SCREENED VERSION OF u,t //
            double tBar=t-xi0g*xi0g*mDSqr;
            double uBar=u-xi0g*xi0g*mDSqr;
            
            double tHat=t-xi0q*xi0q*mQSqr;
            double uHat=u-xi0q*xi0q*mQSqr;
            
            
            // GET PHASE SPACE DISTRIBUTION //
            double f1Q=PhaseSpaceDistribution::fQ(pT1,yP1,0.0,Xi,TEff,qSupp);
            double f2Q=PhaseSpaceDistribution::fQ(pT2,yP2,0.0,Xi,TEff,qSupp);
            double f3Q=PhaseSpaceDistribution::fQ(pT3,yP3,0.0,Xi,TEff,qSupp);
            
            
            double f1G=PhaseSpaceDistribution::fG(pT1,yP1,0.0,Xi,TEff,qSupp);
            double f2G=PhaseSpaceDistribution::fG(pT2,yP2,0.0,Xi,TEff,qSupp);
            double f3G=PhaseSpaceDistribution::fG(pT3,yP3,0.0,Xi,TEff,qSupp);
            
            
            // LEADING ORDER //
            {
                
                // CALCULATE DILEPTON ENERGY AND MOMENTUM //
                double q0=p1+p2;
                double qVec[3]={p1Vec[0]+p2Vec[0],p1Vec[1]+p2Vec[1],p1Vec[2]+p2Vec[2]};
                double q=sqrt(qVec[0]*qVec[0]+qVec[1]*qVec[1]+qVec[2]*qVec[2]);
                
                // SET JACOBIANS FOR NUMERICAL INTEGRATION MEASUERES \int d^3p/((2pi)^3 2p)  //
                double Jacobian=(4.0*M_PI*p1*p1)/(2.0*p1)*(4.0*M_PI*p2*p2)/(2.0*p2)*(p1Max-p1Min)*(p2Max-p2Min)/std::pow(2.0*M_PI,6)*std::pow(2.0*M_PI,4);
                
                // GET INVARIANT MASS QSqr //
                double QSqr=q0*q0-q*q; double Q=sqrt(QSqr);
                
                // GET PRE-FACTOR //
                double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
                
                // MATRIX ELEMENT //
                double M2qqGammaStar=4.0*Nc*QSqr;
                
                // POLARIZATION TENSOR //
                double TracePi=f1Q*f2Q*M2qqGammaStar;
                
                
                // GET BIN AND UPDATE DILEPTON RATE //
                int iQ=int((Q-QMin)/dQ);
                
                if(iQ>=0 && iQ<=NQBins-1){
                    dNlld4xd4Q_LO[iQ]+=Weight*Jacobian*PreFactor*TracePi;
                }
                
            }
            
            // NEXT-TO-LEADING ORDER //
            {
                
                
                // CALCULATE DILEPTON ENERGY AND MOMENTUM //
                double q0=p1+p2-p3;
                double qVec[3]={p1Vec[0]+p2Vec[0]-p3Vec[0],p1Vec[1]+p2Vec[1]-p3Vec[1],p1Vec[2]+p2Vec[2]-p3Vec[2]};
                double q=sqrt(qVec[0]*qVec[0]+qVec[1]*qVec[1]+qVec[2]*qVec[2]);
                
                double Jacobian=(4.0*M_PI*p1*p1)/(2.0*p1)*(4.0*M_PI*p2*p2)/(2.0*p2)*(4.0*M_PI*p3*p3)/(2.0*p3)*(p1Max-p1Min)*(p2Max-p2Min)*(p3Max-p3Min)/std::pow(2.0*M_PI,9)*std::pow(2.0*M_PI,4);
                
                // CHECK IF PHASE-SPACE IS VALID //
                if(q0>q){
                    
                    // GET INVARIANT MASS QSqr //
                    double QSqr=q0*q0-q*q; double Q=sqrt(QSqr);
                    
                    // GET PRE-FACTOR AND CALCULATE DILEPTON RATE //
                    double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
                    
                    // MATRIX ELEMENT //
                    double M2qgGammaStar=8.0*g*g*CF*Nc*(-(u-2.0*QSqr)/s-(s-2.0*QSqr)/uHat-2.0*QSqr*QSqr/(s*uHat));
                    double M2qqGammaStar=8.0*g*g*CF*Nc*((u-2.0*QSqr)/tHat+(t-2.0*QSqr)/uHat+2.0*QSqr*QSqr/(tHat*uHat));
                    
                    if(M2qgGammaStar<0 || M2qqGammaStar<0){
                        std::cerr << "#WARNING NEGATIVE MATRIX ELEMENT" << M2qgGammaStar << " " << M2qqGammaStar << std::endl;
                    }

                    
                    // POLARIZATION TENSOR //
                    double TracePi_A=f1Q*f2G*(1.0-f3Q)*M2qgGammaStar;
                    double TracePi_B=f1Q*f2Q*(1.0+f3G)*M2qqGammaStar;
                    
                    // GET BIN AND UPDATE DILEPTON RATE //
                    int iQ=int((Q-QMin)/dQ);
                    
                    if(iQ>=0 && iQ<=NQBins-1){
                        dNlld4xd4Q_NLO_A[iQ]+=Weight*Jacobian*PreFactor*TracePi_A;
                        dNlld4xd4Q_NLO_B[iQ]+=Weight*Jacobian*PreFactor*TracePi_B;

                    }
                    
                }
                
            }
            
            
        }
        
        
        
    }
    
    
    // CREATE OUTPUT OF dN/d4xdQ //
    void CreatedNdQOutput(double QMin,double QMax,int NQBins,int NSamples,double Xi,double TEff,double qSupp,double wTilde,std::string fname){
        
        
        // SET UP BINNING IN INVARIANT MASS Q //
        double *dNlld4xd4Q_LO=new double[NQBins];
        double *dNlld4xd4Q_NLO_A=new double[NQBins];
        double *dNlld4xd4Q_NLO_B=new double[NQBins];
        
        for(int iQ=0;iQ<NQBins;iQ++){
            dNlld4xd4Q_LO[iQ]=0.0;
            dNlld4xd4Q_NLO_A[iQ]=0.0;
            dNlld4xd4Q_NLO_B[iQ]=0.0;
        }
        // CALCULATE RATE //
        CalculateRate_dNdQ(QMin,QMax,NQBins,NSamples,1.0,Xi,TEff,qSupp,wTilde,dNlld4xd4Q_LO,dNlld4xd4Q_NLO_A,dNlld4xd4Q_NLO_B);
        
        
        // NORMALIZE TO NUMBER OF SAMPLES //
        for(int iQ=0;iQ<NQBins;iQ++){
            
            dNlld4xd4Q_LO[iQ] /=double(NSamples);
            dNlld4xd4Q_NLO_A[iQ]/=double(NSamples);
            dNlld4xd4Q_NLO_B[iQ]/=double(NSamples);
            
            
        }
        
        
        // CREATE OUTPUT //
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        // GET SCREENING MASSES //
        double mDSqr=PhaseSpaceDistribution::mDSqr(Xi,TEff,qSupp);
        double mQSqr=PhaseSpaceDistribution::mQSqr(Xi,TEff,qSupp);
        
        OutStream << "#wTilde=" << wTilde << " #Xi=" << Xi << " #TEff=" << TEff << " #qSupp=" << qSupp << std::endl;
        OutStream << "#mDSqr=" << mDSqr << " #mQSqr=" << mQSqr << std::endl;
        
        double dQ=(QMax-QMin)/double(NQBins);
        
        for(int iQ=0;iQ<NQBins;iQ++){
            
            double Q=QMin+(iQ+0.5)*dQ;
            OutStream << Q << " " << dNlld4xd4Q_LO[iQ]/dQ << " " << GetAnalyticLORate(Q,TEff) << " "  << dNlld4xd4Q_NLO_A[iQ]/dQ << " "  << dNlld4xd4Q_NLO_B[iQ]/dQ << std::endl;
        }
        
        OutStream.close();
        
    }
    
    
    void CalculatedNdQ(double QMin,double QMax,int NQBins,int NSamples,double dNchdEta,double Area,double EtaOverS,double TauMin,double TauMax,int NTau,std::string fname){
        
        
        // SET UP BINNING IN INVARIANT MASS Q //
        double *dNlld4xd4Q_LO=new double[NQBins];
        double *dNlld4xd4Q_NLO_A=new double[NQBins];
        double *dNlld4xd4Q_NLO_B=new double[NQBins];
        
        for(int iQ=0;iQ<NQBins;iQ++){
            dNlld4xd4Q_LO[iQ]=0.0;
            dNlld4xd4Q_NLO_A[iQ]=0.0;
            dNlld4xd4Q_NLO_B[iQ]=0.0;
        }
        
        // PERFORM INTEGRATION OF THE RATE IN TIME //
        double dTau=(TauMax-TauMin)/double(NTau);
        
        for(int iTau=0;iTau<NTau;iTau++){
            
            // GET TAU //
            double Tau=TauMin+(iTau+0.5)*dTau;
            
            // GET HYDRO-ATTRACTOR VALUES //
            double T,wTilde,e,pL,eQOvereG;
            HydroAttractor::GetValues(dNchdEta,Area,EtaOverS,Tau,T,wTilde,e,pL,eQOvereG);
            
            // GET PARAMETERS OF PHASE-SPACE DISTRIBUTION //
            double Xi,TEff,qSupp;
            PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,TEff,qSupp);
            
            // GET WEIGHT //
            double Weight=Tau*dTau*Area;
            
            // CALCULATE RATE //
            CalculateRate_dNdQ(QMin,QMax,NQBins,NSamples,Weight,Xi,TEff,qSupp,wTilde,dNlld4xd4Q_LO,dNlld4xd4Q_NLO_A,dNlld4xd4Q_NLO_B);
            
            std::cerr << "#PROGESS IS " << 100.0*(iTau+1)/double(NTau) << "% Tau=" << Tau  << std::endl;
            
            
        }
        
        
        // NORMALIZE TO NUMBER OF SAMPLES //
        for(int iQ=0;iQ<NQBins;iQ++){
            
            dNlld4xd4Q_LO[iQ] /=double(NSamples);
            dNlld4xd4Q_NLO_A[iQ]/=double(NSamples);
            dNlld4xd4Q_NLO_B[iQ]/=double(NSamples);

            
        }
        
        
        // CREATE OUTPUT //
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        
        double dQ=(QMax-QMin)/double(NQBins);
        
        for(int iQ=0;iQ<NQBins;iQ++){
            
            double Q=QMin+(iQ+0.5)*dQ;
            OutStream << Q << " " << dNlld4xd4Q_LO[iQ]/dQ << " "  << dNlld4xd4Q_NLO_A[iQ]/dQ << " "  << dNlld4xd4Q_NLO_B[iQ]/dQ << std::endl;
        }
        
        OutStream.close();
        
        
    }
    
    
}