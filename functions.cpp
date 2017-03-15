#include "AddEpi.h"
#include "MersenneTwister.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>

using namespace std;
extern MTRand rnd;

void InitializeHaplotypes(int N, int L, Hap * H)
{
    int ForInitializingAFS[N] ;
    int MaxPos,RandPos,temp ;
    int ind, site ; // incrementors
    // initialize ForInitializingAFS array
    for(ind=0; ind < N; ind++){
        ForInitializingAFS[ind] = ind ;
    }
    for(site=0; site<L; site++){
        // To initialize Haplotypes, sample without replacement
        // Durstenfeld-Fisher-Yates shuffle ForInitializingAFS array
        // use shuffled indices to assign first half Individuals 1, second half -1 ;
        // shuffle the array, use it, recycle
        MaxPos=(N-1) ;
        while(MaxPos>0){
            temp = ForInitializingAFS[MaxPos] ;
            RandPos=rnd.rand()*MaxPos ; // [0,MaxPos-1], includes boundaries
            ForInitializingAFS[MaxPos] = ForInitializingAFS[RandPos] ;
            ForInitializingAFS[RandPos] = temp ;
            MaxPos-- ;
        }
        for(ind=0; ind<N; ind++){
            Locus TempLoc ;
            H[ ForInitializingAFS[ind] ].Loci.push_back( TempLoc ) ;
            if( ind<(N/2) ){
                // ForInitializingAFS[ind] is randomly selected individual!
                H[ ForInitializingAFS[ind] ].Loci[site].Eff = 1 ;
            }else{
                H[ ForInitializingAFS[ind] ].Loci[site].Eff = -1 ;
            }
        }
    }
}


void InitializeEpistasisTable(double Vint, int L, double * ET)
{
    int i ;
    double two = 2.0 ;
    double L2 = double(L) ; // to calc rescaling factor, don't wanna do math on doubles and int's
    double Rescale = sqrt((two*Vint)/(L2*(L2-1))) ;
    for(i=0; i<(L*(L-1)/2); i++){
        ET[i] = (gasdev()*Rescale) ; // multiply standard normal by desired standard dev to rescale (sigma*x + mu)
    }
}

void InitializePositions(int * Pos, int L, int SegSize)
{
    int i,j ;
    vector<int> PreviousPos ;
    int temp,flag ;
    
    flag = 0 ;
    for(i=0; i<L; i++){
        if( PreviousPos.empty() ){
            temp = int(rnd.randExc(SegSize)) ;
            Pos[i] = temp  ; // real number in [0,SegSize), int -> includes 0 but not SegSize
            PreviousPos.push_back( temp ) ;
        }else{
            do{
                temp = int(rnd.randExc(SegSize))  ;
                flag = 0 ;
                for(j=0; j<PreviousPos.size(); j++){
                    if( temp == PreviousPos[j] ){
                        flag = 1 ;
                    }
                }
            }while(flag) ;
            Pos[i] = temp ;
            PreviousPos.push_back( temp ) ;
        }
    }
    //sort array
    sort(Pos, Pos + L) ;
}


double MeanFitness(int &I, double * Warray)
{
    double mean, Idbl ;
    int ind ;
    mean = 0.0 ;
    Idbl = double(I) ;
    for(ind=0; ind<I; ind++){
        mean += Warray[ind] ;
    }
    return(mean/Idbl) ;
}


double VarFitness(int &I, double * Warray, double &Wm)
{
    double var, Idbl ;
    int ind ;
    var = 0.0 ;
    Idbl = double(I) ;
    for(ind=0; ind<I; ind++){
        var += pow((Warray[ind]-Wm),2) ;
    }
    return(var/Idbl) ;
}

double SkewFitness(int &I, double * Warray, double &Wm, double &Wv)
{
    // Pearson's moment coefficient of skewness, 3rd standardized moment
    double skew, Idbl ;
    int ind ;
    skew = 0.0 ;
    Idbl = double(I) ;
    for(ind=0; ind<I; ind++){
        skew += pow((Warray[ind]-Wm)/(sqrt(Wv)),3) ;
    }
    return(skew/Idbl) ;
}

double MarginalFitness(Hap * p, int &L, double &VaFit, double * ET)
{
    int site ;
    double w = 0.0 ;
    int i,j,k ;
    
    //Compute Additive component of fitness
    for(site=0; site < L; site++){
        w += (*(p)).Loci[site].Eff  ;
    }

    w = w*VaFit ;

    //Compute epistatic component of fitness
    k=0 ;
    for(i=0; i<(L-1); i++){
        for(j=(i+1); j<L; j++){
            w += ET[k] * (*(p)).Loci[i].Eff * (*(p)).Loci[j].Eff ;
            k++ ;
        }
    }
    return(w) ;
}



void HomoRecomb(Hap * Recipient, Hap * Donor, vector<int> &RecPos)
{
    int j ;
    for(j=0; j<RecPos.size(); j++){
        (*(Recipient)).Loci[RecPos[j]].Eff = (*(Donor)).Loci[RecPos[j]].Eff ;
    }
}

void SexLikeRecomb(Hap * Recipient, Hap * Donor, int &St, int &St2, int &L)
{
    int j ;
    if(St < St2){
        for( j=St ; j<=St2 ; j++){ // replace Loci for each
            //recipient[ j ] = donor[ j ] ;
            (*(Recipient)).Loci[j].Eff = (*(Donor)).Loci[j].Eff ;
        }
    }else{ // recomb tract runs over boundary
        for(j=St; j<=(L-1) ; j++){
            (*(Recipient)).Loci[j].Eff = (*(Donor)).Loci[j].Eff ;
        }
        for(j=0; j<=St2 ; j++){ // if 0, then it went over by 1
            (*(Recipient)).Loci[j].Eff = (*(Donor)).Loci[j].Eff ;
        }
    }
}


void CalculatePairwiseAFS(Hap * p[], int &I, int &L, vector<TwoLocusAFS> * PWAFS)
{
    int i,j,ind,m ;
    int a,b ;
    double aAF, bAF;
    TwoLocusAFS TempAFS ;
    for(i=0; i<(L-1); i++){
        PWAFS[i].clear() ;
    }
    
    for(i=0; i<(L-1); i++){
        m=0 ;
        for(j=(i+1); j<L; j++){
            a=0, b=0 ;
            for(ind=0; ind<I; ind++){
                if( (*(p[ind])).Loci[i].Eff == 1){
                    a++ ;
                }
                if( (*(p[ind])).Loci[j].Eff == 1){
                    b++ ;
                }
            }
            TempAFS.one = a ;
            TempAFS.two = b ;
                   
            PWAFS[m].push_back(TempAFS) ;
            m++ ;
        }
    }
}

void CalcTwoLocusHaplotypeFreqs(Hap * p[], int &I, int &L, vector<XYZUfreqs> * TLHF)
{
    int i,j,ind,m ;
    int xBS,yBS,zBS,uBS ;
    XYZUfreqs TempFreqs ;
    
    for(i=0; i<(L-1); i++){
        TLHF[i].clear() ;
    }
    
    for(i=0; i<(L-1); i++){
        m=0 ;
        for(j=(i+1); j<L; j++){
            xBS=0, yBS=0, zBS=0, uBS=0 ;
            for(ind=0; ind<I; ind++){
                
                if( ((*(p[ind])).Loci[i].Eff == 1) && ((*(p[ind])).Loci[j].Eff == 1) ){
                    xBS++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == -1) && ((*(p[ind])).Loci[j].Eff == 1) ){
                    yBS++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == 1) && ((*(p[ind])).Loci[j].Eff == -1) ){
                    zBS++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == -1) && ((*(p[ind])).Loci[j].Eff == -1) ){
                    uBS++ ;
                }
            }
            TempFreqs.x = xBS ;
            TempFreqs.y = yBS ;
            TempFreqs.z = zBS ;
            TempFreqs.u = uBS ;

            TLHF[m].push_back( TempFreqs ) ;
            m++ ;
        }
    }
}

void CalcChangeHaploFreqsDeltar(Hap * p[], int &I, int &L, vector<XYZUfreqs> * TLHFBS, vector<XYZUfreqs> * TLHFAS,  vector<double> * Deltar)
{
    int i,j,m ;
    double rBS,rAS ;
    double rBSnum, rBSdenom, rASnum, rASdenom ;
    double logrBS, logrAS ;
    
    for(i=0; i<(L-1); i++){
        m=0 ;
        for(j=(i+1); j<L; j++){
            
            rBSnum = double( (TLHFBS[i][m].x)*(TLHFBS[i][m].u) ) ;
            rBSdenom = double( (TLHFBS[i][m].y)*(TLHFBS[i][m].z) ) ;
            rBS = rBSnum/rBSdenom ;
            rASnum = double( (TLHFAS[i][m].x)*(TLHFAS[i][m].u) ) ;
            rASdenom = double( (TLHFAS[i][m].y)*(TLHFAS[i][m].z) ) ;
            rAS = rASnum/rASdenom ;
            logrBS = log(rBS) ;
            logrAS = log(rAS) ;
            Deltar[m].push_back( logrAS - logrBS ) ;

            m++ ;
        }
    }
}

void CalculateLDstart(Hap * p[], int &I, int &L, vector<double> * DprimeBD)
{
    int i,j,ind,m ; // incrememtors i/j go through sites, i/m go through pairwise array or vectors
    int a,b ; // site allele frequencies
    int X ; // haplotype allele frequencies, x=11, y=-11, z=1-1, u=-1-1
    double Pa, Pb, PX ;
    double Dab, Dmax ;
    double Idbl = double(I) ;
    //calculate Dab = Pab - (Pa*Pb)
    //if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
    //if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
    
    for(i=0; i<(L-1); i++){
        DprimeBD[i].clear() ;
    }
    
    for(i=0; i<(L-1); i++){
        m=0 ;
        for(j=(i+1); j<L; j++){
            a=0, b=0, X=0 ;
            for(ind=0; ind<I; ind++){
                if( ((*(p[ind])).Loci[i].Eff == 1) && ((*(p[ind])).Loci[j].Eff == 1) ){
                    X++ ;
                    a++ ;
                    b++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == -1) && ((*(p[ind])).Loci[j].Eff == 1) ){
                    b++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == 1) && ((*(p[ind])).Loci[j].Eff == -1) ){
                    a++ ;
                }
            }
            // ALLELE/HAPLOTYPE FREQUENCY FROM CURRENT GENERATION
            Pa = a/Idbl ;
            Pb = b/Idbl ;
            PX = X/Idbl ;
            //Calculate D' ;
            Dab = PX - (Pa*Pb);
            if( Dab >= 0 ){
                if( Pa*(1-Pb) < Pb*(1-Pa) ){
                    Dmax = Pa*(1-Pb) ;
                }else{
                    Dmax = Pb*(1-Pa) ;
                }
            }else if (Dab < 0){
                if( Pa*Pb < (1-Pa)*(1-Pb) ){
                    Dmax = Pa*Pb ;
                }else{
                    Dmax = (1-Pa)*(1-Pb) ;
                }
            }
            DprimeBD[m].push_back( Dab/Dmax ) ;
            m++ ;
        }
    }
    /*
     NOTES:
     Although each pair of alleles has its own D, the values of different pairs are constrained by the fact that the allele freq's at both loci must sum to 1. If both loci are diallelic, the constraint is strong enough that only one value of D is needed to characterize LD between those loci, i.e. D_AB == D_ab, and the sign of D is arbitrary and depends on which allele you start with.
     */
    
}







void CalculateLDmetrics(Hap * p[], int &I, int &L, vector<XYZUfreqs> * TLHF, vector<double> * DprimeBD, vector<double> * DQLEBD)
{
    int i,j,ind,m ; // incrememtors i/j go through sites, i/m go through pairwise array or vectors
    int a,b ; // site allele frequencies
    int X,Y,Z,U ; // haplotype allele frequencies, x=11, y=-11, z=1-1, u=-1-1
    double Pa, Pb, PX, PY, PZ, PU ;
    double PxBS, PyBS, PzBS, PuBS ;
    double Dab, Dmax, Dqle ;
    double Idbl = double(I) ;
    //calculate Dab = Pab - (Pa*Pb)
    //if Dab positive, Dmax = min( Pa*(1-Pb), Pb*(1-Pa) )
    //if Dab negative, Dmax = min( Pa*(Pb), (1-Pa)*(1-Pb) )
    
    for(i=0; i<(L-1); i++){
        DprimeBD[i].clear() ;
        DQLEBD[i].clear() ;
    }
    
    for(i=0; i<(L-1); i++){
        m=0 ;
        for(j=(i+1); j<L; j++){
            a=0, b=0, X=0, Y=0, Z=0, U=0 ;
            for(ind=0; ind<I; ind++){
                if( ((*(p[ind])).Loci[i].Eff == 1) && ((*(p[ind])).Loci[j].Eff == 1) ){
                    X++ ;
                    a++ ;
                    b++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == -1) && ((*(p[ind])).Loci[j].Eff == 1) ){
                    Y++ ;
                    b++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == 1) && ((*(p[ind])).Loci[j].Eff == -1) ){
                    Z++ ;
                    a++ ;
                }
                if( ((*(p[ind])).Loci[i].Eff == -1) && ((*(p[ind])).Loci[j].Eff == -1) ){
                    U++ ;
                }
            }
            // ALLELE/HAPLOTYPE FREQUENCY FROM CURRENT GENERATION
            Pa = a/Idbl ;
            Pb = b/Idbl ;
            PX = X/Idbl ;
            PY = Y/Idbl ;
            PZ = Z/Idbl ;
            PU = U/Idbl ;
            //Calculate D' ;
            Dab = PX - (Pa*Pb);
            if( Dab >= 0 ){
                if( Pa*(1-Pb) < Pb*(1-Pa) ){
                    Dmax = Pa*(1-Pb) ;
                }else{
                    Dmax = Pb*(1-Pa) ;
                }
            }else if (Dab < 0){
                if( Pa*Pb < (1-Pa)*(1-Pb) ){
                    Dmax = Pa*Pb ;
                }else{
                    Dmax = (1-Pa)*(1-Pb) ;
                }
            }
            DprimeBD[m].push_back( Dab/Dmax ) ;
            
            // See definition of D in Kimura 1965: weigh by fitnesses?
            // D= XU - YZ, calculate from current gen, after selection
            // little x, y, z, u take from previous generation
            
            PxBS = ( TLHF[i][m].x/Idbl ) ;
            PyBS = ( TLHF[i][m].y/Idbl ) ;
            PzBS = ( TLHF[i][m].z/Idbl ) ;
            PuBS = ( TLHF[i][m].u/Idbl ) ;

            // THIS NEEDS TO BE UPDATED! USE HAPLO FREQS FROM PREV GENERATION
            
            Dqle = (0.5)*((PX*PU) - (PY*PZ))*( (1.0/PxBS)+(1.0/PyBS)+(1.0/PzBS)+(1.0/PuBS) ) ;
            DQLEBD[m].push_back( Dqle ) ;
            m++ ;
            

        }
    }
    /*
     NOTES:
     Although each pair of alleles has its own D, the values of different pairs are constrained by the fact that the allele freq's at both loci must sum to 1. If both loci are diallelic, the constraint is strong enough that only one value of D is needed to characterize LD between those loci, i.e. D_AB == D_ab, and the sign of D is arbitrary and depends on which allele you start with.
     */
    
}

double CalculateThetaPi(Hap * p[], int &I, int &L)
{
    int site, ind ;
    double AF ;
    double pi = 0.0 ;
    double Idbl = double(I) ;
    double bc = BinomCoeff(Idbl,2.0) ;
    for(site=0; site<L; site++){
        AF = 0.0 ;
        for(ind=0; ind<I; ind++){
            if ( (*(p[ind])).Loci[site].Eff == 1 ){
                AF+=1.0 ;
            }
        }
        pi += (AF)*(Idbl-AF)/bc ;
    }
    return(pi) ;
}

double BinomCoeff(double n, double k)
{
    if(n>k)
    {
        double r = 1 ;
        do
        {
            r*=(n/(n-k)) ;
            n-- ;
        }while(n>k) ;
        return r ;
    }
    else
    {
        return 0 ;
    }
}


ifstream& GotoLine(std::ifstream& file, unsigned int num)
{
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


void ReadParameterFile(int &Loc, int &Ind, int &Gen, double &RecRate, double &GeoTrLen, double &Vadd, double &Vint, int &GenSize, int &SegSize, int &NmSm)
{
    ifstream fin("Parameters") ;
    GotoLine(fin, 2);

    fin >> Loc  ;
    fin >> Ind  ;
    fin >> Gen  ;
    fin >> RecRate  ;
    fin >> GeoTrLen  ;
    fin >> Vadd  ;
    fin >> Vint  ;
    fin >> GenSize ;
    fin >> SegSize ;
    fin >> NmSm  ;

    fin.close() ;
}






