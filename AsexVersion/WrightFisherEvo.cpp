#include "AddEpi.h"
#include "MersenneTwister.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>

using namespace std;
extern MTRand rnd;


void WrightFisherEvo(int &Loci, int &Individuals, int &Generations, double &RecombinationRate, double &GeomTractLength, double &Va, double &Vi, int &GenomeSize, int &SegSize, int &Sim)
{
    
    int site,ind,i,j,gen,k ; // incrementors
    int RandInd, STOP ;
    double Wmax, ThetaPiStart, ThetaPi, Wmean, Wvar, Wskew ;
    // Recombination related
    int NumRecs, Start, Start2, TractLen, End, RecTotal, Free ;
    // Frequently used constants
    double L = double(Loci) ;
    double AddFitConst = sqrt(Va/L) ;
    // selection coefficient sqrt(Va/L)
    //with additive fitness these need to be small to avoid negative absolute fitnesses
    int Positions[Loci] ;
    vector< vector<int> > RecombinantPositions ;
    vector<int> TempDonorSNPs ;
    
    vector <double> PrintingPoints ;
    double PP = 0.1 ;
    for(i=0; i<9; i++){
        PrintingPoints.push_back(PP) ;
        PP = PP+0.1 ;
    }
    int CurrentPrintingPoint = (PrintingPoints.size() - 1) ;
   
    
    /*============
     DECLARE DYNAMIC VARIABLES
     MAKE SURE TO DELETE! :)
     ==============*/
    Hap * Haplotypes = new Hap [Individuals] ; // array of Hap structs: vector<Locus> Loci ; unsigned int NumHap ;
    Hap ** Pop = new Hap *[Individuals] ; // array of pointers to Hap structs
    Hap ** NextGen = new Hap *[Individuals] ;  // array of pointers to Hap structs
    double * WAbs = new double [Individuals]; // "WAbs" will hold the Absolute fitness of each individual:
    double * EpistasisTable = new double [ ((Loci*(Loci-1))/2) ] ;  // array of pairwise epistatic selection coefficients
    vector<double> * DprimeByDist = new vector<double> [ (Loci-1) ] ; // array of vector<double>
    vector<double> * DqleByDist = new vector<double> [ (Loci-1) ] ;  // array of vector<double>
    vector<double> * DeltarByDist = new vector<double> [ (Loci-1) ] ;  // array of vector<double>
    vector<TwoLocusAFS> * PairwiseAFS = new vector<TwoLocusAFS> [ (Loci-1) ] ; // array of vectors of TwoLocusAFS: int one ; int two ;
    vector<XYZUfreqs> * TwoLocusHaploFreqsBS = new vector<XYZUfreqs> [ (Loci-1) ] ; // BS = before selection
    vector<XYZUfreqs> * TwoLocusHaploFreqsAS = new vector<XYZUfreqs> [ (Loci-1) ] ; // AS = after selection


    string TempSimFileName = "AsexResults" + boost::lexical_cast<string>(Sim);
    ofstream out ;
    out.open(TempSimFileName.c_str()) ;

	out.precision(20) ;
    
    // START CLOCK
    time_t begin, end;
    begin = time(0);
    
    
    /*============
     INITIALIZATION
     ==============*/
    STOP = 0 ;
	RecombinationRate = 0 ; // ASEXUAL SIMULATIONS
    InitializeHaplotypes(Individuals, Loci, Haplotypes) ; //NumHap initialized later
    InitializeEpistasisTable(Loci, EpistasisTable, Sim) ;
    InitializePositions(Positions, Loci, Sim) ;

    
    out << "POSITIONS\n" ;
    for(i=0; i<(Loci); i++){
        out << Positions[i] << " " ;
    }
    out << "\n" ;
    
    out << "AddFitConst: " << AddFitConst << "\n" ;

    out << "EpistasisTable_Ns\n" ;
    // this print epistasis values according to distance, with the top two being adjacent loci
    for(i=0; i<(Loci-1); i++){
        out << Individuals*EpistasisTable[i] << " " ;
        k=0 ;
        for(j=(Loci-1); j>(i+1); j--){
            k+=j ;
            out << Individuals*EpistasisTable[i+k] << " " ;
        }
        out << "\n" ;
    }
    out << "\n" ;
    
    // Assign Individuals to Haplotypes, initiate number in population
    for(ind=0; ind<Individuals; ind++){
        Pop[ind] = &Haplotypes[ind] ;
        (*(Pop[ind])).NumHap = 1 ;
    }

     
    //CALCULATE STARTING METRICS
    CalculateLDstart(Pop, Individuals, Loci, DprimeByDist) ;
    out << "Dprime_BY_DIST_gen0\n" ;
    for(i=0; i<(L-1); i++){
        for(j=0; j<DprimeByDist[i].size() ;j++){
            out << DprimeByDist[i][j] << " " ;
        }
        out << "\n" ;
    }
    out << "\n" ;
    
    ThetaPiStart = CalculateThetaPi(Pop, Individuals, Loci) ;
    out << "ThetaPi_gen0: " << ThetaPiStart << "\n" ;
    
    for(ind=0; ind<Individuals; ind++){
        // 1+MarginalFitness as in Piganeau et al 2001
        WAbs[ind] = 1.0 + MarginalFitness(Pop[ind], Loci, AddFitConst, EpistasisTable) ;
    }
    Wmean = MeanFitness(Individuals, WAbs) ;
    Wvar = VarFitness(Individuals, WAbs, Wmean) ;
    Wskew = SkewFitness(Individuals, WAbs, Wmean, Wvar) ;
    out << "Wmean_gen0: " << Wmean << "\n" ;
    out << "Wvar_gen0: " << Wvar << "\n" ;
    out << "Wskew_gen0: " << Wskew << "\n" ;
    out << "WabsDist_gen0:" << "\n" ;
    for(ind=0; ind<Individuals; ind++){
        out << WAbs[ind] << " " ;
    }
    out << "\n" ;
    /*===========================
     BEGIN WRIGHT-FISHER EVOLUTION
     ===========================*/
    
    for (gen = 1; gen <= Generations; gen++){
        
        Wmax = 0.0 ;
        // CALCULATE FITNESS HERE
        for(ind=0; ind<Individuals; ind++){
            // 1+MarginalFitness as in Piganeau et al 2001
            WAbs[ind] = 1.0 + MarginalFitness(Pop[ind], Loci, AddFitConst, EpistasisTable) ;
            if (Wmax < WAbs[ind]){
                Wmax = WAbs[ind] ;
            }
        }
        
        // CREATE NEXT GENERATION
        for(ind=0; ind<Individuals; ind++){
            do{
                RandInd = int(rnd.randExc(Individuals)) ; // real number in [0,Nv), int -> includes 0 but not Nv
            } while (rnd.rand() > WAbs[RandInd]/Wmax);// Randomly sample individual based on it's relative
            // Need to sample from Pop to simulate drift, not the original population Haplotypes
            NextGen[ind] = Pop[RandInd] ;
        }
        
        // Can this be done in previous loop?
        for(ind=0; ind<Individuals; ind++){
            Haplotypes[ind].NumHap = 0 ;
        }
        
        //Move population to next generation, pointer reassignment
        for(ind=0; ind<Individuals; ind++){
            Pop[ind] = NextGen[ind] ;
            (*(Pop[ind])).NumHap++ ;
        }
        ThetaPi = CalculateThetaPi(Pop, Individuals, Loci) ;
        // Calculate stats right after selection
        if(STOP == 1){
            //Calculate statistics every SampGen
            //LD, Wbar, Wmax, diversity

            CalcTwoLocusHaplotypeFreqs(Pop, Individuals, Loci, TwoLocusHaploFreqsAS) ;
            CalcChangeHaploFreqsDeltar(Pop, Individuals, Loci, TwoLocusHaploFreqsBS, TwoLocusHaploFreqsAS, DeltarByDist) ;

            Wmean = MeanFitness(Individuals, WAbs) ;
            Wvar = VarFitness(Individuals, WAbs, Wmean) ;
            Wskew = SkewFitness(Individuals, WAbs, Wmean, Wvar) ;
            ThetaPi = CalculateThetaPi(Pop, Individuals, Loci) ;
            out << "ThetaPi_gen" << gen << "\t" << "ThetaDecay_" << PrintingPoints[CurrentPrintingPoint+1] << "\t" << ThetaPi << "\n" ;
            out << "Wmean_gen" << gen << "\t" << "ThetaDecay_" << PrintingPoints[CurrentPrintingPoint+1] << "\t" << Wmean << "\n" ;
            out << "Wmax_gen" << gen << "\t" << "ThetaDecay_" << PrintingPoints[CurrentPrintingPoint+1] << "\t" << Wmax << "\n" ;
            out << "Wvar_gen"<< gen << "\t" << "ThetaDecay_" << PrintingPoints[CurrentPrintingPoint+1] << "\t" << Wvar << "\n" ;
            out << "Wskew_gen"<< gen << "\t" << "ThetaDecay_" << PrintingPoints[CurrentPrintingPoint+1] << "\t" << Wskew << "\n" ;
            out << "WabsDist_gen"<< gen << "\t" << "ThetaDecay_" << PrintingPoints[CurrentPrintingPoint+1] << "\n" ;
            for(ind=0; ind<Individuals; ind++){
                out << WAbs[ind] << " " ;
            }
            out << "\n" ;
            CalculatePairwiseAFS(Pop, Individuals, Loci, PairwiseAFS) ;

            CalculateLDmetrics(Pop, Individuals, Loci, TwoLocusHaploFreqsBS, DprimeByDist, DqleByDist) ;
            out << "Dprime_BY_DIST_gen"<< gen << "\t" << "ThetaDecay_" << PrintingPoints[CurrentPrintingPoint+1] << "\n" ;
            for(i=0; i<(L-1); i++){
                for(j=0; j<DprimeByDist[i].size() ;j++){
                    out << DprimeByDist[i][j] << " " ;
                }
                out << "\n" ;
            }
            out << "\n" ;

            if(CurrentPrintingPoint<0){
                break ;
            }else{
                STOP = 0 ;
            }
            
        }
        if(ThetaPi <= (PrintingPoints[CurrentPrintingPoint]*ThetaPiStart)){
            CalcTwoLocusHaplotypeFreqs(Pop, Individuals, Loci, TwoLocusHaploFreqsBS) ;
            STOP = 1 ;
            CurrentPrintingPoint-- ;

        }
        
        
        // RECOMBINE HAPLOTYPES
        for(ind=0; ind<Individuals; ind++){
            NumRecs = int(poisdev(RecombinationRate)) ; // this ind is a recipient of DNA
            if(NumRecs>0){
                RecombinantPositions.clear() ;
                // check how many are effective recombinations
                for(i=0; i<NumRecs; i++){
                    TempDonorSNPs.clear() ;
                    Start = int(rnd.randExc(GenomeSize)) ;
                    do{
                        TractLen = geometric(GeomTractLength) ;
                    }while( TractLen >= GenomeSize-1 ) ;
                    End = (Start+TractLen) ;
                    if( End <= (GenomeSize-1)){
                        j=0 ;
                        for(j=0; j<Loci; j++){
                            // recombination happens to right of Start, so only consider Positions to right
                            if(Start <= Positions[j]){
                                while( (Positions[j] <= End) && (j<Loci) ){
                                    // push j into vector, j++
                                    TempDonorSNPs.push_back( j ) ;
                                    j++ ;
                                }
                                break ;
                            }
                        }
                    }else{
                        j=0 ;
                        for(j=0; j<Loci; j++){
                            // recombination happens to right of Start, so only consider Positions to right
                            if(Start <= Positions[j]){
                                while( (Positions[j] <= (GenomeSize-1)) && (j<Loci) ){ // Positions only include up to GenomeSize-1
                                    // push j into vector, j++
                                    TempDonorSNPs.push_back( j ) ; // record j, not Positions, b/c 0-(L-1) loci
                                    j++ ;
                                }
                                break ;
                            }
                        }
                        k=0 ;
                        while( (Positions[k] <= (End-(GenomeSize-1))) && (k<Loci)){
                            TempDonorSNPs.push_back( k ) ;
                            k++ ;
                        }
                    
                    }
                    if( !TempDonorSNPs.empty() ){
                        RecombinantPositions.push_back( TempDonorSNPs ) ;
                    }
                }
                
                if( !RecombinantPositions.empty() ){
                    (*(Pop[ind])).NumHap-- ; // recipient haplotype taken out of original haplotype category
                    // Find unoccupied haplotype
                    // replace it with this recomb recipient, set NumHap = 1
                    Free=0 ;
                    while( (Haplotypes[Free].NumHap > 0) && (Free < Individuals)){
                        Free++ ;
                    }
                    
                    Haplotypes[Free].Loci = (*(Pop[ind])).Loci ; // Copy recipient haplotype to recombine
                    Pop[ind] = &Haplotypes[Free] ; // reassign Pop pointer to recipient haplotype
                    Haplotypes[Free].NumHap = 1 ; // was previously 0
                    
                    //As it stands, multiple recomb events can have same start position and interfere with eachother
                    for(i=0; i<RecombinantPositions.size(); i++){
                        // RANDOMLY SELECT A DONOR
                        do{
                            RandInd = int(rnd.randExc(Individuals)) ; // real number in [0,Nv), int -> includes 0 but not Nv
                        } while (RandInd == ind) ; // Cant select same individual!!
                        
                        HomoRecomb( Pop[ind], Pop[RandInd], RecombinantPositions[i]) ;

                    }
                }
            }
        }
    }
    /*=========================
     END WRIGHT-FISHER EVOLUTION
     =========================*/
    
    end = time(0);
    int TimeTaken = int(difftime(end,begin)) ;
    out << "hours:min:sec " << TimeTaken/3600 << ":" << (TimeTaken%3600)/60 << ":" << TimeTaken%60 << "\n" ;
    
    
    out.close() ;

    delete [] Haplotypes ;
    delete [] Pop ;
    delete [] NextGen ;
    delete [] WAbs ;
    delete [] EpistasisTable ;
    delete [] DprimeByDist ;
    delete [] DqleByDist ;
    delete [] DeltarByDist ;
    delete [] TwoLocusHaploFreqsBS ;
    delete [] PairwiseAFS ;
    
}

