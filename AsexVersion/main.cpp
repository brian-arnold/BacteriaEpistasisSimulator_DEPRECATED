#include "AddEpi.h"
#include "MersenneTwister.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
MTRand rnd;



/*===================================================
TO DO LIST
    -Rescaling studies (decrease N, increase Va/Vi, do you see similar patterns?)
    -sexual recombinaiton? 2 random positions, reciprocal swapping
 ==================================================*/



int main()
{
    // S appendix indicates starting variable to be passed to WF function
    int LociS, IndividualsS, GenerationsS, GenomeSizeS, SegSizeS, NumSim ;
    double UnscaledRecombRateS, RecombinationRateS, GeomTractLengthS, VaS, ViS ;
    
    ReadParameterFile(LociS, IndividualsS, GenerationsS, UnscaledRecombRateS, GeomTractLengthS, VaS, ViS, GenomeSizeS, SegSizeS, NumSim) ;
    //RecombinationRateS = UnscaledRecombRateS*LociS ;
    RecombinationRateS = UnscaledRecombRateS ;
    int sim ;
    for(sim=1; sim<=NumSim; sim++){
        WrightFisherEvo(LociS, IndividualsS, GenerationsS, RecombinationRateS, GeomTractLengthS, VaS, ViS, GenomeSizeS, SegSizeS, sim) ;
    }
    
    
    return 0 ;
}

