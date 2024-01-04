#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

enum ModelType
{
    AD=0, // adult defence - adults can only defend their own breeding position
    FS=1 // adults and juveniles compete anew every time step over breeding positions
};

class Parameters
{
    public:

        // juvenile dispersal 
        double d{0.5};

        // number of adult breeders in a patch
        unsigned n{3};

        // adult survival rate
        double m_baseline{0.2};

        std::string file_name{"iter_intergenerational"};

        // maximum amount of time before we give up on 
        // converging relatedness coefficients
        unsigned long max_eco_time{10000};

        // maximum amount of time before we give up on 
        // converging evolving variables
        unsigned long max_evo_time{90000};

        // initial value for survival help
        // before we let the trait evolve
        double survival_help_init{0.0};
        // initial value for fecundity help
        // before we let the trait evolve
        double fecundity_help_init{0.0};
};

#endif 
