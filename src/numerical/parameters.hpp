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

        // competition parameter 
        double k{1.0};

        // number of adult breeders in a patch
        int n{3};

        std::string file_name{"iter_intergenerational"};

        // maximum amount of time before we give up on 
        // converging relatedness coefficients
        unsigned long max_eco_time{10000};

        // maximum amount of time before we give up on 
        // converging evolving variables
        unsigned long max_evo_time{90000};

        // threshold difference below we assume a trait
        // or ecological variable is at equilibrium
        double convergence_threshold{1e-07};

        // initial value for survival help
        // before we let the trait evolve
        double x_init{0.0};
        // initial value for fecundity help
        // before we let the trait evolve
        double y_init{0.0};


        // cost and benefit parameters 
        //
        // benefit coefficient of survival help
        double B_surv{0.1};
        
        // benefit coefficient of fecundity help
        double B_fec{0.1};
        
        // cost coefficient of survival help in terms of survival
        double C_x_surv{0.1};

        // cost coefficient of fecundity help in terms of survival
        double C_y_surv{0.1};
        
        // cost coefficient of survival help in terms of fecundity
        double C_x_fec{0.1};

        // cost coefficient of fecundity help in terms of fecundity
        double C_y_fec{0.1};

        // max ceiling to survival (if no one dies no evolution)
        double surv_max{0.95};

        // the baseline number of offspring that are produced
        // before costs and benefits
        double baseline_fecundity{1.0};

        // whether the type of the model is adult defence
        // or full scramble
        ModelType type{FS};

        // the amount of variation in a trait
        // (more variation faster evolution)
        double var_x{0.05};
        double var_y{0};

        unsigned skip_output{10};

        double vanish_relatedness{1e-05};

};

#endif 
