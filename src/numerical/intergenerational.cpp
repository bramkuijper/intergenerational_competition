#include "intergenerational.hpp"
#include "parameters.hpp"

// the class constructor that sets everything up
InterGenerational::InterGenerational(Parameters const &params_arg) :
    params{params_arg}
    ,x{params_arg.survival_help_init}
    ,y{params_arg.fecundity_help_init}
{
}

// focal breeder's survival probability
double InterGenerational::s_a(double const surv_help, double const fec_help)
{
} // end s_a

//probability adult successfully competes for breeding position
void InterGenerational::p_keep_FS()
{
    s_a(x, y) / 
}

// iteratively solve relatedness coefficients
void InterGenerational::solve_relatedness()
{
    double rd_tplus1{0.0}, qa_val{0.0};

    for (time_step_eco = 0; time_step_eco <= params.max_eco_time; ++time_step_eco)
    {
        qa_val = q_a();
        rd_tplus1 = q_a();
    }
}

// the 
void InterGenerational::run()
{
    for (time_step_evo = 0; time_step_evo <= params.max_evo_time; ++time_step_evo)
    {
        solve_relatedness();
    }
}

