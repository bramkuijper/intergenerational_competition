#ifndef _INTERGENERATIONAL_HPP_
#define _INTERGENERATIONAL_HPP_

#include "parameters.hpp"

class InterGenerational
{
    public:
        // the class constructor that sets everything up
        InterGenerational(Parameters const &params);

        // the function that actually runs the thing
        void run();

    private:

        // the actually evolving traits: survival and fecundity help
        double x{0.0}; // survival help
        double y{0.0}; // fecundity help

        // equilibrium relatedness between two distinct breeders
        double rd{0.0};

        // object containin all parameters
        Parameters params{};

        // function to iteratively solve relatedness coefficients
        void solve_relatedness();

        // keep track of the current time step in ecological...
        unsigned long time_step_eco{0};
        // ... and evolutionary time
        unsigned long time_step_evo{0};
};

#endif
