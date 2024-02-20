#ifndef _INTERGENERATIONAL_HPP_
#define _INTERGENERATIONAL_HPP_

#include <fstream>
#include <iostream>
#include <sstream>
#include "parameters.hpp"

class InterGenerational
{
    public:
        // the class constructor that sets everything up
        InterGenerational(Parameters const &params);

        // the function that actually runs the thing
        void run();

    private:

        // object containin all parameters
        Parameters params{};

        // the actually evolving traits: survival and fecundity help
        double x{0.0}; // survival help
        double y{0.0}; // fecundity help

        // equilibrium relatedness between two distinct breeders
        double rd_ad{0.0};

        // function to iteratively solve relatedness coefficients
        void solve_relatedness();

        std::ofstream output_file;

        // keep track of the current time step in ecological...
        unsigned long time_step_eco{0};
        // ... and evolutionary time
        unsigned long time_step_evo{0};

        // output functions
        void write_data_headers();
        void write_data();
        void write_parameters();


        // fitness functions and their derivatives
        

        // total number of newborns at a local patch
        double u(
                double const y_loc 
                ,double const y
                ,double const x_loc
                ,double const x) const;
        
        // derivative of u wrt first argument (y)
        double d_u_d_y1() const;


        // derivative of u wrt third argument (x)
        double d_u_d_x1() const;
        
        // focal breeder's fecundity
        double f(
                double const y_foc, 
                double const y_loc, 
                double const x_foc) const;

        // derivative of fecundity wrt first arg
        double d_f_d_y1() const;

        // derivative of fecundity wrt 2nd arg
        double d_f_d_y2() const;

        // derivative of fecundity wrt x (last arg)
        double d_f_d_x() const;

        // focal breeder's survival probability
        double sa(
                double const x_foc, 
                double const x_loc, 
                double const y_foc) const;

        // derivative of the focal breeder's survival probability
        // wrt to the focal's survival help trait
        double d_sa_d_x1() const;

        // derivative of the focal breeder's survival probability
        // wrt to the focal's fecundity help trait
        double d_sa_d_y() const;

        // derivative of the focal breeder's survival probability
        // wrt to the local breeders' survival help trait
        double d_sa_d_x2() const;
        
        // partial derivative dW_FS / d y_foc
        double d_WFS_d_yfoc() const;

        // partial derivative dW_FS / d x_foc
        double d_WFS_d_xfoc() const;

        // partial derivative dW_FS / d y_loc
        double d_WFS_d_yloc() const;

        // partial derivative dW_FS / d_x_loc
        double d_WFS_d_xloc() const;
       
        // AD model

        // partial derivative dW_FS / d y_foc
        double d_WAD_d_yfoc() const;

        // partial derivative dW_FS / d x_foc
        double d_WAD_d_xfoc() const;

        // partial derivative dW_FS / d y_loc
        double d_WAD_d_yloc() const;

        // partial derivative dW_FS / d_x_loc
        double d_WAD_d_xloc() const;

        double q_a();
        double ell();
};
#endif
