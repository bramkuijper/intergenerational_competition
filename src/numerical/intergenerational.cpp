#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include "intergenerational.hpp"
#include "parameters.hpp"

// the class constructor that sets everything up
InterGenerational::InterGenerational(Parameters const &params_arg) :
    params{params_arg}
    ,x{params_arg.x_init}
    ,y{params_arg.y_init}
    ,output_file{params_arg.file_name}
{
}

void InterGenerational::write_data_headers()
{
    output_file << "t_eco;t_evo;x;y;r_ad" << std::endl;
}

void InterGenerational::write_data()
{
    output_file << time_step_eco << ";" 
        << time_step_evo << ";" 
        << x << ";" 
        << y << ";" 
        << rd_ad << std::endl;
}

void InterGenerational::write_parameters()
{
    output_file << std::endl << std::endl 
        << "type;" << (params.type == AD ? "AD" : "FS") << std::endl
        << "n;" << params.n << std::endl
        << "k;" << params.k << std::endl
        << "d;" << params.d << std::endl
        << "B_surv;" << params.B_surv << std::endl
        << "B_fec;" << params.B_fec << std::endl
        << "C_x_surv;" << params.C_x_surv << std::endl
        << "C_x_fec;" << params.C_x_fec << std::endl
        << "C_y_surv;" << params.C_y_surv << std::endl
        << "C_y_fec;" << params.C_y_fec << std::endl
        << "surv_max;" << params.surv_max << std::endl
        << "baseline_fecundity;" << params.baseline_fecundity << std::endl
        << "var;" << params.var << std::endl
        << "x_init;" << params.x_init << std::endl
        << "y_init;" << params.y_init << std::endl;
}


// focal breeder's fecundity
double InterGenerational::f(
        double const y_foc, 
        double const y_loc, 
        double const x_foc) const
{
    double val{
        params.baseline_fecundity - params.C_y_fec * y_foc 
                                + params.B_fec * y_loc 
                                - params.C_x_fec * x_foc};

    // set lower bound to fecundity of 0
    val = val < 0.0 ? 0.0 : val;

    return(val);
} // end f()

// derivative of fecundity wrt to fecundity help, 1st argument
// of the function f(y,y,x,x)
double InterGenerational::d_f_d_y1() const
{
    return(-params.C_y_fec);
} // end df_d_yfoc

// derivative of fecundity wrt to fecundity help, 2nd argument
// of the function f(y,y,x,x)
double InterGenerational::d_f_d_y2() const
{
    return(params.B_fec);
} // end df_d_yfoc

// derivative of fecundity wrt to survival help, 3rd argument 
// of the function f(y,y,x,x)
double InterGenerational::d_f_d_x() const
{ 
    return(-params.C_x_fec);
} // end df_d_yfoc


// focal breeder's survival probability
double InterGenerational::sa(
        double const x_foc, 
        double const x_loc, 
        double const y_foc) const
{
    // survival = 1.0 - 
    return(std::clamp(
                params.surv_max * (1.0 - params.C_x_surv * x_foc 
                                + params.B_surv * x_loc 
                                - params.C_y_surv * y_foc)
            ,0.0, 1.0) // range between 0 and 1 for survival
        ); 
} // end s_a

// total number of newborns at a local patch
double InterGenerational::u(
        double const y_loc 
        ,double const y
        ,double const x_loc
        ,double const x) const
{
    return( 
             (1.0 - params.d) * f(y_loc, y, x_loc)
             +
             params.d * f(y, y, x)
          );
} // end u()

double InterGenerational::d_u_d_y1() const
{
    return((1.0 - params.d) * (d_f_d_y1() + d_f_d_y2()));
}

double InterGenerational::d_u_d_x1() const
{
    return((1.0 - params.d) * d_f_d_x());
}

// derivative of the focal breeder's survival probability
// wrt to the focal's survival help trait
double InterGenerational::d_sa_d_x1() const
{
    return(- params.surv_max * params.C_x_surv);
}

// derivative of the focal breeder's survival probability
// wrt to the focal's fecundity help trait
double InterGenerational::d_sa_d_y() const
{
    return(- params.surv_max * params.C_y_surv);
}

// derivative of the focal breeder's survival probability
// wrt to the local breeders' survival help trait
double InterGenerational::d_sa_d_x2() const
{
    return(params.surv_max * params.B_surv);
}


// partial derivative dW_AD / d x_loc
double InterGenerational::d_WAD_d_xloc() const
{
    // the denominator of parts 1 and 2 of W_AD (i.e., in case the adult survives
    double denom_ad_survives{1.0 + params.n * params.k * u(y,y,x,x)};
    double denom_ad_dies{params.n * u(y,y,x,x)};

    double part1{
        (d_sa_d_x2() * (1.0 + params.n * params.k * u(y,y,x,x))
        - sa(x,x,y) * params.n * params.k * d_u_d_x1())/
        (denom_ad_survives * denom_ad_survives)
    };

    double part2{
         (params.n * (d_sa_d_x1() + d_sa_d_x2()) * params.k * (1.0 - params.d) * f(y,y,x) * denom_ad_survives
             - params.n * sa(x,x,y) * params.k * (1.0 - params.d) * f(y,y,x) * params.n * params.k * d_u_d_x1()) /
             (denom_ad_survives * denom_ad_survives)};
            
    double part3{
        (params.n * -(d_sa_d_x1() + d_sa_d_x2()) * (1.0 - params.d) * f(y,y,x) / denom_ad_dies
            - params.n * (1.0 - sa(x,x,y)) * (1.0 - params.d) * f(y,y,x) * params.n * d_u_d_x1()) /
            (denom_ad_dies * denom_ad_dies)};

    // derivatives wrt to xloc are 0 for part 4 and part 5

    return(part1  + part2 + part3);
} // end d_WAD_d_xloc


double InterGenerational::d_WAD_d_xfoc() const
{
    double denom_ad_survives{1.0 + params.n * params.k * u(y,y,x,x)};
    double denom_ad_dies{params.n * u(y,y,x,x)};

    double part1{d_sa_d_x1() / denom_ad_survives};

    double part2{params.n * sa(x,x,y) * params.k * (1.0 - params.d) * d_f_d_x() / denom_ad_survives};

    double part3{params.n * (1.0 - sa(x,x,y)) * (1.0 - params.d) * d_f_d_x() / denom_ad_dies};

    double part4{params.d * params.n * sa(x,x,y) * params.k * d_f_d_x() / denom_ad_survives};

    double part5{params.d * params.n * (1.0 - sa(x,x,y)) * d_f_d_x() / denom_ad_dies};

    return{part1 + part2 + part3 + part4 + part5};
} // end d_WAD_d_xfoc()

// partial derivative dW_FS / d x_loc
double InterGenerational::d_WFS_d_xloc() const
{
    double denominator{params.n * sa(x,x,y) + params.k * params.n * u(y,y,x,x)};

    double d_denominator_d_xloc{
        params.n * (d_sa_d_x1() + d_sa_d_x2()) + params.n * params.k * d_u_d_x1()};

    double part1{(d_sa_d_x2() * denominator
        - sa(x,x,y) * d_denominator_d_xloc) / (denominator * denominator)};

    double part2{
        -params.n * params.k * (1.0 - params.d) * 
            f(y,y,x) * d_denominator_d_xloc / (denominator * denominator)};

    // there is no 3rd part wrt xloc
    return(part1 + part2);

} // end d_WFS_d_xloc()

// partial derivative dW_AD / d y_loc
double InterGenerational::d_WAD_d_yloc() const
{
    // the denominator of parts 1 and 2 of W_AD (i.e., in case the adult survives
    double denom_ad_survives{1.0 + params.n * params.k * u(y,y,x,x)};
    double denom_ad_dies{params.n * u(y,y,x,x)};

    double part1{
        -sa(x,x,y) * (params.n * params.k * d_u_d_y1()) / (denom_ad_survives * denom_ad_survives)
    };

    double part2{
        // (derivative of numerator) * denominator
        (
         params.n * d_sa_d_y() * params.k * (1.0 - params.d) * f(y,y,x)
            + params.n * sa(x,x,y) * params.k * (1.0 - params.d) * d_f_d_y2()
            ) / denom_ad_survives
         - params.n * sa(x,x,y) * params.k * (1.0 - params.d) * f(y,y,x) * params.n * params.k * d_u_d_y1() / 
            (denom_ad_survives * denom_ad_survives)
    };

    double part3{
        params.n * -d_sa_d_y() * (1.0 - params.d) * f(y,y,x) / denom_ad_dies
            + params.n * (1.0 - sa(x,x,y)) * (1.0 - params.d) * d_f_d_y2() / denom_ad_dies
            - params.n * (1.0 - sa(x,x,y)) * (1.0 - params.d) * f(y,y,x) * params.n * d_u_d_y1() 
    };

    double part4{
        params.d * params.n * sa(x,x,y) * params.k * d_f_d_y2() / denom_ad_survives
    };

    double part5{
        params.d * params.n * (1.0 - sa(x,x,y)) * d_f_d_y2() / denom_ad_dies
    };

    return(part1 + part2 + part3 + part4 + part5);
} // end d_WAD_d_yloc()

double InterGenerational::d_WAD_d_yfoc() const
{
    // the denominator of parts 1 and 2 of W_AD (i.e., in case the adult survives
    double denom_ad_survives{1.0 + params.n * params.k * u(y,y,x,x)};
    double denom_ad_dies{params.n * u(y,y,x,x)};

    double part1{d_sa_d_y() / denom_ad_survives};

    double part2{params.n * sa(x,x,y) * params.k * (1.0 - params.d) * d_f_d_y1() / denom_ad_survives};

    double part3{params.n * (1.0 - sa(x,x,y)) * (1.0 - params.d) * d_f_d_y1() / denom_ad_dies};

    double part4{params.d * params.n * sa(x,x,y) * params.k * d_f_d_y1() / denom_ad_survives};

    double part5{params.d * params.n * (1.0 - sa(x,x,y)) * d_f_d_y1() / denom_ad_dies};

    return{part1 + part2 + part3 + part4 + part5};
} // end d_WAD_d_yfoc()


// partial derivative dW_FS / d y_loc
double InterGenerational::d_WFS_d_yloc() const
{
    double denominator{params.n * sa(x,x,x) + params.k * params.n * u(y,y,x,x)};

    // derivative of the local denominator wrt y1
    double d_denominator_d_yloc{
        (params.n * d_sa_d_y() + params.n * params.k * d_u_d_y1())};

    double part1{
        -sa(x,x,y) * d_denominator_d_yloc /
            (denominator * denominator)};

    double part2{
        (params.n * params.k * (1.0 - params.d) * d_f_d_y2() * denominator - 
            params.n * params.k * (1.0 - params.d) * f(y,y,x) * d_denominator_d_yloc) / 
            (denominator * denominator)};

    double part3{params.n * params.k * params.d * d_f_d_y2() / denominator};

    return(part1 + part2 + part3);
} // end d_WFS_d_yloc

// partial derivative dW_FS / d x_foc
double InterGenerational::d_WFS_d_yfoc() const
{
    double denominator{params.n * sa(x,x,y) + params.k * params.n * u(y,y,x,x)};

    double val{
        d_sa_d_y() / denominator
            +
            params.n * params.k * (1.0 - params.d) * d_f_d_y1() / denominator
            +
            params.n * params.k * params.d * d_f_d_y1() / denominator};

    return(val);
} // end d_WFS_d_yfoc()

// partial derivative dW_FS / d x_foc
double InterGenerational::d_WFS_d_xfoc() const
{
    double denominator{params.n * sa(x,x,y) + params.k * params.n * u(y,y,x,x)};

    double val{
        d_sa_d_x1() / denominator
            +
            params.n * params.k * (1.0 - params.d) * d_f_d_x() / denominator
            +
            params.n * params.k * params.d * d_f_d_x() / denominator
    };

    return(val);
} // end d_WFS_d_xfoc

// probability an adult secures the breeding position
double InterGenerational::q_a()
{
    if (params.type == AD)
    {
        return(params.n * sa(x,x,y) / 
                (1.0 + params.k * params.n * u(y,y,x,x)));
    }

    // full scramble:
    return(params.n * sa(x,x,y) / 
            (params.n * sa(x,x,y) + params.k * params.n * u(y,y,x,x)));

} // end q_a

// probability a locally born juvenile secures a breeding 
// position, conditional upon it being vacated by an adult
double  InterGenerational::ell()
{
    return(params.n * (1.0 - params.d) * f(y,y,x) / 
            (params.k * params.n * u(y,y,x,x)));
} // end ell()

// iteratively solve relatedness coefficients
void InterGenerational::solve_relatedness()
{
    double rd_ad_tplus1{0.0}, qa_val{0.0}, ell_val{0.0};

    qa_val = q_a();
    assert(qa_val >= 0);
    assert(qa_val <= 1.0);

    ell_val = ell();
    assert(ell_val >= 0);
    assert(ell_val <= 1.0);

    for (time_step_eco = 0; time_step_eco <= params.max_eco_time; ++time_step_eco)
    {
        rd_ad_tplus1 = qa_val * qa_val * rd_ad
            + 2 * qa_val * (1.0 - qa_val) * ell_val * (1.0 / params.n + (params.n - 1.0)/ params.n * rd_ad)
                    + (1.0 - qa_val) * (1.0 - qa_val) * ell_val * ell_val * (1.0 / params.n + (params.n - 1.0)/ params.n * rd_ad);

        if (std::fabs(rd_ad_tplus1 - rd_ad) < params.vanish_relatedness)
        {
            break;
        }

        rd_ad = rd_ad_tplus1;

        assert(rd_ad >= 0);
        assert(rd_ad <= 1.0);
    }
} // end solve_relatedness()

// the heart of the code
void InterGenerational::run()
{
    // aux variable to check convergence 
    bool converged;

    // aux variables to store updated values of both traits
    double xtplus1, ytplus1;

    // write the headers to the output file
    write_data_headers();

    // loop over time
    for (time_step_evo = 0; time_step_evo <= params.max_evo_time; ++time_step_evo)
    {
        // solve relatedness coefficients
        solve_relatedness();

        // update value of x and y
        if (params.type == AD)
        {
            xtplus1 = x + params.var *(d_WAD_d_xfoc() + d_WAD_d_xloc() * rd_ad);
            ytplus1 = y + params.var * (d_WAD_d_yfoc() + d_WAD_d_yloc() * rd_ad);
        }
        else
        {
            xtplus1 = x + params.var *(d_WFS_d_xfoc() + d_WFS_d_xloc() * rd_ad);
            ytplus1 = y + params.var * (d_WFS_d_yfoc() + d_WFS_d_yloc() * rd_ad);
        }

        converged = true;

        // check whether survival help has reached an equilibrium
        if (std::fabs(x - xtplus1) > params.convergence_threshold)
        {
            converged = false;
        }
        
        // check whether fecundity help has reached an equilibrium
        if (std::fabs(y - ytplus1) > params.convergence_threshold)
        {
            converged = false;
        }

        if (time_step_evo % params.skip_output == 0)
        {
            write_data();
        }

        x = xtplus1;
        y = ytplus1;


        if (x < 0)
        {
            x = 0;
        }
        
        if (y < 0)
        {
            y = 0;
        }

        // if x and y are not chang
        if (converged)
        {
            break;
        }
    }

    // output parameters at the end of the output file
    write_parameters();
} // end run();

