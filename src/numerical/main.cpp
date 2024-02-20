#include "parameters.hpp"
#include "intergenerational.hpp"

int main(int argc, char **argv)
{
    Parameters params;

    InterGenerational simulation_object(params);

    simulation_object.run();

    return 0;
}
