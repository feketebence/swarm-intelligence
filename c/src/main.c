

//
//  

#include <stdio.h>
#include <stdlib.h>

#include "initializer.h"
#include "running.h"

int main(int argc, const char * argv[])
{
    printf("Swarm Intelligence Course, 2021\n");
    printf("#1 Assignment Molecular Dynamics (Brownian Dynamics) simulation\n");
    printf("Self Organizing Systemn on Simple Rules\n");
   
    init_simulation();
    init_simulation_box();

    init_particles();
    init_files();
    
    run_simulation();
    
    return 0;
}
