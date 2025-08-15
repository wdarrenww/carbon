#include <iostream>
#include <string>
#include "../include/carbon.hpp"

using namespace carbon;

void run_water_simulation() {
    std::cout << "creating water molecule simulation...\n";
    
    // create water molecule
    System system = molecules::create_water_molecule();
    
    // set up simulation parameters
    SimulationParameters params;
    params.timestep = 0.5;           // 0.5 fs timestep
    params.total_steps = 1000;       // 1000 steps = 0.5 ps
    params.output_frequency = 10;    // save every 10 steps
    params.energy_frequency = 1;     // log energy every step
    params.temperature = 300.0;      // 300 k
    params.trajectory_file = "water.xyz";
    params.energy_file = "water_energy.log";
    params.verbose = true;
    
    // create and run simulation
    Simulation sim(params);
    sim.get_system() = system;
    
    if (sim.initialize()) {
        sim.run();
        std::cout << "water simulation completed. check water.xyz and water_energy.log\n\n";
    } else {
        std::cout << "failed to initialize water simulation\n\n";
    }
}

void run_methane_simulation() {
    std::cout << "creating methane molecule simulation...\n";
    
    // create methane molecule
    System system = molecules::create_methane_molecule();
    
    // set up simulation parameters
    SimulationParameters params;
    params.timestep = 1.0;           // 1.0 fs timestep
    params.total_steps = 500;        // 500 steps = 0.5 ps
    params.output_frequency = 5;     // save every 5 steps
    params.energy_frequency = 1;     // log energy every step
    params.temperature = 300.0;      // 300 k
    params.trajectory_file = "methane.xyz";
    params.energy_file = "methane_energy.log";
    params.verbose = true;
    
    // create and run simulation
    Simulation sim(params);
    sim.get_system() = system;
    
    if (sim.initialize()) {
        sim.run();
        std::cout << "methane simulation completed. check methane.xyz and methane_energy.log\n\n";
    } else {
        std::cout << "failed to initialize methane simulation\n\n";
    }
}

void run_nacl_simulation() {
    std::cout << "creating nacl dimer simulation...\n";
    
    // create nacl dimer
    System system = molecules::create_nacl_dimer();
    
    // set up simulation parameters
    SimulationParameters params;
    params.timestep = 1.0;           // 1.0 fs timestep
    params.total_steps = 1000;       // 1000 steps = 1.0 ps
    params.output_frequency = 10;    // save every 10 steps
    params.energy_frequency = 1;     // log energy every step
    params.temperature = 300.0;      // 300 k
    params.trajectory_file = "nacl.xyz";
    params.energy_file = "nacl_energy.log";
    params.verbose = true;
    
    // create and run simulation
    Simulation sim(params);
    sim.get_system() = system;
    
    if (sim.initialize()) {
        sim.run();
        std::cout << "nacl simulation completed. check nacl.xyz and nacl_energy.log\n\n";
    } else {
        std::cout << "failed to initialize nacl simulation\n\n";
    }
}

void interactive_mode() {
    std::cout << "entering interactive building mode...\n\n";
    
    System system;
    InteractiveBuilder::build_system(system);
    
    if (system.num_atoms() > 0) {
        std::cout << "\nsystem built with " << system.num_atoms() << " atoms and " 
                  << system.num_bonds() << " bonds\n";
        
        std::cout << "run simulation? (y/n): ";
        std::string response;
        std::getline(std::cin, response);
        
        if (response == "y" || response == "yes" || response == "Y") {
            SimulationParameters params;
            params.timestep = 1.0;
            params.total_steps = 1000;
            params.output_frequency = 10;
            params.energy_frequency = 1;
            params.temperature = 300.0;
            params.trajectory_file = "custom.xyz";
            params.energy_file = "custom_energy.log";
            params.verbose = true;
            
            Simulation sim(params);
            sim.get_system() = system;
            
            if (sim.initialize()) {
                sim.run();
                std::cout << "custom simulation completed. check custom.xyz and custom_energy.log\n";
            }
        }
    }
}

void print_usage() {
    std::cout << "usage: carbon [mode]\n";
    std::cout << "modes:\n";
    std::cout << "  water     - run water molecule simulation\n";
    std::cout << "  methane   - run methane molecule simulation\n";
    std::cout << "  nacl      - run nacl dimer simulation\n";
    std::cout << "  build     - interactive molecule building\n";
    std::cout << "  demo      - run all demo simulations\n";
    std::cout << "  help      - show this help\n";
}

int main(int argc, char* argv[]) {
    print_banner();
    
    std::string mode = "demo";
    if (argc > 1) {
        mode = argv[1];
    }
    
    if (mode == "help" || mode == "-h" || mode == "--help") {
        print_usage();
    } else if (mode == "water") {
        run_water_simulation();
    } else if (mode == "methane") {
        run_methane_simulation();
    } else if (mode == "nacl") {
        run_nacl_simulation();
    } else if (mode == "build") {
        interactive_mode();
    } else if (mode == "demo") {
        std::cout << "running demo simulations...\n\n";
        run_water_simulation();
        run_methane_simulation();
        run_nacl_simulation();
        std::cout << "all demo simulations completed!\n";
    } else {
        std::cout << "unknown mode: " << mode << "\n";
        print_usage();
        return 1;
    }
    
    return 0;
}