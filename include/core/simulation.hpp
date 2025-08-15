#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "system.hpp"
#include "integrator.hpp"
#include "forces.hpp"
#include "../io/xyz_io.hpp"

namespace carbon {

struct SimulationParameters {
    double timestep = 1.0;           // integration timestep (fs)
    int total_steps = 1000;          // total simulation steps
    int output_frequency = 10;       // write trajectory every n steps
    int energy_frequency = 1;        // write energy every n steps
    double temperature = 300.0;      // target temperature (k)
    bool constant_temperature = false; // enable temperature coupling
    std::string trajectory_file = "trajectory.xyz";
    std::string energy_file = "energy.log";
    bool verbose = true;             // print progress
};

class Simulation {
private:
    System system;
    VelocityVerletIntegrator integrator;
    SimulationParameters params;
    
    std::ofstream energy_log;
    TrajectoryWriter trajectory;
    
    int current_step;
    double current_time;
    
public:
    explicit Simulation(const SimulationParameters& parameters = SimulationParameters())
        : params(parameters), trajectory(params.trajectory_file), 
          current_step(0), current_time(0.0) {
        integrator.set_timestep(params.timestep);
    }
    
    System& get_system() { return system; }
    const System& get_system() const { return system; }
    
    VelocityVerletIntegrator& get_integrator() { return integrator; }
    
    void set_parameters(const SimulationParameters& parameters) {
        params = parameters;
        integrator.set_timestep(params.timestep);
    }
    
    bool initialize() {
        if (system.num_atoms() == 0) {
            std::cerr << "error: no atoms in system\n";
            return false;
        }
        
        // open trajectory file
        if (!trajectory.open()) {
            std::cerr << "error: cannot open trajectory file " << params.trajectory_file << "\n";
            return false;
        }
        
        // open energy log
        energy_log.open(params.energy_file);
        if (!energy_log.is_open()) {
            std::cerr << "error: cannot open energy file " << params.energy_file << "\n";
            return false;
        }
        
        // write energy header
        energy_log << "# molecular dynamics energy log\n";
        energy_log << "# columns: step time(fs) kinetic(kcal/mol) potential(kcal/mol) total(kcal/mol) temperature(k)\n";
        energy_log << std::fixed << std::setprecision(6);
        
        // initialize velocities if needed
        bool has_velocities = false;
        for (const auto& atom : system.atoms) {
            if (atom.v.magnitude() > 1e-12) {
                has_velocities = true;
                break;
            }
        }
        
        if (!has_velocities && params.temperature > 0.0) {
            integrator.initialize_velocities(system, params.temperature);
            if (params.verbose) {
                std::cout << "initialized velocities for temperature " << params.temperature << " k\n";
            }
        }
        
        // calculate initial forces
        integrator.get_force_calculator().calculate_forces(system);
        
        current_step = 0;
        current_time = 0.0;
        
        // write initial frame
        write_trajectory_frame();
        write_energy_log();
        
        if (params.verbose) {
            std::cout << "simulation initialized\n";
            std::cout << "  atoms: " << system.num_atoms() << "\n";
            std::cout << "  bonds: " << system.num_bonds() << "\n";
            std::cout << "  timestep: " << params.timestep << " fs\n";
            std::cout << "  total steps: " << params.total_steps << "\n";
            std::cout << "  total time: " << (params.total_steps * params.timestep) << " fs\n";
            std::cout << "  trajectory file: " << params.trajectory_file << "\n";
            std::cout << "  energy file: " << params.energy_file << "\n\n";
        }
        
        return true;
    }
    
    void run() {
        if (params.verbose) {
            std::cout << "starting molecular dynamics simulation...\n";
        }
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        for (int step = 1; step <= params.total_steps; ++step) {
            current_step = step;
            current_time = step * params.timestep;
            
            // temperature coupling
            if (params.constant_temperature && params.temperature > 0.0) {
                integrator.set_temperature(system, params.temperature);
            }
            
            // integration step
            integrator.step(system);
            
            // output
            if (step % params.output_frequency == 0) {
                write_trajectory_frame();
            }
            
            if (step % params.energy_frequency == 0) {
                write_energy_log();
            }
            
            // progress reporting
            if (params.verbose && step % (params.total_steps / 10) == 0) {
                double progress = 100.0 * step / params.total_steps;
                double temperature = integrator.calculate_temperature(system);
                double total_energy = get_total_energy();
                
                std::cout << "step " << step << "/" << params.total_steps 
                          << " (" << std::fixed << std::setprecision(1) << progress << "%) "
                          << "t=" << std::setprecision(3) << current_time << " fs, "
                          << "temp=" << std::setprecision(1) << temperature << " k, "
                          << "energy=" << std::setprecision(3) << total_energy << " kcal/mol\n";
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        if (params.verbose) {
            std::cout << "\nsimulation completed in " << duration.count() << " ms\n";
            std::cout << "performance: " << std::fixed << std::setprecision(2) 
                      << (params.total_steps * 1000.0 / duration.count()) << " steps/second\n";
            
            print_final_statistics();
        }
        
        finalize();
    }
    
    void run_single_step() {
        current_step++;
        current_time = current_step * params.timestep;
        
        if (params.constant_temperature && params.temperature > 0.0) {
            integrator.set_temperature(system, params.temperature);
        }
        
        integrator.step(system);
        
        if (current_step % params.output_frequency == 0) {
            write_trajectory_frame();
        }
        
        if (current_step % params.energy_frequency == 0) {
            write_energy_log();
        }
    }
    
    double get_kinetic_energy() const {
        return system.total_kinetic_energy();
    }
    
    double get_potential_energy() const {
        return integrator.get_force_calculator().calculate_total_potential_energy(system);
    }
    
    double get_total_energy() const {
        return get_kinetic_energy() + get_potential_energy();
    }
    
    double get_temperature() const {
        return integrator.calculate_temperature(system);
    }
    
    int get_current_step() const {
        return current_step;
    }
    
    double get_current_time() const {
        return current_time;
    }
    
    void finalize() {
        trajectory.close();
        energy_log.close();
        
        if (params.verbose) {
            std::cout << "simulation files closed\n";
        }
    }

private:
    void write_trajectory_frame() {
        double total_energy = get_total_energy();
        trajectory.write_frame(system, current_time, total_energy);
    }
    
    void write_energy_log() {
        double ke = get_kinetic_energy();
        double pe = get_potential_energy();
        double total = ke + pe;
        double temp = get_temperature();
        
        energy_log << std::setw(8) << current_step << " "
                   << std::setw(12) << current_time << " "
                   << std::setw(15) << ke << " "
                   << std::setw(15) << pe << " "
                   << std::setw(15) << total << " "
                   << std::setw(12) << temp << "\n";
        energy_log.flush();
    }
    
    void print_final_statistics() {
        double final_ke = get_kinetic_energy();
        double final_pe = get_potential_energy();
        double final_total = final_ke + final_pe;
        double final_temp = get_temperature();
        
        std::cout << "\nfinal system state:\n";
        std::cout << "  kinetic energy: " << std::fixed << std::setprecision(6) << final_ke << " kcal/mol\n";
        std::cout << "  potential energy: " << final_pe << " kcal/mol\n";
        std::cout << "  total energy: " << final_total << " kcal/mol\n";
        std::cout << "  temperature: " << std::setprecision(1) << final_temp << " k\n";
        std::cout << "  trajectory frames: " << trajectory.get_frame_count() << "\n";
        
        Vec3 com = system.center_of_mass();
        std::cout << "  center of mass: (" << std::setprecision(3) 
                  << com.x << ", " << com.y << ", " << com.z << ")\n";
    }
};

}