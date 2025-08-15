#pragma once

#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "system.hpp"
#include "forces.hpp"
#include "vec3.hpp"

namespace carbon {

class VelocityVerletIntegrator {
private:
    double dt;  // timestep in femtoseconds
    ForceCalculator force_calc;
    
    // conversion factors
    static constexpr double AMU_TO_KG = 1.66054e-27;        // amu to kg
    static constexpr double FS_TO_S = 1e-15;                // fs to s
    static constexpr double ANGSTROM_TO_M = 1e-10;          // angstrom to m
    static constexpr double KCAL_TO_J = 4184.0;             // kcal to j
    static constexpr double AVOGADRO = 6.02214076e23;       // avogadro's number
    
    // conversion factor for acceleration: (kcal/mol/angstrom) / amu -> angstrom/fs^2
    static constexpr double FORCE_TO_ACCEL = KCAL_TO_J / (AVOGADRO * AMU_TO_KG * ANGSTROM_TO_M) * (FS_TO_S * FS_TO_S);

public:
    explicit VelocityVerletIntegrator(double timestep = 1.0) : dt(timestep) {}
    
    void set_timestep(double timestep) {
        dt = timestep;
    }
    
    double get_timestep() const {
        return dt;
    }
    
    // perform one integration step
    void step(System& system) {
        const int n_atoms = system.num_atoms();
        
        // store previous accelerations
        std::vector<Vec3> prev_accel(n_atoms);
        for (int i = 0; i < n_atoms; ++i) {
            prev_accel[i] = system.atoms[i].f / system.atoms[i].mass * FORCE_TO_ACCEL;
        }
        
        // update positions: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
        for (int i = 0; i < n_atoms; ++i) {
            Atom& atom = system.atoms[i];
            atom.r += atom.v * dt + 0.5 * prev_accel[i] * (dt * dt);
        }
        
        // apply periodic boundary conditions if enabled
        system.apply_pbc();
        
        // calculate new forces
        force_calc.calculate_forces(system);
        
        // calculate new accelerations
        std::vector<Vec3> new_accel(n_atoms);
        for (int i = 0; i < n_atoms; ++i) {
            new_accel[i] = system.atoms[i].f / system.atoms[i].mass * FORCE_TO_ACCEL;
        }
        
        // update velocities: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
        for (int i = 0; i < n_atoms; ++i) {
            Atom& atom = system.atoms[i];
            atom.v += 0.5 * (prev_accel[i] + new_accel[i]) * dt;
        }
    }
    
    // perform multiple integration steps
    void run(System& system, int n_steps) {
        for (int step = 0; step < n_steps; ++step) {
            this->step(system);
        }
    }
    
    // calculate system temperature from kinetic energy
    double calculate_temperature(const System& system) const {
        if (system.num_atoms() == 0) return 0.0;
        
        double total_ke = system.total_kinetic_energy();
        int degrees_of_freedom = 3 * system.num_atoms();  // 3n dof for n atoms
        
        // temperature from equipartition theorem: 0.5 * k_B * T per dof
        // convert from kcal/mol to kelvin
        constexpr double R_GAS_CONSTANT = 1.987e-3;  // kcal/(mol*k)
        
        return (2.0 * total_ke) / (degrees_of_freedom * R_GAS_CONSTANT);
    }
    
    // scale velocities to set temperature
    void set_temperature(System& system, double target_temp) {
        if (target_temp <= 0.0) return;
        
        double current_temp = calculate_temperature(system);
        if (current_temp <= 0.0) return;
        
        double scale_factor = std::sqrt(target_temp / current_temp);
        
        for (auto& atom : system.atoms) {
            atom.v *= scale_factor;
        }
    }
    
    // initialize random velocities based on maxwell-boltzmann distribution
    void initialize_velocities(System& system, double temperature, unsigned int seed = 0) {
        if (seed == 0) {
            seed = static_cast<unsigned int>(std::time(nullptr));
        }
        std::srand(seed);
        
        constexpr double R_GAS_CONSTANT = 1.987e-3;  // kcal/(mol*k)
        
        for (auto& atom : system.atoms) {
            // maxwell-boltzmann velocity distribution
            // <v^2> = k_B*T/m = R*T/m (in our units)
            double variance = R_GAS_CONSTANT * temperature / atom.mass;
            double sigma = std::sqrt(variance);
            
            // box-muller transform for gaussian random numbers
            double u1 = (std::rand() + 1.0) / (RAND_MAX + 2.0);
            double u2 = (std::rand() + 1.0) / (RAND_MAX + 2.0);
            double z0 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            
            u1 = (std::rand() + 1.0) / (RAND_MAX + 2.0);
            u2 = (std::rand() + 1.0) / (RAND_MAX + 2.0);
            double z1 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            
            u1 = (std::rand() + 1.0) / (RAND_MAX + 2.0);
            u2 = (std::rand() + 1.0) / (RAND_MAX + 2.0);
            double z2 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            
            atom.v = Vec3(z0 * sigma, z1 * sigma, z2 * sigma);
        }
        
        // remove net momentum
        Vec3 total_momentum;
        double total_mass = 0.0;
        for (const auto& atom : system.atoms) {
            total_momentum += atom.v * atom.mass;
            total_mass += atom.mass;
        }
        
        Vec3 cm_velocity = total_momentum / total_mass;
        for (auto& atom : system.atoms) {
            atom.v -= cm_velocity;
        }
    }
    
    ForceCalculator& get_force_calculator() {
        return force_calc;
    }
    
    const ForceCalculator& get_force_calculator() const {
        return force_calc;
    }
};

}