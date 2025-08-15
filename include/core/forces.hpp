#pragma once

#include "system.hpp"
#include "vec3.hpp"
#include <cmath>

namespace carbon {

class ForceCalculator {
public:
    // physical constants
    static constexpr double COULOMB_CONSTANT = 332.0637;  // kcal*angstrom/(mol*e^2)
    
    ForceCalculator() = default;
    
    // calculate all forces and update atom forces
    void calculate_forces(System& system) {
        system.clear_forces();
        
        calculate_bond_forces(system);
        calculate_angle_forces(system);
        calculate_nonbonded_forces(system);
    }
    
    // calculate harmonic bond forces
    void calculate_bond_forces(System& system) {
        for (const auto& bond : system.bonds) {
            if (bond.i >= system.num_atoms() || bond.j >= system.num_atoms()) {
                continue;
            }
            
            Atom& atom_i = system.atoms[bond.i];
            Atom& atom_j = system.atoms[bond.j];
            
            Vec3 dr = system.minimum_image_distance(atom_i.r, atom_j.r);
            double r = dr.magnitude();
            
            if (r < 1e-12) continue;  // avoid division by zero
            
            // harmonic potential: U = k * (r - r0)^2
            // force magnitude: F = -dU/dr = -2k * (r - r0)
            double k_eff = bond.effective_spring_constant();
            double r0_eff = bond.effective_equilibrium_length();
            double force_magnitude = -2.0 * k_eff * (r - r0_eff);
            
            Vec3 force_direction = dr.normalized();
            Vec3 force = force_magnitude * force_direction;
            
            atom_i.add_force(-force);  // newton's third law
            atom_j.add_force(force);
        }
    }
    
    // calculate harmonic angle forces
    void calculate_angle_forces(System& system) {
        for (const auto& angle : system.angles) {
            if (angle.i >= system.num_atoms() || angle.j >= system.num_atoms() || 
                angle.k >= system.num_atoms()) {
                continue;
            }
            
            Atom& atom_i = system.atoms[angle.i];
            Atom& atom_j = system.atoms[angle.j];  // central atom
            Atom& atom_k = system.atoms[angle.k];
            
            Vec3 rij = system.minimum_image_distance(atom_j.r, atom_i.r);
            Vec3 rkj = system.minimum_image_distance(atom_j.r, atom_k.r);
            
            double rij_mag = rij.magnitude();
            double rkj_mag = rkj.magnitude();
            
            if (rij_mag < 1e-12 || rkj_mag < 1e-12) continue;
            
            // calculate current angle
            double cos_theta = rij.dot(rkj) / (rij_mag * rkj_mag);
            cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
            double theta = std::acos(cos_theta);
            
            // harmonic angle potential: U = ka * (theta - theta0)^2
            // force magnitude: F = -dU/dtheta = -2 * ka * (theta - theta0)
            double force_magnitude = -2.0 * angle.ka * (theta - angle.theta0);
            
            // avoid singularity at theta = 0 or pi
            double sin_theta = std::sin(theta);
            if (std::abs(sin_theta) < 1e-8) continue;
            
            // force directions
            Vec3 ui = (rij / rij_mag - (rij.dot(rkj) / (rij_mag * rkj_mag * rkj_mag)) * rkj) / sin_theta;
            Vec3 uk = (rkj / rkj_mag - (rij.dot(rkj) / (rij_mag * rij_mag * rkj_mag)) * rij) / sin_theta;
            
            Vec3 force_i = force_magnitude / rij_mag * ui;
            Vec3 force_k = force_magnitude / rkj_mag * uk;
            Vec3 force_j = -(force_i + force_k);
            
            atom_i.add_force(force_i);
            atom_j.add_force(force_j);
            atom_k.add_force(force_k);
        }
    }
    
    // calculate lennard-jones and coulomb forces
    void calculate_nonbonded_forces(System& system) {
        for (int i = 0; i < system.num_atoms(); ++i) {
            for (int j = i + 1; j < system.num_atoms(); ++j) {
                if (are_bonded(system, i, j)) {
                    continue;  // skip bonded pairs
                }
                
                calculate_lj_force(system, i, j);
                calculate_coulomb_force(system, i, j);
            }
        }
    }
    
    // calculate lennard-jones force between two atoms
    void calculate_lj_force(System& system, int i, int j) {
        Atom& atom_i = system.atoms[i];
        Atom& atom_j = system.atoms[j];
        
        Vec3 dr = system.minimum_image_distance(atom_i.r, atom_j.r);
        double r = dr.magnitude();
        
        if (r < 1e-12) return;  // avoid division by zero
        
        // combine lj parameters using lorentz-berthelot rules
        double sigma_ij = 0.5 * (atom_i.sigma + atom_j.sigma);
        double epsilon_ij = std::sqrt(atom_i.epsilon * atom_j.epsilon);
        
        if (epsilon_ij < 1e-12) return;  // no interaction
        
        double sigma_over_r = sigma_ij / r;
        double sigma6 = std::pow(sigma_over_r, 6);
        double sigma12 = sigma6 * sigma6;
        
        // lj potential: U = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
        // force magnitude: F = -dU/dr = 24*epsilon/r * [2*(sigma/r)^12 - (sigma/r)^6]
        double force_magnitude = 24.0 * epsilon_ij / r * (2.0 * sigma12 - sigma6);
        
        Vec3 force_direction = dr.normalized();
        Vec3 force = force_magnitude * force_direction;
        
        atom_i.add_force(-force);
        atom_j.add_force(force);
    }
    
    // calculate coulomb force between two atoms
    void calculate_coulomb_force(System& system, int i, int j) {
        Atom& atom_i = system.atoms[i];
        Atom& atom_j = system.atoms[j];
        
        if (std::abs(atom_i.charge) < 1e-12 || std::abs(atom_j.charge) < 1e-12) {
            return;  // no charge
        }
        
        Vec3 dr = system.minimum_image_distance(atom_i.r, atom_j.r);
        double r = dr.magnitude();
        
        if (r < 1e-12) return;  // avoid division by zero
        
        // coulomb potential: U = k * q1 * q2 / r
        // force magnitude: F = -dU/dr = k * q1 * q2 / r^2
        double force_magnitude = COULOMB_CONSTANT * atom_i.charge * atom_j.charge / (r * r);
        
        Vec3 force_direction = dr.normalized();
        Vec3 force = force_magnitude * force_direction;
        
        atom_i.add_force(-force);
        atom_j.add_force(force);
    }
    
    // energy calculations
    double calculate_bond_energy(const System& system) const {
        double energy = 0.0;
        
        for (const auto& bond : system.bonds) {
            if (bond.i >= system.num_atoms() || bond.j >= system.num_atoms()) {
                continue;
            }
            
            const Atom& atom_i = system.atoms[bond.i];
            const Atom& atom_j = system.atoms[bond.j];
            
            Vec3 dr = system.minimum_image_distance(atom_i.r, atom_j.r);
            double r = dr.magnitude();
            
            double k_eff = bond.effective_spring_constant();
            double r0_eff = bond.effective_equilibrium_length();
            double dr_eq = r - r0_eff;
            energy += k_eff * dr_eq * dr_eq;
        }
        
        return energy;
    }
    
    double calculate_angle_energy(const System& system) const {
        double energy = 0.0;
        
        for (const auto& angle : system.angles) {
            if (angle.i >= system.num_atoms() || angle.j >= system.num_atoms() || 
                angle.k >= system.num_atoms()) {
                continue;
            }
            
            double current_angle = system.angle_between_atoms(angle.i, angle.j, angle.k);
            double delta_theta = current_angle - angle.theta0;
            energy += angle.ka * delta_theta * delta_theta;
        }
        
        return energy;
    }
    
    double calculate_lj_energy(const System& system) const {
        double energy = 0.0;
        
        for (int i = 0; i < system.num_atoms(); ++i) {
            for (int j = i + 1; j < system.num_atoms(); ++j) {
                if (are_bonded(system, i, j)) {
                    continue;
                }
                
                const Atom& atom_i = system.atoms[i];
                const Atom& atom_j = system.atoms[j];
                
                Vec3 dr = system.minimum_image_distance(atom_i.r, atom_j.r);
                double r = dr.magnitude();
                
                if (r < 1e-12) continue;
                
                double sigma_ij = 0.5 * (atom_i.sigma + atom_j.sigma);
                double epsilon_ij = std::sqrt(atom_i.epsilon * atom_j.epsilon);
                
                if (epsilon_ij < 1e-12) continue;
                
                double sigma_over_r = sigma_ij / r;
                double sigma6 = std::pow(sigma_over_r, 6);
                double sigma12 = sigma6 * sigma6;
                
                energy += 4.0 * epsilon_ij * (sigma12 - sigma6);
            }
        }
        
        return energy;
    }
    
    double calculate_coulomb_energy(const System& system) const {
        double energy = 0.0;
        
        for (int i = 0; i < system.num_atoms(); ++i) {
            for (int j = i + 1; j < system.num_atoms(); ++j) {
                if (are_bonded(system, i, j)) {
                    continue;
                }
                
                const Atom& atom_i = system.atoms[i];
                const Atom& atom_j = system.atoms[j];
                
                if (std::abs(atom_i.charge) < 1e-12 || std::abs(atom_j.charge) < 1e-12) {
                    continue;
                }
                
                Vec3 dr = system.minimum_image_distance(atom_i.r, atom_j.r);
                double r = dr.magnitude();
                
                if (r < 1e-12) continue;
                
                energy += COULOMB_CONSTANT * atom_i.charge * atom_j.charge / r;
            }
        }
        
        return energy;
    }
    
    double calculate_total_potential_energy(const System& system) const {
        return calculate_bond_energy(system) + 
               calculate_angle_energy(system) +
               calculate_lj_energy(system) + 
               calculate_coulomb_energy(system);
    }

private:
    bool are_bonded(const System& system, int i, int j) const {
        for (const auto& bond : system.bonds) {
            if ((bond.i == i && bond.j == j) || (bond.i == j && bond.j == i)) {
                return true;
            }
        }
        return false;
    }
};

}