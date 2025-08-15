#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include "atom.hpp"
#include "bond.hpp"
#include "angle.hpp"
#include "vec3.hpp"

namespace carbon {

class System {
public:
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    
    // simulation box info
    Vec3 box_size;
    bool periodic;
    
    System() : box_size(50.0, 50.0, 50.0), periodic(false) {}
    
    explicit System(const Vec3& box) : box_size(box), periodic(true) {}
    
    // atom management
    int add_atom(const Atom& atom) {
        atoms.push_back(atom);
        return atoms.size() - 1;
    }
    
    int add_atom(int Z, double mass, const Vec3& position = Vec3()) {
        atoms.emplace_back(Z, mass, 0.0, 0.0, 0.0, position);
        return atoms.size() - 1;
    }
    
    void remove_atom(int idx) {
        if (idx >= 0 && idx < static_cast<int>(atoms.size())) {
            atoms.erase(atoms.begin() + idx);
            // update bond and angle indices
            remove_bonds_involving_atom(idx);
            remove_angles_involving_atom(idx);
            update_bond_indices_after_removal(idx);
            update_angle_indices_after_removal(idx);
        }
    }
    
    // bond management
    void add_bond(int i, int j, double k, double r0) {
        if (i >= 0 && i < static_cast<int>(atoms.size()) && 
            j >= 0 && j < static_cast<int>(atoms.size()) && i != j) {
            bonds.emplace_back(i, j, k, r0);
        }
    }
    
    void add_bond(const Bond& bond) {
        bonds.push_back(bond);
    }
    
    void remove_bond(int i, int j) {
        bonds.erase(
            std::remove_if(bonds.begin(), bonds.end(),
                [i, j](const Bond& bond) {
                    return (bond.i == i && bond.j == j) || (bond.i == j && bond.j == i);
                }),
            bonds.end()
        );
    }
    
    // angle management
    void add_angle(int i, int j, int k, double ka, double theta0) {
        if (i >= 0 && i < static_cast<int>(atoms.size()) && 
            j >= 0 && j < static_cast<int>(atoms.size()) &&
            k >= 0 && k < static_cast<int>(atoms.size()) &&
            i != j && j != k && i != k) {
            angles.emplace_back(i, j, k, ka, theta0);
        }
    }
    
    void add_angle(const Angle& angle) {
        if (angle.is_valid()) {
            angles.push_back(angle);
        }
    }
    
    void remove_angle(int i, int j, int k) {
        angles.erase(
            std::remove_if(angles.begin(), angles.end(),
                [i, j, k](const Angle& angle) {
                    return (angle.i == i && angle.j == j && angle.k == k) ||
                           (angle.i == k && angle.j == j && angle.k == i);
                }),
            angles.end()
        );
    }
    
    // system queries
    int num_atoms() const { return static_cast<int>(atoms.size()); }
    int num_bonds() const { return static_cast<int>(bonds.size()); }
    int num_angles() const { return static_cast<int>(angles.size()); }
    
    bool is_empty() const { return atoms.empty(); }
    
    Vec3 center_of_mass() const {
        if (atoms.empty()) return Vec3();
        
        Vec3 com;
        double total_mass = 0.0;
        
        for (const auto& atom : atoms) {
            com += atom.r * atom.mass;
            total_mass += atom.mass;
        }
        
        return com / total_mass;
    }
    
    double total_kinetic_energy() const {
        double ke = 0.0;
        for (const auto& atom : atoms) {
            ke += atom.kinetic_energy();
        }
        return ke;
    }
    
    // force clearing
    void clear_forces() {
        for (auto& atom : atoms) {
            atom.clear_forces();
        }
    }
    
    // periodic boundary conditions
    void apply_pbc() {
        if (!periodic) return;
        
        for (auto& atom : atoms) {
            atom.r = apply_pbc(atom.r, box_size);
        }
    }
    
    Vec3 minimum_image_distance(const Vec3& r1, const Vec3& r2) const {
        if (periodic) {
            return minimum_image_vector(r1, r2, box_size);
        } else {
            return r2 - r1;
        }
    }
    
    double distance_between_atoms(int i, int j) const {
        if (i < 0 || i >= num_atoms() || j < 0 || j >= num_atoms()) {
            return 0.0;
        }
        
        Vec3 dr = minimum_image_distance(atoms[i].r, atoms[j].r);
        return dr.magnitude();
    }
    
    double angle_between_atoms(int i, int j, int k) const {
        if (i < 0 || i >= num_atoms() || j < 0 || j >= num_atoms() || 
            k < 0 || k >= num_atoms()) {
            return 0.0;
        }
        
        Vec3 rij = minimum_image_distance(atoms[j].r, atoms[i].r);
        Vec3 rkj = minimum_image_distance(atoms[j].r, atoms[k].r);
        
        double dot_product = rij.dot(rkj);
        double mag_i = rij.magnitude();
        double mag_k = rkj.magnitude();
        
        if (mag_i < 1e-12 || mag_k < 1e-12) return 0.0;
        
        double cos_theta = dot_product / (mag_i * mag_k);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));  // clamp to [-1,1]
        
        return std::acos(cos_theta);
    }
    
    // system information
    void print_info() const {
        std::cout << "System info:\n";
        std::cout << "  Atoms: " << num_atoms() << "\n";
        std::cout << "  Bonds: " << num_bonds() << "\n";
        std::cout << "  Angles: " << num_angles() << "\n";
        std::cout << "  Box size: " << box_size << "\n";
        std::cout << "  Periodic: " << (periodic ? "yes" : "no") << "\n";
        std::cout << "  Total kinetic energy: " << total_kinetic_energy() << " kcal/mol\n";
    }

private:
    void remove_bonds_involving_atom(int atom_idx) {
        bonds.erase(
            std::remove_if(bonds.begin(), bonds.end(),
                [atom_idx](const Bond& bond) {
                    return bond.involves_atom(atom_idx);
                }),
            bonds.end()
        );
    }
    
    void update_bond_indices_after_removal(int removed_idx) {
        for (auto& bond : bonds) {
            if (bond.i > removed_idx) bond.i--;
            if (bond.j > removed_idx) bond.j--;
        }
    }
    
    void remove_angles_involving_atom(int atom_idx) {
        angles.erase(
            std::remove_if(angles.begin(), angles.end(),
                [atom_idx](const Angle& angle) {
                    return angle.involves_atom(atom_idx);
                }),
            angles.end()
        );
    }
    
    void update_angle_indices_after_removal(int removed_idx) {
        for (auto& angle : angles) {
            if (angle.i > removed_idx) angle.i--;
            if (angle.j > removed_idx) angle.j--;
            if (angle.k > removed_idx) angle.k--;
        }
    }
};

}