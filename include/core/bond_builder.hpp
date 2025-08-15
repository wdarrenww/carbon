#pragma once

#include "system.hpp"
#include "elements.hpp"
#include <vector>

namespace carbon {

class BondBuilder {
public:
    struct BondParameters {
        double bond_tolerance = 1.2;    // factor to multiply sum of covalent radii
        double max_bond_length = 3.0;   // maximum bond length in angstroms
        bool auto_angles = true;        // automatically generate angle terms
        double angle_tolerance = 0.1;   // tolerance for angle detection
        
        // default bond parameters by element pairs
        struct DefaultBondParams {
            double k = 300.0;           // spring constant (kcal/mol/angstrom^2)
            double order = 1.0;         // bond order
        };
        
        DefaultBondParams get_default_params(int z1, int z2) const {
            DefaultBondParams params;
            
            // c-h bonds
            if ((z1 == 6 && z2 == 1) || (z1 == 1 && z2 == 6)) {
                params.k = 340.0;
            }
            // o-h bonds  
            else if ((z1 == 8 && z2 == 1) || (z1 == 1 && z2 == 8)) {
                params.k = 450.0;
            }
            // n-h bonds
            else if ((z1 == 7 && z2 == 1) || (z1 == 1 && z2 == 7)) {
                params.k = 410.0;
            }
            // c-c bonds
            else if (z1 == 6 && z2 == 6) {
                params.k = 320.0;
            }
            // c-o bonds
            else if ((z1 == 6 && z2 == 8) || (z1 == 8 && z2 == 6)) {
                params.k = 360.0;
            }
            // c-n bonds
            else if ((z1 == 6 && z2 == 7) || (z1 == 7 && z2 == 6)) {
                params.k = 350.0;
            }
            
            return params;
        }
    };
    
    static void detect_bonds(System& system, const BondParameters& params = BondParameters()) {
        const ElementTable& table = get_element_table();
        
        // clear existing bonds
        system.bonds.clear();
        
        // detect bonds between all atom pairs
        for (int i = 0; i < system.num_atoms(); ++i) {
            for (int j = i + 1; j < system.num_atoms(); ++j) {
                const Atom& atom_i = system.atoms[i];
                const Atom& atom_j = system.atoms[j];
                
                // get covalent radii
                double r_cov_i = table.get_covalent_radius(atom_i.Z);
                double r_cov_j = table.get_covalent_radius(atom_j.Z);
                
                if (r_cov_i < 1e-12 || r_cov_j < 1e-12) continue;
                
                // calculate distance
                double distance = system.distance_between_atoms(i, j);
                
                // check if atoms are bonded
                double bond_cutoff = params.bond_tolerance * (r_cov_i + r_cov_j);
                
                if (distance > 0.1 && distance < bond_cutoff && distance < params.max_bond_length) {
                    // get bond parameters
                    auto bond_params = params.get_default_params(atom_i.Z, atom_j.Z);
                    
                    // add bond
                    system.add_bond(i, j, bond_params.k, distance);
                }
            }
        }
        
        // detect angles if requested
        if (params.auto_angles) {
            detect_angles(system, params);
        }
    }
    
    static void detect_angles(System& system, const BondParameters& params = BondParameters()) {
        // clear existing angles
        system.angles.clear();
        
        // find all angles formed by bonds
        for (int j = 0; j < system.num_atoms(); ++j) {  // central atom
            std::vector<int> bonded_atoms;
            
            // find all atoms bonded to j
            for (const auto& bond : system.bonds) {
                if (bond.i == j) {
                    bonded_atoms.push_back(bond.j);
                } else if (bond.j == j) {
                    bonded_atoms.push_back(bond.i);
                }
            }
            
            // create angles for all pairs of bonded atoms
            for (size_t a = 0; a < bonded_atoms.size(); ++a) {
                for (size_t b = a + 1; b < bonded_atoms.size(); ++b) {
                    int i = bonded_atoms[a];
                    int k = bonded_atoms[b];
                    
                    // calculate current angle
                    double current_angle = system.angle_between_atoms(i, j, k);
                    
                    // get default angle parameters
                    double ka = get_default_angle_constant(system.atoms[i].Z, 
                                                           system.atoms[j].Z, 
                                                           system.atoms[k].Z);
                    
                    // add angle
                    system.add_angle(i, j, k, ka, current_angle);
                }
            }
        }
    }
    
    static void generate_bonds_from_template(System& system, const std::vector<std::pair<int, int>>& template_bonds,
                                           const BondParameters& params = BondParameters()) {
        system.bonds.clear();
        
        for (const auto& bond_pair : template_bonds) {
            int i = bond_pair.first;
            int j = bond_pair.second;
            
            if (i >= 0 && i < system.num_atoms() && j >= 0 && j < system.num_atoms()) {
                double distance = system.distance_between_atoms(i, j);
                auto bond_params = params.get_default_params(system.atoms[i].Z, system.atoms[j].Z);
                
                system.add_bond(i, j, bond_params.k, distance);
            }
        }
        
        if (params.auto_angles) {
            detect_angles(system, params);
        }
    }

private:
    static double get_default_angle_constant(int z_i, int z_j, int z_k) {
        // default angle spring constants based on central atom type
        switch (z_j) {
            case 1:  // hydrogen (rarely central)
                return 35.0;
            case 6:  // carbon
                return 40.0;  // sp3 carbon
            case 7:  // nitrogen
                return 50.0;
            case 8:  // oxygen
                return 55.0;
            case 15: // phosphorus
                return 45.0;
            case 16: // sulfur
                return 45.0;
            default:
                return 40.0;
        }
    }
};

}