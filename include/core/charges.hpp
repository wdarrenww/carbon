#pragma once

#include "system.hpp"
#include "elements.hpp"
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace carbon {

class ChargeAssignment {
public:
    enum class Method {
        GASTEIGER,          // gasteiger partial charges
        SIMPLE_ELECTRONEGATIVITY, // based on electronegativity differences
        FORMAL_CHARGES,     // formal oxidation state charges
        ZERO,              // all charges set to zero
        CUSTOM             // user-defined charges
    };
    
    struct ChargeParameters {
        Method method = Method::SIMPLE_ELECTRONEGATIVITY;
        int max_iterations = 50;      // for iterative methods
        double convergence = 1e-6;    // convergence criterion
        double damping = 0.5;         // damping factor for iterative methods
        bool preserve_total_charge = true;  // maintain overall system charge
        double target_total_charge = 0.0;   // desired total charge
    };
    
    static void assign_charges(System& system, const ChargeParameters& params = ChargeParameters()) {
        switch (params.method) {
            case Method::GASTEIGER:
                assign_gasteiger_charges(system, params);
                break;
            case Method::SIMPLE_ELECTRONEGATIVITY:
                assign_electronegativity_charges(system, params);
                break;
            case Method::FORMAL_CHARGES:
                assign_formal_charges(system, params);
                break;
            case Method::ZERO:
                assign_zero_charges(system);
                break;
            default:
                assign_zero_charges(system);
        }
        
        if (params.preserve_total_charge) {
            normalize_total_charge(system, params.target_total_charge);
        }
    }
    
private:
    static void assign_zero_charges(System& system) {
        for (auto& atom : system.atoms) {
            atom.charge = 0.0;
        }
    }
    
    static void assign_formal_charges(System& system, const ChargeParameters& params) {
        // assign formal charges based on common oxidation states
        for (auto& atom : system.atoms) {
            switch (atom.Z) {
                case 1:  // hydrogen
                    atom.charge = 1.0;   // h+ in most compounds
                    break;
                case 6:  // carbon
                    atom.charge = 0.0;   // usually neutral in organics
                    break;
                case 7:  // nitrogen
                    atom.charge = 0.0;   // depends on bonding, default neutral
                    break;
                case 8:  // oxygen
                    atom.charge = 0.0;   // depends on bonding, default neutral
                    break;
                case 9:  // fluorine
                    atom.charge = -1.0;  // f- in ionic compounds
                    break;
                case 11: // sodium
                    atom.charge = 1.0;   // na+
                    break;
                case 12: // magnesium
                    atom.charge = 2.0;   // mg2+
                    break;
                case 15: // phosphorus
                    atom.charge = 0.0;   // depends on oxidation state
                    break;
                case 16: // sulfur
                    atom.charge = 0.0;   // depends on oxidation state
                    break;
                case 17: // chlorine
                    atom.charge = -1.0;  // cl-
                    break;
                case 19: // potassium
                    atom.charge = 1.0;   // k+
                    break;
                case 20: // calcium
                    atom.charge = 2.0;   // ca2+
                    break;
                case 35: // bromine
                    atom.charge = -1.0;  // br-
                    break;
                case 53: // iodine
                    atom.charge = -1.0;  // i-
                    break;
                default:
                    atom.charge = 0.0;
            }
        }
    }
    
    static void assign_electronegativity_charges(System& system, const ChargeParameters& params) {
        // simple electronegativity-based charge assignment
        std::unordered_map<int, double> electronegativity = {
            {1, 2.20},   // H
            {6, 2.55},   // C
            {7, 3.04},   // N
            {8, 3.44},   // O
            {9, 3.98},   // F
            {11, 0.93},  // Na
            {12, 1.31},  // Mg
            {15, 2.19},  // P
            {16, 2.58},  // S
            {17, 3.16},  // Cl
            {19, 0.82},  // K
            {20, 1.00},  // Ca
            {35, 2.96},  // Br
            {53, 2.66}   // I
        };
        
        // initialize charges to zero
        for (auto& atom : system.atoms) {
            atom.charge = 0.0;
        }
        
        // assign charges based on electronegativity differences in bonds
        for (const auto& bond : system.bonds) {
            if (bond.i >= system.num_atoms() || bond.j >= system.num_atoms()) continue;
            
            Atom& atom_i = system.atoms[bond.i];
            Atom& atom_j = system.atoms[bond.j];
            
            double chi_i = electronegativity.count(atom_i.Z) ? electronegativity[atom_i.Z] : 2.0;
            double chi_j = electronegativity.count(atom_j.Z) ? electronegativity[atom_j.Z] : 2.0;
            
            double delta_chi = chi_j - chi_i;
            double charge_transfer = 0.1 * delta_chi * bond.order;  // empirical scaling
            
            atom_i.charge += charge_transfer;
            atom_j.charge -= charge_transfer;
        }
        
        // scale charges to reasonable range
        for (auto& atom : system.atoms) {
            atom.charge = std::max(-2.0, std::min(2.0, atom.charge));
        }
    }
    
    static void assign_gasteiger_charges(System& system, const ChargeParameters& params) {
        // simplified gasteiger electronegativity equalization
        std::unordered_map<int, double> electronegativity = {
            {1, 2.20}, {6, 2.55}, {7, 3.04}, {8, 3.44}, {9, 3.98},
            {11, 0.93}, {12, 1.31}, {15, 2.19}, {16, 2.58}, {17, 3.16},
            {19, 0.82}, {20, 1.00}, {35, 2.96}, {53, 2.66}
        };
        
        std::unordered_map<int, double> hardness = {
            {1, 20.02}, {6, 10.39}, {7, 11.54}, {8, 13.36}, {9, 14.66},
            {11, 6.11}, {12, 3.15}, {15, 8.90}, {16, 10.14}, {17, 9.38},
            {19, 4.34}, {20, 2.37}, {35, 7.59}, {53, 6.76}
        };
        
        // initialize charges
        for (auto& atom : system.atoms) {
            atom.charge = 0.0;
        }
        
        // iterative charge equilibration
        for (int iter = 0; iter < params.max_iterations; ++iter) {
            std::vector<double> new_charges(system.num_atoms());
            bool converged = true;
            
            for (int i = 0; i < system.num_atoms(); ++i) {
                const Atom& atom = system.atoms[i];
                
                double chi_i = electronegativity.count(atom.Z) ? electronegativity[atom.Z] : 2.0;
                double eta_i = hardness.count(atom.Z) ? hardness[atom.Z] : 10.0;
                
                double sum_chi = 0.0;
                double sum_eta_inv = 1.0 / eta_i;  // self term
                int coordination = 1;  // self
                
                // sum over bonded neighbors
                for (const auto& bond : system.bonds) {
                    int j = -1;
                    if (bond.i == i) j = bond.j;
                    else if (bond.j == i) j = bond.i;
                    
                    if (j >= 0 && j < system.num_atoms()) {
                        const Atom& neighbor = system.atoms[j];
                        double chi_j = electronegativity.count(neighbor.Z) ? electronegativity[neighbor.Z] : 2.0;
                        double eta_j = hardness.count(neighbor.Z) ? hardness[neighbor.Z] : 10.0;
                        
                        sum_chi += chi_j;
                        sum_eta_inv += 1.0 / eta_j;
                        coordination++;
                    }
                }
                
                // calculate new charge
                double avg_chi = (chi_i + sum_chi) / coordination;
                double effective_eta = coordination / sum_eta_inv;
                
                new_charges[i] = params.damping * (avg_chi - chi_i) / effective_eta + 
                                (1.0 - params.damping) * atom.charge;
                
                if (std::abs(new_charges[i] - atom.charge) > params.convergence) {
                    converged = false;
                }
            }
            
            // update charges
            for (int i = 0; i < system.num_atoms(); ++i) {
                system.atoms[i].charge = new_charges[i];
            }
            
            if (converged) break;
        }
        
        // scale to reasonable range
        for (auto& atom : system.atoms) {
            atom.charge = std::max(-2.0, std::min(2.0, atom.charge));
        }
    }
    
    static void normalize_total_charge(System& system, double target_charge) {
        double total_charge = 0.0;
        for (const auto& atom : system.atoms) {
            total_charge += atom.charge;
        }
        
        double correction = (target_charge - total_charge) / system.num_atoms();
        
        for (auto& atom : system.atoms) {
            atom.charge += correction;
        }
    }

public:
    // utility functions for analyzing charges
    static double calculate_total_charge(const System& system) {
        double total = 0.0;
        for (const auto& atom : system.atoms) {
            total += atom.charge;
        }
        return total;
    }
    
    static double calculate_dipole_moment(const System& system) {
        Vec3 dipole;
        for (const auto& atom : system.atoms) {
            dipole += atom.r * atom.charge;
        }
        return dipole.magnitude();
    }
    
    static void print_charge_summary(const System& system) {
        double total_charge = calculate_total_charge(system);
        double dipole = calculate_dipole_moment(system);
        
        std::cout << "charge analysis:\n";
        std::cout << "  total charge: " << std::fixed << std::setprecision(3) << total_charge << " e\n";
        std::cout << "  dipole moment: " << std::setprecision(3) << dipole << " e*angstrom\n";
        std::cout << "  individual atomic charges:\n";
        
        const ElementTable& table = get_element_table();
        for (int i = 0; i < system.num_atoms(); ++i) {
            const Atom& atom = system.atoms[i];
            std::string symbol = table.get_symbol(atom.Z);
            std::cout << "    " << i << " " << symbol << ": " 
                      << std::setprecision(4) << atom.charge << " e\n";
        }
    }
    
    // specialized charge assignment for common functional groups
    static void assign_amino_acid_charges(System& system) {
        // simple amino acid charge assignment
        for (auto& atom : system.atoms) {
            switch (atom.Z) {
                case 1:  // hydrogen
                    atom.charge = 0.3;   // slightly positive
                    break;
                case 6:  // carbon
                    if (is_carbonyl_carbon(system, atom)) {
                        atom.charge = 0.5;   // carbonyl carbon
                    } else {
                        atom.charge = 0.0;   // aliphatic carbon
                    }
                    break;
                case 7:  // nitrogen
                    atom.charge = -0.3;  // amino nitrogen
                    break;
                case 8:  // oxygen
                    if (is_carbonyl_oxygen(system, atom)) {
                        atom.charge = -0.5;  // carbonyl oxygen
                    } else {
                        atom.charge = -0.6;  // hydroxyl oxygen
                    }
                    break;
                default:
                    atom.charge = 0.0;
            }
        }
        
        normalize_total_charge(system, 0.0);
    }
    
private:
    static bool is_carbonyl_carbon(const System& system, const Atom& atom) {
        // check if carbon is double-bonded to oxygen
        for (const auto& bond : system.bonds) {
            int other_idx = -1;
            if (bond.i < system.num_atoms() && system.atoms[bond.i].r == atom.r) {
                other_idx = bond.j;
            } else if (bond.j < system.num_atoms() && system.atoms[bond.j].r == atom.r) {
                other_idx = bond.i;
            }
            
            if (other_idx >= 0 && other_idx < system.num_atoms() &&
                system.atoms[other_idx].Z == 8 && bond.is_double()) {
                return true;
            }
        }
        return false;
    }
    
    static bool is_carbonyl_oxygen(const System& system, const Atom& atom) {
        // check if oxygen is double-bonded to carbon
        for (const auto& bond : system.bonds) {
            int other_idx = -1;
            if (bond.i < system.num_atoms() && system.atoms[bond.i].r == atom.r) {
                other_idx = bond.j;
            } else if (bond.j < system.num_atoms() && system.atoms[bond.j].r == atom.r) {
                other_idx = bond.i;
            }
            
            if (other_idx >= 0 && other_idx < system.num_atoms() &&
                system.atoms[other_idx].Z == 6 && bond.is_double()) {
                return true;
            }
        }
        return false;
    }
};

}