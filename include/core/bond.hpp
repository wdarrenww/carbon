#pragma once

#include <string>

namespace carbon {

enum class BondType {
    SINGLE = 1,
    DOUBLE = 2,
    TRIPLE = 3,
    AROMATIC = 4,
    AMIDE = 5,
    IONIC = 6
};

struct Bond {
    int i, j;            // atom indices in the system
    double k;            // spring constant (kcal/mol/angstrom^2)
    double r0;           // equilibrium bond length (angstroms)
    BondType type;       // bond type
    double order;        // bond order (1.0, 2.0, 3.0, 1.5 for aromatic, etc.)
    bool is_rigid;       // constrained bond (e.g., for shake)
    
    Bond() : i(0), j(0), k(0.0), r0(0.0), type(BondType::SINGLE), order(1.0), is_rigid(false) {}
    
    Bond(int i_, int j_, double k_, double r0_, BondType type_ = BondType::SINGLE, double order_ = 1.0) 
        : i(i_), j(j_), k(k_), r0(r0_), type(type_), order(order_), is_rigid(false) {}
    
    bool involves_atom(int atom_idx) const {
        return i == atom_idx || j == atom_idx;
    }
    
    int other_atom(int atom_idx) const {
        if (i == atom_idx) return j;
        if (j == atom_idx) return i;
        return -1; // atom not in this bond
    }
    
    bool is_single() const { return type == BondType::SINGLE && order <= 1.1; }
    bool is_double() const { return type == BondType::DOUBLE || (order >= 1.8 && order <= 2.2); }
    bool is_triple() const { return type == BondType::TRIPLE || order >= 2.8; }
    bool is_aromatic() const { return type == BondType::AROMATIC || (order >= 1.3 && order <= 1.7); }
    
    std::string type_string() const {
        switch (type) {
            case BondType::SINGLE: return "single";
            case BondType::DOUBLE: return "double"; 
            case BondType::TRIPLE: return "triple";
            case BondType::AROMATIC: return "aromatic";
            case BondType::AMIDE: return "amide";
            case BondType::IONIC: return "ionic";
            default: return "unknown";
        }
    }
    
    // adjust spring constant based on bond order
    double effective_spring_constant() const {
        if (is_rigid) return 1e6;  // very stiff for constraints
        
        double base_k = k;
        
        // stronger bonds have higher spring constants
        if (is_triple()) {
            base_k *= 2.5;
        } else if (is_double()) {
            base_k *= 1.8;
        } else if (is_aromatic()) {
            base_k *= 1.4;
        }
        
        return base_k;
    }
    
    // adjust equilibrium length based on bond order
    double effective_equilibrium_length() const {
        double base_r0 = r0;
        
        // multiple bonds are shorter
        if (is_triple()) {
            base_r0 *= 0.87;
        } else if (is_double()) {
            base_r0 *= 0.93;
        } else if (is_aromatic()) {
            base_r0 *= 0.96;
        }
        
        return base_r0;
    }
};

}