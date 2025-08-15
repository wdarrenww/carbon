#pragma once

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

// core data structures
#include "core/vec3.hpp"
#include "core/atom.hpp"
#include "core/bond.hpp"
#include "core/angle.hpp"
#include "core/system.hpp"
#include "core/elements.hpp"
#include "core/bond_builder.hpp"
#include "core/fragments.hpp"
#include "core/charges.hpp"
#include "core/polymers.hpp"
#include "core/lattice_builder.hpp"

// physics engine
#include "core/forces.hpp"
#include "core/integrator.hpp"
#include "core/simulation.hpp"

// input/output
#include "io/xyz_io.hpp"
#include "io/builder.hpp"

namespace carbon {

// version information
constexpr int MAJOR_VERSION = 1;
constexpr int MINOR_VERSION = 0;
constexpr int PATCH_VERSION = 0;

inline std::string get_version_string() {
    return std::to_string(MAJOR_VERSION) + "." + 
           std::to_string(MINOR_VERSION) + "." + 
           std::to_string(PATCH_VERSION);
}

inline void print_banner() {
    std::cout << "carbon molecular simulation engine v" << get_version_string() << "\n";
    std::cout << "a c++ simulation engine for molecular dynamics\n";
    std::cout << "build date: " << __DATE__ << " " << __TIME__ << "\n\n";
}

// convenience functions for creating common molecules
namespace molecules {

inline System create_water_molecule() {
    System system;
    
    // create atoms
    int o_idx = system.add_atom(create_oxygen(Vec3(0.0, 0.0, 0.0)));
    int h1_idx = system.add_atom(create_hydrogen(Vec3(0.757, 0.586, 0.0)));
    int h2_idx = system.add_atom(create_hydrogen(Vec3(-0.757, 0.586, 0.0)));
    
    // add bonds (typical water parameters)
    system.add_bond(o_idx, h1_idx, 450.0, 0.957);  // o-h bond
    system.add_bond(o_idx, h2_idx, 450.0, 0.957);  // o-h bond
    
    return system;
}

inline System create_methane_molecule() {
    System system;
    
    // create carbon at center
    int c_idx = system.add_atom(create_carbon(Vec3(0.0, 0.0, 0.0)));
    
    // create hydrogens in tetrahedral geometry
    double bond_length = 1.09;  // c-h bond length
    double angle = 109.5 * M_PI / 180.0;  // tetrahedral angle
    
    int h1_idx = system.add_atom(create_hydrogen(Vec3(bond_length, 0.0, 0.0)));
    int h2_idx = system.add_atom(create_hydrogen(Vec3(-bond_length/3.0, bond_length*2.0/3.0*std::sqrt(2), 0.0)));
    int h3_idx = system.add_atom(create_hydrogen(Vec3(-bond_length/3.0, -bond_length/3.0*std::sqrt(2), bond_length*std::sqrt(6)/3.0)));
    int h4_idx = system.add_atom(create_hydrogen(Vec3(-bond_length/3.0, -bond_length/3.0*std::sqrt(2), -bond_length*std::sqrt(6)/3.0)));
    
    // add c-h bonds
    double k_ch = 340.0;  // typical c-h spring constant
    system.add_bond(c_idx, h1_idx, k_ch, bond_length);
    system.add_bond(c_idx, h2_idx, k_ch, bond_length);
    system.add_bond(c_idx, h3_idx, k_ch, bond_length);
    system.add_bond(c_idx, h4_idx, k_ch, bond_length);
    
    return system;
}

inline System create_nacl_dimer() {
    System system;
    
    // create ions
    Atom na = create_sodium(Vec3(-1.4, 0.0, 0.0));
    Atom cl = create_chlorine(Vec3(1.4, 0.0, 0.0));
    
    // set charges
    na.charge = 1.0;
    cl.charge = -1.0;
    
    system.add_atom(na);
    system.add_atom(cl);
    
    return system;
}

inline System create_ethanol() {
    System system;
    
    // create carbon backbone: CH3-CH2-OH
    int c1 = system.add_atom(create_carbon(Vec3(-1.269, -0.378, 0.000)));  // CH3
    int c2 = system.add_atom(create_carbon(Vec3(0.000, 0.518, 0.000)));    // CH2
    int o = system.add_atom(create_oxygen(Vec3(1.269, -0.378, 0.000)));     // OH
    
    // add hydrogens
    int h1 = system.add_atom(create_hydrogen(Vec3(-2.157, 0.256, 0.000)));
    int h2 = system.add_atom(create_hydrogen(Vec3(-1.269, -1.009, 0.890)));
    int h3 = system.add_atom(create_hydrogen(Vec3(-1.269, -1.009, -0.890)));
    int h4 = system.add_atom(create_hydrogen(Vec3(0.000, 1.149, 0.890)));
    int h5 = system.add_atom(create_hydrogen(Vec3(0.000, 1.149, -0.890)));
    int h6 = system.add_atom(create_hydrogen(Vec3(2.157, 0.256, 0.000)));
    
    // auto-detect bonds and angles
    BondBuilder::detect_bonds(system);
    
    return system;
}

inline System create_glycine() {
    System system;
    
    // NH2-CH2-COOH (glycine amino acid)
    int n = system.add_atom(create_nitrogen(Vec3(-1.896, 0.326, 0.000)));    // NH2
    int ca = system.add_atom(create_carbon(Vec3(-0.495, -0.174, 0.000)));    // alpha carbon
    int c = system.add_atom(create_carbon(Vec3(0.495, 1.025, 0.000)));      // carbonyl carbon
    int o1 = system.add_atom(create_oxygen(Vec3(0.000, 2.155, 0.000)));     // carbonyl oxygen
    int o2 = system.add_atom(create_oxygen(Vec3(1.726, 0.826, 0.000)));     // hydroxyl oxygen
    
    // add hydrogens
    int h1 = system.add_atom(create_hydrogen(Vec3(-2.464, -0.505, 0.000)));  // NH2
    int h2 = system.add_atom(create_hydrogen(Vec3(-2.175, 0.845, 0.816)));   // NH2
    int h3 = system.add_atom(create_hydrogen(Vec3(-0.495, -0.805, 0.890)));  // CH2
    int h4 = system.add_atom(create_hydrogen(Vec3(-0.495, -0.805, -0.890))); // CH2
    int h5 = system.add_atom(create_hydrogen(Vec3(2.265, 1.659, 0.000)));    // COOH
    
    // set partial charges for amino acid
    system.atoms[n].charge = -0.3;   // nitrogen
    system.atoms[c].charge = 0.5;    // carbonyl carbon  
    system.atoms[o1].charge = -0.5;  // carbonyl oxygen
    system.atoms[o2].charge = -0.6;  // hydroxyl oxygen
    system.atoms[h1].charge = 0.3;   // NH
    system.atoms[h2].charge = 0.3;   // NH
    system.atoms[h5].charge = 0.4;   // OH
    
    // auto-detect bonds and angles
    BondBuilder::detect_bonds(system);
    
    return system;
}

inline System create_benzene() {
    System system;
    
    // create benzene ring - hexagon with side length ~1.39 angstrom
    double r = 1.39;  // c-c bond length in benzene
    double angle_step = 2.0 * M_PI / 6.0;  // 60 degrees
    
    std::vector<int> carbons;
    for (int i = 0; i < 6; ++i) {
        double angle = i * angle_step;
        double x = r * std::cos(angle);
        double y = r * std::sin(angle);
        int c = system.add_atom(create_carbon(Vec3(x, y, 0.0)));
        carbons.push_back(c);
    }
    
    // add hydrogens
    for (int i = 0; i < 6; ++i) {
        double angle = i * angle_step;
        double x = (r + 1.08) * std::cos(angle);  // c-h bond length
        double y = (r + 1.08) * std::sin(angle);
        system.add_atom(create_hydrogen(Vec3(x, y, 0.0)));
    }
    
    // create aromatic bonds manually
    for (int i = 0; i < 6; ++i) {
        int j = (i + 1) % 6;
        system.add_bond(carbons[i], carbons[j], 350.0, 1.39, BondType::AROMATIC, 1.5);
    }
    
    // add c-h bonds
    for (int i = 0; i < 6; ++i) {
        system.add_bond(carbons[i], 6 + i, 340.0, 1.08);  // c-h bonds
    }
    
    // add angles
    BondBuilder::detect_angles(system);
    
    return system;
}

inline System create_acetone() {
    System system;
    
    // CH3-CO-CH3
    int c1 = system.add_atom(create_carbon(Vec3(-1.269, 0.000, 0.000)));   // left CH3
    int c2 = system.add_atom(create_carbon(Vec3(0.000, 0.000, 0.000)));    // carbonyl C
    int c3 = system.add_atom(create_carbon(Vec3(1.269, 0.000, 0.000)));    // right CH3
    int o = system.add_atom(create_oxygen(Vec3(0.000, 1.229, 0.000)));      // carbonyl O
    
    // add hydrogens to methyl groups
    system.add_atom(create_hydrogen(Vec3(-1.269, -0.693, 0.890)));
    system.add_atom(create_hydrogen(Vec3(-1.269, -0.693, -0.890)));
    system.add_atom(create_hydrogen(Vec3(-2.157, 0.634, 0.000)));
    system.add_atom(create_hydrogen(Vec3(1.269, -0.693, 0.890)));
    system.add_atom(create_hydrogen(Vec3(1.269, -0.693, -0.890)));
    system.add_atom(create_hydrogen(Vec3(2.157, 0.634, 0.000)));
    
    // set partial charges
    system.atoms[c2].charge = 0.5;   // carbonyl carbon
    system.atoms[o].charge = -0.5;   // carbonyl oxygen
    
    // create bonds manually with proper types
    system.add_bond(c1, c2, 320.0, 1.50);                    // c-c single
    system.add_bond(c2, c3, 320.0, 1.50);                    // c-c single  
    system.add_bond(c2, o, 570.0, 1.23, BondType::DOUBLE);   // c=o double
    
    // add c-h bonds
    for (int i = 4; i < 10; ++i) {
        int c_idx = (i < 7) ? c1 : c3;
        system.add_bond(c_idx, i, 340.0, 1.09);
    }
    
    BondBuilder::detect_angles(system);
    
    return system;
}

inline System create_phosphate_ion() {
    System system;
    
    // PO4^3- tetrahedral
    int p = system.add_atom(create_phosphorus(Vec3(0.0, 0.0, 0.0)));
    
    // tetrahedral geometry for oxygens
    double r = 1.54;  // p-o bond length
    system.add_atom(create_oxygen(Vec3(r, 0.0, 0.0)));
    system.add_atom(create_oxygen(Vec3(-r/3.0, r*2.0/3.0*std::sqrt(2), 0.0)));
    system.add_atom(create_oxygen(Vec3(-r/3.0, -r/3.0*std::sqrt(2), r*std::sqrt(6)/3.0)));
    system.add_atom(create_oxygen(Vec3(-r/3.0, -r/3.0*std::sqrt(2), -r*std::sqrt(6)/3.0)));
    
    // set charges
    system.atoms[p].charge = 1.2;     // phosphorus
    for (int i = 1; i < 5; ++i) {
        system.atoms[i].charge = -1.05;  // oxygens
    }
    
    // create p-o bonds  
    for (int i = 1; i < 5; ++i) {
        system.add_bond(p, i, 400.0, 1.54);
    }
    
    BondBuilder::detect_angles(system);
    
    return system;
}

} // namespace molecules

} // namespace carbon