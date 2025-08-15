#pragma once

#include "system.hpp"
#include "elements.hpp"
#include "bond_builder.hpp"
#include <vector>

namespace carbon {

class LatticeBuilder {
public:
    enum class LatticeType {
        SIMPLE_CUBIC,
        BODY_CENTERED_CUBIC,
        FACE_CENTERED_CUBIC,
        HEXAGONAL_CLOSE_PACKED,
        DIAMOND,
        SODIUM_CHLORIDE,
        CESIUM_CHLORIDE,
        FLUORITE
    };
    
    struct LatticeParameters {
        double a = 5.0;         // lattice parameter a (angstroms)
        double b = 0.0;         // lattice parameter b (if 0, uses a)
        double c = 0.0;         // lattice parameter c (if 0, uses a)
        int nx = 3;             // number of unit cells in x
        int ny = 3;             // number of unit cells in y
        int nz = 3;             // number of unit cells in z
        bool periodic = true;   // set up periodic boundary conditions
        bool auto_bonds = true; // automatically detect bonds
    };
    
    static System build_simple_cubic(int element_z, const LatticeParameters& params = LatticeParameters()) {
        System system;
        
        double a = params.a;
        double b = (params.b > 0) ? params.b : a;
        double c = (params.c > 0) ? params.c : a;
        
        // create atoms at lattice points
        for (int i = 0; i < params.nx; ++i) {
            for (int j = 0; j < params.ny; ++j) {
                for (int k = 0; k < params.nz; ++k) {
                    double x = i * a;
                    double y = j * b;
                    double z = k * c;
                    
                    Atom atom = get_element_table().create_atom(element_z, Vec3(x, y, z));
                    system.add_atom(atom);
                }
            }
        }
        
        if (params.periodic) {
            system.periodic = true;
            system.box_size = Vec3(params.nx * a, params.ny * b, params.nz * c);
        }
        
        if (params.auto_bonds) {
            BondBuilder::detect_bonds(system);
        }
        
        return system;
    }
    
    static System build_fcc(int element_z, const LatticeParameters& params = LatticeParameters()) {
        System system;
        double a = params.a;
        
        // fcc has atoms at: (0,0,0), (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5)
        std::vector<Vec3> basis = {
            Vec3(0.0, 0.0, 0.0),
            Vec3(0.5, 0.5, 0.0),
            Vec3(0.5, 0.0, 0.5),
            Vec3(0.0, 0.5, 0.5)
        };
        
        for (int i = 0; i < params.nx; ++i) {
            for (int j = 0; j < params.ny; ++j) {
                for (int k = 0; k < params.nz; ++k) {
                    for (const auto& offset : basis) {
                        double x = (i + offset.x) * a;
                        double y = (j + offset.y) * a;
                        double z = (k + offset.z) * a;
                        
                        Atom atom = get_element_table().create_atom(element_z, Vec3(x, y, z));
                        system.add_atom(atom);
                    }
                }
            }
        }
        
        if (params.periodic) {
            system.periodic = true;
            system.box_size = Vec3(params.nx * a, params.ny * a, params.nz * a);
        }
        
        if (params.auto_bonds) {
            BondBuilder::detect_bonds(system);
        }
        
        return system;
    }
    
    static System build_diamond(const LatticeParameters& params = LatticeParameters()) {
        System system;
        double a = params.a;
        
        // diamond structure: two interpenetrating fcc lattices
        std::vector<Vec3> basis = {
            Vec3(0.0, 0.0, 0.0),
            Vec3(0.5, 0.5, 0.0),
            Vec3(0.5, 0.0, 0.5),
            Vec3(0.0, 0.5, 0.5),
            Vec3(0.25, 0.25, 0.25),
            Vec3(0.75, 0.75, 0.25),
            Vec3(0.75, 0.25, 0.75),
            Vec3(0.25, 0.75, 0.75)
        };
        
        for (int i = 0; i < params.nx; ++i) {
            for (int j = 0; j < params.ny; ++j) {
                for (int k = 0; k < params.nz; ++k) {
                    for (const auto& offset : basis) {
                        double x = (i + offset.x) * a;
                        double y = (j + offset.y) * a;
                        double z = (k + offset.z) * a;
                        
                        Atom atom = create_carbon(Vec3(x, y, z));
                        system.add_atom(atom);
                    }
                }
            }
        }
        
        if (params.periodic) {
            system.periodic = true;
            system.box_size = Vec3(params.nx * a, params.ny * a, params.nz * a);
        }
        
        if (params.auto_bonds) {
            BondBuilder::BondParameters bond_params;
            bond_params.bond_tolerance = 1.1;  // tighter tolerance for diamond
            BondBuilder::detect_bonds(system, bond_params);
        }
        
        return system;
    }
    
    static System build_nacl(const LatticeParameters& params = LatticeParameters()) {
        System system;
        double a = params.a;
        
        // nacl structure: na+ at (0,0,0) type positions, cl- at (0.5,0,0) type
        for (int i = 0; i < params.nx; ++i) {
            for (int j = 0; j < params.ny; ++j) {
                for (int k = 0; k < params.nz; ++k) {
                    // sodium positions
                    {
                        double x = i * a;
                        double y = j * a;
                        double z = k * a;
                        
                        Atom na = create_sodium(Vec3(x, y, z));
                        na.charge = 1.0;
                        system.add_atom(na);
                    }
                    
                    // chlorine positions  
                    {
                        double x = (i + 0.5) * a;
                        double y = j * a;
                        double z = k * a;
                        
                        Atom cl = create_chlorine(Vec3(x, y, z));
                        cl.charge = -1.0;
                        system.add_atom(cl);
                    }
                    
                    // additional cl positions for complete structure
                    {
                        double x = i * a;
                        double y = (j + 0.5) * a;
                        double z = k * a;
                        
                        Atom cl = create_chlorine(Vec3(x, y, z));
                        cl.charge = -1.0;
                        system.add_atom(cl);
                    }
                    
                    {
                        double x = i * a;
                        double y = j * a;
                        double z = (k + 0.5) * a;
                        
                        Atom cl = create_chlorine(Vec3(x, y, z));
                        cl.charge = -1.0;
                        system.add_atom(cl);
                    }
                }
            }
        }
        
        if (params.periodic) {
            system.periodic = true;
            system.box_size = Vec3(params.nx * a, params.ny * a, params.nz * a);
        }
        
        return system;  // no bonds for ionic crystals
    }
    
    static System build_water_ice(const LatticeParameters& params = LatticeParameters()) {
        System system;
        double a = params.a;
        
        // simple hexagonal ice structure (ice ih approximation)
        for (int i = 0; i < params.nx; ++i) {
            for (int j = 0; j < params.ny; ++j) {
                for (int k = 0; k < params.nz; ++k) {
                    double x = i * a;
                    double y = j * a;
                    double z = k * a;
                    
                    // create water molecule
                    System water = molecules::create_water_molecule();
                    
                    // translate and add to lattice
                    for (auto& atom : water.atoms) {
                        atom.r += Vec3(x, y, z);
                        system.add_atom(atom);
                    }
                    
                    // add bonds
                    int base_idx = system.num_atoms() - 3;  // last 3 atoms added
                    if (base_idx >= 0) {
                        system.add_bond(base_idx, base_idx + 1, 450.0, 0.957);  // o-h
                        system.add_bond(base_idx, base_idx + 2, 450.0, 0.957);  // o-h
                    }
                }
            }
        }
        
        if (params.periodic) {
            system.periodic = true;
            system.box_size = Vec3(params.nx * a, params.ny * a, params.nz * a);
        }
        
        return system;
    }
    
    static System build_graphite_layer(const LatticeParameters& params = LatticeParameters()) {
        System system;
        double a = params.a;  // in-plane lattice constant
        double bond_length = 1.42;  // c-c bond length in graphene
        
        // hexagonal lattice for graphene layer
        for (int i = 0; i < params.nx; ++i) {
            for (int j = 0; j < params.ny; ++j) {
                // two atoms per unit cell in graphene
                double x1 = i * a;
                double y1 = j * a * std::sqrt(3.0);
                
                double x2 = x1 + a * 0.5;
                double y2 = y1 + a * std::sqrt(3.0) / 6.0;
                
                system.add_atom(create_carbon(Vec3(x1, y1, 0.0)));
                system.add_atom(create_carbon(Vec3(x2, y2, 0.0)));
            }
        }
        
        if (params.periodic) {
            system.periodic = true;
            system.box_size = Vec3(params.nx * a, params.ny * a * std::sqrt(3.0), 10.0);
        }
        
        if (params.auto_bonds) {
            BondBuilder::BondParameters bond_params;
            bond_params.bond_tolerance = 1.1;
            bond_params.max_bond_length = 1.6;
            BondBuilder::detect_bonds(system, bond_params);
        }
        
        return system;
    }
};

}