#pragma once

#include "system.hpp"
#include "fragments.hpp"
#include "elements.hpp"
#include "bond_builder.hpp"

namespace carbon {

class PolymerBuilder {
public:
    enum class PolymerType {
        LINEAR,
        BRANCHED,
        CYCLIC,
        NETWORK
    };
    
    struct PolymerParameters {
        int chain_length = 10;           // number of repeat units
        double bond_length = 1.54;       // backbone bond length
        double bond_angle = 109.5;       // backbone bond angle (degrees)
        bool add_hydrogens = true;       // add terminal hydrogens
        bool auto_detect_bonds = true;   // automatically detect bonds
        Vec3 chain_direction = Vec3(1.0, 0.0, 0.0);  // initial chain direction
        double randomness = 0.0;         // conformational randomness (0-1)
    };
    
    // build linear alkane chain
    static System build_alkane_chain(int n_carbons, const PolymerParameters& params = PolymerParameters()) {
        System system;
        
        if (n_carbons < 1) return system;
        
        std::vector<int> carbon_indices;
        double bond_length = params.bond_length;
        double angle_rad = params.bond_angle * M_PI / 180.0;
        
        // add first carbon
        carbon_indices.push_back(system.add_atom(create_carbon(Vec3(0.0, 0.0, 0.0))));
        
        Vec3 direction = params.chain_direction.normalized();
        Vec3 current_pos = Vec3(0.0, 0.0, 0.0);
        
        // add subsequent carbons
        for (int i = 1; i < n_carbons; ++i) {
            current_pos += direction * bond_length;
            
            // add some randomness if requested
            if (params.randomness > 0.0) {
                Vec3 random_offset = Vec3(
                    (std::rand() / double(RAND_MAX) - 0.5) * params.randomness,
                    (std::rand() / double(RAND_MAX) - 0.5) * params.randomness,
                    (std::rand() / double(RAND_MAX) - 0.5) * params.randomness
                );
                current_pos += random_offset;
            }
            
            carbon_indices.push_back(system.add_atom(create_carbon(current_pos)));
            
            // update direction for next carbon (slight rotation around z-axis)
            double rotation = angle_rad - M_PI;  // supplement to tetrahedral angle
            direction = Vec3(
                direction.x * std::cos(rotation) - direction.y * std::sin(rotation),
                direction.x * std::sin(rotation) + direction.y * std::cos(rotation),
                direction.z
            );
        }
        
        // add c-c bonds
        for (int i = 0; i < n_carbons - 1; ++i) {
            system.add_bond(carbon_indices[i], carbon_indices[i+1], 320.0, bond_length);
        }
        
        // add hydrogens if requested
        if (params.add_hydrogens) {
            add_alkane_hydrogens(system, carbon_indices);
        }
        
        if (params.auto_detect_bonds && params.add_hydrogens) {
            BondBuilder::detect_angles(system);
        }
        
        return system;
    }
    
    // build polyethylene chain
    static System build_polyethylene(int n_units, const PolymerParameters& params = PolymerParameters()) {
        System system;
        
        std::vector<int> carbon_indices;
        double bond_length = params.bond_length;
        
        Vec3 current_pos = Vec3(0.0, 0.0, 0.0);
        Vec3 direction = params.chain_direction.normalized();
        
        for (int i = 0; i < n_units; ++i) {
            // add two carbons for each ethylene unit: -CH2-CH2-
            int c1 = system.add_atom(create_carbon(current_pos));
            current_pos += direction * bond_length;
            int c2 = system.add_atom(create_carbon(current_pos));
            current_pos += direction * bond_length;
            
            carbon_indices.push_back(c1);
            carbon_indices.push_back(c2);
            
            // add c-c bond within unit
            system.add_bond(c1, c2, 320.0, bond_length);
            
            // add bond to previous unit
            if (i > 0) {
                system.add_bond(carbon_indices[2*i-1], c1, 320.0, bond_length);
            }
        }
        
        // add hydrogens
        if (params.add_hydrogens) {
            add_alkane_hydrogens(system, carbon_indices);
        }
        
        if (params.auto_detect_bonds && params.add_hydrogens) {
            BondBuilder::detect_angles(system);
        }
        
        return system;
    }
    
    // build polystyrene chain
    static System build_polystyrene(int n_units, const PolymerParameters& params = PolymerParameters()) {
        System system;
        
        std::vector<int> backbone_carbons;
        double bond_length = params.bond_length;
        
        Vec3 current_pos = Vec3(0.0, 0.0, 0.0);
        Vec3 direction = params.chain_direction.normalized();
        Vec3 side_direction = Vec3(0.0, 1.0, 0.0);  // perpendicular for side chains
        
        for (int i = 0; i < n_units; ++i) {
            // backbone: -CH2-CH(phenyl)-
            int c1 = system.add_atom(create_carbon(current_pos));
            current_pos += direction * bond_length;
            int c2 = system.add_atom(create_carbon(current_pos));  // carbon with phenyl
            current_pos += direction * bond_length;
            
            backbone_carbons.push_back(c1);
            backbone_carbons.push_back(c2);
            
            // add backbone bond
            system.add_bond(c1, c2, 320.0, bond_length);
            if (i > 0) {
                system.add_bond(backbone_carbons[2*i-1], c1, 320.0, bond_length);
            }
            
            // add phenyl ring as side chain
            Vec3 phenyl_pos = system.atoms[c2].r + side_direction * 1.54;
            int phenyl_start = FragmentLibrary::add_fragment_to_system(system, "phenyl", phenyl_pos);
            
            // connect phenyl to backbone
            if (phenyl_start >= 0) {
                system.add_bond(c2, phenyl_start, 320.0, 1.54);
            }
        }
        
        // add backbone hydrogens
        if (params.add_hydrogens) {
            add_backbone_hydrogens(system, backbone_carbons);
        }
        
        return system;
    }
    
    // build protein backbone
    static System build_protein_backbone(const std::vector<std::string>& amino_acids, 
                                       const PolymerParameters& params = PolymerParameters()) {
        System system;
        
        Vec3 current_pos = Vec3(0.0, 0.0, 0.0);
        Vec3 direction = params.chain_direction.normalized();
        double peptide_length = 1.33;  // peptide bond length
        double ca_length = 1.54;       // ca-c/n bond length
        
        for (size_t i = 0; i < amino_acids.size(); ++i) {
            // add n-ca-c backbone for each amino acid
            int n = system.add_atom(create_nitrogen(current_pos));
            current_pos += direction * ca_length;
            
            int ca = system.add_atom(create_carbon(current_pos));
            current_pos += direction * ca_length;
            
            int c = system.add_atom(create_carbon(current_pos));
            current_pos += direction * peptide_length;
            
            // add backbone bonds
            system.add_bond(n, ca, 350.0, ca_length);
            system.add_bond(ca, c, 320.0, ca_length);
            
            // add peptide bond to previous residue
            if (i > 0) {
                int prev_c = system.num_atoms() - 6;  // c from previous residue
                system.add_bond(prev_c, n, 400.0, peptide_length, BondType::AMIDE);
            }
            
            // add carbonyl oxygen
            Vec3 o_pos = system.atoms[c].r + Vec3(0.0, 1.23, 0.0);
            int o = system.add_atom(create_oxygen(o_pos));
            system.add_bond(c, o, 570.0, 1.23, BondType::DOUBLE);
            
            // could add side chain based on amino acid type
            // (simplified - just add basic side chains)
            if (amino_acids[i] == "ALA") {
                // alanine: add methyl side chain
                Vec3 cb_pos = system.atoms[ca].r + Vec3(0.0, -1.54, 0.0);
                int cb = system.add_atom(create_carbon(cb_pos));
                system.add_bond(ca, cb, 320.0, 1.54);
                
                if (params.add_hydrogens) {
                    // add methyl hydrogens
                    system.add_atom(create_hydrogen(cb_pos + Vec3(1.09, 0.0, 0.0)));
                    system.add_atom(create_hydrogen(cb_pos + Vec3(-0.545, 0.943, 0.0)));
                    system.add_atom(create_hydrogen(cb_pos + Vec3(-0.545, -0.943, 0.0)));
                }
            }
        }
        
        if (params.auto_detect_bonds) {
            BondBuilder::detect_angles(system);
        }
        
        return system;
    }
    
    // build cyclic polymer (e.g., cyclohexane)
    static System build_cyclic_alkane(int n_carbons, const PolymerParameters& params = PolymerParameters()) {
        System system;
        
        if (n_carbons < 3) return system;
        
        std::vector<int> carbon_indices;
        double radius = params.bond_length / (2.0 * std::sin(M_PI / n_carbons));
        
        // create ring in xy-plane
        for (int i = 0; i < n_carbons; ++i) {
            double angle = 2.0 * M_PI * i / n_carbons;
            double x = radius * std::cos(angle);
            double y = radius * std::sin(angle);
            
            carbon_indices.push_back(system.add_atom(create_carbon(Vec3(x, y, 0.0))));
        }
        
        // add ring bonds
        for (int i = 0; i < n_carbons; ++i) {
            int j = (i + 1) % n_carbons;
            double distance = system.distance_between_atoms(carbon_indices[i], carbon_indices[j]);
            system.add_bond(carbon_indices[i], carbon_indices[j], 320.0, distance);
        }
        
        // add hydrogens
        if (params.add_hydrogens) {
            add_cyclic_hydrogens(system, carbon_indices, radius);
        }
        
        if (params.auto_detect_bonds && params.add_hydrogens) {
            BondBuilder::detect_angles(system);
        }
        
        return system;
    }

private:
    static void add_alkane_hydrogens(System& system, const std::vector<int>& carbon_indices) {
        for (size_t i = 0; i < carbon_indices.size(); ++i) {
            const Vec3& c_pos = system.atoms[carbon_indices[i]].r;
            
            // determine how many hydrogens to add based on bonding
            int num_bonds = 0;
            for (const auto& bond : system.bonds) {
                if (bond.involves_atom(carbon_indices[i])) {
                    num_bonds++;
                }
            }
            
            int num_hydrogens = 4 - num_bonds;  // tetrahedral carbon
            
            // add hydrogens in tetrahedral positions
            for (int h = 0; h < num_hydrogens; ++h) {
                double angle = 2.0 * M_PI * h / num_hydrogens;
                Vec3 h_pos = c_pos + Vec3(
                    1.09 * std::cos(angle),
                    1.09 * std::sin(angle),
                    1.09 * ((h % 2 == 0) ? 1.0 : -1.0)
                );
                
                int h_idx = system.add_atom(create_hydrogen(h_pos));
                system.add_bond(carbon_indices[i], h_idx, 340.0, 1.09);
            }
        }
    }
    
    static void add_backbone_hydrogens(System& system, const std::vector<int>& carbon_indices) {
        // simplified hydrogen addition for polymer backbones
        for (int c_idx : carbon_indices) {
            const Vec3& c_pos = system.atoms[c_idx].r;
            
            // add one or two hydrogens per carbon
            Vec3 h1_pos = c_pos + Vec3(0.0, 1.09, 0.0);
            Vec3 h2_pos = c_pos + Vec3(0.0, -1.09, 0.0);
            
            int h1 = system.add_atom(create_hydrogen(h1_pos));
            int h2 = system.add_atom(create_hydrogen(h2_pos));
            
            system.add_bond(c_idx, h1, 340.0, 1.09);
            system.add_bond(c_idx, h2, 340.0, 1.09);
        }
    }
    
    static void add_cyclic_hydrogens(System& system, const std::vector<int>& carbon_indices, double radius) {
        for (int c_idx : carbon_indices) {
            const Vec3& c_pos = system.atoms[c_idx].r;
            
            // add hydrogens pointing inward and outward from ring
            Vec3 radial = c_pos.normalized();
            Vec3 h1_pos = c_pos + radial * 1.09;        // outward
            Vec3 h2_pos = c_pos - radial * 1.09;        // inward
            
            int h1 = system.add_atom(create_hydrogen(h1_pos));
            int h2 = system.add_atom(create_hydrogen(h2_pos));
            
            system.add_bond(c_idx, h1, 340.0, 1.09);
            system.add_bond(c_idx, h2, 340.0, 1.09);
        }
    }
};

}