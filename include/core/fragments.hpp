#pragma once

#include "system.hpp"
#include "elements.hpp"
#include "bond_builder.hpp"
#include <unordered_map>
#include <string>

namespace carbon {

struct Fragment {
    std::vector<Atom> atoms;
    std::vector<std::pair<int, int>> bonds;  // atom index pairs
    std::vector<std::tuple<int, int, int>> angles;  // i-j-k angle triplets
    std::vector<int> connection_points;  // atoms that can connect to other fragments
    std::string name;
    std::string description;
    
    Fragment(const std::string& name_ = "", const std::string& desc_ = "")
        : name(name_), description(desc_) {}
    
    int num_atoms() const { return static_cast<int>(atoms.size()); }
    int num_bonds() const { return static_cast<int>(bonds.size()); }
    int num_connection_points() const { return static_cast<int>(connection_points.size()); }
};

class FragmentLibrary {
private:
    std::unordered_map<std::string, Fragment> fragments;
    
    void initialize_common_fragments() {
        // methyl group -CH3
        {
            Fragment methyl("methyl", "methyl group -CH3");
            methyl.atoms.push_back(create_carbon(Vec3(0.0, 0.0, 0.0)));
            methyl.atoms.push_back(create_hydrogen(Vec3(0.0, 1.09, 0.0)));
            methyl.atoms.push_back(create_hydrogen(Vec3(0.943, -0.545, 0.0)));
            methyl.atoms.push_back(create_hydrogen(Vec3(-0.943, -0.545, 0.0)));
            
            methyl.bonds = {{0, 1}, {0, 2}, {0, 3}};
            methyl.connection_points = {0};  // carbon can connect
            fragments["methyl"] = methyl;
        }
        
        // hydroxyl group -OH
        {
            Fragment hydroxyl("hydroxyl", "hydroxyl group -OH");
            hydroxyl.atoms.push_back(create_oxygen(Vec3(0.0, 0.0, 0.0)));
            hydroxyl.atoms.push_back(create_hydrogen(Vec3(0.957, 0.0, 0.0)));
            
            hydroxyl.bonds = {{0, 1}};
            hydroxyl.connection_points = {0};  // oxygen can connect
            fragments["hydroxyl"] = hydroxyl;
        }
        
        // amino group -NH2
        {
            Fragment amino("amino", "amino group -NH2");
            amino.atoms.push_back(create_nitrogen(Vec3(0.0, 0.0, 0.0)));
            amino.atoms.push_back(create_hydrogen(Vec3(0.0, 1.014, 0.0)));
            amino.atoms.push_back(create_hydrogen(Vec3(0.878, -0.507, 0.0)));
            
            amino.bonds = {{0, 1}, {0, 2}};
            amino.connection_points = {0};  // nitrogen can connect
            fragments["amino"] = amino;
        }
        
        // carboxyl group -COOH
        {
            Fragment carboxyl("carboxyl", "carboxyl group -COOH");
            carboxyl.atoms.push_back(create_carbon(Vec3(0.0, 0.0, 0.0)));      // C
            carboxyl.atoms.push_back(create_oxygen(Vec3(1.23, 0.0, 0.0)));     // =O
            carboxyl.atoms.push_back(create_oxygen(Vec3(-0.615, 1.066, 0.0)));  // -OH
            carboxyl.atoms.push_back(create_hydrogen(Vec3(-1.572, 1.066, 0.0))); // H
            
            // set partial charges
            carboxyl.atoms[0].charge = 0.5;   // C
            carboxyl.atoms[1].charge = -0.5;  // =O
            carboxyl.atoms[2].charge = -0.6;  // -O
            carboxyl.atoms[3].charge = 0.4;   // H
            
            carboxyl.bonds = {{0, 1}, {0, 2}, {2, 3}};
            carboxyl.connection_points = {0};  // carbon can connect
            fragments["carboxyl"] = carboxyl;
        }
        
        // phenyl ring -C6H5
        {
            Fragment phenyl("phenyl", "phenyl ring -C6H5");
            
            // benzene ring geometry
            double r = 1.39;
            double angle_step = 2.0 * M_PI / 6.0;
            
            // add carbons
            for (int i = 0; i < 6; ++i) {
                double angle = i * angle_step;
                double x = r * std::cos(angle);
                double y = r * std::sin(angle);
                phenyl.atoms.push_back(create_carbon(Vec3(x, y, 0.0)));
            }
            
            // add hydrogens (except for connection point)
            for (int i = 1; i < 6; ++i) {
                double angle = i * angle_step;
                double x = (r + 1.08) * std::cos(angle);
                double y = (r + 1.08) * std::sin(angle);
                phenyl.atoms.push_back(create_hydrogen(Vec3(x, y, 0.0)));
            }
            
            // ring bonds
            for (int i = 0; i < 6; ++i) {
                int j = (i + 1) % 6;
                phenyl.bonds.push_back({i, j});
            }
            
            // c-h bonds (skip first carbon - connection point)
            for (int i = 1; i < 6; ++i) {
                phenyl.bonds.push_back({i, 5 + i});
            }
            
            phenyl.connection_points = {0};  // first carbon can connect
            fragments["phenyl"] = phenyl;
        }
        
        // methylene group -CH2-
        {
            Fragment methylene("methylene", "methylene group -CH2-");
            methylene.atoms.push_back(create_carbon(Vec3(0.0, 0.0, 0.0)));
            methylene.atoms.push_back(create_hydrogen(Vec3(0.0, 1.09, 0.0)));
            methylene.atoms.push_back(create_hydrogen(Vec3(0.0, -1.09, 0.0)));
            
            methylene.bonds = {{0, 1}, {0, 2}};
            methylene.connection_points = {0};  // carbon can connect on both sides
            fragments["methylene"] = methylene;
        }
        
        // carbonyl group =C=O
        {
            Fragment carbonyl("carbonyl", "carbonyl group =C=O");
            carbonyl.atoms.push_back(create_carbon(Vec3(0.0, 0.0, 0.0)));
            carbonyl.atoms.push_back(create_oxygen(Vec3(1.23, 0.0, 0.0)));
            
            carbonyl.atoms[0].charge = 0.5;
            carbonyl.atoms[1].charge = -0.5;
            
            carbonyl.bonds = {{0, 1}};
            carbonyl.connection_points = {0};  // carbon can connect
            fragments["carbonyl"] = carbonyl;
        }
        
        // phosphate group -PO4
        {
            Fragment phosphate("phosphate", "phosphate group -PO4");
            phosphate.atoms.push_back(create_phosphorus(Vec3(0.0, 0.0, 0.0)));
            
            // tetrahedral geometry
            double r = 1.54;
            phosphate.atoms.push_back(create_oxygen(Vec3(r, 0.0, 0.0)));
            phosphate.atoms.push_back(create_oxygen(Vec3(-r/3.0, r*2.0/3.0*std::sqrt(2), 0.0)));
            phosphate.atoms.push_back(create_oxygen(Vec3(-r/3.0, -r/3.0*std::sqrt(2), r*std::sqrt(6)/3.0)));
            phosphate.atoms.push_back(create_oxygen(Vec3(-r/3.0, -r/3.0*std::sqrt(2), -r*std::sqrt(6)/3.0)));
            
            // set charges
            phosphate.atoms[0].charge = 1.2;
            for (int i = 1; i < 5; ++i) {
                phosphate.atoms[i].charge = -0.8;
            }
            
            phosphate.bonds = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
            phosphate.connection_points = {1, 2, 3, 4};  // oxygens can connect
            fragments["phosphate"] = phosphate;
        }
        
        // sulfate group -SO4
        {
            Fragment sulfate("sulfate", "sulfate group -SO4");
            sulfate.atoms.push_back(create_sulfur(Vec3(0.0, 0.0, 0.0)));
            
            // tetrahedral geometry
            double r = 1.47;  // s-o bond length
            sulfate.atoms.push_back(create_oxygen(Vec3(r, 0.0, 0.0)));
            sulfate.atoms.push_back(create_oxygen(Vec3(-r/3.0, r*2.0/3.0*std::sqrt(2), 0.0)));
            sulfate.atoms.push_back(create_oxygen(Vec3(-r/3.0, -r/3.0*std::sqrt(2), r*std::sqrt(6)/3.0)));
            sulfate.atoms.push_back(create_oxygen(Vec3(-r/3.0, -r/3.0*std::sqrt(2), -r*std::sqrt(6)/3.0)));
            
            // set charges
            sulfate.atoms[0].charge = 1.4;
            for (int i = 1; i < 5; ++i) {
                sulfate.atoms[i].charge = -0.85;
            }
            
            sulfate.bonds = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
            sulfate.connection_points = {1, 2, 3, 4};  // oxygens can connect
            fragments["sulfate"] = sulfate;
        }
    }

public:
    FragmentLibrary() {
        initialize_common_fragments();
    }
    
    static FragmentLibrary& instance() {
        static FragmentLibrary library;
        return library;
    }
    
    bool has_fragment(const std::string& name) const {
        return fragments.find(name) != fragments.end();
    }
    
    const Fragment& get_fragment(const std::string& name) const {
        auto it = fragments.find(name);
        if (it != fragments.end()) {
            return it->second;
        }
        static Fragment empty;
        return empty;
    }
    
    void add_fragment(const std::string& name, const Fragment& fragment) {
        fragments[name] = fragment;
    }
    
    std::vector<std::string> list_fragments() const {
        std::vector<std::string> names;
        for (const auto& pair : fragments) {
            names.push_back(pair.first);
        }
        return names;
    }
    
    // add fragment to system at specified position
    static int add_fragment_to_system(System& system, const std::string& fragment_name, 
                                     const Vec3& position = Vec3(), 
                                     const Vec3& rotation = Vec3()) {
        const Fragment& frag = instance().get_fragment(fragment_name);
        if (frag.num_atoms() == 0) return -1;
        
        int start_idx = system.num_atoms();
        
        // add atoms
        for (const auto& atom : frag.atoms) {
            Atom new_atom = atom;
            new_atom.r = rotate_point(new_atom.r, rotation) + position;
            system.add_atom(new_atom);
        }
        
        // add bonds
        for (const auto& bond : frag.bonds) {
            int i = start_idx + bond.first;
            int j = start_idx + bond.second;
            
            if (i < system.num_atoms() && j < system.num_atoms()) {
                double distance = system.distance_between_atoms(i, j);
                auto bond_params = BondBuilder::BondParameters().get_default_params(
                    system.atoms[i].Z, system.atoms[j].Z);
                system.add_bond(i, j, bond_params.k, distance);
            }
        }
        
        // add angles
        for (const auto& angle : frag.angles) {
            int i = start_idx + std::get<0>(angle);
            int j = start_idx + std::get<1>(angle);
            int k = start_idx + std::get<2>(angle);
            
            if (i < system.num_atoms() && j < system.num_atoms() && k < system.num_atoms()) {
                double current_angle = system.angle_between_atoms(i, j, k);
                system.add_angle(i, j, k, 40.0, current_angle);  // default angle constant
            }
        }
        
        return start_idx;
    }
    
    static void print_available_fragments() {
        const auto& lib = instance();
        std::cout << "Available molecular fragments:\n";
        for (const auto& name : lib.list_fragments()) {
            const auto& frag = lib.get_fragment(name);
            std::cout << "  " << name << ": " << frag.description 
                      << " (" << frag.num_atoms() << " atoms, " 
                      << frag.num_connection_points() << " connection points)\n";
        }
    }

private:
    static Vec3 rotate_point(const Vec3& point, const Vec3& rotation) {
        // simple euler rotation - x, then y, then z
        Vec3 result = point;
        
        // rotate around x-axis
        if (std::abs(rotation.x) > 1e-12) {
            double cos_x = std::cos(rotation.x);
            double sin_x = std::sin(rotation.x);
            double y = result.y * cos_x - result.z * sin_x;
            double z = result.y * sin_x + result.z * cos_x;
            result.y = y;
            result.z = z;
        }
        
        // rotate around y-axis
        if (std::abs(rotation.y) > 1e-12) {
            double cos_y = std::cos(rotation.y);
            double sin_y = std::sin(rotation.y);
            double x = result.x * cos_y + result.z * sin_y;
            double z = -result.x * sin_y + result.z * cos_y;
            result.x = x;
            result.z = z;
        }
        
        // rotate around z-axis
        if (std::abs(rotation.z) > 1e-12) {
            double cos_z = std::cos(rotation.z);
            double sin_z = std::sin(rotation.z);
            double x = result.x * cos_z - result.y * sin_z;
            double y = result.x * sin_z + result.y * cos_z;
            result.x = x;
            result.y = y;
        }
        
        return result;
    }
};

}