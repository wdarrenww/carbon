#pragma once

#include <unordered_map>
#include <string>
#include "atom.hpp"

namespace carbon {

struct ElementData {
    int Z;                    // atomic number
    double mass;              // atomic mass (amu)
    double covalent_radius;   // covalent radius (angstroms)
    double vdw_radius;        // van der waals radius (angstroms)
    double sigma;             // lj sigma parameter (angstroms)
    double epsilon;           // lj epsilon parameter (kcal/mol)
    std::string symbol;       // element symbol
    std::string name;         // element name
    
    ElementData() : Z(0), mass(0.0), covalent_radius(0.0), vdw_radius(0.0),
                    sigma(0.0), epsilon(0.0), symbol(""), name("") {}
    
    ElementData(int Z_, double mass_, double cov_r, double vdw_r,
                double sig, double eps, const std::string& sym, const std::string& nm)
        : Z(Z_), mass(mass_), covalent_radius(cov_r), vdw_radius(vdw_r),
          sigma(sig), epsilon(eps), symbol(sym), name(nm) {}
};

class ElementTable {
private:
    std::unordered_map<int, ElementData> elements_by_z;
    std::unordered_map<std::string, ElementData> elements_by_symbol;
    
    void initialize_default_elements() {
        // hydrogen
        add_element(1, 1.008, 0.31, 1.20, 2.51, 0.016, "H", "Hydrogen");
        
        // carbon
        add_element(6, 12.011, 0.76, 1.70, 3.40, 0.086, "C", "Carbon");
        
        // nitrogen
        add_element(7, 14.007, 0.71, 1.55, 3.25, 0.170, "N", "Nitrogen");
        
        // oxygen
        add_element(8, 15.999, 0.66, 1.52, 3.12, 0.210, "O", "Oxygen");
        
        // fluorine
        add_element(9, 18.998, 0.57, 1.47, 2.94, 0.061, "F", "Fluorine");
        
        // sodium
        add_element(11, 22.990, 1.66, 2.27, 2.58, 0.470, "Na", "Sodium");
        
        // magnesium
        add_element(12, 24.305, 1.41, 1.73, 2.69, 0.870, "Mg", "Magnesium");
        
        // phosphorus
        add_element(15, 30.974, 1.07, 1.80, 3.74, 0.200, "P", "Phosphorus");
        
        // sulfur
        add_element(16, 32.065, 1.05, 1.80, 3.60, 0.250, "S", "Sulfur");
        
        // chlorine
        add_element(17, 35.453, 0.99, 1.75, 3.52, 0.265, "Cl", "Chlorine");
        
        // potassium
        add_element(19, 39.098, 2.03, 2.75, 3.40, 0.650, "K", "Potassium");
        
        // calcium
        add_element(20, 40.078, 1.76, 2.31, 3.03, 1.200, "Ca", "Calcium");
        
        // bromine
        add_element(35, 79.904, 1.20, 1.85, 3.97, 0.320, "Br", "Bromine");
        
        // iodine
        add_element(53, 126.904, 1.39, 1.98, 4.45, 0.400, "I", "Iodine");
    }
    
    void add_element(int Z, double mass, double cov_r, double vdw_r,
                     double sig, double eps, const std::string& sym, const std::string& name) {
        ElementData data(Z, mass, cov_r, vdw_r, sig, eps, sym, name);
        elements_by_z[Z] = data;
        elements_by_symbol[sym] = data;
    }

public:
    ElementTable() {
        initialize_default_elements();
    }
    
    static ElementTable& instance() {
        static ElementTable table;
        return table;
    }
    
    bool has_element(int Z) const {
        return elements_by_z.find(Z) != elements_by_z.end();
    }
    
    bool has_element(const std::string& symbol) const {
        return elements_by_symbol.find(symbol) != elements_by_symbol.end();
    }
    
    const ElementData& get_element(int Z) const {
        auto it = elements_by_z.find(Z);
        if (it != elements_by_z.end()) {
            return it->second;
        }
        static ElementData empty;
        return empty;
    }
    
    const ElementData& get_element(const std::string& symbol) const {
        auto it = elements_by_symbol.find(symbol);
        if (it != elements_by_symbol.end()) {
            return it->second;
        }
        static ElementData empty;
        return empty;
    }
    
    Atom create_atom(int Z, const Vec3& position = Vec3(), const Vec3& velocity = Vec3()) const {
        const ElementData& data = get_element(Z);
        return Atom(data.Z, data.mass, 0.0, data.sigma, data.epsilon, position, velocity);
    }
    
    Atom create_atom(const std::string& symbol, const Vec3& position = Vec3(), const Vec3& velocity = Vec3()) const {
        const ElementData& data = get_element(symbol);
        return Atom(data.Z, data.mass, 0.0, data.sigma, data.epsilon, position, velocity);
    }
    
    double get_mass(int Z) const {
        return get_element(Z).mass;
    }
    
    double get_mass(const std::string& symbol) const {
        return get_element(symbol).mass;
    }
    
    double get_covalent_radius(int Z) const {
        return get_element(Z).covalent_radius;
    }
    
    double get_covalent_radius(const std::string& symbol) const {
        return get_element(symbol).covalent_radius;
    }
    
    std::string get_symbol(int Z) const {
        return get_element(Z).symbol;
    }
    
    void print_elements() const {
        std::cout << "Available elements:\n";
        for (const auto& pair : elements_by_z) {
            const ElementData& data = pair.second;
            std::cout << "  " << data.symbol << " (" << data.Z << "): " 
                      << data.name << ", mass=" << data.mass 
                      << ", sigma=" << data.sigma << ", epsilon=" << data.epsilon << "\n";
        }
    }
};

// convenience functions
inline const ElementTable& get_element_table() {
    return ElementTable::instance();
}

inline Atom create_hydrogen(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(1, pos);
}

inline Atom create_carbon(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(6, pos);
}

inline Atom create_nitrogen(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(7, pos);
}

inline Atom create_oxygen(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(8, pos);
}

inline Atom create_fluorine(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(9, pos);
}

inline Atom create_sodium(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(11, pos);
}

inline Atom create_magnesium(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(12, pos);
}

inline Atom create_phosphorus(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(15, pos);
}

inline Atom create_sulfur(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(16, pos);
}

inline Atom create_chlorine(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(17, pos);
}

inline Atom create_potassium(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(19, pos);
}

inline Atom create_calcium(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(20, pos);
}

inline Atom create_bromine(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(35, pos);
}

inline Atom create_iodine(const Vec3& pos = Vec3()) {
    return get_element_table().create_atom(53, pos);
}

}