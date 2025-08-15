#pragma once

#include "vec3.hpp"

namespace carbon {

struct Atom {
    // atomic properties
    int Z;           // atomic number
    double mass;     // atomic mass (amu)
    double charge;   // partial charge (e)
    
    // lennard-jones parameters
    double sigma;    // lj diameter parameter (angstroms)
    double epsilon;  // lj well depth parameter (kcal/mol)
    
    // state variables
    Vec3 r;          // position (angstroms)
    Vec3 v;          // velocity (angstroms/fs)
    Vec3 f;          // force (kcal/mol/angstrom)
    
    Atom() : Z(0), mass(0.0), charge(0.0), sigma(0.0), epsilon(0.0),
             r(), v(), f() {}
    
    Atom(int Z_, double mass_, double charge_ = 0.0, 
         double sigma_ = 0.0, double epsilon_ = 0.0,
         const Vec3& r_ = Vec3(), const Vec3& v_ = Vec3())
        : Z(Z_), mass(mass_), charge(charge_), sigma(sigma_), epsilon(epsilon_),
          r(r_), v(v_), f() {}
    
    void clear_forces() {
        f = Vec3();
    }
    
    void add_force(const Vec3& force) {
        f += force;
    }
    
    double kinetic_energy() const {
        return 0.5 * mass * v.magnitude_squared();
    }
    
    bool is_hydrogen() const { return Z == 1; }
    bool is_carbon() const { return Z == 6; }
    bool is_nitrogen() const { return Z == 7; }
    bool is_oxygen() const { return Z == 8; }
    bool is_sodium() const { return Z == 11; }
    bool is_chlorine() const { return Z == 17; }
};

}