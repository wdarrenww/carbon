#pragma once

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace carbon {

struct Angle {
    int i, j, k;     // atom indices: i-j-k angle with j as central atom
    double ka;       // angle spring constant (kcal/mol/rad^2)
    double theta0;   // equilibrium angle (radians)
    
    Angle() : i(0), j(0), k(0), ka(0.0), theta0(0.0) {}
    
    Angle(int i_, int j_, int k_, double ka_, double theta0_)
        : i(i_), j(j_), k(k_), ka(ka_), theta0(theta0_) {}
    
    bool involves_atom(int atom_idx) const {
        return i == atom_idx || j == atom_idx || k == atom_idx;
    }
    
    int get_central_atom() const {
        return j;
    }
    
    bool is_valid() const {
        return i != j && j != k && i != k;
    }
    
    // convert angle to degrees for display
    double theta0_degrees() const {
        return theta0 * 180.0 / M_PI;
    }
    
    static double degrees_to_radians(double degrees) {
        return degrees * M_PI / 180.0;
    }
    
    static double radians_to_degrees(double radians) {
        return radians * 180.0 / M_PI;
    }
};

}