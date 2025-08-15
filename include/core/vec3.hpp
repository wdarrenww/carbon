#pragma once

#include <cmath>
#include <ostream>

namespace carbon {

struct Vec3 {
    double x, y, z;
    
    Vec3() : x(0.0), y(0.0), z(0.0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }
    
    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }
    
    Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }
    
    Vec3 operator/(double scalar) const {
        return Vec3(x / scalar, y / scalar, z / scalar);
    }
    
    Vec3& operator+=(const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    
    Vec3& operator-=(const Vec3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    
    Vec3& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }
    
    double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }
    
    double magnitude_squared() const {
        return x * x + y * y + z * z;
    }
    
    Vec3 normalized() const {
        double mag = magnitude();
        if (mag > 1e-12) {
            return *this / mag;
        }
        return Vec3();
    }
    
    void normalize() {
        double mag = magnitude();
        if (mag > 1e-12) {
            *this /= mag;
        }
    }
};

inline Vec3 operator*(double scalar, const Vec3& vec) {
    return vec * scalar;
}

inline std::ostream& operator<<(std::ostream& os, const Vec3& vec) {
    os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return os;
}

// distance functions
inline double distance(const Vec3& a, const Vec3& b) {
    return (a - b).magnitude();
}

inline double distance_squared(const Vec3& a, const Vec3& b) {
    return (a - b).magnitude_squared();
}

// periodic boundary condition utilities
inline Vec3 apply_pbc(const Vec3& pos, const Vec3& box_size) {
    Vec3 result = pos;
    result.x = pos.x - box_size.x * std::floor(pos.x / box_size.x);
    result.y = pos.y - box_size.y * std::floor(pos.y / box_size.y);
    result.z = pos.z - box_size.z * std::floor(pos.z / box_size.z);
    return result;
}

inline Vec3 minimum_image_vector(const Vec3& r1, const Vec3& r2, const Vec3& box_size) {
    Vec3 dr = r2 - r1;
    dr.x = dr.x - box_size.x * std::round(dr.x / box_size.x);
    dr.y = dr.y - box_size.y * std::round(dr.y / box_size.y);
    dr.z = dr.z - box_size.z * std::round(dr.z / box_size.z);
    return dr;
}

}