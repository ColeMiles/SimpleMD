#ifndef VEC3_H
#define VEC3_H

#include <cmath>

struct vec3 {
    double x;
    double y;
    double z;

    vec3(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}

    inline vec3& operator=(const vec3& other) {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }

    inline vec3 operator+(const vec3& other) const {
        return vec3(x + other.x, y + other.y, z + other.z);
    }

    inline vec3& operator+=(const vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // Negation
    inline vec3 operator-() const{
        return vec3(-x, -y, -z);
    }

    // Subtraction
    inline vec3 operator-(const vec3& other) const {
        return vec3(x - other.x, y - other.y, z - other.z);
    }

    inline vec3& operator-=(const vec3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    inline vec3 operator*(const double scale) const {
        return vec3(x * scale, y * scale, z * scale);
    }

    inline vec3& operator*=(const double scale) {
        x *= scale;
        y *= scale;
        z *= scale;
        return *this;
    }

    inline vec3 operator/(const double scale) const {
        return vec3(x / scale, y / scale, z / scale);
    }

    inline vec3& operator/=(const double scale) {
        x /= scale;
        y /= scale;
        z /= scale;
        return *this;
    }

    inline double square_mag() const {
        return x * x + y * y + z * z;
    }

    inline double mag() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    inline vec3 normalize() const {
        double magnitude = mag();
        return *this / magnitude;
    }

    inline double dot(const vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    inline vec3 project(const vec3 & other) const {
        return other.normalize() * (*this).dot(other);
    }
};

// Left multiply override
inline vec3 operator*(const double scale, const vec3& vec) {
    return vec * scale;
}

// I feel like this could be handled better, but I just want these
//   sectioned away in case another library defines a "dot" or 
//   "cross" function
namespace vec3ops {
    inline double dot(const vec3& a, const vec3& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    inline vec3 cross(const vec3& a, const vec3& b) {
        return vec3(a.y * b.z - a.z * b.y, 
                    a.z * b.x - a.x * b.z,
                    a.x * b.y - a.y * b.x);
    }
}

#endif