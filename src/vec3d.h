#ifndef VEC3D_H
#define VEC3D_H

#include <iostream>
#include <cmath>
#include <optional>

/* Vec3D represents a 3-dimensional vector, or equivalently, a point in 3D space. */
struct Vec3D {
    /* Three `double` components, all 0 by default.
    
    Note: Later, we may no longer set the components to 0 by default with a member default
    value, so that `Vec3D` can be trivially constructible. */
    double x = 0, y = 0, z = 0;

    /* Mathematical negation, +=, -=, *=, /= operators */

    /* Element-wise negation for `Vec3D`s */
    auto operator-() const {return Vec3D{-x, -y, -z};}

    /* Element-wise addition assignment operator for `Vec3D`s */
    auto& operator+= (const Vec3D &rhs) {x += rhs.x; y += rhs.y; z += rhs.z; return *this;}
    /* Element-wise subtraction assignment operator for `Vec3D`s */
    auto& operator-= (const Vec3D &rhs) {x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this;}
    /* Element-wise multiplication assignment operator for `Vec3D`s */
    auto& operator*= (double d) {x *= d; y *= d; z *= d; return *this;}
    /* Element-wise division assignment operator for `Vec3D`s */
    auto& operator/= (double d) {return *this *= (1 / d);}  /* Multiply by 1/d for less divisions */

    /* Compute magnitude (length) of this vector */
    auto mag() const {return std::sqrt(x * x + y * y + z * z);}
    /* Compute squared magnitude (squared length) of this vector */
    auto mag_squared() const {return x * x + y * y + z * z;}

    /* Compute unit vector (forward declared) since it requires operator/, which has
    not yet been defined */

    Vec3D unit_vector() const;
};

/* Math utility functions; vector addition/subtraction, multiplication/division by a scalar */

/* Performs element-wise addition on two `Vec3D`s */
auto operator+ (const Vec3D &a, const Vec3D &b) {auto ret = a; ret += b; return ret;}
/* Performs element-wise subtraction on two `Vec3D`s */
auto operator- (const Vec3D &a, const Vec3D &b) {auto ret = a; ret -= b; return ret;}
/* Performs element-wise multiplication by `d` on a `Vec3D` */
auto operator* (const Vec3D &a, double d) {auto ret = a; ret *= d; return ret;}
/* Performs element-wise multiplication by `d` on a `Vec3D` */
auto operator* (double d, const Vec3D &a) {return a * d;}
/* Performs element-wise division by `d` on a `Vec3D` */
auto operator/ (const Vec3D &a, double d) {auto ret = a; ret /= d; return ret;}

/* Dot and cross product of two vectors. Note that these are not static because
in OOP, static functions ought to not depend on the values of the member variables,
or on the existence of actual instances of the class which they are a static member of.
See https://softwareengineering.stackexchange.com/a/113034/426687. */

/* Computes the dot product of `a` and `b` */
auto dot(const Vec3D &a, const Vec3D &b) {return a.x * b.x + a.y * b.y + a.z * b.z;}
/* Computes the cross product of `a` and `b` */
auto cross(const Vec3D &a, const Vec3D &b) {
    return Vec3D{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

/* Overload operator<< to allow printing `Vec3D`s to output streams */
std::ostream& operator<< (std::ostream& os, const Vec3D &v) {
    os << "(" << v.x << " " << v.y << " " << v.z << ") ";
    return os;
}

/* Returns the unit vector of this `Vec3D`. */
Vec3D Vec3D::unit_vector() const {
    /* The unit vector is found by dividing the vector by its length/magnitude */
    return *this / this->mag();
}

/* Computes the linear interpolation between two `Vec3D`s `a` and `b`, with interpolation
parameter equal to `t`. */
auto vec3d_lerp(const Vec3D &a, const Vec3D &b, double t) {
    return a + t * (b - a);
}

/* `Point3D` is a type alias for `Vec3D`, declared to improve clarity in the code */
using Point3D = Vec3D;

#endif