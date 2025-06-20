// vec3.hpp (CPU-only version)

#pragma once

#include <cmath>

struct Vec3 {
  double x, y, z;

  // Constructors
  Vec3() : x(0.0), y(0.0), z(0.0) {}
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  // Arithmetic operators
  Vec3 operator+(const Vec3& rhs) const {
    return Vec3(x + rhs.x, y + rhs.y, z + rhs.z);
  }

  Vec3& operator+=(const Vec3& rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }

  Vec3 operator-(const Vec3& rhs) const {
    return Vec3(x - rhs.x, y - rhs.y, z - rhs.z);
  }

  Vec3& operator-=(const Vec3& rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;
    return *this;
  }

  Vec3 operator*(double scalar) const {
    return Vec3(x * scalar, y * scalar, z * scalar);
  }

  Vec3 operator/(double scalar) const {
    return Vec3(x / scalar, y / scalar, z / scalar);
  }

  Vec3& operator*=(double scalar) {
    x *= scalar; y *= scalar; z *= scalar;
    return *this;
  }

  Vec3& operator/=(double scalar) {
    x /= scalar; y /= scalar; z /= scalar;
    return *this;
  }

  // Dot product
  double dot(const Vec3& rhs) const {
    return x * rhs.x + y * rhs.y + z * rhs.z;
  }

  // Cross product
  Vec3 cross(const Vec3& rhs) const {
    return Vec3(
		y * rhs.z - z * rhs.y,
		z * rhs.x - x * rhs.z,
		x * rhs.y - y * rhs.x
		);
  }

  // Norm (length)
  double norm() const {
    return std::sqrt(x * x + y * y + z * z);
  }

  // Return a normalized copy
  Vec3 normalized() const {
    double n = norm();
    return (n > 1e-14) ? (*this / n) : Vec3(0, 0, 0);
  }

  // Index access
  double& operator[](int i) { return *(&x + i); }
  const double& operator[](int i) const { return *(&x + i); }

  // Named accessors
  double& x_() { return x; }
  double& y_() { return y; }
  double& z_() { return z; }

  const double& x_() const { return x; }
  const double& y_() const { return y; }
  const double& z_() const { return z; }
};

// Global scalar * Vec3
inline Vec3 operator*(double a, const Vec3& v) {
  return Vec3(a * v.x, a * v.y, a * v.z);
}
