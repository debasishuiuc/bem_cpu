// geometry_utility.hpp

#pragma once

#include "vec3.hpp"

// Interpolates position, tangents, normal, and geometry for a 6-node triangle
void interp_p(const Vec3& p1, const Vec3& p2, const Vec3& p3,
              const Vec3& p4, const Vec3& p5, const Vec3& p6,
              double al, double be, double ga,
              double xi, double eta,
              Vec3& x,             // interpolated position
              Vec3& dx_dxi,        // tangent vector in xi-direction
              Vec3& dx_deta,       // tangent vector in eta-direction
              Vec3& n,             // unit normal vector
              double& hxi,         // |dx_dxi|
              double& het,         // |dx_deta|
              double& hs,          // |n| (area element)
              int Ichoose);        // 1 = position only, >1 = include tangents and normal
