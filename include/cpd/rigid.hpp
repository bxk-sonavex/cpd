// cpd - Coherent Point Drift
// Copyright (C) 2017 Pete Gadomski <pete.gadomski@gmail.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

/// \file
///
/// Rigid coherent point drift.

#pragma once

#include <cpd/transform.hpp>

namespace cpd {

/// Should rigid registrations allow reflections by default?
const bool DEFAULT_REFLECTIONS = false;

/// The result of a rigid coherent point drift run.
struct RigidResult: public Result {
	/// The rotation component of the transformation.
	Matrix rotation;
	/// The translation component of the transformation.
	Vector translation;
	/// The scaling component of the transformation.
	double scale;
	/// Returns a single matrix that contains all the transformation
	/// information.
	Matrix matrix() const;
};

/// Rigid coherent point drift.
///
/// Scaling and reflections can be turned on and off.
class Rigid: public Transform<RigidResult> {
public:
	Rigid();

	/// Sets whether this rigid transform allows reflections.
	Rigid& setReflections(bool reflections);

	/// Computes one iteration of the rigid transformation.
	RigidResult computeMStep(const Matrix& fixed, const Matrix& moving,
													 const Probabilities& probabilities,
													 double sigma2) const;

private:
	bool m_reflections;
};

/// Runs a rigid registration on two matrices.
RigidResult rigid(const Matrix& fixed, const Matrix& moving);

}  // namespace cpd
