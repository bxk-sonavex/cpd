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
/// Basic typedefs for our flavors of Eigen matrices.

#pragma once

#include <Eigen/Dense>

namespace cpd {

/// Our base matrix class.
typedef Eigen::MatrixXf Matrix;

/// Typedef for our specific type of vector.
typedef Eigen::VectorXf Vector;

/// Typedef for our specific type of array.
typedef Eigen::ArrayXf Array;

/// Apply a transformation matrix to a set of points.
///
/// The transformation matrix should be one column wider than the point matrix.
Matrix applyTransformation(Matrix points, const Matrix& transform);

/// Loads a matrix from a delimited text file.
Matrix loadMatrixFromFile(const std::string& path, const char delimiter=' ');

}  // namespace cpd
