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

#include <cpd/transform.hpp>

namespace cpd {

double computeSigma2(const Matrix& fixed, const Matrix& moving) {
	return ((moving.rows() * (fixed.transpose() * fixed).trace())
					+ (fixed.rows() * (moving.transpose() * moving).trace())
					- 2 * fixed.colwise().sum() * moving.colwise().sum().transpose())
					/ (fixed.rows() * moving.rows() * fixed.cols());
}

}  // namespace cpd
