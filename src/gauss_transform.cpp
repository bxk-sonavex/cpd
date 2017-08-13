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

#define _USE_MATH_DEFINES

#include <cmath>
#include <cpd/gauss_transform.hpp>

namespace cpd {

GaussTransform::~GaussTransform() {
}

std::unique_ptr<GaussTransform> GaussTransform::makeDefault() {
	return std::unique_ptr < GaussTransform > (new GaussTransformDirect());
}

Probabilities GaussTransformDirect::computeEStep(const Matrix& fixed,
																								 const Matrix& moving,
																								 double sigma2, double outliers) const {
	double ksig = -2.0 * sigma2;
	double truncate = log(1e-3);
	size_t cols = fixed.cols();
	double outlier_tmp = (outliers * moving.rows()
												* std::pow(-ksig * M_PI, 0.5 * cols))
											 / ((1 - outliers) * fixed.rows());
	Vector p = Vector::Zero(moving.rows());
	Vector p1 = Vector::Zero(moving.rows());
	Vector pt1 = Vector::Zero(fixed.rows());
	Matrix px = Matrix::Zero(moving.rows(), cols);
	double l = 0.0;

	for (Matrix::Index i = 0; i < fixed.rows(); ++i) {
		double sp = 0;
		for (Matrix::Index j = 0; j < moving.rows(); ++j) {
			double razn = (fixed.row(i) - moving.row(j)).array().pow(2).sum();
			razn = razn / ksig;
			if (razn < truncate) {
				p(j) = 0;
			} else {
				p(j) = std::exp(razn);
				sp += p(j);
			}
		}

		sp += outlier_tmp;
		pt1(i) = 1 - outlier_tmp / sp;
		for (Matrix::Index j = 0; j < moving.rows(); ++j) {
			p1(j) += p(j) / sp;
			px.row(j) += fixed.row(i) * p(j) / sp;
		}
		l += -std::log(sp);
	}
	l += cols * fixed.rows() * std::log(sigma2) / 2;

	return {p1, pt1, px, l};
}

}  // namespace cpd
