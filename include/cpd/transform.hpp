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
/// Generic coherent point drift transform.
///
/// Downstreams shouldn't need to include this file directly, use a realization
/// of a transform (e.g. `Rigid`) instead.

#pragma once

#include <iostream>
#include <chrono>
#include <memory>
#include <cpd/gauss_transform.hpp>
#include <cpd/matrix.hpp>

namespace cpd {

/// The default number of iterations allowed.
const size_t DEFAULT_MAX_ITERATIONS = 100;
/// The default outlier weight.
const double DEFAULT_OUTLIERS = 0.5;
/// The default tolerance.
const double DEFAULT_TOLERANCE = 1e-5;
/// The default sigma2.
const double DEFAULT_SIGMA2 = 0.0;

/// The result of a generic transform run.
struct Result {
	/// The final moved points.
	Matrix points;
	/// The final sigma2 value.
	double sigma2;
	/// The runtime.
	std::chrono::microseconds runtime;
	/// The number of iterations.
	size_t iterations;
};

/// Computes the default sigma2 for the given matrices.
double computeSigma2(const Matrix& fixed, const Matrix& moving);

/// Generic coherent point drift transform.
///
/// An abstract base class for real transforms, e.g. `Rigid` or `Nonrigid`.
template<typename Result>
class Transform {
public:
	virtual ~Transform() {
	}

	Transform()
			: m_GaussTransform(GaussTransform::makeDefault()),
				m_MaxIterations(DEFAULT_MAX_ITERATIONS), m_outliers(DEFAULT_OUTLIERS),
				m_sigma2(DEFAULT_SIGMA2), m_tolerance(DEFAULT_TOLERANCE) {
	}

	/// Sets the gauss transform.
	Transform& gauss_transform(std::unique_ptr<GaussTransform> gauss_transform) {
		this->m_GaussTransform = std::move(gauss_transform);
		return *this;
	}

	/// Sets the max iterations for this transform.
	Transform& max_iterations(double max_iterations) {
		this->m_MaxIterations = max_iterations;
		return *this;
	}

	/// Sets the outlier tolerance.
	Transform& outliers(double outliers) {
		this->m_outliers = outliers;
		return *this;
	}

	/// Sets the sigma2 value for this transform.
	Transform& sigma2(double sigma2) {
		this->m_sigma2 = sigma2;
		return *this;
	}

	/// Sets the final tolerance.
	Transform& tolerance(double tolerance) {
		this->m_tolerance = tolerance;
		return *this;
	}

	/// Runs this transform for the provided matrices.
	Result run(Matrix fixed, Matrix moving) {
		auto tic = std::chrono::high_resolution_clock::now();

		this->init(fixed, moving);

		Result result;
		result.points = moving;
		if (this->m_sigma2 <= 0.0) {
			result.sigma2 = cpd::computeSigma2(fixed, moving);
		} else {
			result.sigma2 = this->m_sigma2;
		}

		size_t iter = 0;
		double ntol = this->m_tolerance + 10.0;
		double l = 0.;
		while (iter < this->m_MaxIterations && ntol > this->m_tolerance
					 && result.sigma2 > 10 * std::numeric_limits<double>::epsilon()) {

//			auto ticEM = std::chrono::high_resolution_clock::now();
			clock_t clockStart = clock();

			Probabilities P = this->m_GaussTransform->computeEStep(fixed, result.points,
																														 result.sigma2,
																														 this->m_outliers);
			this->modifyProbabilities(P);

			ntol = std::abs((P.l - l) / P.l);
			l = P.l;

			result = this->computeMStep(fixed, moving, P, result.sigma2);

//			auto tocEM = std::chrono::high_resolution_clock::now();
//			// fractional duration: no duration_cast needed
//			std::chrono::duration<double, std::milli> runetime_ms = tocEM - ticEM;
			double execTime = (double) (clock() - clockStart) / CLOCKS_PER_SEC;
			++iter;

			std::cout << iter << " of " << this->m_MaxIterations << ": dL = " << ntol
								<< ", sigma2 = " << result.sigma2 << ", " << execTime << " sec (CPU)"
								<< std::endl;
//								<< ", sigma2 = " << result.sigma2 << ", " << runetime_ms.count()
//								<< " ms (CPU)" << std::endl;
		}

		auto toc = std::chrono::high_resolution_clock::now();
		// integral duration: requires duration_cast
		result.runtime = std::chrono::duration_cast < std::chrono::microseconds
										 > (toc - tic);
		result.iterations = iter;
		return result;
	}

	/// Initialize this transform for the provided matrices.
	///
	/// This is called before beginning each run, but after normalization. In
	/// general, transforms shouldn't need to be initialized, but the nonrigid
	/// does.
	virtual void init(const Matrix& fixed, const Matrix& moving) {
	}

	/// Modifies `Probabilities` in some way.
	///
	/// Some types of transform need to tweak the probabilities before moving on
	/// with an interation. The default behavior is to do nothing.
	virtual void modifyProbabilities(Probabilities& probabilities) const {
	}

	/// Computes one iteration of the transform.
	virtual Result computeMStep(const Matrix& fixed, const Matrix& moving,
															const Probabilities& probabilities,
															double sigma2) const = 0;

private:
	std::unique_ptr<GaussTransform> m_GaussTransform;
	size_t m_MaxIterations;
	double m_outliers;
	double m_sigma2;
	double m_tolerance;
};
}  // namespace cpd
