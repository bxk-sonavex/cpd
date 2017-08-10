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

// #include <cpd/jsoncpp.hpp>
#include <cpd/rigid.hpp>
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
	double execTime;
	clock_t clockStart;

	if (!(argc == 4 || argc == 5)) {
		std::cout << "ERROR: invalid usage" << std::endl;
		std::cout << "USAGE: cpd-rigid <sigma> <fixed> <moving> [outfile]"
							<< std::endl;
		return 1;
	}

	clockStart = clock();
	float sigma = std::atof(argv[1]);
	cpd::Matrix fixed = cpd::matrix_from_path(argv[2]);
	cpd::Matrix moving = cpd::matrix_from_path(argv[3]);

	std::cout << argv[2] << ": " << fixed.rows() << std::endl;
	std::cout << argv[3] << ": " << moving.rows() << std::endl;
  std::cout << "sigma = " << sigma << std::endl;
  
  execTime = (double) (clock() - clockStart) / CLOCKS_PER_SEC;
  std::cout << "Loading data took " << execTime << " sec (CPU)" << std::endl;

	clockStart = clock();
	cpd::Rigid rigid;
	rigid.sigma2(sigma);
	cpd::RigidResult result = rigid.run(fixed, moving);
	// std::cout << cpd::to_json(result) << std::endl;
	execTime = (double) (clock() - clockStart) / CLOCKS_PER_SEC;
	std::cout << "Registration took " << execTime << " sec (CPU)" << std::endl;

  cpd::Matrix transform = result.matrix();
  std::cout << "Transformation Matrix" << std::endl;
  std::cout << transform << std::endl;

	if (argc == 5) {
    std::cout << "Save result to ... " << argv[4];
		std::ofstream outfile(argv[4]);
		outfile << result.points << std::endl;
		outfile.close();
    std::cout << "done" << std::endl;
	}

	return 0;
}
