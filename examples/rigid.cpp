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

#include <pcl/io/pcd_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/console/print.h>
#include <pcl/console/parse.h>
#include <pcl/console/time.h>

float default_sigma = 0.0f;
float default_leaf_size = 0.1f;
std::string default_field("z");
double default_filter_min = -std::numeric_limits<double>::max();
double default_filter_max = std::numeric_limits<double>::max();

void printHelp(int, char **argv) {
	pcl::console::print_error("Syntax is: %s fixed.pcd moving.pcd outfile.pcd <options>\n",
														argv[0]);
	pcl::console::print_info("  options:\n");
	pcl::console::print_info("    -sigma X     = the VoxelGrid leaf size (default: ");
	pcl::console::print_value("%f", default_sigma);
	pcl::console::print_info(")\n");
	pcl::console::print_info("    -leaf x,y,z  = the VoxelGrid leaf size (default: ");
	pcl::console::print_value("%f, %f, %f", default_leaf_size, default_leaf_size,
														default_leaf_size);
	pcl::console::print_info(")\n");
	pcl::console::print_info("    -field X     = filter data along this field "
													 "name (default: ");
	pcl::console::print_value("%s", default_field.c_str());
	pcl::console::print_info(")\n");
	pcl::console::print_info("    -fmin  X     = filter all data with values "
													 "along the specified field smaller than this value "
													 "(default: ");
	pcl::console::print_value("-inf");
	pcl::console::print_info(")\n");
	pcl::console::print_info("    -fmax  X     = filter all data with values along "
													 "the specified field larger than this value (default: ");
	pcl::console::print_value("inf");
	pcl::console::print_info(")\n");
}

void downsample(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr &input,
								pcl::PointCloud<pcl::PointXYZ> &output,
								float leaf_x, float leaf_y, float leaf_z, const std::string &field,
								double fmin,
								double fmax) {
	pcl::console::TicToc tt;
	tt.tic();

	pcl::console::print_highlight("Downsampling ");

	pcl::VoxelGrid<pcl::PointXYZ> grid;
	grid.setInputCloud(input);
	grid.setFilterFieldName(field);
	grid.setFilterLimits(fmin, fmax);
	grid.setLeafSize(leaf_x, leaf_y, leaf_z);
	grid.filter(output);

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms : ");
	pcl::console::print_value("%d", output.width * output.height);
	pcl::console::print_info(" points]\n");

	std::cout << "width: " << output.width << ", height: " << output.height << std::endl;
}

bool loadCloud(const std::string &filename,
							 pcl::PointCloud<pcl::PointXYZ> &cloud) {
	pcl::console::TicToc tt;
	pcl::console::print_highlight("Loading ");
	pcl::console::print_value("%s ", filename.c_str());

	tt.tic();
	if (pcl::io::loadPCDFile(filename, cloud) < 0) {
		return false;
	}

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms : ");
	pcl::console::print_value("%d", cloud.width * cloud.height);
	pcl::console::print_info(" points]\n");
	pcl::console::print_info("Available dimensions: ");
	pcl::console::print_value("%s\n", pcl::getFieldsList(cloud).c_str());

	std::cout << "width: " << cloud.width << ", height: " << cloud.height << std::endl;

	return true;
}

void saveCloud(const std::string &filename,
							 const pcl::PointCloud<pcl::PointXYZ> &output) {
	pcl::console::TicToc tt;
	tt.tic();

	pcl::console::print_highlight("Saving ");
	pcl::console::print_value("%s ", filename.c_str());

	pcl::PCDWriter w;
	w.writeBinaryCompressed(filename, output);

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms : ");
	pcl::console::print_value("%d", output.width * output.height);
	pcl::console::print_info(" points]\n");
}

int main(int argc, char** argv) {
	double execTime;
	clock_t clockStart;

	if (argc < 4) {
		printHelp(argc, argv);
		return -1;
	}

	// Parse the command line arguments for .pcd files
	std::vector<int> idxFiles;
	idxFiles = pcl::console::parse_file_extension_argument(argc, argv, ".pcd");
	if (idxFiles.size() != 3) {
		pcl::console::print_error("Need two input files and one output file.\n");
		return -1;
	}

	// Command line parsing
	float sigma = default_sigma;
	pcl::console::parse_argument(argc, argv, "-sigma", sigma);
	pcl::console::print_info("Using a sigma of: ");
	pcl::console::print_value("%f\n", sigma);

	float leaf_x = default_leaf_size,
			leaf_y = default_leaf_size,
			leaf_z = default_leaf_size;

	std::vector<double> values;
	pcl::console::parse_x_arguments(argc, argv, "-leaf", values);
	if (values.size() == 1) {
		leaf_x = static_cast<float>(values[0]);
		leaf_y = static_cast<float>(values[0]);
		leaf_z = static_cast<float>(values[0]);
	}
	else if (values.size() == 3) {
		leaf_x = static_cast<float>(values[0]);
		leaf_y = static_cast<float>(values[1]);
		leaf_z = static_cast<float>(values[2]);
	}
	else {
		pcl::console::print_error(
				"Leaf size must be specified with either 1 or 3 numbers (%lu given).\n",
				values.size());
	}
	pcl::console::print_info("Using a leaf size of: ");
	pcl::console::print_value("%f, %f, %f\n", leaf_x, leaf_y, leaf_z);

	std::string field(default_field);
	pcl::console::parse_argument(argc, argv, "-field", field);
	double fmin = default_filter_min,
			fmax = default_filter_max;
	pcl::console::parse_argument(argc, argv, "-fmin", fmin);
	pcl::console::parse_argument(argc, argv, "-fmax", fmax);
	pcl::console::print_info("Filtering data on field: ");
	pcl::console::print_value("%s", field.c_str());
	pcl::console::print_info(" between: ");
	if (fmin == -std::numeric_limits<double>::max()) {
		pcl::console::print_value("-inf ->");
	}
	else {
		pcl::console::print_value("%f ->", fmin);
	}

	if (fmax == std::numeric_limits<double>::max()) {
		pcl::console::print_value("inf\n");
	}
	else {
		pcl::console::print_value("%f\n", fmax);
	}

	// Load the first file (fixed point cloud)
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcFixed(new pcl::PointCloud<pcl::PointXYZ>);
	if (!loadCloud(argv[idxFiles[0]], *pcFixed)) {
		return -1;
	}
	// Load the second file (moving point cloud)
	pcl::PointCloud<pcl::PointXYZ> pcMoving;
	if (!loadCloud(argv[idxFiles[1]], pcMoving)) {
		return -1;
	}

	// Apply the voxel grid downsampling
	pcl::PointCloud<pcl::PointXYZ> pcFixedDownsampled;
	downsample(pcFixed, pcFixedDownsampled, leaf_x, leaf_y, leaf_z, field, fmin, fmax);

	// Save into the third file
//	saveCloud(argv[idxFiles[2]], pcFixedDownsampled);

// Get Eigen matrix
	cpd::Matrix ptsMoving = pcMoving.getMatrixXfMap().transpose();
	ptsMoving.conservativeResize(ptsMoving.rows(), ptsMoving.cols() - 1);
//	std::cout << "(Eigen) #row : " << ptsMoving.rows()
//						<< " #col : " << ptsMoving.cols() << std::endl;
	cpd::Matrix ptsFixed = pcFixedDownsampled.getMatrixXfMap().transpose();
	ptsFixed.conservativeResize(ptsFixed.rows(), ptsFixed.cols() - 1);

	// Perform rigid registration
	clockStart = clock();
	cpd::Rigid rigid;
	rigid.sigma2(sigma);
	cpd::RigidResult result = rigid.run(ptsFixed, ptsMoving);
	// std::cout << cpd::to_json(result) << std::endl;
	execTime = (double) (clock() - clockStart) / CLOCKS_PER_SEC;
	std::cout << "Registration took " << execTime << " sec (CPU)" << std::endl;

	cpd::Matrix transform = result.matrix();
	std::cout << "Transformation Matrix" << std::endl;
	std::cout << transform << std::endl;

	// Save into the third file
//	saveCloud(argv[idxFiles[2]], pcFixedDownsampled);
	std::cout << "Save result to ... [" << argv[idxFiles[2]];
	std::ofstream outfile(argv[idxFiles[2]]);
	outfile << result.points << std::endl;
	outfile.close();
	std::cout << "] done" << std::endl;

	return 0;
}

// EOF
