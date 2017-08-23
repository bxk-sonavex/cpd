#include <boost/thread/thread.hpp>
#include <Eigen/Dense>
#include <cpd/rigid.hpp>
#include <pcl/io/ply_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/search/kdtree.h>
#include <pcl/console/print.h>
#include <pcl/console/parse.h>
#include <pcl/console/time.h>
#include <pcl/visualization/pcl_visualizer.h>

float default_sigma = 0.0f;
float default_leaf_size = 0.1f;
std::string default_field("z");
double default_filter_min = -std::numeric_limits<double>::max();
double default_filter_max = std::numeric_limits<double>::max();

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp) {
	// The machine epsilon has to be scaled to the magnitude of the values used
	// and multiplied by the desired precision in ULPs (units in the last place)
	return std::abs(x - y) < std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
				 // unless the result is subnormal
			|| std::abs(x - y) < std::numeric_limits<T>::min();
}

/***
 * Checks if a matrix is a valid rotation matrix.
 * @param
 */
bool isRotationMatrix(const cpd::Matrix &R) {
	cpd::Matrix Rt = R.transpose();
	cpd::Matrix shouldBeIdentity = Rt * R;
	cpd::Matrix diff = shouldBeIdentity - Eigen::Matrix<float, 3, 3>::Identity();

//	return diff.norm() < 1e-6;
	return diff.norm() < std::numeric_limits<float>::epsilon();
}

/***
 * Function: Get from a cv::Mat to a pcl::PointCloud<pcl::PointXYZ>
 * @param[in] OpencVPointCloud Matrix of OpenCV cv::Mat format
 * @param[out] PointCloud of pcl::PointCloud<pcl::PointXYZ> format
 */
//pcl::PointCloud<pcl::PointXYZ>::Ptr MatToPoinXYZ(cv::Mat OpencVPointCloud) {
//	pcl::PointCloud<pcl::PointXYZ>::Ptr pcXYZ(new pcl::PointCloud<pcl::PointXYZ>);
//
//	for (int i = 0; i < OpencVPointCloud.cols; i++) {
//		pcl::PointXYZ point;
//		point.x = OpencVPointCloud.at<float>(0, i);
//		point.y = OpencVPointCloud.at<float>(1, i);
//		point.z = OpencVPointCloud.at<float>(2, i);
//		pcXYZ->points.push_back(point);
//	}
//	pcXYZ->width = (int) pcXYZ->points.size();
//	pcXYZ->height = 1;
//
//	return pcXYZ;
//}

//cv::Mat PoinXYZToMat(pcl::PointCloud<pcl::PointXYZ>::ConstPtr &pcXYZ) {
//	cv::Mat OpenCVPointCloud(3, pcXYZ.size(), CV_64FC1);
//	for (int i = 0; i < point_cloud_ptr->points.size(); i++) {
//		OpenCVPointCloud.at<double>(0, i) = pcXYZ->points.at(i).x;
//		OpenCVPointCloud.at<double>(1, i) = pcXYZ->points.at(i).y;
//		OpenCVPointCloud.at<double>(2, i) = pcXYZ->points.at(i).z;
//	}
//
//	return OpenCVPointCloud;
//}

/***
 * Calculates rotation matrix to Euler angles
 */
Eigen::Vector3d convertRotationMatrixToEulerAngles(const cpd::Matrix &R) {
	pcl::console::TicToc tt;
	tt.tic();

	assert(isRotationMatrix(R));

	pcl::console::print_highlight(stderr, "Converting to Euler angles ");
	double sy = sqrt((double) (R(0, 0) * R(0, 0) + R(1, 0) * R(1, 0)));
	bool singular = sy < std::numeric_limits<float>::epsilon();

	double x, y, z;
	if (!singular) {
		x = std::atan2((double) R(2, 1), (double)R(2, 2));
		y = std::atan2((double)-R(2, 0), sy);
		z = std::atan2((double) R(1, 0), (double)R(0, 0));
	}
	else {
		x = std::atan2((double)-R(1, 2), (double)R(1, 1));
		y = std::atan2((double)-R(2, 0), sy);
		z = 0;
	}

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms : Euler angles: ");
	pcl::console::print_value("%f, %f, %f", x, y, z);
	pcl::console::print_info(" ]\n");

	return Eigen::Vector3d(x, y, z);
}

float computeDistanceAB(pcl::PointCloud<pcl::PointXYZ> &pcA,
												pcl::PointCloud<pcl::PointXYZ> &pcB) {
	float max_dist = -std::numeric_limits<float>::max();
	pcl::search::KdTree<pcl::PointXYZ> tree_b;
	tree_b.setInputCloud(pcB.makeShared());
	for (size_t i = 0; i < pcA.points.size(); ++i) {
		std::vector<int> indices(1);
		std::vector<float> sqr_distances(1);
		tree_b.nearestKSearch(pcA.points[i], 1, indices, sqr_distances);
		if (sqr_distances[0] > max_dist) {
			max_dist = sqr_distances[0];
		}
	}

	return std::sqrt(max_dist);
}

float computeHausdorffDistance(pcl::PointCloud<pcl::PointXYZ> &pcA,
															 pcl::PointCloud<pcl::PointXYZ> &pcB) {
	// Estimate
	pcl::console::TicToc tt;
	tt.tic();

	pcl::console::print_highlight(stderr, "Computing Hausdorff distance ");
	float max_dist_a = computeDistanceAB(pcA, pcB);  // compare A to B
	float max_dist_b = computeDistanceAB(pcB, pcA);  // compare B to A
	float dist = std::max(max_dist_a, max_dist_b);

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms : ");
	pcl::console::print_info("A->B: ");
	pcl::console::print_value("%f", max_dist_a);
	pcl::console::print_info(", B->A: ");
	pcl::console::print_value("%f", max_dist_b);
	pcl::console::print_info(", Hausdorff Distance: ");
	pcl::console::print_value("%f", dist);
	pcl::console::print_info(" ]\n");

	return dist;
}

void printHelp(int, char **argv) {
	pcl::console::print_error("Syntax is: %s fixed.ply moving.ply outfile.ply <options>\n",
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

void downsampleCloud(const pcl::PointCloud<pcl::PointXYZ>::ConstPtr &input,
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
}

bool loadCloud(const std::string &filename,
							 pcl::PointCloud<pcl::PointXYZ> &cloud) {
	pcl::console::TicToc tt;
	pcl::console::print_highlight("Loading ");
	pcl::console::print_value("%s ", filename.c_str());

	tt.tic();
	if (pcl::io::loadPLYFile(filename, cloud) < 0) {
		return (false);
	}

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms : ");
	pcl::console::print_value("%d", cloud.width * cloud.height);
	pcl::console::print_info(" points]\n");
	pcl::console::print_info("Available dimensions: ");
	pcl::console::print_value("%s\n", pcl::getFieldsList(cloud).c_str());

	return true;
}

void saveCloud(const std::string &filename,
							 const pcl::PointCloud<pcl::PointXYZ> &cloud) {
	pcl::console::TicToc tt;
	tt.tic();

	pcl::console::print_highlight("Saving ");
	pcl::console::print_value("%s ", filename.c_str());

	pcl::io::savePLYFileBinary(filename, cloud);

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms : ");
	pcl::console::print_value("%d", cloud.width * cloud.height);
	pcl::console::print_info(" points]\n");
}

void addCloudViz(pcl::visualization::PCLVisualizer &viewer,
								 pcl::PointCloud<pcl::PointXYZ>::ConstPtr pc,
								 const std::string& name,
								 const float r,
								 const float g, const float b) {
	viewer.addPointCloud<pcl::PointXYZ>(pc, name);
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, r, g,
																					b,
																					name);
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE,
																					1,
																					name);
}

bool convertMatrixToPointCloud(const cpd::Matrix &matrix,
															pcl::PointCloud<pcl::PointXYZ>::Ptr cloud) {
	try {
		pcl::PointXYZ point;
		for (unsigned int i = 0; i < matrix.rows(); i++) {
			point.x = matrix(i, 0);
			point.y = matrix(i, 1);
			point.z = matrix(i, 2);
			cloud->push_back(point);
		}
		cloud->width = (int) cloud->points.size();
		cloud->height = 1;

		return true;
	} catch (const std::exception& e) {
		std::cout << " a standard exception was caught, with message '" << e.what() << "'\n";
		return false;
	}
}

/***
 *
 */
void visResult(pcl::PointCloud<pcl::PointXYZ>::ConstPtr pcFixed,
							 pcl::PointCloud<pcl::PointXYZ>::ConstPtr pcMoving,
							 pcl::PointCloud<pcl::PointXYZ>::ConstPtr pcRegistered,
							 const cpd::Matrix &roi) {
	pcl::visualization::PCLVisualizer viewer("Marker Registration");

	addCloudViz(viewer, pcFixed, "Fixed", 1.0f, 1.0f, 1.0f);
	addCloudViz(viewer, pcMoving, "Moving", 1.0f, 1.0f, 0.0f);
	addCloudViz(viewer, pcRegistered, "Registered", 0.0f, 0.0f, 1.0f);

	// Display Doppler ROI as a 3D plan cutting through the center peaks of the registered marker
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcROI(new pcl::PointCloud<pcl::PointXYZ>);
	convertMatrixToPointCloud(roi, pcROI);
	viewer.addPolygon<pcl::PointXYZ>(pcROI, 0.0f, 0.0f, 1.0f, "DopplerROI");
	viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 6,
																		 "DopplerROI");

	viewer.setBackgroundColor(0, 0, 0);
	viewer.addCoordinateSystem(1.0, "first");
	viewer.setRepresentationToSurfaceForAllActors();
	viewer.initCameraParameters();
//	Clipping plane [near,far] 78.3397, 268.894
//	Focal point [x,y,z] 24.7589, 16.1316, 29.89
//	Position [x,y,z] -57.9468, -29.629, -102.325
//	View up [x,y,z] 0.172914, -0.959171, 0.223812
//	Camera view angle [degrees] 30
//	Window size [x,y] 918, 771
//	Window position [x,y] 277, 84
	while (!viewer.wasStopped()) {
		viewer.spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
}

int main(int argc, char** argv) {
	double execTime;
	clock_t clockStart;

	pcl::console::print_highlight("Parsing inputs\n");
	if (argc < 4) {
		printHelp(argc, argv);
		return -1;
	}

	// Parse the command line arguments for .ply files
	std::vector<int> idxFiles;
	idxFiles = pcl::console::parse_file_extension_argument(argc, argv, ".ply");
	if (idxFiles.size() != 3) {
		pcl::console::print_error("Need two input files and one output file.\n");
		return -1;
	}

	// Load the first file (fixed point cloud)
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcFixed(new pcl::PointCloud<pcl::PointXYZ>);
	if (!loadCloud(argv[idxFiles[0]], *pcFixed)) {
		return -1;
	}
	// Load the second file (moving point cloud)
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcMoving(new pcl::PointCloud<pcl::PointXYZ>);
	if (!loadCloud(argv[idxFiles[1]], *pcMoving)) {
		return -1;
	}

	// Parse command line options
	float sigma = default_sigma;
	pcl::console::parse_argument(argc, argv, "-sigma", sigma);

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

	bool performDownsample = ((leaf_x * leaf_y * leaf_z) == 0) ? false : true;
	pcl::PointCloud<pcl::PointXYZ> pcFixedDownsampled;
	if (performDownsample) {
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

		// Downsample using the voxel grid average
		downsampleCloud(pcFixed, pcFixedDownsampled, leaf_x, leaf_y, leaf_z, field, fmin,
										fmax);
	}
	else {
		pcFixedDownsampled = *pcFixed;
	}

	if (sigma == 0) {
		sigma = computeHausdorffDistance(*pcFixed, *pcMoving);
	}
	else if (sigma == -1) {
		sigma = computeHausdorffDistance(pcFixedDownsampled, *pcMoving);
	}
	else {
	}
	pcl::console::print_info("Using a sigma of: ");
	pcl::console::print_value("%f\n", sigma);

	// Get Eigen matrix
	cpd::Matrix ptsMoving = pcMoving->getMatrixXfMap().transpose();
	assert(ptsMoving.cols() >= 3);
	ptsMoving.conservativeResize(ptsMoving.rows(), 3);
	cpd::Matrix ptsFixed = pcFixedDownsampled.getMatrixXfMap().transpose();
	assert(ptsFixed.cols() >= 3);
	ptsFixed.conservativeResize(ptsFixed.rows(), 3);

	// Perform rigid registration
	pcl::console::TicToc tt;
	tt.tic();

	pcl::console::print_highlight("Performing registration\n");
	cpd::Rigid rigid;
	rigid.setSigma2(sigma);
	cpd::RigidResult result = rigid.run(ptsFixed, ptsMoving);

	pcl::console::print_info("[done, ");
	pcl::console::print_value("%g", tt.toc());
	pcl::console::print_info(" ms]\n");

	// Compute motor position (center of the marker)
	cpd::Matrix maxCoeff = result.points.colwise().maxCoeff(),
			minCoeff = result.points.colwise().minCoeff();
	cpd::Matrix aveCoeff = (maxCoeff + minCoeff) / 2.0f;
	float motorPosition = aveCoeff(2);
	float ratio = 2.4 / 4;	// channel
	float offset = (maxCoeff - minCoeff)(0) * ratio / 2;

	cpd::Matrix roi = cpd::Matrix::Constant(4, 3, motorPosition);
	roi.col(0) << aveCoeff(0) - offset, aveCoeff(0) + offset,
			aveCoeff(0) + offset, aveCoeff(0) - offset;
	roi.col(1) << minCoeff(1), minCoeff(1), maxCoeff(1), maxCoeff(1);

	std::cout << "Motor position: " << motorPosition << std::endl;
	std::cout << "Doppler ROI: (" << minCoeff(1) << ", "
	<< aveCoeff(0) - offset
	<< ", " << maxCoeff(1) << ", "
	<< aveCoeff(0) + offset
	<< ")" << std::endl;
	std::cout << "Doppler ROI: \n" << roi << std::endl;

	// Convert transformation matrix to Euler angles
	cpd::Matrix mtxTransform = result.matrix();
	Eigen::Vector3d angles = convertRotationMatrixToEulerAngles(mtxTransform);
	angles = angles * 180 / M_PI;

	pcl::console::print_info("Transformation matrix:\n");
	std::cout << mtxTransform << std::endl;
	pcl::console::print_info("Euler angles (deg.): ");
	pcl::console::print_value("%f, %f, %f\n", angles(0), angles(1), angles(2));

	// Save into the third file
//	std::cout << result.points.rows() << ", " << result.points.cols() << std::endl;
	assert(result.points.cols() == 3);
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcRegistered(new pcl::PointCloud<pcl::PointXYZ>);
	convertMatrixToPointCloud(result.points, pcRegistered);
	saveCloud(argv[idxFiles[2]], *pcRegistered);

	// viz initialization and result
	visResult(pcFixed, pcMoving, pcRegistered, roi);

	return 0;
}

// EOF
