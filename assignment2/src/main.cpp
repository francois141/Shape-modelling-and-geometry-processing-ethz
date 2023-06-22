#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>

#include <algorithm>
#include <iostream>
#include <chrono>
#include <chrono>
#include <chrono>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

#define LARGE_POSITIVE_VALUE 1000

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Normals evaluated via PCA method, #P x3
Eigen::MatrixXd NP;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int degree_polynomial = 1;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: Smooth parameter wendland function
double smooth_parameter = 100;

// Parameter: grid resolution
int resolutionX = 20;
int resolutionY = 20;
int resolutionZ = 20;

// Parameter: acceleration resolution
int acceleration_resolution = 10;

// Parameter: enlarge factor
double enlarge_grid_factor = 1.10;

// Parameter: epsilon factor
double epsilon_factor = 0.01;

// Parameter : Will compute optimal BB if > 0
int automatic_rotation = 0;

// Parameter : Output path of the mesh (to save it)
std::string output_path = "output.off";

// Parameter : Search with brute force if > 0
int brute_force = 0;

// Parameter : Check that the base constrained point is the nearest one
int check_closest_point = 0;


// Intermediate result: grid points, at which the implicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

// Intermediate result
vector<vector<int>> values_index;

// Functions
void createGrid();
void evaluateImplicitFunc();
void evaluateImplicitFunc_PolygonSoup();
void getLines();
void pcaNormal();
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);

inline double weightFunction(double radius) {
    return pow(1 - radius / smooth_parameter, 4) * (4 * radius / smooth_parameter + 1);
}

template <typename T>
T clip(const T& n, const T& lower, const T& upper) {
    return std::max(lower, std::min(n, upper));
}


// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid()
{
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Enlarge a bit the grid
    bb_min = enlarge_grid_factor * bb_min;
    bb_max = enlarge_grid_factor * bb_max;

    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolutionX - 1);
    const double dy = dim[1] / (double)(resolutionY - 1);
    const double dz = dim[2] / (double)(resolutionZ - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolutionX * resolutionY * resolutionZ, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolutionX; ++x)
    {
        for (unsigned int y = 0; y < resolutionY; ++y)
        {
            for (unsigned int z = 0; z < resolutionZ; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionX * (y + resolutionY * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
    // Reset the value of course
    grid_values.resize(pow(acceleration_resolution,3) );
    grid_values.setZero(pow(acceleration_resolution,3));
}

void rotateFigure() {
    if(automatic_rotation == 0) {
        return;
    }

    Eigen::MatrixXd points = P;

    Eigen::MatrixXd centered = points.rowwise() - points.colwise().mean();
    Eigen::MatrixXd covariance_matrix = (centered.adjoint().eval() * centered) / double(points.rows() - 1);

    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(covariance_matrix);
    Eigen::MatrixXd rotationMatrix = eigensolver.pseudoEigenvectors();

    P = P * rotationMatrix;
    N = N * rotationMatrix;
}
// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    // Sphere center
    Eigen::MatrixXd bb_min = P.colwise().minCoeff();
    Eigen::MatrixXd bb_max = P.colwise().maxCoeff();

    Eigen::RowVector3d dim = bb_max - bb_min;
    const double diagonal_bounding_box_size = (bb_max - bb_min).lpNorm<2>();
    double h_radius = wendlandRadius * diagonal_bounding_box_size;

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolutionX * resolutionY * resolutionZ);

    // Initialize system
    unsigned int size_col_a;

    switch(degree_polynomial) {
        case 0:
            size_col_a = 1;
            break;
        case 1:
            size_col_a = 4;
            break;
        case 2:
            size_col_a = 10;
            break;
        default:
            return;
    }

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolutionX; ++x)
    {
        cout << x <<  ":" << resolutionX << endl;
        for (unsigned int y = 0; y < resolutionY; ++y)
        {
            for (unsigned int z = 0; z < resolutionZ; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionX * (y + resolutionY * z);

                // Value at (x,y,z) = implicit function for the sphere
                Eigen::RowVector3d current_point = grid_points.row(index);

                // Find points
                vector<Eigen::RowVector3d> points(0);
                vector<double> values_c(0);
                vector<double> distances(0);

                if(brute_force) {
                    for(int i = 0; i < constrained_points.rows();i++) {
                        Eigen::RowVector3d constrained_point = constrained_points.row(i);
                        const double current_radius = (current_point - constrained_point).lpNorm<2>();
                        if(current_radius < h_radius) {
                            points.push_back(constrained_points.row(i));
                            values_c.push_back(constrained_values(i,0));
                            distances.push_back(current_radius);
                        }
                    }
                }
                // We use the acceleration grid
                else {
                    // Trick ==> We set wendland radius a ratio of the bounding box and then it is really easy
                    int indexX = (current_point(0) - bb_min(0)) / (bb_max(0) - bb_min(0)) * acceleration_resolution;
                    int indexY = (current_point(1) - bb_min(1)) / (bb_max(1) - bb_min(1)) * acceleration_resolution;
                    int indexZ = (current_point(2) - bb_min(2)) / (bb_max(2) - bb_min(2)) * acceleration_resolution;

                    for(int x_accelerated_grid = max(indexX-1,0); x_accelerated_grid <= min(indexX + 1, acceleration_resolution - 1); x_accelerated_grid++) {
                        for(int y_accelerated_grid = max(indexY-1,0); y_accelerated_grid <= min(indexY + 1, acceleration_resolution - 1); y_accelerated_grid++) {
                            for(int z_accelerated_grid = max(indexZ-1,0); z_accelerated_grid <= min(indexZ + 1, acceleration_resolution - 1); z_accelerated_grid++) {
                                const int acceleration_grid_idx =  x_accelerated_grid + acceleration_resolution * (y_accelerated_grid + acceleration_resolution * z_accelerated_grid);
                                for(int i = 0; i < values_index[acceleration_grid_idx].size(); i++) {
                                    const int currentIdx = values_index[acceleration_grid_idx][i];
                                    Eigen::RowVector3d constrained_point = constrained_points.row(currentIdx);
                                    const double current_radius = (current_point - constrained_point).lpNorm<2>();
                                    if(current_radius < h_radius) {
                                        points.push_back(constrained_points.row(currentIdx));
                                        values_c.push_back(constrained_values(currentIdx,0));
                                        distances.push_back(current_radius);
                                    }
                                }
                            }
                        }
                    }
                }

                Eigen::MatrixXd A = Eigen::MatrixXd(points.size(), size_col_a);
                Eigen::MatrixXd f = Eigen::MatrixXd(values_c.size(), 1);
                Eigen::MatrixXd w = Eigen::MatrixXd(points.size(), points.size());

                for(int i = 0; i < A.rows();i++) {
                    A(i,0) = 1.0;
                    if(degree_polynomial >= 1) {
                        A(i,1) = points[i].x();
                        A(i,2) = points[i].y();
                        A(i,3) = points[i].z();
                    }
                    if(degree_polynomial >= 2) {
                        A(i,4) = points[i].x()*points[i].x();
                        A(i,5) = points[i].y()*points[i].y();
                        A(i,6) = points[i].z()*points[i].z();

                        A(i,7) = points[i].x()*points[i].y();
                        A(i,8) = points[i].y()*points[i].z();
                        A(i,9) = points[i].x()*points[i].z();
                    }
                }

                for(int i = 0; i < points.size();i++) {
                    f(i,0) = values_c[i];
                }

                w.setZero();

                for(int i = 0; i < w.rows();i++) {
                    w(i,i) = weightFunction(distances[i]);
                }

                // Solve the equations
                if(A.rows() > 0) {

                    Eigen::MatrixXd A_w = w * A;
                    Eigen::MatrixXd f_w = w * f;

                    Eigen::RowVectorXd c;

                    c.resize(1, size_col_a);
                    c = A_w.colPivHouseholderQr().solve(f_w).transpose();

                    Eigen::RowVectorXd input_point_values = Eigen::RowVectorXd(1,size_col_a);

                    input_point_values(0) = 1;

                    if(degree_polynomial >= 1) {
                        input_point_values(1) = current_point(0);
                        input_point_values(2) = current_point(1);
                        input_point_values(3) = current_point(2);
                    }
                    if(degree_polynomial >= 2) {
                        input_point_values(4) = current_point(0) * current_point(0);
                        input_point_values(5) = current_point(1) * current_point(1);
                        input_point_values(6) = current_point(2) * current_point(2);

                        input_point_values(7) = current_point(0) * current_point(1);
                        input_point_values(8) = current_point(1) * current_point(2);
                        input_point_values(9) = current_point(0) * current_point(2);
                    }

                    grid_values[index] = input_point_values.dot(c);
                } else {
                    grid_values[index] = LARGE_POSITIVE_VALUE;
                }
            }
        }
    }

    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end - start;

    std::cout << "elapsed time for the normal generation: " << elapsed_seconds.count() << "s\n";
}


void evaluateImplicitFunc_PolygonSoup()
{
    // Sphere center
    Eigen::MatrixXd bb_min = P.colwise().minCoeff();
    Eigen::MatrixXd bb_max = P.colwise().maxCoeff();

    Eigen::RowVector3d dim = bb_max - bb_min;
    const double diagonal_bounding_box_size = (bb_max - bb_min).lpNorm<2>();
    double h_radius = wendlandRadius * diagonal_bounding_box_size;

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolutionX * resolutionY * resolutionZ);

    // Initialize system
    unsigned int size_col_a;

    switch(degree_polynomial) {
        case 0:
            size_col_a = 1;
            break;
        case 1:
            size_col_a = 4;
            break;
        case 2:
            size_col_a = 10;
            break;
        default:
            return;
    }

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolutionX; ++x)
    {
        cout << x <<  ":" << resolutionX << endl;
        for (unsigned int y = 0; y < resolutionY; ++y)
        {
            for (unsigned int z = 0; z < resolutionY; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionX * (y + resolutionY * z);

                // Value at (x,y,z) = implicit function for the sphere
                Eigen::RowVector3d current_point = grid_points.row(index);

                // Find points
                vector<Eigen::RowVector3d> points(0);
                vector<double> values_c(0);
                vector<double> distances(0);

                // Trick ==> We set wendland radius a ratio of the bounding box and then it is really easy
                int indexX = (current_point(0) - bb_min(0)) / (bb_max(0) - bb_min(0)) * acceleration_resolution;
                int indexY = (current_point(1) - bb_min(1)) / (bb_max(1) - bb_min(1)) * acceleration_resolution;
                int indexZ = (current_point(2) - bb_min(2)) / (bb_max(2) - bb_min(2)) * acceleration_resolution;

                for(int x_accelerated_grid = max(indexX-1,0); x_accelerated_grid <= min(indexX + 1, acceleration_resolution - 1); x_accelerated_grid++) {
                    for(int y_accelerated_grid = max(indexY-1,0); y_accelerated_grid <= min(indexY + 1, acceleration_resolution - 1); y_accelerated_grid++) {
                        for(int z_accelerated_grid = max(indexZ-1,0); z_accelerated_grid <= min(indexZ + 1, acceleration_resolution - 1); z_accelerated_grid++) {
                            const int acceleration_grid_idx =  x_accelerated_grid + acceleration_resolution * (y_accelerated_grid + acceleration_resolution * z_accelerated_grid);
                            for(int i = 0; i < values_index[acceleration_grid_idx].size(); i++) {
                                const int currentIdx = values_index[acceleration_grid_idx][i];

                                if(currentIdx % 3 != 0) {
                                    continue; // We only want O(N) points which lies on the 0 level set
                                }

                                Eigen::RowVector3d constrained_point = constrained_points.row(currentIdx);
                                const double current_radius = (current_point - constrained_point).lpNorm<2>();
                                if(current_radius < h_radius) {
                                    points.push_back(constrained_points.row(currentIdx));
                                    distances.push_back(current_radius);
                                    double current_value = (current_point - constrained_point).dot(N.row(currentIdx/3));
                                    values_c.push_back(current_value + constrained_values(currentIdx,0));
                                }
                            }
                        }
                    }
                }

                Eigen::MatrixXd A = Eigen::MatrixXd(points.size(), size_col_a);
                Eigen::MatrixXd f = Eigen::MatrixXd(values_c.size(), 1);
                Eigen::MatrixXd w = Eigen::MatrixXd(points.size(), points.size());

                for(int i = 0; i < A.rows();i++) {
                    A(i,0) = 1.0;
                    if(degree_polynomial >= 1) {
                        A(i,1) = points[i].x();
                        A(i,2) = points[i].y();
                        A(i,3) = points[i].z();
                    }
                    if(degree_polynomial >= 2) {
                        A(i,4) = points[i].x()*points[i].x();
                        A(i,5) = points[i].y()*points[i].y();
                        A(i,6) = points[i].z()*points[i].z();

                        A(i,7) = points[i].x()*points[i].y();
                        A(i,8) = points[i].y()*points[i].z();
                        A(i,9) = points[i].x()*points[i].z();
                    }
                }

                for(int i = 0; i < points.size();i++) {
                    f(i,0) = values_c[i];
                }

                w.setZero();

                for(int i = 0; i < w.rows();i++) {
                    w(i,i) = weightFunction(distances[i]);
                }

                // Solve the equations
                if(A.rows() > 0) {

                    Eigen::MatrixXd A_w = w * A;
                    Eigen::MatrixXd f_w = w * f;

                    Eigen::RowVectorXd c;

                    c.resize(1, size_col_a);
                    c = A_w.colPivHouseholderQr().solve(f_w).transpose();

                    Eigen::RowVectorXd input_point_values = Eigen::RowVectorXd(1,size_col_a);

                    input_point_values(0) = 1;

                    if(degree_polynomial >= 1) {
                        input_point_values(1) = current_point(0);
                        input_point_values(2) = current_point(1);
                        input_point_values(3) = current_point(2);
                    }
                    if(degree_polynomial >= 2) {
                        input_point_values(4) = current_point(0) * current_point(0);
                        input_point_values(5) = current_point(1) * current_point(1);
                        input_point_values(6) = current_point(2) * current_point(2);

                        input_point_values(7) = current_point(0) * current_point(1);
                        input_point_values(8) = current_point(1) * current_point(2);
                        input_point_values(9) = current_point(0) * current_point(2);
                    }

                    grid_values[index] = input_point_values.dot(c);
                } else {
                    grid_values[index] = LARGE_POSITIVE_VALUE;
                }
            }
        }
    }
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines()
{
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x < resolutionX; ++x)
    {
        for (unsigned int y = 0; y < resolutionY; ++y)
        {
            for (unsigned int z = 0; z < resolutionZ; ++z)
            {
                int index = x + resolutionX * (y + resolutionY * z);
                if (x < resolutionX - 1)
                {
                    int index1 = (x + 1) + y * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolutionY - 1)
                {
                    int index1 = x + (y + 1) * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolutionZ - 1)
                {
                    int index1 = x + y * resolutionX + (z + 1) * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

// Estimation of the normals via PCA.
void pcaNormal()
{
    NP = N;
    // Sphere center
    Eigen::MatrixXd bb_min = P.colwise().minCoeff();
    Eigen::MatrixXd bb_max = P.colwise().maxCoeff();

    Eigen::RowVector3d dim = bb_max - bb_min;
    const double diagonal_bounding_box_size = (bb_max - bb_min).lpNorm<2>();
    double h_radius = 0.2 * diagonal_bounding_box_size;

    for(int i = 0; i < P.rows();i++) {
        Eigen::MatrixXd points(0,P.cols());

        for(int j = 0; j < P.rows();j++) {
            double distance = (P.row(i) - P.row(j)).lpNorm<2>();
            if(distance < h_radius) {
                points.conservativeResize(points.rows()+1,points.cols());
                points.row(points.rows()-1) = P.row(j);
            }
        }

        Eigen::MatrixXd centered = points.rowwise() - points.colwise().mean();
        Eigen::MatrixXd covariance_matrix = (centered.adjoint().eval() * centered) / double(points.rows() - 1);

        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(covariance_matrix);

        int bestIndex = -1;
        double minValue = numeric_limits<double>::max();

        for(int k = 0; k < 3;k++) {
            double absEigenValue = abs(eigensolver.eigenvalues().row(k).value());
            if(absEigenValue < minValue) {
                minValue = absEigenValue;
                bestIndex = k;
            }
        }

        if(eigensolver.eigenvectors().col(bestIndex).real().dot(N.row(i).transpose()) < 0)  {
            NP.row(i) = -eigensolver.eigenvectors().col(bestIndex).real();
        } else {
            NP.row(i) = eigensolver.eigenvectors().col(bestIndex).real();
        }
    }
}




bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        // Show imported points
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }

    if (key == '2')
    {
        // Rotate figure first
        rotateFigure();
        // Show all constraints
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add your code for computing auxiliary constraint points here
        // Add code for displaying all points, as above

        std::chrono::time_point<std::chrono::system_clock> start, end, end_force;
        start = std::chrono::system_clock::now();

        // Step 1) Fix epsilon value
        Eigen::RowVector3d bounding_box_min_corner = P.colwise().minCoeff();
        Eigen::RowVector3d bounding_box_max_corner = P.colwise().maxCoeff();
        Eigen::RowVector3d bounding_box_differences = bounding_box_max_corner - bounding_box_min_corner;

        const double epsilon = epsilon_factor * bounding_box_differences.lpNorm<2>();

        // Step 2) Resize the data structures
        constrained_points.resize(3 * P.rows(), 3);
        constrained_values.resize(3 * P.rows(), 1);
        Eigen::MatrixXd constrained_colors = Eigen::MatrixXd(3*P.rows(),3);
        constrained_colors.setZero();

        // Step 3) Create the points
        auto store_value = [&](int idx, Eigen::RowVector3d point, double eps, Eigen::RowVector3d color) -> void {
            constrained_points.row(idx) = point;
            constrained_values(idx) = eps;
            constrained_colors.row(idx) = color;
        };

        const Eigen::RowVector3d COLOR_BLUE = Eigen::RowVector3d(0,0,1);
        const Eigen::RowVector3d COLOR_RED = Eigen::RowVector3d(1,0,0);
        const Eigen::RowVector3d COLOR_GREEN = Eigen::RowVector3d(0,1,0);

        for(int i = 0; i < P.rows();i++) {
            // Normals can be unormalized in the input
            N.row(i) = N.row(i).eval() / N.row(i).lpNorm<2>();
            store_value(3*i, P.row(i), 0, COLOR_BLUE);
            store_value(3*i+1, P.row(i) + epsilon*N.row(i), epsilon, COLOR_RED);
            store_value(3*i+2, P.row(i) - epsilon*N.row(i), -epsilon, COLOR_GREEN);
        }

        end_force = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds_force = end_force - start;

        std::cout << "elapsed time for the normal generation with brute force: " << elapsed_seconds_force.count() << "s\n";

        // Step 4) Fill the points in the grid
        values_index = vector<vector<int>>(pow(acceleration_resolution,3), vector<int>(0));
        for(int i = 0; i < constrained_points.rows();i++) {
            int indexX = (constrained_points(i,0) - bounding_box_min_corner(0)) / (bounding_box_max_corner(0) - bounding_box_min_corner(0)) * acceleration_resolution;
            int indexY = (constrained_points(i,1) - bounding_box_min_corner(1)) / (bounding_box_max_corner(1) - bounding_box_min_corner(1)) * acceleration_resolution;
            int indexZ = (constrained_points(i,2) - bounding_box_min_corner(2)) / (bounding_box_max_corner(2) - bounding_box_min_corner(2)) * acceleration_resolution;

            indexX = clip(indexX,0, acceleration_resolution-1);
            indexY = clip(indexY,0, acceleration_resolution-1);
            indexZ = clip(indexZ,0, acceleration_resolution-1);

            int index = indexX + acceleration_resolution * (indexY + acceleration_resolution * indexZ);

            assert(index >= 0);
            assert(index < acceleration_resolution*acceleration_resolution*acceleration_resolution);

            values_index[index].push_back(i);
        }

        end = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds = end - start;

        std::cout << "elapsed time for the normal generation: " << elapsed_seconds.count() << "s\n";

        // Step 5) Check that the nearest point of each constrained point is the nearest and warn user when it isn't the case

        if(check_closest_point > 0){

            Eigen::MatrixXd bb_min = P.colwise().minCoeff();
            Eigen::MatrixXd bb_max = P.colwise().maxCoeff();

            Eigen::RowVector3d dim = bb_max - bb_min;

            unsigned int count_not_nearest = 0;

            for(auto current_point_index = 0; current_point_index < constrained_points.rows();current_point_index++) {

                int best_index = -1;
                double best_norm = std::numeric_limits<double>::max();

                if (current_point_index % 3 == 0) {
                    continue; // We check only the generated constraints
                }

                auto current_point = constrained_points.row(current_point_index);

                int indexX = (current_point(0) - bb_min(0)) / (bb_max(0) - bb_min(0)) * acceleration_resolution;
                int indexY = (current_point(1) - bb_min(1)) / (bb_max(1) - bb_min(1)) * acceleration_resolution;
                int indexZ = (current_point(2) - bb_min(2)) / (bb_max(2) - bb_min(2)) * acceleration_resolution;

                for (int x_accelerated_grid = max(indexX - 1, 0);
                x_accelerated_grid <= min(indexX + 1, acceleration_resolution - 1); x_accelerated_grid++) {
                    for (int y_accelerated_grid = max(indexY - 1, 0);
                    y_accelerated_grid <= min(indexY + 1, acceleration_resolution - 1); y_accelerated_grid++) {
                        for (int z_accelerated_grid = max(indexZ - 1, 0); z_accelerated_grid <= min(indexZ + 1,
                                                                                                                acceleration_resolution -
                                                                                                                1); z_accelerated_grid++) {
                            const int acceleration_grid_idx = x_accelerated_grid + acceleration_resolution *
                                                                                               (y_accelerated_grid +
                                                                                                acceleration_resolution *
                                                                                                z_accelerated_grid);

                            for (int i = 0; i < values_index[acceleration_grid_idx].size(); i++) {
                                const int current_index = values_index[acceleration_grid_idx][i];
                                Eigen::RowVector3d constrained_point = constrained_points.row(current_index);
                                const double current_radius = (current_point - constrained_point).lpNorm<2>();
                                if (current_index != current_point_index && current_radius < best_norm) {
                                    best_norm = current_radius;
                                    best_index = current_index;
                                }
                            }
                        }
                    }
                }

                if(best_index != (current_point_index - (current_point_index % 3))) {
                    count_not_nearest++;
                }
            }

            cout << "Ratio of non-nearest point is : " << (double)count_not_nearest / (2.0*P.rows()) << endl;
            cout << "The absolute value is : " << count_not_nearest << endl;
        }

        viewer.data().point_size = 11;
        viewer.data().add_points(constrained_points, constrained_colors);
    }

    if (key == '3')
    {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/
    }

    if (key == '4')
    {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolutionX, resolutionY, resolutionZ, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);

        igl::writeOFF(output_path, V, F);
    }

    if (key == '5')
    {
        // Use the structure for key=='3' but replace the function evaluateImplicitFunc();
        // with a function performing the approximation of the implicit surface from polygon soup
        // Ref: Chen Shen, James F. O’Brien, and Jonathan Richard Shewchuk. Interpolating and approximating implicit surfaces from polygon soup.

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function --> Function to be modified here
        evaluateImplicitFunc_PolygonSoup();

        // get grid lines
        getLines();

        // Display the reconstruction
        callback_key_down(viewer, '4', modifiers);
    }

    if (key == '6' || key == '7' || key == '8')
    {
        pcaNormal();
        // To use the normals estimated via PCA instead of the input normals and then restaurate the input normals
        Eigen::MatrixXd N_tmp = N;
        N = NP;

        switch (key)
        {
        case '6':
            // Implement PCA Normal Estimation --> Function to be modified here
            pcaNormal();
            callback_key_down(viewer, '2', modifiers);
            break;
        case '7':
            callback_key_down(viewer, '3', modifiers);
            break;
        case '8':
            callback_key_down(viewer, '2', modifiers);
            callback_key_down(viewer, '3', modifiers);
            callback_key_down(viewer, '4', modifiers);
            break;
        default:
            break;
        }

        // Restore input normals
        N = N_tmp;
    }

    return true;
}

bool callback_load_mesh(Viewer &viewer, string filename)
{
    igl::readOFF(filename, P, F, N);
    callback_key_down(viewer, '1', 0);
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage ex2_bin <mesh.off>" << endl;
        igl::readOFF("../data/cat.off", P, F, N);
    }
    else
    {
        // Read points and normals
        igl::readOFF(argv[1], P, F, N);
    }

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            ImGui::InputInt("Resolution X", &resolutionX, 0, 0);
            ImGui::InputInt("Resolution Y", &resolutionY, 0, 0);
            ImGui::InputInt("Resolution Z", &resolutionZ, 0, 0);
            ImGui::InputInt("Resolution acceleration grid", &acceleration_resolution, 0, 0);
            ImGui::InputInt("Polynomial degree", &degree_polynomial, 0, 0);

            ImGui::InputDouble("Wendland Radius", &wendlandRadius, 0,0);
            ImGui::InputDouble("Wendland smooth parameter", &smooth_parameter, 0, 0);
            ImGui::InputDouble("Epsilon", &epsilon_factor);
            ImGui::InputDouble("Enlarge factor", &enlarge_grid_factor);
            ImGui::InputInt("Brute force", &brute_force);
            ImGui::InputInt("Check closest point",&check_closest_point);

            ImGui::InputInt("Automatic rotation", &automatic_rotation);
            ImGui::InputText("Ouput path", output_path);

            if (ImGui::Button("Reset Grid", ImVec2(-1, 0)))
            {
                std::cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer, '3', 0);
            }
        }
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
