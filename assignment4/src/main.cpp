#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <algorithm>

/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

// color per face
Eigen::MatrixXd C;

bool showingUV = false;
bool freeBoundary = false;
bool fastBoundary = false;
bool verticesFromBoundary = false;
bool outputCoutOfParametrization = false;

double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

enum distortionTypes {
    CONFORMAL_ANGLE_PRESEVING = 0,
    ISOMETRIC_LENGTH_PRESERVING = 1,
    AUTHALIC_AREA_PRESERVING = 2,
};


// Parameter: Choose with distortion type
distortionTypes distortionType = CONFORMAL_ANGLE_PRESEVING;

// Parameter: Distortion type
int distortionTypeInput = 0;

double white_color_threshold = 0.0;
double red_color_threshold = 1.0;

bool visualizeResult = false;

inline Eigen::RowVector3d lerp(const double& alpha, const Eigen::RowVector3d& lower, const Eigen::RowVector3d& upper) {
    return alpha * upper + (1 - alpha) * lower;
}

void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

    if(UV.size() != 0 && !visualizeResult)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}

    if(C.size() != 0 && visualizeResult)
    {
        viewer.data().set_colors(C);
    }
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;

}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}


static inline void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
    const unsigned int nbPositions = positions.rows();
    const unsigned int nbTotalPoints = V.rows();
    const unsigned int nbIndices = indices.rows();

    // Step 1: Fill C
    C.resize(2*nbPositions, 2*nbTotalPoints);
    vector<Eigen::Triplet<double>> cTriplets;
    for(int i = 0; i < indices.rows(); i++) {
        cTriplets.push_back(Eigen::Triplet<double>(i, indices[i], 1));
        cTriplets.push_back(Eigen::Triplet<double>(nbIndices + i, nbTotalPoints + indices[i], 1));
    }
    C.setFromTriplets(cTriplets.begin(),cTriplets.end());

    // step 2: Fill d
    d.setZero(2*nbPositions);
    for(int i = 0; i < nbPositions;i++) {
        d[i] = positions(i,0);
        d[nbPositions + i] = positions(i,1);
    }
}

static Eigen::SparseMatrix<double> GetAreaMatrix() {
    Eigen::VectorXd areas;
    igl::doublearea(V,F, areas);

    Eigen::SparseMatrix<double> areasDiagonal;
    areasDiagonal.resize(areas.rows(), areas.rows());

    vector<Eigen::Triplet<double>> areaTriplets;
    for(int i = 0; i < areas.size();i++) {
        areaTriplets.push_back(Eigen::Triplet<double>(i,i,0.5*areas(i)));
    }

    areasDiagonal.setFromTriplets(areaTriplets.begin(), areaTriplets.end());

    return areasDiagonal;
}

static void FillAWithL(Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &L) {
    vector<Eigen::Triplet<double>> aTriplets;
    for(int i = 0; i < L.outerSize();i++) {
        for(SparseMatrix<double>::InnerIterator it(L,i); it; ++it) {
            aTriplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
            aTriplets.push_back(Eigen::Triplet<double>(V.rows() + it.row(), V.rows() + it.col(), it.value()));
        }
    }
    A.setFromTriplets(aTriplets.begin(), aTriplets.end());
}

static inline Eigen::VectorXd SolveLinearSystem(Eigen::SparseMatrix<double> &A,Eigen::SparseMatrix<double> &C,Eigen::VectorXd &b,Eigen::VectorXd &d)
{
    // Setup left part
    Eigen::SparseMatrix<double> Zero;
    Zero.resize(C.rows(), C.rows());

    Eigen::SparseMatrix<double> C_t = C.transpose();

    Eigen::SparseMatrix<double> row1, row2, A_system;
    // Concatenate first row
    igl::cat(1, A, C, row1);
    // Concatenate second row
    igl::cat(1, C_t, Zero, row2);
    // Concatenate the matrices together
    igl::cat(2, row1, row2, A_system);

    // Setup right part
    Eigen::VectorXd b_system;
    b_system.resize(b.rows() + d.rows());
    b_system << b,d;

    // Solve the system
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
    solver.analyzePattern(A_system);
    solver.factorize(A_system);
    return solver.solve(b_system);
}


void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;

	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary) {
        // The boundary vertices should be fixed to positions on the unit disc. Find these position and
        // save them in the #V x 2 matrix fixed_UV_position.
        igl::boundary_loop(F, fixed_UV_indices);
        igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
    } else {
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
        vector<double> maxDistanceFromPoint;
        vector<int> maxDistanceFromPointIndex;

        vector<vector<int> > graph;
        igl::adjacency_list(F, graph);

        if(verticesFromBoundary){
            igl::boundary_loop(F, fixed_UV_indices);
            const unsigned int nbIndices = fixed_UV_indices.rows();
            for (int i = 0; i < nbIndices; i++) {
                VectorXd distances;
                VectorXi previous_vertex;
                igl::dijkstra(fixed_UV_indices(i,0), set<int>(), graph, distances, previous_vertex);

                Eigen::VectorXd::Index index;
                maxDistanceFromPoint.push_back(distances.maxCoeff(&index));
                maxDistanceFromPointIndex.push_back(index);
            }
        } else if(fastBoundary) {
            VectorXd distances;
            VectorXi previous_vertex;
            igl::dijkstra(0, set<int>(), graph, distances, previous_vertex);

            Eigen::VectorXd::Index index;
            maxDistanceFromPoint.push_back(distances.maxCoeff(&index));
            maxDistanceFromPointIndex.push_back(index);
        } else {
                const unsigned int nbVertices = V.rows();
                for (int i = 0; i < nbVertices; i++) {
                    VectorXd distances;
                    VectorXi previous_vertex;
                    igl::dijkstra(i, set<int>(), graph, distances, previous_vertex);

                    Eigen::VectorXd::Index index;
                    maxDistanceFromPoint.push_back(distances.maxCoeff(&index));
                    maxDistanceFromPointIndex.push_back(index);
                }
        }

        int index1 = max_element(maxDistanceFromPoint.begin(), maxDistanceFromPoint.end()) - maxDistanceFromPoint.begin();
        int index2 = maxDistanceFromPointIndex[index1];

        fixed_UV_indices.resize(2, 1);
        fixed_UV_indices << index1, index2;

        Eigen::RowVector2d point1, point2;
        point1 << 1,0;
        point2 << -1,0;

        fixed_UV_positions.resize(2,2);
        fixed_UV_positions.row(0) = point1;
        fixed_UV_positions.row(1) = point2;
	}

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

    const unsigned int nbPoints = V.rows();
    A.resize(2*nbPoints, 2*nbPoints);
    b.setZero(2 * nbPoints);

    Eigen::SparseMatrix<double> Dx,Dy, Dx_t, Dy_t;
    if(type != '1') {
        computeSurfaceGradientMatrix(Dx, Dy);
        Dx_t = Dx.transpose();
        Dy_t = Dy.transpose();
    }

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh
        Eigen::SparseMatrix<double> adjacencyList;
        igl::adjacency_matrix(F, adjacencyList);

        SparseVector<double> nbNeighbours;
        igl::sum(adjacencyList, 1, nbNeighbours);

        SparseMatrix<double> nbNeiboursDiagonal;
        igl::diag(nbNeighbours, nbNeiboursDiagonal);

        SparseMatrix<double> L = adjacencyList - nbNeiboursDiagonal;

        FillAWithL(A,L);
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~
        Eigen::SparseMatrix<double> areaMatrix = GetAreaMatrix();
        Eigen::SparseMatrix<double> L = (Dx_t * areaMatrix * Dx + Dy_t * areaMatrix * Dy);

        FillAWithL(A,L);
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
        Eigen::SparseMatrix<double> areaMatrix = GetAreaMatrix();

        Eigen::SparseMatrix<double> A11 = Dx_t * areaMatrix * Dx + Dy_t * areaMatrix * Dy;
        Eigen::SparseMatrix<double> A12 = -Dy_t * areaMatrix * Dx + Dx_t * areaMatrix * Dy;
        Eigen::SparseMatrix<double> A21 = Dy_t * areaMatrix * Dx - Dx_t * areaMatrix * Dy;
        Eigen::SparseMatrix<double> A22 = Dx_t * areaMatrix * Dx + Dy_t * areaMatrix * Dy;

        Eigen::SparseMatrix<double> A1, A2;
        igl::cat(1, A11, A12, A1);
        igl::cat(1, A21, A22, A2);

        igl::cat(2, A1, A2, A);
	}

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
        Eigen::VectorXd D11 = Dx * UV.col(0);
        Eigen::VectorXd D12 = Dy * UV.col(0);
        Eigen::VectorXd D21 = Dx * UV.col(1);
        Eigen::VectorXd D22 = Dy * UV.col(1);

        Eigen::SparseMatrix<double> areaMatrix = GetAreaMatrix();
        Eigen::SparseMatrix<double> L = (Dx_t * areaMatrix * Dx + Dy_t * areaMatrix * Dy);

        FillAWithL(A,L);

        const unsigned int nbFaces = F.rows();

        Eigen::MatrixXd RotationMatrix;
        RotationMatrix.resize(nbFaces,4);

        for(int i = 0; i < nbFaces;i++) {
            Eigen::Matrix2d jacobian;
            jacobian << D11[i], D12[i], D21[i], D22[i];

            Eigen::Matrix2d U,sigma,V;
            SSVD2x2(jacobian,U,sigma,V);

            Eigen::MatrixXd vt = V.transpose();

            Eigen::Matrix2d correctionMatrix;
            correctionMatrix << 1,0,0,(U * V.transpose()).determinant();

            Eigen::Matrix2d rotationMatrixFace = U * correctionMatrix * vt;

            RotationMatrix.row(i) << rotationMatrixFace(0,0), rotationMatrixFace(0,1), rotationMatrixFace(1,0) ,rotationMatrixFace(1,1);
        }

        Eigen::VectorXd ru = Dx_t*areaMatrix*RotationMatrix.col(0) + Dy_t * areaMatrix * RotationMatrix.col(1);
        Eigen::VectorXd rv = Dx_t*areaMatrix*RotationMatrix.col(2) + Dy_t * areaMatrix * RotationMatrix.col(3);

        igl::cat(1,ru, rv, b);
	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail
    Eigen::VectorXd x = SolveLinearSystem(A,C,b,d);

	// The solver will output a vector
	UV.resize(V.rows(), 2);
	UV.col(0) = x.block(0,0, V.rows(), 1);
	UV.col(1) = x.block(V.rows(),0, V.rows(), 1);

    if(outputCoutOfParametrization) {
        for(int i = 0; i < UV.rows();i++) {
            cout << UV(i,0) << " " << UV(i,1) << endl;
        }
    }
}

void computeDistortion() {
    Eigen::SparseMatrix<double> Dx,Dy;
    computeSurfaceGradientMatrix(Dx, Dy);

    Eigen::VectorXd D11 = Dx * UV.col(0);
    Eigen::VectorXd D12 = Dy * UV.col(0);
    Eigen::VectorXd D21 = Dx * UV.col(1);
    Eigen::VectorXd D22 = Dy * UV.col(1);

    const unsigned int nbFaces = F.rows();

    Eigen::VectorXd faceDistortion;
    faceDistortion.resize(nbFaces);

    for(int i = 0; i < nbFaces;i++) {

        Eigen::Matrix2d jacobian;
        jacobian << D11[i], D12[i], D21[i], D22[i];

        // We have to declare the variables outside the switch statement
        Eigen::Matrix2d jacobian_transpose;
        Eigen::Matrix2d identity;

        Eigen::Matrix2d U,sigma,V;
        Eigen::Matrix2d rotationMatrix;
        Eigen::Matrix2d correctionMatrix;

        double distortion;

        switch (distortionType) {
            case CONFORMAL_ANGLE_PRESEVING:
                jacobian_transpose = jacobian.transpose();
                identity.setIdentity();

                distortion = (jacobian + jacobian_transpose - identity * jacobian.trace()).squaredNorm();
                break;
            case ISOMETRIC_LENGTH_PRESERVING:
                SSVD2x2(jacobian, U, sigma, V);

                correctionMatrix << 1,0,0, (signbit((U * V.transpose()).determinant())) ? -1 : 1;
                rotationMatrix = U * correctionMatrix * V.transpose();

                distortion = (jacobian - rotationMatrix).squaredNorm();
                break;
            case AUTHALIC_AREA_PRESERVING:
                distortion = pow(jacobian.determinant() - 1,2);
                break;
        }

        // Finally save the distortion
        faceDistortion.row(i) << distortion;
    }

    // Normalize distortion
    assert(white_color_threshold <= red_color_threshold);

    double minFaceDistortion = faceDistortion.minCoeff();
    double maxFaceDistortion = faceDistortion.maxCoeff();
    for(int i = 0; i < nbFaces;i++) {
        faceDistortion(i) = (faceDistortion(i) - minFaceDistortion) / (maxFaceDistortion - minFaceDistortion);
        faceDistortion(i) = std::clamp(faceDistortion(i), white_color_threshold, red_color_threshold);
    }

    for(int i = 0; i < nbFaces;i++) {
        faceDistortion(i) = (faceDistortion(i) - white_color_threshold) / (red_color_threshold - white_color_threshold);
    }

    // Array to store the color of the faces
    C.resize(nbFaces, 3);

    // Save the actual distortion
    const Eigen::RowVector3d white = Eigen::RowVector3d(1,1,1);
    const Eigen::RowVector3d red = Eigen::RowVector3d(1,0,0);

    for(int i = 0; i < nbFaces;i++) {
        double blending = faceDistortion(i,0);
        C.row(i) << lerp(blending, white, red);
    }
}

void detectAndDisplayFlippedTriangles() {
    Eigen::SparseMatrix<double> Dx,Dy;
    computeSurfaceGradientMatrix(Dx, Dy);

    Eigen::VectorXd D11 = Dx * UV.col(0);
    Eigen::VectorXd D12 = Dy * UV.col(0);
    Eigen::VectorXd D21 = Dx * UV.col(1);
    Eigen::VectorXd D22 = Dy * UV.col(1);

    const unsigned int nbFaces = F.rows();

    Eigen::VectorXd faceDistortion;
    faceDistortion.resize(nbFaces);

    const Eigen::RowVector3d white = Eigen::RowVector3d(1,1,1);
    const Eigen::RowVector3d red = Eigen::RowVector3d(1,0,0);

    C.resize(nbFaces, 3);

    for(int i = 0; i < nbFaces;i++) {
        Eigen::Matrix2d jacobian;
        jacobian << D11[i], D12[i], D21[i], D22[i];

        double jacobianDeterminant = jacobian.determinant();
        if(jacobianDeterminant > 0) {
            C.row(i) << white;
        } else if(jacobianDeterminant < 0) {
            C.row(i) << red;
        } else {
            cout << "There is a determinant equal 0. This is not very good." << endl;
        }
    }
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
        visualizeResult = false;
		computeParameterization(key);
		break;
	case '5':
			// Add your code for detecting and displaying flipped triangles in the
			// UV domain here
        visualizeResult = true;
        detectAndDisplayFlippedTriangles();
		break;
    case '6':
        if(distortionTypeInput == 0) distortionType = CONFORMAL_ANGLE_PRESEVING;
        if(distortionTypeInput == 1) distortionType = ISOMETRIC_LENGTH_PRESERVING;
        if(distortionTypeInput == 2) distortionType = AUTHALIC_AREA_PRESERVING;
            visualizeResult = true;

        computeDistortion();
        break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core();
      viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core();
        viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core().align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
  if(argc != 2) {
    cout << "Usage ex4_bin <mesh.off/obj>" << endl;
    load_mesh("../data/cathead.obj");
  }
  else
  {
    // Read points and normals
    load_mesh(argv[1]);
  }

	igl::opengl::glfw::imgui::ImGuiPlugin plugin;
  viewer.plugins.push_back(&plugin);
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  plugin.widgets.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);
			ImGui::Checkbox("Fast boundary", &fastBoundary);
			ImGui::Checkbox("Print parametrization in cout", &outputCoutOfParametrization);
            ImGui::InputInt("Distortion type (1: Conformal, 2: Isometric, 3: Authalic", &distortionTypeInput);
            ImGui::Checkbox("Select fixed vertices from boundary", &verticesFromBoundary);
            ImGui::InputDouble("White value threshold", &white_color_threshold);
            ImGui::InputDouble("Red value threshold", &red_color_threshold);
		}
	};

  viewer.callback_key_pressed = callback_key_pressed;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_init = callback_init;

  viewer.launch();
}