#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>
#include <bits/stdc++.h>
#include <igl/readTGF.h>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/sum.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/adjacency_matrix.h>
#include <igl/normalize_row_sums.h>
#include <igl/column_to_quats.h>
#include <igl/dqs.h>
#include <igl/diag.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/grad.h>
#include <igl/polar_dec.h>
#include <chrono>
#include <cstdlib>

#include "util/Quaternion.h"

using namespace std;

// Data from the skeleton
Eigen::MatrixXd C;
Eigen::MatrixXi E;
Eigen::VectorXi P;
Eigen::MatrixXi BE;

// Transformed point for the skeleton
Eigen::MatrixXd TransformedPoint;

// Data from the mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

// Choose the type of input we get
enum MESH {
    HAND,
    BODY_VERTEX,
    BODY_FACE,
};

// Transformed point for the vertices
Eigen::MatrixXd TransformedVertices;

// Input rotations and quaternions
Eigen::MatrixXd rotations;
Eigen::MatrixXd quaternions;
int animationFrames;
int sizeRotationsInput;

// Used in the forward kinematics
vector<bool> visited;
vector<Eigen::MatrixXd> parentTransformation;

// Weights
Eigen::MatrixXd weights;

// Parameters for the computation of the weights
Eigen::VectorXi handles;
Eigen::VectorXi handlesFromFile;
bool load_weights_from_file = true;
int handles_to_select = 50;

// Parameters for the visualization of the weights
int weight_to_display = 0;
bool display_only_selected_handles = false;

// Visualization of the unposed images
bool visualize_unpose = false;

// Parameters for the animation
// Display options
bool display_skeleton = true;
bool display_weights = false;

// Rotations in the skeleton
bool interpolate_with_quaternions = false;

// Type of animation
bool linear_blend_skinning = true;
bool dual_quaternion_skinning = false;
bool per_face_skinning = false;

// Per face skinning
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> choleskySolver;
Eigen::SparseMatrix<double> Gradient;
Eigen::SparseMatrix<double> areaMatrix;
Eigen::SparseMatrix<double> A;
Eigen::SparseMatrix<double> A_free_to_free, A_fixed_to_free;
Eigen::VectorXi free_vertices, handle_vertices;

// Used to load poses in unpose
Eigen::MatrixXd currentPose;
Eigen::MatrixXi currentPoseFace;

// Unpose & displacements
Eigen::MatrixXd weightsPerPosePerDisplacement;
vector<vector<Eigen::MatrixXd>> referencePoses;
vector<Eigen::MatrixXd> displacementsPerPose;

// Visualization parameters
bool add_displacement = true;
double smooth_displacement = 0.05;

// Fixed indices for the bone
const unsigned int FIXED_BONE_HAND = 10;
const unsigned int FIXED_BONE_BODY = 0;

// Per face displacement
bool per_face_displacement = true;
bool rotation_and_skew = false;

int getEdgeIndex(int currentVertexIndex) {
    for (int i = 0; i < E.rows(); i++) {
        if (E(i, 1) == currentVertexIndex) {
            return i;
        }
    }

    cerr << "There is an error in the code. This section should not be reachable" << endl;
    exit(0);

    return -1;
}

int getParentIndex(int currentIdx) {
    for (int i = 0; i < E.rows(); i++) {
        if (E(i, 1) == currentIdx) {
            return E(i, 0);
        }
    }

    return -1;
}

vector<Eigen::MatrixXd> loadTransformations(int idx) {
    vector<Eigen::MatrixXd> transformationMatrices;
    // Add transformations
    const unsigned int start = idx * (sizeRotationsInput / animationFrames);
    for (int i = 0; i < E.rows(); i++) {
        transformationMatrices.push_back(rotations.block(start + 3 * i, 0, 3, 3));
    }
    return transformationMatrices;
}

void loadHand() {
    // Inline mesh of a cube
    bool ok = igl::readTGF("data/hand/hand.tgf", C, E);
    if (!ok) {
        cout << "Fail read TGF" << endl;
    }
    ok = igl::readDMAT("data/hand/hand-pose_matrot.dmat", rotations);
    if (!ok) {
        cout << "Fail read dmat" << endl;
    }
    ok = igl::readDMAT("data/hand/hand-pose_quat.dmat", quaternions);
    if (!ok) {
        cout << "Fail read dmat2" << endl;
    }
    ok = igl::readOFF("data/hand/hand.off", V, F);
    if (!ok) {
        cout << "Fail read off" << endl;
    }
    ok = igl::readDMAT("data/hand/hand-handles.dmat", handlesFromFile);
    if (!ok) {
        cout << "Fail read TGF" << endl;
    }

    // Safety assertion for the moment
    assert(C.rows() > 0 && C.cols() > 0);
    assert(E.rows() > 0 && E.cols() > 0);
    assert(rotations.rows() > 0 && rotations.cols() > 0);
    assert(quaternions.rows() > 0 && quaternions.cols() > 0);
    assert(V.rows() > 0 && V.cols() > 0);
    assert(F.rows() > 0 && F.cols() > 0);

    // Set size for the rotated point
    TransformedPoint.setZero(C.rows(), 3);

    // Resize parent transformation
    parentTransformation = vector<Eigen::MatrixXd>(C.rows());

    // Set number of animation frames
    animationFrames = 34;
    sizeRotationsInput = 2040;
}

void loadBody() {
    // Inline mesh of a cube
    bool ok = igl::readTGF("data/context-aware/skeleton.tgf", C, E);
    if (!ok) {
        cout << "Fail read TGF" << endl;
    }
    ok = igl::readDMAT("data/context-aware/all_frames.dmat", rotations);
    if (!ok) {
        cout << "Fail read dmat" << endl;
    }
    ok = igl::readOBJ("data/context-aware/reference.obj", V, F);
    if (!ok) {
        cout << "Fail read off" << endl;
    }
    ok = igl::readDMAT("data/context-aware/handles.dmat", handlesFromFile);
    if (!ok) {
        cout << "Fail read TGF" << endl;
    }

    // Safety assertion for the moment
    assert(C.rows() > 0 && C.cols() > 0);
    assert(E.rows() > 0 && E.cols() > 0);
    assert(rotations.rows() > 0 && rotations.cols() > 0);
    //assert(quaternions.rows() > 0 && quaternions.cols() > 0);
    assert(V.rows() > 0 && V.cols() > 0);
    assert(F.rows() > 0 && F.cols() > 0);

    // Set size for the rotated point
    TransformedPoint.setZero(C.rows(), 3);

    // Resize parent transformation
    parentTransformation = vector<Eigen::MatrixXd>(C.rows());

    // Set number of animation frames
    animationFrames = 394;
    sizeRotationsInput = 33096;
}

void displayWeights(igl::opengl::glfw::Viewer &viewer) {
    Eigen::MatrixXd colorVertices;
    colorVertices.setZero(V.rows(), 3);

    // We can display only the selected weights or display also the computed weights
    for (int i = 0; i < V.rows(); i++) {
        if (display_only_selected_handles && weights(i, weight_to_display) == 1.0) {
            colorVertices(i, 0) = weights(i, weight_to_display);
            colorVertices(i, 1) = weights(i, weight_to_display);
            colorVertices(i, 2) = weights(i, weight_to_display);
        } else if (!display_only_selected_handles) {
            colorVertices(i, 0) = weights(i, weight_to_display);
            colorVertices(i, 1) = weights(i, weight_to_display);
            colorVertices(i, 2) = weights(i, weight_to_display);
        }

    }

    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(colorVertices);
    viewer.data().show_lines = true;
    viewer.data().show_faces = true;
}

void displaySkeleton(igl::opengl::glfw::Viewer &viewer) {
    viewer.data().show_overlay_depth = false;
    viewer.data().set_points(TransformedPoint, Eigen::RowVector3d(1, 0.8, 0.2));
    viewer.data().set_edges(TransformedPoint, E, Eigen::RowVector3d(0.0,0.0,0.0));
}

void displayMesh(igl::opengl::glfw::Viewer &viewer) {
    viewer.data().set_mesh(TransformedVertices, F);
    viewer.data().show_lines = true;
    viewer.data().show_faces = true;
}

///////////////////////////////////////////////////

// Skeletal animation - 3

///////////////////////////////////////////////////

void rotateBone(int currentIndex, const vector<Eigen::MatrixXd> &localTransformations) {
    // Parent index
    const unsigned int parentIndex = getParentIndex(currentIndex);

    // Base case in the recursion
    if (parentIndex == -1) {
        parentTransformation[currentIndex] = Eigen::MatrixXd::Identity(4, 4);
        TransformedPoint.row(currentIndex) = C.row(currentIndex);
        return;
    }

    // Rotate parent first
    if (!visited[parentIndex]) {
        rotateBone(parentIndex, localTransformations);
    }
    visited[currentIndex] = true;

    // Compute our current point
    Eigen::Vector4d homogeneousPoint;
    homogeneousPoint << C.row(currentIndex).transpose(), 1;

    // Compute translation first
    Eigen::Matrix4d translationMatrix = Eigen::Matrix4d::Identity();
    translationMatrix.block(0, 3, 3, 1) = (-C.row(parentIndex)).transpose();

    // Compute rotation matrix
    Eigen::Matrix4d rotationMatrix = Eigen::Matrix4d::Identity();
    rotationMatrix.block(0, 0, 3, 3) = localTransformations[getEdgeIndex(currentIndex)];

    // Compute translation second
    Eigen::Matrix4d translationMatrixSecond = Eigen::Matrix4d::Identity();
    translationMatrixSecond.block(0, 3, 3, 1) = (C.row(parentIndex)).transpose();

    // Compute global transformation
    parentTransformation[currentIndex] =
            parentTransformation[parentIndex] * translationMatrixSecond * rotationMatrix * translationMatrix;

    // Update the point
    TransformedPoint.row(currentIndex) = (parentTransformation[currentIndex] * homogeneousPoint).block(0, 0, 3,
                                                                                                       1).transpose();
}

void getInterpolatedTransformation(std::chrono::time_point<std::chrono::system_clock> start) {
    // Get two indices
    auto current = chrono::system_clock::now();
    chrono::milliseconds s = chrono::duration_cast<chrono::milliseconds>(current - start);
    const unsigned int idxT1 = int(floor(s.count() / 1000.0));
    const unsigned int idxT2 = int(ceil(s.count() / 1000.0));
    double t = (idxT2 - (s.count() / 1000.0));

    // Get rotation matrices
    vector<Eigen::MatrixXd> firstTransformation = loadTransformations(idxT1);
    vector<Eigen::MatrixXd> secondTransformation = loadTransformations(idxT2);

    // Linearly interpolate the transformations
    vector<Eigen::MatrixXd> interpolatedTransformation(firstTransformation.size());
    for (int i = 0; i < firstTransformation.size(); i++) {
        interpolatedTransformation[i] = t * firstTransformation[i] + (1.0 - t) * secondTransformation[i];
    }

    // Draw the transformations
    visited = vector<bool>(C.rows(), false);
    for (int i = 0; i < C.rows(); i++) {
        rotateBone(i, interpolatedTransformation);
    }
}

void rotateQuaternion(std::chrono::time_point<std::chrono::system_clock> start) {
    // Get two indices
    auto current = chrono::system_clock::now();
    chrono::milliseconds s = chrono::duration_cast<chrono::milliseconds>(current - start);
    double t = s.count() / (animationFrames * 1000.0);

    // Absolute transformation matrices
    Eigen::MatrixXd T;

    // Create quaternions
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> dQ;
    igl::column_to_quats(quaternions, dQ);

    // Interpolate quaternions
    Eigen::Quaterniond idQuaternion(1, 0, 0, 0);
    for (int i = 0; i < E.rows(); i++) {
        dQ[i] = idQuaternion.slerp(t, dQ[i]);
    }

    // Perform forward kinematics
    igl::forward_kinematics(C, E, P, dQ, T);

    // Apply the transformations
    vector<Eigen::MatrixXd> transformedMatrices;
    for (int i = 0; i < E.rows(); i++) {
        Eigen::MatrixXd currentTransformation = Eigen::Matrix4d::Identity();
        currentTransformation.block(0, 0, 3, 4) = T.block(4 * i, 0, 4, 3).transpose();
        transformedMatrices.push_back(currentTransformation);
    }

    // Perform lbs
    for (int i = 0; i < C.rows(); i++) {
        Eigen::Vector4d homogeneousPoint;
        homogeneousPoint << C.row(i).transpose(), 1;
        // Update the parentTransformation
        if (getParentIndex(i) != -1) {
            Eigen::MatrixXd transformation = transformedMatrices[getEdgeIndex(i)] * homogeneousPoint;
            TransformedPoint.row(i) = transformation.block(0, 0, 3, 1).transpose();
            parentTransformation[i] = transformedMatrices[getEdgeIndex(i)];
        } else {
            TransformedPoint.row(i) = C.row(i);
            parentTransformation[i] = Eigen::Matrix4d::Identity();
        }
    }
}

///////////////////////////////////////////////////

// Weights computation - 4

///////////////////////////////////////////////////

Eigen::SparseMatrix<double> getLaplacianMatrix() {
    // Get adjacency matrix
    Eigen::SparseMatrix<double> adjacencyList;
    igl::adjacency_matrix(F, adjacencyList);

    // Valence degree of each vertex
    Eigen::SparseVector<double> nbNeighbours;
    igl::sum(adjacencyList, 1, nbNeighbours);

    // Create a diagonal matrix
    Eigen::SparseMatrix<double> nbNeiboursDiagonal;
    igl::diag(nbNeighbours, nbNeiboursDiagonal);

    // Finally build the laplacian
    Eigen::SparseMatrix<double> L = adjacencyList - nbNeiboursDiagonal;
    return L;
}

Eigen::MatrixXd ComputeWeights() {

    // Reset handles
    handles.setConstant(V.rows(), -1);

    // Choose handles
    for (int i = 0; i < E.rows(); i++) {
        Eigen::RowVector3d point = (C.row(E(i, 0)) + C.row(E(i, 1))) / 2;
        vector<pair<Eigen::RowVector3d, int>> distances;
        for (int j = 0; j < V.rows(); j++) {
            distances.push_back(make_pair(V.row(j), j));
        }

        const unsigned int nbNearest = min(handles_to_select, (int) distances.size());
        // Nth element is faster than sorting. In this case we only need a partial sort
        std::nth_element(distances.begin(), distances.begin() + nbNearest, distances.end(),
                         [&](const pair<Eigen::RowVector3d, int> &v1, const pair<Eigen::RowVector3d, int> &v2) -> bool {
                             return (v1.first - point).lpNorm<2>() < (v2.first - point).lpNorm<2>();
                         });
        // Set handles
        for (int j = 0; j < nbNearest; j++) {
            handles(distances[j].second) = i;
        }
    }

    // load from file
    if (load_weights_from_file) {
        handles = handlesFromFile;
    }

    // Get laplacian
    Eigen::SparseMatrix<double> A = getLaplacianMatrix();

    // Generate list of fixed and free indices
    int numFree = (handles.array() == -1).cast<int>().sum();
    int numFixed = V.rows() - numFree;

    // Create the handles
    Eigen::VectorXi handleFree;
    Eigen::VectorXi handleFixed;
    Eigen::MatrixXd weightsConstraints;

    handleFree.setZero(numFree);
    handleFixed.setZero(numFixed);
    weightsConstraints.setZero(V.rows(), E.rows());

    int count = 0;
    int countFree = 0;
    for (long vi = 0; vi < V.rows(); ++vi) {
        if (handles[vi] >= 0) {
            handleFixed[count++] = vi;
            weightsConstraints(vi, handles[vi]) = 1;
        } else {
            handleFree[countFree++] = vi;
        }
    }

    // Compute A
    Eigen::SparseMatrix<double> A_free_to_free;
    igl::slice(A, handleFree, handleFree, A_free_to_free);

    // Compute b
    Eigen::SparseMatrix<double> A_fixed_to_free;
    igl::slice(A, handleFree, handleFixed, A_fixed_to_free);
    Eigen::MatrixXd V_fixed = igl::slice(weightsConstraints, handleFixed, 1);
    Eigen::MatrixXd b = -A_fixed_to_free * V_fixed;

    // Solve the linear system
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> choleskySolver;
    choleskySolver.compute(A_free_to_free);
    Eigen::MatrixXd computedWeights = choleskySolver.solve(b);

    // Rescale the matrices
    igl::normalize_row_sums(computedWeights, computedWeights);

    // Add the weights back to the solution
    igl::slice_into(computedWeights, handleFree, 1, weightsConstraints);

    // Sanity check
    for (int i = 0; i < computedWeights.rows(); i++) {
        double sum = 0;
        for (int j = 0; j < computedWeights.cols(); j++) {
            sum += computedWeights(i, j);
        }
        assert(abs(sum - 1.0) <= 1e-3);
    }

    // Return the values
    return weightsConstraints;
}

///////////////////////////////////////////////////

// Linear blend skinning - 5

///////////////////////////////////////////////////

void animateLBS() {
    // Initialize transformed vertices
    TransformedVertices.setZero(V.rows(), 3);
    assert(weights.rows() == V.rows());

    // Apply linear skinning
    for (int i = 0; i < weights.rows(); i++) {
        // Create homogeneous point
        Eigen::Vector4d homogeneousPoint;
        homogeneousPoint << V.row(i).transpose(), 1;
        for (int j = 0; j < weights.cols(); j++) {
            // Transform the point
            Eigen::RowVector4d transformed = (parentTransformation[E(j, 1)] * homogeneousPoint).transpose();
            // Add it
            TransformedVertices.row(i) = TransformedVertices.row(i) + weights(i, j) * transformed.block(0, 0, 1, 3);
        }
    }
}

///////////////////////////////////////////////////

// Dual quaternion skinning - 6

///////////////////////////////////////////////////

void animateDqs() {
    // Initialize transformed vertices
    TransformedVertices.setZero(V.rows(), 3);
    assert(weights.rows() == V.rows());

    // Apply dqs
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> vQ;
    vector<Eigen::Vector3d> vT;

    for (int i = 0; i < E.rows(); i++) {
        // Get transformation matrix
        Eigen::Matrix4d currentTransformationMatrix = parentTransformation[E(i, 1)];
        // Decompose it
        Eigen::Quaterniond currentRotationQuaternion((Eigen::Matrix3d) currentTransformationMatrix.block(0, 0, 3, 3));
        Eigen::Vector3d currentTranslation(currentTransformationMatrix.block(0, 3, 3, 1));

        // Push it to the vectors
        vQ.push_back(currentRotationQuaternion);
        vT.push_back(currentTranslation);
    }

    // Finally run dqs
    igl::dqs(V, weights, vQ, vT, TransformedVertices);
}

///////////////////////////////////////////////////

// Animation per face - 7

///////////////////////////////////////////////////

void computeAreaMatrix(Eigen::SparseMatrix<double> &areaMatrix, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
    const unsigned int m = F.rows();
    Eigen::VectorXd areaVector;
    igl::doublearea(V, F, areaVector);

    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < F.rows(); i++) {
        triplets.emplace_back(i, i, 0.5 * areaVector(i));
        triplets.emplace_back(i + F.rows(), i + F.rows(), 0.5 * areaVector(i));
        triplets.emplace_back(i + 2 * F.rows(), i + 2 * F.rows(), 0.5 * areaVector(i));
    }

    areaMatrix.resize(3 * m, 3 * m);
    areaMatrix.setFromTriplets(triplets.begin(), triplets.end());
}

void computeGradientMatrix(Eigen::SparseMatrix<double> &G, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
    const unsigned int m = F.rows();

    Eigen::SparseMatrix<double> Gradient;
    igl::grad(V, F, Gradient);

    G.resize(3 * F.rows(), V.rows());
    G.reserve(Gradient.nonZeros());

    int idx = 0;
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < m; i++, idx += 3) {
        for (int j = 0; j < Gradient.cols(); ++j) {
            if (Gradient.coeff(i, j) != 0) triplets.emplace_back(idx, j, Gradient.coeff(i, j));
            if (Gradient.coeff(i + m, j) != 0) triplets.emplace_back(idx + 1, j, Gradient.coeff(i + m, j));
            if (Gradient.coeff(i + 2 * m, j) != 0) triplets.emplace_back(idx + 2, j, Gradient.coeff(i + 2 * m, j));
        }
    }

    G.setFromTriplets(triplets.begin(), triplets.end());
}

vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> getRotationsPerFace(vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> &dQ) {
    // Perform forward kinematics
    vector<Eigen::Vector3d> vT;
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> vQ;
    igl::forward_kinematics(C, E, P, dQ, vQ, vT);

    // Compute weights
    Eigen::MatrixXd h;
    h.setZero(F.rows(), weights.cols());
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < weights.cols(); j++) {
            double w = 0.0;
            w += weights(F(i, 0), j);
            w += weights(F(i, 1), j);
            w += weights(F(i, 2), j);
            h(i, j) = w / 3.0;
        }
    }

    // Perform rotation using SLERP
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> qF;
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Quaterniond quat(0, 0, 0, 0);
        for (int j = 0; j < h.cols(); j++) {
            quat = add(quat, mult(h(i, j), log(vQ[j])));
        }
        qF.push_back(exp(quat));
    }

    return qF;
}

void precomputePerFace(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const unsigned int FIXED_INDEX) {

    // Compute gradient matrix
    computeGradientMatrix(Gradient, V, F);

    // Compute area matrix
    computeAreaMatrix(areaMatrix, V, F);

    // Perform substitution
    vector<int> free, fixed;
    for (int i = 0; i < weights.rows(); i++) {
        if (weights(i, FIXED_INDEX) == 1.0) {
            fixed.push_back(i);
        } else {
            free.push_back(i);
        }
    }

    // Convert in eigen vector
    free_vertices.setZero(free.size());
    handle_vertices.setZero(fixed.size());
    for (int i = 0; i < free.size(); i++)
        free_vertices(i, 0) = free[i];
    for (int i = 0; i < handle_vertices.size(); i++)
        handle_vertices(i, 0) = fixed[i];

    // Create matrix a
    A = Gradient.transpose() * areaMatrix * Gradient;
    igl::slice(A, free_vertices, free_vertices, A_free_to_free);
    igl::slice(A, free_vertices, handle_vertices, A_fixed_to_free);

    // Solve system
    choleskySolver.compute(A_free_to_free);
}

void animatePerFace(std::chrono::time_point<std::chrono::system_clock> start) {

    // Get timeframe
    auto current = chrono::system_clock::now();
    chrono::milliseconds s = chrono::duration_cast<chrono::milliseconds>(current - start);
    double t = s.count() / (animationFrames * 1000.0);

    // Create a quaternion vector
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> dQ;
    igl::column_to_quats(quaternions, dQ);

    // Interpolate quaternions
    Eigen::Quaterniond idQuaternion(1, 0, 0, 0);
    for (int i = 0; i < E.rows(); i++) {
        dQ[i] = idQuaternion.slerp(t, dQ[i]);
    }

    // Get rotations
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> rotationCurrent = getRotationsPerFace(dQ);

    // Compute right side
    Eigen::MatrixXd rotatedPoints;
    rotatedPoints.setZero(3 * F.rows(), 3);
    for (int i = 0; i < F.rows(); i++) {
        Eigen::RowVector3d trans11 = V.row(F(i, 0));
        Eigen::RowVector3d trans12 = V.row(F(i, 1));
        Eigen::RowVector3d trans13 = V.row(F(i, 2));
        Eigen::RowVector3d trans21 = rotationCurrent[i] * V.row(F(i, 0));
        Eigen::RowVector3d trans22 = rotationCurrent[i] * V.row(F(i, 1));
        Eigen::RowVector3d trans23 = rotationCurrent[i] * V.row(F(i, 2));

        Eigen::MatrixXd m1, m2;

        m1.setZero(3, 3);
        m2.setZero(3, 3);

        m1.col(0) = trans21 - trans23;
        m1.col(1) = trans22 - trans23;
        m1.col(2) = trans21.cross(trans22);

        m2.col(0) = trans11 - trans13;
        m2.col(1) = trans12 - trans13;
        m2.col(2) = trans11.cross(trans12);

        rotatedPoints.block(3 * i, 0, 3, 3) = (m1 * m2.inverse()).transpose();
    }

    // Compute b
    Eigen::MatrixXd transformedPoints = Gradient.transpose() * areaMatrix * rotatedPoints;

    // Create matrix b
    Eigen::MatrixXd rotatedFixed = igl::slice(transformedPoints, free_vertices, 1);
    Eigen::MatrixXd V_fixed = igl::slice(V, handle_vertices, 1);
    Eigen::MatrixXd b_system = rotatedFixed - A_fixed_to_free * V_fixed;

    // Solve
    Eigen::MatrixXd solution = choleskySolver.solve(b_system);

    // Store solution
    TransformedVertices = V;
    igl::slice_into(solution, free_vertices, 1, TransformedVertices);
}

///////////////////////////////////////////////////

// Unpose per vertex - 8.1

///////////////////////////////////////////////////

double RBF(const vector<Eigen::MatrixXd> &pose1, const vector<Eigen::MatrixXd> &pose2) {
    assert(pose1.size() == pose2.size());
    double distance = 0.0;

    // Add each rotation individually
    // I also tried with translation
    for(int i = 0; i < pose1.size();i++) {
        distance += (pose1[i].block(0,0,3,3) - pose2[i].block(0,0,3,3)).norm();
    }

    // Return the exponential of the distance
    return exp(-smooth_displacement*distance*distance);
}

void getUnposeTransformation(vector<Eigen::MatrixXd> &transformations) {
    // Draw the transformations
    visited = vector<bool>(C.rows(), false);
    for (int i = 0; i < C.rows(); i++) {
        rotateBone(i, transformations);
    }
}

void precomputeUnposedExamples() {
    // Clear data structures
    referencePoses = vector<vector<Eigen::MatrixXd>>(0);
    displacementsPerPose = vector<Eigen::MatrixXd>(0);

    // Load and unpose all examples
    for (int currentIdx = 0; currentIdx < 4; currentIdx++) {
        // Load
        bool ok = igl::readOBJ("data/context-aware/eg" + to_string(currentIdx) + ".obj", currentPose, currentPoseFace);
        if (!ok) {
            cout << "Fail read obj" << endl;
        }
        Eigen::MatrixXd currentPoseRotation;
        ok = igl::readDMAT("data/context-aware/eg" + to_string(currentIdx) + ".dmat", currentPoseRotation);
        if (!ok) {
            cout << "Fail read off" << endl;
        }

        // Visualization of the unpose
        if (visualize_unpose) {
            igl::opengl::glfw::Viewer viewer;
            viewer.data().set_mesh(currentPose, currentPoseFace);
            viewer.data().show_lines = true;
            viewer.data().show_faces = true;
            viewer.launch();
        }

        // Add unposeTransformations
        vector<Eigen::MatrixXd> unposeTransformation;
        for (int i = 0; i < E.rows(); i++) {
            unposeTransformation.push_back(currentPoseRotation.block(3 * i, 0, 3, 3));
        }

        // Get transformation
        getUnposeTransformation(unposeTransformation);

        // Perform the unpose
        Eigen::MatrixXd trans;
        trans.setZero(V.rows(), 3);
        for (int i = 0; i < weights.rows(); i++) {
            // Create homogeneous point
            Eigen::Vector4d homogeneousPoint;
            homogeneousPoint << currentPose.row(i).transpose(), 1;

            Eigen::Matrix4d unpose;
            unpose.setZero();
            for (int j = 0; j < weights.cols(); j++) {
                unpose += weights(i, j) * parentTransformation[E(j, 1)];
            }

            Eigen::RowVector4d transformed = (unpose.inverse() * homogeneousPoint).transpose();
            trans.row(i) = transformed.block(0, 0, 1, 3);
        }

        // Add parentTransformation to the reference pose
        referencePoses.push_back(parentTransformation);

        // Visualization of the unpose
        if (visualize_unpose) {
            igl::opengl::glfw::Viewer viewer;
            viewer.data().set_mesh(trans, F);
            viewer.data().show_lines = true;
            viewer.data().show_faces = true;
            viewer.launch();
        }

        // Get displacement
        Eigen::MatrixXd currentDisplacement;
        currentDisplacement.setZero(V.rows(), 3);
        for (int i = 0; i < V.rows(); i++) {
            currentDisplacement.row(i) = trans.row(i) - V.row(i);
        }
        displacementsPerPose.push_back(currentDisplacement);
    }

    // Second step compute parameters for the displacements
    Eigen::MatrixXd A;
    A.setZero(referencePoses.size(), referencePoses.size());

    for(int i = 0; i < A.rows();i++) {
        for(int j = 0; j < A.cols();j++) {
            A(i,j) = RBF(referencePoses[i], referencePoses[j]);
        }
    }

    Eigen::MatrixXd b;
    b.setZero(referencePoses.size(), referencePoses.size());
    b.setIdentity();

    weightsPerPosePerDisplacement = A.inverse() * b;
}

void displayBodyPose(igl::opengl::glfw::Viewer &viewer, std::chrono::time_point<std::chrono::system_clock> start) {

    // Get two indices
    auto current = chrono::system_clock::now();
    chrono::milliseconds s = chrono::duration_cast<chrono::milliseconds>(current - start);
    const unsigned int idxT1 = int(floor(s.count() / 1000.0));
    const unsigned int idxT2 = int(ceil(s.count() / 1000.0));
    double t = (idxT2 - (s.count() / 1000.0));

    // Get rotation matrices
    vector<Eigen::MatrixXd> firstTransformation = loadTransformations(idxT1);
    vector<Eigen::MatrixXd> secondTransformation = loadTransformations(idxT2);

    // Linearly interpolate the transformations
    vector<Eigen::MatrixXd> transformations(firstTransformation.size());
    for (int i = 0; i < firstTransformation.size(); i++) {
        transformations[i] = t * firstTransformation[i] + (1.0 - t) * secondTransformation[i];
    }

    // Create parent transformation
    visited = vector<bool>(C.rows(), false);
    for (int i = 0; i < C.rows(); i++) {
        rotateBone(i, transformations);
    }
    transformations = parentTransformation;

    // Prepare transformed point
    Eigen::MatrixXd trans;
    trans.setZero(V.rows(), 3);

    // Compute displacements
    Eigen::MatrixXd poseDistances;
    poseDistances.setZero(4, 4);
    for (int i = 0; i < referencePoses.size(); i++) {
        for (int j = 0; j < referencePoses.size(); j++) {
            poseDistances(i, j) = RBF(transformations, referencePoses[j]);
        }
    }

    // Create weights
    Eigen::VectorXd poseWeights;
    poseWeights.setZero(referencePoses.size());
    for (int i = 0; i < referencePoses.size(); i++) {
        for (int j = 0; j < referencePoses.size(); j++) {
            poseWeights(i) += weightsPerPosePerDisplacement(i, j) * RBF(transformations, referencePoses[j]);
        }
    }

    // Apply partition of unity on the weights
    double sum = 0.0;
    for (int i = 0; i < referencePoses.size(); i++) {
        sum += poseWeights(i, 0);
    }
    for (int i = 0; i < referencePoses.size(); i++) {
        poseWeights(i, 0) /= sum;
    }

    // Compute the actual transformation with the displacements
    for (int i = 0; i < V.rows(); i++) {

        // Compute displacements
        Eigen::Vector3d currentDisplacement;
        currentDisplacement.setZero();
        for (int j = 0; j < referencePoses.size(); j++) {
            currentDisplacement += poseWeights(j) * displacementsPerPose[j].row(i);
        }

        // Create homogeneous point
        Eigen::Vector4d homogeneousPoint;
        homogeneousPoint << V.row(i).transpose(), 1;

        // Add the displacement
        if (add_displacement) {
            homogeneousPoint(0) += currentDisplacement(0);
            homogeneousPoint(1) += currentDisplacement(1);
            homogeneousPoint(2) += currentDisplacement(2);
        }

        // Perform LBS
        for (int j = 0; j < weights.cols(); j++) {
            // Transform the point
            Eigen::RowVector4d transformed = (transformations[E(j, 1)] * homogeneousPoint).transpose();
            // Add it
            trans.row(i) = trans.row(i) + weights(i, j) * transformed.block(0, 0, 1, 3);
        }
    }

    // Display mesh
    viewer.data().set_mesh(trans, F);
    viewer.data().show_lines = true;
    viewer.data().show_faces = true;
}

///////////////////////////////////////////////////

// Unpose per face - 8.2

///////////////////////////////////////////////////

void precomputeUnposedExamplesPerFace() {
    // Clear data structures
    referencePoses = vector<vector<Eigen::MatrixXd>>(0);
    displacementsPerPose = vector<Eigen::MatrixXd>(0);

    // Load them again
    for (int currentIdx = 0; currentIdx < 4; currentIdx++) {

        // One iteration takes time, we print to ack the progress
        cout << "Current iteration : [" << currentIdx + 1 << ":4]" << endl;

        // Load
        bool ok = igl::readOBJ("data/context-aware/eg" + to_string(currentIdx) + ".obj", currentPose, currentPoseFace);
        if (!ok) {
            cout << "Fail read obj" << endl;
        }
        Eigen::MatrixXd currentPoseRotation;
        ok = igl::readDMAT("data/context-aware/eg" + to_string(currentIdx) + ".dmat", currentPoseRotation);
        if (!ok) {
            cout << "Fail read off" << endl;
        }

        // Precompute
        precomputePerFace(currentPose, currentPoseFace, FIXED_BONE_BODY);

        // Visualization of the unpose
        if (visualize_unpose) {
            igl::opengl::glfw::Viewer viewer;
            viewer.data().set_mesh(currentPose, currentPoseFace);
            viewer.data().show_lines = true;
            viewer.data().show_faces = true;
            viewer.launch();
        }

        // Add unposeTransformations
        vector<Eigen::MatrixXd> unposeTransformation;
        for (int i = 0; i < E.rows(); i++) {
            unposeTransformation.push_back(currentPoseRotation.block(3 * i, 0, 3, 3));
        }

        // Get transformation
        getUnposeTransformation(unposeTransformation);

        vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> dQ;
        for(int i = 0; i < unposeTransformation.size();i++) {
            Eigen::Matrix3d currentRotation = unposeTransformation[i].block(0,0,3,3);
            dQ.push_back(Eigen::Quaterniond(currentRotation));
        }

        // Get rotations
        vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> rotationCurrent = getRotationsPerFace(dQ);

        // Perform the unpose
        Eigen::MatrixXd trans;
        trans.setZero(V.rows(), 3);

        // Compute right side
        Eigen::MatrixXd rotatedPoints;
        rotatedPoints.setZero(3 * F.rows(), 3);
        for (int i = 0; i < F.rows(); i++) {
            Eigen::RowVector3d trans11 = currentPose.row(currentPoseFace(i, 0));
            Eigen::RowVector3d trans12 = currentPose.row(currentPoseFace(i, 1));
            Eigen::RowVector3d trans13 = currentPose.row(currentPoseFace(i, 2));
            Eigen::RowVector3d trans21 = rotationCurrent[i].inverse() * currentPose.row(currentPoseFace(i, 0));
            Eigen::RowVector3d trans22 = rotationCurrent[i].inverse() * currentPose.row(currentPoseFace(i, 1));
            Eigen::RowVector3d trans23 = rotationCurrent[i].inverse() * currentPose.row(currentPoseFace(i, 2));

            Eigen::MatrixXd m1, m2;

            m1.setZero(3, 3);
            m2.setZero(3, 3);

            m1.col(0) = trans21 - trans23;
            m1.col(1) = trans22 - trans23;
            m1.col(2) = trans21.cross(trans22);

            m2.col(0) = trans11 - trans13;
            m2.col(1) = trans12 - trans13;
            m2.col(2) = trans11.cross(trans12);

            rotatedPoints.block(3 * i, 0, 3, 3) = (m1 * m2.inverse()).transpose();
        }

        // Compute b
        Eigen::MatrixXd transformedPoints = Gradient.transpose() * areaMatrix * rotatedPoints;

        // Create matrix b
        Eigen::MatrixXd rotatedFixed = igl::slice(transformedPoints, free_vertices, 1);
        Eigen::MatrixXd V_fixed = igl::slice(currentPose, handle_vertices, 1);
        Eigen::MatrixXd b_system = rotatedFixed - A_fixed_to_free * V_fixed;

        // Solve
        Eigen::MatrixXd solution = choleskySolver.solve(b_system);

        // Store solution
        trans = currentPose;
        igl::slice_into(solution, free_vertices, 1, trans);

        // Visualization of the unpose
        if (visualize_unpose) {
            igl::opengl::glfw::Viewer viewer;
            viewer.data().set_mesh(trans, F);
            viewer.data().show_lines = true;
            viewer.data().show_faces = true;
            viewer.launch();
        }

        // Add parentTransformation to the reference pose
        referencePoses.push_back(parentTransformation);

        // Get displacement
        Eigen::MatrixXd currentDisplacement;
        currentDisplacement.setZero(3 * F.rows(), 3);

        // Compute normals from
        Eigen::MatrixXd Nfrom;
        igl::per_face_normals(V, F, Nfrom);

        // Compute normals to
        Eigen::MatrixXd Nto;
        igl::per_face_normals(trans, currentPoseFace, Nto);

        // Compute displacement now
        for (int i = 0; i < F.rows(); i++) {
            Eigen::RowVector3d trans11 = V.row(F(i, 0));
            Eigen::RowVector3d trans12 = V.row(F(i, 1));
            Eigen::RowVector3d trans13 = V.row(F(i, 2));
            Eigen::RowVector3d trans21 = trans.row(currentPoseFace(i, 0));
            Eigen::RowVector3d trans22 = trans.row(currentPoseFace(i, 1));
            Eigen::RowVector3d trans23 = trans.row(currentPoseFace(i, 2));

            Eigen::MatrixXd m1, m2;

            m1.setZero(3, 3);
            m2.setZero(3, 3);

            m1.col(0) = trans21 - trans23;
            m1.col(1) = trans22 - trans23;
            m1.col(2) = Nto.row(i);

            m2.col(0) = trans11 - trans13;
            m2.col(1) = trans12 - trans13;
            m2.col(2)  = Nfrom.row(i);

            currentDisplacement.block(3 * i, 0, 3, 3) = m1 * m2.inverse();
        }
        displacementsPerPose.push_back(currentDisplacement);
    }

    // Second step compute parameters for the displacements
    Eigen::MatrixXd A;
    A.setZero(referencePoses.size(), referencePoses.size());

    for(int i = 0; i < A.rows();i++) {
        for(int j = 0; j < A.cols();j++) {
            A(i,j) = RBF(referencePoses[i], referencePoses[j]);
        }
    }

    Eigen::MatrixXd b;
    b.setZero(referencePoses.size(), referencePoses.size());
    b.setIdentity();

    weightsPerPosePerDisplacement = A.inverse() * b;

    // Finally we can precompute the base mesh
    precomputePerFace(V, F, FIXED_BONE_BODY);
}

void displayBodyPosePerFace(igl::opengl::glfw::Viewer &viewer, std::chrono::time_point<std::chrono::system_clock> start) {

    // Get two indices
    auto current = chrono::system_clock::now();
    chrono::milliseconds s = chrono::duration_cast<chrono::milliseconds>(current - start);
    const unsigned int idxT1 = int(floor(s.count() / 1000.0));
    const unsigned int idxT2 = int(ceil(s.count() / 1000.0));
    double t = (idxT2 - (s.count() / 1000.0));

    // Get rotation matrices
    vector<Eigen::MatrixXd> firstTransformation = loadTransformations(idxT1);
    vector<Eigen::MatrixXd> secondTransformation = loadTransformations(idxT2);

    // Create a quaternion vector
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> dQ;

    // Linearly interpolate the transformations and convert into a quaternion
    vector<Eigen::MatrixXd> transformations(firstTransformation.size());
    for (int i = 0; i < firstTransformation.size(); i++) {
        transformations[i] = (t * firstTransformation[i] + (1.0 - t) * secondTransformation[i]);
        Eigen::Matrix3d localRotation = transformations[i].block(0,0,3,3);
        dQ.push_back(Eigen::Quaterniond(localRotation));

    }

    // Create parent transformation
    visited = vector<bool>(C.rows(), false);
    for (int i = 0; i < C.rows(); i++) {
        rotateBone(i, transformations);
    }
    transformations = parentTransformation;

    // Compute displacements
    Eigen::MatrixXd poseDistances;
    poseDistances.setZero(4, 4);
    for (int i = 0; i < referencePoses.size(); i++) {
        for (int j = 0; j < referencePoses.size(); j++) {
            poseDistances(i, j) = RBF(transformations, referencePoses[j]);
        }
    }

    // Create weights
    Eigen::VectorXd poseWeights;
    poseWeights.setZero(referencePoses.size());
    for (int i = 0; i < referencePoses.size(); i++) {
        for (int j = 0; j < referencePoses.size(); j++) {
            poseWeights(i) += weightsPerPosePerDisplacement(i, j) * RBF(transformations, referencePoses[j]);
        }
    }

    // Apply partition of unity on the weights
    double sum = 0.0;
    for (int i = 0; i < referencePoses.size(); i++) {
        sum += poseWeights(i, 0);
    }
    for (int i = 0; i < referencePoses.size(); i++) {
        poseWeights(i, 0) /= sum;
    }

    // Get rotations
    vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>> rotationCurrent = getRotationsPerFace(dQ);

    // Compute right side
    Eigen::MatrixXd rotatedPoints;
    rotatedPoints.setZero(3 * F.rows(), 3);
    for (int i = 0; i < F.rows(); i++) {
        Eigen::RowVector3d trans11 = V.row(F(i, 0));
        Eigen::RowVector3d trans12 = V.row(F(i, 1));
        Eigen::RowVector3d trans13 = V.row(F(i, 2));
        Eigen::RowVector3d trans21 = rotationCurrent[i] * V.row(F(i, 0));
        Eigen::RowVector3d trans22 = rotationCurrent[i] * V.row(F(i, 1));
        Eigen::RowVector3d trans23 = rotationCurrent[i] * V.row(F(i, 2));

        Eigen::MatrixXd m1, m2;

        m1.setZero(3, 3);
        m2.setZero(3, 3);

        m1.col(0) = trans21 - trans23;
        m1.col(1) = trans22 - trans23;
        m1.col(2) = trans21.cross(trans22);

        m2.col(0) = trans11 - trans13;
        m2.col(1) = trans12 - trans13;
        m2.col(2) = trans11.cross(trans12);

        Eigen::MatrixXd currentTransformation = m1 * m2.inverse();

        if(per_face_displacement) {
            Eigen::Matrix3d displacementRotation;
            displacementRotation.setZero();

            Eigen::Quaterniond Qsum(0,0,0,0);
            Eigen::Matrix3d Ssum;
            Ssum.setZero();

            for(int j = 0; j < referencePoses.size();j++) {
                // Improve this pose then
                Eigen::Matrix3d currentDisplacement = displacementsPerPose[j].block(3*i,0,3,3);

                if(rotation_and_skew) {
                    // We perform the equation (10) here
                    Eigen::Matrix3d Q,S;
                    igl::polar_dec(currentDisplacement, Q, S);

                    Qsum = add(Qsum, mult(poseWeights(j), log(Eigen::Quaterniond(Q))));
                    Ssum += poseWeights(j) * S;
                } else {
                    // Simple decomposition - equation (9)
                    displacementRotation += poseWeights(j) * currentDisplacement;
                }
            }

            if(rotation_and_skew) {
                // Apply equation (10)
                currentTransformation = currentTransformation * Qsum * Ssum;
            } else {
                // Apply equation (9)
                currentTransformation = currentTransformation * displacementRotation;
            }
        }

        rotatedPoints.block(3 * i, 0, 3, 3) = currentTransformation.transpose();
    }

    // Compute b
    Eigen::MatrixXd transformedPoints = Gradient.transpose() * areaMatrix * rotatedPoints;

    // Create matrix b
    Eigen::MatrixXd rotatedFixed = igl::slice(transformedPoints, free_vertices, 1);
    Eigen::MatrixXd V_fixed = igl::slice(V, handle_vertices, 1);
    Eigen::MatrixXd b_system = rotatedFixed - A_fixed_to_free * V_fixed;

    // Solve
    Eigen::MatrixXd solution = choleskySolver.solve(b_system);

    // Store solution
    Eigen::MatrixXd trans = V;
    igl::slice_into(solution, free_vertices, 1, trans);

    // Display mesh
    viewer.data().set_mesh(trans, F);
    viewer.data().show_lines = true;
    viewer.data().show_faces = true;
}

int main(int argc, char *argv[]) {

    // Set current mesh to hand and change if we load the body
    MESH currentMesh;

    // Load
    if (argc > 1 && (strcmp(argv[1], "body") == 0)) {
        cout << "Load body with vertex displacement" << endl;
        currentMesh = BODY_VERTEX;
        loadBody();
    } else if(argc > 1 && (strcmp(argv[1], "face") == 0)) {
        cout << "Load body with face displacement" << endl;
        currentMesh = BODY_FACE;
        loadBody();
    } else {
        cout << "Load hand" << endl;
        currentMesh = HAND;
        loadHand();
    }

    // Create child - parent graph
    igl::directed_edge_parents(E, P);

    // Compute weights
    weights = ComputeWeights();

    // Viewer
    igl::opengl::glfw::Viewer viewer;
    viewer.core().is_animating = true;
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);

    // Callback functions
    function<bool(igl::opengl::glfw::Viewer& viewer)> animationCallbackFunction;
    function<void(void)> menuCallbackFunction;

    // Start chrono
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    // Display data
    if (currentMesh == HAND) {
        // Precompute for the per face skinning
        precomputePerFace(V, F, FIXED_BONE_HAND);

        // Bind callback functions
        animationCallbackFunction = [&](igl::opengl::glfw::Viewer &) {
            // Clear
            viewer.data().clear();

            // Check counter
            auto current = chrono::system_clock::now();
            chrono::seconds s = chrono::duration_cast<chrono::seconds>(current - start);
            if (s.count() == animationFrames) {
                start = std::chrono::system_clock::now();
            }

            // Interpolate
            vector<Eigen::MatrixXd> transformations;
            if (!interpolate_with_quaternions) {
                getInterpolatedTransformation(start);
            } else {
                rotateQuaternion(start);
            }

            // Draw
            if (display_skeleton) {
                displaySkeleton(viewer);
            } else if (display_weights) {
                displayWeights(viewer);
            } else {
                if (linear_blend_skinning) {
                    animateLBS();
                } else if (dual_quaternion_skinning) {
                    animateDqs();
                } else if (per_face_skinning) {
                    animatePerFace(start);
                } else {
                    cout << "nothing selected" << endl;
                }

                // Display it
                displayMesh(viewer);
            }

            return false;
        };

        menuCallbackFunction = [&]() {
            // Draw parent menu content
            menu.draw_viewer_menu();

            // Draw actual options
            if (ImGui::CollapsingHeader("Skeletal animation options", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::Checkbox("Interpolate with quaterions", &interpolate_with_quaternions);
            }
            if (ImGui::CollapsingHeader("Display options", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::Checkbox("Display skeleton", &display_skeleton);
                ImGui::Checkbox("Display weights", &display_weights);
                ImGui::InputInt("Weight to display", &weight_to_display);
                ImGui::Checkbox("Display only selected handles", &display_only_selected_handles);
            }
            if (ImGui::CollapsingHeader("Rotation weights options", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::Checkbox("Load weights from file", &load_weights_from_file);
                ImGui::InputInt("Number of handles to select per bone", &handles_to_select);
            }
            if (ImGui::CollapsingHeader("Animation options", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::Checkbox("Linear blend skinning", &linear_blend_skinning);
                ImGui::Checkbox("Dual quaternion skinning", &dual_quaternion_skinning);
                ImGui::Checkbox("Per face skinning", &per_face_skinning);
            }
            if (ImGui::CollapsingHeader("Animation", ImGuiTreeNodeFlags_DefaultOpen)) {
                if (ImGui::Button("Restart", ImVec2(-1, 0)))
                {
                    cout << "restart" << endl;

                    // Compute weights
                    weights = ComputeWeights();

                    // Precompute for the per face skinning
                    precomputePerFace(V,F,FIXED_BONE_HAND);

                    // Restart chrono
                    start = std::chrono::system_clock::now();
                }
            }
        };
    } else if(currentMesh == BODY_VERTEX) {
        // Precompute the unposed example
        precomputeUnposedExamples();

        // Bind callback functions
        animationCallbackFunction = [&](igl::opengl::glfw::Viewer &) {
            // Clear
            viewer.data().clear();

            // Check counter
            auto current = chrono::system_clock::now();
            chrono::seconds s = chrono::duration_cast<chrono::seconds>(current - start);
            if (s.count() == animationFrames) {
                start = std::chrono::system_clock::now();
            }

            // Perform animation
            displayBodyPose(viewer, start);
            return false;
        };

        menuCallbackFunction = [&]() {
            // Draw parent menu content
            menu.draw_viewer_menu();

            if (ImGui::CollapsingHeader("Add displacement", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::Checkbox("Display displacement", &add_displacement);
                ImGui::InputDouble("Smooth displacement", &smooth_displacement);
            }
            if (ImGui::CollapsingHeader("Animation", ImGuiTreeNodeFlags_DefaultOpen)) {
                if (ImGui::Button("Restart", ImVec2(-1, 0)))
                {
                    cout << "restart" << endl;

                    // Precompute the unposed example
                    precomputeUnposedExamples();

                    // Restart chrono
                    start = std::chrono::system_clock::now();
                }
            }
        };
    } else if(currentMesh == BODY_FACE) {
        // Load unposed example
        precomputeUnposedExamplesPerFace();

        // Bind callback functions
        animationCallbackFunction = [&](igl::opengl::glfw::Viewer &) {
            // Clear
            viewer.data().clear();

            // Check counter
            auto current = chrono::system_clock::now();
            chrono::seconds s = chrono::duration_cast<chrono::seconds>(current - start);
            if (s.count() == animationFrames) {
                start = std::chrono::system_clock::now();
            }

            // Perform animation
            displayBodyPosePerFace(viewer, start);
            return false;
        };

        menuCallbackFunction = [&]() {
            // Draw parent menu content
            menu.draw_viewer_menu();

            if (ImGui::CollapsingHeader("Add displacement", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::Checkbox("Display displacement", &per_face_displacement);
                ImGui::Checkbox("Divide rotation and skew while computing displacements", &rotation_and_skew);
            }
            if (ImGui::CollapsingHeader("Animation", ImGuiTreeNodeFlags_DefaultOpen)) {
                if (ImGui::Button("Restart", ImVec2(-1, 0)))
                {
                    cout << "restart" << endl;

                    // Restart chrono
                    start = std::chrono::system_clock::now();
                }
            }
        };
    } else {
        cerr << "Unknown mesh" << endl;
        exit(0);
    }

    // Bind functions and launch
    menu.callback_draw_viewer_menu = menuCallbackFunction;
    viewer.callback_pre_draw = animationCallbackFunction;
    viewer.launch();

    return 0;
}