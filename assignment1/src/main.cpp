#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui.h>

/*** insert any libigl headers here ***/
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/edge_topology.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

void subdivide_sqrt3(const Eigen::MatrixXd &V,
					 const Eigen::MatrixXi &F,
					 Eigen::MatrixXd &Vout,
					 Eigen::MatrixXi &Fout){

    // Step 1 : Compute bc coordinates
    Eigen::MatrixXd BC;
    igl::barycenter(V,F,BC);

    // Step 2 : Setup new DS
    Vout.setZero(V.rows() + F.rows(), V.cols());
    Fout.setZero(3*F.rows(), F.cols());
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(V, F,EV,FE,EF);

    // Step 3 : Fill DS
    Vout << V,BC;

    unsigned int currentIndex = 0;
    set<unsigned int> cornersVetices;

    for(int i = 0; i < EV.rows();i++) {
        int vertex1 = EV(i,0);
        int vertex2 = EV(i,1);
        int face1 = EF(i,0) + V.rows();
        int face2 = EF(i,1) + V.rows();

        if(face1 == V.rows()-1) {
            Fout.row(currentIndex++) = Eigen::Vector3i(face2, vertex2, vertex1);
            cornersVetices.insert(vertex1);
        } else if(face2 == V.rows()-1) {
            Fout.row(currentIndex++) = Eigen::Vector3i(face1, vertex1, vertex2);
            cornersVetices.insert(vertex2);
        } else {
            Fout.row(currentIndex++) = Eigen::Vector3i(face2,face1,vertex1);
            Fout.row(currentIndex++) = Eigen::Vector3i(face1,face2,vertex2);
        }
    }

    // Step 4 : Modify position of previous vertices
    auto compute_an = [](const int n) -> double {
        return (4.0 - 2.0 * cos(2.0  * M_PI / n)) / 9.0;
    };

    auto transform_point = [](const Eigen::Vector3d& v, const Eigen::Vector3d& vi, const double a_i, const int n) -> Eigen::Vector3d  {
        return (1.0 - a_i) * v + a_i / n * vi;
    };

    igl::adjacency_list(F,VV);

    for(int i = 0; i < V.rows();i++) {
        if(cornersVetices.find(i) != cornersVetices.end())  {
            continue;
        }

        Eigen::Vector3d vi = Eigen::Vector3d::Zero();

        const unsigned int neighbours_size = VV[i].size();
        for(auto neighbour : VV[i]) {
            vi += V.row(neighbour);
        }

        double a_n = compute_an(neighbours_size);
        Vout.row(i) = transform_point(V.row(i),vi, a_n , neighbours_size);
    }
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        igl::vertex_triangle_adjacency(V.rows(),F,VF,VFi);

        for(int i = 0; i < VF.size();i++) {
            cout << "| Vertex : " << i << " | ";
            cout << "Connected faces : ";
            for(const int &face: VF[i]) {
                cout << face << " ";
            }
            cout << " |\n";
        }
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
        igl::adjacency_list(F, VV);

        for(int i = 0; i < VV.size();i++) {
            cout << "| Vertex : " << i << " | ";
            cout << "Connected vertices : ";
            for(const int &vertex: VV[i]) {
                cout << vertex << " ";
            }
            cout << " |\n";
        }
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
        igl::per_face_normals(V,F,FN);

        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V,F,VN);

        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
        const double THRESHOLD_CORNER_VALUE = 85;
        igl::per_corner_normals(V,F,THRESHOLD_CORNER_VALUE,CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        igl::facet_components(F,cid);

        // Print some information about the components
        vector<int> sizesOfComp(cid.maxCoeff()+1,0);
        for(int i = 0; i < cid.size();i++) {
            sizesOfComp[cid(i)]++;
        }

        for(int i = 0; i < sizesOfComp.size();i++) {
            cout << "Component " << i << " has " << sizesOfComp[i] << " faces." << "\n";
        }

        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        igl::jet(cid, true, component_colors_per_face);

        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
		Eigen::MatrixXd Vout;
		Eigen::MatrixXi Fout;

        // Fill the subdivide_sqrt3() function with your code for sqrt(3) subdivision.
        subdivide_sqrt3(V,F,Vout,Fout);

        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V,F);
  viewer.data().compute_normals();
  viewer.core().align_camera_center(V, F);
  return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else {
        filename = std::string("../data/plane.off"); // Default mesh
    }
	
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiPlugin plugin;
	viewer.plugins.push_back(&plugin);
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    plugin.widgets.push_back(&menu);
    
    viewer.launch();
}
