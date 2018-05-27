#include <iostream>
#include "gui_utils.h"
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>

#include <SimulationMesh.h>
#include <FEMElementAnalytic.h>
#include <FEMElementFD.h>
#include <FEMElementAutoDiff.h>

// Viewer mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd colors_per_vertex;
Eigen::MatrixXd V_pinned;
Eigen::MatrixXd V_pinned_color;
Eigen::MatrixXd V_selected;
Eigen::MatrixXd V_selected_color;

Eigen::Vector2d V_com;

bool IsMouseDown = false;
int MouseDownMod = -1;

// sim parameters
double desiredFrameRate = 30;
bool isSimPlaying = false;
double maxStress = 0.02;
double bulkModulus = 50;
double shearModulus = 50;

// mouse interaction
enum MouseMode { SET_SUPPORT, DEFORM_NODES};
MouseMode mouse_mode = SET_SUPPORT;
int down_mouse_x = -1, down_mouse_y = -1;
int selected_v = -1;
int over_v = -1;

enum ElementType { SPRING, FEM_FD, FEM_ANALYTIC, FEM_AUTODIFF };
ElementType element_type = FEM_ANALYTIC;

enum CoM_ConstraintType { CENTER_OF_SUPPORT, SOMEWHERE_IN_SUPPORT };
CoM_ConstraintType com_constraint_type = CENTER_OF_SUPPORT;

// adapted from: https://stackoverflow.com/a/7811134
typedef struct { double r,g,b;} Color;
Color jetColorFromScalar(double v, double vmin, double vmax)
{
   Color c = {1.0,1.0,1.0}; // white
   double dv;

   if (v < vmin)
	  v = vmin;
   if (v > vmax)
	  v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
	  c.r = 0;
	  c.g = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
	  c.r = 0;
	  c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
	  c.r = 4 * (v - vmin - 0.5 * dv) / dv;
	  c.b = 0;
   } else {
	  c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
	  c.b = 0;
   }

   return(c);
}

SimulationMesh simMesh;

void loadTriangleMesh(const SimulationMesh::TriangleMesh &triMesh) {
	simMesh.loadTriangleMesh(triMesh);

	if(element_type == SPRING) {
		simMesh.createSpringsFromTriangles();
	}
	else if(element_type == FEM_FD) {
		simMesh.createElementsFromTriangles<FEMElementFD>();
	}
	else if(element_type == FEM_ANALYTIC) {
		simMesh.createElementsFromTriangles<FEMElementAnalytic>();
	}
	else if(element_type == FEM_AUTODIFF ) {
		simMesh.createElementsFromTriangles<FEMElementAutoDiff>();
	}

//	simMesh.addGravityForces(Vector2d(0, -9.8));
	//set some boundary conditions...
//	for (int i = 0; i < nCols;i++)
//	{
//		Vector2d p = simMesh.x.segment<2>(2*i);
//		simMesh.setPinnedNode(i, p);
//	}
//	for (int i = 0; i < nRows;i++)
//	{
//		int ii = nCols*i;
//		Vector2d p = simMesh.x.segment<2>(2*ii);
//		simMesh.setPinnedNode(ii, p);
//	}
}

void makeSimMesh() {

	int nRows = 10; //10 DEBUG CHANGE HERE
	int nCols = 10; //10 DEBUG CHANGE HERE
	SimulationMesh::TriangleMesh triMesh;
	SimulationMesh::generateSquareTriMesh(triMesh, -1, 0, 0.1, 0.1, nRows, nCols);

	loadTriangleMesh(triMesh);
}

void updateViewerVertices() {
	V.resize(simMesh.numNodes, 3);
	for (int i = 0; i < simMesh.numNodes; ++i) {
		for (int j = 0; j < 2; ++j) {
			V(i,j) = simMesh.x[2*i+j];
		}
//		V(i,2) = 0;
	}
}

void updateViewerVerticesColors() {
	// color vertices according to deformation
	VectorXd defEnergyPerNode;
	simMesh.computeDefoEnergyPerNode(defEnergyPerNode);
	maxStress = std::max(1e-5, defEnergyPerNode.maxCoeff());
	colors_per_vertex.resize(defEnergyPerNode.size(), 3);
	for (int i = 0; i < defEnergyPerNode.size(); ++i) {
		double v = defEnergyPerNode[i];
		Color color = jetColorFromScalar(v, 0, maxStress);
		colors_per_vertex.row(i) << color.r, color.g, color.b;
	}
}

void updateGlobalPinnedVertices() {

	int nVertsPinned = simMesh.pinnedNodeElements.size();
	V_pinned.resize(nVertsPinned, 3);
	V_pinned_color.resize(nVertsPinned, 3);

	int i = 0;
	for(Element* pinnedNode : simMesh.pinnedNodeElements) {
		int nodeIndex = pinnedNode->getNodeIndex(0);
		for (int j = 0; j < 2; ++j)
			V_pinned(i,j) = simMesh.x[2*nodeIndex+j];
		V_pinned(i, 2) = 0;
		V_pinned_color.row(i) << 0, 0, 0;
		i++;
	}

}

void updateViewerFaces() {
	int nTriangles = simMesh.triangles.size();
	F.resize(nTriangles, 3);
	for (int i = 0; i < nTriangles; ++i) {
		for (int j = 0; j < 3; ++j) {
			F(i,j) = simMesh.triangles[i][2-j];
		}
	}
}

void update_display(igl::viewer::Viewer& viewer) {
	viewer.data.clear();
	viewer.data.set_mesh(V, F);
	viewer.data.compute_normals();
	viewer.data.set_colors(colors_per_vertex);

	if(over_v!=-1){
		Eigen::MatrixXd points(1, 3);
		points.row(0) << simMesh.x[2*over_v+0], simMesh.x[2*over_v+1], 0;
		Eigen::MatrixXd colors(1, 3);
		colors.row(0) << 1, 1, 0;
		viewer.data.add_points(points, colors);
	}

	Eigen::MatrixXd nodes(simMesh.constraintNodes.size(),3);
	Eigen::MatrixXd nodesColor(simMesh.constraintNodes.size(),3);
	int i=0;
	for (const auto &f : simMesh.constraintNodes) {
		nodes.row(i) << f.second[0], f.second[1], 0;
		if(simMesh.supportNodes.find(f.first) != simMesh.supportNodes.end())
			nodesColor.row(i) << 0, 0, 0;
		else
			nodesColor.row(i) << 1, 0, 1;
		i++;
	}
	viewer.data.add_points(nodes, nodesColor);
}

void solve() {
	simMesh.solveShapeOpt(com_constraint_type);
}

bool callback_pre_draw(igl::viewer::Viewer& viewer) {
	if(isSimPlaying){
		solve();
	}
	updateViewerVertices();
	updateViewerVerticesColors();
	updateGlobalPinnedVertices();
	update_display(viewer);
	return false;
}

bool callback_mouse_down(igl::viewer::Viewer& viewer, int button, int modifier) {
	down_mouse_x = viewer.current_mouse_x;
	down_mouse_y = viewer.current_mouse_y;
	IsMouseDown = true;
	MouseDownMod = modifier;

	if (mouse_mode == SET_SUPPORT) {
		int v = pick_vertex(viewer, down_mouse_x, down_mouse_y, V, F);
		if (v !=-1) {
			if(MouseDownMod == GLFW_MOD_CONTROL){
				simMesh.constraintNodes.erase(v);
				simMesh.supportNodes.erase(v);
			}
			else {
				Vector2d pos(simMesh.x[2*v+0], simMesh.x[2*v+1]);
				simMesh.constraintNodes[v] = pos;
				simMesh.supportNodes.insert(v);
				selected_v = v;
			}
			update_display(viewer);
			return true;
		}
	}
	else if (mouse_mode == DEFORM_NODES) {
		int v = pick_vertex(viewer, down_mouse_x, down_mouse_y, V, F);
		if (v !=-1) {
			if(MouseDownMod == GLFW_MOD_CONTROL){
				simMesh.constraintNodes.erase(v);
				simMesh.supportNodes.erase(v);
			}
			else {
				Vector2d pos(simMesh.x[2*v+0], simMesh.x[2*v+1]);
				simMesh.constraintNodes[v] = pos;
				simMesh.supportNodes.erase(v);
				selected_v = v;
			}
			update_display(viewer);
			return true;
		}
	}

	return false;
}

bool callback_mouse_move(igl::viewer::Viewer& viewer, int mouse_x, int mouse_y) {

	// mark point hovering over
	over_v = pick_vertex(viewer, mouse_x, mouse_y, V, F);

	// move constrained node
	if(IsMouseDown && MouseDownMod == 0 && selected_v != -1) {
		Eigen::Vector3f pos = computePosition(viewer, mouse_x, mouse_y, V.row(selected_v));
		simMesh.constraintNodes[selected_v] = Vector2d(pos[0], pos[1]);
		return true;
	}

	// bulk add/remove
	if (IsMouseDown) {
		int v = pick_vertex(viewer, mouse_x, mouse_y, V, F);
		if (v !=-1) {
			if(MouseDownMod == GLFW_MOD_CONTROL) {
				simMesh.constraintNodes.erase(v);
				simMesh.supportNodes.erase(v);
			}
			else if(MouseDownMod == GLFW_MOD_SHIFT) {
				simMesh.constraintNodes[v] = Vector2d(simMesh.x[2*v+0], simMesh.x[2*v+1]);
				if(mouse_mode == SET_SUPPORT)
					simMesh.supportNodes.insert(v);
				else
					simMesh.supportNodes.erase(v);
			}
			update_display(viewer);
			return true;
		}
	}

	return false;
}

bool callback_mouse_up(igl::viewer::Viewer& viewer, int button, int modifier) {

	IsMouseDown = false;
	selected_v = -1;

	return false;
}

bool callback_init(igl::viewer::Viewer& viewer) {

	viewer.ngui->addGroup("Simulation");

	nanogui::ComboBox *widgetElement = viewer.ngui->addVariable("Element Type", element_type, true);
	widgetElement->setItems({ "spring", "FEM finite differences", "FEM analytic", "FEM AutoDiff" });
	widgetElement->setCallback([](int i){
		if(i == 0) {
			element_type = SPRING;
			simMesh.createSpringsFromTriangles();
			simMesh.x = simMesh.X;
		}
		else if(i == 1) {
			element_type = FEM_FD;
			simMesh.createElementsFromTriangles<FEMElementFD>();
			simMesh.x = simMesh.X;
		}
		else if(i == 2){
			element_type = FEM_ANALYTIC;
			simMesh.createElementsFromTriangles<FEMElementAnalytic>();
			simMesh.x = simMesh.X;
		}
		else {
			element_type = FEM_AUTODIFF;
			simMesh.createElementsFromTriangles<FEMElementAutoDiff>();
			simMesh.x = simMesh.X;
		}
	});

	viewer.ngui->addVariable("CoM Constraint", com_constraint_type, true)->setItems({ "center of support", "inside the support" });

	// play button
	nanogui::Button *buttonPlay = viewer.ngui->addButton("Play",[]{});
	buttonPlay->setFlags(nanogui::Button::Flags::ToggleButton);
	buttonPlay->setChangeCallback([](bool pushed){ isSimPlaying = pushed;});
	// step button
	nanogui::Button *buttonStep = viewer.ngui->addButton("Step",[]{});
	buttonStep->setCallback([](){solve();});

	viewer.ngui->addVariable("Color Scale", maxStress);

	viewer.ngui->addVariable("Mouse Mode", mouse_mode, true)->setItems({ "Set Support", "Deform handles"});

	viewer.screen->setSize({1200, 1000});
	viewer.screen->performLayout();
	return true;
}

bool load_mesh(igl::viewer::Viewer &viewer, std::string filename)
{
	igl::read_triangle_mesh(filename,V,F);

	SimulationMesh::TriangleMesh triMesh;
	triMesh.v.resize(V.rows());
	for (int i = 0; i < V.rows(); ++i) {
		triMesh.v[i] = {V(i,0), V(i,1)};
	}
	triMesh.tri.resize(F.rows());
	for (int i = 0; i < F.rows(); ++i) {
		triMesh.tri[i] = {F(i,0), F(i,1), F(i,2)};
	}
	loadTriangleMesh(triMesh);

	selected_v = -1;
	colors_per_vertex.resize(V.size(), 3);
	colors_per_vertex.setOnes();
	update_display(viewer);
	viewer.core.align_camera_position(V);

	return true;
}

bool callback_load_mesh(igl::viewer::Viewer &viewer, std::string filename)
{
  load_mesh(viewer, filename);

  return true;
}

int main(int argc, char *argv[]) {

	makeSimMesh();

	updateViewerVertices();
	updateViewerFaces();
	updateViewerVerticesColors();

	using Viewer = igl::viewer::Viewer;

	Viewer viewer;
	viewer.data.point_size = 10.0;
	viewer.data.line_width = 3.0;
	viewer.core.align_camera_center(V,F);
	viewer.core.is_animating = true;
	viewer.callback_init = callback_init;
	viewer.callback_load_mesh = callback_load_mesh;
	viewer.callback_pre_draw = callback_pre_draw;

	viewer.callback_mouse_down = callback_mouse_down;
	viewer.callback_mouse_up = callback_mouse_up;
	viewer.callback_mouse_move = callback_mouse_move;

	viewer.launch();
}
