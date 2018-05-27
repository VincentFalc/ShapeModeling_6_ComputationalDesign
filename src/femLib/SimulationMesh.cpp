#include "SimulationMesh.h"

#include <NewtonFunctionMinimizer.h>
#include <GradientDescentFunctionMinimizer.h>
#include <SQPFunctionMinimizer.h>

#include <TotalEnergyFunction.h>
#include <EqualityOverload.h>
#include <EqualityOverloadwithInequality.h>

#include "FEMElementAnalytic.h"
#include "Spring.h"

SimulationMesh::SimulationMesh(){
	energyFunction = new TotalEnergyFunction(this);
}

SimulationMesh::~SimulationMesh(){
	delete energyFunction;
}

void SimulationMesh::generateSquareTriMesh(SimulationMesh::TriangleMesh &triMesh, double startX, double startY, double dX, double dY, int xSize, int ySize) {
	triMesh.v.resize(xSize*ySize);
	for (int i=0; i<xSize; i++)
		for (int j=0; j<ySize;j++)
			triMesh.v[i*ySize+j] = {startX + (double)i * dX, startY + (double)j * dY};

	triMesh.tri.resize(2*(xSize-1)*(ySize-1));
	for (int i=0; i<xSize-1; i++)
		for (int j=0; j<ySize-1; j++){
			triMesh.tri[2*(i*(ySize-1) + j) + 0] = {{i * ySize + j, i * ySize + (j + 1), (i+1)*ySize + j}};
			triMesh.tri[2*(i*(ySize-1) + j) + 1] = {{i * ySize + (j+1), (i+1) * ySize + (j + 1), (i+1)*ySize + j}};
		}
}

void SimulationMesh::loadTriangleMesh(const TriangleMesh &triMesh)
{
	clear();

	numNodes = triMesh.v.size();
	numElements = triMesh.tri.size();

	x.resize(2 * numNodes);
	X.resize(2 * numNodes);
	f_ext.resize(2 * numNodes);
	m.resize(2 * numNodes);

	for (int i=0; i<numNodes; i++){
		Vector2d p = triMesh.v[i];
		x[2 * i + 0] = p[0]; x[2 * i + 1] = p[1];
		X[2 * i + 0] = p[0]; X[2 * i + 1] = p[1];
		f_ext[2 * i + 0] = 0; f_ext[2 * i + 1] = 0;
		//the masses for the node are obtained by lumping together/distributing the mass of the elements that share the nodes...
		m[2 * i + 0] = 0; m[2 * i + 1] = 0;
	}

	triangles = triMesh.tri;
}

void SimulationMesh::createSpringsFromTriangles()
{
	energyFunction->elements.clear();
	for (int i=0; i<numElements; i++){

		for (int j = 0; j < 3; ++j) {
			std::array<int, 2> nodes;
			nodes[0] = triangles[i][j];
			int j2 = (j<2) ? j+1 : 0;
			nodes[1] = triangles[i][j2];
			Spring* spring = new Spring(nodes, X);
			energyFunction->elements.push_back(spring);
		}
	}

	Vector2d pMin = {HUGE_VAL, HUGE_VAL};
	Vector2d pMax = {-HUGE_VAL, -HUGE_VAL};
	for (int i = 0; i < X.size()/2; ++i) {
		for (int j = 0; j < 2; ++j) {
			pMin[j] = std::min(pMin[j], X[2*i+j]);
			pMax[j] = std::max(pMax[j], X[2*i+j]);
		}
	}

	m.setZero();
	double massDensity = 1.0;
	double totalMass = massDensity*((pMax[0]-pMin[0])*(pMax[1]-pMin[1]));
	double massPerNode = totalMass / (X.size()/2);
	for (int j=0; j<m.size(); j++)
		this->m[j] = massPerNode;
}

void SimulationMesh::addGravityForces(const Vector2d &g){
	for (size_t i=0;i<numNodes;i++)
		for (size_t j=0;j<2;j++)
			f_ext[2*i + j] = g[j] * m[2*i + j];
}

void SimulationMesh::applyForceAt(int i, const Vector2d &force)
{
	f_ext[2*i+0] = force(0);
	f_ext[2*i+1] = force(1);
}

void SimulationMesh::clearExtForces()
{
	f_ext.setZero();
}

void SimulationMesh::solveGradientDescent()
{
	solve_statics<GradientDescentFunctionMinimizer>();
}

void SimulationMesh::solveNewtonsMethod()
{
	solve_statics<NewtonFunctionMinimizer>();
}

void SimulationMesh::solveShapeOpt(int COMconstraintType)
{
	/* Ex 2.1, 2.2: solve shape optimization problem using the SQPFunctionMinimizer
	 * make sure you use the correct constraint type as set in GUI by using COMconstraintType
	 */
	ObjectiveFunction* objective;
	//FunctionConstraints* constraints;
	//Ex2Constraint* constraints;
	Ex2ConstraintIneque* constraints;

	objective = energyFunction; // Or TotalEnergyFunction(this) ? //Directly an instance of ObjectiveFunction with Energy calculation ? Seems correct ? 
	constraints = new Ex2ConstraintIneque;//Ex2Constraint;

	//We give a ref of the current class to the constraints instance
	constraints->setRef(this);

	if(COMconstraintType==0){ // use the constraints from 2.1 = EQUALITY
		constraints->setEquality();
	} else if (COMconstraintType==1){  // use the constraints from 2.2 = INEQUALITY
		constraints->setINEquality();
	}

	SQPFunctionMinimizer minimizer;
	minimizer.minimize(objective, constraints, x);

	constraints->printMass();
/*
	std::cout << "x              = " << x.transpose() << std::endl;
	std::cout << "f(x)           = " << objective->computeValue(x) << std::endl;
	std::cout << "c(x)           = " << constraints->getEqualityConstraintValues(x) << std::endl;
	std::cout << "# iterations   = " << minimizer.lastNumberIterations << std::endl;
*/
}

void SimulationMesh::togglePinnedNode(int i)
{
	auto it = std::find_if(pinnedNodeElements.begin(), pinnedNodeElements.end(),
						[i](Element *e){
							return e->getNodeIndex(0) == i;
						});
	if(it == pinnedNodeElements.end())
		setPinnedNode(i, { x[2*i+0], x[2*i+1] });
	else
		unpinNode(i);
}

void SimulationMesh::setPinnedNode(int nodeIndex, const Vector2d &point)
{
	for (auto it = pinnedNodeElements.begin(); it != pinnedNodeElements.end(); ++it) {
		FixedPointElement* fps = dynamic_cast<FixedPointElement*>(*it);
		if (fps->getNodeIndex(0) == nodeIndex) {
			fps->setTargetPosition(point);
			return;
		}
	}
	pinnedNodeElements.push_back(new FixedPointElement(nodeIndex, point));
}

void SimulationMesh::unpinNode(int nodeIndex)
{
	for (auto it = pinnedNodeElements.begin(); it != pinnedNodeElements.end(); ++it) {
		FixedPointElement* fps = dynamic_cast<FixedPointElement*>(*it);
		if (fps->getNodeIndex(0) == nodeIndex) {
			pinnedNodeElements.erase(it);
			break;
		}
	}
}

void SimulationMesh::computeDefoEnergyPerNode(VectorXd &defEnergyPerNode)
{
	defEnergyPerNode.resize(numNodes);
	defEnergyPerNode.setZero();

	for (auto element : energyFunction->elements) {
		double energy = element->getEnergy(x, X);
		for (int i = 0; i < element->getNumNodes(); ++i) {
			defEnergyPerNode[element->getNodeIndex(i)] += energy / (double)element->getNumNodes();
		}
	}
}

void SimulationMesh::clear(){
	numNodes = 0;
	energyFunction->elements.clear();
	pinnedNodeElements.clear();
	x.resize(0);

}

