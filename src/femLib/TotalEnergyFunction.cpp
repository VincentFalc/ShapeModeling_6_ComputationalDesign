#include "TotalEnergyFunction.h"
#include "SimulationMesh.h"

#include "FixedPointElement.h"

TotalEnergyFunction::TotalEnergyFunction(SimulationMesh *simMesh)
	: simMesh(simMesh) {

}

TotalEnergyFunction::~TotalEnergyFunction(void){
}

//The net energy is: 1/2 a'M a + E + x'F, where E is the potential energy stored in the various elements
double TotalEnergyFunction::computeValue(const VectorXd& x){
	double totalEnergy = 0;

	for (size_t i=0;i<elements.size();i++)
		totalEnergy += elements[i]->getEnergy(x, simMesh->X);

	for (size_t i=0;i<simMesh->pinnedNodeElements.size();i++)
		totalEnergy += simMesh->pinnedNodeElements[i]->getEnergy(x, simMesh->X);

	totalEnergy -= x.dot(simMesh->f_ext);

	return totalEnergy;
}

void TotalEnergyFunction::addGradientTo(VectorXd& grad, const VectorXd& x) {
	if (grad.size() != x.size())
		resize(grad, x.size());

	//take into account the gradient of the deformation energy
	for (size_t i=0;i<elements.size();i++)
		elements[i]->addEnergyGradientTo(x, simMesh->X, grad);

	for (size_t i=0;i<simMesh->pinnedNodeElements.size();i++)
		simMesh->pinnedNodeElements[i]->addEnergyGradientTo(x, simMesh->X, grad);

	grad -= simMesh->f_ext;
}


void TotalEnergyFunction::addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {
	for (size_t i = 0;i < elements.size();i++)
		elements[i]->addEnergyHessianTo(x, simMesh->X, hessianEntries);

	for (size_t i = 0;i < simMesh->pinnedNodeElements.size();i++)
		simMesh->pinnedNodeElements[i]->addEnergyHessianTo(x, simMesh->X, hessianEntries);
}

//void TotalEnergyFunction::init()
//{
//	Vector2d left = {HUGE_VAL, HUGE_VAL};
//	Vector2d right = {-HUGE_VAL, -HUGE_VAL};
//	for(FixedPointElement* fe : simMesh->pinnedNodeElements) {
//		if(fe->getTargetPosition()[0] > right[0])
//			right = fe->getTargetPosition();
//		if(fe->getTargetPosition()[0] < left[0])
//			left = fe->getTargetPosition();
//	}
//	comObj.bounds[0] = left;
//	comObj.bounds[1] = right;
//}
