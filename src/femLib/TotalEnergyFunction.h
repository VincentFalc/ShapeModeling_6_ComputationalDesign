#pragma once

#include <ObjectiveFunction.h>

#include "Element.h"

class SimulationMesh;

class TotalEnergyFunction : public ObjectiveFunction {
public:
	TotalEnergyFunction(SimulationMesh* simMesh);
	virtual ~TotalEnergyFunction(void);

	virtual double computeValue(const VectorXd& s);
	virtual void addGradientTo(VectorXd& grad, const VectorXd& s);

	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& s);

	std::vector<Element*> elements;

public:
	SimulationMesh* simMesh;
};
