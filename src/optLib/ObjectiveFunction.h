#pragma once

#include "MathHelper.h"

class ObjectiveFunction{
public:
	// this should always return the current value of the objective function
	virtual double computeValue(const VectorXd& x) = 0;
	virtual void addGradientTo(VectorXd& grad, const VectorXd& x) = 0;
	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) = 0;
};
