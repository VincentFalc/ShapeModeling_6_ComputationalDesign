#pragma once

#include <math.h>
#include <MathHelper.h>

/**
	This class implements the interface for an elementary energy unit. As a function of deformed, undeformed, 
	and other parameters, such as boundary conditions, each class that extends this one will define a potential energy.
	The deformed energy depends on a number of nodes.
*/
class Element{

public:
	Element() {}
	virtual ~Element() {}

	// Returns the number of nodes this unit depends on
	virtual int getNumNodes() const = 0;
	// Returns the global index of node `i`
	virtual int getNodeIndex(int i) const = 0;

	Vector2d getNodePos(int i, const VectorXd &x) const {
		return x.segment<2>(2*getNodeIndex(i));
	}

	// Returns the element's mass
	virtual double getMass() const = 0;

	// Returns the energy value given deformed `x` and undeformed `X` state
	virtual double getEnergy(const VectorXd& x, const VectorXd& X) = 0;
	// Adds the gradient to `grad` given deformed `x` and undeformed `X` state
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) = 0;
	// Adds the hessian entries to `hesEntries` given deformed `x` and undeformed `X` state
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) = 0;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
