#pragma once

#include "Element.h"

/**
	This class implements the interface for an elementary energy unit. As a function of deformed, undeformed,
	and other parameters, such as boundary conditions, each class that extends this one will define a potential energy.
	The deformed energy depends on a number of nodes.
*/
class Spring : public Element {

public:
	Spring(const std::array<int, 2> &nodeIndices, const VectorXd &X)
		: nodeIndices(nodeIndices) {
	}
	virtual ~Spring() {}

	// Returns the number of nodes this unit depends on
	virtual int getNumNodes() const {
		return 2;
	}
	// Returns the global index of node `i`
	virtual int getNodeIndex(int i) const {
		return nodeIndices[i];
	}

	// Returns the element's mass
	virtual double getMass() const {
		return 0;
	}

	// Returns the energy value given deformed `x` and undeformed `X` state
	virtual double getEnergy(const VectorXd& x, const VectorXd& X) {

		// EXERCISE 1.2

		Vector2d V = getNodePos(0, X) - getNodePos(1, X);
		Vector2d v = getNodePos(0, x) - getNodePos(1, x);
		double l = v.norm(); double L = V.norm();

		return 0.5 * k * std::pow(l/L-1,2.0)*L;
	}

	// Adds the gradient to `grad` given deformed `x` and undeformed `X` state
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {

		// EXERCISE 1.2

		Vector2d dEdx[2];

		Vector2d V = getNodePos(0, X) - getNodePos(1, X);
		Vector2d v = getNodePos(0, x) - getNodePos(1, x);
		double l = v.norm(); double L = V.norm();

		Vector2d f = -k*(l/L-1)/l * v;

		dEdx[0] = -f;
		dEdx[1] = f;

		// Add the local gradient `dEdx` to the global gradient `grad`
		for (int i = 0;i<2;i++)
			for (int j = 0;j<2;j++)
				grad[2*nodeIndices[i] + j] += dEdx[i][j];
	}

	// Adds the hessian entries to `hesEntries` given deformed `x` and undeformed `X` state
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {

		// EXERCISE 1.4

		Matrix2d ddEdxdx[2][2];

		Vector2d V = getNodePos(0, X) - getNodePos(1, X);
		Vector2d v = getNodePos(0, x) - getNodePos(1, x);
		double l = v.norm(); double L = V.norm();

		double epsilon = l/L-1;
		Matrix2d vv = v*v.transpose()/(v.dot(v));
		Matrix2d dfdx1 = -k*(1/L*vv + epsilon/l*(Matrix2d::Identity()-vv));

		ddEdxdx[0][0] = -dfdx1;
		ddEdxdx[0][1] = dfdx1;

		ddEdxdx[1][0] = dfdx1;
		ddEdxdx[1][1] = -dfdx1;

		// Add the local hessians `ddEdxdx` to the global hessian using Eigen::Triplets
		for (int i = 0; i<2; i++)
			for (int j = 0; j< 2;j++)
				addSparseMatrixDenseBlockToTriplet(hesEntries, 2*nodeIndices[i], 2*nodeIndices[j], ddEdxdx[i][j], true);
	}

protected:
	// the collection of nodes that define the triangle element
	std::array<int, 2> nodeIndices;
	// spring stiffness
	double k = 20.0;
};
