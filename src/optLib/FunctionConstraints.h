#pragma once

#include <math.h>
#include "MathHelper.h"

/*!
	A multi-dimensional function that expresses linear equality and inequality constraints applicable to an objective function:

	Equality constrains:    A(p) = b
	Inequality constraint:  d <= C(p) <= f
	Bound constraint:       l <= p <= u
*/
class FunctionConstraints {
public:
	/*! Constructor
	* Fill the constant variables in the constructor of the derived class:
	*  b, d, f, l, u
	*/
	FunctionConstraints() {}

	virtual ~FunctionConstraints() {}

	// Returns the number of equality constraints.
	virtual int getEqualityConstraintCount() {
		return (int)getEqualityConstraintsTargetValues().size();
	}

	// Returns b of A(p) = b.
	virtual const VectorXd& getEqualityConstraintsTargetValues() {
		return b;
	}

	// Returns A(p) of A(p) = b.
	// Derive from this to compute A(p). Fill `eqConstraintsVals` and return it.
	virtual const VectorXd& getEqualityConstraintValues(const VectorXd& p) {
		return eqConstraintVals;
	}

	// Computes the Jacobian dA/dp of the equality constraints A.
	virtual void addEqualityConstraintsJacobianEntriesTo(std::vector<Tripletd>& jacobianEntries, const VectorXd& p) {

	}

	// Returns the number of inequality constraints.
	virtual int getInequalityConstraintCount() {
		return (int)getInequalityConstraintsMinValues().size();
	}

	// Returns the value of the inequality constraint C(p).
	virtual const VectorXd& getInequalityConstraintValues(const VectorXd& p) {
		return ineqConstraintVals;
	}

	// Returns d of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMinValues() {
		return d;
	}

	// Returns f of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMaxValues() {
		return f;
	}

	// Computes the Jacobian dA/dp of the inequality constraints C.
	virtual void addInequalityConstraintsJacobianEntriesTo(std::vector<Tripletd>& jacobianEntries, const VectorXd& p) {

	}

	// Returns l of constraint l <= p <= u
	virtual const VectorXd& getBoundConstraintsMinValues() {
		return l;
	}

	// Returns u of constraint l <= p <= u
	virtual const VectorXd& getBoundConstraintsMaxValues() {
		return u;
	}

protected:
	VectorXd b;
	VectorXd d;
	VectorXd f;
	VectorXd l;
	VectorXd u;
	VectorXd eqConstraintVals;		// A(p)
	VectorXd ineqConstraintVals;	// C(p)
};
