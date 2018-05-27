#pragma once

#include "ObjectiveFunction.h"
#include "FunctionConstraints.h"

/**
	Use the Sequential Quadratic Programming method to optimize a function, subject to constraints.
	
	Task: Find p that minimize f(p), such that Ap = b and d <= Cp <= f
*/
class SQPFunctionMinimizer {
 public:
  /*!
   *
   * @param maxIterations               maximum number of iterations
   * @param solveResidual               abortion criterium
   * @param maxLineSearchIterations     maximum number of line search iterations
   */
	SQPFunctionMinimizer(int maxIterations = 10, double solveResidual=0.0001, int maxLineSearchIterations = 20);
	virtual ~SQPFunctionMinimizer();

	/**
		min f(p) subject to the constraints...
	*/
	bool minimize(ObjectiveFunction *objective, FunctionConstraints* constraints, VectorXd &p);

protected:
	virtual void computeGradient(ObjectiveFunction *objective, const VectorXd& pi);
	virtual void computeHessian(ObjectiveFunction *objective, const VectorXd& pi);
	virtual void computeConstraintsAndJacobians(FunctionConstraints *constraints, const VectorXd& pi);

private:

	void computeSearchDirection(const SparseMatrixd& hessian,
							  const VectorXd& gradient,
							  const VectorXd &p,
							  VectorXd &dp,
							  const SparseMatrixd& A,
							  const VectorXd& b,
							  const VectorXd& d,
							  const SparseMatrixd& C,
							  const VectorXd& f,
							  const VectorXd& minBounds,
							  const VectorXd& maxBounds);

	double doLineSearch(ObjectiveFunction *function, FunctionConstraints *constraints, VectorXd &p, const VectorXd &dp, int maxSteps);

	static double computeMerritFunction(ObjectiveFunction *objective, FunctionConstraints *constraints, const Eigen::VectorXd &x, double mu_k);

	double computeMu(const Eigen::VectorXd &x) const;

	static void computeConstraints(FunctionConstraints *constraints, const Eigen::VectorXd &x, VectorXd &ck);

public:
	int maxLineSearchIterations;
	int maxIterations;
	double solveResidual;
	int lastNumberIterations = -1;

	double mu_min = 1e-10;
	double mu_max = 1e10;
	mutable double mu_k = 1e-10;

	SparseMatrixd H, A, C;
	std::vector<Tripletd> smEntries;
	VectorXd gradient, ck;
};

