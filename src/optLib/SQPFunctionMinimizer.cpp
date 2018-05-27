#include "SQPFunctionMinimizer.h"

//TODO: all jacobians and constraints should be represented in terms of the triplets...
//#include "MathHelper.h"

#include "ObjectiveFunction.h"

#include "OoqpEigenInterface.h"
#include "ooqpei_assert_macros.h"

#include <fstream>
#include <iostream>

SQPFunctionMinimizer::SQPFunctionMinimizer(int maxIterations, double solveResidual, int maxLineSearchIterations)
	:	maxIterations(maxIterations),
		solveResidual(solveResidual),
		maxLineSearchIterations(maxLineSearchIterations)
{
}

SQPFunctionMinimizer::~SQPFunctionMinimizer(){
}

void SQPFunctionMinimizer::computeGradient(ObjectiveFunction *objective, const VectorXd& pi) {
	int nParameters = (int)pi.size();
	resize(gradient, nParameters);
	gradient.setZero();
	objective->addGradientTo(gradient, pi);
}

void SQPFunctionMinimizer::computeHessian(ObjectiveFunction *objective, const VectorXd& pi) {
	int nParameters = (int)pi.size();
	resize(H, nParameters, nParameters);
	smEntries.clear();
	objective->addHessianEntriesTo(smEntries, pi);
	H.setFromTriplets(smEntries.begin(), smEntries.end());
}

void SQPFunctionMinimizer::computeConstraintsAndJacobians(FunctionConstraints* constraints, const VectorXd& pi) {
	int nParameters = (int)pi.size();

	smEntries.clear();
	resize(A, constraints->getEqualityConstraintCount(), nParameters);
	constraints->addEqualityConstraintsJacobianEntriesTo(smEntries, pi);
	A.setFromTriplets(smEntries.begin(), smEntries.end());

	smEntries.clear();
	resize(C, constraints->getInequalityConstraintCount(), nParameters);
	constraints->addInequalityConstraintsJacobianEntriesTo(smEntries, pi);
	C.setFromTriplets(smEntries.begin(), smEntries.end());

	// compute constraints
	computeConstraints(constraints, pi, ck);
}


/**
	min f(p) subject to the constraints...
*/
bool SQPFunctionMinimizer::minimize(ObjectiveFunction *objective, FunctionConstraints* constraints, VectorXd &p){

	const int nParameters = (int)p.size();
	VectorXd dp = VectorXd::Zero(nParameters);
	VectorXd pi = p;

	// Iterate - like Newton
	bool optimizationConverged = false;
	int i;
	for (i=0; i<maxIterations; i++) {

		computeGradient(objective, pi);

		computeHessian(objective, pi);

		computeConstraintsAndJacobians(constraints, pi);

		// Find the direction to step at
		computeSearchDirection(H,
		                       gradient,
		                       pi,
		                       dp,
		                       A,
							   constraints->getEqualityConstraintsTargetValues(),
							   constraints->getInequalityConstraintsMinValues(),
		                       C,
							   constraints->getInequalityConstraintsMaxValues(),
							   constraints->getBoundConstraintsMinValues(),
							   constraints->getBoundConstraintsMaxValues());
		

		if (dp.norm() < solveResidual)	{
			lastNumberIterations = i;
			optimizationConverged = true;
			break;
		}

		// Do a line search
		double alpha = doLineSearch(objective, constraints, pi, dp, maxLineSearchIterations);

		pi += alpha*dp;
	}

	lastNumberIterations = i;
	p = pi;
	return optimizationConverged;

}

void SQPFunctionMinimizer::computeSearchDirection(const SparseMatrixd& hessian,
                                                  const VectorXd& gradient,
                                                  const VectorXd &p,
                                                  VectorXd &dp,
												  const SparseMatrixd& A,
                                                  const VectorXd& b,
                                                  const VectorXd& d,
												  const SparseMatrixd& C,
                                                  const VectorXd& f,
                                                  const VectorXd& l,
                                                  const VectorXd& u) {

	/** We want dp to minimize: F(p+dp) ~ F(p) + dp' grad + 1/2 dp' H dp
	* while maintaining the constraints A*(p+dp) = b and d <= C*(p+dp) <= f and l <= p+dp <= u.
	* re-writing the constraints as
	* A*dp = b - A*p
	* d-C*p <= C*dp <= f-C*p
	* l-p <= dp <= u-p
	* we can get a canonical QP form:
	* min 1/2 x' Q x + c' x s. t. A x = b, d <= Cx <= f, and l <= x <= u
	* where:
	* x = dp
	* Q = hessian
	* c = gradient
	* A = A
	* b = bMinusAp
	* C = C
	* d = dMinusCp
	* f = fMinusCp
	* l = minMinusp
	* u = maxMinusp
	*/

	int np = p.size();

	VectorXd fMinusCp;
	VectorXd dMinusCp;

	if (C.size() > 0) {
		const VectorXd Cp = C*p;
		fMinusCp = f - Cp;
		dMinusCp = d - Cp;
	}

	VectorXd bMinusAp;
	if (A.size() > 0) {
		bMinusAp = b - A*p;
	}

	VectorXd lMinusdp = l - p;
	VectorXd uMinusdp = u - p;

	bool success = ooqpei::OoqpEigenInterface::solve(hessian,
										   gradient,
										   A,
										   bMinusAp,
										   C,
										   dMinusCp,
										   fMinusCp,
										   lMinusdp,
										   uMinusdp,
										   dp);
}

double SQPFunctionMinimizer::doLineSearch(ObjectiveFunction *objective, FunctionConstraints *constraints, VectorXd &p, const VectorXd &dp, int maxSteps){

	double alpha = 1.0;
	const double mu = computeMu(p);
	const double initialValue = computeMerritFunction(objective, constraints, p, mu);
	VectorXd pc;

	for (int j=0; j<maxSteps; j++)
	{
		pc = p + alpha*dp;

//		const double mu = computeMu(p);
		const double newLineSearchValue = computeMerritFunction(objective, constraints, pc, mu);

		if (!std::isfinite(newLineSearchValue) || newLineSearchValue >= initialValue) {
			alpha /= 2.0;
		}
		else {
			return alpha;
		}
	}

	return alpha;
}

double SQPFunctionMinimizer::computeMerritFunction(ObjectiveFunction *objective, FunctionConstraints *constraints, const VectorXd &x, double mu) {
	double f = objective->computeValue(x);

	VectorXd c_at_x;
	computeConstraints(constraints, x, c_at_x);
	double c_at_x_norm1 = 0;
	for (int i = 0; i < c_at_x.size(); ++i)
		c_at_x_norm1 += abs(c_at_x[i]);

	return f + mu*c_at_x_norm1;
}

double SQPFunctionMinimizer::computeMu(const VectorXd &x) const
{
	double ck_norm1 = 0;
	for (int i = 0; i < ck.size(); ++i)
		ck_norm1 += abs(ck[i]);

	double xTHx = x.transpose()*H*x;
	double sigma = (xTHx > 0) ? 1.0 : 0.0;
	double rho = 0.1;
	double mu_new = (abs(gradient.dot(x) + 0.5*sigma*xTHx)/((1.0-rho)*ck_norm1));

	mu_new = std::max(mu_max, std::min(mu_min, mu_new));

	mu_k = std::max(mu_k, mu_new);

//	mu_k = mu_new;
	return mu_k;
}

void SQPFunctionMinimizer::computeConstraints(FunctionConstraints *constraints, const VectorXd &x, VectorXd &c)
{
	VectorXd c_ineq = constraints->getInequalityConstraintValues(x);
	VectorXd c_min = constraints->getInequalityConstraintsMinValues();
	VectorXd c_max = constraints->getInequalityConstraintsMaxValues();

	VectorXd c_eq = constraints->getEqualityConstraintValues(x);
	VectorXd c_val = constraints->getEqualityConstraintsTargetValues();

	c.resize(c_ineq.size() + c_eq.size());
	c.setZero();
	for (int i = 0; i < c_ineq.size(); ++i) {
		if(c_ineq[i] < c_min[i])
			c[i] = c_min[i]-c_ineq[i];
		if(c_ineq[i] > c_max[i])
			c[i] = c_ineq[i]-c_max[i];
	}
	for (int i = 0; i < c_eq.size(); ++i)
		c[c_ineq.size()+i] = c_eq[i] - c_val[i];
}
