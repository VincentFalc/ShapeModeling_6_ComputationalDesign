#include <iostream>

#include <FunctionConstraints.h>
#include <SQPFunctionMinimizer.h>
//#define DEBUG

class Ex1Objective : public ObjectiveFunction {
public:
	// Ex 1
	virtual double computeValue(const VectorXd& x){
#ifdef DEBUG
		std::cout << "> SQP - computeValue - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
#endif
		const double &x1 = x[0];
		const double &x2 = x[1];
		return std::pow(x1,2.0) + std::pow(x2, 2.0);
	}

	virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {
#ifdef DEBUG
		std::cout << "> SQP - addGradientTo - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> SQP - addGradientTo - size of grad (output) : " << grad.rows() << " by " << grad.cols() << std::endl;
#endif
		// Ex 1.1
		const double &x1 = x[0];
		const double &x2 = x[1];
		//grad.resize(2);
		grad[0] += 2*x1;
		grad[1] += 2*x2;
	}


	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {
#ifdef DEBUG
		std::cout << "> SQP - addHessianEntriesTo - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> SQP - addGradientTo - size of hessian (output) : " << hessianEntries.size() << std::endl;
#endif
		// Ex 1.2
		const double &x1 = x[0];
		const double &x2 = x[1];
		hessianEntries.push_back(Tripletd(0, 0, 2));
		hessianEntries.push_back(Tripletd(1, 1, 2));
	}

};

class Ex1Constraint : public FunctionConstraints {
public:
	// Ex 1
	
	// 1 2 | x1 = 3
	//     | x2 = 

	// Returns b of A(p) = b.
	virtual const VectorXd& getEqualityConstraintsTargetValues() {
		b.resize(1);
		b[0] = 3;

#ifdef DEBUG
		std::cout << "> SQP - getEqualityConstraintsTargetValues - size of b (input) : " << b.rows() << " by " << b.cols() << std::endl;
#endif
		return b;
	}

	// Returns A(p) of A(p) = b.
	// Derive from this to compute A(p). Fill `eqConstraintsVals` and return it.
	virtual const VectorXd& getEqualityConstraintValues(const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		eqConstraintVals.resize(1);
		eqConstraintVals[0] = x1 + 2 * x2;
		
#ifdef DEBUG
		std::cout << "> SQP - getEqualityConstraintValues - size of eqConstraintVals (input) : " << eqConstraintVals.rows() << " by " << eqConstraintVals.cols() << std::endl;
		std::cout << "> SQP - getEqualityConstraintValues - eqConstraintVals : " << eqConstraintVals << std::endl;
#endif
		return eqConstraintVals;
	}

	// Computes the Jacobian dA/dp of the equality constraints A.
	virtual void addEqualityConstraintsJacobianEntriesTo(std::vector<Tripletd>& jacobianEntries, const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		jacobianEntries.push_back(Tripletd(0, 0, 1)); // 2*x1 ? 
		jacobianEntries.push_back(Tripletd(0, 1, 2)); // 2*x2 ? 
	}

/* NOT USEFUL HERE :) 
	// Returns the value of the inequality constraint C(p).
	virtual const VectorXd& getInequalityConstraintValues(const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		ineqConstraintVals.resize(1);
		ineqConstraintVals[0] = x1 + 2 * x2;

		return ineqConstraintVals;
	}

	// Returns d of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMinValues() {
		d.resize(1);
		d[0] = 3;
		return d;
	}

	// Returns f of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMaxValues() {
		f.resize(1);
		f[0] = 3;
		return f;
	}

	// Computes the Jacobian dA/dp of the inequality constraints C.
	virtual void addInequalityConstraintsJacobianEntriesTo(std::vector<Tripletd>& jacobianEntries, const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		jacobianEntries.push_back(Tripletd(0, 0, 1));
		jacobianEntries.push_back(Tripletd(0, 1, 2));
	}
*/
	// Returns l of constraint l <= p <= u
	virtual const VectorXd& getBoundConstraintsMinValues() {
		l.resize(2);
		l[0] = -100; // ??? 
		l[1] = -100; // ??? 
		return l;
	}

	// Returns u of constraint l <= p <= u
	virtual const VectorXd& getBoundConstraintsMaxValues() {
		u.resize(2);
		u[0] = 100; // ??? 
		u[1] = 100; // ??? 
		return u;
	}

};

int main() {

	ObjectiveFunction* objective;
	FunctionConstraints* constraints;

	// Ex 1: uncomment these lines
	objective = new Ex1Objective;
	constraints = new Ex1Constraint;

	SQPFunctionMinimizer minimizer;
	VectorXd x(2); x << 100, -100;
	minimizer.minimize(objective, constraints, x);

	std::cout << "x              = " << x.transpose() << std::endl;
	std::cout << "f(x)           = " << objective->computeValue(x) << std::endl;
	std::cout << "c(x)           = " << constraints->getEqualityConstraintValues(x) << std::endl;
	std::cout << "# iterations   = " << minimizer.lastNumberIterations << std::endl;
}
