#pragma once

#include "ObjectiveFunction.h"

class RosenbrockFunction : public ObjectiveFunction {
public:

	RosenbrockFunction() {
	a = 1; b = 100;
	}

	virtual double computeValue(const VectorXd& x) {

		// Ex 1.1

		const double &x1 = x[0];
		const double &x2 = x[1];
		return std::pow(a-x1,2.0) + b*std::pow(x2-x1*x1, 2.0);
	}

	virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {

		// Ex 1.1

		const double &x1 = x[0];
		const double &x2 = x[1];
		grad[0] += -2*(a-x1) - 4*b*(x2-x1*x1)*x1;
		grad[1] += 2*b*(x2-x1*x1);
	}


	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {

		// Ex 1.2

		const double &x1 = x[0];
		const double &x2 = x[1];
		hessianEntries.push_back(Tripletd(0, 0, 2 + 12*b*x1*x1 - 4*b*x2));
		hessianEntries.push_back(Tripletd(1, 0, -4*b*x1));
		hessianEntries.push_back(Tripletd(1, 1, 2*b));
	}

	double a, b;
};
