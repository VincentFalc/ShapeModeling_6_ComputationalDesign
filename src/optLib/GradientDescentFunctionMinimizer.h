#pragma once

#include "ObjectiveFunction.h"

#include <string>
#include <iostream>
#include <cmath>
#include <cfloat>

class GradientDescentFunctionMinimizer{
public:
	GradientDescentFunctionMinimizer(int maxIterations=100, double solveResidual=1e-5, int maxLineSearchIterations=15)
		: maxIterations(maxIterations), solveResidual(solveResidual), maxLineSearchIterations(maxLineSearchIterations){
	}

	virtual ~GradientDescentFunctionMinimizer(){}

	int getLastIterations() { return lastIterations; }

	virtual bool minimize(ObjectiveFunction *function, VectorXd &x){

		//number of parameters...
		int N = (int) x.size();
		resize(xi, N);
		resize(dx, N);
		resize(gradient, N);

		xi = x;

		bool optimizationConverged = false;

		int i=0;
		for(; i < maxIterations; i++) {
			computeSearchDirection(function, xi, dx);

			if (dx.norm() < solveResidual){
				optimizationConverged = true;
				break;
			}

			doLineSearch(function, dx, xi);
		}

		lastIterations = i;

		//p now holds the parameter values at the start of the iteration...
		x = xi;

		//and done!
		return optimizationConverged;
	}

protected:
	// Since the gradient of a function gives the direction of steepest descent, all one needs to do is go in that direction...
	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd& dx) {

		// Ex. 1.1

		dx.setZero();
		function->addGradientTo(dx, x);
	}

	virtual void doLineSearch(ObjectiveFunction *function, const VectorXd& dx, VectorXd& xi)
	{
		// Ex. 1.1

		// line search now...
		double alpha = 1.0;
		VectorXd xc(xi);
		double initialValue = function->computeValue(xc);

		for(int j = 0; j < maxLineSearchIterations; j++) {
			// try a new solution
			xi = xc - dx * alpha;

			// now check the new function value at this point...
			double newLineSearchValue = function->computeValue(xi);

			if((!std::isfinite(newLineSearchValue) || newLineSearchValue > initialValue)
					&& j < maxLineSearchIterations -1)
				// restore and try again...
				alpha /= 2.0;
			else
				// found a better solution!
				return;
		}

		// couldn't find a good value. Return what we now have and hope for the best...
		std::cout << "line search failed." << std::endl;
	}

protected:
	double solveResidual = 1e-5;
	int maxIterations = 100;
	int maxLineSearchIterations = 15;

	VectorXd xi, dx, gradient;

	// some stats about the last time `minimize` was called
	int lastIterations = -1;
};
