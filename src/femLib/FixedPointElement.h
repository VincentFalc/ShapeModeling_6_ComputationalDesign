#pragma once

#include "Element.h"

// zero rest length spring connected to a target position (could be the mouse, or something else...)
class FixedPointElement : public Element{
public:
	FixedPointElement(int nodeIndex, const Vector2d &targetPosition, double K = 10)
		: nodeIndex(nodeIndex), targetPosition(targetPosition), k(K) {
	}
	~FixedPointElement(){
	}

	virtual int getNumNodes() const {
		return 1;
	}
	virtual int getNodeIndex(int i) const {
		assert(i == 0);
		return nodeIndex;
	}

	virtual double getMass() const {
		return 0;
	}

	virtual double getEnergy(const VectorXd& x, const VectorXd& X){
		Vector2d p = x.segment<2>(2*nodeIndex);

		double v0 = k;
		double v3 = 0.500000;
		double v4 = v3 * v0;
		double v5 = p[0];
		double v6 = targetPosition[0];
		double v7 = v5 - v6;
		double v8 = v7 * v7;
		double v9 = p[1];
		double v10 = targetPosition[1];
		double v11 = v9 - v10;
		double v12 = v11 * v11;
		double v13 = v8 + v12;
		double v14 = v4 * v13;
		return v14;
	}
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {
		Vector2d p = x.segment<2>(2*nodeIndex);

		double v2 = k;
		double v3 = p[0];
		double v4 = targetPosition[0];
		double v5 = v3 - v4;
		double v6 = v5 + v5;
		double v7 = 0.500000;
		double v8 = v7 * v2;
		double v9 = v6 * v8;
		double v10 = p[1];
		double v11 = targetPosition[1];
		double v12 = v10 - v11;
		double v13 = v12 + v12;
		double v14 = v13 * v8;
		grad[2*nodeIndex + 0] += v9;
		grad[2*nodeIndex + 1] += v14;
	}
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries){
		double v0 = k;
		double v1 = 2.000000;
		double v2 = 0.500000;
		double v3 = v2 * v0;
		double v4 = v1 * v3;
		double v5 = 0.000000;
		hesEntries.push_back(Tripletd(2*nodeIndex+0,2*nodeIndex+0,v4));
		hesEntries.push_back(Tripletd(2*nodeIndex+0,2*nodeIndex+1,v5));
		hesEntries.push_back(Tripletd(2*nodeIndex+1,2*nodeIndex+0,v5));
		hesEntries.push_back(Tripletd(2*nodeIndex+1,2*nodeIndex+1,v4));
	}

	void setTargetPosition(const Vector2d &p) {
		targetPosition = p;
	}

	Vector2d getTargetPosition() const {
		return targetPosition;
	}

private:
	double k;
	int nodeIndex;
	Vector2d targetPosition;
};
