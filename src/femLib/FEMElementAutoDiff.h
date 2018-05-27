#pragma once

#include "FEMElement.h"

#include <AutoDiff.h>

template <class S> using Vector2 = Eigen::Matrix<S, 2, 1>;
template <class S> using Matrix2 = Eigen::Matrix<S, 2, 2>;

typedef AutoDiff<double, double> AD;
typedef AutoDiff<AD, AD> ADD;

/* This is the code from the tutorial using this AutoDiff implementation:
 *
 *		// df/dx1
		AD x1 = M_PI; x1.deriv() = 1;
		AD x2 = 3; x2.deriv() = 0;
		AD f = sin(x1*x2);
		double f_val = f.value();
		double dfdx1 = f.deriv();

		// ddf/dx1dx2
		ADD x1 = M_PI;
		x1.value().deriv() = 1; x1.deriv().value() = 0;
		ADD x2 = 3;
		x2.value().deriv() = 0; x1.deriv().value() = 1;
		ADD f = sin(x1*x2);
		double f_val = f.value().value();
		double ddfdx1dx2 = f.deriv().deriv();

		// ddf/dx2dx2
		ADD x1 = M_PI;
		x1.value().deriv() = 0; x1.deriv().value() = 0;
		ADD x2 = 3;
		x2.value().deriv() = 1; x1.deriv().value() = 1;
		ADD f = sin(x1*x2);
		double f_val = f.value().value();
		double ddfdx1dx2 = f.deriv().deriv();
*/

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElementAutoDiff : public FEMElement {
public:
	FEMElementAutoDiff(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: FEMElement(nodeIndices, X) {
	}
	~FEMElementAutoDiff(){
	}

	virtual double getEnergy(const VectorXd& x, const VectorXd& X){

		// Ex. 2.3

		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		return computeEnergyT({x0, x1, x2});
	}

	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {

		// Ex. 2.3

		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		Vector2<AD> xd[3] = {
			{x0[0], x0[1]},
			{x1[0], x1[1]},
			{x2[0], x2[1]}
		};

		for (int i = 0; i < 3; ++i) {
			for (int j= 0; j < 2; ++j) {
				xd[i][j].deriv() = 1;
				grad[2*getNodeIndex(i) + j] += computeEnergyT(xd).deriv();
				xd[i][j].deriv() = 0;
			}
		}
	}

	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {

		// Ex. 2.3

		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		Vector2<ADD> xdd[3] = {
			{x0[0], x0[1]},
			{x1[0], x1[1]},
			{x2[0], x2[1]}
		};

		for (int i1 = 0; i1 < 3; ++i1) {
			for (int j1 = 0; j1 < 2; ++j1) {
				xdd[i1][j1].deriv() = 1;
				for (int i2 = 0; i2 < 3; ++i2) {
					for (int j2 = 0; j2 < 2; ++j2) {
						xdd[i2][j2].value().deriv() = 1;
						hesEntries.push_back(Tripletd(2*getNodeIndex(i1)+j1, 2*getNodeIndex(i2)+j2, computeEnergyT(xdd).deriv().deriv()));
						xdd[i2][j2].value().deriv() = 0;
					}
				}
				xdd[i1][j1].deriv() = 0;
			}
		}

	}

private:
	// This is the same function as FEMElement::computeGradient,
	// however with a template parameter T.
	template<class T>
	void computeDeformationGradientT(const Vector2<T> (&x)[3], Matrix2<T>& dxdX) const {
		//edge vectors
		Vector2<T> v1 = x[1] - x[0];
		Vector2<T> v2 = x[2] - x[0];
		Matrix2<T> dx;
		dx << v1[0], v2[0],
				v1[1], v2[1];
		Matrix2<T> dxInvT;
		for (int i = 0; i < 4; ++i)
			dxInvT.data()[i] = dXInv.data()[i];
		dxdX = dx * dxInvT;
	}

	template<class T>
	T computeEnergyT(const Vector2<T> (&x)[3]) const {
		Matrix2<T> F;
		computeDeformationGradientT(x, F);

		T normF2 = F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0) + F(1,1)*F(1,1);
		T detF = F.determinant();
		T energyDensity = shearModulus/2 * (normF2-2) - shearModulus * log(detF) + bulkModulus/2 * log(detF) * log(detF);

		return energyDensity * restShapeArea;
	}
};
