#pragma once

#include "FEMElement.h"

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElementFD : public FEMElement {
public:
	FEMElementFD(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: FEMElement(nodeIndices, X) {
	}
	~FEMElementFD(){
	}

	virtual double getEnergy(const VectorXd& x, const VectorXd& X)
	{
		// Ex 2.1

		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		return computeEnergy({x0, x1, x2});
	}
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad)
	{
		// Ex 2.1

		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		Vector2d dEdx[3];
		computeGradientComponents({x0, x1, x2}, dEdx);

		for (int i = 0;i<3;i++)
			for (int j = 0;j<2;j++)
				grad[2*nodeIndices[i] + j] += dEdx[i][j];
	}

	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries)
	{
		// Ex 2.2

		// get vertices of this element
		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		//compute the hessian blocks and write in global hessian
		Matrix2d ddEdxdx[3][3];
		computeHessianComponents({x0, x1, x2}, ddEdxdx);
		for (int i = 0;i<3;i++)
			for (int j = 0;j < 3;j++)
				addSparseMatrixDenseBlockToTriplet(hesEntries, 2*nodeIndices[i], 2*nodeIndices[j], ddEdxdx[i][j], true);
	}

private:

	double computeEnergy(const Vector2d (&x)[3]) const
	{
		Matrix2d F;
		computeDeformationGradient(x, F);

		double normF2 = F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0) + F(1,1)*F(1,1);
		double detF = F.determinant();

		double energyDensity = shearModulus/2 * (normF2-2) - shearModulus * log(detF) + bulkModulus/2 * log(detF) * log(detF);

		return energyDensity * restShapeArea;
	}

	void computeGradientComponents(const Vector2d (&x)[3], Vector2d (&dEdx)[3]) const
	{
		// loop over 3 vertices
		for (int i = 0; i < 3; ++i) {
			// loop over dimensions
			for (int j = 0; j < 2; ++j) {
				Vector2d xp[3] = {x[0], x[1], x[2]};
				Vector2d xm[3] = {x[0], x[1], x[2]};
				xp[i][j] += h;
				xm[i][j] -= h;
				double ep = computeEnergy(xp);
				double em = computeEnergy(xm);
				dEdx[i][j] = (ep-em)/(2*h);
			}
		}
	}

	void computeHessianComponents(const Vector2d (&x)[3], Matrix2d (&ddEdxdx)[3][3]) const
	{
		// loop over 3 vertices
		for (int i = 0; i < 3; ++i) {

			// loop over dimensions
			for (int j = 0; j < 2; ++j) {
				Vector2d xp[3] = {x[0], x[1], x[2]};
				Vector2d xm[3] = {x[0], x[1], x[2]};
				xp[i][j] += h;
				xm[i][j] -= h;

				Vector2d dEdx_p[3];
				computeGradientComponents(xp, dEdx_p);

				Vector2d dEdx_m[3];
				computeGradientComponents(xm, dEdx_m);

				for (int k = 0; k < 3; ++k) {
					for (int l = 0; l < 2; ++l) {
						ddEdxdx[i][k](j,l) = (dEdx_p[k](l) - dEdx_m[k](l)) / (2*h);
					}
				}
			}
		}
	}

private:
	double h = 1e-8; // step size
};
