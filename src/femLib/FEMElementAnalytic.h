#pragma once

#include "FEMElement.h"

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElementAnalytic : public FEMElement {
public:
	FEMElementAnalytic(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: FEMElement(nodeIndices, X) {
	}
	~FEMElementAnalytic(){
	}

	virtual double getEnergy(const VectorXd& x, const VectorXd& X){
		//compute the deformation gradient
		Matrix2d F;
		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);
		computeDeformationGradient({x0, x1, x2}, F);

		double energyDensity = 0;

		double normF2 = F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0) + F(1,1)*F(1,1);
		double detF = F.determinant();
		energyDensity += shearModulus/2 * (normF2-2) - shearModulus * log(detF) + bulkModulus/2 * log(detF) * log(detF);

		return energyDensity * restShapeArea;
	}

	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad) {
		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		//compute the local gradient, and write it into global gradient
		Vector2d dEdx[3];
		computeGradientComponents({x0,x1,x2}, dEdx);
		for (int i = 0;i<3;i++)
			for (int j = 0;j<2;j++)
				grad[2*nodeIndices[i] + j] += dEdx[i][j];
	}

	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries) {
		Vector2d x0 = getNodePos(0, x);
		Vector2d x1 = getNodePos(1, x);
		Vector2d x2 = getNodePos(2, x);

		//compute the hessian blocks and write in global hessian
		Matrix2d ddEdxdx[3][3];
		computeHessianComponents({x0,x1,x2}, ddEdxdx);
		for (int i = 0;i<3;i++)
			for (int j = 0;j < 3;j++)
				addSparseMatrixDenseBlockToTriplet(hesEntries, 2*nodeIndices[i], 2*nodeIndices[j], ddEdxdx[i][j], true);
	}

private:
	void computeGradientComponents(const Vector2d (&x)[3], Vector2d (&dEdx)[3]){
		//compute the gradient of the energy using the chain rule: dE/dx = dE/dF * dF/dx. dE/dF is the first Piola-Kirchoff stress sensor, for which nice expressions exist.

		//compute the deformation gradient
		Matrix2d F;
		computeDeformationGradient(x, F);

		Matrix2d dEdF;
		Matrix2d FinvT = F.inverse().transpose();

		double normF2 = F(0,0)*F(0,0) + F(0,1)*F(0,1) + F(1,0)*F(1,0) + F(1,1)*F(1,1);
		double detF = F.determinant();
		dEdF = F * shearModulus + FinvT * (-shearModulus + bulkModulus*log(detF));

		//dF/dx is going to be some +/- Xinv terms. The forces on nodes 1,2 can be writen as: dE/dF * XInv', while the force on node 0 is -f1-f2;
		dEdx[1] = Vector2d(dEdF(0,0) * dXInv(0,0) + dEdF(0,1) * dXInv(0,1), dEdF(1,0) * dXInv(0,0) + dEdF(1,1) * dXInv(0,1)) * restShapeArea;
		dEdx[2] = Vector2d(dEdF(0,0) * dXInv(1,0) + dEdF(0,1) * dXInv(1,1), dEdF(1,0) * dXInv(1,0) + dEdF(1,1) * dXInv(1,1)) * restShapeArea;
		dEdx[0] = -dEdx[1]-dEdx[2];
	}
	void computeHessianComponents(const Vector2d (&x)[3], Matrix2d (&ddEdxdx)[3][3]) {

		Matrix2d F;
		computeDeformationGradient(x, F);

		Matrix2d Finv = F.inverse();
		Matrix2d FinvT = Finv.transpose();
		Matrix2d dF, dP, tmpM, dH;
		const double dDs[6][4] = { { -1,-1,0,0 },{ 0,0,-1,-1 },{ 1,0,0,0 },{ 0,0,1,0 },{ 0,1,0,0 },{ 0,0,0,1 } };
		for (int i = 0; i < 6; ++i)
		{
			for (int x = 0; x < 4; ++x)
				dF(x / 2, x % 2) = dDs[i][x];
			dF = dF * dXInv;
			dP = shearModulus * dF;
			double J = F.determinant();
			dP = dP + (shearModulus - bulkModulus * log(J)) * FinvT * dF.transpose() * FinvT;
			tmpM = Finv * dF;
			dP = dP + bulkModulus * (tmpM(0, 0) + tmpM(1, 1)) * FinvT;
			dH = restShapeArea * dP * dXInv.transpose();
			int row = i / 2, subrow = i % 2;
			ddEdxdx[row][1](subrow, 0) = dH(0, 0); ddEdxdx[row][1](subrow, 1) = dH(1, 0);
			ddEdxdx[row][2](subrow, 0) = dH(0, 1); ddEdxdx[row][2](subrow, 1) = dH(1, 1);
			ddEdxdx[row][0](subrow, 0) = -(dH(0, 0) + dH(0, 1));
			ddEdxdx[row][0](subrow, 1) = -(dH(1, 0) + dH(1, 1));
		}
	}
};
