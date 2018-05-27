#pragma once

#include "Element.h"

#include <array>

enum MaterialModel2D {
	MM_LINEAR_ISOTROPIC=0,
	MM_STVK,
	MM_NEO_HOOKEAN
};

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElement : public Element {
public:
	FEMElement(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: nodeIndices(nodeIndices) {

		assert(nodeIndices.size() == 3);

		precomputeFromUndeformedState(X);
	}
	~FEMElement(){
	}

	virtual int getNumNodes() const {
		return 3;
	}
	virtual int getNodeIndex(int i) const {
		return nodeIndices[i];
	}
	virtual double getMass() const {
		return restShapeArea * massDensity;
	}

protected:
	// sets important properties of the rest shape using the set of points passed in as parameters
	void precomputeFromUndeformedState(const VectorXd &X){
		//edge vectors
		Vector2d v1 = getNodePos(1, X) - getNodePos(0, X);
		Vector2d v2 = getNodePos(2, X) - getNodePos(0, X);

		//matrix that holds three edge vectors
		Matrix2d dX;
		dX << v1[0], v2[0],
				v1[1], v2[1];

		dXInv = dX.inverse();

		//compute the area of the element...
		restShapeArea = 1 / 2.0 * fabs(cross2d(v1, v2));
	}

	// as a deformation measure, we need to compute the deformation gradient F. F maps deformed vectors dx to undeformed coords dX: dx = F*dX.
	// for linear basis functions, an easy way to compute it is by looking at the matrix that maps deformed traingle/tet edges to their underformed counterparts (F = dx * inv(dX)).
	void computeDeformationGradient(const Vector2d (&x)[3], Matrix2d& dxdX) const {
		//edge vectors
		Vector2d v1 = x[1] - x[0];
		Vector2d v2 = x[2] - x[0];
		Matrix2d dx;
		dx << v1[0], v2[0],
				v1[1], v2[1];
		dxdX = dx * dXInv;
	}

//protected:
public:
	// the collection of nodes that define the triangle element
	std::array<int, 3> nodeIndices;
	// material parameters
	double shearModulus = 50, bulkModulus = 500;
	// relates area/volume to the mass of the element
	double massDensity = 1;
	// precomputed values
	double restShapeArea = 0;
	Matrix2d dXInv;
};
