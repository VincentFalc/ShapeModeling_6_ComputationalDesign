#pragma once

//#pragma warning( disable : 4996)

#define EPSILON 1e-10
#define TINY 0.0000001
#define IS_ZERO(x) (fabs(x)<EPSILON)

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Eigen::VectorXd;
using Eigen::MatrixXd;
typedef Eigen::SparseMatrix<double> SparseMatrixd;
typedef Eigen::Triplet<double> Tripletd;

#ifdef _DEBUG // DEBUG
typedef Eigen::Matrix<double, 2, 2, Eigen::DontAlign> Matrix2d;
typedef Eigen::Matrix<double, 2, 1, Eigen::DontAlign> Vector2d;
#else // RELEASE
using Eigen::Matrix2d;
using Eigen::Vector2d;
#endif // _DEBUG

inline void resize(SparseMatrixd& sm, int rows, int cols) {
	if (sm.rows() != rows || sm.cols() != cols)
		sm.resize(rows, cols);
	sm.setZero();
}

inline void resize(VectorXd& v, int n) {
	if (v.size() != n)
		v.resize(n);
	v.setZero();
}

inline void resize(MatrixXd& m, int rows, int cols) {
	if (m.rows() != rows || m.cols() != cols)
		m.resize(rows, cols);
	m.setZero();
}

template<class MATType>
void addSparseMatrixDenseBlockToTriplet(std::vector<Tripletd>& triplets, int startX, int startY, const MATType& block, bool writeOnlyLowerDiagonalValues = false) {
    for (int i = 0; i < block.rows(); i++)
        for (int j = 0; j < block.cols(); j++)
        if (startX + i >= startY + j || !writeOnlyLowerDiagonalValues)
			triplets.push_back(Tripletd(startX + i, startY + j, block(i, j)));
}

inline double cross2d(const Vector2d &a, const Vector2d &b) {
	return a.x()*b.y() - b.x()*a.y();
}
