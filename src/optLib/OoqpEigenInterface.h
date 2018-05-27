#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "MathHelper.h"

namespace ooqpei {

class OoqpEigenInterface{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*!
   * Solve min 1/2 x' Q x + c' x, such that A x = b, d <= Cx <= f, and l <= x <= u.
   * @param [in] Q a symmetric positive semidefinite matrix (nxn)
   * @param [in] c a vector (nx1)
   * @param [in] A a (possibly null) matrices (m_axn)
   * @param [in] b a vector (m_ax1)
   * @param [in] C a (possibly null) matrices (m_cxn)
   * @param [in] d a vector (m_cx1)
   * @param [in] f a vector (m_cx1)
   * @param [in] l a vector (nx1)
   * @param [in] u a vector (nx1)
   * @param [out] x a vector of variables (nx1)
   * @return true if successful
   */
  static bool solve(const SparseMatrixd& Q,
                    const VectorXd& c,
                    const SparseMatrixd& A,
                    const VectorXd& b,
                    const SparseMatrixd& C,
                    const VectorXd& d, const VectorXd& f,
                    const VectorXd& l, const VectorXd& u,
                    VectorXd& x);

  /*!
   * Solve min 1/2 x' Q x + c' x, such that A x = b, and d <= Cx <= f
   * @param [in] Q a symmetric positive semidefinite matrix (nxn)
   * @param [in] c a vector (nx1)
   * @param [in] A a (possibly null) matrices (m_axn)
   * @param [in] b a vector (m_ax1)
   * @param [in] C a (possibly null) matrices (m_cxn)
   * @param [in] d a vector (m_cx1)
   * @param [in] f a vector (m_cx1)
   * @param [out] x a vector of variables (nx1)
   * @return true if successful
   */
  static bool solve(const SparseMatrixd& Q,
                    VectorXd& c,
                    const SparseMatrixd& A,
                    VectorXd& b,
                    const SparseMatrixd& C,
                    VectorXd& d, VectorXd& f,
                    VectorXd& x);

  /*!
   * Solve min 1/2 x' Q x + c' x, such that A x = b, and l <= x <= u.
   * @param [in] Q a symmetric positive semidefinite matrix (nxn)
   * @param [in] c a vector (nx1)
   * @param [in] A a (possibly null) matrices (m_axn)
   * @param [in] b a vector (m_ax1)
   * @param [in] l a vector (nx1)
   * @param [in] u a vector (nx1)
   * @param [out] x a vector of variables (nx1)
   * @return true if successful
   */
  static bool solve(const SparseMatrixd& Q,
                    VectorXd& c,
                    const SparseMatrixd& A,
                    VectorXd& b,
                    VectorXd& l, VectorXd& u,
                    VectorXd& x);

  /*!
   * Solve min 1/2 x' Q x + c' x, such that Cx <= f
   * @param [in] Q a symmetric positive semidefinite matrix (nxn)
   * @param [in] c a vector (nx1)
   * @param [in] C a (possibly null) matrices (m_cxn)
   * @param [in] f a vector (m_cx1)
   * @param [out] x a vector of variables (nx1)
   * @return true if successful
   */
  static bool solve(const SparseMatrixd& Q,
                    VectorXd& c,
                    const SparseMatrixd& C,
                    VectorXd& f,
                    VectorXd& x);

  /*!
   * Solve min 1/2 x' Q x + c' x
   * @param [in] Q a symmetric positive semidefinite matrix (nxn)
   * @param [in] c a vector (nx1)
   * @param [out] x a vector of variables (nx1)
   * @return true if successful
   */
  static bool solve(const SparseMatrixd& Q,
                    VectorXd& c,
                    VectorXd& x);

  /*!
   * Change to true to print debug information.
   * @return true if in debug mode
   */
  static bool isInDebugMode() { return isInDebugMode_; };
  static void setIsInDebugMode(bool isInDebugMode) {
    isInDebugMode_ = isInDebugMode;
  }

 private:
  /*!
   * Determine which limits are active and which are not.
   * @param [in]  l
   * @param [in]  u
   * @param [out] useLowerLimit
   * @param [out] useUpperLimit
   * @param [out] lowerLimit
   * @param [out] upperLimit
   */
  static void generateLimits(const VectorXd& l, const VectorXd& u,
                      Eigen::Matrix<char, Eigen::Dynamic, 1>& useLowerLimit,
                      Eigen::Matrix<char, Eigen::Dynamic, 1>& useUpperLimit,
                      VectorXd& lowerLimit, VectorXd& upperLimit, bool ignoreEqualLowerAndUpperLimits);

  static void printProblemFormulation(
      const SparseMatrixd& Q, const VectorXd& c,
      const SparseMatrixd& A, const VectorXd& b,
      const SparseMatrixd& C, const VectorXd& d, const VectorXd& f,
          const VectorXd& l, const VectorXd& u);

  static void printLimits(Eigen::Matrix<char, Eigen::Dynamic, 1>& useLowerLimit,
                          Eigen::Matrix<char, Eigen::Dynamic, 1>& useUpperLimit,
                          VectorXd& lowerLimit,
                          VectorXd& upperLimit);

  static void printSolution(int& status, VectorXd& x);

 private:
  static bool isInDebugMode_;
};

} /* namespace ooqpei */
