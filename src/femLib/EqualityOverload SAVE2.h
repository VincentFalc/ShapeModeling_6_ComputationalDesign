#include <iostream>

#include <./../optLib/FunctionConstraints.h>
#include <./../optLib/SQPFunctionMinimizer.h>

#define DEBUG

class Ex2Constraint : public FunctionConstraints
{

  protected:
	SimulationMesh *simParent;
	const int NbLineToDisplay = 100;

  public:
	//Function to get the parent (Added)
	virtual const void setRef(SimulationMesh *simParent)
	{
		this->simParent = simParent;
	}

	virtual const double getAvgCoordsSupport()
	{
		//Calculus of xMax XMin
		double maxX = 0;
		double minX = 0;
		bool firstIt = true;

		//Get the min max of the support points
		for (std::set<int>::iterator it = this->simParent->supportNodes.begin(); it != this->simParent->supportNodes.end(); ++it)
		{
			int indice = *it;

#ifdef DEBUG
		std::cout << " Current support point : << " << indice << " of coords << " << this->simParent->x[2*indice] << " " << this->simParent->x[2*indice+1] << " of mass " <<  this->simParent->m[2*indice] << std::endl;
#endif
			//Vector2d currentSupportpoint = this->simParent->constraintNodes.find(indice)->second;
			double xPos = this->simParent->x[2*indice];

			//Instanciation for the first iteration
			if (firstIt)
			{
				maxX = xPos;
				minX = xPos;
				firstIt = false;
			}

			//Get the max and the min
			if (xPos > maxX)
			{
				maxX = xPos;
			}
			if (xPos < minX)
			{
				minX = xPos;
			}
		}
#ifdef DEBUG
		std::cout << " Max X coord : " << maxX << std::endl;
		std::cout << " Min X coord : " << minX << std::endl;
#endif

		double averageCoords = 0.5 * (maxX + minX);
#ifdef DEBUG
		std::cout << " Calculated average coords " << averageCoords << std::endl;
#endif
		//Final calculus of the right hand side constraint
		return averageCoords;
	}

	virtual const double getTotalMass()
	{
		double totalMass = 0;
		int nbMeshPoints = this->simParent->x.rows()/2;

#ifdef DEBUG
		std::cout << " mass matrix size : " <<  this->simParent->m.rows() << " rows and col = " << this->simParent->m.cols() << std::endl;
#endif

		//Get the min max of the support points
		for (int i = 0; i < nbMeshPoints; i++)
		{
			totalMass += this->simParent->m[2*i]; // mass of the point

#ifdef DEBUG
		std::cout << " Current support point : << " << i << " of coords " << this->simParent->x[2*i] << " \t " << this->simParent->x[2*i+1] << " \t\tof mass " <<  this->simParent->m[2*i] << " \tand " <<  this->simParent->m[2*i+1]  << std::endl;
#endif
		}

		//Sanity check
		assert(totalMass != 0); //Prevent division by zero
		return totalMass;
	}

	// Returns b of A(p) = b.
	virtual const VectorXd &getEqualityConstraintsTargetValues()
	{
		//Hypothèse 1 : Il faut un B qui soit aussi grand que deux fois la liste des points, pour avoir tous les X et tous les Y contraints.
		//Hypotèhse 2 : Il faut avoir un B qui soit aussi grand que la liste des points contraints *2 pour avoir tous les X et tous les Y.
		//Note : c'est stocké en (x1,y1),(x2,y2),...

		// Initialize
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		int nbMeshPoints = this->simParent->x.rows()/2;
		std::vector<Tripletd> coefficients;
		int sizeB = 2 * nbMeshPoints + 1;

		//Resize vector to good sizes
		b.resize(sizeB);

		//ADD "x = xi_traget" constraints
		for (std::map<int, Vector2d>::iterator it = this->simParent->constraintNodes.begin(); it != this->simParent->constraintNodes.end(); ++it)
		{
			int indice = it->first;

			coefficients.push_back(Tripletd(2 * indice, 0, it->second[0]));		//Store the x constrained value
			coefficients.push_back(Tripletd(2 * indice + 1, 0, it->second[1])); //Store the y constrained value
		}

		// ADD 0.5*(Xmax-xMin) constraint
		double averageCoords = getAvgCoordsSupport();

		//Storage of the last constraint
		coefficients.push_back(Tripletd(sizeB - 1, 0, averageCoords)); //Store the y constrained value in the last cell

		Eigen::SparseMatrix<double> Btmp(sizeB, 1);
		Eigen::VectorXd B(sizeB, 1);
		//Transform triplet into SparseMatrix
		Btmp.setFromTriplets(coefficients.begin(), coefficients.end());
		B = VectorXd(Btmp);

#ifdef DEBUG
		std::cout << " B Overview : " << std::endl;
		for (int i = 0; i < B.rows(); i++) //NbLineToDisplay first and last lines
		{
			if (i < NbLineToDisplay || i > B.rows() - NbLineToDisplay)
			{
				for (int j = 0; j < B.row(0).size(); j++)
				{
					std::cout << B.row(i)[j] << " ";
				}
				std::cout << "\n";
			}
		}
		std::cout << "> End of B Overview" << std::endl;
#endif
		return b;
	}

	virtual const std::vector<Tripletd> getATriplet(std::vector<Tripletd> &coefficients, const VectorXd &p)
	{

		//Sanity check
		assert(p.rows() % 2 == 0); //All points are 2D (X followed by Y component)
		coefficients.clear();

		// Initialize
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		int nbMeshPoints = p.rows() / 2;

		//ADD "x = xi_traget" left values
		for (std::map<int, Vector2d>::iterator it = this->simParent->constraintNodes.begin(); it != this->simParent->constraintNodes.end(); ++it)
		{
			int indice = it->first;
			coefficients.push_back(Tripletd(2 * indice, 2 * indice, 1));		 //xi/dxi = 1 in the "good" cell
			coefficients.push_back(Tripletd(2 * indice + 1, 2 * indice + 1, 1)); //xi/dxi = 1 in the "good" cell
		}

		//Get the Total mass of the system (for division)
		double totalMass = getTotalMass();

#ifdef DEBUG
		std::cout << " Mass of points :  " << std::endl;
#endif
		//ADD XCom(x)
		for (int i = 0; i < nbMeshPoints; i++)
		{
			coefficients.push_back(Tripletd(nbMeshPoints * 2, 2 * i, this->simParent->m[2*i] / totalMass)); //(mi)/M attributed to x components
#ifdef DEBUG
		std::cout << this->simParent->m[2*i] ;
#endif
		}
#ifdef DEBUG
		std::cout << std::endl;
#endif
		return coefficients;
	}

	virtual const Eigen::SparseMatrix<double> getA(const VectorXd &p)
	{

		//Sanity check
		assert(p.rows() % 2 == 0); //All points are 2D (X followed by Y component)

		// Initialize
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		int nbMeshPoints = p.rows() / 2;
		int currIndice = 0;
		std::vector<Tripletd> coefficients;
		Eigen::SparseMatrix<double> A(nbMeshPoints * 2 + 1, nbMeshPoints * 2); //One line more for the mass constraint

		//Transform triplet into SparseMatrix
		getATriplet(coefficients, p); //Get the A List of triplet
		A.setFromTriplets(coefficients.begin(), coefficients.end());

		return A;
	}

	// Returns A(p) of A(p) = b.
	// Derive from this to compute A(p). Fill `eqConstraintsVals` and return it.
	virtual const VectorXd &getEqualityConstraintValues(const VectorXd &p)
	{

		Eigen::SparseMatrix<double> A = getA(p);
		eqConstraintVals.resize(A.rows());
		eqConstraintVals = A * p; //Compute A(p)

#ifdef DEBUG

		Eigen::MatrixXd printableA = MatrixXd(A);

		std::cout << " System A * p = equalityConstraints Overview : " << std::endl;
		for (int i = 0; i < eqConstraintVals.rows(); i++) //NbLineToDisplay first and last lines
		{
			if (i < NbLineToDisplay || i > eqConstraintVals.rows() - NbLineToDisplay)
			{
				for (int j = 0; j < printableA.row(0).size(); j++)
				{
					std::cout << printableA.row(i)[j] << " ";
				}
				std::cout << " * ";
				if(i < p.rows()){
					for (int j = 0; j < p.row(0).size(); j++)
					{
						std::cout << p.row(i)[j] << " ";
					}
				}
				std::cout << "\t = ";
				for (int j = 0; j < eqConstraintVals.row(0).size(); j++)
				{
					std::cout << eqConstraintVals.row(i)[j] << " ";
				}
				std::cout << "\n";
			}
		}
		std::cout << "> End of eqConstraintVals Overview" << std::endl;
#endif

		return eqConstraintVals;
	}

	// Computes the Jacobian dA/dp of the equality constraints A.
	virtual void addEqualityConstraintsJacobianEntriesTo(std::vector<Tripletd> &jacobianEntries, const VectorXd &p)
	{
		getATriplet(jacobianEntries, p);

	}

	// Returns l of constraint l <= p <= u
	virtual const VectorXd &getBoundConstraintsMinValues()
	{
		assert(this->simParent->x.rows() % 2 == 0); //There is 2 component per node, xi followed by yi

		//ALL POINTS
		int nbMeshPoints = this->simParent->x.rows() / 2;
		l.resize(nbMeshPoints * 2);

		for (int i = 0; i < nbMeshPoints; i++)
		{
			l[2 * i] = this->simParent->x[2 * i] - 1000;
			l[2 * i + 1] = this->simParent->x[2 * i + 1] - 1000;
		}

		return l;
	}

	// Returns u of constraint l <= p <= u
	virtual const VectorXd &getBoundConstraintsMaxValues()
	{
		assert(this->simParent->x.rows() % 2 == 0); //There is 2 component per node, xi followed by yi

		//ALL POINTS
		int nbMeshPoints = this->simParent->x.rows() / 2;
		u.resize(nbMeshPoints * 2);

		std::cout << nbMeshPoints << " nbMeshPoints " << std::endl;

		for (int i = 0; i < nbMeshPoints; i++)
		{
			u[2 * i] = this->simParent->x[2 * i] + 1000;
			u[2 * i + 1] = this->simParent->x[2 * i + 1] + 1000;
		}

		return u;
	}
};