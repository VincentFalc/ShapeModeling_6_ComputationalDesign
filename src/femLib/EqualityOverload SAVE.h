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
			Vector2d currentSupportpoint = this->simParent->constraintNodes.find(indice)->second;

			//Instanciation for the first iteration
			if (firstIt)
			{
				maxX = currentSupportpoint[0];
				minX = currentSupportpoint[0];
				firstIt = false;
			}

			//Get the max and the min
			if (currentSupportpoint[0] > maxX)
			{
				maxX = currentSupportpoint[0];
			}
			if (currentSupportpoint[0] < minX)
			{
				minX = currentSupportpoint[0];
			}
		}

		//Final calculus of the right hand side constraint
		return 0.5 * (maxX - minX);
	}

	virtual const double getTotalMass()
	{
		double totalMass = 0;
		int nbMeshPoints = this->simParent->x.rows();

		//Get the min max of the support points
		for (int i = 0; i < nbMeshPoints; i++)
		{
			totalMass += this->simParent->m[i]; // mass of the point
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
/*
		// Initialize
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		int currIndice = 0 ;

		//Resize vector to good sizes
		b.resize(nbConstrainedPoints*2+1);

		//ADD "x = xi_traget" constraints
		for (std::map<int,Vector2d>::iterator it=this->simParent->constraintNodes.begin(); it!=this->simParent->constraintNodes.end(); ++it){
			int indice = it->first;
			
			b[currIndice] = it->second[0];
			currIndice++;
			b[currIndice] = it->second[1];
			currIndice++;
		}

		// ADD 0.5*(Xmax-xMin) constraint
		//Calculus of xMax XMin
		double maxX = 0;
		double minX = 0;
		bool firstIt = true;

		//Get the min max of the support points
		for (std::set<int>::iterator it=this->simParent->supportNodes.begin(); it!=this->simParent->supportNodes.end(); ++it){
			int indice = *it;
			Vector2d currentSupportpoint = this->simParent->constraintNodes.find(indice)->second;

			//Instanciation for the first iteration
			if(firstIt){
				maxX = currentSupportpoint[0];
				minX = currentSupportpoint[0];
				firstIt = false;
			}

			//Get the max and the min
			if(currentSupportpoint[0]>maxX){
				maxX = currentSupportpoint[0];
			}
			if(currentSupportpoint[0]<minX){
				minX = currentSupportpoint[0];
			}

		}

		//Final calculus of the right hand side constraint
		double averageCoords = 0.5 * (maxX - minX);

		//Storage of the last constraint
		b[currIndice] = averageCoords;

		//Sanity check
		assert(nbConstrainedPoints*2 == currIndice); //Not +1 because we start at 0. So the last cell is size-1
*/
#ifdef DEBUG
		std::cout << "> SQP - getEqualityConstraintsTargetValues - size of b (input) : " << b.rows() << " by " << b.cols() << std::endl;
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

		//ADD XCom(x)
		for (int i = 0; i < nbMeshPoints; i++)
		{
			coefficients.push_back(Tripletd(nbMeshPoints * 2, 2 * i, this->simParent->m[i] / totalMass)); //(mi)/M attributed to x components
		}

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
		std::cout << " eqConstraintVals Overview : " << std::endl;
		for (int i = 0; i < eqConstraintVals.rows(); i++) //NbLineToDisplay first and last lines
		{
			if (i < NbLineToDisplay || i > eqConstraintVals.rows() - NbLineToDisplay)
			{
				for (int j = 0; j < eqConstraintVals.row(0).size(); j++)
				{
					std::cout << eqConstraintVals.row(i)[j] << " ";
				}
				std::cout << "\n";
			}
		}
		std::cout << "> End of eqConstraintVals Overview" << std::endl;
#endif

		/*
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		int nbMeshPoints = p.rows()/2;
		eqConstraintVals.resize(nbConstrainedPoints*2+1);
		int currIndice = 0 ;

		//ADD "x = xi_traget" left values
		for (std::map<int,Vector2d>::iterator it=this->simParent->constraintNodes.begin(); it!=this->simParent->constraintNodes.end(); ++it){
			int indice = it->first;
			assert(indice+nbMeshPoints < p.rows());

			eqConstraintVals[currIndice] = p[indice];
			eqConstraintVals[currIndice+nbConstrainedPoints] = p[indice+nbMeshPoints];
			currIndice++;
		}

		currIndice = (currIndice) * 2;

		//ADD XCom(x) ? 
		double totalMass = 0;
		double sum = 0;

		//Get the min max of the support points
		for (int i = 0 ; i < nbMeshPoints ; i ++){
			totalMass += this->simParent->m[i]; // mass of the point
			sum += this->simParent->m[i] * p[i]; // mass of the point * coordX
		}

		double xCOM = (1/totalMass) * sum;
		//currIndice++;
		eqConstraintVals[currIndice] = xCOM;

		assert(nbConstrainedPoints*2 == currIndice);

	#ifdef DEBUG
		std::cout << "> SQP - getEqualityConstraintValues - size of eqConstraintVals (input) : " << eqConstraintVals.rows() << " by " << eqConstraintVals.cols() << std::endl;
		std::cout << "> SQP - getEqualityConstraintValues - eqConstraintVals : " << eqConstraintVals << std::endl;
		std::cout << "> SQP - getEqualityConstraintValues - size of p (input) : " << p.rows() << " by " << p.cols() << std::endl;
		std::cout << "> SQP - getEqualityConstraintValues - p : " << p << std::endl;
#endif
*/
		return eqConstraintVals;
	}

	// Computes the Jacobian dA/dp of the equality constraints A.
	virtual void addEqualityConstraintsJacobianEntriesTo(std::vector<Tripletd> &jacobianEntries, const VectorXd &p)
	{
		getATriplet(jacobianEntries, p);

#ifdef DEBUG
		int nbMeshPoints = p.rows() / 2;
		Eigen::SparseMatrix<double> A(nbMeshPoints * 2 + 1, nbMeshPoints * 2); //One line more for the mass constraint
		A.setFromTriplets(jacobianEntries.begin(), jacobianEntries.end());

		std::cout << " jacobianEntries Overview : " << std::endl;
		Eigen::MatrixXd printableM = MatrixXd(A);
		for (int i = 0; i < printableM.rows(); i++) //NbLineToDisplay first and last lines
		{
			if (i < NbLineToDisplay || i > printableM.rows() - NbLineToDisplay)
			{
				for (int j = 0; j < printableM.row(0).size(); j++)
				{
					std::cout << printableM.row(i)[j] << " ";
				}
				std::cout << "\n";
			}
		}
		std::cout << "> End of jacobianEntries Overview" << std::endl;
#endif

		/*
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		int nbMeshPoints = p.rows()/2;
		int currIndice = 0 ;

		//ADD "x = xi_traget" left values
		for (std::map<int,Vector2d>::iterator it=this->simParent->constraintNodes.begin(); it!=this->simParent->constraintNodes.end(); ++it){
			int indice = it->first;
			jacobianEntries.push_back(Tripletd(currIndice, indice, 1)); //xi/dxi
			currIndice++;
		}

		currIndice = (currIndice) * 2;
		//currIndice++;

		double totalMass = 0;

		//Get the min max of the support points
		for (int i = 0 ; i < nbMeshPoints ; i ++){
			totalMass += this->simParent->m[i]; // mass of the point
		}

		//ADD XCom(x)
		for (int i = 0 ; i < nbMeshPoints ; i ++){
			jacobianEntries.push_back(Tripletd(currIndice, i, this->simParent->m[i]/totalMass)); //(mi)/M
		}

		assert(nbConstrainedPoints*2 == currIndice);
	*/
	}
	/*
	// Returns the value of the inequality constraint C(p).
	virtual const VectorXd& getInequalityConstraintValues(const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		ineqConstraintVals.resize(1);
		ineqConstraintVals[0] = x1 + 2 * x2;

		return ineqConstraintVals;
	}

	// Returns d of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMinValues() {
		d.resize(1);
		d[0] = 3;
		return d;
	}

	// Returns f of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMaxValues() {
		f.resize(1);
		f[0] = 3;
		return f;
	}

	// Computes the Jacobian dA/dp of the inequality constraints C.
	virtual void addInequalityConstraintsJacobianEntriesTo(std::vector<Tripletd>& jacobianEntries, const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		jacobianEntries.push_back(Tripletd(0, 0, 1));
		jacobianEntries.push_back(Tripletd(0, 1, 2));
	}
*/
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

		/*
	// CONSTRAINTED POINTS ONLY
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		l.resize(nbConstrainedPoints*2);
		int currIndice = 0 ;

		for (std::map<int,Vector2d>::iterator it=this->simParent->constraintNodes.begin(); it!=this->simParent->constraintNodes.end(); ++it){
			int indice = it->first;

			l[currIndice] = it->second[0]-1000;
			l[currIndice+nbConstrainedPoints] = it->second[1]-1000;
			currIndice++;
		}
*/
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

		/*
		// CONSTRAINTED POINTS ONLY
		int nbConstrainedPoints = this->simParent->constraintNodes.size();
		u.resize(nbConstrainedPoints*2);
		int currIndice = 0 ;

		for (std::map<int,Vector2d>::iterator it=this->simParent->constraintNodes.begin(); it!=this->simParent->constraintNodes.end(); ++it){
			int indice = it->first;

			u[currIndice] = it->second[0]+1000;
			u[currIndice+nbConstrainedPoints] = it->second[1]+1000;
			currIndice++;
		} */
		return u;
	}
};

/*
class Ex2Objective : public ObjectiveFunction {
public:

	virtual double computeValue(const VectorXd& x){
#ifdef DEBUG
		std::cout << "> SQP - computeValue - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
#endif
		const double &x1 = x[0];
		const double &x2 = x[1];
		return std::pow(x1,2.0) + std::pow(x2, 2.0);
	}

	virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {
#ifdef DEBUG
		std::cout << "> SQP - addGradientTo - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> SQP - addGradientTo - size of grad (output) : " << grad.rows() << " by " << grad.cols() << std::endl;
#endif
		// Ex 1.1
		const double &x1 = x[0];
		const double &x2 = x[1];
		//grad.resize(2);
		grad[0] += 2*x1;
		grad[1] += 2*x2;
	}


	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {
#ifdef DEBUG
		std::cout << "> SQP - addHessianEntriesTo - size of x (input) : " << x.rows() << " by " << x.cols() << std::endl;
		std::cout << "> SQP - addGradientTo - size of hessian (output) : " << hessianEntries.size() << std::endl;
#endif
		// Ex 1.2
		const double &x1 = x[0];
		const double &x2 = x[1];
		hessianEntries.push_back(Tripletd(0, 0, 2));
		hessianEntries.push_back(Tripletd(1, 1, 2));
	}

};
*/