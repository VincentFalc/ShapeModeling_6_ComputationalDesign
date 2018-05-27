---
author: VF
layout: post-full
type: image
featimg: woody-high.gif
title: Shape Modeling Shape Optimization
tags: [Shape Modeling, IGL, Eigen, Shape, Optimization]
category: [Shape Modeling]
---

##### TL;DR
This project was linked to a lecture about Shape Modeling. This project presents a method to optimize shapes by automatic deformation, according to arbitrary constraints.
<br/>
The followed constraint is to make the input mesh to stand in equilibrium on some chosen vertices, as if we wanted to 3D-print it, and make it stand by itself.
<br/>
We would like to add constraints, for example if the user want to impose the position of some vertices of the mesh in space.
<br/>

##### Table of Contents  
[1. Sequential Quadratic Programming](#Sequential)  <br/>
[1.1. Minimizer Implementeation](#Minimizer)  <br/>
[2. ”Make It Stand” - Editor](#Editor)  <br/>
[2.1. Equality Constraint](#Equality)  <br/>
[2.2. Inequality Constraint](#Inequality)  <br/>
[2.3. The most amazing balancing object!](#amazing)  <br/>

<a name="Sequential"/>

###" 1. Sequential Quadratic Programming

<a name="Minimizer"/>

#### 1.1. Minimizer Implementation

Again, the problem is an optimization problem. We are trying to fit to an objective function and to fit to some constraint.
<br/>
The whole problem can be solved by constructing a constrained system, and so by building a system-equivalent matrix.
As a first step, we can try our function on a simple problem : We **minimize the sum of the square of 2 values**, whereas we want the two values to validate the equation **x1 + 2*x2 = 3**
<br/>
We use SQPMinimizer, a child of the class ObjectiveFunction to build the related system.

###### Usage
We compute the value of expression to minimze, its gradient and its hessian anatically, and evaluate it on the provided coordinates.
<br/>
We consider the constraints system as following :
* B is 3 (left hand side of the constraints system).
* eqConstraintVals is the evaluation of the left hand side of the system, at given coordinates.
* jacobianEntries is the evaluation of the gradient of the left hand side of the system, at given coordinates.
* l and u are not constrained in this configuration, and so are set at low/high values respectively.

```C++

    (...)

	virtual double computeValue(const VectorXd& x){
		const double &x1 = x[0];
		const double &x2 = x[1];
		return std::pow(x1,2.0) + std::pow(x2, 2.0);
	}

	virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {
		const double &x1 = x[0];
		const double &x2 = x[1];
		grad[0] += 2*x1;
		grad[1] += 2*x2;
	}


	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {
		const double &x1 = x[0];
		const double &x2 = x[1];
		hessianEntries.push_back(Tripletd(0, 0, 2));
		hessianEntries.push_back(Tripletd(1, 1, 2));
	}

    (...)

	// Returns b of A(p) = b.
	virtual const VectorXd& getEqualityConstraintsTargetValues() {
		b.resize(1);
		b[0] = 3;
		return b;
	}

	// Returns A(p) of A(p) = b.
	// Derive from this to compute A(p). Fill `eqConstraintsVals` and return it.
	virtual const VectorXd& getEqualityConstraintValues(const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		eqConstraintVals.resize(1);
		eqConstraintVals[0] = x1 + 2 * x2;

		return eqConstraintVals;
	}

	// Computes the Jacobian dA/dp of the equality constraints A.
	virtual void addEqualityConstraintsJacobianEntriesTo(std::vector<Tripletd>& jacobianEntries, const VectorXd& p) {
		const double &x1 = p[0];
		const double &x2 = p[1];
		jacobianEntries.push_back(Tripletd(0, 0, 1)); // 2*x1 ? 
		jacobianEntries.push_back(Tripletd(0, 1, 2)); // 2*x2 ? 
	}

	// Returns l of constraint l <= p <= u
	virtual const VectorXd& getBoundConstraintsMinValues() {
		l.resize(2);
		l[0] = -100; // ??? 
		l[1] = -100; // ??? 
		return l;
	}

	// Returns u of constraint l <= p <= u
	virtual const VectorXd& getBoundConstraintsMaxValues() {
		u.resize(2);
		u[0] = 100; // ??? 
		u[1] = 100; // ??? 
		return u;
	}
```

###### Result

The result converge in only one iteration : We are in a Quadratic programming situation.
Therefore, grad_f(x) is linear to x ; grad_C(x) is constant ; C(x) is linear to x. So the global system has an order of 1, and can be "directly" solved.
<br/>
The solution found effectively validate the objective function (minimization of the square of the two function, noted f(x)) and validate the constraint (affine equation equal 3, noted c(x)).
<br/>

```
x              = 0.6 1.2
f(x)           = 1.8
c(x)           = 3
# iterations   = 1
```

__________________________________

<a name="Editor"/>

### 2. ”Make It Stand” - Editor

<a name="Equality"/>

#### 2.1. Equality Constraint

We want to maximize the stability of the mesh, which can be interpreted as placing the center of mass of the mesh exactly in the middle of the vertices acting as support.
<br/>
These vertices, the support, are fixed in space. The user can change the shape of the object by fixing the position of other vertices, without being part of the support.

###### Usage
We want to compute the different parts of the constraint system.
<br/>
We compute the position of the middle of the support points. This is where we expect the center of mass to be. This is calculated by a simple min/max search and a mean computation.
We get the total mass by simply looping on each vertice of the mesh and summing the local mass values.
<br/>
Then we construct b (left side of the system), by putting the values of the constrained positions in the mesh points order. We finnally add one constraint at the end of the vector : the coordinates of the middle of the support points previously calculated.
<br/>
We construct A to allow to calcul A(p) and the jacobian easily. A is constructed as a diagonal of 1, on the coordinates we want to constrain, and a final row, composed of the local mass of points divided by the total mass of all points.
<br/>

```C++

	//Get the expected position of the Center Of Mass
	getAvgCoordsSupport()
	{
		(...)
		//Get the min max of the support points
		for (std::set<int>::iterator it = this->simParent->supportNodes.begin(); it != this->simParent->supportNodes.end(); ++it)
		{
			int indice = *it;
			double xPos = this->simParent->x[2*indice];

			(...)

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

		double averageCoords = 0.5 * (maxX + minX);

		//Final calculus of the right hand side constraint
		return averageCoords;
	}

	virtual const double getTotalMass()
	{
		//Get the min max of the support points
		for (int i = 0; i < nbMeshPoints; i++)
		{
			totalMass += this->simParent->m[2*i]; // mass of the point
		}

		(...)
		return totalMass;
	}

	// Returns b of A(p) = b.
	virtual const VectorXd &getEqualityConstraintsTargetValues()
	{
		(...)

		//ADD "x = xi_traget" constraints
		for (std::map<int, Vector2d>::iterator it = this->simParent->constraintNodes.begin(); it != this->simParent->constraintNodes.end(); ++it)
		{
			int indice = it->first;

			coefficients.push_back(Tripletd(2 * indice, 0, it->second[0]));		//Store the x constrained value
			coefficients.push_back(Tripletd(2 * indice + 1, 0, it->second[1])); //Store the y constrained value
		}

		// ADD Right side of the constraint system (0.5 * (xMax+XMin))
		double averageCoords = getAvgCoordsSupport();

		//Storage of the last constraint
		coefficients.push_back(Tripletd(sizeB - 1, 0, averageCoords)); //Store the y constrained value in the last cell

		(...)
		//Transform triplet into SparseMatrix
		b.setFromTriplets(coefficients.begin(), coefficients.end());
		(...)

		return b;
	}

	virtual const std::vector<Tripletd> getATriplet(std::vector<Tripletd> &coefficients, const VectorXd &p)
	{
		(...)

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
			coefficients.push_back(Tripletd(nbMeshPoints * 2, 2 * i, this->simParent->m[2*i] / totalMass)); //(mi)/M attributed to x components
		}

		return coefficients;
	}

	virtual const Eigen::SparseMatrix<double> getA(const VectorXd &p)
	{
		(...)
		//Transform triplet into SparseMatrix
		getATriplet(coefficients, p); //Get the A List of triplet
		A.setFromTriplets(coefficients.begin(), coefficients.end());

		return A;
	}

	// Returns A(p) of A(p) = b.
	virtual const VectorXd &getEqualityConstraintValues(const VectorXd &p)
	{
		Eigen::SparseMatrix<double> A = getA(p);
		eqConstraintVals.resize(A.rows());
		eqConstraintVals = A * p; //Compute A(p)
		return eqConstraintVals;
	}

	// Computes the Jacobian dA/dp of the equality constraints A.
	virtual void addEqualityConstraintsJacobianEntriesTo(std::vector<Tripletd> &jacobianEntries, const VectorXd &p)
	{
		getATriplet(jacobianEntries, p);
	}

```

###### Result
<div style="text-align:center">
<p align="center">
  <img style="height : 300px;" src="./results/woody-low.png"> 
  <img style="height : 300px;" src="./results/woody-low-1iter.png"> <br/>
  Left : Woody-low original | Right : Woody-low after one iteration<br/>
  <img style="height : 300px;" src="./results/woody-low-withConstraints.png">
  <img style="height : 300px;" src="./results/woody-high.gif"><br/>
  Left : Woody-low after iterations with constraints | Right : Woody-high evolution without and then with constraints (gif)<br/>
  <img style="height : 300px;" src="./results/square-1iter.png">
  <img style="height : 300px;" src="./results/square-withConstraints.png"><br/>
  Left : Square after one iteration | Right : Square with some constraints <br/>
  <img style="height : 300px;" src="./results/square-withConstraints-lot.png">
  <img style="height : 300px;" src="./results/square.gif"><br/>
  Left : Square with a lot of constraints | Right : Square evolution with and then with constraints (gif)<br/>
</p>
</div>

We see that the center of mass position (xCOM) is aligned with the center of the support (AVG).

```
Square :
AVG : -0.75
Max : -0.5
Min : -1
xCOM : -0.75

Woody-low : 
AVG : 112.61
Max : 135.514
Min : 89.7069
xCOM : 112.61

```
__________________________________

<a name="Inequality"/>

#### 2.2. Inequality Constraint

Forcing the center of mass to be exactly aligned with the center of the support may be too strong. We would like the center of mass to be, at least "inside" the support, which is enough to guarantee the stand of the object.
<br/>
We can translate it with inequality constraints.

###### Usage

If the inequality is selected, the matrix A is  computed in a slightly different manner : without the last line, corresponding to the equality constraint.
<br/>
Then, we add the inequality constraint as 2 "equality constraint" with d and f. These vector will be handled in a way (out of the code presneted below) that will correspond to an inequality constraint.
<br/>

```C++

	// Returns d of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMinValues() {
		(...)

		// ADD Right side of the constraint system (XMin)
		double minCoords = getMinSupport();

		//Storage of the last constraint
		coefficients.push_back(Tripletd(sizeD - 1, 0, minCoords)); //Store the y constrained value in the last cell

		Eigen::SparseMatrix<double> Dtmp(sizeD, 1);
		d.resize(sizeD);
		//Transform triplet into SparseMatrix
		Dtmp.setFromTriplets(coefficients.begin(), coefficients.end());
		d = VectorXd(Dtmp);

		return d;
	}

	// Returns f of d <= C(p) <= f
	virtual const VectorXd& getInequalityConstraintsMaxValues() {
		(...)

		// ADD Right side of the constraint system (XMin)
		double maxCoords = getMaxSupport();

		//Storage of the last constraint
		coefficients.push_back(Tripletd(sizeF - 1, 0, maxCoords)); //Store the y constrained value in the last cell

		Eigen::SparseMatrix<double> Ftmp(sizeF, 1);
		f.resize(sizeF);
		//Transform triplet into SparseMatrix
		Ftmp.setFromTriplets(coefficients.begin(), coefficients.end());
		f = VectorXd(Ftmp);

		return f;
	}

	// Computes the Jacobian dA/dp of the inequality constraints C.
	virtual void addInequalityConstraintsJacobianEntriesTo(std::vector<Tripletd>& jacobianEntries, const VectorXd& p) {
		(...)
		
		//ADD XCom(x)
		for (int i = 0; i < nbMeshPoints; i++)
		{
			jacobianEntries.push_back(Tripletd(0, 2 * i, this->simParent->m[2*i] / totalMass)); //(mi)/M attributed to x components
		}
	}

```

###### Result
<div style="text-align:center">
<p align="center">
  <img style="height : 300px;" src="./results/square-inequal-iter1.png">
  <img style="height : 300px;" src="./results/square-inequal-stable.png"><br/>
  Left : Square after 1 iteration | Right : Square in stable state<br/>
  <img style="height : 300px;" src="./results/woody-low-inequality-iter1.png">
  <img style="height : 300px;" src="./results/woody-low-inequality-constraints.png"><br/>
  Left : Woody after 1 iteration (stable) | Right : Woody stable with constraints<br/>
</p>
</div>

We see that the center of mass position (xCOM) is "inside" the support (between Max and Min).

```
Square : 
AVG : -0.8
Max : -0.6
Min : -1
xCOM : -0.6

Woody-low :
AVG : 112.61
Max : 135.514
Min : 89.7069
xCOM : 135.514

```

__________________________________

<a name="amazing"/>

#### 2.3. The most amazing balancing object! 

The whole code is effective for 2D Meshes. But we can use it with 3D meshes, if the mesh is not "too thick". 
<br/>
For example, an octopod would not works, as it is heavily dispersed in the 3D space. Whereas a mesh of a character, approximatively flat, could be deformed successfully.
<br/>

###### Result
<div style="text-align:center">
<p align="center">
  <img style="height : 600px;" src="./results/Danseuse2.png">
  <img style="height : 600px;" src="./results/printed.jpg">
</p>
</div>

The second picture is a 3D print it, but it was too small to really stay on his foot as expected. However, we can feel that it does has a equilbrium point.

__________________________________


For more information, see : 
[-> Repository](https://github.com/VincentFalc/ShapeModeling_6_ComputationalDesign)

Note that the code is not optimized by thoughtful laziness: this code is not made to be reused as is in other applications, and if this is the case, the code will be reviewed and optimized. Thanks for your understanding.
