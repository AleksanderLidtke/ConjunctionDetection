#include "VectorOperations.h"

namespace VectorOperations{
	void linspace(std::vector<double>* outputVecPtr, double start, double stop, int* noPoints){
		/* Fill the vector with noPoints evenly-spaced numbers between start and stop (inclusive).
		@param outputVectPtr - vector which will be filled with numbers.
		@param start, stop - beginning and end numbers of the sequence. 
		@param noPoints - number of points that will be produced, has to be equal to length of outputVecPtr.
		*/
		double spacing = (stop-start)/( (double) *noPoints - 1.0 ); // Spacing between the points.
		for(std::vector<int>::size_type  i=0; i<outputVecPtr->size(); i++){
			outputVecPtr->at(i) = start + ((int) i )*spacing;
		};
	};

	void centralLinspace(std::vector<double>* outputVecPtr, double centre, double spacing, int* noPoints){
		/* Given a number at the centre of an interval create a set of evenly spaced numbers centred on the given value. 
		Fill a vector with those values.
		@param outputVectPtr - vector which will be filled with numbers.
		@param centre - gives the centre of the interval.
		@param spacing - the separation between the numbers in the interval.
		@param noPoints - number of points that will be produced, has to be equal to length of outputVecPtr and be even.
		@throws FatalError if the number of points desired is not even and hence cannot be centred on a location with an even number of points on each side.
		*/
		if( *noPoints%2 != 0 ){
			std::string msg1 = "Requested an odd number of points to be put on both sides of a value, which is impossible.";
			throw FatalError(msg1,__FILE__,__LINE__);
		}else{
			linspace( outputVecPtr, centre-(*noPoints-1)/2.0*spacing, centre+(*noPoints-1)/2.0*spacing, noPoints);
		}
	};

	void vectorDifference(std::vector<double>* outputVecPtr, std::vector<double>* vect1Ptr, std::vector<double>* vec2Ptr){
		/* Compute the difference of the two vectors (vector2 - vector1) of any length.
		@param outputVecPtr - pointer to a vector where the difference will be stored.
		@param vec1Ptr, vec2Ptr - pointers to both vectors difference between which is to be computed. Must be the same size as outputVecPtr.
		*/
		for(std::vector<int>::size_type  i=0; i<vec2Ptr->size(); i++){
			outputVecPtr->at(i) = vec2Ptr->at(i)-vect1Ptr->at(i);
		};
	};

	double vectorMagnitudeSquared(std::vector<double>* vecPtr){
		/* Find the magnitude of a vector squared (avoid taking the square root for speed).
		@param vecPtr - pointer to the vector whose magnitude will be computed, can be any length.
		@return - sum of the squares of all the vecPtr's entries.
		*/
		double result = 0.;
		for(std::vector<int>::size_type  i=0; i<vecPtr->size(); i++){
			result = result + vecPtr->at(i)*vecPtr->at(i);
		};
		return result;
	};

	double vectorMagnitude(std::vector<double>* vecPtr){
		/* Find the magnitude of a vector.
		@param vecPtr - pointer to the vector whose magnitude will be computed, can be any length.
		@return - sum of the squares of all the vecPtr's entries.
		*/
		double result = 0.;
		for(std::vector<int>::size_type  i=0; i<vecPtr->size(); i++){
			result = result + vecPtr->at(i)*vecPtr->at(i);
		};
		return std::sqrt( result );
	};

	void unitVector(std::vector<double>* vecPtr){
		/* Scale the vector to a unit vector (divide by its two-norm).
		@param vecPtr - vector to be made unit.
		*/
		double twoNorm = sqrt( vectorMagnitudeSquared(vecPtr) ); // Compute the two-norm.

		for(std::vector<int>::size_type  i=0; i<vecPtr->size(); i++){
			vecPtr->at(i) = vecPtr->at(i)/twoNorm; // Divide each entry.
		};
	};

	double dotProduct(std::vector<double>* vec1Ptr, std::vector<double>* vec2Ptr){
		/* Compute a dot (inner) product of two vectors.
		@param vec1Ptr, vec2Ptr - pointers to the vectors to be dorred together, must be the same length.
		@return square root of the sum of each dimension's components multiplied together, e.g. sqrt(a1*a2 + b1*b2 + c1*c2).
		*/
		double result =0.;
		for(std::vector<int>::size_type  i=0; i<vec1Ptr->size(); i++){
			result = result + vec1Ptr->at(i)*vec2Ptr->at(i);
		};
		return result;
	};

	void crossProduct(std::vector<double>* outputVecPtr, std::vector<double>* vec1Ptr, std::vector<double>* vec2Ptr){
		/* Compute the cross (outer) product of two vectors of lenght 3.
		@param outputVecPtr - the result will be stored in this vector, size 3.
		@param vec1Ptr, vec2Ptr - vectors to be crossed together, size 3.
		*/

		outputVecPtr->at(0) = vec1Ptr->at(1)*vec2Ptr->at(2) - vec1Ptr->at(2)*vec2Ptr->at(1);
		outputVecPtr->at(1) = vec1Ptr->at(2)*vec2Ptr->at(0) - vec1Ptr->at(0)*vec2Ptr->at(2);
		outputVecPtr->at(2) = vec1Ptr->at(0)*vec2Ptr->at(1) - vec1Ptr->at(1)*vec2Ptr->at(0);
	};

	std::vector<double> vectorMultiplyByScalar(std::vector<double>* vecPtr, double coefficient){
		/* Multiply a vector by a scalar.
		@param vecPtr - vector to be multiplied, any length.
		@param coefficient - the scalar by which the entries of the vector will be multiplied by.
		@return - vector with the same shape as vecPtr but with each of the entries multiplied by the coefficient.
		*/
		std::vector<double> result = std::vector<double>(vecPtr->size(), 0.0);
		for(std::vector<int>::size_type  i=0; i<vecPtr->size(); i++){
			result.at(i) = coefficient * vecPtr->at(i);
		};

		return result;
	};

	std::vector< std::vector<double> > matrixMultiplyByScalar(std::vector< std::vector<double> >* matPtr, double coefficient){
		/* Multiply a matrix by a scalar.
		@param matPtr - 2D matrix to be multiplied, any length of both axes.
		@param coefficient - the scalar by which the entries of the matrix will be multiplied by.
		@return - matrix with the same shape as matPtr but with each of the entries multiplied by the coefficient.
		*/
		std::vector< std::vector<double> > result = std::vector< std::vector<double> >( matPtr->size(), std::vector<double>(matPtr->at(0).size(), 0.0) );
		for(std::vector<int>::size_type  i=0; i<matPtr->size(); i++){ // Outer loop - go through all the rows.
			for(std::vector<int>::size_type  j=0; j<matPtr->at(0).size(); j++){ // Inner loop - go through all the columns in every row.
				result.at(i).at(j) = coefficient * matPtr->at(i).at(j);
			};
		};
		return result;
	};

	std::vector<double> vectorMultiplyByMatrix(std::vector<double>* vecPtr, std::vector< std::vector<double> >* matPtr){
		/* Multiply a square matrix by a vector. It is assumed that the shapes of both will be suitable for such multiplication.
		@param vacPtr - vector v of size N.
		@param matPtr - 2D matrix A of size NxN.
		@result - vector of size N that is the result of Av.
		@throws - std::out_of_range if the shapes of the matrix and the vector are not suitable to perform the multiplication.
		*/
		std::vector<double> result = std::vector<double>(vecPtr->size(), 0.0);
		for(std::vector<int>::size_type  i=0; i<matPtr->size(); i++){ // Outer loop - go through all the rows of the matrix.
			for(std::vector<int>::size_type  j=0; j<matPtr->at(0).size(); j++){ // Inner loop - go through all the columns of the current row.
				result.at(i) = result.at(i) + matPtr->at(i).at(j)*vecPtr->at(j);
			};
		};
		return result;
	};

	std::vector< std::vector<double> > matrixMultiplyByMatrix(std::vector< std::vector<double> >* mat1Ptr, std::vector< std::vector<double> >* mat2Ptr){
		/* Multiply two matrices together in the following order: mat1Ptr*mat2Ptr (premultiply matrix2 by matrix1).
		@param mat1Ptr - NxM 2D matrix.
		@param mat2Ptr - MxP 2D matrix.
		@return - NxP 2D matrix.
		@throws - std::out_of_range if the shapes of the matrices are not suitable to perform the multiplication.
		*/
		double sum; // Temporary variable.
		//int N = (int) mat1Ptr->size();
		//int M = (int) mat1Ptr->at(0).size();
		//int P = (int) (int) mat2Ptr->at(0).size();

		std::vector< std::vector<double> > result = std::vector< std::vector<double> >( mat1Ptr->size(), std::vector<double>(mat2Ptr->at(0).size(), 0.0) );
		for(std::vector<int>::size_type  i=0; i<result.size(); ++i){ // Go through result rows - N.
	 		for (std::vector<int>::size_type j=0; j<result.at(0).size(); ++j){ // And collumns - P.
				sum = 0; // Reset this every time a new element of the result is to be computed.
			
				for (std::vector<int>::size_type k=0; k<mat2Ptr->size(); ++k){ // Go through M and sum element-wise multiplication of the ith row of mat1 with the jth column of mat2.
					sum += mat1Ptr->at(i).at(k) * mat2Ptr->at(k).at(j);
				}
			
				result.at(i).at(j) = sum;
			};
		};
		/* OLD VERSION */
		/*double sum; // Temporary variable.
		std::vector< std::vector<double> > result = std::vector< std::vector<double> >( mat1Ptr->size(), std::vector<double>(mat2Ptr->at(0).size(), 0.0) );
		for(std::vector<int>::size_type  i=0; i<mat1Ptr->size(); ++i){
	 		for (std::vector<int>::size_type j=0; j<mat2Ptr->size(); ++j){
				sum = 0; // Reset this.
				for (std::vector<int>::size_type k=0; k<mat2Ptr->at(0).size(); ++k){
					sum += mat1Ptr->at(i).at(k) * mat2Ptr->at(k).at(j);
				}
				result.at(i).at(j) = sum;
			};
		};*/
		return result;
	};

	double twoByTwoDeterminant( std::vector< std::vector<double> >* matPtr ){
		/* Compute the determinant of a 2x2 matrix. Use the fact that for this type of matrix an easy analytical solution is well known.
		@param - 2x2 matrix.
		*/
		return matPtr->at(0).at(0)*matPtr->at(1).at(1) - matPtr->at(0).at(1)*matPtr->at(1).at(0);
	};

	double threeByThreeDeterminant( std::vector< std::vector<double> >* matPtr ){
		/* Compute the determinant of a 3x3 matrix. Use the fact that for this type of matrix an easy analytical solution is well known.
		@param - 3x3 matrix.
		*/
		return matPtr->at(0).at(0)*( matPtr->at(1).at(1)*matPtr->at(2).at(2) - matPtr->at(1).at(2)*matPtr->at(2).at(1) ) - matPtr->at(0).at(1)*( matPtr->at(1).at(0)*matPtr->at(2).at(2) - matPtr->at(1).at(2)*matPtr->at(2).at(0) ) + matPtr->at(0).at(2)*( matPtr->at(1).at(0)*matPtr->at(2).at(1) - matPtr->at(1).at(1)*matPtr->at(2).at(0) );
	};

	std::vector< std::vector<double> > inverseTwoByTwo( std::vector< std::vector<double> >* matPtr ){
		/* Inverse a 2x2 matrix. Use the fact that for this type of matrix an easy analytical solution is well known.
		@param - a 2x2 matrix.
		@throws - MathWarning when the matrix cannot be inversed, i.e. its determinant is 0.
		*/
		std::vector< std::vector<double> > result = std::vector< std::vector<double> >( 2, std::vector<double>(2,0.0) );
		result.at(0).at(0) = matPtr->at(1).at(1); // 1st row of the inverted matrix.
		result.at(0).at(1) = -1.0*matPtr->at(0).at(1);
		result.at(1).at(1) = matPtr->at(0).at(0); // Second row.
		result.at(1).at(0) = -1.0*matPtr->at(1).at(0);

		double det = twoByTwoDeterminant( matPtr );
		if( (det<DBL_MIN) && (det>-1.0*DBL_MIN) ){ // Deteriminant is very, very small.
			std::string msg = "Attempted to inverse a matrix with determinant equal to 0.";
			throw MathWarning(msg,__FILE__,__LINE__);
		}
		return matrixMultiplyByScalar(&result, 1.0/det); // Multiply through by the reciprocal of the determinant et viola!
	};

	std::vector< std::vector<double> > inverseThreeByThree( std::vector< std::vector<double> >* matPtr ){
		/* Inverse a 3x3 matrix. Use the fact that for this type of matrix an easy analytical solution is well known.
		@param - a 3x3 matrix.
		@throws - MathWarning when the matrix cannot be inversed, i.e. its determinant is 0.
		*/
		std::vector< std::vector<double> > result = std::vector< std::vector<double> >( 3, std::vector<double>(3,0.0) );
		result.at(0).at(0) = matPtr->at(1).at(1)*matPtr->at(2).at(2) - matPtr->at(1).at(2)*matPtr->at(2).at(1); // 1st row of the inverted matrix.
		result.at(0).at(1) = matPtr->at(0).at(2)*matPtr->at(2).at(1) - matPtr->at(0).at(1)*matPtr->at(2).at(2);
		result.at(0).at(2) = matPtr->at(0).at(1)*matPtr->at(1).at(2) - matPtr->at(0).at(2)*matPtr->at(1).at(1);
		result.at(1).at(0) = matPtr->at(1).at(2)*matPtr->at(2).at(0) - matPtr->at(1).at(0)*matPtr->at(2).at(2); // Second row.
		result.at(1).at(1) = matPtr->at(0).at(0)*matPtr->at(2).at(2) - matPtr->at(0).at(2)*matPtr->at(2).at(0);
		result.at(1).at(2) = matPtr->at(0).at(2)*matPtr->at(1).at(0) - matPtr->at(0).at(0)*matPtr->at(1).at(2);
		result.at(2).at(0) = matPtr->at(1).at(0)*matPtr->at(2).at(1) - matPtr->at(1).at(1)*matPtr->at(2).at(0); // Third row.
		result.at(2).at(1) = matPtr->at(0).at(1)*matPtr->at(2).at(0) - matPtr->at(0).at(0)*matPtr->at(2).at(1);
		result.at(2).at(2) = matPtr->at(0).at(0)*matPtr->at(1).at(1) - matPtr->at(0).at(1)*matPtr->at(1).at(0);
		double det = threeByThreeDeterminant( matPtr );
		if( (det<DBL_MIN) && (det>-1.0*DBL_MIN) ){ // Deteriminant is very, very small.
			std::string msg = "Attempted to inverse a matrix with determinant equal to 0.";
			throw MathWarning(msg,__FILE__,__LINE__);
		}
		return matrixMultiplyByScalar(&result, 1.0/det); // Multiply through by the reciprocal of the determinant et viola!
	};

	std::vector< std::vector<double> > transposeMatrix( std::vector< std::vector<double> >* originalPtr ){
		/* Find the transpose of a square matrix. */
		std::vector< std::vector<double> >  output = std::vector< std::vector<double> >( originalPtr->at(0).size(), std::vector<double>(originalPtr->size(), 0.0) ); // Initialise the output with swapped dimensions c.f. original matrix.
	
		for(size_t i=0; i<originalPtr->size(); i++){ // Go through rows of the original matrix.
			for(size_t j=0; j<originalPtr->at(i).size(); j++){ // They become collumns of the transposed matrix.
				output.at(j).at(i) = originalPtr->at(i).at(j);
			}
		}

		return output;
	};

	std::vector< std::vector<double> > addMatrices(std::vector< std::vector<double> >* mat1Ptr, std::vector< std::vector<double> >* mat2Ptr){
		/* Add two matrices together, i.e. element-wise addition, and return the result.
		@param pointers to two 2D matrices of arbitrary but identical size.
		@return matrix elements of which are sums of the elements of mat1 and mat2, same size as mat1 and mat2.
		*/
		std::vector< std::vector<double> > result = std::vector< std::vector<double> >( mat1Ptr->size(), std::vector<double>(mat1Ptr->at(0).size(), 0.0) );
		for(std::vector<int>::size_type i=0; i<mat1Ptr->size(); i++){ // Go through the rows.
			for(std::vector<int>::size_type j=0; j<mat1Ptr->size(); j++){ // Go through the columns.
				result.at(i).at(j) = mat1Ptr->at(i).at(j) + mat2Ptr->at(i).at(j);
			}
		}
		return result;
	};

	double mean(std::vector<double>* inputVecPtr){
		/* Compute the arithmetic mean of the entries in a vector.
		@param - pointer to a vector for which the arithmetic mean of the entries is to be computed.
		*/
		double result = 0.0;
		for(std::vector<int>::size_type i=0; i<inputVecPtr->size(); i++){ // Compute mean of the residuals.
			result += inputVecPtr->at(i);
		}
		result = result/(double) inputVecPtr->size();
		return result;
	};

	double standardDeviation(std::vector<double>* inputVecPtr){
		/* Compute the standard deviation of the entries in a vector.
		@param - pointer to a vector for which the standard deviation is to be computed.
		*/
		double result = 0.0;
		double average = mean( inputVecPtr ); // Need the mean of the vector to compute the standard deviation.

		for(std::vector<int>::size_type i=0; i<inputVecPtr->size(); i++){ // Compute the standard deviation.
			result += (inputVecPtr->at(i)-average)*(inputVecPtr->at(i)-average);
		}
	
		result = result/(double) inputVecPtr->size(); // Final operations on the final result, i.e. on the sum of squares of differences between the terms and the mean.
		result = std::sqrt(result);
		return result;
	};
}//END namespace VectorOperations

namespace EquationsSolving{
	std::vector< std::vector<double> > luDecomposition(std::vector< std::vector<double> > A, std::vector<int> &index, double &d){
		/* Decompose the coefficients' matrix A into upper- and lower-triangular parts (LU decomposition)
		and doing row swappings, if necessary. Use the Crout's method to do the LU decomposition.
		@param A - matrix of coefficients, NxN.
		@param index - stores information about the row swaps that have been performed.
		@param d - keeps count of how many row swaps have been made (will be +ve if the number was even and -ve if it was odd).
		@returns - NxN matrix that is a combination of the matrices L and U.
		@raises FatalError - raise this exception if the provided matrix A is singular.
		@reference Numerical Recipes in C++, Second Edition
		*/
		const double smallNum = 1.0e-20; // Used to help properly condition the matrix (used instead of 0.0 to avoid singularities).
		double bigNum = 0.0, sum = 0.0; // Intermediate variables.
		double dum = 0.0;	// A dummy storage variable.
		int imax = 0;	// Location of the maximum element for pivoting.

		int n = int(A.size()); // No. rows of the matrix.
		std::vector<int> indx(n); // Values of this vector are transferred to index input -> avoid out of range memory problems.
		std::vector<double> vv(n); // Stores the implicit scaling factors of each row.
		d = 1.0; // No rows interchanged so far -> starting value is +ve.

		for(int i = 0; i < n; i++){	// Do the scaling of each row.
			double temp;
			// Loop over the matrix to find the scaling value to use implicit pivoting, i.e. scale all biggest row values to range [0,1].
			for(int j = 0; j < n; j++)
				if((temp = fabs(A[i][j])) > bigNum) bigNum = temp;
			if(bigNum == 0.0) {
				char msg[] = "Singular matrix encountered when performing the LU decomposition: no non-zero largest element in a given row.\n";
				throw FatalError(msg,__FILE__,__LINE__);
			}
			vv[i] = 1.0/bigNum;
		}

		for(int j = 0; j < n; j++){	// Loop over the columns to use the Crout's method.
			for (int i = 0; i < j; i++){
				// beta_ij = A_ij - sum alpha_ik beta_kj   for i!=j (off-diagonal)
				sum = A[i][j];
				for (int k = 0; k < i; k++) sum -= A[i][k]*A[k][j];
				A[i][j] = sum;
			}
			bigNum = 0.0;
			for (int i = j; i < n; i++){ // Same equation as above but for i=j (diagonal).
				sum = A[i][j];
				for (int k = 0; k < j; k++) sum -= A[i][k]*A[k][j];
				A[i][j] = sum;
				if ((dum = vv[i]*fabs(sum)) >= bigNum){	// Check if this is a better pivot than the previous ones.
					bigNum = dum;
					imax = i;
				}
			}
			if (j != imax){	// Check if rows need to be swapped.
				for (int k = 0; k < n; k++){ // Do the swap by iterating over each element of the imax'th row.
					dum = A[imax][k];
					A[imax][k] = A[j][k];
					A[j][k] = dum;
				}
				d = -d; // Change the pairity of d -> indicates if rows are being swapped.
				vv[imax] = vv[j]; // Swap the scale factor as well.
			}
			indx[j] = imax; // Store the row permutation index.
			// Reduce singularity of the matrix that could be cuased by machine precision by substituting a small number for the zero.
			if (A[j][j] == 0.0) A[j][j] = smallNum;
			if (j != n-1){
				dum = 1.0/A[j][j];	// Divide by the pivot element.
				for (int i = j+1; i < n; i++) A[i][j] *= dum;
			}
		}
		index = indx;
		return A;	// The copy of the initially supplied matrix has been LU decomposed now.
	};

	void luSubst(std::vector < std::vector<double> >* A, std::vector<int>* index, std::vector<double>* b){
		/* Solve a system of equations Ax=b for the unknown vector x by using a LU decomposed coefficients' matrix A and the
		information about the row swaps that were performed when decomposing it.
		@param A - decomposed matrix of coefficients, NxN, as returned by @see luDecomposition.
		@param index - stores information about the row swaps that have been performed as returned by @see luDecomposition.
		@param b - right hand side vector, length N. The solution of the set of equations, i.e. the vector "x" will be written to this RHS vector.
		@reference Numerical Recipes in C++, Second Edition
		*/

		int N = int(A->size());	// No. rows in the matrix i.e. no. equations and variables (assume a well-conditioned set of equations).
		int ii = 0;
		for (int i = 0; i < N; i++){ // Iterate over rows.
			int ip = index->at(i); // Unscramble the row permutations -> shift rows accorind to the index.
			double sum = b->at(ip); // Use sum as a temporary variable.
			b->at(ip) = b->at(i);
			if (ii != 0) // +ve ii becomes the index of the first non-vanishing element of b.
				for (int j = ii-1; j < i; j++) sum -= A->at(i).at(j)*b->at(j); // Do the forward substitution.
			else if (sum!= 0.0)
				ii = i+1;
			b->at(i) = sum;
			//b[i] = sum;
		}
		for (int i = N-1; i >= 0; i--){	// Do the back sbstitution.
			double sum = b->at(i);
			for (int j = i+1; j < N; j++) sum -= A->at(i).at(j)*b->at(j);
			b->at(i) = sum/A->at(i).at(i);
		}
	};
}//END namespace EquationsSolving
