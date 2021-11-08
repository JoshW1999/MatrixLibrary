import java.util.Random;
import java.lang.Math;
import java.util.ArrayList; 
import java.util.Stack;
import java.util.Arrays; 


public class Matrix{ 

	
	int nRows;
	int nCols; 

	double[][] mat;
	
	private static final double EPSILON = 0.00001;

	Matrix(int rows, int cols){
		this.nRows = rows; 
		this.nCols = cols; 
		this.mat = new double[rows][cols];
	}

	Matrix(double[][] m){

		this.nRows = m.length; 
		this.nCols = m[0].length; 
		this.mat = m; 
	}	

	Matrix(double[] m){
		this.nRows = 1; 
		this.nCols = m.length; 

		this.mat = new double[1][m.length];
		this.mat[0] = m;

	}

	// Copy Constructor
	Matrix(Matrix m){
		this.mat = deepCopy(m.mat);
		this.nRows = m.nRows; 
		this.nCols = m.nCols;  
	}

	public double[][] deepCopy(double[][] matrix) {
   		return java.util.Arrays.stream(matrix).map(el -> el.clone()).toArray($ -> matrix.clone());
	}

	public double[][] toArray(){
		return this.mat; 
	}

	public int[] shape(){
		int[] s = {this.nRows, this.nCols};
		return s;
	}

	public static boolean equals(Matrix m1, Matrix m2){
		return Arrays.deepEquals(m1.mat, m2.mat);
	}

	public void printMatrix(){

		
		for(int i = 0; i < nRows; i++){
			System.out.print("[ ");
			for(int j = 0; j < this.nCols; j++){
				System.out.print(this.mat[i][j] + " ");
			}
			System.out.println(" ]");
		}
		System.out.println();
	}

	// Fills Matrix with Ones
	public static Matrix ones(Matrix m){

		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				m.mat[i][j] = (double)1;

		return m; 
	}

	public static Matrix ones(int rows, int cols){

		Matrix m = new Matrix(rows, cols);

		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				m.mat[i][j] = (double)1;

		return m; 
	}

	// Matrix.zeros?

	// Returns Identity Matrix of n-by-n dimensionality
	// Note: Identity matricies are necessarily square matrices
	public static Matrix eye(int n){
		Matrix m = new Matrix(n,n);

		for(int i = 0; i < m.nRows; i++)
			m.mat[i][i] = (double)1;
		
		return m; 
	}

	// Returns Matrix with randomly initialized elements between 0 and 1
	public static Matrix rand(int rows, int cols){
		Matrix m = new Matrix(rows, cols);

		for(int i = 0; i < rows; i++)
			for(int j = 0; j < cols; j++)
				m.mat[i][j] = Math.random();
		
		return m; 
	}

	// Returns Matrix with randomly initialized integers between min and max
	public static Matrix rand(int rows, int cols, int min, int max){
		Matrix m = new Matrix(rows, cols);

		for(int i = 0; i < rows; i++)
			for(int j = 0; j < cols; j++)
				m.mat[i][j] = (int)(Math.random()*(max-min)) + min;

		return m; 
	}

	// Returns diagonal of Matrix
	// Note: m must be a square matrix
	public static Matrix diag(Matrix m) throws Exception{
		
		if(m.nRows != m.nCols)
			throw new NonSquareMatrixException("Matrix must have equal dimensions: (" + m.nRows + "," + m.nCols + ").");


		Matrix diag_m = new Matrix(1, m.nRows);

		for(int i = 0; i < m.nRows; i++)
			diag_m.mat[0][i] = m.mat[i][i]; 

		return diag_m; 
	}
	
	// Methods for getting column from Matrix
	public static Matrix getColumn(Matrix m, int index){
		return new Matrix(Matrix.getColumn(m.mat, index));
	}

	public static double[] getColumn(double[][] array, int index){
   	 	double[] column = new double[array[0].length]; // Here I assume a rectangular 2D array! 
    	
    	for(int i=0; i<column.length; i++)
       		column[i] = array[i][index];
    
   	 	return column;
	}


	// Returns Column Vector of specified column index
	public static Matrix getColVec(Matrix m, int index){
		Matrix ret = new Matrix(m.mat[0].length,1);

		for(int i = 0; i < ret.nRows; i++)
			ret.mat[i][0] = m.mat[i][index];

		return ret; 
	}


	// Returns Reduced Row Echelon Form of Matrix in format A|b
	public static Matrix REF(Matrix m) throws Exception{
		
		if(m.nRows != m.nCols-1)
			throw new NonSquareMatrixException("Submatrix A is not square: (" + m.nRows + "," + (m.nCols-1) + ").");

		int N = m.nRows;

		for(int i = 0; i < N-1; i++){ // -1 since last column is b

			int max = i; 

			for(int j = i+1; j < N; j++){

				if(Math.abs(m.mat[j][i]) > Math.abs(m.mat[max][i]))
					max = j; 
			}

			// Exchanging Rows such that max is at position i
			double[] temp = m.mat[i];
			m.mat[i] = m.mat[max];
			m.mat[max] = temp;



			if(Math.abs(m.mat[i][i]) <= EPSILON)
				throw new SingularMatrixException("Matrix is singular.");

			for(int k = i+1; k < N; k++){
				double alpha = m.mat[k][i] / m.mat[i][i];

				m.mat[k][N] -= alpha*m.mat[i][N];

				for(int j = i; j < N; j++){
					m.mat[k][j] -= (alpha * m.mat[i][j]);
				}
			}
		}

		return m; 
	}


	// Returns x in equation Ax = b
	// E.g. solution to system of equations.
	public static Matrix RREF(Matrix m) throws Exception{
		
		int N = m.nRows;

		m = Matrix.REF(m);
		
		// RREF Operation
		Matrix x = new Matrix(1,N);

		for(int i = N-1; i >= 0; i--){
			double sum = 0.0;

			for(int j = i+1; j < N; j++){
				sum += m.mat[i][j]*x.mat[0][j];
			}

			x.mat[0][i] = (m.mat[i][N] - sum) / m.mat[i][i]; 
		}

		//x.printMatrix();
		return x; 
	}
	
	// Returns linear product of all matrix arguments
	public static Matrix mul(Matrix ...args) throws Exception{
		Matrix ret = args[0];

		for(int m = 1; m < args.length; m++){
			
			if(ret.nCols != args[m].nRows)
				throw new IncompatibleDimensionsException("Incompatible Dimensions: (" + ret.nRows + ", " + ret.nCols + ") and (" + args[m].nRows + ", " + args[m].nCols + ")" );

			Matrix res = new Matrix(ret.nRows, args[m].nCols);

			for(int i = 0; i < ret.nRows; i++){

				for(int j = 0; j < args[m].nCols; j++){

					for(int k = 0; k < args[m].nRows; k++){

						res.mat[i][j] += ret.mat[i][k] * args[m].mat[k][j];

					}

				}

			}

			ret = res; 

		}

		return ret; 
	}


	// PLU Decommposition of Matrix m
	// where m = PLU
	// P: Permutation Matrix (Tracks row swaps)
	// L: Lower Triangular Matrix
	// U: Upper Triangular Matrix
	public static Matrix[] PLU(Matrix m) throws Exception{
		
		if(m.nRows != m.nCols)
			throw new NonSquareMatrixException("The matrix m (" + m.nRows + "," + m.nCols + ") must be square!");

		int N = m.nRows;

		double subAmt; 
		double focus; 

		Matrix lower = Matrix.eye(N);
		Matrix upper = new Matrix(m); 
		Matrix p = Matrix.eye(N);

		Stack<int[]>swapList = new Stack<int[]>();
				
		
		// Returns permutation matrix at index 0, 
		// lower matrix at index 1,
		// upper matrix at index 2.  
		Matrix[] ret = new Matrix[3];

		// Essentially just performing REF operations here
		// But we are logging alpha in the lower triangular matrix
		// as well as swapping lower triangular matrix rows when necessary.
		for(int i = 0; i < N; i++){

			int max = i; 

			for(int j = i+1; j < N; j++){

				if(Math.abs(upper.mat[j][i]) > Math.abs(upper.mat[max][i])){
					max = j; 	
				}

			}

			// Exchanging Rows such that max is at position i in lower triangular matrix
			if(max != i && i > 0){
				for(int k = 0; k < i; k++){
					
					double temp2 = lower.mat[i][k];
					lower.mat[i][k] = lower.mat[max][k];
					lower.mat[max][k] = temp2; 
					
				}	
			}
			

			if(i != max){
				int[] swapPair = {i, max};
				swapList.push(swapPair);
			}

			double[] temp = upper.mat[i];
			upper.mat[i] = upper.mat[max];
			upper.mat[max] = temp;
			
			if(Math.abs(upper.mat[i][i]) <= EPSILON)
				throw new SingularMatrixException("Matrix is singular");

			for(int k = i+1; k < N; k++){
				
				double alpha = upper.mat[k][i] / upper.mat[i][i];
				
				// Logging alpha into Lower Matrix
				lower.mat[k][i] = alpha;  
				

				for(int j = i; j < N; j++){
					upper.mat[k][j] -= (alpha * upper.mat[i][j]);
				}

				
			}
		}

		ret[1] = lower; 
		ret[2] = upper;
		
		// Applying stored swapping operations
		// to permutation matrix such that it encodes
		// all swaps applied to upper and lower. 

		while(!swapList.empty()){
			int[] sle = swapList.pop();

			double[] t = p.mat[sle[0]];
			p.mat[sle[0]] = p.mat[sle[1]];
			p.mat[sle[1]] = t;
		}


		ret[0] = p;

		return ret; 
	}


	// PLU, but with return values altered
	// to return det(p) as well, which is 1 if an even
	// number of swaps are performed on p and -1 otherwise.
	private static Matrix[] detPLU(Matrix m) throws Exception{
		if(m.nRows != m.nCols)
			throw new NonSquareMatrixException("The matrix m (" + m.nRows + "," + m.nCols + ") must be square!");

		int N = m.nRows;

		double subAmt; 
		double focus; 

		Matrix lower = Matrix.eye(N);
		Matrix upper = new Matrix(m); 
		Matrix p = Matrix.eye(N);

		Stack<int[]>swapList = new Stack<int[]>();
				
		
		// Returns permutation matrix at index 0, 
		// lower matrix at index 1,
		// upper matrix at index 2.
		// numSwaps at index 3.   
		Matrix[] ret = new Matrix[4];

		// Essentially just performing REF operations here
		// But we are logging alpha in the lower triangular matrix
		// as well as swapping lower triangular matrix rows when necessary.
		for(int i = 0; i < N; i++){

			int max = i; 

			for(int j = i+1; j < N; j++){

				if(Math.abs(upper.mat[j][i]) > Math.abs(upper.mat[max][i])){
					max = j; 	
				}

			}

			// Exchanging Rows such that max is at position i in lower triangular matrix
			if(max != i && i > 0){
				for(int k = 0; k < i; k++){
					
					double temp2 = lower.mat[i][k];
					lower.mat[i][k] = lower.mat[max][k];
					lower.mat[max][k] = temp2; 
					
				}	
			}
			

			if(i != max){
				int[] swapPair = {i, max};
				swapList.push(swapPair);
			}

			double[] temp = upper.mat[i];
			upper.mat[i] = upper.mat[max];
			upper.mat[max] = temp;
			
			if(Math.abs(upper.mat[i][i]) <= EPSILON)
				throw new SingularMatrixException("Matrix is singular.");

			for(int k = i+1; k < N; k++){
				
				double alpha = upper.mat[k][i] / upper.mat[i][i];
				
				// Logging alpha into Lower Matrix
				lower.mat[k][i] = alpha;  
				

				for(int j = i; j < N; j++){
					upper.mat[k][j] -= (alpha * upper.mat[i][j]);
				}

				
			}
		}

		ret[1] = lower; 
		ret[2] = upper;
		ret[3] = new Matrix(1,1);
		ret[3].mat[0][0] = swapList.size()%2 == 0 ? 1 : -1; // Determinant of Identity matrix with swapped rows is dependent on
															// # of row swaps. Even # of row swaps -> -1, odd -> 1

		// Applying stored swapping operations
		// to permutation matrix such that it encodes
		// all swaps applied to upper and lower. 
		while(!swapList.empty()){
			int[] sle = swapList.pop();

			double[] t = p.mat[sle[0]];
			p.mat[sle[0]] = p.mat[sle[1]];
			p.mat[sle[1]] = t;
		}

		//Matrix.mul(Matrix.mul(p, lower), upper).printMatrix();

		ret[0] = p;

		return ret; 
	}

	// Returns Determinant of Matrix
	// Note: Uses fact that det(m) = det(P)*det(L)*det(U)
	public static double det(Matrix m) throws Exception{

		if(m.nRows != m.nCols)
			throw new NonSquareMatrixException("The matrix m (" + m.nRows + "," + m.nCols + ") must be square!");
			

		Matrix[] plu; 
		Matrix upperTri; 
		double prod; 

		try{ 
			
			plu = detPLU(m);
			upperTri = new Matrix(plu[2]);	
			prod = plu[3].mat[0][0]; // detPLU[3] is 1x1 matrix containing determinant of P. 
		
		} catch(SingularMatrixException e){
			return 0; // Singular Matrix always has determinant = 0
		}

		for(int i = 0; i < m.nRows; i++)
			prod *= upperTri.mat[i][i];

		return prod; 
	}

	// Returns trace of Matrix 
	// (sum of diagonal of square matrix)
	public static double trace(Matrix m) throws Exception{

		if(m.nRows != m.nCols)
			throw new NonSquareMatrixException("The matrix m (" + m.nRows + "," + m.nCols + ") must be square!");

		double trace = 0; 

		for(int i = 0; i < m.nRows; i++)
				trace += m.mat[i][i];

		return trace; 

	}

	// Returns flattened matrix e.g. 1D Vector
	// of all elements in m.
	public static Matrix flatten(Matrix m){
		Matrix ret = new Matrix(1, m.nRows*m.nCols);

		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				ret.mat[0][i*m.nCols + j] = m.mat[i][j];

		return ret; 
	}

	// Returns reshaped Matrix via flattening it and linearly
	// appending elements to new matrix array.
	public static Matrix reshape(Matrix m, int rows, int cols) throws Exception{
			
		if(m.nRows * m.nCols != rows*cols)
			throw new IncompatibleDimensionsException("Matrix must have equal number of elements as desired reshape: (" + m.nRows + "," + m.nCols + ") = " + m.nRows * m.nCols + ", (" + rows + "," + cols + ") = " + rows*cols);

		Matrix ret = new Matrix(rows, cols);
		Matrix mFlat = Matrix.flatten(m);

		int k = 0; 
		for(int i = 0; i < rows; i++){

			for(int j = 0; j < cols; j++){

				ret.mat[i][j] = mFlat.mat[0][k];
				k++; 

			}
		}

		return ret; 

	}

	// Returns the transpose of a matrix
	// e.g. flipped i,j values
	public static Matrix transpose(Matrix m){
		
		Matrix tran = new Matrix(m.nCols, m.nRows);

		for(int i = 0; i < m.nRows; i++){

			for(int j = 0; j < m.nCols; j++){

				tran.mat[j][i] = m.mat[i][j]; 

			}
		}

		return tran;
	}

	// Returns minor of Matrix e.g. matrix m ignoring specified row and col.
	private static Matrix minor(Matrix m, int row, int col){

		Matrix min = new Matrix(m.nRows-1, m.nCols-1);
		int k = 0, l = 0; 

		for(int i = 0; i < m.nRows; i++){

			if (i == row) continue; 

			for(int j = 0; j < m.nCols; j++){

				if (j == col) continue; 

				min.mat[k][l] = m.mat[i][j];

				l = ++l % min.nCols;
				
				if(l == 0) k++; // l becomes 0 when we've filled all columns in a row.  

			}

		}

		return min; 

	}

	// Returns cofactor of Matrix m 
	// E.g. Matrix containing determinants of all minors in m
	public static Matrix cofactor(Matrix m) throws Exception{

		Matrix c = new Matrix(m.nRows, m.nCols);

		for(int i = 0; i < m.nRows; i++){

			for(int j = 0; j < m.nCols; j++){

				c.mat[i][j] = Math.pow(-1, i+j) * Matrix.det(Matrix.minor(m,i,j));

			}
		}

		return c; 
	}

	// Returns transpose of the cofactor matrix. 
	public static Matrix tranCof(Matrix m) throws Exception{ 
		Matrix c = new Matrix(m.nRows, m.nCols);

		for(int i = 0; i < m.nRows; i++){

			for(int j = 0; j < m.nCols; j++){

				c.mat[j][i] = Math.pow(-1, i+j) * Matrix.det(Matrix.minor(m,i,j));

			}
		}

		return c; 
	}

	// Element-wise multiplication to m
	public static Matrix mul(double a, Matrix m){
		
		Matrix ret = new Matrix(m);

		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				ret.mat[i][j] *= a; 

		return ret;
	}

	// Element-wise divison to m
	public static Matrix div(double a, Matrix m){
		
		Matrix ret = new Matrix(m);

		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				ret.mat[i][j] /= a; 

		return ret;
	}

	// Element-wise addition to m
	public static Matrix add(double a, Matrix m){
		
		Matrix ret = new Matrix(m);

		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				ret.mat[i][j] += a; 

		return ret;
	}

	// Element-wise subtraction to m
	public static Matrix sub(double a, Matrix m){
		
		Matrix ret = new Matrix(m);

		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				ret.mat[i][j] -= a; 

		return ret;
	}

	// Adding two matrices to eachother
	public static Matrix add(Matrix ...args) throws Exception{

		Matrix res = new Matrix(args[0]);

		for(int m = 1; m < args.length; m++){
			
			if(res.nRows != args[m].nRows || res.nCols != args[m].nCols)
				throw new IncompatibleDimensionsException("Incompatible Dimensions: (" + res.nRows + ", " + res.nCols + ") and (" + args[m].nRows + ", " + args[m].nCols + ")" );
				

			for(int i = 0; i < res.nRows; i++)
				for(int j = 0; j < res.nCols; j++)
					res.mat[i][j] += args[m].mat[i][j];
		}
		

		return res; 

	}

	// Subtract two or more matrices from first matrix argument.
	public static Matrix sub(Matrix ...args) throws Exception{

		Matrix res = new Matrix(args[0]);

		for(int m = 1; m < args.length; m++){
			
			if(res.nRows != args[m].nRows || res.nCols != args[m].nCols)
				throw new IncompatibleDimensionsException("Incompatible Dimensions: (" + res.nRows + ", " + res.nCols + ") and (" + args[m].nRows + ", " + args[m].nCols + ")" );
				

			for(int i = 0; i < res.nRows; i++)
				for(int j = 0; j < res.nCols; j++)
					res.mat[i][j] -= args[m].mat[i][j];
		}
		

		return res; 

	}

	// Returns the inverse of a square matrix M
	// Use the adjugate method: M^-1 = 1/det(M) * C^T
	// Where C^T is the transpose of the cofactor matrix of M. 
	public static Matrix inv(Matrix m) throws Exception{
		return Matrix.mul(1/Matrix.det(m), Matrix.tranCof(m));
	}

	// Returns the p-norm of a Vector m
	public static double norm(Matrix m, double p){

		// Vector Norm
		double res = 0; 

		// Row Vector
		if(m.nRows == 1){

			for(int i = 0; i < m.nCols; i++)
				res += Math.pow(Math.abs(m.mat[0][i]), p);

			res = Math.pow(res, 1/p);

		}

		// Column Vector
		else if(m.nCols == 1){

			for(int i = 0; i < m.nRows; i++)
				res += Math.pow(Math.abs(m.mat[i][0]), p);

			res = Math.pow(res, 1/p);
		}

		else{ 
			throw new RuntimeException("P-Norm is not valid for multi-dimensional matrices. Did you want Matrix.norm(Matrix)?");
		}

		return res; 	 
	}

	// Returns the Frobenius Norm
	// (2-norm for matrices)
	public static double norm(Matrix m){

		double res = 0; 
		for(int i = 0; i < m.nRows; i++)
			for(int j = 0; j < m.nCols; j++)
				res += Math.pow(m.mat[i][j],2);

		return Math.sqrt(res);
	}


	// Matrix Chain Multiplication
	// Efficiently Multiply several compatible matrices in a single call.
	// More efficient for repeated products of large matrices.
	public static Matrix MCM(Matrix ...args) throws Exception{

		// Compute Parenthesizations
		int[] p = new int[args.length+1];

		p[0] = args[0].nRows;
		p[1] = args[0].nCols; 		
		for(int i = 1; i < args.length; i++){

			if(args[i].nRows != p[i])
				throw new IncompatibleDimensionsException("Incompatible Dimensions: (" + p[0] + ", " + p[i] + ") and (" + args[i].nRows + ", " + args[i].nCols + ")" );

			p[i+1] = args[i].nCols;
		}

		int[][] parens = Matrix.MatrixChainOrder(p);
		
		// Perform Product with optimal parenthesizations
		Matrix res = Matrix.MCMHelper(parens[1][args.length], parens, args);

		return res; 


	}

	// Recursively iterates through parenthesization table to seperate Matrices into their most
	// efficient products.
	private static Matrix MCMHelper(int paren, int[][] parens, Matrix ...args) throws Exception{

		// Base Case 1: Best way to multiply two matrices 
		// is to just multiply them together.
		if(args.length == 2)
			return Matrix.mul(args[0], args[1]);
		
		// Base Case 2: If there's only one matrix, return it.
		if(args.length == 1)
			return args[0];
		
		// Ensures that paren value is transformed for smaller argument array.
		if(paren > args.length)
			paren -= args.length; 

		// 0 isn't a valid parenthesization, replace with 1. 
		else if(paren == args.length)
			paren = 1; 

		// Compute optimal parenthesizations for left and right products.
		Matrix[] leftP = new Matrix[paren];
		Matrix[] rightP = new Matrix[args.length-paren];

		for(int i = 0; i < leftP.length; i++){
			leftP[i] = args[i];
		}

		for(int i = 0; i < rightP.length; i++){
			rightP[i] = args[i+leftP.length];
		}


		Matrix left = MCMHelper(parens[1][paren],parens,leftP); // Left of Parenthesization
		Matrix right = MCMHelper(parens[paren+1][args.length],parens,rightP); // Right of Parenthesization

		return Matrix.mul(left,right);

	}

	// Produces table of most efficient matrix products
	// Matrix Mi has dimensions pi-1, pi
	private static int[][] MatrixChainOrder(int[] p){

		int N = p.length; 

		int[][] m = new int[N][N];
		int[][] s = new int[N][N];

		int j = 0;

		for(int l = 2; l < N; l++){ // Chain Length

			for(int i = 1; i < N-l+1; i++){

				j = i+l-1;
				m[i][j] = Integer.MAX_VALUE;

				for(int k = i; k <= j-1; k++){

					int q = m[i][k] + m[k+1][j] + p[i-1] * p[k] * p[j];

					if(q < m[i][j]){
						m[i][j] = q;
						s[i][j] = k; 	
					}  
				} 

			}
		}

		return s;

	}

}
