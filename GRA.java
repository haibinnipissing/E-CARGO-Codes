/**
 * An implementation of Group Role Assignment;
 * Adapted from the implementation of the classic hungarian algorithm for the assignment problem.
 *  by Gary Baker (GPL v3), 2007
 * Please cite:
 * [1]	H. Zhu, M.C. Zhou, and R. Alkins, “Group Role Assignment via a Kuhn-Munkres Algorithm-based Solution”, IEEE Trans. on Systems, Man, and Cybernetics, Part A: Systems and Humans, vol. 42, no. 3, May 2012, pp. 739-750.
 * [2]  H. Zhu, and M. Zhou, “Role-Based Collaboration and its Kernel Mechanisms,” IEEE Trans. on Systems, Man, and Cybernetics, Part C: Applications and Reviews, vol. 36, no. 4, 2006, pp. 578-589.  
 */

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.lang.Integer;
 class RatedMunkres {
    static public int[][] computeAssignments(double[][] matrix) {
        // subtract minumum value from rows and columns to create lots of zeroes
        reduceMatrix(matrix);

        // non negative values are the index of the starred or primed zero in the row or column
        int[] starsByRow = new int[matrix.length]; Arrays.fill(starsByRow,-1);
        int[] starsByCol = new int[matrix[0].length]; Arrays.fill(starsByCol,-1);
        int[] primesByRow = new int[matrix.length]; Arrays.fill(primesByRow,-1);

        // 1s mean covered, 0s mean not covered
        int[] coveredRows = new int[matrix.length];
        int[] coveredCols = new int[matrix[0].length];

        // star any zero that has no other starred zero in the same row or column
        initStars(matrix, starsByRow, starsByCol);
        coverColumnsOfStarredZeroes(starsByCol,coveredCols);

        while (!allAreCovered(coveredCols)) {

            int[] primedZero = primeSomeUncoveredZero(matrix, primesByRow, coveredRows, coveredCols);

            while (primedZero == null) {
                // keep making more zeroes until we find something that we can prime (i.e. a zero that is uncovered)
                makeMoreZeroes(matrix,coveredRows,coveredCols);
                primedZero = primeSomeUncoveredZero(matrix, primesByRow, coveredRows, coveredCols);
            }

            // check if there is a starred zero in the primed zero's row
            int columnIndex = starsByRow[primedZero[0]];
            if (-1 == columnIndex){

                // if not, then we need to increment the zeroes and start over
                incrementSetOfStarredZeroes(primedZero, starsByRow, starsByCol, primesByRow);
                Arrays.fill(primesByRow,-1);
                Arrays.fill(coveredRows,0);
                Arrays.fill(coveredCols,0);
                coverColumnsOfStarredZeroes(starsByCol,coveredCols);
            } else {

                // cover the row of the primed zero and uncover the column of the starred zero in the same row
                coveredRows[primedZero[0]] = 1;
                coveredCols[columnIndex] = 0;
            }
        }

        // ok now we should have assigned everything
        // take the starred zeroes in each column as the correct assignments

        int[][] retval = new int[matrix.length][];
        for (int i = 0; i < starsByCol.length;  i++) {
            retval[i] = new int[]{starsByCol[i],i};
        }
        return retval;
    }

    static private boolean allAreCovered(int[] coveredCols) {
        for (int covered : coveredCols) {
            if (0 == covered) return false;
        }
        return true;
    }


    /**
     * the first step of the hungarian algorithm
     * is to find the smallest element in each row
     * and subtract it's values from all elements
     * in that row
     *
     * @return the next step to perform
     */
    static private void reduceMatrix(double[][] matrix) {

        for (int i = 0; i < matrix.length; i++) {

            // find the min value in the row
            double minValInRow = Double.MAX_VALUE;
            for (int j = 0; j < matrix[i].length; j++) {
                if (minValInRow > matrix[i][j]) {
                    minValInRow = matrix[i][j];
                }
            }

            // subtract it from all values in the row
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] -= minValInRow;
            }
        }

        for (int i = 0; i < matrix[0].length; i++) {
            double minValInCol = Double.MAX_VALUE;
            for (int j = 0; j < matrix.length; j++) {
                if (minValInCol > matrix[j][i]) {
                    minValInCol = matrix[j][i];
                }
            }

            for (int j = 0; j < matrix.length; j++) {
                matrix[j][i] -= minValInCol;
            }

        }

    }

    /**
     * init starred zeroes
     *
     * for each column find the first zero
     * if there is no other starred zero in that row
     * then star the zero, cover the column and row and
     * go onto the next column
     *
     * @param costMatrix
     * @param starredZeroes
     * @param coveredRows
     * @param coveredCols
     * @return the next step to perform
     */
    static private void initStars(double costMatrix[][], int[] starsByRow, int[] starsByCol) {


        int [] rowHasStarredZero = new int[costMatrix.length];
        int [] colHasStarredZero = new int[costMatrix[0].length];

        for (int i = 0; i < costMatrix.length; i++) {
            for (int j = 0; j < costMatrix[i].length; j++) {
                if (0 == costMatrix[i][j] && 0 == rowHasStarredZero[i] && 0 == colHasStarredZero[j]) {
                    starsByRow[i] = j;
                    starsByCol[j] = i;
                    rowHasStarredZero[i] = 1;
                    colHasStarredZero[j] = 1;
                    break; // move onto the next row
                }
            }
        }
    }


    /**
     * just marke the columns covered for any coluimn containing a starred zero
     * @param starsByCol
     * @param coveredCols
     */
    static private void coverColumnsOfStarredZeroes(int[] starsByCol, int[] coveredCols) {
        for (int i = 0; i < starsByCol.length; i++) {
            coveredCols[i] = -1 == starsByCol[i] ? 0 : 1;
        }
    }


    /**
     * finds some uncovered zero and primes it
     * @param matrix
     * @param primesByRow
     * @param coveredRows
     * @param coveredCols
     * @return
     */
    static  private int[] primeSomeUncoveredZero(double matrix[][], int[] primesByRow,
                                       int[] coveredRows, int[] coveredCols) {


        // find an uncovered zero and prime it
        for (int i = 0; i < matrix.length; i++) {
            if (1 == coveredRows[i]) continue;
            for (int j = 0; j < matrix[i].length; j++) {
                // if it's a zero and the column is not covered
                if (0 == matrix[i][j] && 0 == coveredCols[j]) {

                    // ok this is an unstarred zero
                    // prime it
                    primesByRow[i] = j;
                    return new int[]{i,j};
                }
            }
        }
        return null;

    }

    /**
     *
     * @param unpairedZeroPrime
     * @param starsByRow
     * @param starsByCol
     * @param primesByRow
     */
    static  private void incrementSetOfStarredZeroes(int[] unpairedZeroPrime, int[] starsByRow, int[] starsByCol, int[] primesByRow) {

        // build the alternating zero sequence (prime, star, prime, star, etc)
        int i, j = unpairedZeroPrime[1];

        Set<int[]> zeroSequence = new LinkedHashSet<int[]>();
        zeroSequence.add(unpairedZeroPrime);
        boolean paired = false;
        do {
            i = starsByCol[j];
            paired = -1 != i && zeroSequence.add(new int[]{i,j});
            if (!paired) break;

            j = primesByRow[i];
            paired = -1 != j && zeroSequence.add(new int[]{ i, j });

        } while (paired);


        // unstar each starred zero of the sequence
        // and star each primed zero of the sequence
        for (int[] zero : zeroSequence) {
            if (starsByCol[zero[1]] == zero[0]) {
                starsByCol[zero[1]] = -1;
                starsByRow[zero[0]] = -1;
            }
            if (primesByRow[zero[0]] == zero[1]) {
                starsByRow[zero[0]] = zero[1];
                starsByCol[zero[1]] = zero[0];
            }
        }

    }


    static  private void makeMoreZeroes(double[][] matrix, int[] coveredRows, int[] coveredCols) {

        // find the minimum uncovered value
        double minUncoveredValue = Double.MAX_VALUE;
        for (int i = 0; i < matrix.length; i++) {
            if (0 == coveredRows[i]) {
                for (int j = 0; j < matrix[i].length; j++) {
                    if (0 == coveredCols[j] && matrix[i][j] < minUncoveredValue) {
                        minUncoveredValue = matrix[i][j];
                    }
                }
            }
        }

        // add the min value to all covered rows
        for (int i = 0; i < coveredRows.length; i++) {
            if (1 == coveredRows[i]) {
                for (int j = 0; j < matrix[i].length; j++) {
                        matrix[i][j] += minUncoveredValue;
                }
            }
        }

        // subtract the min value from all uncovered columns
        for (int i = 0; i < coveredCols.length; i++) {
            if (0 == coveredCols[i]) {
                for (int j = 0; j < matrix.length; j++) {
                    matrix[j][i] -= minUncoveredValue;
                }
            }
        }
    }
	public static double RatedAssign(int [] L,  double [][] Q, int [][] T, int m, int n, double th) {
		double v=0.0f;
		int cnt=0, LL[]=new int [m];
		double Q1[][]= new double [m][n];
		//Check if it can be a square matrix.
		for (int i = 0; i<n; i++) cnt +=L[i];
		if (cnt > m) return 0.0f;//Not enough agents.
		//Adjust Q with the threshold.
		for (int i = 0; i < m; i++)
			for (int j =0; j< n; j++) 
				if (Q[i][j]<=th) Q1[i][j]=-m*m;
				else Q1[i][j]=Q[i][j];
		double CC[]= new double [n]; //		CC[] is the numbers of qualified agents for roles.
		int D[]= new int [n];		//D is a vector for the difference between the number required current agents and the number of actual current agents for each role.   
		for (int j =0; j< n; j++)
		{	for (int i = 0; i < m; i++)			{		if (Q1[i][j]>th)	 CC[j]=CC[j]+1;			}
			D[j]=(int)( L[j]-CC[j]);		if (D[j]>0) return 0.0f;//One role has not enough agents..
		}
		//Create the index vector.
		int index=0;
		for (int j = 0; j<n; j++)
				for (int k = 0; k<L[j]; k++) LL[index++] =j;
		for (int k = index; k < m; k++)//fill the left columns (roles).
			LL[k]=n;

		double [][] M = new double [m][m];
		for (int i = 0; i<m; i++)
		{ 	index =0;
			for (int j = 0; j<n; j++)
			{
				for (int k = 0; k<L[j]; k++)
					M[i][index++]=1-Q1[i][j];
			}
		}
		int [][] N = computeAssignments(M);
		//Obtaing the matrix T.
		for (int i = 0; i<m; i++)
		{ 	for (int j = 0; j<n; j++)
			{					T[i][j]=0;
			}
		}		
		for (int i = 0; i < N.length; i++)
			if (LL[N[i][1]]< n)
			{
				if (Q1[N[i][0]][LL[N[i][1]]]>0.0f) T[N[i][0]][LL[N[i][1]]]=1;
				else return -1.0f;
			}
		for (int i = 0; i<m; i++)
		 	for (int j = 0; j<n; j++)
				v += Q[i][j]*T[i][j];
		return v;
	}
}


public class GRA {
	 public static  void printDMatrix (double [][]x, int m, int n){
			DecimalFormat tw = new DecimalFormat("0.00");
			for (int i = 0; i < m; i++)
			{	for (int j =0; j< n; j++)
				{
				System.out.print (tw.format(x[i][j]));		System.out.print (" ");
				}
			System.out.println ();
			}
			System.out.println ();
		}	
		 public static  void printIMatrix (int [][]x, int m, int n){
			DecimalFormat tw = new DecimalFormat("0");
			for (int i = 0; i < m; i++)
			{	for (int j =0; j< n; j++)
				{
				System.out.print (tw.format(x[i][j]));		System.out.print (" ");
				}
			System.out.println ();
			}
			System.out.println ();
		}	

public static int sigmaL(int []L){
	int total=0;
	for(int j=0; j<L.length; j++) 
		total+=L[j];
	return total;
}
	
/*	public static void main(String[] args)
	{
		Random generator = new Random();
		DecimalFormat df = new DecimalFormat("0.00");
		int m = 1000, n = 120; 
		int []L=new int [n];
		double [][]Q=new double [m][n];		
		try 
		{
			BufferedWriter out = new BufferedWriter(new FileWriter(("Result"), true));
			out.write("Q: \n");	// Random Q
			for(int r=0; r<m; r++)
			{
				for(int c=0; c<n; c++)
				{
					Q[r][c] = generator.nextDouble();
					out.write(df.format(Q[r][c]) + "	");
				}
				out.write("\n");
			}			
			out.write("\n");	
			out.write("\nL: \n");	//Random L
			for (int i =0; i<n; i++) 	
			{ 	
				L[i] = generator.nextInt(m/n)+1;
				out.write(L[i] + "	");
			}
			out.write("\n");
			out.close();
		}
		catch (IOException e) {System.out.println ("Error in writing into a file!");}
		//TEST parameters:
		int[][] T = new int[m][n];
		long t11 = System.nanoTime();
		double v1=RatedMunkres.RatedAssign(L, Q, T, m, n, 0);
		long t12 = System.nanoTime();
		double diff1 = (double)(t12-t11)/1000000;
		printDMatrix (Q, m, n);
		printIMatrix (T, m, n);
		System.out.print("L=[");	
		for (int j=0; j<n; j++) {System.out.print(L[j]+" ");}	System.out.println("]");
		System.out.println ("Total GRA ="+v1+" "+"Time = "+diff1+"ms");
		System.out.println();
		return;
	}	*/
public static void main(String[] args)
{
	Random generator = new Random();
	DecimalFormat df = new DecimalFormat("0.00");
	int m =6;
	int n =3;
	int L[]={2, 1, 2};
//	int L[]={1, 2, 2, 1};
/*	double [][]Q={
			{0.56,0.35,0.80,0.62},
			{0.49,0.09,0.33,0.58},
			{0.91,0.31,0.90,0.15},
			{0.58,0.64,0.24,0.0},
			{0.38,0.54,0.72,0.20},
			{0.85,0.34,0.43,0.18},};
*/	
	double [][]Q={
			{	0.56, 0.35, 0.50}, 
	{0.92, 0.89, 0.30}, 
	{0.09, 0.26, 0.20}, 
	{0.14, 0.32, 0.45}, 
	{0.42, 0.18, 0.49}, 
	{0.86, 0.20, 1.12}};
	int[][] T = new int[m][n];
	long t11 = System.nanoTime();
	double v1=RatedMunkres.RatedAssign(L, Q, T, m, n, 0);
	long t12 = System.nanoTime();
	double diff1 = (double)(t12-t11)/1000000;
	printDMatrix (Q, m, n);
	printIMatrix (T, m, n);
	double W[] = {1, 1, 1};
	double [][]Q1=new double [m][n];
	int[][] T1 = new int[m][n];
	for (int i =0; i< m; i++)
		for (int j =0; j <n; j++)
			Q1[i][j]=Q[i][j]*W[j];
	double v2=RatedMunkres.RatedAssign(L, Q1, T1, m, n, 0);
	printDMatrix (Q1, m, n);
	printIMatrix (T1, m, n);
	System.out.print("L=[");	
	for (int j=0; j<n; j++) {System.out.print(L[j]+" ");}	System.out.println("]");
	System.out.println ("Total GRA ="+v1+" "+v2+" "+"Time = "+diff1+"ms");
	System.out.println();
	return;
}	

/*//		int m = 13;
int m =12;
int n =4;
int L[]={1,2,4,2};
double [][]Q={
		{0.18,0.82,0.29,0.01},
		{0.35,0.80,0.58,0.35},
		{0.84,0.85,0.86,0.36},
		{0.96,0.51,0.45,0.64},
//		{0.92,0,0,0},
		{0.22,0.33,0.68,0.33},
		{0.96,0.50,0.10,0.73},
		{0.25,0.18,0.23,0.39},
		{0.56,0.35,0.80,0.62},
		{0.49,0.09,0.33,0.58},
		{0.38,0.54,0.72,0.20},
		{0.91,0.31,0.34,0.15},
		{0.85,0.34,0.43,0.18},};
//		{0.44,0.06,0.66,0.37}};
//		{0.93,0.06,0.66,0.37}};
*/

}
