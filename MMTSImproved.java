/*
 * This program is for Temporal Role Transfer with Strong Restriction with KM algorithm.
Please cite:
[1]H. Zhu, and M. Zhou, “M–M Role-Transfer Problems and Their Solutions,” IEEE Trans. on Systems, Man, and Cybernetics, Part A: Systems and Humans, vol. 39, no. 2, Mar. 2009, pp. 448-459. 
[2]H. Zhu, and M. Zhou, “Role Transfer Problems and Algorithms,” IEEE Trans. on Systems, Man, and Cybernetics, Part A: Systems and Humans, vol. 38, no. 6, Nov. 2008, pp. 1442-1450. 
[3]H. Zhu, and M.C. Zhou, “Efficient Role Transfer Based on Kuhn–Munkres Algorithm”, IEEE Trans. on Systems, Man, and Cybernetics, Part A: Systems and Humans, vol. 42, no.2, 2012, pp. 491 - 496. 
*/

import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.Random;
import java.util.Set;

class TestResult
{
	public int scale;
	public double time;
	public TestResult()
	{
		scale=0;
		time =0;
	}
}
class Scheme
{
public int [][] SS;
public int v;
public Scheme (int n)
{
	SS= new int [n][n];
	v = 0;
}
}
class SequencedScheme
{
public int [][][]SS;
public int [] SI;
public int v;
public SequencedScheme (int n)
{
	SS= new int [n][][];
	SI = new int [n];
	v = 0;
}
}

public class MMTSImproved {
    static public int[][] computeAssignments(float[][] matrix) {


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
    static private void reduceMatrix(float[][] matrix) {

        for (int i = 0; i < matrix.length; i++) {

            // find the min value in the row
            float minValInRow = Float.MAX_VALUE;
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
            float minValInCol = Float.MAX_VALUE;
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
    static private void initStars(float costMatrix[][], int[] starsByRow, int[] starsByCol) {


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
    static  private int[] primeSomeUncoveredZero(float matrix[][], int[] primesByRow,
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


    static  private void makeMoreZeroes(float[][] matrix, int[] coveredRows, int[] coveredCols) {

        // find the minimum uncovered value
        float minUncoveredValue = Float.MAX_VALUE;
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
	public static float SimpleAssign(int [] L,  int [][] Q, int [][] T, int m, int n, float th) {
		float v=0.0f;
		int cnt=0; int [] LL=new int [m];
		float Q1[][]= new float [m][n];
		//Check if it can be a sqaure matrix.
		for (int i = 0; i<n; i++) cnt +=L[i];
		if (cnt > m) return 0.0f;//Not enough agents.
		//Adjust Q with the threshold.
		for (int i = 0; i < m; i++)
			for (int j =0; j< n; j++) 
				if (Q[i][j]<th) Q1[i][j]=-m*m;
				else Q1[i][j]=Q[i][j];
		float CC[]= new float [n]; //		CC[] is the numbers of qualified agents for roles.
		int D[]= new int [n];		//D is a vector for the difference between the number required current agents and the number of actual current agents for each role.   
		for (int j =0; j< n; j++)
		{	for (int i = 0; i < m; i++)			{		if (Q1[i][j]>=th)	 CC[j]=CC[j]+1;			}
			D[j]=(int)( L[j]-CC[j]);		if (D[j]>0) return 0.0f;//One role has not enough agents..
		}
		//Create the index vector.
		int index=0;
		for (int j = 0; j<n; j++)
				for (int k = 0; k<L[j]; k++) LL[index++] =j;
		for (int k = index; k < m; k++)//fill the left columns (roles).
			LL[k]=n;

		float [][] M = new float [m][m];
		for (int i = 0; i<m; i++)
		{ 	index =0;
			for (int j = 0; j<n; j++)
			{
				for (int k = 0; k<L[j]; k++)
					M[i][index++]=1-Q1[i][j];
			}
			//for (int k = index; k < m; k++)//fill the left columns (roles).
				//M[i][k]=1;
		}

/*		
	for (int i = 0; i < m; i++)
			{System.out.print (LL[i]);System.out.print (" ");	}
		System.out.println ();
		
		for (int i = 0; i < m; i++)
		{	for (int j =0; j< m; j++)
			{	
			System.out.print (M[i][j]);	System.out.print (" ");	}System.out.println ();
		}
*/			int [][] N = computeAssignments(M);
/*
		for (int i = 0; i < N.length; i++)
		{	for (int j =0; j< N[0].length; j++)
			{	
			System.out.print (N[i][j]);	System.out.print (" ");	}System.out.println ();
		}
*/
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
//		System.out.println ("v="+v);
		//check if it is a workable group:
/*		System.out.println ("T:");
		for (int i = 0; i<m; i++)
		{ 	for (int j = 0; j<n; j++)
				{System.out.print(T[i][j]+" ");}
		System.out.println ();
		}
		System.out.println ();
		for (int j =0; j< n; j++)
		{	CC[j]=0;
			for (int i = 0; i < m; i++)			{			 if (Q[i][j]>th) CC[j]=CC[j]+T[i][j];			}
			D[j]=(int)(L[j]-CC[j]);		
			if (D[j]>0)  return -1.0f;//One role has not enough agents..
		}
*/
		//		System.out.println ();
		for (int i = 0; i<m; i++)
		 	for (int j = 0; j<n; j++)
				v += Q[i][j]*T[i][j];
		return v;
	}

	public static Scheme SplittingScheme(int [] L, int m, int n)
	{
		int l = L.length;
		int [] LTAG= new int [l];
		for (int i = 0; i< l; i++) LTAG[i]=0;
		Scheme S= new Scheme(n);
		for (int i=0; i<n; i++)
			for (int j=0; j<n; j++) S.SS[i][j]=-1;

		int s = 0, acc=0; 
	
		while (s<n)
		{  int v=0;
			int u=0;
			while ((u<n)&&(LTAG[u] !=0)) u++;
			if (u<n)
			{	if (v<n) S.SS[s][v++]=u; LTAG[u]=1;
				acc=L[u];
				int k = u+1;
				if (k<n)
					while (k<n)
					{	if ((LTAG[k]!=1) && ((acc+ L[k])<=m))
				  		{
							acc+= L[k];
							if (v<n) S.SS[s][v++]=k;
							LTAG[k]=1;//Put a tag to mean this role has been considered.
				  		}
						k++;
					}
				}
			else break;
		s++;
		}
		S.v = s;
		return  S;
	   	}

	public static SequencedScheme splitting (int [] L, int [][] M, int s,  int [][]S, int m, int n)
	{
		SequencedScheme  SM= new SequencedScheme (s);
		for (int i = 0; i<=s-1; i++) 
		{
		    int sl=0; 
			for (int r=0; r<=n-1; r++) 
				if (S[i][r]==-1) { sl = r; break;} 
			SM.SS[i]=new int[m][sl];
			SM.SI[i]=sl;
			for (int j = 0; j<=sl-1; j++)
			{	 
					for (int l =0; l <=m-1; l ++)
						SM.SS[i][l][j]= M[l][S[i][j]];				
			}
		}
		SM.v = s;
		return SM;
		}

	public static int NewTMMPS (int [] L, int [][] Mc, int [][]Mp, int [][][]RSM, int [][][]RSQ, int m, int n)
	{
	int [][] Mr = new int [m][n];
	for (int i =0; i<m; i++)
		for (int j =0; j<n; j++)
			Mr [i][j] = Mc[i][j] + Mp[i][j];

	int [] CC = new int [n];
	int [] D = new int [n];
	
	for (int j =0; j< n; j++) CC[j]=0;
			//Step 1: Preparation
	for (int j =0; j< n; j++)
	{	for (int i = 0; i < m; i++)
		{		CC[j]=CC[j]+Mr[i][j];
		}
	D[j]=L[j]-CC[j];
	if (D[j]>0) return 0;
	}
	int msplit =m;
	float res = 0;
	SequencedScheme SM; 
	do {
		Scheme Scheme = SplittingScheme(L, msplit, n);
		int s = Scheme.v;
		SM = splitting (L, Mr, s, Scheme.SS, m, n); 
		int [][] T= new int [m][m];
		for (int i = 0; i < s; i++)
		{	int []LN= new int [SM.SI[i]];
			for (int j =0; j<SM.SI[i]; j++) LN[j]=L[Scheme.SS[i][j]];
			res =SimpleAssign (LN, SM.SS[i], T, m, SM.SI[i], 1.0f);
			if (res<0)  {msplit--; break;};
			for (int j = 0; j < m; j++) 				for (int k = 0; k<n; k++) 				RSM[i][j][k]=0;
			for (int j = 0; j < m; j++)
			{	for (int k = 0; k<SM.SI[i]; k++)
				{	RSM[i][j][Scheme.SS[i][k]]= T[j][k];
					RSQ[i][j][Scheme.SS[i][k]] = Mr[j][k]-T[j][k];
				}
			for (int p = SM.SI[i]; p < n; p++)
									RSQ[i][j][p] = Mr[j][p]-T[j][p];				
			}
		}
	} while (res<0);
	 return 	SM.v;
	}
	
/*			public static void main(String[] args) {
				int C[][] =new int [10][5];
				int Q[][] =new int [10][5];
				//	......................
/*				int L[]={3, 3, 4, 4, 2};
			    C[0][0]=1;		C[1][0]=1;		C[2][2]=1;		C[3][0]=1; C[4][1]=1; C[5][2]=1;C[6][1]=1;C[7][1]=1;C[8][4]=1;C[9][2]=1;
			    Q[0][1]=1; Q[0][3]=1; Q[0][4]=1; Q[1][3]=1;	Q[2][0]=1; Q[2][1]=1; Q[2][3]=1; Q[2][4]=1; 
			    Q[3][1]=1;	Q[4][0]=1; Q[5][0]=1;Q[5][3]=1; Q[6][2]=1;
				int n1, n2;
				n1 = 10; n2 = 5;
  
 
				int L[]={4, 3, 4, 3};
			    C[0][0]=1;		C[1][2]=1;		C[2][0]=1;		C[3][0]=1; C[4][0]=1; C[5][1]=1;C[6][1]=1;C[7][1]=1;
			    
			    Q[0][1]=1; 
			    Q[1][0]=1; 
			    Q[2][2]=1; Q[3][1]=1;	Q[3][3]=1; 
			    Q[4][2]=1; Q[4][3]=1; 
			    Q[5][2]=1; 
			    Q[7][0]=1; Q[7][2]=1;Q[7][3]=1; 
				int n1, n2;
				n1 = 8; n2 = 4;
 				System.out.println ("Before Transferring...");
			System.out.println ("The Current Role Matrix:");
			for (int i = 0; i < n1; i++)
			{	for (int j =0; j< n2; j++)
				{
				System.out.print (C[i][j]);
				System.out.print (" ");
				}
			System.out.println ();
			}
			System.out.println ();
			System.out.println ("The Qualified Role Matrix:");
		for (int i = 0; i < n1; i++)
		{	for (int j =0; j< n2; j++)
			{
			System.out.print (Q[i][j]);
			System.out.print (" ");
			}
		System.out.println ();
		}
		System.out.println ();
//		int s=n2;
		int [][][]CL= new int [n2][n1][n2]; 
		int [][][]QL= new int [n2][n1][n2]; 
		int	res = NewTMMPS(L, C,  Q, CL, QL, n1, n2);
		if (res!=0) 
				System.out.println ("Success....Scale = "+res+"\n After Transferring...");
		else 
			System.out.println ("Fail...\n");
		System.out.println ("The Current Role Matrix List:");
		for (int k = 0; k < res; k++)
			{for (int i = 0; i < n1; i++)				
			  {	for (int j =0; j< n2; j++){	System.out.print (CL[k][i][j]);	System.out.print (" ");	} 		System.out.println ();		}	System.out.println ();
			}
      }
}
*/
	public static TestResult RandomTest() {
		Random generator = new Random();
		int n1, n2;
		n1 = 10; n2 = 5;
		n1 = 20; n2 = 10;
		n1 = 30; n2 = 15;
		n1 = 40; n2 = 20;
		n1 = 50; n2 = 25;
		int L[]= new int [n2];
		int C0[][] =new int [n1][n2];
		int Q[][] =new int [n1][n2];
		//	......................
		int sum =0;
		do {
		for (int i =0; i<n2; i++)
		{ 	L[i]=generator.nextInt(3)+2;	
			sum +=L[i];
		}
		}while (sum <n1);
		sum =0;
		int a[]=new int [n1];
		int r[]=new int [n2];
		
		do {
		for (int i = 0; i < n1; i++)			
			{if (a[i]==0)
				for (int j =0; j< n2; j++)
				{	if (r[j]<L[j]) C0[i][j]=generator.nextInt(2);				if (C0[i][j]==1) {a[i]=1; r[j] ++; sum ++; break;}			}			
			}
		}while (sum <n1);
//------
/*	System.out.print (sum+" ");
		for (int i = 0; i < n1; i++)
			System.out.print (a[i]+" ");
*/
		for (int j =0; j< n2; j++)
		{ int cnt=0;
			for (int i = 0; i < n1; i++)
			{	if (C0[i][j]==1) cnt++;
				if (cnt==L[j])
				{	for (int k = i+1; k < n1; k++)	C0[k][j]=0;
				break;
				}
			}
		}
//-------------
	for (int i = 0; i < n1; i++)
		{	for (int j =0; j< n2; j++)
			{
			if (C0[i][j]==0)			Q[i][j]=generator.nextInt(2);
			}
		}
		for (int j =0; j< n2; j++)
		{ int cnt=0;
			for (int i = 0; i < n1; i++)
			{	if (Q[i][j]==1) cnt++;
				if (cnt>=L[j])
				{	for (int k = i+1; k < n1; k++)	Q[k][j]=0;
				break;
				}
			}
		}

	System.out.println ("Before Transferring...");
		System.out.println ("The Current Role Matrix:");
		
		for (int i =0; i<n2; i++)
			System.out.print (L[i]+" ");
		System.out.println ();

		for (int i = 0; i < n1; i++)
		{	for (int j =0; j< n2; j++)
			{
			if (C0[i][j]==1||Q[i][j]==1) System.out.print (1);
			else System.out.print (0);
			System.out.print (" ");
			}
		System.out.println ();
		}
		System.out.println ();
/*		System.out.println ("The Qualified Role Matrix:");
	for (int i = 0; i < n1; i++)
	{	for (int j =0; j< n2; j++)
		{
		System.out.print (Q[i][j]);
		System.out.print (" ");
		}
	System.out.println ();
	}
	System.out.println ();
*/
		TestResult tr = new TestResult();
	int [][][]CL= new int [n2][n1][n2]; 
	int [][][]QL= new int [n2][n1][n2]; 
	long t1 = System.nanoTime();
	tr.scale = NewTMMPS(L, C0,  Q, CL, QL, n1, n2);
	long t2 = System.nanoTime();
	double diff = (double)(t2-t1)/1000000;
	tr.time = diff;
	System.out.println ("Time = "+diff+"ms");
	System.out.println ("Scale ="+tr.scale);
	System.out.println ("The Current Role Matrix List:");
	for (int k = 0; k < tr.scale; k++)
		{for (int i = 0; i < n1; i++)
			{	for (int j =0; j< n2; j++)						{						System.out.print (CL[k][i][j]);						System.out.print (" ");						}
			System.out.println ();
			}
			System.out.println ();
		}
	return tr;
	}
	public static void main(String[] args) {
		int num = 300;
		TestResult []TR=new TestResult [num]; 
		for (int i = 0; i<num; i++)
		{	TR [i]=RandomTest();
		}
		for (int i = 0; i<num; i++)
		{System.out.print (TR[i].scale);
		System.out.println (" "+TR[i].time);
		}
	}
}
