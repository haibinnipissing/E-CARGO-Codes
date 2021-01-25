/**
 * An implementation of the KMB algorithm to solve Group Multi-Role Assignment;
 * Adapted from the implementation of the classic hungarian algorithm for the assignment problem.
*  by Gary Baker (GPL v3), 2007
 * @author Dongning Liu, Siqin Zhang, and Haibin Zhu 2009.
 * Please cite:
 *[1] H. Zhu, D. Liu, S. Zhang, S. Teng, and Y. Zhu, “Solving the Group Multi-Role Assignment Problem by Improving the ILOG Approach”, IEEE Trans. on Systems, Man, and Cybernetics: Systems, vol. 47, no. 12, Dec. 2017, pp. 3418-3424. 
 *[2] H. Zhu, D. Liu, S. Zhang, Y. Zhu, L. Teng, S. Teng, “Solving the Many to Many Assignment Problem by Improving the Kuhn-Munkres Algorithm with Backtracking”, Theoretical Computer Science, vol. 618, March 2016, pp. 30-41. 
 */
import java.text.DecimalFormat;
import java.util.*;


class TestResult
{
	public boolean TF;
	public double time;
	public float count;
	public TestResult()
	{
		TF=false;
		time =0;
		count=0.0f;
	}
}
public class KMB {

    static public int[][] computeAssignments(float[][] matrix,int [] RC,int [] LL,int n,int[] L,int[] BA) {

        reduceMatrix(matrix);
        

        int[] starsByRow = new int[matrix.length]; Arrays.fill(starsByRow,-1);
        int[] starsByCol = new int[matrix[0].length]; Arrays.fill(starsByCol,-1);
        int[] primesByRow = new int[matrix.length]; Arrays.fill(primesByRow,-1);
        int[][] infinity=new int[matrix.length][matrix[0].length];
        int[] coveredRows = new int[matrix.length];
        int[] coveredCols = new int[matrix[0].length];

        initStars(matrix, starsByRow, starsByCol,RC,LL,infinity,n,L,BA);   
        coverColumnsOfStarredZeroes(starsByCol,coveredCols);

        while (!allAreCovered(coveredCols)) {
  
            int[] primedZero = primeSomeUncoveredZero(matrix, primesByRow, coveredRows, coveredCols,starsByRow,RC,LL,infinity);

            while (primedZero == null) {
                makeMoreZeroes(matrix,coveredRows,coveredCols,infinity);
                primedZero = primeSomeUncoveredZero(matrix, primesByRow, coveredRows, coveredCols,starsByRow,RC,LL,infinity);
            }
          
            int i=primedZero[0];
            int j=primedZero[1];
            if(LL[j]<n&&L[LL[j]]>1&&BA[RC[i]]>1){
                int t=j;
                int s=i;
                while(t>=0){
                	if(LL[t]==LL[j])
                		t--;
                	else break;
                }
                while(s>=0){
                	if(RC[s]==RC[i])
                		s--;
                	else break;
                }
                for(int k=t+1;LL[k]<=LL[j];)
                {
                 if(LL[k]==LL[j]&&k!=j)
                 for(int h=s+1;RC[h]<=RC[i];)
                 {
                    if(RC[h]==RC[i]&&h!=i&&infinity[h][k]==0)
                    	infinity[h][k]=-1;
                    h++;
                    if(h==matrix.length) break;
                 }
                 k++;
                 if(k==matrix.length) break;
                }
            }
      
            int columnIndex = starsByRow[primedZero[0]];
            if (-1 == columnIndex){

                incrementSetOfStarredZeroes(primedZero, starsByRow, starsByCol, primesByRow,RC,LL,infinity,matrix,n,L,BA);
                
               for(int f=0;f<matrix.length;f++)
                if(primesByRow[f]!=-1&&LL[primesByRow[f]]<n&&L[LL[primesByRow[f]]]>1&&BA[RC[f]]>1){
                    int t=primesByRow[f];
                    int s=f;
                    while(t>=0){
                    	if(LL[t]==LL[primesByRow[f]])
                    		t--;
                    	else break;
                    }
                    while(s>=0){
                    	if(RC[s]==RC[f])
                    		s--;
                    	else break;
                    }
                   for(int k=t+1;LL[k]<=LL[primesByRow[f]];)
                   {
                   if(LL[k]==LL[primesByRow[f]]&&k!=primesByRow[f])
                     for(int h=s+1;RC[h]<=RC[f];)
                     {
                      if(RC[h]==RC[f]&&h!=f&&infinity[h][k]==-1)
                         infinity[h][k]=0;
                      h++;
                      if(h==matrix.length) break;
                     }
                   k++;
                   if(k==matrix.length) break;
                   }
                }           
                      
                Arrays.fill(primesByRow,-1);
                Arrays.fill(coveredRows,0);
                Arrays.fill(coveredCols,0);
                coverColumnsOfStarredZeroes(starsByCol,coveredCols);
            } else {

                coveredRows[primedZero[0]] = 1;
                coveredCols[columnIndex] = 0;
            }
        }

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



    static private void reduceMatrix(float[][] matrix) {

        for (int i = 0; i < matrix.length; i++) {

            float minValInRow = Float.MAX_VALUE;
            for (int j = 0; j < matrix[i].length; j++) {
                if (minValInRow > matrix[i][j]) {
                    minValInRow = matrix[i][j];
                }
            }

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

    static private void initStars(float costMatrix[][], int[] starsByRow, int[] starsByCol,int [] RC,int [] LL,int[][] infinity,int n,int[] L,int[] BA) {


        int [] rowHasStarredZero = new int[costMatrix.length];
        int [] colHasStarredZero = new int[costMatrix[0].length];
        for (int i = 0; i < costMatrix.length; i++) {
            for (int j = 0; j < costMatrix[i].length; j++) {
                if (0 == costMatrix[i][j] && 0 == rowHasStarredZero[i] && 0 == colHasStarredZero[j]&&infinity[i][j]==0) {
                    starsByRow[i] = j;
                    starsByCol[j] = i;
                    rowHasStarredZero[i] = 1;
                    colHasStarredZero[j] = 1;
                    
                    if(LL[j]<n&&L[LL[j]]>1&&BA[RC[i]]>1){
                        int t=j;
                        int s=i;
                        while(t>=0){
                        	if(LL[t]==LL[j])
                        		t--;
                        	else break;
                        }
                        while(s>=0){
                        	if(RC[s]==RC[i])
                        		s--;
                        	else break;
                        }
                     for(int k=t+1;LL[k]<=LL[j];)
                      {
                      if(LL[k]==LL[j]&&k!=j)
                       for(int h=s+1;RC[h]<=RC[i];)
                       {
                    	   if(RC[h]==RC[i]&&h!=i)
                             infinity[h][k]=1;
                    	   h++;
                    	   if(h==costMatrix.length) break;
                       }
                       k++;
                       if(k==costMatrix.length) break;
                      }
                    }
                   
                    break; 
                }
                     
            }
        }
    }

    static private void coverColumnsOfStarredZeroes(int[] starsByCol, int[] coveredCols) {
        for (int i = 0; i < starsByCol.length; i++) {
            coveredCols[i] = -1 == starsByCol[i] ? 0 : 1;
        }
    }


    static  private int[] primeSomeUncoveredZero(float matrix[][], int[] primesByRow,
                                       int[] coveredRows, int[] coveredCols, int[] starsByRow,int [] RC,int [] LL,int[][] infinity) {


        for (int i = 0; i < matrix.length; i++) {
            if (1 == coveredRows[i]) continue;
            for (int j = 0; j < matrix[i].length; j++) {
                if (0 == matrix[i][j] && 0 == coveredCols[j]&&infinity[i][j]==0) {
                    primesByRow[i] = j;
 
                    return new int[]{i,j};
                }
            }
        }
        return null;

    }

    static  private void incrementSetOfStarredZeroes(int[] unpairedZeroPrime, int[] starsByRow, int[] starsByCol, int[] primesByRow,int []RC,int []LL,int [][] infinity,float[][] matrix,int n,int[] L,int[] BA) {

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

        for (int[] zero : zeroSequence) {

            if(starsByRow[zero[0]]==zero[1]){
            if(LL[zero[1]]<n&&L[LL[zero[1]]]>1&&BA[RC[zero[0]]]>1){
                int t=zero[1];
                int s=zero[0];
                while(t>=0){
                	if(LL[t]==LL[zero[1]])
                		t--;
                	else break;
                }
                while(s>=0){
                	if(RC[s]==RC[zero[0]])
                		s--;
                	else break;
                }
              for(int k=t+1;LL[k]<=LL[zero[1]];)
              { 
                if(LL[k]==LL[zero[1]]&&k!=zero[1])
                 for(int h=s+1;RC[h]<=RC[zero[0]];)
                 { 
                   if(RC[h]==RC[zero[0]]&&h!=zero[0])
                       infinity[h][k]=0;
                   h++;
                   if(h==matrix.length) break;
                 }
                 k++;
                 if(k==matrix.length) break;
            	}
              }
            }
            
            if (primesByRow[zero[0]] == zero[1]) {
                starsByRow[zero[0]] = zero[1];
                starsByCol[zero[1]] = zero[0];
               primesByRow[zero[0]]=-1;
             }
            
        }
        
         for (int[] zero : zeroSequence){
        	  if(starsByRow[zero[0]]==zero[1]){
        		if(LL[zero[1]]<n&&L[LL[zero[1]]]>1&&BA[RC[zero[0]]]>1){
                    int t=zero[1];
                    int s=zero[0];
                    while(t>=0){
                    	if(LL[t]==LL[zero[1]])
                    		t--;
                    	else break;
                    }
                    while(s>=0){
                    	if(RC[s]==RC[zero[0]])
                    		s--;
                    	else break;
                    }
                  for(int k=t+1;LL[k]<=LL[zero[1]];)
                  {
                      if(LL[k]==LL[zero[1]]&&k!=zero[1])
                       for(int h=s+1;RC[h]<=RC[zero[0]];)
                       {
                         if(RC[h]==RC[zero[0]]&&h!=zero[0])
                             infinity[h][k]=1;
                         h++;
                         if(h==matrix.length) break;
                       }
                      k++;
                      if(k==matrix.length) break;
                  }
        		}
        	 }
          }

    }


    static  private void makeMoreZeroes(float[][] matrix, int[] coveredRows, int[] coveredCols,int [][] infinity) {

        float minUncoveredValue = Float.MAX_VALUE;
        for (int i = 0; i < matrix.length; i++) {
            if (0 == coveredRows[i]) {
                for (int j = 0; j < matrix[i].length; j++) {
                    if (0 == coveredCols[j] && matrix[i][j] < minUncoveredValue&&infinity[i][j]==0) {
                        minUncoveredValue = matrix[i][j];
                    }
                }
            }
        }

        for (int i = 0; i < coveredRows.length; i++) {
            if (1 == coveredRows[i]) {
                for (int j = 0; j < matrix[i].length; j++) {
                        matrix[i][j] += minUncoveredValue;
                }
            }
        }

        for (int i = 0; i < coveredCols.length; i++) {
            if (0 == coveredCols[i]) {
                for (int j = 0; j < matrix.length; j++) {
                    matrix[j][i] -= minUncoveredValue;
                }
            }
        }
    }
	public static float RatedAssign(int [] L, int[] BA, float [][] Q, int [][] T, int m, int n) {
		float v=0.0f;
		int cnt=0, BAcnt=0;
		for (int i = 0; i<m; i++) {if(i<n) cnt+=L[i];BAcnt +=BA[i];}
		int LL[]=new int [BAcnt];
		int RC[]=new int[BAcnt];
		float Q1[][]= new float [BAcnt][n];
 
		if (cnt >BAcnt) return 0.0f;
		
		int[] arrL=Arrays.copyOf(L, n);
		int[] arrBA=Arrays.copyOf(BA, m);
		int len=arrBA.length;
		int s=0,zero=0;;
		Arrays.sort(arrBA);
		for(int i=arrL.length-1;i>=0;i--)
			if(arrL[i]<=len){
			  for(int j=arrBA.length-1;s<arrL[i];j--){
			        arrBA[j]--;
			        s++;
			        if(arrBA[j]==0)
			        	zero++;
			        }
			s=0;
			len=arrBA.length-zero;
			Arrays.sort(arrBA);
			}else
			  return 0.0f;
		
		int index=0;
		for (int j = 0; j<n; j++)
				for (int k = 0; k<L[j]; k++) LL[index++] =j;
		int t=n;
		for (int k = index; k < BAcnt; k++)
			{LL[k]=t;
			t++;}

		index=0;
		for (int j = 0; j<m; j++)
				for (int k = 0; k<BA[j]; k++) RC[index++] =j;
		
		for (int j = 0; j<n; j++)
		{ 	index =0;
			for (int i = 0; i<m; i++)
			{
				for (int k = 0; k<BA[i]; k++)
					Q1[index++][j]=Q[i][j];
			}
		}
		
		float [][] M = new float [BAcnt][BAcnt];
		for (int i = 0; i<BAcnt; i++)
		{ 	index =0;
			for (int j = 0; j<n; j++)
			{
				for (int k = 0; k<L[j]; k++)
					M[i][index++]=1-Q1[i][j];
			}
		}
               
		System.out.print("Matrix=");
		System.out.print(BAcnt);
		System.out.print("*");
		System.out.println(BAcnt);
		
		int [][] N = computeAssignments(M,RC,LL,n,L,BA);
				
		for (int i = 0; i < N.length; i++)
			if (LL[N[i][1]]< n)
			  T[RC[N[i][0]]][LL[N[i][1]]]=1;

		for (int i = 0; i<m; i++)
		 	for (int j = 0; j<n; j++)
				v += Q[i][j]*T[i][j];
		return v;
	}

	public static TestResult RandomTest(List<Double> ltime,List<Float> lmax) {
		Random generator = new Random();
		int m, n;

		m=5;n=5;

		int[] L=new int[n];
		int[] BA=new int[m];
		float[][] Q=new float[m][n];
		int T[][] =new int [m][n];
        int temp=0;
		       
		for (int i =0; i<m; i++){ 
			temp=generator.nextInt(6);
				if(temp==0) BA[i]=temp+1;
				else 
				BA[i]=temp;
				}
		for (int i =0; i<n; i++){
			temp=generator.nextInt(6);
				if(temp==0) L[i]=temp+1;
				else
				L[i]=temp;
			}
        
		for (int i = 0; i < m; i++) 		for (int j =0; j< n; j++)		
 		{ if (generator.nextInt(n)<n) {Q[i][j]=((float)generator.nextInt(100)+1) / 100; }}


		DecimalFormat tw = new DecimalFormat("0.000");

	System.out.println ("Before Transferring...");
	System.out.print ("m=");
	System.out.print (m);
	System.out.print (",");
	System.out.print ("n=");
	System.out.println (n);
	System.out.println ("The Role Range Vector:");
	for (int i =0; i<n; i++)
		System.out.print (L[i]+" ");
	System.out.println ("");
	System.out.println ("The Agent Range Vector:");
	for (int i =0; i<m; i++)
		System.out.print (BA[i]+" ");
	System.out.println ();
	System.out.println ("The Qualified Role Matrix:");
	for (int i = 0; i < m; i++)
	{	for (int j =0; j< n; j++)
		{
		System.out.print (tw.format(Q[i][j]));		System.out.print (" ");
		}
	System.out.println ();
	}
	
	TestResult tr = new TestResult();
	long t1 = System.nanoTime();
	double v=RatedAssign(L,BA,  Q, T, m, n);
	long t2 = System.nanoTime();
	double diff = (double)(t2-t1)/1000000;
	tr.time = diff;
	System.out.println ("Time = "+diff+"ms");
	
	
	System.out.println ();
	for (int i =0; i<n; i++)
		System.out.print (L[i]+" ");		System.out.println ();	
System.out.println ("The Assignment Matrix:");
float count = 0;
	for (int i = 0; i < m; i++)
		{	for (int j =0; j< n; j++)
			{	
			count += Q[i][j]*T[i][j];
			System.out.print (T[i][j]);	System.out.print (" ");	
			} System.out.println ();
			}
	System.out.print ("count = "+count);	System.out.println ();

		if (v>0) 
		{	System.out.println ("Success....\n After Assignment");
		System.out.println ("Max G-Q = "+v);
		System.out.println ();
		tr.TF=true;
		tr.count=count;
		}
		else 
			{
			System.out.println ("Fail...\n");
			tr.TF = false;
			}
				       
	return tr;
	}

	public static void main(String[] args) {

	List<Double> ltime=new ArrayList<Double>();
	List<Float>  lmax=new ArrayList<Float>();
	
	int num =100;
	TestResult []TR=new TestResult [num];
	int stop=num;
	for (int i = 0; i<num; i++)
	{	TR [i]=RandomTest(ltime,lmax);
	}
	for (int i = 0; i<stop; i++)
	{	if (TR[i].TF) System.out.print (1);
	else System.out.print (0);
	System.out.println (" "+TR[i].time);
	}
  }
}


