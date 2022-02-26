/*
 * This short program is to solve a simplified GMEO problem. 
 * Created by Haibin Zhu, May 30, 2017
 * Please cite:
 [1] H. Zhu, “Avoiding Critical Members in a Team by Redundant Assignment,” IEEE Trans. on Systems, Man, and Cybernetics: Systems, vol. 50, no. 7, July 2020, pp. 2729-2740.
 [2] H. Zhu, E-CARGO and Role-Based Collaboration: Modeling and Solving Problems in the Complex World, Wiley-IEEE Press, NJ, USA, Dec. 2021. 
 [3] H. Zhu, M.C. Zhou, and R. Alkins, “Group Role Assignment via a Kuhn-Munkres Algorithm-based Solution”, IEEE Trans. on Systems, Man, and Cybernetics, Part A: Systems and Humans, vol. 42, no. 3, May 2012, pp. 739-750.
 [4] H. Zhu, and M. Zhou, “Role-Based Collaboration and its Kernel Mechanisms,” IEEE Trans. on Systems, Man, and Cybernetics, Part C: Applications and Reviews, vol. 36, no. 4, July. 2006, pp. 578-589.

*/


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.Random;
class TestResult
{
	public long GP1;
	public int na, nt, nd;
	public TestResult()
	{
		GP1=0;//Group Performance GRA
		na= 0;
		nd = 0;
		nt =0;
	}
};

public class GMEO_SS {
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
	private static String filename="Result.txt";
	public static void main(String[] args)
	{	
		Random generator = new Random();
		DecimalFormat df = new DecimalFormat("0.00");
		
		int m = 30, n = 6;
		int T[][]= new int [m-1][n];
//		try 
//		{BufferedWriter out = new BufferedWriter(new FileWriter(("GMEO.txt"), true));	
		m=600;n = 50;
		int L[]=new int [n];
		TestResult []testR=new TestResult[100];
		for (int kk=0; kk< 100; kk++){
			testR[kk]=new TestResult();
//		double Q[][]=new double [m][n];
		int na = 0;
		for (int j =0; j<n; j++) { L[j]=generator.nextInt(m/n)+1; na+=L[j];};	
		 int QD[][]=new int [m][n];//Initialize all the elements of   with 0s.
		long t1 = System.nanoTime();
		int st = 0, nt=0;
		for (int i =0; i<m; i++) for (int j =0; j<n; j++) QD[i][j]=0;
	 	for (int j =0; j<n; j++){	
			for (int k =st; k<=st+L[j]; k++){ 
				 QD[k][j]=1;
				 nt++;
			}
			 st=st+L[j];
		}		
		long t2 = System.nanoTime();
		long diff = (t2-t1)/1000000;
		testR[kk].GP1=diff;
		testR[kk].na=na;
		testR[kk].nt=nt;
		int actualm=0;
		for(int r=0; r<m; r++)	for(int c=0; c<n; c++)	{if (QD[r][c]==1) {actualm++; break;}}
		testR[kk].nd=actualm;
		}
		System.out.println("na	nt	nd	Time (ms)");
		for (int kk=0; kk< 100; kk++)
			System.out.println(""+testR[kk].na+" "+testR[kk].nt+" "+testR[kk].nd+" "+testR[kk].GP1);
		
/*		System.out.println("na= "+na+" nt= "+nt);
		int actualm=0;
		for(int r=0; r<m; r++)	for(int c=0; c<n; c++)	{if (QD[r][c]==1) {actualm++; break;}}
		System.out.println("real_m= "+actualm);
		System.out.print("L=[");
		for(int c=0; c<n; c++) System.out.print(L[c]+"  ");
		System.out.print("]\n");
*/		
/*		int QDD[][]=new int [m-1][n];
		int i1 =6, j1=2;
		for (int i =0; i<i1; i++) for (int j =0; j<n; j++) QDD[i][j]=QD[i][j];
		for (int i =i1; i<m-1; i++) for (int j =0; j<n; j++) QDD[i][j]=QD[i+1][j];		
		for (int i =0; i<m-1; i++) for (int j =0; j<n; j++) T[i][j]=0;
		 st=0;
		for (int j =0; j<j1; j++){
			st += L[j];
			for (int i =0; i<i1; i++){				
				if (i!=st) T[i][j]=QDD[i][j];
			}
		}
		for (int j =j1; j<n; j++) for (int i =0; i<i1; i++) T[i][j]=QDD[i][j];
		for (int j =0; j<=j1; j++) for (int i =i1; i<m-1; i++) T[i][j]=QDD[i][j];
		for (int j =j1+1; j<n; j++){
			st += L[j-1];
			for (int i =i1; i<m-1; i++){				
				if (i!=st-1) T[i][j]=QDD[i][j];
			}
		}
		printIMatrix(T, m-1, n);
*/		//LOG result:		
	}
}

