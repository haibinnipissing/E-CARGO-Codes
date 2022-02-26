
/**
 * Using Simulated Annealing and GRA to solve the Agent Classification Problem
 * Please cite: 
 [1] H. Zhu, “Agent Categorization with Group Role Assignment with Constraints (GRA+) and Simulated Annealing (SA),” IEEE Trans. on Computational Social Systems, vol. 7, no. 5, Oct. 2020, pp. 1234-1245.  
 [2] H. Zhu, E-CARGO and Role-Based Collaboration: Modeling and Solving Problems in the Complex World, Wiley-IEEE Press, NJ, USA, Dec. 2021. 
 [3] H. Zhu, M.C. Zhou, and R. Alkins, “Group Role Assignment via a Kuhn-Munkres Algorithm-based Solution”, IEEE Trans. on Systems, Man, and Cybernetics, Part A: Systems and Humans, vol. 42, no. 3, May 2012, pp. 739-750.
 [4] H. Zhu, and M. Zhou, “Role-Based Collaboration and its Kernel Mechanisms,” IEEE Trans. on Systems, Man, and Cybernetics, Part C: Applications and Reviews, vol. 36, no. 4, July. 2006, pp. 578-589.
 * An implementation of Agent Categorization with Simulated Annealing Algorithm. 
 * @author Haibin Zhu, Jan. 13, 2020
*/

import java.text.DecimalFormat;
import java.util.Random;
import java.util.ArrayList;
import java.util.Collections;
import ilog.concert.*;
import ilog.cplex.*;

class GRA_ILOG {
	private int m;	//number of agents
	private int n;	//number of roles
	private double[] Q;	//Qualification matrix
	private int[] L;	//Requirement array
	private int[][] A;  //Assignment array
	
	DecimalFormat df = new DecimalFormat("0.00");
	
	double optimized_result = 0;
	boolean bILOG_result;
	
	public GRA_ILOG(int nagent, int nrole, double[][] QM, int[]RA)
	{
		m = nagent;
		n = nrole;
		
		Q = new double[m*n];
		for(int i=0, r=0; r<m; r++) for (int c=0; c<n; c++){Q[i] = QM[r][c]; i++; }
		
		L = new int[n];
		L = RA;
		
		A = new int[m][n];
		for(int r=0; r<m; r++) for (int c=0; c<n; c++) A[r][c] = 0;
		
		//LOG:
		System.out.println("Qualification Matrix: ");
		for (int i=0;i<m*n;i++)
		{
			System.out.print(df.format(Q[i])+"	");
			if ((i+1)%(n) == 0) System.out.print("\n");
			
		}
		System.out.print("\n");	
		
		System.out.println("Requirement Array: ");
		for(int i=0; i<n; i++)
		{
			System.out.print(L[i]+"	");
		}
		System.out.print("\n");	
	}
	
	public double resolve(int[][]TR)
	{
		try
		{
			//Creat cplex obj
			IloCplex cplex = new IloCplex();	//initialize the cplex object
			
			IloIntVar[]x = cplex.intVarArray(m*n, 0, 1);	//initialize the variables array under cplex.
			
			//cplex.addMinimize(cplex.scalProd(x, Q));	//add the optimize objective to cplex.
			cplex.addMaximize(cplex.scalProd(x, Q));	//add the optimize objective to cplex.
			
			//Add Constrains:
			
			//Constrain type 1: unique constrains here, one person can only be assigned on one role at one time, 
			//thus there are number of 'm' constrains here need to be inserted into the cplex obj.
			for(int i=0; i<m; i++)
			{
				IloLinearNumExpr exprUniConstrain = cplex.linearNumExpr();
				for(int j = 0; j<n; j++)
				{
					exprUniConstrain.addTerm(1, x[n*i+j]);
				}
				cplex.addLe(exprUniConstrain, 1.0);
				
			}
			//Constrain type 2: Add role requirement constrains, 
			//the number of people assigned on each role should meet the requirement on that role.
			//Hence, n constrains will be added.
			
			for (int i = 0; i<n; i++)
			{
				IloLinearNumExpr exprReqConstrain = cplex.linearNumExpr();
				for (int j = 0; j<m; j++)
				{
				exprReqConstrain.addTerm(1, x[i+j*n]);
				}
				cplex.addGe(exprReqConstrain, L[i]);				
			}
			
			//Solve LP
			//long t1 = System.nanoTime();
			if (cplex.solve()) 
			{
				bILOG_result = true;
				optimized_result = cplex.getObjValue();
				
				double[] val = cplex.getValues(x);
				int ncols = cplex.getNcols();
				//cplex.output().println("Num COL: " + ncols);
				
				cplex.output().println("Result Table: " );
				for (int j=0; j<ncols; j++)
				{
					A[j/n][j%n] = (int)val[j]; 
					System.out.print(A[j/n][j%n] + " ");
					TR[j/n][j%n] = A[j/n][j%n];
					//System.out.print(val[j]+ "	");
					if ((j+1)%(n) == 0) {System.out.print("\n");}	
				}
				cplex.end();
			}
			else
			{
				cplex.end();
				bILOG_result = true;
			}	
		}
		catch (IloException e){System.err.println("Concert exception" + e + " caught");}
		
		
		return(optimized_result);
	}
	
	public double getOptimizedResult()
	{
		return optimized_result;
		
	}
}


class Agent {
    private int id;
    private String AgentName;            

    //Constructor
    //creates a Agent given its name and (id,y) location
	public Agent(String AgentName, int id){
    	this.AgentName = AgentName;
        this.id = id;
    }            
        
    /**
	 * @return the id
	 */
	public int getid() {
		return id;
	}
	/**
	 * @param id the id to set
	 */
	public void setid(int id) {
		this.id = id;
	}

	/**
	 * @return the AgentName
	 */
	public String getAgentName() {
		return AgentName;
	}

	/**
	 * @param AgentName the AgentName to set
	 */
	public void setAgentName(String AgentName) {
		this.AgentName = AgentName;
	}
	
}
class Roles{
	static int m=24, n =6;
	double Q [][];
	int T[][];
	static double Sim [][]={{100,30,50,60,40,60,30,30,30,50,60,30,60,50,40,20,70,60,60,50,70,60,60,60},
        	{30,100,30,30,50,20,80,80,80,10,20,70,20,10,40,40,20,40,30,40,30,20,20,20},
        	{50,30,100,60,20,40,30,30,30,40,50,40,50,50,30,20,40,40,40,30,60,80,80,80},
        	{60,30,60,100,0,90,20,20,20,30,40,30,80,60,50,0,50,50,50,40,60,60,60,60},
        	{40,50,20,0,100,10,30,20,30,20,10,40,0,0,40,90,0,20,20,60,20,40,50,60},
        	{60,20,40,90,10,100,20,20,0,40,40,30,50,60,0,10,50,40,60,50,80,30,30,30},
        	{30,80,30,20,30,20,100,90,90,20,20,90,0,0,40,30,0,20,0,0,0,0,0,0},
        	{30,80,30,20,20,20,90,100,90,10,10,85,0,0,40,30,0,20,0,0,0,0,0,0},
        	{30,80,30,20,30,0,90,90,100,10,10,90,0,0,40,30,0,20,0,0,0,0,0,0},
        	{50,10,40,30,20,40,20,10,10,100,40,0,30,30,0,0,20,30,30,30,40,30,30,30},
        	{60,20,50,40,10,40,20,10,10,40,100,0,40,20,0,0,30,40,30,30,40,85,85,85},
        	{30,70,40,30,40,30,90,85,90,0,0,100,0,0,40,30,0,20,0,0,0,30,30,30},
        	{60,20,50,80,0,50,0,0,0,30,40,0,100,40,0,0,50,40,50,20,60,50,50,50},
        	{50,10,50,60,0,60,0,0,0,30,20,0,40,100,0,0,40,40,10,0,50,40,40,40},
        	{40,40,30,50,40,0,40,40,40,0,0,40,0,0,100,20,30,20,40,40,0,40,40,40},
        	{20,40,20,0,90,10,30,30,30,0,0,30,0,0,20,100,0,30,20,30,0,30,30,30},
        	{70,20,40,50,0,50,0,0,0,20,30,0,50,40,30,0,100,30,20,40,50,40,40,40},
        	{60,40,40,50,20,40,20,20,20,30,40,20,40,40,20,30,30,100,50,40,50,50,50,50},
        	{60,30,40,50,20,60,0,0,0,30,30,0,50,10,40,20,20,50,100,40,40,70,70,70},
        	{50,40,30,40,60,50,0,0,0,30,30,0,20,0,40,30,40,40,40,100,40,60,60,60},
        	{70,30,60,60,20,80,0,0,0,40,40,0,60,50,0,0,50,50,40,40,100,50,50,50},
        	{60,20,80,60,40,30,0,0,0,30,85,30,50,40,40,30,40,50,70,60,50,100,90,90},
        	{60,20,80,60,50,30,0,0,0,30,85,30,50,40,40,30,40,50,70,60,50,90,100,90},
        	{60,20,80,60,60,30,0,0,0,30,85,30,50,40,40,30,40,50,70,60,50,90,90,100}};
    //to hold a Roles of cities
    private ArrayList<Agent> Roles = new ArrayList<Agent>();
    
    //we assume initial value of Similarity is 0 
    private int Similarity = 0;
    
    //Constructor
    //starts an empty Roles
    public Roles(){
        for (int i = 0; i < n; i++) {
            Roles.add(null);
        }
    }
    
    //another Constructor
    //starts a Roles from another Roles
    @SuppressWarnings("unchecked")
	public Roles(ArrayList<Agent> Roles){
        this.Roles = (ArrayList<Agent>) Roles.clone();
    }
    
    /**
      Returns Roles information
      @return currentRoles
     */
    public ArrayList<Agent> getRoles(){
        return Roles;
    }
     
    /**
     * Creates a random Roles (i.e. individual or candidate solution)
     */
    public void generateIndividual() {
        // Loop through all our destination cities and add them to our Roles
        for (int AgentIndex = 0; AgentIndex < n; AgentIndex++) {
          setAgent(AgentIndex, RolesManager.getAgent(AgentIndex));
        }
        // Randomly reorder the Roles
        Collections.shuffle(Roles);
    }

    /**
     * Returns a Agent from the Roles given the Agent's index
     * @param index
     * @return Agent at that index
     */
    public Agent getAgent(int index) {
        return Roles.get(index);
    }

    /**
     * Sets a Agent in a certain position within a Roles
     * @param index
     * @param Agent
     */
    public void setAgent(int index, Agent Agent) {
        Roles.set(index, Agent);
        // If the Roles has been altered we need to reset the fitness and Similarity
        Similarity = 0;
    }
    
    /**
     * Computes and returns the total Similarity of the Roles
     * @return Similarity total Similarity of the Roles
     */
    public double getTotalSimilarity(){

    	T=new int [m][n];
    	Q=new double [m][n];
    	for (int i =0; i < m; i++)
        	for (int j =0; j < n; j++)
        		Q[i][j]=Sim[i][Roles.get(j).getid()];
    	int L[] = {3,3,3,3,3,3};
  //  	int L[] = {2,2,2,2,2,2};
		GRA_ILOG ILOG = new GRA_ILOG(m, n, Q, L);
		double v1 = ILOG.resolve(T);//ILOG.resolve(TR, time);
        return v1;
    }

    /**
     * Get number of cities on our Roles
     * @return number how many cities there are in the Roles!
     */
    public int RolesSize() {
        return Roles.size();
    }
    
    @Override
    /**
     * To print out a list of all the cities in the Roles
     */
    public String toString() {
        String s = getAgent(0).getAgentName();
        for (int i = 1; i < RolesSize(); i++) {
            s += ", " + getAgent(i).getAgentName();
        }
        return s;
    }
}

class RolesManager {

    // Holds our cities
    private static ArrayList<Agent> destinationCities = new ArrayList<Agent>();

    /**
	 * Adds a destination Agent
	 * @param Agent
	 */
	public static void addAgent(Agent Agent) {
		destinationCities.add(Agent);
	}

	/**
	 * returns a Agent given its index
	 * @param index
	 * @return Agent the Agent at index 
	 */
	public static Agent getAgent(int index){
		return (Agent)destinationCities.get(index);
	}

	/**
	 * Returns the number of destination cities 
	 * @return size the number of destination cities
	 */
	public static int numberOfCities(){
		return destinationCities.size();
	}
    
}

class Utility {


	/**
	 * Computes and returns the Euclidean Similarity between two cities
	 * @param Agent1 the first Agent
	 * @param Agent2 the second Agent
	 * @return Similarity the dist between Agent1 and Agent2
	 */
	/**
	 * Calculates the acceptance probability
	 * @param currentSimilarity the total Similarity of the current Roles
	 * @param newSimilarity the total Similarity of the new Roles
	 * @param temperature the current temperature
	 * @return value the probability of whether to accept the new Roles
	 */
	public static double acceptanceProbability(double currentSimilarity, double newSimilarity, double temperature) {
		// If the new solution is better, accept it
		if (newSimilarity > currentSimilarity) {
			return 1.0;
		}
		// If the new solution is worse, calculate an acceptance probability
		return Math.exp((currentSimilarity - newSimilarity) / temperature);
	}

	/**
	 * this method returns a random number n such that
	 * 0.0 <= n <= 1.0
	 * @return random such that 0.0 <= random <= 1.0
	 */
	static double randomDouble()
	{
		Random r = new Random();
		return r.nextInt(1000) / 1000.0;
	}
	
	/**
	 * returns a random int value within a given range
	 * min inclusive .. max not inclusive
	 * @param min the minimum value of the required range (int)
	 * @param max the maximum value of the required range (int)
	 * @return rand a random int value between min and max [min,max)
	 */ 
	public static int randomInt(int min , int max) {
		Random r = new Random();
		double d = min + r.nextDouble() * (max - min);
		return (int)d;
	}
}

public class GitHub {
    
    public static void main(String[] args) {
        // Create and add our cities
    	int m =24, n =6;
        Agent Agent = new Agent("Anthropology",0);
        RolesManager.addAgent(Agent);
        Agent Agent2 = new Agent("Biology", 1);
        RolesManager.addAgent(Agent2);
        Agent Agent3 = new Agent("Child and Family Studies",2);
        RolesManager.addAgent(Agent3);
        Agent Agent4 = new Agent("Classical Studies",3);
        RolesManager.addAgent(Agent4);
        Agent Agent5 = new Agent("Computer Science",4 );
        RolesManager.addAgent(Agent5);
        Agent Agent6 = new Agent("English Studies",5 );
        RolesManager.addAgent(Agent6);
        Agent Agent7 = new Agent("Environmental and Physical Geography",6);
        RolesManager.addAgent(Agent7);
        Agent Agent8 = new Agent("Environmental Biology and Technology",7);
        RolesManager.addAgent(Agent8);
        Agent Agent9 = new Agent("Environmental Geography",8);
        RolesManager.addAgent(Agent9);
        Agent Agent10 = new Agent("Fine Arts",9);
        RolesManager.addAgent(Agent10);        
        Agent Agent11 = new Agent("Gender Eqality and Social Justice",10);
        RolesManager.addAgent(Agent11);
        Agent Agent12 = new Agent("Geography",11);
        RolesManager.addAgent(Agent12);
        Agent Agent13 = new Agent("History",12);
        RolesManager.addAgent(Agent13);
        Agent Agent14 = new Agent("Liberal Arts",13);
        RolesManager.addAgent(Agent14);
        Agent Agent15 = new Agent("Liberal Science",14);
        RolesManager.addAgent(Agent15);
        Agent Agent16 = new Agent("Mathematics",15);
        RolesManager.addAgent(Agent16);
        Agent Agent17 = new Agent("Native Studies",16);
        RolesManager.addAgent(Agent17);
        Agent Agent18 = new Agent("Philosophy",17);
        RolesManager.addAgent(Agent18);
        Agent Agent19 = new Agent("Political Science",18);
        RolesManager.addAgent(Agent19);
        Agent Agent20 = new Agent("Psychology",19);
        RolesManager.addAgent(Agent20);
        Agent Agent21 = new Agent("Religions and Cultures",20);
        RolesManager.addAgent(Agent21);
        Agent Agent22 = new Agent("Social Welfare and Social Develpment",21);
        RolesManager.addAgent(Agent22);
        Agent Agent23 = new Agent("Social Work",22);
        RolesManager.addAgent(Agent23);
        Agent Agent24 = new Agent("Social Work",23);
        RolesManager.addAgent(Agent24);

        //Set initial temp
        double temp = 100000;

        //Cooling rate
 //       double coolingRate = 0.003;
        double coolingRate = 0.005;
    //      double coolingRate = 0.007;
  //          double coolingRate = 0.01;
  //      double coolingRate = 0.015; 
     
   //     double coolingRate = 0.02; 
//        double coolingRate = 0.03; //(No!)
        
 //       double coolingRate = 0.05; //(No!)
           
        //create random intial solution
        Roles currentSolution = new Roles();
        currentSolution.generateIndividual();
        
        System.out.println("Total Similarity of initial solution: " + currentSolution.getTotalSimilarity());
        System.out.println("Roles: " + currentSolution);
		long t1 = System.nanoTime();
        // We would like to keep track if the best solution
        // Assume best solution is the current solution
        Roles best = new Roles(currentSolution.getRoles());
        int x=0, y=0;
        // Loop until system has cooled
        while (temp > 1) {
            // Create new neighbour Roles
            Roles newSolution = new Roles(currentSolution.getRoles());
            ArrayList<Agent> NewlyUsedCities = new ArrayList<Agent>();
            for (int i =0; i<n; i++) NewlyUsedCities.add(currentSolution.getAgent(i));

            // Get random positions in the Roles
            int RolesPos1 = Utility.randomInt(0 , n);
          //  int RolesPos2 = Utility.randomInt(n , newSolution.m);
            int RolesPos2;
            RolesPos2 = (n+x++)%m; 
            //to make sure that RolesPos1 and RolesPos2 are different
//    		while(newSolution.getRoles().get(RolesPos1).getX() == newSolution.getRoles().get(RolesPos2).getX()) {RolesPos2 = Utility.randomInt(0 , newSolution.m);}
    		while(NewlyUsedCities.contains(RolesManager.getAgent(RolesPos2))) 
    		{
    			RolesPos2 = (n+x++)%m;
    		}

            // Get the cities at selected positions in the Roles
            Agent AgentSwap2 = RolesManager.getAgent(RolesPos2);
            // Change a Agent
            newSolution.setAgent(RolesPos1, AgentSwap2);
            // Get energy of solutions
            double currentSimilarity   = currentSolution.getTotalSimilarity();
            double neighbourSimilarity = newSolution.getTotalSimilarity();

            // Decide if we should accept the neighbour
            double rand = Utility.randomDouble();
            if (Utility.acceptanceProbability(currentSimilarity, neighbourSimilarity, temp) > rand) {
                currentSolution = new Roles(newSolution.getRoles());
            }

            // Keep track of the best solution found
            if (currentSolution.getTotalSimilarity() > best.getTotalSimilarity()) {
                best = new Roles(currentSolution.getRoles());
            }
            
            // Cool system
            temp *= 1 - coolingRate;
            y++;
        }
		long t2 = System.nanoTime();

        System.out.println("Maximized Assinged Similarity Score: " + best.getTotalSimilarity());
        double total=0;
        for (int j =0; j< m; j++) 
        	for (int i =0; i< m; i++)
        		total+=Roles.Sim[j][i];
    	System.out.println("Total values: " + total +".");
		System.out.println("Iteration: " + y +".");
		double diff = (double)(t2-t1)/1000000;
		System.out.println("Time: " + diff +"ms");
        System.out.println("Roles: " + best);
        for (int j =0; j< n; j++) {
        	for (int i =0; i< m; i++)
        		if (1==best.T[i][j]) System.out.print(RolesManager.getAgent(i).getAgentName()+",");
            System.out.print(";");
        }
    }
}
