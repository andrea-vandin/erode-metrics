package it.imt.erode.metrics.algorithm.affine;

import java.io.BufferedWriter;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.eclipse.ui.console.MessageConsoleStream;

import com.google.ortools.Loader;
import com.google.ortools.linearsolver.MPConstraint;
import com.google.ortools.linearsolver.MPObjective;
import com.google.ortools.linearsolver.MPSolver;
import com.google.ortools.linearsolver.MPVariable;

import it.imt.erode.commandline.CRNReducerCommandLine;
import it.imt.erode.crn.interfaces.ICRN;
import it.imt.erode.crn.interfaces.ICRNReaction;
import it.imt.erode.crn.interfaces.ISpecies;
import it.imt.erode.metrics.data.Matrix;
import it.imt.erode.metrics.data.SolutionAndOptimalValue;
import it.imt.erode.onthefly.Pair;
import it.imt.erode.partition.interfaces.IBlock;
import it.imt.erode.partition.interfaces.IPartition;
import it.sssa.erode.or_tools.BasicExample;

public class SimplePolicyIteration {
	//TODO: in computing all upper-bounds, I should use bigdecimals...

	private static final double INFINITY = java.lang.Double.POSITIVE_INFINITY;
	private static final double NEGATIVE_INFINITY = java.lang.Double.NEGATIVE_INFINITY;
	protected static final double TOLERANCE = 1e-6;
	
	private static final HashMap<ISpecies, BigDecimal> EMPTYROW=new LinkedHashMap<>(0);
	
	private MessageConsoleStream out;
	private BufferedWriter bwOut;
	private ICRN affineCRN;
	//private LinkedHashMap<ISpecies, ArrayList<ICRNReaction>> Arows;
	private LinkedHashMap<ISpecies, HashMap<ISpecies,BigDecimal>> Arows;
	private BigDecimal[] bentries;
	private static final String I="I"; 

	public ICRN getAffineCRN() {
		return affineCRN;
	}

	public static void main(String[] args) {
		SimplePolicyIteration.basicExampleOrTools(null,null);
	}

	public SimplePolicyIteration(MessageConsoleStream out, BufferedWriter bwOut, ICRN crn) {
		super();
		this.out = out;
		this.bwOut = bwOut;
		this.affineCRN=crn;

		//Drop unsupported models
		if(!(affineCRN.isElementary()&&isUnaryReag_BinaryProd(affineCRN))) {
			CRNReducerCommandLine.print(out,bwOut,"Model not supported. I return\n");
			throw new UnsupportedOperationException("Unsupported model");
		}

		//Init B
		bentries=new BigDecimal[affineCRN.getSpeciesSize()];
		for(int s=0;s<affineCRN.getSpeciesSize();s++){
			bentries[s]=BigDecimal.ZERO;
		}
		
		//Compute back the row A_i for every species i, and Bi for every species i
		Arows = new LinkedHashMap<>(crn.getSpeciesSize());
		for(ICRNReaction r : crn.getReactions()) {
			Pair ij = getIJ(r);
			if(isBentry(ij)) {
				addToBEntry((ISpecies)ij.getFirst(),r);
			}
			else {
				addToRowOfA((ISpecies)ij.getFirst(),r);
			}
			
		}
				
	}

	private boolean isBentry(Pair ij) {
		ISpecies spi = (ISpecies)ij.getSecond();
		if(spi.getName().equals(I)) {
			return true;
		}
		else {
			return false;
		}
		
	}

	/**
	 * We assume to have a matrix A_{(i,j)}. This method adds the contribution of the reaction r to the row of i (Ai)
	 * 
	 * An entry b_{row}     represents "\dot{x_{row}} = b_{row}" =CRN=> I -{b_{row}}-> I + x_{row}
	 * An entry A_{row,col} represents "\dot{x_{row}} = A_{row,col} x_{col}" =CRN=> x_{col} -{A_{row,col}}-> x_{col} + x_{row}
	 * @param i
	 * @param r
	 */
	private void addToRowOfA(ISpecies i, ICRNReaction r) {
		HashMap<ISpecies, BigDecimal> Ai = Arows.get(i);
		if(Ai == null) {
			Ai = new LinkedHashMap<>();
			Arows.put(i, Ai);
		}
		
		Pair ij = getIJ(r);
		ISpecies j = (ISpecies)ij.getSecond();
		BigDecimal Aij = Ai.get(j);
		if(Aij==null) {
			Aij=r.getRate();
		}
		else {
			Aij=Aij.add(r.getRate());
		}
		Ai.put(j, Aij);
//		ArrayList<ICRNReaction> Ai = Arows.get(i);
//		if(Ai == null) {
//			Ai = new ArrayList<>();
//			Arows.put(i, Ai);
//		}
//		Ai.add(r);
	}
	private void addToBEntry(ISpecies i, ICRNReaction r) {
		BigDecimal current = bentries[i.getID()];
		if(current==null) {
			bentries[i.getID()]=r.getRate();
		}
		else {
			bentries[i.getID()]=current.add(r.getRate());
		}
	}

	protected HashMap<ISpecies, BigDecimal> getRowOfA(ISpecies i){
		HashMap<ISpecies, BigDecimal> Ai = Arows.get(i);
		if(Ai!=null) {
			return Ai;
		}
		else {
			//we return an empty list - WE DO NOT ADD IT TO A! It shall not be modified
			return new LinkedHashMap<>(0);
		}
	}

	/**
	 * An entry b_{row}     represents "\dot{x_{row}} = b_{row}"             =CRN=> I       -{b_{row}    }-> I       + x_{row}
	 * An entry A_{row,col} represents "\dot{x_{row}} = A_{row,col} x_{col}" =CRN=> x_{col} -{A_{row,col}}-> x_{col} + x_{row}
	 * Therefore, 
	 * - i is the 'other' product (the product that is not also a reagent
	 * - j is the reagent (I for b entries)
	 *  
	 * @param r
	 * @return
	 */
	private static Pair getIJ(ICRNReaction r) {
		ISpecies j = r.getReagents().getFirstReagent();
		ISpecies i = r.getProducts().getFirstReagent();
		if(i.equals(j)) {
			i = r.getProducts().getSecondReagent();
		}
		return new Pair(i,j);
	}

	public Matrix computeMetrics(IPartition R, double lambda) {

		if(lambda<=0) {
			CRNReducerCommandLine.print(out,bwOut,"Lambda must be positive. I return\n");
			return null;
		}
		
		Loader.loadNativeLibraries();

		Matrix Dr = buildCostMatrixOfR(R);

		HashMap<Pair,Matrix> pi = constructInitialPolicy(R,lambda,Dr);
		
		
		//Compute delta_lambda^pi
		SolutionAndOptimalValue solAndMatrix = computeInvariantImprovement(lambda,pi);
		Matrix current_invariant= solAndMatrix.getSolution();
		
		//while exists (i,j) s.t. ...
		current_invariant = iterativeStep(lambda, pi, current_invariant);
		

		Matrix D = current_invariant;
		return D;
	}

	protected void println(String msg) {
		CRNReducerCommandLine.println(out, bwOut,msg);
	}
	
	public Matrix iterativeStep(double lambda, HashMap<Pair, Matrix> pi, Matrix current_invariant) {
		boolean existsAPair = false;
		int pairsUpdated=0;
		int iteration=1;
		do {
			CRNReducerCommandLine.println(out, bwOut,"Iteration "+iteration+"...");
			existsAPair = false;
			pairsUpdated=0;
			for(int i=0;i<affineCRN.getSpeciesSize();i++){
				ISpecies spi=affineCRN.getSpecies().get(i);
				for(int j=0;j<affineCRN.getSpeciesSize();j++){
					ISpecies spj=affineCRN.getSpecies().get(j);

					SolutionAndOptimalValue solutionAndValue = computeTransportScheduleMap(lambda,current_invariant,0,getRowOfA(spi),getRowOfA(spj));
					double Delta_current_invariant = solutionAndValue.getValue() + absbiminusbj(i, j);
					double prev =getOrDefaultValueIfHMNull(current_invariant, i, j,0);
					if(Delta_current_invariant < prev && prev - Delta_current_invariant > TOLERANCE )
					{
						//System.out.println("      delta="+Delta_current_invariant);
						//System.out.println("       prev="+prev);
						pairsUpdated++;
						existsAPair=true;
						Pair ij = new Pair(spi, spj);
						//Update the policy at (i,j)
						Matrix pi_ij = solutionAndValue.getSolution(); //computeTransportScheduleMap(lambda,current_policy,getRowOfA(spi),getRowOfA(spj));
						pi.put(ij, pi_ij);

						//Prepare for next iteration by updating the invariat
						SolutionAndOptimalValue solAndM = computeInvariantImprovement(lambda,pi);
						current_invariant= solAndM.getSolution();
						
						//double prev_new =getOrDefaultValueIfHMNull(current_invariant, i, j,0);
						//CRNReducerCommandLine.println(out, bwOut,"  Value for "+ij+" changed from "+prev+" to "+prev_new);
					}
				}
			}
			CRNReducerCommandLine.println(out, bwOut,"  Considered "+pairsUpdated+" pairs");
			iteration++;
		}while(existsAPair);
		return current_invariant;
	}

	/**
	 * Builds the cost matrix D_R of the equivalence/partition R
	 * - 0 if ij belong to R
	 * - 1 otherwise
	 * @param R
	 * @return
	 */
	private Matrix buildCostMatrixOfR(IPartition R) {
		Matrix Dr = new Matrix(1.0,affineCRN.getSpecies());
		IBlock current = R.getFirstBlock();
		while(current!=null) {
			for(ISpecies i : current.getSpecies()) {
				for(ISpecies j : current.getSpecies()) {
					Dr.setValue(i,j,0.0);
				}
			}
			current=current.getNext();
		}
		return Dr;
	}


	protected HashMap<Pair, Matrix> constructInitialPolicy(IPartition R,double lambda, Matrix Dr) {
		HashMap<Pair, Matrix> pi = new LinkedHashMap<Pair, Matrix>();

		//We explicitly compute the matrix Klambda(Dr)(Ai,Aj) for every pair in R, 
		//We do not explicitly compute the matrix 0 for pair not in R
		IBlock current = R.getFirstBlock();
		while(current!=null) {
			for(ISpecies i : current.getSpecies()) {
				for(ISpecies j : current.getSpecies()) {
					Pair ij = new Pair(i, j);
					HashMap<ISpecies, BigDecimal> Ai = getRowOfA(i);
					HashMap<ISpecies, BigDecimal> Aj = getRowOfA(j);
					SolutionAndOptimalValue solAndVal = computeTransportScheduleMap(lambda,Dr,1,Ai,Aj);
					Matrix pi_ij = solAndVal.getSolution();
					pi.put(ij, pi_ij);
				}
			}
			current=current.getNext();
		}

		return pi;
	}

	/**
	 * Here we should use or-tools to compute K(D)(c,d) tosolve a transport schedule problem
	 * We shall use a solver that provides solutions that are 'vertex'
	 * 	- This is enforced, e.g., by any solver using the simplex method 
	 * @param lambda
	 * @param D
	 * @param c
	 * @param d
	 * @return
	 */
	protected SolutionAndOptimalValue computeTransportScheduleMap(double lambda, Matrix D, double defValD, HashMap<ISpecies,BigDecimal> c, HashMap<ISpecies,BigDecimal> d) {
		//TODO: are we getting vertex solutions!? 
		// GLOP should provide vertex solutions. Not sure about SCIP 

		//The following code creates the data for the problem.
		int numWorkers = affineCRN.getSpeciesSize();
		int numTasks = affineCRN.getSpeciesSize();
		//The costs are contained in D


		// Create the linear solver with the solver backend.
		//String backend="SCIP";
		String backend="GLOP";
		MPSolver solver = MPSolver.createSolver(backend);
		if (solver == null) {
			System.out.println("Could not create solver "+backend);
			return null;
		}


		//Create the variables
		//w[i][j] the matrix to be returned 
		MPVariable[][] w = new MPVariable[numWorkers][numTasks];
		//s[i] and os[i] are the variables s and overline{s}. 
		MPVariable[] s = new MPVariable[numWorkers];
		MPVariable[] os = new MPVariable[numTasks];
		for (int i = 0; i < numWorkers; ++i) {
			s[i] = solver.makeNumVar(0, INFINITY, "s["+i+"]");
			os[i] = solver.makeNumVar(0, INFINITY, "os["+i+"]");
			for (int j = 0; j < numTasks; ++j) {
				//x[i][j] = solver.makeIntVar(0, 1, "");
				w[i][j] = solver.makeNumVar(0, INFINITY, "w["+i+"]["+j+"]");
			}
		}


		//Create the constraints
		//First set of constraints for each i
		for (int i = 0; i < numWorkers; ++i) {
			//I compute the value c_i^+ + d_i^- on the rhs of the constraints
			ISpecies spi=affineCRN.getSpecies().get(i);
			BigDecimal cip = c.get(spi);
			if(cip==null || cip.compareTo(BigDecimal.ZERO)<0) {
				cip=BigDecimal.ZERO;
			}
			BigDecimal dim = d.get(spi);
			if(dim==null || dim.compareTo(BigDecimal.ZERO)>0) {
				dim=BigDecimal.ZERO;
			}
			else {
				dim=dim.abs();
			}
			double cip_dim=cip.add(dim).doubleValue();
			
			//I provide cipdim as lower and upper, making it: ... = cipdim
			MPConstraint constraint = solver.makeConstraint(cip_dim, cip_dim, "constraint_i_"+i);
			for (int j = 0; j < numTasks; ++j) {
				constraint.setCoefficient(w[i][j], 1);
			}
			constraint.setCoefficient(s[i], 1);
		}
		//Second set of constraints for each j
		for (int j = 0; j < numTasks; ++j) {
			//I compute the value c_j^- + d_j^+ on the rhs of the constraints
			ISpecies spj=affineCRN.getSpecies().get(j);
			BigDecimal cjm = c.get(spj);
			if(cjm==null || cjm.compareTo(BigDecimal.ZERO)>0) {
				cjm=BigDecimal.ZERO;
			}
			else {
				cjm=cjm.abs();
			}
			
			BigDecimal djp = d.get(spj);
			if(djp==null || djp.compareTo(BigDecimal.ZERO)<0) {
				djp=BigDecimal.ZERO;
			}
			double cjm_djp=cjm.add(djp).doubleValue();
			
			//I provide cjm_djp as lower and upper, making it: ... = cjm_djp
			MPConstraint constraint = solver.makeConstraint(cjm_djp, cjm_djp, "constraint_j_"+j);
			for (int i = 0; i < numWorkers; ++i) {
				constraint.setCoefficient(w[i][j], 1);
			}
			constraint.setCoefficient(os[j], 1);
		}
		//The constraints >= 0 are obtained implicitly when declaring the variables


		//Create the objective function
		MPObjective objective = solver.objective();
		for (int i = 0; i < numWorkers; ++i) {
			//s[i]
			objective.setCoefficient(s[i], lambda);
			
			//wij
			for (int j = 0; j < numTasks; ++j) {
				double Dij=getOrDefaultValueIfHMNull(D, i, j,defValD);
				objective.setCoefficient(w[i][j], Dij);
			}
		}
		for (int j = 0; j < numTasks; ++j) {
			objective.setCoefficient(os[j], lambda);
		}
		objective.setMinimization();

		
		//Invoke the solver
		MPSolver.ResultStatus resultStatus = solver.solve();

		
		// Check that the problem has a feasible solution.
		if (resultStatus == MPSolver.ResultStatus.OPTIMAL
		    || resultStatus == MPSolver.ResultStatus.FEASIBLE) {
//		  System.out.println("Total cost: " + objective.value());
//		  for (int i = 0; i < numWorkers; ++i) {
//		    for (int j = 0; j < numTasks; ++j) {
//		    	System.out.println("Worker " + i + " assigned to task " + j + " with weight "+w[i][j].solutionValue()+".  Cost = " + getOrDefaultValueIfHMNull(D, i, j,defValD));
//		    }
//		  }
		} else {
		  System.err.println("No solution found.");
		  throw new UnsupportedOperationException("Something went wrong with the transport schedule map optimization");
		}
		
		Matrix KlDcd = new Matrix(0,affineCRN.getSpecies());
		for (int i = 0; i < numWorkers; ++i) {
			ISpecies spi = affineCRN.getSpecies().get(i);
			for (int j = 0; j < numTasks; ++j) {
				ISpecies spj = affineCRN.getSpecies().get(j);
				double v = w[i][j].solutionValue();
				KlDcd.setValue(spi, spj, v);
			}
		}
		return new SolutionAndOptimalValue(KlDcd,objective.value());
	}

	/**
	 * Computes the policy delta_lambda^pi
	 * We don't necessarily need to compute vertex solutions
	 * @param lambda
	 * @param pi
	 * @return the policy
	 */
	protected SolutionAndOptimalValue computeInvariantImprovement(double lambda, HashMap<Pair, Matrix> pi) {
		//The following code creates the data for the problem.
		int numWorkers = affineCRN.getSpeciesSize();
		int numTasks = affineCRN.getSpeciesSize();
		//Not clear what the costs are hereThe costs are contained in D


		// Create the linear solver with the backend.
		String backend="SCIP";
		//String backend="GLOP";
		MPSolver solver = MPSolver.createSolver(backend);
		if (solver == null) {
			System.out.println("Could not create solver "+backend);
			return null;
		}


		//Create the variables
		//D[i][j] the matrix to be returned
		MPVariable[][] D = new MPVariable[numWorkers][numTasks];
		for (int i = 0; i < numWorkers; ++i) {
			for (int j = 0; j < numTasks; ++j) {
				D[i][j] = solver.makeNumVar(0, INFINITY, "D["+i+"]["+j+"]");
			}
		}


		//Create the constraints
		//The constraints >= 0 are obtained implicitly when declaring the variables
		//Second set of constraints
		for (int i = 0; i < numWorkers; ++i) {
			ISpecies spi=affineCRN.getSpecies().get(i);
			HashMap<ISpecies, BigDecimal> Ai = getRowOrEmptyOne(Arows,spi);
			for (int j = 0; j < numTasks; ++j) {
				ISpecies spj=affineCRN.getSpecies().get(j);
				HashMap<ISpecies, BigDecimal> Aj = getRowOrEmptyOne(Arows,spj);
				Pair ij = new Pair(spi, spj);
				Matrix pi_ij = pi.get(ij);

				//Compute the upper bound
				double sumhk = 0;
				for (int h = 0; h < numTasks; ++h) {
					for (int k = 0; k < numTasks; ++k) {
						double pi_ij_hk = getOrDefaultValueIfHMNull(pi_ij, h, k,0);
						sumhk+=pi_ij_hk;
					}
				}
				
				double bibj = absbiminusbj(i, j);
				
				BigDecimal sumkBD=BigDecimal.ZERO;
				for (int k = 0; k < numTasks; ++k) {
					ISpecies spk = affineCRN.getSpecies().get(k);
					//double pi_ij_hk = pi_ij.getValue(h, k);
					//sumhk+=pi_ij_hk;
					BigDecimal Aik = getBDOrZeroIfNull(Ai, spk).abs();
					BigDecimal Ajk = getBDOrZeroIfNull(Aj, spk).abs();
					sumkBD = sumkBD.add(Aik.add(Ajk));
				}
				double sumk=sumkBD.doubleValue();
				
				double ub_ij=2*lambda*sumhk - bibj - lambda*sumk;
				
				//Now we are ready to build the constraints for ij
				MPConstraint constraint = solver.makeConstraint(NEGATIVE_INFINITY, ub_ij, "constraint_"+i+"_"+j);
				//sum across all pairs h,k in n such that h,k! i,j
				for (int h = 0; h < numTasks; ++h) {
					for (int k = 0; k < numTasks; ++k) {
						if(i==h && j==k) {
							//we skip this pair (h,k)==(i,j)
						}
						else {
							double pi_ij_hk = getOrDefaultValueIfHMNull(pi_ij, h, k,0);
							constraint.setCoefficient(D[h][k], pi_ij_hk);
						}
					}
				}
				double pi_ij_ij = getOrDefaultValueIfHMNull(pi_ij, i, j,0);
				constraint.setCoefficient(D[i][j], pi_ij_ij - 1);
			}
		}
		
		
		//Create the objective function
		MPObjective objective = solver.objective();
		for (int i = 0; i < numWorkers; ++i) {
			for (int j = 0; j < numTasks; ++j) {
				objective.setCoefficient(D[i][j], 1);
			}
		}
		objective.setMinimization();

		
		
		//Invoke the solver
		MPSolver.ResultStatus resultStatus = solver.solve();
		// Check that the problem has a feasible solution.
		if (resultStatus == MPSolver.ResultStatus.OPTIMAL
				|| resultStatus == MPSolver.ResultStatus.FEASIBLE) {
//			

		} else {
			throw new UnsupportedOperationException("Something went wrong with the policy improvement LP."); //System.err.println("No solution found.");
		}

		double optValue = objective.value();
		
		//Compute the matrix to be returned
		Matrix deltalp = new Matrix(0,affineCRN.getSpecies());
		for (int i = 0; i < numWorkers; ++i) {
			ISpecies spi = affineCRN.getSpecies().get(i);
			for (int j = 0; j < numTasks; ++j) {
				ISpecies spj = affineCRN.getSpecies().get(j);
				deltalp.setValue(spi, spj, D[i][j].solutionValue());
			}
		}

		return new SolutionAndOptimalValue(deltalp,optValue);
	}

	private HashMap<ISpecies, BigDecimal> getRowOrEmptyOne(
			LinkedHashMap<ISpecies, HashMap<ISpecies, BigDecimal>> arows2, ISpecies spi) {
		HashMap<ISpecies, BigDecimal> ret = Arows.get(spi);
		if(ret==null) {
			return EMPTYROW;
		}
		return ret;
	}

	public double getOrDefaultValueIfHMNull(Matrix pi_ij, int h, int k, double defVal) {
		if(pi_ij==null) {

			return defVal;
		}
		else {
			if(defVal!=pi_ij.getDefault()) {
				throw new UnsupportedOperationException("matrix does not have the expected default value");
			}
			double ret=pi_ij.getValue(h, k);
			return ret;
		}
		
	}

	protected double absbiminusbj(int i, int j) {
		return (bentries[i].subtract(bentries[j])).abs().doubleValue();
	}

	private BigDecimal getBDOrZeroIfNull(HashMap<ISpecies, BigDecimal> Ai, ISpecies spk) {
		if(Ai==null) {
			return BigDecimal.ZERO;
		}
		BigDecimal Aik=Ai.get(spk);
		if(Aik==null) {
			Aik=BigDecimal.ZERO;
		}
		return Aik;
	}

	private boolean isUnaryReag_BinaryProd(ICRN affineCRN) {
		for(ICRNReaction r : affineCRN.getReactions()) {
			if(!r.getReagents().isUnary()) {
				return false;
			}
			if(!r.getProducts().isBinary()) {
				return false;
			}
		}
		return true;
	}
	
	public static void basicExampleOrTools(MessageConsoleStream out, BufferedWriter bwOut) {
		String msg=BasicExample.myBasicExample();
		CRNReducerCommandLine.print(out,bwOut,"Basic example or-tools:\n");
		CRNReducerCommandLine.print(out,bwOut,msg);
	}

}
