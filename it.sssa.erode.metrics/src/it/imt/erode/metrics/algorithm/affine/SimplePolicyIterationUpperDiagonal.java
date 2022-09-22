package it.imt.erode.metrics.algorithm.affine;

import java.io.BufferedWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.eclipse.ui.console.MessageConsoleStream;


import it.imt.erode.crn.interfaces.ICRN;
import it.imt.erode.crn.interfaces.ISpecies;
import it.imt.erode.metrics.data.Matrix;
import it.imt.erode.metrics.data.SolutionAndOptimalValue;
import it.imt.erode.onthefly.Pair;
import it.imt.erode.partition.interfaces.IBlock;
import it.imt.erode.partition.interfaces.IPartition;

/**
 * DOES NOT WORK!
 * For hermann5 we get different results v
 * @author andrea
 *
 */
public class SimplePolicyIterationUpperDiagonal extends SimplePolicyIteration {
 
	public static void main(String[] args) {
		SimplePolicyIterationUpperDiagonal.basicExampleOrTools(null,null);
	}

	public SimplePolicyIterationUpperDiagonal(MessageConsoleStream out, BufferedWriter bwOut, ICRN crn) {
		super(out,bwOut,crn);				
	}

	@Override
	public Matrix iterativeStep(double lambda, HashMap<Pair, Matrix> pi, Matrix current_invariant) {
		boolean existsAPair = false;
		int pairsUpdated=0;
		int iteration=1;
		do {
			println("Iteration "+iteration+"...");
			existsAPair = false;
			pairsUpdated=0;
			for(int i=0;i<getAffineCRN().getSpeciesSize();i++){
				ISpecies spi=getAffineCRN().getSpecies().get(i);
				for(int j=i+1;j<getAffineCRN().getSpeciesSize();j++){
					ISpecies spj=getAffineCRN().getSpecies().get(j);

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
						
						Pair ji = new Pair(spj, spi);
						pi.put(ji, pi_ij);

						//Prepare for next iteration by updating the invariat
						SolutionAndOptimalValue solAndM = computeInvariantImprovement(lambda,pi);
						current_invariant= solAndM.getSolution();
						
						//double prev_new =getOrDefaultValueIfHMNull(current_invariant, i, j,0);
						//CRNReducerCommandLine.println(out, bwOut,"  Value for "+ij+" changed from "+prev+" to "+prev_new);
					}
				}
			}
			println("  Considered "+pairsUpdated+" pairs");
			iteration++;
		}while(existsAPair);
		return current_invariant;
	}
	
	@Override
	protected HashMap<Pair, Matrix> constructInitialPolicy(IPartition R,double lambda, Matrix Dr) {
		HashMap<Pair, Matrix> pi = new LinkedHashMap<Pair, Matrix>();

		//We explicitly compute the matrix Klambda(Dr)(Ai,Aj) for every pair in R, 
		//We do not explicitly compute the matrix 0 for pair not in R
		IBlock current = R.getFirstBlock();
		while(current!=null) {
			ArrayList<ISpecies> currentList=new ArrayList<>(current.getSpecies());
			for(int _i=0;_i<currentList.size();_i++) {
				ISpecies i = currentList.get(_i);
				for(int _j=_i;_j<currentList.size();_j++) {
					ISpecies j = currentList.get(_j);
					Pair ij = new Pair(i, j);
					HashMap<ISpecies, BigDecimal> Ai = getRowOfA(i);
					HashMap<ISpecies, BigDecimal> Aj = getRowOfA(j);
					SolutionAndOptimalValue solAndVal = computeTransportScheduleMap(lambda,Dr,1,Ai,Aj);
					Matrix pi_ij = solAndVal.getSolution();
					pi.put(ij, pi_ij);
					
					Pair ji = new Pair(j, i);
					pi.put(ji, pi_ij);
				}
			}
			/*
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
			*/
			current=current.getNext();
		}

		return pi;
	}

}
