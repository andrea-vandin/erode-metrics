package it.imt.erode.metrics.data;

public class SolutionAndOptimalValue {

	private Matrix solution;
	private double optimalValue;
	
	public SolutionAndOptimalValue(Matrix m, double value){
		this.solution=m;
		this.optimalValue=value;
	}

	public Matrix getSolution() {
		return solution;
	}

	public double getValue() {
		return optimalValue;
	}
	
	@Override
	public String toString() {
		return optimalValue + "\n" + solution.toString(); 
	}
	
}
