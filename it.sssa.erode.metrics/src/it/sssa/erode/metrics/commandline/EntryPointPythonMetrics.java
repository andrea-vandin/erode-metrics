package it.sssa.erode.metrics.commandline;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import org.eclipse.ui.console.MessageConsoleStream;

import com.microsoft.z3.Z3Exception;

import it.imt.erode.auxiliarydatastructures.IPartitionAndBoolean;
import it.imt.erode.commandline.CRNReducerCommandLine;
import it.imt.erode.commandline.CommandsReader;
import it.imt.erode.commandline.EntryPointForPython;
import it.imt.erode.crn.interfaces.ICRN;
import it.imt.erode.crn.interfaces.ISpecies;
import it.imt.erode.importing.UnsupportedFormatException;
import it.imt.erode.onthefly.Pair;
import it.imt.erode.partition.interfaces.IPartition;
import it.sssa.erode.or_tools.StringAndPairs;
import py4j.GatewayServer;

public class EntryPointPythonMetrics {

	private CRNReducerCommandLine erode;
	private MessageConsoleStream out = null;
	private BufferedWriter bwOut=null;
	private boolean printPartitions;
	private boolean printModels;
	private LinkedHashMap<String, ISpecies> nameToSpecies;
	private ISpecies[] idToSpecies;
	private String[] idToSpeciesNames;
	private boolean modelLoaded=false;
	private IPartition defaultPartition;
	
	public static void main(String[] args) {
		EntryPointPythonMetrics entry = new EntryPointPythonMetrics(true, true);
		
		int port=-1;
		if(args!=null && args.length>0) {
			port = Integer.valueOf(args[0]);
		}
		GatewayServer gatewayServer = null;
		if(port==-1){
			gatewayServer=new GatewayServer(entry);
		}
		else{
			gatewayServer=new GatewayServer(entry,port);
		}
				
		gatewayServer.start();
		
		CRNReducerCommandLine.println(entry.out,entry.bwOut,"Gateway server started on port "+gatewayServer.getPort()+ " (while pythonPort is "+gatewayServer.getPythonPort()+", and pythonAddress is "+gatewayServer.getPythonAddress().toString()+")");

		/*
		entry.loadModel("distr/AM.ode");
		try {
			int[] obtained = entry.computeBB();
			System.out.println("The partition"+Arrays.toString(obtained));
		} catch (Z3Exception | UnsupportedFormatException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		 
		try {
			HashMap<Pair, Double> pairs = entry.computeMetrics();
			System.out.println("The pairs\n"+pairs);
		} catch (UnsupportedFormatException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		*/
	}
	
	private void checkModelLoaded() throws UnsupportedOperationException {
		if(!modelLoaded) {
			CRNReducerCommandLine.println(out,bwOut,"Please first load a model.");
			throw new UnsupportedOperationException("Please first load a model.");
		}
	}
	
	public EntryPointPythonMetrics(boolean printPartitions, boolean printModels) {
		this.printPartitions=printPartitions;
		this.printModels=printModels;
		erode=new CRNReducerCommandLine(new CommandsReader(new ArrayList<String>(0),out,bwOut));
		
		CRNReducerCommandLine.println(out,bwOut,"ERODE instantiated");
	}
	
	public String getModelString() {
		checkModelLoaded();
		return erode.getCRN().toString();
	}
	
	public void printPartition(int[] partitionArray){
		CRNReducerCommandLine.println(out,bwOut,getPartitionString(partitionArray));
	}
	public String getPartitionString(int[] partitionArray){
		checkModelLoaded();
		boolean numbersAreIDOfRepresentativeSpecies=false;
		IPartition partition = EntryPointForPython.importPartition(idToSpecies, partitionArray,numbersAreIDOfRepresentativeSpecies);
		return partition.toString();
	}
	
	public void populateAuxiliarySpeciesDataStructures(ICRN crn) {
		nameToSpecies=new LinkedHashMap<>(crn.getSpecies().size());
		idToSpecies=new ISpecies[crn.getSpecies().size()];
		idToSpeciesNames=new String[crn.getSpecies().size()];
		int i=0;
		for(ISpecies sp : crn.getSpecies()) {
			nameToSpecies.put(sp.getName(), sp);
			idToSpecies[i]=sp;
			idToSpeciesNames[i]=sp.getName();
			i++;
		}
	}
	
	private int completeImporting() {
		modelLoaded=true;
		CRNReducerCommandLine.println(out,bwOut);
		populateAuxiliarySpeciesDataStructures(erode.getCRN());
		defaultPartition=erode.getPartition();
		
//		idToSpecies=new ISpecies[erode.getCRN().getSpecies().size()];
//		idToSpeciesNames=new String[idToSpecies.length];
//		int i=0;
//		for (ISpecies species : erode.getCRN().getSpecies()) {
//			idToSpecies[i]=species;
//			idToSpeciesNames[i]=species.getName();
//			i++;
//		}

		/*for(i=0;i<idToSpecies.length;i++){
			ISpecies species = idToSpecies[i];
			species.setInitialConcentration(BigDecimal.ONE, "1.0");
		}*/

		if(printModels){
			CRNReducerCommandLine.println(out,bwOut,erode.getCRN());
		}

		return erode.getCRN().getSpecies().size();
	}
	
	public int loadModel(String fileName){
		int ret=-1;
		if(fileName.endsWith(".ode")){
			ret= importERODE(fileName);
		}
		else if(fileName.endsWith(".net")){
			ret= importBNG(fileName);
		}
		else {
			throw new UnsupportedOperationException("Either .crn or .net formats are supported");
		}
		return ret;
	}
	
	private int importBNG(String fileName){
		erode.handleImportBNGCommand("importBNG({fileIn=>"+fileName+"})",out,bwOut);
		return completeImporting();
	}
	private int importERODE(String fileName){
		erode.handleLoadCommand("load({fileIn=>"+fileName+"})",false,out,bwOut);
		return completeImporting();
	}
	
//	private int completeImporting() {
//		CRNReducerCommandLine.println(out,bwOut);
//		idToSpecies=new ISpecies[erode.getCRN().getSpecies().size()];
//		idToSpeciesNames=new String[idToSpecies.length];
//		int i=0;
//		for (ISpecies species : erode.getCRN().getSpecies()) {
//			idToSpecies[i]=species;
//			idToSpeciesNames[i]=species.getName();
//			i++;
//		}
//
//		/*for(i=0;i<idToSpecies.length;i++){
//			ISpecies species = idToSpecies[i];
//			species.setInitialConcentration(BigDecimal.ONE, "1.0");
//		}*/
//
//		if(printModels){
//			CRNReducerCommandLine.println(out,bwOut,erode.getCRN());
//		}
//
//		return erode.getCRN().getSpecies().size();
//	}
	
	public String[]/*ArrayList<String>*//*void*//*String*/ computeMetrics() throws UnsupportedFormatException, IOException{
		checkModelLoaded();
		StringAndPairs ret = erode.handleMetrics("", out, bwOut);
		HashMap<Pair,Double> pairs=ret.getPairs();
		String[] pairsString=new String[pairs.size()+1];
		StringBuffer sb = new StringBuffer();
		pairsString[0]=ret.getMsg();//"(test,test)=4";
		int i=1;
		for(Entry<Pair, Double> entry : pairs.entrySet()) {
			String s = entry.getKey().toString()+"="+entry.getValue()+"\n";
			pairsString[i]=s;
			sb.append(s);
			i++;
		}
		return pairsString;
		//return sb.toString();
	}
	
	public int[] computeBB() throws UnsupportedFormatException, Z3Exception, IOException{
		int[] initialPartitionArray = new int[idToSpecies.length];
		Arrays.fill(initialPartitionArray, 1);
		return computeBB(initialPartitionArray,false);
	}
	
	public int[] computeBB(int[] initialPartitionArray) throws UnsupportedFormatException, Z3Exception, IOException{
		return computeBB(initialPartitionArray,false);
	}
	
	public int[] computeBB(int[] initialPartitionArray, boolean numbersAreIDOfRepresentativeSpecies) throws UnsupportedFormatException, Z3Exception, IOException{
		checkModelLoaded();
		CRNReducerCommandLine.println(out,bwOut,"Computing BB reduction");
		
		IPartition initialPartition = EntryPointForPython.importPartition(idToSpecies, initialPartitionArray,numbersAreIDOfRepresentativeSpecies);
		erode.setPartition(initialPartition);
		if(printPartitions){
			CRNReducerCommandLine.println(out,bwOut,"Initial partition:\n"+initialPartition);
		}
		
		


		IPartitionAndBoolean obtainedPartitionAndBool = erode.handleReduceCommand("reduceBE({computeOnlyPartition=>true,print=>false})",false,"be",out,bwOut);
		IPartition obtainedPartition = obtainedPartitionAndBool.getPartition();
		int[] obtainedPartitionToExport=null;
		try {
		//IPartition obtainedPartition = crnreducer.handleReduceCommand("reduceEFL({computeOnlyPartition=>true,print=>false})",false,"EFL");
		if(printPartitions){
			CRNReducerCommandLine.println(out,bwOut,"Obtained partition:\n"+obtainedPartition);
		}
	
		obtainedPartitionToExport = EntryPointForPython.exportPartition(idToSpecies,obtainedPartition);
		CRNReducerCommandLine.println(out,bwOut,"BBE reduction completed");
		}finally {
			erode.setPartition(defaultPartition);
		}
		return obtainedPartitionToExport;
	}
	
}
