package org.omg;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.omg.tools.*;

public class PMG{
	/**Output File containing the list of graph. */
	static BufferedWriter outFile, matrixFile, rejectedFile, CDKFile;
	static AtomicInteger availThreads;
	final static AtomicLong molCounter = new AtomicLong(0);
	final static AtomicLong rejectedByCDK = new AtomicLong(0);
	final static AtomicLong pendingTasks = new AtomicLong(0);
	final static AtomicLong startedTasks = new AtomicLong(0);
	final static SaturationChecker satCheck = new SaturationChecker();
	final static LinkedBlockingQueue<Runnable> taskQueue = new LinkedBlockingQueue<Runnable>();
	static ThreadPoolExecutor executor;
	static int executorCount;
	static boolean wFile;
	static String formula = null;
	static int method = MolProcessor.MIN_CAN;
	static boolean hashMap = false;
	static boolean cdk = false;
		
	private static void startup () {
		ArrayList<String> atomSymbols = Util.parseFormula(formula);
		if (atomSymbols == null) System.exit(1);
		MolProcessor mp = new MolProcessor(atomSymbols, formula, method, hashMap, cdk);
		startedTasks.getAndIncrement();
		pendingTasks.getAndIncrement();
		executor.execute(mp);
	}

	public static void main(String[] args) throws IOException{
		wFile = false;
		String out = "default_out.sdf";
		executorCount = Runtime.getRuntime().availableProcessors();
		if (args.length == 0) {
			usage();
		}
		formula = args[0];
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-p")){
				executorCount = Integer.parseInt(args[++i]);
				if (executorCount < 1) executorCount = 1;
			}
			else if(args[i].equals("-m")){
				method = interpretMethod(args[++i]);
			}
			else if(args[i].equals("-hashmap")){
				hashMap = true;
			}
			else if(args[i].equals("-cdk")){
				cdk = true;
			}
			else if(args[i].equals("-o")){
				out = args[++i];
				wFile = true;
			}
			else if(args[i].equals("-fr")){
				System.err.println("Fragments are not supported in this version.");
			}
		}
		

		if (wFile) {
			outFile = new BufferedWriter(new FileWriter(out));
			rejectedFile =  new BufferedWriter(new FileWriter("Rejected_"+out));
			CDKFile = new BufferedWriter(new FileWriter("CDK_"+out));
			matrixFile = new BufferedWriter(new FileWriter(out+".txt"));
		}
		executor = new ThreadPoolExecutor(executorCount, executorCount, 0L, TimeUnit.MILLISECONDS, taskQueue);
		availThreads = new AtomicInteger(executorCount-1);

		System.out.println("CDK-free "+formula);
		switch(method){
		case MolProcessor.BRT_FRC: System.out.println("Using brute-force."); break;
		case MolProcessor.CAN_AUG: System.out.println("Using canonical augmentation with bliss as canonizer."); break;
		case MolProcessor.MIN_CAN: System.out.println("Using only minimization."); break;
		case MolProcessor.SEM_CAN: System.out.println("Using semi-canonization and "+(hashMap?"hash map.":"minimization.")); break;
		}
		if (executorCount > 1)    System.out.println("Parallel execution with "+executorCount+" threads.");
		long before = System.currentTimeMillis();
		startup(); 	
		wait2Finish();	// wait for all tasks to finish, close the output file and return the final count
		executor.shutdown();	
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		System.out.println("Unique molecular graphs:  " + (molCounter.get()-MolProcessor.duplicate.get()));
		if (method == MolProcessor.SEM_CAN) System.out.println("Duplicates removed in the end: "+MolProcessor.duplicate.get());
		if (cdk) System.out.println("Final molecule count: "+(molCounter.get()-MolProcessor.duplicate.get()-rejectedByCDK.get())+" after rejecting "+rejectedByCDK.get()+" by CDK.");
		System.out.println("Duration: " + (after - before) + " milliseconds");
		System.out.println("Started Tasks: "+startedTasks.get());
	}

	private static int interpretMethod(String methodStr) {
		int m=0;
		try{
			m = Integer.parseInt(methodStr);
		} catch (NumberFormatException nfe) {
			System.err.println("Invalid method.");
			usage();
		}
		return m;
	}

	private static void usage() {
		System.out.println("Usage: PMG <formula> [options]");
		System.out.println("Providing a formula for the elemental compositoin is obligatory, e.g., C4H7NO3.");
		System.out.println("You can further specify the following options.");
		System.out.println("\t-m \tThe method to use: 0=semi-canonicity; 1=minimization; 2=canonical-augmentation; 3=brute-force");
		System.out.println("\t-p  \tThe number of parallel threads to use; by default will use all available cores");
		System.out.println("\t-o  \tThe name of the output file");
		System.out.println("\t-hashmap \tEnables using a hashmap with semi-canonicity instead of the minimizer");
		System.out.println("\t-cdk \tEnables using CDK for removing unacceptable molecular structures in the end.");
		System.exit(0);
	}

	private static void wait2Finish() {
		while (0 < pendingTasks.get()){
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		try {
			if (wFile) {outFile.close(); matrixFile.close(); rejectedFile.close(); CDKFile.close(); }
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
