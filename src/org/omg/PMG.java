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
	static BufferedWriter outFile, matrixFile;
	static AtomicInteger availThreads;
	final static AtomicLong molCounter = new AtomicLong(0);
	final static AtomicLong pendingTasks = new AtomicLong(0);
	final static AtomicLong startedTasks = new AtomicLong(0);
	final static SaturationChecker satCheck = new SaturationChecker();
	final static LinkedBlockingQueue<Runnable> taskQueue = new LinkedBlockingQueue<Runnable>();
	static ThreadPoolExecutor executor;
	static int executorCount;
	static boolean wFile;
	static String formula = null;
		
	private static void startup (String formula) {
		ArrayList<String> atomSymbols = Util.parseFormula(formula);
		if (atomSymbols == null) System.exit(1);
		MolProcessor mp = new MolProcessor(atomSymbols);
		startedTasks.getAndIncrement();
		pendingTasks.getAndIncrement();
		executor.execute(mp);
	}

	public static void main(String[] args) throws IOException{
		wFile = false;
		String out = "default_out.sdf";
		executorCount = Runtime.getRuntime().availableProcessors();
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-p")){
				executorCount = Integer.parseInt(args[++i]);
				if (executorCount < 1) executorCount = 1;
			}
			else if(args[i].equals("-mf")){
				formula = args[++i];
			}
			else if(args[i].equals("-o")){
				out = args[++i];
				wFile = true;
			}
			else if(args[i].equals("-fr")){
				System.err.println("Fragments are not supported in this version.");
			}
		}
		
		if (formula == null) {
			help();
			System.exit(0);
		}
		
		if (wFile) {
			outFile = new BufferedWriter(new FileWriter(out));
			matrixFile = new BufferedWriter(new FileWriter(out+".txt"));
		}
		executor = new ThreadPoolExecutor(executorCount, executorCount, 0L, TimeUnit.MILLISECONDS, taskQueue);
		availThreads = new AtomicInteger(executorCount-1);

		System.out.println("CDK-free "+formula);
		if (MolProcessor.canAug)  System.out.println("Using canonical augmentation with bliss as canonizer.");
		if (MolProcessor.semiCan) System.out.println("Using semi-canonization and "+(MolProcessor.hashMap?"hash map.":"minimization."));
		if (executorCount > 1)    System.out.println("Parallel execution with "+executorCount+" threads.");
		long before = System.currentTimeMillis();
		startup(formula); 	
		wait2Finish();	// wait for all tasks to finish, close the output file and return the final count
		executor.shutdown();	
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		System.out.println("molecules:  " + molCounter.get());
		System.out.println("duplicates: "+MolProcessor.duplicate.get()+"; non-duplicates: "+(molCounter.get()-MolProcessor.duplicate.get()));
		System.out.println("Duration: " + (after - before) + " milliseconds");
		System.out.println("Started Tasks: "+startedTasks.get());
		System.in.read();
	}

	private static void help() {
		System.out.println("You can use the following options.");
		System.out.println("Specifying a formula is obligatory.");
		System.out.println("\t-mf \tA formula - the elemental composition");
		System.out.println("\t-p  \tThe paralellism degree (number of parallel threads)");
		System.out.println("\t-o  \tThe name of the output file");
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
			if (wFile) {outFile.close(); matrixFile.close();}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
