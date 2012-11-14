package org.omg;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import org.omg.tools.Util;
import org.openscience.cdk.tools.SaturationChecker;

public class PMG{
	/**Output File containing the list of graph. */
	static BufferedWriter outFile, rejectedFile, CDKFile;
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
	static int method = MolProcessor.OPTIMAL;
	static boolean hashMap = false;
	static boolean cdk = false;
	static boolean checkBad = true;
	private static boolean java7 = false;
		
	public static void main(String[] args) throws IOException{
		wFile = false;
		String out = "default_out.sdf";
		ArrayList<String> fragFiles = new ArrayList<>();
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
			else if(args[i].equals("-noBadSub")){
				checkBad = false;
			}
			else if(args[i].equals("-o")){
				out = args[++i];
				wFile = true;
			}
			else if(args[i].equals("-forkjoin")){
				java7 = true;
			}
			else if(args[i].equals("-fr")){
				if (fragFiles.size() == 0) {
					fragFiles.add(args[++i]);
				} else {
					System.out.println("Disregarding extra fragment: "+args[++i]);
				}				
			}
		}
		

		if (wFile) {
			outFile = new BufferedWriter(new FileWriter(out));
			rejectedFile =  new BufferedWriter(new FileWriter("Rejected_"+out));
			CDKFile = new BufferedWriter(new FileWriter("CDK_"+out));
		}
		startupMessages();
		long before = System.currentTimeMillis();
		ArrayList<String> atomSymbols = Util.parseFormula(formula);
		if (atomSymbols == null) System.exit(1);
		MolProcessor mp = new MolProcessor(atomSymbols, formula, method, (method==MolProcessor.SEM_CAN) && hashMap, cdk, checkBad);
		for (String fileName : fragFiles) {
//			mp.addFragment(fileName);
			mp.useFragment(fileName);
		}

		long finalCount;
		if (java7){
			ForkJoinPool forkJoinPool = new ForkJoinPool();
			finalCount = forkJoinPool.invoke(mp);
		} else {
			executor = new ThreadPoolExecutor(executorCount, executorCount, 0L, TimeUnit.MILLISECONDS, taskQueue);
			availThreads = new AtomicInteger(executorCount-1);
			startedTasks.getAndIncrement();
			pendingTasks.getAndIncrement();
			executor.execute(mp);
			wait2Finish();	// wait for all tasks to finish, close the output file and return the final count
			executor.shutdown();	
			finalCount = molCounter.get();
		}
		
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		System.out.println("Unique molecular graphs:  " + (finalCount-MolProcessor.duplicate.get()));
		if (method == MolProcessor.SEM_CAN || mp.frag) System.out.println("Duplicates removed in the end: "+MolProcessor.duplicate.get());
		if (cdk) System.out.println("Final molecule count: "+(finalCount-MolProcessor.duplicate.get()-rejectedByCDK.get())+" after rejecting "+rejectedByCDK.get()+" by CDK.");
		System.out.println("Duration: " + (after - before) + " milliseconds");
		System.out.println("Started Tasks: "+startedTasks.get());
	}

	private static void startupMessages() {
		System.out.println("Processing "+formula);
		switch(method){
		case MolProcessor.BRT_FRC: System.out.println("Using brute-force."); break;
		case MolProcessor.CAN_AUG: System.out.println("Using canonical augmentation with bliss as canonizer."); break;
		case MolProcessor.MIN_CAN: System.out.println("Using only minimization."); break;
		case MolProcessor.SEM_CAN: System.out.println("Using semi-canonization and "+(hashMap?"hash map.":"minimization in the end.")); break;
		case MolProcessor.OPTIMAL: System.out.println("Using semi-canonization and minimization at each step."); break;
		}
		if (executorCount > 1)    System.out.println("Parallel execution with "+executorCount+" threads.");
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
		System.out.println("You can further specify the following options. Note that if you don't specify a method (-m) then an optimal method will be used, which is a mix of semi-canonicity and minimization.");
		System.out.println("\t-m \tThe method to use: 0=semi-canonicity; 1=minimization; 2=canonical-augmentation; 3=brute-force");
		System.out.println("\t-p  \tThe number of parallel threads to use; by default will use all available cores");
		System.out.println("\t-o  \tThe name of the output file");
		System.out.println("\t-hashmap \tEnables using a hashmap with semi-canonicity instead of the minimizer");
		System.out.println("\t-cdk \tEnables using CDK for removing unacceptable molecular structures in the end.");
		System.out.println("\t-noBadSub \tDisables checking for bad substructures during the generation.");
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
			if (wFile) {outFile.close(); rejectedFile.close(); CDKFile.close(); }
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
