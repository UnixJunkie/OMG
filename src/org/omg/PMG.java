package org.omg;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import org.omg.tools.Util;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.tools.SaturationChecker;

public class PMG{
	/**Output File containing the list of graph. */
	static BufferedWriter outFile;
	static AtomicInteger availThreads;
	final static AtomicLong molCounter = new AtomicLong(0);
	final static AtomicLong rejectedByCDK = new AtomicLong(0);
	final static AtomicLong pendingTasks = new AtomicLong(0);
	final static AtomicLong startedTasks = new AtomicLong(0);
	final static LinkedBlockingQueue<Runnable> taskQueue = new LinkedBlockingQueue<Runnable>();
	static ThreadPoolExecutor executor;
	static ExecutorService fileWriterExecutor;
	static int executorCount=0;
	static boolean wFile = false;
	static String formula = null;
	static int method = MolProcessor.OPTIMAL;
	static boolean hashMap = false;
	static boolean cdk = true;
	static boolean checkBad = false;
	private static boolean verbose = false;
	private static String goodlist = null; //prescribed fragments to use as filter after generation process
    private static String badlist = null; //badlist of fragments to use as filter after generation process. Molecules should not contain them 
    private static String fragFile = null;
    private static MolProcessor mp;
    
	public static void main(String[] args) throws IOException, CDKException{
		if (args.length == 0) {
			usage();
			System.exit(0);
		}
		String out = parseArguments(args);
		if (executorCount < 1) executorCount = Runtime.getRuntime().availableProcessors();	// default: use all cores available

		if (wFile) {
			outFile = new BufferedWriter(new FileWriter(out));
			fileWriterExecutor = Executors.newSingleThreadExecutor();	// use an independent thread to write to file
		}
		availThreads = new AtomicInteger(executorCount);	// number of tasks waiting in the queue; or, to generate in parallel
		startupMessages();
		
		long before = System.currentTimeMillis();
		ArrayList<String> atomSymbols = Util.parseFormula(formula);
		if (atomSymbols == null) System.exit(1);	// exit if it failed to parse the formula
		mp = new MolProcessor(atomSymbols, formula, method, (method==MolProcessor.SEM_CAN) && hashMap, cdk, checkBad, fragFile != null, goodlist, badlist);
		if (fragFile != null) mp.useFragment(fragFile);
		executor = new ThreadPoolExecutor(executorCount, executorCount, 0L, TimeUnit.MILLISECONDS, taskQueue);
		startedTasks.getAndIncrement();
		pendingTasks.getAndIncrement();
		executor.execute(mp);	// start the generation
		wait2Finish();			// wait for all tasks to finish (pendingTasks == 0), and close the output file 
		executor.shutdown();	// stop the executor --> kill the threads
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		report(after - before);
	}

	private static String parseArguments(String[] args){
		int i=0;
		String out=null;
		try{
			formula = args[0];
			for(i = 0; i < args.length; i++){
				if(args[i].equals("-p")){
					executorCount = Integer.parseInt(args[++i]);
				}
				else if(args[i].startsWith("-p")){
					executorCount = Integer.parseInt(args[i].substring(2));
				}
				else if(args[i].equals("-m")){
					method = Integer.parseInt(args[++i]);
				}
				else if(args[i].startsWith("-m")){
					method = Integer.parseInt(args[i].substring(2));
				}
				else if(args[i].equals("-hashmap")){
					hashMap = true;
				}
				else if(args[i].equals("-nocdk")){
					cdk = false;
				}
				else if(args[i].equals("-filter")){
					checkBad = true;
				}
				else if(args[i].equals("-o")){
					out = args[++i];
					wFile = true;
				}
				else if(args[i].equals("-v")){
					verbose = true;
				}
				else if(args[i].equals("-fr")){
					if (fragFile == null) {
						fragFile = args[++i];
					} else {
						System.out.println("Disregarding extra fragment: "+args[++i]+". Consider using them with -goodlist.");
					}				
				}
				else if(args[i].equals("-goodlist")){
					if (goodlist == null)
						goodlist = args[++i];
					else 
						System.out.println("Disregarding the extra goodlist: "+args[++i]);
				}
				else if(args[i].equals("-badlist")){
					if (badlist == null)
						badlist = args[++i];
					else 
						System.out.println("Disregarding the extra badlist: "+args[++i]);
				} else 
					System.out.println("Warning: Invalid arguments are ignored. Run the program with no arguments to see the usage.");
			}
		}catch(ArrayIndexOutOfBoundsException ae){
			System.err.println("Insufficient parameters provided. Run the program with no arguments to see the usage.");
			System.exit(1);
		}catch(NumberFormatException ne){
			System.err.println("Invalid number given: "+args[i]);
			System.exit(1);
		}
		return out;
	}

	private static void usage() {
		System.out.println("Usage: PMG <formula> [options]\n");
		System.out.println("Providing a formula for the elemental compositoin is obligatory, e.g., C4H7NO3.");
		System.out.println("You can further specify the following options.");
		System.out.println("\t-filter \tFilter bad substructures in the molecular structure.");
		System.out.println("\t-fr \tA file containing one substructure to use as initial structure for generation");
		System.out.println("\t-goodlist\tA file containing required substructures of the molecule (checked in the end) - only active if cdk is used");
		System.out.println("\t-badlist\tA file containing forbidden substructures (checked in the end) - only active if cdk is used");
		System.out.println("\t-hashmap\tEnables using a hashmap with semi-canonicity instead of the minimizer");
		System.out.println("\t-nocdk \tDisables using CDK for removing unacceptable molecular structures in the end.");
		System.out.println("\t-p  \tThe number of parallel threads to use; by default will use all available cores");
		System.out.println("\t-o  \tThe name of the output file");
		System.out.println("\t-v  \tverbose mode");
		System.out.println("\t-m \tThe method to use: \n\t\t\t0=semi-canonicity; \n\t\t\t1=minimization; \n\t\t\t2=canonical-augmentation;");
		System.out.println("Note that if you don't specify a method (-m) then an optimal method will be used, which is a mix of semi-canonicity and minimization.");
	}

	private static void startupMessages() {
		System.out.print("Processing "+formula+" using");
		switch(method){
		case MolProcessor.CAN_AUG: System.out.println(" canonical augmentation with bliss as canonizer"); break;
		case MolProcessor.MIN_CAN: System.out.println(" only minimization"); break;
		case MolProcessor.SEM_CAN: System.out.println(" semi-canonization and "+(hashMap?"hash map":"minimization in the end")); break;
		case MolProcessor.OPTIMAL: System.out.println(" semi-canonization and minimization at each step"); break;
		}
		System.out.print("Selected options: ");
		if (cdk) System.out.print("cdk"+(goodlist==null?"":" with good-list")+(badlist==null?"":(goodlist==null?" with":" and")+" bad-list")+", ");
		if (checkBad) System.out.print("bad-sub filter, ");
		if (wFile) System.out.print("with output to file, ");
		System.out.println((executorCount > 1)? ""+executorCount+" threads." : "sequential execution.");
		if (!cdk && (goodlist!= null || badlist != null)) System.out.println("Warning: The provided good- and/or bad-list will not be used because CDK is disabled.");
	}

	private static void report(long duration) {
		if (verbose) {
			if (method == MolProcessor.SEM_CAN || mp.frag) System.out.println("Duplicates removed in the end: "+MolProcessor.duplicate.get());
			if (cdk) System.out.println("Rejecting "+rejectedByCDK.get()+" by CDK.");
			System.out.println("Started Tasks: "+startedTasks.get());
		}
		System.out.println("Unique molecule count:  " + molCounter.get());
		System.out.println("Duration: " + duration + " milliseconds");
	}

	private static void wait2Finish() {
		int timer = 60;
		while (0 < pendingTasks.get()){
			try {
				Thread.sleep(1000);
				if (timer > 0)
					timer--;
				else {
					System.out.println("Molecules generated so far: "+molCounter.get());
					timer = 60;
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		try {
			if (wFile) { 
				fileWriterExecutor.shutdown(); 
				while (!fileWriterExecutor.awaitTermination(60, TimeUnit.SECONDS));
				outFile.close(); 
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			System.err.println("Failed while waiting for file outputs to complete.");
			e.printStackTrace();
		}
	}

}
