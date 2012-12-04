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
	final static SaturationChecker satCheck = new SaturationChecker();
	final static LinkedBlockingQueue<Runnable> taskQueue = new LinkedBlockingQueue<Runnable>();
	static ThreadPoolExecutor executor;
	static ExecutorService fileWriterExecutor;
	static int executorCount;
	static boolean wFile;
	static String formula = null;
	static int method = MolProcessor.OPTIMAL;
	static boolean hashMap = false;
	static boolean cdk = true;
	static boolean checkBad = false;
	private static boolean java7 = false, verbose = false;
	private static String goodlist = null; //prescribed fragments to use as filter after generation process
    private static String badlist = null; //badlist of fragments to use as filter after generation process. Molecules should not contain them 
    private static MolProcessor mp;
    
	public static void main(String[] args) throws IOException, CDKException{
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
			else if(args[i].equals("-forkjoin")){
				java7 = true;
			}
			else if(args[i].equals("-v")){
				verbose = true;
			}
			else if(args[i].equals("-fr")){
				if (fragFiles.size() == 0) {
					fragFiles.add(args[++i]);
				} else {
					System.out.println("Disregarding extra fragment: "+args[++i]);
				}				
			}
			else if(args[i].equals("-goodlist")){
                try {
                    goodlist = args[++i];
                } catch (Exception e) {
                    System.err.println("No file with posterior substructures provided");
                    System.exit(1);
                }                   
            }
            else if(args[i].equals("-badlist")){
                try {
                    badlist = args[++i];
                } catch (Exception e) {
                    System.err.println("No file with posterior substructures provided");
                    System.exit(1);
                }                   
            }
		}
		

		if (wFile) {
			outFile = new BufferedWriter(new FileWriter(out));
			//executorCount *= 2;		// in case of IO-bound, we increase the number of threads
			fileWriterExecutor = Executors.newSingleThreadExecutor();
		}
		availThreads = new AtomicInteger(executorCount);	// number of tasks waiting in the queue; or, to generate in parallel
		startupMessages();
		long before = System.currentTimeMillis();
		ArrayList<String> atomSymbols = Util.parseFormula(formula);
		if (atomSymbols == null) System.exit(1);
		mp = new MolProcessor(atomSymbols, formula, method, (method==MolProcessor.SEM_CAN) && hashMap, cdk, checkBad, fragFiles.size() != 0, goodlist, badlist);
		for (String fileName : fragFiles) {
//			mp.addFragment(fileName);
			mp.useFragment(fileName);
		}

		long finalCount=0;
		if (java7){
			System.out.println("Fork/Join not supported.");
//			ForkJoinPool forkJoinPool = new ForkJoinPool();
//			finalCount = forkJoinPool.invoke(mp);
		} else {
			executor = new ThreadPoolExecutor(executorCount, executorCount, 0L, TimeUnit.MILLISECONDS, taskQueue);
			startedTasks.getAndIncrement();
			pendingTasks.getAndIncrement();
			executor.execute(mp);
			wait2Finish();	// wait for all tasks to finish, close the output file and return the final count
			executor.shutdown();	
			finalCount = molCounter.get();
		}
		
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		System.out.println("Unique molecule count:  " + (finalCount));
		if (verbose && (method == MolProcessor.SEM_CAN || mp.frag)) System.out.println("Duplicates removed in the end: "+MolProcessor.duplicate.get());
		if (cdk) {
			if (verbose) System.out.println("Rejecting "+rejectedByCDK.get()+" by CDK.");
		}
		System.out.println("Duration: " + (after - before) + " milliseconds");
		if (verbose) System.out.println("Started Tasks: "+startedTasks.get());
	}

	private static void startupMessages() {
		System.out.print("Processing "+formula+" using");
		switch(method){
//		case MolProcessor.BRT_FRC: System.out.println(" brute-force"); break;
		case MolProcessor.CAN_AUG: System.out.println(" canonical augmentation with bliss as canonizer"); break;
		case MolProcessor.MIN_CAN: System.out.println(" only minimization"); break;
		case MolProcessor.SEM_CAN: System.out.println(" semi-canonization and "+(hashMap?"hash map":"minimization in the end")); break;
		case MolProcessor.OPTIMAL: System.out.println(" semi-canonization and minimization at each step"); break;
		}
		System.out.print("Selected options: ");
		if (cdk) System.out.print("cdk"+(goodlist==null?"":" with good-list")+(badlist==null?"":(goodlist==null?" with":" and")+" bad-list")+", ");
		if (checkBad) System.out.print("bad-sub filter, ");
		if (executorCount > 1) System.out.print(""+executorCount+" threads, ");
		if (wFile) System.out.print("with output to file, ");
		System.out.print("\b\b.\n");
		if (!cdk && (goodlist!= null || badlist != null)) System.out.println("Warning: The provided good- and/or bad-list will not be used because CDK is disabled.");
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
		System.out.println("\t-badlist\tA file containing forbidden substructures (checked in the end) - only active if cdk is used");
		System.out.println("\t-filter \tFilter bad substructures in the molecular structure.");
		System.out.println("\t-fr \tA file containing one substructure to use as initial structure for generation");
		System.out.println("\t-goodlist\tA file containing requried substructures of the molecule (checked in the end) - only active if cdk is used");
		System.out.println("\t-hashmap\tEnables using a hashmap with semi-canonicity instead of the minimizer");
		System.out.println("\t-m \tThe method to use: 0=semi-canonicity; 1=minimization; 2=canonical-augmentation; 3=brute-force");
		System.out.println("\t-nocdk \tDisables using CDK for removing unacceptable molecular structures in the end.");
		System.out.println("\t-p  \tThe number of parallel threads to use; by default will use all available cores");
		System.out.println("\t-o  \tThe name of the output file");
		System.out.println("\t-v  \tverbose mode");
		System.exit(0);
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
