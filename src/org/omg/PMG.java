package org.omg;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
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

/**
 * Open Molecule Generation
 * The main class collecting parameters and setting global objects
 * 
 * @author julio
 */
public class PMG{
	/**Output File containing the list of graph. */
	BufferedWriter outFile;

	AtomicLong mol_counter;
	SaturationChecker satCheck = new SaturationChecker();
	private int atomCount;
	boolean parallelExecution = true;

	static int executorCount=6;
	AtomicLong startedTasks;
	MultiCoreExecutor executor;

	int nH;
	static boolean wFile = false;

	// TODO: This is only for checking duplicates. This should be unnecessary in a good version.
	private Set<String> molSet = Collections.synchronizedSet(new HashSet<String>());

	public PMG(){
		mol_counter = new AtomicLong(0);
		startedTasks = new AtomicLong(0);
	}
	
	boolean generateParallelTask(MolHelper2 molecule) {
		startedTasks.getAndIncrement();
		try {
			executor.execute(new Generator(this, molecule));
			return true;
		} catch (RejectedExecutionException e) {	// if it failed to execute in parallel (due to overload), then continue sequentially
			return false;
		}
	}

	private MolHelper2 initialize(String formula, String fragments, String output) {
		MolHelper2 mol;
		try {
			if (wFile) outFile = new BufferedWriter(new FileWriter(output));

			mol = new MolHelper2();

			System.out.println("PMG: Parallel processing of "+formula+ " started (using bliss as canonizer and with "+executorCount+" threads).");
			System.out.print("Current atom order is: ");
//			while (true) {
			if (fragments == null)
				nH = mol.initialize(formula);
			else 
				nH = mol.initialize(formula, fragments);
			for (IAtom atom:mol.acontainer.atoms()) System.out.print(atom.getSymbol());
			System.out.println();
//			if (formula.equals("C4H7NO3")) {
//				if (!mol.acontainer.getAtom(0).getSymbol().equals("C")) continue;
//				if (!mol.acontainer.getAtom(4).getSymbol().equals("N")) continue;
//			}
//			break;
//			}


			atomCount = mol.atomCount;
			return mol;
		} catch (IOException e) {
			System.err.println("Could not open "+output+" for writing. Continuting without file output...");
			wFile = false;
			e.printStackTrace();
		} catch (CDKException e) {
			e.printStackTrace();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return null;	// upon unsuccessful initialization
	}

	public static void main(String[] args) throws IOException{
		// parse the command-line arguments
		String formula = "C4H10";
		String fragments = null;
		String out = "default_out.sdf";
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-p")){
				executorCount = Integer.parseInt(args[++i]);
			}
			else if(args[i].equals("-mf")){
				formula = args[++i];
			}
			else if(args[i].equals("-o")){
				out = args[++i];
				wFile = true;
			}
			else if(args[i].equals("-fr")){
				System.err.println("No fragments currently supported in the parallel MG.");
			}
		}
		
		PMG pmg = new PMG();
		// TODO
		pmg.executor = new MultiExecutor(executorCount);
		MolHelper2 mol = pmg.initialize(formula, fragments, out);
		
		// do the real processing
		long before = System.currentTimeMillis();
		pmg.generateParallelTask(mol);
		pmg.wait2Finish();	// wait for all tasks to finish, and close the output files
		pmg.shutdown();	// shutdown the executor service(s)
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		System.out.println("molecules " + pmg.getFinalCount());
		System.out.println("Duration: " + (after - before) + " miliseconds\n");
	}



	private void wait2Finish() {
		int time = 0, min = 0;
		while (0 < startedTasks.get()){
			try {
				Thread.sleep(1000);
				if (60 == time++) {
					time = 0;
					System.out.println("Almost "+ ++min +" minute(s) passed and so far the count of molecules = "+mol_counter.get());
				}
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		try {
			if (wFile) outFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void shutdown() {
		executor.shutdown();
	}

	public long getFinalCount() {
		while (0 < startedTasks.get()){
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		return mol_counter.get();
	}
}
