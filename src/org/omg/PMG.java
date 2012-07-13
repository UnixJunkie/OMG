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
	boolean parallelExecution = true;

	final int executorCount;
	MultiCoreExecutor executor;

	int nH;
	static boolean wFile = false;

	// TODO: This is only for checking duplicates. This should be unnecessary in a good version.
	private Set<String> molSet = Collections.synchronizedSet(new HashSet<String>());

	public PMG(int ec){
		mol_counter = new AtomicLong(0);
		executorCount = ec;
	}
	
	boolean generateParallelTask(Runnable molecule, boolean force) {
		try {
			executor.execute(molecule, force);
			return true;
		} catch (RejectedExecutionException e) {	// if it failed to execute in parallel (due to overload), then continue sequentially
			return false;
		}
	}
		
	boolean generateParallelTask(MolHelper2 molecule, boolean force) {
		return generateParallelTask(new Generator(this, molecule), force);
	}

	private void startup (String formula) {
		ArrayList<String> atomSymbols = Util.parseFormula(formula);
		if (atomSymbols == null) System.exit(1);
		MolProcessor mp = new MolProcessor(atomSymbols, this);
		executor.execute(mp);
	}
/*	
	private void addAtomAndStart(Atom[] atomsInMolecule, ArrayList<String> atomSymbols, int i) {
		List<Atom> atomPossibilities = Atom.create(atomSymbols.get(i));
		for (Atom a:atomPossibilities) {
			atomsInMolecule[i] = a;
			if (i+1<atomSymbols.size()) {
				addAtomAndStart(atomsInMolecule, atomSymbols, i+1);
			} else {
				MolProcessor mp = new MolProcessor(atomsInMolecule, Util.nH);
				startedTasks.getAndIncrement();
				executor.execute(mp);
			}
		}
	}
*/
	private MolHelper2 initialize(String formula, String fragments, String output) {
		MolHelper2 mol;
		try {
			mol = new MolHelper2();

			System.out.println("PMG: Parallel processing of "+formula+ " started (using bliss as canonizer and with "+executorCount+" threads).");
			System.out.print("Current atom order is: ");
			if (fragments == null)
				nH = mol.initialize(formula);
			else 
				nH = mol.initialize(formula, fragments);
			for (IAtom atom:mol.acontainer.atoms()) System.out.print(atom.getSymbol());
			System.out.println();
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
		String formula = null;
		String fragments = null;
		String out = "default_out.sdf";
		int eCount = Runtime.getRuntime().availableProcessors();
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-p")){
				eCount = Integer.parseInt(args[++i]);
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
		
		if (formula == null) {
			help();
			System.exit(0);
		}
		
		new PMG(eCount).calculate(formula, out);
	}

	/**
	 * @param formula
	 * @param out
	 * @param pmg
	 * @throws IOException
	 */
	private void calculate(String formula, String out)
			throws IOException {
		// TODO: Enable using MultiExecutor
		executor = new SingleExecutor(executorCount);
//		pmg.executor = new MultiExecutor(executorCount);			
		
		if (wFile) outFile = new BufferedWriter(new FileWriter(out));

		startup(formula); 
		
//		MolHelper2 mol = pmg.initialize(formula, fragments, out);
		
		// do the real processing
		long before = System.currentTimeMillis();
//		pmg.startedTasks.getAndIncrement();
//		pmg.generateParallelTask(mol, true);
		long finalCount = wait2Finish();	// wait for all tasks to finish, close the output file and return the final count
		shutdown();	// shutdown the executor service(s)
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		System.out.println("molecules " + finalCount);
//		System.out.println("N=3 " + N3.get());
		System.out.println("Duration: " + (after - before) + " milliseconds\n");
	}



	private static void help() {
		// TODO Auto-generated method stub
		System.out.println("You can use the following options.");
		System.out.println("Specifying a formula is obligatory.");
		System.out.println("\t-mf \tA formula - the elemental composition");
		System.out.println("\t-p  \tThe paralellism degree (number of parallel threads)");
		System.out.println("\t-o  \tThe name of the output file");
	}

	private long wait2Finish() {
		int time = 0, min = 0;
		while (executor.busy()){
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
		return mol_counter.get();
	}

	private void shutdown() {
		executor.shutdown();
	}

}
