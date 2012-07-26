package org.omg;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.tools.SaturationChecker;
/**
 * Open Molecule Generation
 * The main class collecting parameters and setting global objects
 * 
 * @author mahdi
 */
public class PMG7{
	/**Output File containing the list of molecules. */
	static BufferedWriter outFile;
	private final static ForkJoinPool forkJoinPool = new ForkJoinPool();
	
	private static class Generator extends RecursiveTask<Long> {
		final MolHelper2 mol;
		final static SaturationChecker satCheck = new SaturationChecker();
		final boolean wFile;
		final int nH;

		public Generator(final MolHelper2 mol, final boolean wf, final int nh) {
			super();
			this.mol = mol;
			this.wFile = wf;
			this.nH = nh;
		}

		@Override
		protected Long compute() {
			long count = 0;
			try {
				if (mol.isComplete(satCheck, nH) && mol.isConnected()) {
//					if (!molSet.add(mol.canString)) System.err.println("Duplicate");
					count = 1;
					if(wFile){
						mol.writeMol(outFile);
					}
				}
				ArrayList<MolHelper2> extMolList = mol.addOneBond();
				List<Generator> tasks = new ArrayList<>();
				for (final MolHelper2 child : extMolList) {
					tasks.add(new Generator(child, wFile, nH)); 
				}
				for(final Generator task : invokeAll(tasks)) { 
					count += task.join(); 
				}
			} catch (CloneNotSupportedException e){
				e.printStackTrace();
			} catch (CDKException e) {
				e.printStackTrace();
			}		
			return count;			
		}
		
	}
	

	// TODO: This is only for checking duplicates. This should be unnecessary in a good version.
	private Set<String> molSet = Collections.newSetFromMap(new ConcurrentHashMap<String,Boolean>());

	static int hCount;
	private static MolHelper2 initialize(String formula, String fragments) {
		MolHelper2 mol;
		try {
			mol = new MolHelper2();

			System.out.println("PMG7: Fork-join processing of "+formula+ " started:");
			hCount = mol.initialize(formula, fragments);
			for (IAtom atom:mol.acontainer.atoms()) System.out.print(atom.getSymbol());
			System.out.println();
			return mol;
		} catch (IOException e) {
			System.err.println("Could not read the fragments.");
			e.printStackTrace();
		} catch (CDKException e) {
			e.printStackTrace();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return null;	// upon unsuccessful initialization
	}

	public static void main(String[] args){
		// parse the command-line arguments
		String formula = null;
		String fragments = null;
		boolean wFile = false;
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-mf")){
				formula = args[++i];
			}
			else if(args[i].equals("-o")){
				String out = args[++i];
				try {
					outFile = new BufferedWriter(new FileWriter(out));
				} catch (IOException e) {
					System.err.println("Could not open "+out+" for writing.");
					e.printStackTrace();
					System.exit(10);
				}
				wFile  = true;
			}
			else if(args[i].equals("-fr")){
				System.err.println("No fragments currently supported in the parallel MG.");
			}
		}
		
		if (formula == null) {
			help();
			System.exit(0);
		}
		
		MolHelper2 mol = initialize(formula, fragments);
		
		// do the real processing
		long before = System.currentTimeMillis();
		long finalCount = forkJoinPool.invoke(new Generator(mol, wFile, hCount));
		long after = System.currentTimeMillis();		

		// Report the number of generated molecules
		System.out.println("molecules " + finalCount);
		System.out.println("Duration: " + (after - before) + " milliseconds\n");
	}


	private static void help() {
		// TODO Auto-generated method stub
		System.out.println("You can use the following options.");
		System.out.println("Specifying a formula is obligatory.");
		System.out.println("\t-mf \tA formula - the elemental composition");
		System.out.println("\t-o  \tThe name of the output file");
	}
}
