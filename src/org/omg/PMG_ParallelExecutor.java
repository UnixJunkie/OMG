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
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
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
public class PMG_ParallelExecutor{
	/**Output File containing the list of graph. */
//	BufferedWriter outFile;

	AtomicLong mol_counter;
	private SaturationChecker satCheck = new SaturationChecker();
	
	private static final int executorCount = 6;
	ExecutorService executor[];
	AtomicLong startedTasks;

	private int nH;
//	private static boolean wfile = false;
	
	public PMG_ParallelExecutor(){
		mol_counter = new AtomicLong(0);
		startedTasks = new AtomicLong(0);
		executor = new ExecutorService[executorCount];
		for (int i=0; i<executorCount; i++){
			executor[i] = Executors.newSingleThreadExecutor();
		}
	}

	public class AtomicModInteger {
	    private AtomicInteger value;
	    private final int max;

	    public AtomicModInteger(int start, int max) {
	        this.value = new AtomicInteger(start);
	        this.max = max;
	    }

	    public int get() {
	        return value.get();
	    }

	    /* Simple modification of AtomicInteger.incrementAndGet() */
	    public int incrementAndGet() {
	        for (;;) {
	            int current = get();
	            int next = (current + 1) % max;
	            if (value.compareAndSet(current, next))
	                return next;
	        }
	    }
	}

	AtomicModInteger loadBalance = new AtomicModInteger(0, executorCount);
	private void generateTaskMol(MoleculeGraph molecule) {
		startedTasks.getAndIncrement();
		executor[loadBalance.incrementAndGet()].execute(new Generator(molecule));
	}

	public static void main(String[] args) throws IOException{

		// parse the command-line arguments
		String formula = "C4H10";
		String fragments = null;
		String out = "default_out.sdf";
		
		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-mf")){
				formula = args[i+1];
			}
			else if(args[i].equals("-o")){
				System.err.println("No output file currently supported in the parallel MG.");
			}
			else if(args[i].equals("-fr")){
				System.err.println("No fragments currently supported in the parallel MG.");
			}
		}

//		try {
//		outFile = new BufferedWriter(new FileWriter(output));
//	} catch (IOException e) {
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	}
		// start the process of generating structures
		long before = System.currentTimeMillis();
		System.out.println(formula);
		
		MoleculeGraph mol;
		PMG_ParallelExecutor gen = new PMG_ParallelExecutor();
		try {
			mol = new MoleculeGraph();
			if (fragments == null)
				gen.nH = mol.initialize(formula);
			else 
				gen.nH = mol.initialize(formula, fragments);
			
			gen.generateTaskMol(mol);
		} catch (CDKException e) {
			e.printStackTrace();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
//	try {
//		outFile.close();
//	} catch (IOException e) {
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	}
		
		System.out.println("molecules " + gen.getFinalCount());

		long after = System.currentTimeMillis();		
		System.out.println("Duration: " + (after - before) + " miliseconds\n");
		
		gen.shutdown();
	}

	
	private void shutdown() {
		for (ExecutorService e : executor){
			e.shutdown();
		}
	}


	private class Generator implements Runnable {
		MoleculeGraph mol;
		
		public Generator(MoleculeGraph mol) {
			super();
			this.mol = mol;
		}

		@Override
		public void run() {
			try {
				if (mol.isComplete(satCheck, nH)){
					if (mol.isConnected()) {
						mol_counter.incrementAndGet();
//						if(wfile){
//							StringWriter writer = new StringWriter();
//							MDLV2000Writer mdlWriter = new MDLV2000Writer(writer);
//							mdlWriter.write(acprotonate);
//							outFile.write(writer.toString());
//							outFile.write("> <Id>\n"+(mol_counter+1)+"\n\n> <can_string>\n"+canstr2+"\n\n$$$$\n");
//						}
					}	
//				return;
				}
				else{
					// get all possible ways to add one bond to the molecule
					ArrayList<MoleculeGraph> extMolList = mol.addOneBond();
					
					// recursively process all extended molecules
					for (MoleculeGraph  molecule : extMolList) {
						generateTaskMol(molecule);	
					}
//				return;				
				}
			} catch (CloneNotSupportedException e){
				e.printStackTrace();
			} catch (CDKException e) {
				e.printStackTrace();
			}
			startedTasks.decrementAndGet();
		}
		
	}

	public long getFinalCount() {
		// TODO make sure the computation is finished before returning the final count
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
