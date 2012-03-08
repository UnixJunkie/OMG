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

import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

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
 * @author julio, Mahdi
 */
public class OMGPar{
	/**Output File containing the list of graph. */
	BufferedWriter outFile;
	HashMap<String, Byte> globalmap;
	AtomicInteger mol_counter;
	private int nH;
	private static boolean wfile = false;
	private SaturationChecker satCheck;
//	private final static int threadCount = 20;
//	private ExecutorService threadPool;
	
    private static final int poolSize = 20;
    
//    int maxPoolSize = 20;
// 
//    long keepAliveTime = 100;
 
    ExecutorService threadPool = null;
    
 
//    final ArrayBlockingQueue<Runnable> queue = new ArrayBlockingQueue<Runnable>(500);
 
 
	public OMGPar() { 
		threadPool = Executors.newFixedThreadPool(poolSize);
//	    threadPool = new ThreadPoolExecutor(poolSize, maxPoolSize, keepAliveTime, TimeUnit.SECONDS, queue);
	}

	private int threadCount=0;
	private void generateMol(IAtomContainer atom, String canon) {
		try {
			threadPool.execute(new Generator(atom, canon));
			threadCount++;
			System.out.println("Threads successful: "+threadCount);
		} catch (RejectedExecutionException ree){
			System.exit(20);
		}
	}

	public static void main(String[] args) throws IOException{

		OMGPar gen = new OMGPar();
		String formula = "C4H10";
		String fragments = null;
		String out = "default_out.sdf";


		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-mf")){
				formula = args[i+1];
			}
			else if(args[i].equals("-o")){
				System.out.println("Outputing to file is currently disabled.");
//				out = args[i+1];
//				wfile = true;
			}
			else if(args[i].equals("-fr")){
				System.out.println("Fragments not supported for the moment.");
				// fragments = args[i+1];
			}
		}
		try {
			gen.initializeMolecule(formula,fragments, out);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	void initializeMolecule(String formula, String fragments, String output) throws CDKException, FileNotFoundException, CloneNotSupportedException {
		long before = System.currentTimeMillis();

		IAtomContainer acontainer = MolecularFormulaManipulator.getAtomContainer(
				MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));

		nH = 0;
		List<IAtom> listcont = new ArrayList<IAtom>();
		int atom_counter = 0;
		for(IAtom atom: acontainer.atoms()){
			if(atom.getSymbol().equals("H")){
				nH++;
				listcont.add(atom);
			}
			else{

				atom.setID("a"+atom_counter);
				atom.setFlag(1, false);
				atom_counter++;
			}
		}
		for(IAtom atom: listcont){
			acontainer.removeAtom(atom);
		}
		if(fragments != null){
			InputStream ins = new BufferedInputStream(new FileInputStream(fragments));
			MDLV2000Reader reader = new MDLV2000Reader(ins);
			ChemFile fileContents = (ChemFile)reader.read(new ChemFile());		        
			IChemSequence sequence = fileContents.getChemSequence(0);

			for (int i=0; i<sequence.getChemModelCount(); i++) {
				IAtomContainer frag = sequence.getChemModel(i).getMoleculeSet().getAtomContainer(0);

				try {
					acontainer = MolManipulator.buildFromFragment(acontainer,frag);
				} catch (CloneNotSupportedException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}

			acontainer = MolManipulator.getcanonical(acontainer);
		}		




		mol_counter  = new AtomicInteger(0);
		try {
			outFile = new BufferedWriter(new FileWriter(output));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			satCheck = new SaturationChecker();

			//			}
			globalmap = new HashMap<String, Byte>();

			generateMol(acontainer, null);
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		long after = System.currentTimeMillis();
		try {
			outFile.close();
			System.out.println(formula);
			System.out.println("molecules " + mol_counter);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.out.println("Duration: " + (after - before) + " miliseconds\n");
		
//		threadPool.shutdown(); TODO: where to put?
	}





	class Generator implements Runnable {
		IAtomContainer acontainer;
		String canstr2;
		
		public Generator(IAtomContainer acontainer, String canstr2) {
			super();
			this.acontainer = acontainer;
			this.canstr2 = canstr2;
		}

		// previously known as generateMol()
		public void run() { // throws CloneNotSupportedException, CDKException, IOException {
			/*We check that the atoms are connected in the same molecules
			 * this is, they are not separated fragments*/
			/*We add hydrogens in order to check if the molecule is saturated.
			 * We will accept the molecule if the number of hydrogens necessary to saturate 
			 * is the same as the hydrogens in the original formula*/
			try{
				IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();
		
				for (IAtom atom : acprotonate.atoms()) {
					IAtomType type = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder()).findMatchingAtomType(acprotonate, atom);
					
					AtomTypeManipulator.configure(atom, type);
				}
				CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
				hAdder.addImplicitHydrogens(acprotonate);
				// TODO: Is it safe to use satCheck concurrently?
				if(satCheck.isSaturated(acprotonate)&&(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH)){
					if(ConnectivityChecker.partitionIntoMolecules(acontainer).getAtomContainerCount() == 1){
						mol_counter.getAndIncrement();
	
						// Not needed if we know there are no duplicates, e.g., when no initial fragments are given
	//					if(!globalmap.containsKey(canstr2)){
	//						globalmap.put(canstr2, null);
	//						if(wfile){
	//							StringWriter writer = new StringWriter();
	//							MDLV2000Writer mdlWriter = new MDLV2000Writer(writer);
	//							mdlWriter.write(acprotonate);
	//							outFile.write(writer.toString());
	//							outFile.write("> <Id>\n"+(mol_counter+1)+"\n\n> <can_string>\n"+canstr2+"\n\n$$$$\n");
	//						}
	//					}					}
					}	
		
					hAdder = null;
					acprotonate = null;
					acontainer = null;
					return;
				}
				else{
		
					ArrayList<int[]> extBondlist = MolManipulator.extendMol(acontainer);
		
					HashMap<String, Byte> map = new HashMap<String, Byte>();
					
					for(int[] bond : extBondlist){
						IAtomContainer m_ext = (IAtomContainer) acontainer.clone();
		
						//Only one object bond is used, and recycled at every time we use it
						IBond bondAdd = m_ext.getBond(m_ext.getAtom(bond[0]), m_ext.getAtom(bond[1]));
						if(bondAdd == null){					
							m_ext.addBond(bond[0], bond[1], IBond.Order.SINGLE);
						}
						else if(bondAdd.getOrder() == IBond.Order.SINGLE){
							m_ext.getBond(m_ext.getAtom(bond[0]), m_ext.getAtom(bond[1])).setOrder(IBond.Order.DOUBLE);
						}
						else if(bondAdd.getOrder() == IBond.Order.DOUBLE){
							m_ext.getBond(m_ext.getAtom(bond[0]), m_ext.getAtom(bond[1])).setOrder(IBond.Order.TRIPLE);
						}
						else if(bondAdd.getOrder() == IBond.Order.TRIPLE){
							m_ext.getBond(m_ext.getAtom(bond[0]), m_ext.getAtom(bond[1])).setOrder(IBond.Order.QUADRUPLE);					
						}
		
						// end add bond
						IAtomContainer canonM_ext = m_ext; // TODO: MolManipulator.getcanonical(m_ext);
		
				        
						String canstr =  MolManipulator.array2string(MolManipulator.mol2array(canonM_ext));
		
						if(!map.containsKey(canstr)){
							map.put(canstr, null);				
							////					check for canonical augmentation
							boolean lastNotInFrag = false;
							String lastBondID[] = new String[2];
							int i = 0;
							while((!lastNotInFrag)&&(canonM_ext.getBondCount() > i)){
								i++;
								lastBondID[0] = canonM_ext.getBond(canonM_ext.getBondCount()-i).getAtom(0).getID();
								lastBondID[1] = canonM_ext.getBond(canonM_ext.getBondCount()-i).getAtom(1).getID();
								//					we remove select the last bond if it was not in the fragment, or if it was in the fragment
								//					but its desgree augmented
								bondAdd = canonM_ext.getBond(canonM_ext.getBondCount()-i);
								if(!(bondAdd.getProperty("BondINfrag")!= null)){
									lastNotInFrag = true;
								} 
								else if((bondAdd.getProperty("BondINfrag") == IBond.Order.SINGLE)&&
										(bondAdd.getOrder() == IBond.Order.DOUBLE)){
									lastNotInFrag = true;
								}
								else if((bondAdd.getProperty("BondINfrag") == IBond.Order.DOUBLE)&&
										(bondAdd.getOrder() == IBond.Order.TRIPLE)){
									lastNotInFrag = true;
								}
								else if((bondAdd.getProperty("BondINfrag") == IBond.Order.TRIPLE)&&
										(bondAdd.getOrder() == IBond.Order.QUADRUPLE)){
									lastNotInFrag = true;
								}	
							}
		
							int[] lastBond = new int[2];
							for(IAtom atom : m_ext.atoms()){
								if(atom.getID().equals(lastBondID[0]))
									lastBond[0]=m_ext.getAtomNumber(atom);
								if(atom.getID().equals(lastBondID[1]))
									lastBond[1] = m_ext.getAtomNumber(atom);						
							}
							bondAdd = m_ext.getBond(m_ext.getAtom(lastBond[0]),m_ext.getAtom(lastBond[1]));
							if(bondAdd.getOrder() == IBond.Order.SINGLE){
								m_ext.removeBond(bondAdd);
							}
							else if(bondAdd.getOrder() == IBond.Order.DOUBLE){
								bondAdd.setOrder(IBond.Order.SINGLE);
							}
							else if(bondAdd.getOrder() == IBond.Order.TRIPLE){
								bondAdd.setOrder(IBond.Order.DOUBLE);
							}
							else if(bondAdd.getOrder() == IBond.Order.QUADRUPLE){
							}

							if(MolManipulator.aresame(acontainer, m_ext)||(acontainer.getBondCount()==0)){
//							if(MolManipulator.aresame(acontainer, MolManipulator.getcanonical(m_ext))||(acontainer.getBondCount()==0)){
								// generate a parallel task instead of a recursive call
								generateMol(canonM_ext, canstr);	
	
							}	
						}
					}		
			//		return;				
				}
				
			} catch (CloneNotSupportedException e) {
				e.printStackTrace();
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
		}

	}

	public int getFinalCount() {
		// TODO check to make sure the count is finished!
		return mol_counter.get();
	}
}
