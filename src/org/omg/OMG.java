package org.omg;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

/**
 * Open Molecule Generation
 * The main class collecting parameters and setting global objects
 * 
 * @author Julio Peironcely
 */
public class OMG{
	SDFWriter outFile;
	HashMap<String, Byte> globalmap;
	public static int mol_counter;
	int maxOpenings;
	private int nH;
	private static boolean wfile = false;
	private static String fragments = null;
	private CDKHydrogenAdder hAdder;
	private CDKAtomTypeMatcher matcher;
	private SaturationChecker satCheck;
	private static ConnectivityChecker conCheck;
	private static Map<String, Double> valencetable; 
	static {

		// initialize the table
		valencetable = new HashMap<String, Double>();
		// TODO: read atom symbols from CDK
		valencetable.put("C", new Double(4));
		valencetable.put("N", new Double(5));
		valencetable.put("O", new Double(2));
		valencetable.put("S", new Double(6));
		valencetable.put("P", new Double(5));
		valencetable.put("F", new Double(1));
		valencetable.put("I", new Double(1));
		valencetable.put("Cl", new Double(1));
	}
	public static void main(String[] args) throws IOException{

		OMG gen = new OMG();
		String formula = null;
		String out = "default_out.sdf";
		
		if (args.length > 0) {
			
			for(int i = 0; i < args.length; i++){
				if(args[i].equals("-h")){
					System.out.println("OMG generates chemical structures");
					System.out.println("");
					System.out.println("Usage: java -jar OMG.jar -ec <elemental_composition> [-o <out_file.sdf>, -fr <in_fragments.sdf>]");
					System.out.println("");
					System.out.println("Required Parameters");
					System.out.println("-ec:  elemental composition of the molecules to be generated.");
					System.out.println("");
					System.out.println("Optional Parameters");
					System.out.println("-o:   SDF file where to store the molecules. ");
					System.out.println("-fr:  SDF file containing prescribed one or multiple substructures. In the case");
					System.out.println("         of multiple substructures, they have to be non-overlapping. ");
					System.out.println("");
					System.out.println("");
					System.out.println("Examples:");
					System.out.println("java -jar OMG.jar -ec C6H6");
					System.out.println("");
					System.out.println("java -jar OMG.jar -ec C6H6 -o out_C6H6.sdf");
					System.out.println("");
					System.out.println("java -jar OMG.jar -ec C2H5NO2 -fr fragment_CO2.sdf");
					System.out.println("");

					System.exit(1);
				}
				if(args[i].equals("-ec")){
					try {
						formula = args[i+1];
				    } catch (Exception e) {
				        System.err.println("No formula provided");
				        System.exit(1);
				    }
				}
				else if(args[i].equals("-o")){
					try {
						out = args[i+1];
						wfile = true;
				    } catch (Exception e) {
				        System.err.println("No output file provided");
				        System.exit(1);
				    }						
				}
				else if(args[i].equals("-fr")){
					try {
						fragments = args[i+1];
				    } catch (Exception e) {
				        System.err.println("No file with prescribed substructures provided");
				        System.exit(1);
				    }					
				}
			}
		}
		else{
			System.err.println("Provide at least an elemental composition. Type OMG.jar -h for help");
			System.exit(1);
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
	public OMG(){
		matcher = CDKAtomTypeMatcher.getInstance(SilentChemObjectBuilder.getInstance());
		hAdder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
	}
	void initializeMolecule(String formula, String fragments, String output) throws CDKException, FileNotFoundException, CloneNotSupportedException {
		long before = System.currentTimeMillis();
		IAtomContainer acontainer = null;
//		do{
		acontainer = MolecularFormulaManipulator.getAtomContainer(
				MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));
//		}while(!acontainer.getAtom(0).getSymbol().equals("C"));
		nH = 0;
		List<IAtom> listcont = new ArrayList<IAtom>();
		int atom_counter = 0;
		String symbol = null;
		for(IAtom atom: acontainer.atoms()){
			symbol = atom.getSymbol();
			if(symbol.equals("H")){
				nH++;
				maxOpenings--;
				listcont.add(atom);
			}
			else{
				maxOpenings += valencetable.get(symbol);
				atom.setID("a"+atom_counter);
				atom.setFlag(1, false);
				atom_counter++;
				System.out.print(atom.getSymbol());
			}
		}
		System.out.println();
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
			acontainer = MolManipulator.getcanonical(acontainer, fragments);
		}		

		mol_counter = 0;
		try {
			outFile = new SDFWriter(new FileWriter(output));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {

			satCheck = new SaturationChecker();
			conCheck = new ConnectivityChecker();
			globalmap = new HashMap<String, Byte>();
			generateMol(acontainer, null, false);
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
	}

	private void generateMol(IAtomContainer acontainer, String canstr2, boolean extraAtoms) throws CloneNotSupportedException, CDKException, IOException {
		/*We check that the atoms are connected in the same molecules
		 * this is, they are not separated fragments*/
		/*We add hydrogens in order to check if the molecule is saturated.
		 * We will accept the molecule if the number of hydrogens necessary to saturate 
		 * is the same as the hydrogens in the original formula*/
		IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();
		int bondCount = 0;
		boolean isComplete = false;
		try {
			for (IAtom atom : acprotonate.atoms()) {
				IAtomType type = matcher.findMatchingAtomType(acprotonate, atom);
				AtomTypeManipulator.configure(atom, type);
			}
			hAdder.addImplicitHydrogens(acprotonate);
			if(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH){
				isComplete = satCheck.isSaturated(acprotonate);
			}
		} catch (IllegalArgumentException iae){
			isComplete = false;
		}
		if(isComplete && !extraAtoms){
			if(conCheck.partitionIntoMolecules(acprotonate).getAtomContainerCount() == 1){
				if(fragments != null){
					if(!globalmap.containsKey(canstr2)){
						globalmap.put(canstr2, null);
						mol_counter++;
						if(wfile){
							acprotonate.setProperty("Id",mol_counter);
							acprotonate.setProperty("can_string",canstr2);
							outFile.write(acprotonate);
						}		
					}				
				}else{
					mol_counter++;
					if(wfile){
						acprotonate.setProperty("Id",mol_counter);
						acprotonate.setProperty("can_string",canstr2);
						outFile.write(acprotonate);
					}
				}
			}	
			for (IBond b:acprotonate.bonds()) {
				bondCount += b.getOrder().ordinal()+1;
			}
			if (maxOpenings>bondCount*2) {
				generateMol(acontainer, canstr2, true);				
			}
		}
		else{
			ArrayList<int[]> extBondlist = MolManipulator.extendMol(acontainer);
			HashMap<String, Byte> map = new HashMap<String, Byte>();
			IBond bondAdd = null;
			IAtomContainer m_ext = null;
			Order ord = null;
			for(int[] bond : extBondlist){
				ord = null;
				m_ext = (IAtomContainer) acontainer.clone();

				//Only one object bond is used, and recycled at every time we use it
				bondAdd = m_ext.getBond(m_ext.getAtom(bond[0]), m_ext.getAtom(bond[1]));
				if(bondAdd == null){					
					m_ext.addBond(bond[0], bond[1], IBond.Order.SINGLE);
				}
				else{
					ord = bondAdd.getOrder();
					if(ord == IBond.Order.SINGLE){
						m_ext.getBond(m_ext.getAtom(bond[0]), m_ext.getAtom(bond[1])).setOrder(IBond.Order.DOUBLE);
					}
					else if(ord == IBond.Order.DOUBLE){
						m_ext.getBond(m_ext.getAtom(bond[0]), m_ext.getAtom(bond[1])).setOrder(IBond.Order.TRIPLE);
					}
				}
				// end add bond
				IAtomContainer canonM_ext = MolManipulator.getcanonical(m_ext,fragments);
			
				String canstr =  MolManipulator.array2string(MolManipulator.mol2array(canonM_ext));

				if(!map.containsKey(canstr)){
					map.put(canstr, null);				
					////					check for canonical augmentation
					boolean lastNotInFrag = false;
					String lastBondID[] = new String[2];
					int i = 0;
					int[] lastBond = new int[2];
					while((!lastNotInFrag)&&(canonM_ext.getBondCount() > i)){
						i++;
						lastBondID[0] = canonM_ext.getBond(canonM_ext.getBondCount()-i).getAtom(0).getID();
						lastBondID[1] = canonM_ext.getBond(canonM_ext.getBondCount()-i).getAtom(1).getID();
						//					we remove the last bond if it was not in the fragment, or if it was in the fragment
						//					but its degree augmented
						bondAdd = canonM_ext.getBond(canonM_ext.getBondCount()-i);
						if(!(bondAdd.getProperty("BondINfrag")!= null)){
							lastNotInFrag = true;
						} 
						else if((bondAdd.getProperty("BondINfrag") == IBond.Order.SINGLE)&&
								(bondAdd.getOrder() == IBond.Order.DOUBLE)){
							lastNotInFrag = false;
							canonM_ext.getBond(canonM_ext.getBondCount()-i).setProperty("BondINfrag", ""+canonM_ext.getBond(canonM_ext.getBondCount()-i).getOrder());
						}
						else if((bondAdd.getProperty("BondINfrag") == IBond.Order.DOUBLE)&&
								(bondAdd.getOrder() == IBond.Order.TRIPLE)){
							lastNotInFrag = false;
							canonM_ext.getBond(canonM_ext.getBondCount()-i).setProperty("BondINfrag", ""+canonM_ext.getBond(canonM_ext.getBondCount()-i).getOrder());

						}
					}
					lastBond[0] = m_ext.getAtomNumber(AtomContainerManipulator.getAtomById(m_ext, lastBondID[0]));
					lastBond[1] = m_ext.getAtomNumber(AtomContainerManipulator.getAtomById(m_ext, lastBondID[1]));
					bondAdd = m_ext.getBond(m_ext.getAtom(lastBond[0]),m_ext.getAtom(lastBond[1]));
					ord = bondAdd.getOrder();
					if(ord == IBond.Order.SINGLE){
						m_ext.removeBond(bondAdd);
					}
					else if(ord == IBond.Order.DOUBLE){
						bondAdd.setOrder(IBond.Order.SINGLE);
					}
					else if(ord == IBond.Order.TRIPLE){
						bondAdd.setOrder(IBond.Order.DOUBLE);
					}				
					if(MolManipulator.aresame(acontainer, MolManipulator.getcanonical(m_ext, fragments))||(acontainer.getBondCount()==0)){			
						generateMol(canonM_ext, canstr, false);	
					}	
				}
			}		
		}
	}

	public int getFinalCount() {
		// TODO Auto-generated method stub
		return mol_counter;
	}
}
