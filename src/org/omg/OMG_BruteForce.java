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
public class OMG_BruteForce{
	/**Output File containing the list of graph. */
	BufferedWriter outFile;
	HashMap<String, Byte> globalmap;
	int mol_counter;
	private int nH;
	private static boolean wfile = false;
	private SaturationChecker satCheck;
	public static void main(String[] args) throws IOException{

		OMG_BruteForce gen = new OMG_BruteForce();
		String formula = "C4H10";
		String fragments = null;
		String out = "default_out.sdf";


		for(int i = 0; i < args.length; i++){
			if(args[i].equals("-mf")){
				formula = args[i+1];
			}
			else if(args[i].equals("-o")){
				out = args[i+1];
				wfile = true;
			}
			else if(args[i].equals("-fr")){
				fragments = args[i+1];
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
	public OMG_BruteForce(){

	}
	void initializeMolecule(String formula, String fragments, String output) throws CDKException, FileNotFoundException, CloneNotSupportedException {
		long before = System.currentTimeMillis();
		System.out.println("OMG: Sequntial processing of "+formula+ " started (using nauty as canonizer).");
		System.out.print("Current atom order is: ");

		IAtomContainer acontainer;
		while(true) {
		acontainer = MolecularFormulaManipulator.getAtomContainer(
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
		for (IAtom atom:acontainer.atoms()) System.out.print(atom.getSymbol());
		System.out.println();
		if (formula.equals("C4H7NO3")) {
			if (!acontainer.getAtom(0).getSymbol().equals("C")) continue;
			if (!acontainer.getAtom(7).getSymbol().equals("N")) continue;
		}
		break;
		}

		mol_counter = 0;
		try {
			if (wfile){
				outFile = new BufferedWriter(new FileWriter(output));
//				for (IAtom atom:acontainer.atoms()) outFile.write(atom.getSymbol());
//				outFile.write("\n");
			}

			satCheck = new SaturationChecker();
			globalmap = new HashMap<String, Byte>();
			generateMol(acontainer, null);

			long after = System.currentTimeMillis();

			outFile.close();
			System.out.println("molecules " + mol_counter);
			System.out.println("Duration: " + (after - before) + " miliseconds\n");
		}  catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}






	private void generateMol(IAtomContainer acontainer, String canstr2) throws CloneNotSupportedException, CDKException, IOException {
		/*We check that the atoms are connected in the same molecules
		 * this is, they are not separated fragments*/
		/*We add hydrogens in order to check if the molecule is saturated.
		 * We will accept the molecule if the number of hydrogens necessary to saturate 
		 * is the same as the hydrogens in the original formula*/
		IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();
		
		boolean isComplete ;
		CDKHydrogenAdder hAdder;
		try {
			for (IAtom atom : acprotonate.atoms()) {
				IAtomType type = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder()).findMatchingAtomType(acprotonate, atom);

				AtomTypeManipulator.configure(atom, type);
			}
			hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
			hAdder.addImplicitHydrogens(acprotonate);
			isComplete = satCheck.isSaturated(acprotonate)&&(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH);
		} catch (IllegalArgumentException iae){
			isComplete = false;
//			System.err.println("incomplete at: "+mol_counter);
		}
		if(isComplete){
			if(ConnectivityChecker.partitionIntoMolecules(acontainer).getAtomContainerCount() == 1){

					// Not needed if we know there are no duplicates, e.g., when no initial fragments are given
				if(!globalmap.containsKey(canstr2)){
					globalmap.put(canstr2, null);
					mol_counter++;
					System.out.println(mol_counter);
					if(wfile){
						StringWriter writer = new StringWriter();
						MDLV2000Writer mdlWriter = new MDLV2000Writer(writer);
						mdlWriter.write(acprotonate);
						outFile.write(writer.toString());
						outFile.write("> <Id>\n"+(mol_counter)+"\n\n> <can_string>\n"+canstr2+"\n\n$$$$\n");
					}
				}
			}	

			hAdder = null;
			acprotonate = null;
			acontainer = null;
			return;
		}
		else{
			if(!globalmap.containsKey(canstr2)){
				globalmap.put(canstr2, null);

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
					IAtomContainer canonM_ext = MolManipulator.getcanonical(m_ext);
	
			        
					String canstr =  MolManipulator.array2string(MolManipulator.mol2array(canonM_ext));
	
//					if(!map.containsKey(canstr)){
//						map.put(canstr, null);				
						generateMol(canonM_ext, canstr);	
//					}
				}
			}		
			return;				
		}

	}

	public int getFinalCount() {
		// TODO Auto-generated method stub
		return mol_counter;
	}
}
