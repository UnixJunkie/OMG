package org.omg;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.management.RuntimeErrorException;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import fi.tkk.ics.jbliss.Graph;


/**
 * Changing a multi-graph to a single graph (that can be used by existing graph isomorphism checkers like: bliss or nauty
 * To do so, each "edge" in the multi-graph is converted to a vertex in the normal graph. All edges will have the same color
 * which will be different from the colors of other vertices.
 * 
 * @author mmajid
 *
 */
public class MolHelper {
	IAtomContainer acontainer;
	int atomCount=0;
	
	public MolHelper() {
	}


	public MolHelper(IAtomContainer acontainer) {
		this.acontainer = acontainer;
		atomCount = acontainer.getAtomCount();
	}


	public int initialize (String formula) throws CloneNotSupportedException{
		acontainer = MolecularFormulaManipulator.getAtomContainer(
				MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));

		// Count and remove the Hydrogens
		int nH = 0;
		List<IAtom> listcont = new ArrayList<IAtom>();
		for(IAtom atom: acontainer.atoms()){
			String symbol = atom.getSymbol();
			if(symbol.equals("H")){
				nH++;
				listcont.add(atom);
			} else {
				atom.setID(""+atomCount++);
				atom.setFlag(1, false);
			}
		}
		for(IAtom atom: listcont){
			acontainer.removeAtom(atom);
		}
		return nH;
	}
	
	public int initialize (String formula, String fragments) throws CloneNotSupportedException, CDKException, FileNotFoundException{
		throw new RuntimeErrorException(new Error("The fragments are not yet implemented."));
		/*
		// generate the atom as in the above constructor
		int nH = this.initialize(formula);

		// Assuming that fragments is not null, add the fragments to the atom structure and canonicalize
		InputStream ins = new BufferedInputStream(new FileInputStream(fragments));
		MDLV2000Reader reader = new MDLV2000Reader(ins);
		ChemFile fileContents = (ChemFile)reader.read(new ChemFile());		        
		IChemSequence sequence = fileContents.getChemSequence(0);
		
		for (int i=0; i<sequence.getChemModelCount(); i++) {
			IAtomContainer frag = sequence.getChemModel(i).getMoleculeSet().getAtomContainer(0);
			
			try {
				acontainer = MolManipulator.buildFromFragment(acontainer,frag);
			} catch (CloneNotSupportedException e1) {
				e1.printStackTrace();
			}
		}
		// TODO: make it canonical?
		
		acontainer = MolManipulator.getcanonical(acontainer);
		return nH;
		*/
	}
		
	/**
	 * Checks if a molecule is complete, including checking the number of hydrogens against nH, and saturation.
	 * @param satCheck
	 * @param nH
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	public boolean isComplete(SaturationChecker satCheck, int nH) throws CloneNotSupportedException, CDKException{
		IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();

		for (IAtom atom : acprotonate.atoms()) {
			IAtomType type = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder()).findMatchingAtomType(acprotonate, atom);
			
			AtomTypeManipulator.configure(atom, type);
		}
		CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
		hAdder.addImplicitHydrogens(acprotonate);
		return (satCheck.isSaturated(acprotonate)&&(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH));
	}
	
	
	/**
	 * Checks if all atoms in the molecule are connected.
	 * @return
	 */
	public boolean isConnected() {
		return (ConnectivityChecker.partitionIntoMolecules(acontainer).getAtomContainerCount() == 1);
	}
	
	
	
	/**
	 * It checks all possible ways to add one bond to a molecule, while considering the canonical augmentation of the corresponding graph.
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	ArrayList<MolHelper> addOneBond() throws CloneNotSupportedException, CDKException{
    	Set<String> visited = new HashSet<>();
		ArrayList<MolHelper> extMolList = new ArrayList<MolHelper>();
		
		ArrayList<int[]> extBondlist = MolManipulator.extendMol(acontainer);
		
		for(int[] bond : extBondlist){
			// add one bond and canonize 
			IAtomContainer copyMol = (IAtomContainer) acontainer.clone();
			incBond(bond[0], bond[1], copyMol);
			int[] perm1 = new Graph().canonize(copyMol);
			IAtomContainer canExtMol = Graph.relabel(copyMol, perm1);
			if (visited.add(molString(canExtMol)) == false) continue;	

			if (acontainer.getBondCount()==0){ // no need to check canonical augmentation
				extMolList.add(new MolHelper(canExtMol)); 
				continue;
			}

			// remove the last bond and canonize again (to check for canonical augmentation)
			int maxBond = 0; 
			for (int i=1; i<perm1.length; i++) if (perm1[maxBond] < perm1[i]) maxBond = i;
			
			copyMol = (IAtomContainer) canExtMol.clone();
			Iterator<IBond> bonds = copyMol.bonds().iterator();
			IBond lastBond;
			do {
				lastBond = bonds.next();
			} while (Integer.parseInt(lastBond.getID()) < maxBond);
			decBond(lastBond, copyMol);	// remove a bond ....
			int[] perm = new Graph().canonize(copyMol);
			copyMol = Graph.relabel(copyMol, perm);
			
			if (aresame(acontainer, copyMol)){
				extMolList.add(new MolHelper(canExtMol)); 
			}
		}
		return extMolList;
	}

	/**
	 * It checks all possible ways to add one bond to a molecule, while considering the canonical augmentation of the corresponding graph.
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	ArrayList<MolHelper> addOneBond2(int a) throws CloneNotSupportedException, CDKException{
		ArrayList<MolHelper> extMolList = new ArrayList<MolHelper>();
		Set<String> visited = new HashSet<>();
//		int vCount = acontainer.getAtomCount();

    	for (int j = 0; j < a; j++){
    		IAtomContainer copyMol = (IAtomContainer) acontainer.clone();
    		while (incBond(a, j, copyMol)) {	
    			CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(copyMol.getBuilder());
    			IAtomType type1 = matcher.findMatchingAtomType(copyMol,copyMol.getAtom(a));
    			IAtomType type2 = matcher.findMatchingAtomType(copyMol,copyMol.getAtom(j));
    			if(type1 == null || type2 == null) continue;

    			// canonize 
    			int[] perm1 = new Graph().canonize(copyMol);
    			IAtomContainer canExtMol = Graph.relabel(copyMol, perm1);

    			if (visited.add(molString(canExtMol)) == false) continue;	

    			if (acontainer.getBondCount()==0){ // no need to check canonical augmentation
    				extMolList.add(new MolHelper(canExtMol)); 
    				continue;
    			}

    			// remove the last bond and canonize again (to check for canonical augmentation)
    			int maxBond = 0; 
    			for (int i=1; i<perm1.length; i++) if (perm1[maxBond] < perm1[i]) maxBond = i;

    			copyMol = (IAtomContainer) canExtMol.clone();
    			Iterator<IBond> bonds = copyMol.bonds().iterator();
    			IBond lastBond;
    			do {
    				lastBond = bonds.next();
    			} while (Integer.parseInt(lastBond.getID()) < maxBond);
    			decBond(lastBond, copyMol);	// remove a bond ....
    			int[] perm = new Graph().canonize(copyMol);
    			copyMol = Graph.relabel(copyMol, perm);

    			if (aresame(acontainer, copyMol)){
    				extMolList.add(new MolHelper(canExtMol)); 
    			}
    		}
		}
		return extMolList;
	}


	
	/**
	 * 
	 * @param ac
	 * @return
	 */
	private static String molString(IAtomContainer ac) {
		StringBuffer s = new StringBuffer();
		int aCount = ac.getAtomCount();
		int[][] ar = molArray(ac);	
		// Turn the adjacency graph into a string by taking only the lower half.
		for (int i=0; i<aCount; i++)
			for (int j=0; j<i; j++){
				s.append(ar[i][j]);
			}

		return s.toString();
	}

	public static int orderNumber(Order o){
		if (o == null)
			return 0;
		if(o == IBond.Order.SINGLE){
			return 1;
		}
		if(o == IBond.Order.DOUBLE){
			return 2;
		}
		if(o == IBond.Order.TRIPLE){
			return 3;
		}
		return 4;
	}

	/**
	 * @param ac
	 * @param aCount
	 * @return
	 * @throws NumberFormatException
	 */
	private static int[][] molArray(IAtomContainer ac)
			throws NumberFormatException {
		int aCount = ac.getAtomCount();
		int[][] ar = new int[aCount][aCount];
		for(IBond bond : ac.bonds()){
			//we read only the bonds and store the bond degree as ar[left][right] and ar[right][left]
			int left = Integer.parseInt(bond.getAtom(0).getID());
			int right = Integer.parseInt(bond.getAtom(1).getID());
			ar[right][left] = ar[left][right] = orderNumber(bond.getOrder());
		}
		return ar;
	}

	public static boolean aresame(IAtomContainer ac1, IAtomContainer ac2) {
		return Arrays.deepEquals(molArray(ac1),molArray(ac2));
	}

	private boolean incBond(int leftAtom, int rightAtom, IAtomContainer m_ext) {
		//Only one object bond is used, and recycled at every time we use it
		IBond bondAdd = m_ext.getBond(m_ext.getAtom(leftAtom), m_ext.getAtom(rightAtom));
		if(bondAdd == null){					
			m_ext.addBond(leftAtom, rightAtom, IBond.Order.SINGLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.SINGLE){
			bondAdd.setOrder(IBond.Order.DOUBLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.DOUBLE){
			bondAdd.setOrder(IBond.Order.TRIPLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.TRIPLE){
			bondAdd.setOrder(IBond.Order.QUADRUPLE);					
		} 
		else if(bondAdd.getOrder() == IBond.Order.QUADRUPLE){
			return false;
		}
		return true;
	}
	
	private boolean decBond(IBond bondDel, IAtomContainer m_ext) {
		// = m_ext.getBond(m_ext.getBondCount()-1);
//		for (IBond bond : m_ext.bonds()) {
//			if (Integer.parseInt(bond.getAtom(0).getID()) < Integer.parseInt(bondDel.getAtom(0).getID()) && 
//				Integer.parseInt(bond.getAtom(1).getID()) < Integer.parseInt(bondDel.getAtom(1).getID())   )
//				bondDel = bond;
//		}
//		IBond bondDel = m_ext.getBond(m_ext.getAtom(leftAtom), m_ext.getAtom(rightAtom));
		if (bondDel == null) return false;
		if(bondDel.getOrder() == IBond.Order.SINGLE){
			m_ext.removeBond(bondDel);
		}
		else if(bondDel.getOrder() == IBond.Order.DOUBLE){
			bondDel.setOrder(IBond.Order.SINGLE);
		}
		else if(bondDel.getOrder() == IBond.Order.TRIPLE){
			bondDel.setOrder(IBond.Order.DOUBLE);
		}
		else if(bondDel.getOrder() == IBond.Order.QUADRUPLE){
			bondDel.setOrder(IBond.Order.TRIPLE);
		}
		return true;
	}
}
