package org.omg;

import java.io.BufferedInputStream;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;

import javax.management.RuntimeErrorException;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
//import org.openscience.cdk.atomtype.GeneratorAtomTypeMatcher;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
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
public class MolHelper2 {
	IAtomContainer acontainer;
	int atomCount=0;
	String canString="";
	private static Map<String, List<Integer>> valenceTable; 
//	private static GeneratorAtomTypeMatcher matcher;
	
	int [] rep;
	private IAtomContainer acprotonate;
	
	public MolHelper2() {
		AtomTypeFactory factory = AtomTypeFactory.getInstance("org/openscience/cdk/dict/data/cdk-atom-types.owl", 
				new ChemObject().getBuilder()
				);
//		IAtomType[] types = factory.getAllAtomTypes();
//		List<IAtomType> typeList = new ArrayList<IAtomType>();
//		for (IAtomType type : types) {
//			typeList.add(type);
//		}
//		matcher = new GeneratorAtomTypeMatcher(typeList);
	}


	public MolHelper2(IAtomContainer acontainer, int[] orbit, String molString) {
		this.acontainer = acontainer;
		atomCount = acontainer.getAtomCount();
		this.rep = orbit.clone();
		this.canString = molString;
//		System.out.print("Mol Orbit:");
//		for (int i:rep) System.out.print(" "+rep[i]);
//		System.out.println();
	}


	public int initialize (String formula) throws CloneNotSupportedException{
		acontainer = MolecularFormulaManipulator.getAtomContainer(
				MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));

		atomCount = 0;
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
		rep = new int[atomCount];
		int r = -1;
		String prev="";
		for (int i=0; i<atomCount; i++){
			String symbol = acontainer.getAtom(i).getSymbol();
			if (!symbol.equals(prev)) {
				prev = symbol;
				r=i;
			}
			rep[i] = r;	// initially all atoms of a type are symmetric
		}
		canString = molString(acontainer);
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
		try{
			acprotonate = (IAtomContainer) acontainer.clone();

			for (IAtom atom : acprotonate.atoms()) {
				IAtomType type = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder()).findMatchingAtomType(acprotonate, atom);

				AtomTypeManipulator.configure(atom, type);
			}
			CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
			hAdder.addImplicitHydrogens(acprotonate);
			return (satCheck.isSaturated(acprotonate)&&(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH));
		} catch (IllegalArgumentException iae) {
			return false;
		}
	}
	

	private int openings = 0;
	private int [] bondCounts;
	private int[] perm1;
	private boolean addUpOpenings(int atomNum, int nH) {
		if (atomCount == atomNum) 
			return openings == nH;
		else {
			IAtom atom = acontainer.getAtom(atomNum);
			String symbol = atom.getSymbol();
			int bondSum = 0;
			for (IBond b:acontainer.getConnectedBondsList(atom))
				bondSum += b.getOrder().ordinal()+1;
			for (Integer valence:valenceTable.get(symbol)) {
				if (valence >= bondSum) {
					openings += (valence - bondSum);
					if (openings <= nH && addUpOpenings(atomNum+1, nH)) {
						return true;
					}
					openings -= (valence - bondSum);
				}
			}
		}
		return false;
	}
	
	public boolean isComplete(int nH) {
		return addUpOpenings (0, nH);
	}
	
	
	/**
	 * Checks if all atoms in the molecule are connected.
	 * @return
	 */
	public boolean isConnected() {
		return (ConnectivityChecker.partitionIntoMolecules(acontainer).getAtomContainerCount() == 1);
	}
	
	/**
	 * Must be called only on "complete" molecules, i.e., after a call has been made to isComplete() method.
	 * @param outFile
	 * @param mol_counter
	 * @throws CDKException
	 * @throws CloneNotSupportedException 
	 */
	public synchronized boolean writeTo(BufferedWriter outFile, long mol_counter) throws CDKException, CloneNotSupportedException{
		StringWriter writer = new StringWriter();
		MDLV2000Writer mdlWriter = new MDLV2000Writer(writer);
		try {
			acprotonate = (IAtomContainer) acontainer.clone();

			for (IAtom atom : acprotonate.atoms()) {
				IAtomType type = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder()).findMatchingAtomType(acprotonate, atom);

				AtomTypeManipulator.configure(atom, type);
			}
			CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
			hAdder.addImplicitHydrogens(acprotonate);

			mdlWriter.write(acprotonate);
			writer.append("> <Id>\n"+(mol_counter)+"\n\n> <can_string>\n"+canString+"\n\n$$$$\n");
			outFile.write(writer.toString());
		} catch (IOException e) {
			System.err.println("Could not write molecule to output file.");
			e.printStackTrace();
		} catch (IllegalArgumentException ila) {
			return false;
		}
		return true;
	}
	
	/**
	 * Writes a molecule to a file in SDF format. This is more or less equivalent to writeTo method but it does not depend on CDK. 
	 * The good is that it does not have the restrictions of having a "valid" atom container and is expected to wrok faster.
	 * The bad side is that it has less features and cannot handle all complicated chemical structures. 
	 * @param outFile
	 * @param mol_counter
	 */
	public synchronized void writeMol(BufferedWriter outFile, long mol_counter) {
		StringWriter writer = new StringWriter();
		writer.write(String.format("%n  PMG%n%n%3d%3d  0  0  0  0  0  0  0  0  0%n",acontainer.getAtomCount(), acontainer.getBondCount()));
		for (IAtom atom : acontainer.atoms()) {
			writer.write(String.format("    0.0000    0.0000    0.0000 %-3s 0  0  0  0  0  0%n", atom.getSymbol()));
		}
		for (IBond bond : acontainer.bonds()){
			writer.write(String.format("%3d%3d%3d%3d%3d%3d%n",acontainer.getAtomNumber(bond.getAtom(0))+1, (acontainer.getAtomNumber(bond.getAtom(1))+1), orderNumber(bond.getOrder()), 0, 0, 0));
		}
		writer.write("M  END\n> <Id>\n"+mol_counter+"\n\n> <can_string>\n"+canString+"\n\n$$$$\n");
		try {
			outFile.write(writer.toString());
		} catch (IOException e) {
			System.err.println("Could not write molecule to output file.");
			e.printStackTrace();
		}
	}
	
	/**
	 * It checks all possible ways to add one bond to a molecule, while considering the canonical augmentation of the corresponding graph.
	 * canAug cannot be true if bigStep is true (this is why this method is made private)
	 * 
	 * @param canAug specifies whether or not to use canonical augmentation.
	 * @param singleStep specifies whether to add only single bonds in each steps or consider also double and triple bonds
	 * @return
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 */
	private ArrayList<MolHelper2> addBond(boolean canAug, boolean bigStep) throws CloneNotSupportedException, CDKException{
    	Set<String> visited = new HashSet<>();
		ArrayList<MolHelper2> extMolList = new ArrayList<MolHelper2>();
		int vCount = acontainer.getAtomCount();

		bondCounts = new int [atomCount];
		for (IBond b:acontainer.bonds()) {
			bondCounts[acontainer.getAtomNumber(b.getAtom(0))] += orderNumber(b.getOrder());
			bondCounts[acontainer.getAtomNumber(b.getAtom(1))] += orderNumber(b.getOrder());
		}
		
		
		// Note that the representative of an atom never has a bigger ID
		for (int left = 0; left < vCount; left++){
			IAtom leftAtom = acontainer.getAtom(left);
//			if (valenceTable.get(leftAtom.getSymbol()).get(0) == bondCounts[left]) continue;
//			if (left>0 && rep[left] <= rep[left-1]) continue;	// make sure each orbit is considered only once
			for (int right = left+1; right < vCount; right++){
				IAtom rightAtom = acontainer.getAtom(right);
//				if (valenceTable.get(rightAtom.getSymbol()).get(0) == bondCounts[right]) continue;
//				if (right>left+1 && rep[right] <= rep[right-1]) continue;	// make sure each orbit is considered only once
				// For the first iteration (in inner loop), we may consider the same orbit as "left"
				
				if (bigStep && acontainer.getBond(leftAtom, rightAtom) != null) continue;
				
				IAtomContainer copyMol = (IAtomContainer) acontainer.clone();
				do {
					if (!incBond(left, right, copyMol)) break;	

					// canonize 
					Graph graph = new Graph();
					int[] perm1 = graph.canonize(copyMol, true);	// ask for the automorphisms to be reported back
					IAtomContainer canExtMol = Graph.relabel(copyMol, perm1);
					int[] orbit = graph.orbitRep;
					String molString = molString(canExtMol);
					if (visited.add(molString) == false) break;	

					if (!canAug || acontainer.getBondCount()==0){ // no need to check canonical augmentation
						extMolList.add(new MolHelper2(canExtMol, orbit, molString)); 
						continue;
					}


					// remove the last bond and canonize again (to check for canonical augmentation)
					int maxBond = 0; 
					for (int p=1; p<perm1.length; p++) if (perm1[maxBond] < perm1[p]) maxBond = p;

					copyMol = (IAtomContainer) canExtMol.clone();
					Iterator<IBond> bonds = copyMol.bonds().iterator();
					IBond lastBond;
					do {
						lastBond = bonds.next();
					} while (Integer.parseInt(lastBond.getID()) < maxBond);
					decBond(lastBond, copyMol);	// remove a bond ....

					int[] perm = graph.canonize(copyMol, true);
					copyMol = Graph.relabel(copyMol, perm);
					String parentString = molString(copyMol);

					if (this.canString.equals(parentString)){
						//				if (aresame(acontainer, copyMol)){
						//				if (aresame(Graph.relabel(acontainer, perm1), copyMol)){	// Is it possible to avoid the second canonization?
						extMolList.add(new MolHelper2(canExtMol, orbit, molString)); 
					}
				} while (bigStep);
			}
		}
		return extMolList;
	}

	ArrayList<MolHelper2> addOneBond() throws CloneNotSupportedException, CDKException{
		return addBond(true, false);
	}

	ArrayList<MolHelper2> addOneBondNoCheck() throws CloneNotSupportedException, CDKException{
		return addBond(false, false);
	}

	ArrayList<MolHelper2> addBigBond() throws CloneNotSupportedException, CDKException {
		return addBond(false, true);
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
		return o.ordinal()+1;
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


	private boolean incBond(int left, int right, IAtomContainer m_ext) {
		IAtom leftAtom = m_ext.getAtom(left);
		IAtom rightAtom = m_ext.getAtom(right);
		
		if (isFull(leftAtom, m_ext))  return false;
		if (isFull(rightAtom, m_ext)) return false;

		//Only one object bond is used, and recycled at every time we use it
		IBond bondAdd = m_ext.getBond(leftAtom, rightAtom);
		if(bondAdd == null){					
			m_ext.addBond(left, right, IBond.Order.SINGLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.SINGLE){
			bondAdd.setOrder(IBond.Order.DOUBLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.DOUBLE){
			bondAdd.setOrder(IBond.Order.TRIPLE);
		}
		else if(bondAdd.getOrder() == IBond.Order.TRIPLE){
			return false;	// disallow quadruple bonds!
//			bondAdd.setOrder(IBond.Order.QUADRUPLE);					
		} 
		else if(bondAdd.getOrder() == IBond.Order.QUADRUPLE){
			return false;
		}
		return true;
	}


	/**
	 * @param atom
	 */
	private boolean isFull(IAtom atom, IAtomContainer m_ext) {
		int bondSum;
		bondSum = 0;
		for (IBond b:m_ext.getConnectedBondsList(atom))
			bondSum += b.getOrder().ordinal()+1;
		return (valenceTable.get(atom.getSymbol()).get(0) <= bondSum);
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
	static {

		// initialize the table
		valenceTable = new HashMap<>();
		// TODO: read atom symbols from CDK
		valenceTable.put("C", Arrays.asList(4) );
		valenceTable.put("N", Arrays.asList(5,3));	// Remember to put the biggest number first!
		valenceTable.put("O", Arrays.asList(2));
		valenceTable.put("S", Arrays.asList(6,4,2));
		valenceTable.put("P", Arrays.asList(5,3));
	}
}
