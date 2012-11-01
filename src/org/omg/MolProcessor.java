package org.omg;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

import org.omg.tools.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.ICDKObject;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.SaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import fi.tkk.ics.jbliss.Graph;

public class MolProcessor implements Runnable{
	static final int SEM_CAN = 0;
	static final int MIN_CAN = 1;
	static final int CAN_AUG = 2;
	static final int BRT_FRC = 3;
	final int method;
	final boolean hashMap;	// otherwise, minimality check (with semiCan)
	final boolean cdkCheck;
	
	public final Atom[] atoms;
	final int nH;
	int maxOpenings;
	final String canString;
	final Graph graph;
	final int[][] adjacency;
	final static AtomicLong duplicate = new AtomicLong(0);
	final int startLeft, startRight;
	private final IAtomContainer acontainer;
	
	private static final Set<String> molSet = Collections.synchronizedSet(new HashSet<String>());
	
	
	public MolProcessor(final ArrayList<String> atomSymbols, String formula,
			final int method, final boolean hm, final boolean cdk){
		this.method = method;
		this.hashMap = hm;
		this.cdkCheck = cdk;
		int nH=0;
		int maxOpenings=0;
		int atomCount=0;
		for (String symbol:atomSymbols){
			if (symbol.equals("H")) {
				nH++;
				maxOpenings--;
			}
		}
		atoms = new Atom[atomSymbols.size()-nH];
		this.nH = nH;
		for (String symbol:atomSymbols){
			if (symbol.equals("H")) continue;
			atoms[atomCount] = new Atom(symbol,atomCount);
			maxOpenings += atoms[atomCount].maxValence;
			atomCount++;
		}
		adjacency = new int [atomCount][atomCount];
		graph = new Graph(atomCount);
		this.maxOpenings = maxOpenings;
		canString = "";
		startLeft=0;
		startRight=1;
		IAtomContainer lcontainer;
		do {	// this stupid loop is here only because CDK is not stable with respect to the order of the atoms it returns
			lcontainer = MolecularFormulaManipulator.getAtomContainer(
				MolecularFormulaManipulator.getMolecularFormula(formula, DefaultChemObjectBuilder.getInstance()));
		}while(!inRightOrder(atomSymbols, lcontainer));	// make sure the order is as we want
		acontainer = lcontainer;
	}

	private boolean inRightOrder(ArrayList<String> atomSymbols, IAtomContainer lcontainer) {
		int atom = 0;
		for (String symbol:atomSymbols){
			if (symbol.equals("H")) continue;	// skip hydrogens
			while (lcontainer.getAtom(atom).getSymbol().equals("H")) { lcontainer.removeAtom(atom); } // remove hydrogens
			String symbol2 = lcontainer.getAtom(atom++).getSymbol();
			if (!symbol2.equals(symbol)) return false;	// compare other atoms
		}
		return true;
	}

	/**
	 * Used with the semi-canonization method, which does not keep 
	 * a canonical form, but a semi-canonical form of the molecule.
	 * 
	 * @param atoms
	 * @param nH
	 * @param maxOpenings
	 * @param adjacency
	 * @param gr
	 */
	public MolProcessor(final Atom[] atoms, final int nH, final int maxOpenings, 
			            final int[][] adjacency, final Graph gr, 
			            final int stL, final int stR, final IAtomContainer acontainer,
			            final int method, final boolean hm, final boolean cdk) {
		this.method = method;
		this.hashMap = hm;
		this.cdkCheck = cdk;
		this.atoms = atoms;
		this.nH = nH;
		this.maxOpenings = maxOpenings;
		graph = gr; 
		this.adjacency = new int [atoms.length][atoms.length];
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++){
				this.adjacency[i][j] = adjacency[i][j];
			}
		this.canString = "";
		this.startLeft = stL;
		this.startRight = stR;
		this.acontainer = acontainer;
	}
		 	
	public MolProcessor(final Atom[] atoms, final int nH, final int maxOpenings, 
			            final int[][] adjacency, final Graph gr, final int[] canPerm,
			            final IAtomContainer acontainer, 
			            final int method, final boolean hm, final boolean cdk) {
		this.method = method;
		this.hashMap = hm;
		this.cdkCheck = cdk;
		this.atoms = atoms;
		this.nH = nH;
		this.maxOpenings = maxOpenings;
		graph = gr; 
		this.adjacency = new int [atoms.length][atoms.length];
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++){
				this.adjacency[canPerm[i]][canPerm[j]] = adjacency[i][j];
			}
		this.canString = molString();
		startLeft=0;
		startRight=1;
		this.acontainer = acontainer;
	}
	
	private String molString(){
		StringBuilder buf = new StringBuilder();
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<i; j++){
				buf.append(adjacency[i][j]);
			}
		return buf.toString();
	}

	public int getBondOrder (final int l, final int r){
		return adjacency[l][r];
	}

	private int openings;
	private boolean addUpOpenings(final int atomNum) {
		if (atoms.length == atomNum) 
			return openings == nH;
		else {
			int bondSum = 0;
			for (int i=0;i<atoms.length;i++)
				bondSum += adjacency[atomNum][i];
			for (Integer valence:atoms[atomNum].valenceList) {
				if (valence >= bondSum) {
					openings += (valence - bondSum);
					if (openings <= nH && addUpOpenings(atomNum+1)) {
						return true;
					}
					openings -= (valence - bondSum);
				}
			}
		}
		return false;
	}
	
	public boolean isComplete() {
		if (maxOpenings == 0) return true;
		openings = 0;
		return addUpOpenings (0);
	}

	private boolean[] seen ;
	private boolean isConnectedDFS(){
		seen = new boolean[atoms.length];
		dfs(0);
		boolean connected = true;
		for (boolean b:seen) connected &= b;
		return connected;
	}
	
	private void dfs(int num) {
		seen[num] = true;
		for (int i=0; i<atoms.length; i++){
			if (adjacency[i][num]>0 && !seen[i]) dfs(i);
		}
	}
	
	/**
	 * Writes a molecule to a file in SDF format. This is more or less equivalent to writeTo method but it does not depend on CDK. 
	 * The good is that it does not have the restrictions of having a "valid" atom container and is expected to work faster.
	 * The bad side is that it has less features and cannot handle all complicated chemical structures. 
	 * @param outFile
	 * @param mol_counter
	 */
	public synchronized void writeMol(final BufferedWriter outFile, final long mol_counter) {
		StringWriter writer = new StringWriter();
		int bondsCount = 0;
		
		for (Atom atom : atoms) {
			writer.write(String.format("    0.0000    0.0000    0.0000 %-3s 0  0  0  0  0  0%n", atom.symbol));
		}
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<i; j++) {
				int order = adjacency[i][j];
				if (order > 0){
					bondsCount++;
					writer.write(String.format("%3d%3d%3d%3d%3d%3d%n",i+1, j+1, order, 0, 0, 0));
				}
			}
		
		writer.write("M  END\n> <Id>\n"+mol_counter+"\n\n> <can_string>\n"+canString+"\n\n$$$$\n");
		try {
			// For the sake of calculating the bondsCount, we write the first line last!
			outFile.write(String.format("%n  PMG%n%n%3d%3d  0  0  0  0  0  0  0  0  0%n",atoms.length, bondsCount)
					      + writer.toString());
		} catch (IOException e) {
			System.err.println("Could not write molecule to output file.");
			e.printStackTrace();
		} 
	}
	
	private boolean incBondSemiCanBig(final int left, final int right){
		if (right<atoms.length-1 && atoms[right].symbol.equals(atoms[right+1].symbol) && 
			adjacency[left][right]==adjacency[left][right+1]) {
			int row;
			for (row=left; row>0; row--)  {
				if (adjacency[row-1][right]!=adjacency[row-1][right+1]) break;
			}
			if (row == 0)
				return false;
		}
		return incBond(left, right);
	}
	
	private boolean incBondSemiCan(final int left, final int right){
		if (right>left+1 && atoms[right].symbol.equals(atoms[right-1].symbol) && 
			adjacency[left][right]==adjacency[left][right-1]) {
			int row;
			for (row=left; row>0; row--)  {
				if (adjacency[row-1][right]!=adjacency[row-1][right-1]) break;
			}
			if (row == 0)
				return false;
		}
		return incBond(left, right);
	}
	
	private boolean incBond(final int left, final int right) {
		if (isFull(left))  return false;
		if (isFull(right)) return false;

		if (adjacency[left][right] > 2) return false;	// no more than Triple bonds
		adjacency[left][right]++;
		adjacency[right][left]++;
		return true;
	}

	private boolean decBond(final int left, final int right) {
		if (adjacency[left][right] <= 0) return false;	// no more than Triple bonds
		adjacency[left][right]--;
		adjacency[right][left]--;
		return true;
	}

	/**
	 * @param atom
	 */
	private boolean isFull(final int atom) {
		int bondSum = 0;
		for (int i=0;i<atoms.length;i++)
			bondSum += adjacency[atom][i];
		return (atoms[atom].maxValence <= bondSum);
	}
		
	//---------------------------------------------------------------
	class SortCompare {
		int maxRow = atoms.length - 1;
		boolean equal;	
		int[] iPerm;
		
		SortCompare(){}
		SortCompare (int max){
			maxRow = max;
		}
		
		boolean isSmallerThan (int i, int base) {
			equal = false;
			int []tPerm = iPerm.clone();
			for (int xCol=base+1; xCol<atoms.length; xCol++){
				int maxIndex=-1;
				int maxVal = -1;
				for (int yCol=atoms.length-1; yCol>=0; yCol--) {
					if (maxVal <= adjacency[i][yCol] && tPerm[yCol] == -1 && atoms[xCol].symbol.equals(atoms[yCol].symbol) && areCompatible(yCol, xCol, base)) {
						maxIndex = yCol;
						maxVal = adjacency[i][yCol];
					}
				}
//				if (maxIndex == -1) return false;	must not happen
				if (adjacency[base][xCol] < maxVal) return true;
				if (adjacency[base][xCol] > maxVal) return false;
				tPerm[maxIndex] = xCol;	// in order to avoid using the same value twice
			}
			equal = true;
			return false;
		}
		
		@SuppressWarnings("unchecked")
		boolean isMinimal(){
			iPerm = new int[atoms.length];
			Arrays.fill(iPerm, -1);
			return sortCheckMin2(0);
		}
		
		private boolean sortCheckMin2(int base) {
			if (base == maxRow) return true;
			for (int row=0; row<=maxRow; row++){	// even check row == base
				// check what happens if base <- row ?
				if (iPerm[row] != -1 || !atoms[base].symbol.equals(atoms[row].symbol)) continue;
				if (!areCompatible(row, base, base)) continue;
				iPerm[row] = base;
				if (isSmallerThan(row, base)) return false;
				if (equal && !sortCheckMin2(base+1)) return false;
				iPerm[row] = -1;
			}
			return true;
		}
		
		private boolean areCompatible(int src, int dst, int base) {
			for (int above=0; above<atoms.length; above++) {
				if (iPerm[above] != -1 && iPerm[above]!=base && adjacency[above][src] != adjacency[iPerm[above]][dst]) return false;
			}
			return true;
		}
	}

	//----------------------------------------------------------------------------------------

	private void outputMatrix(BufferedWriter out) {
		StringWriter writer = new StringWriter();
		writer.write(PMG.formula+"\n");
		writer.write("-----------"+(PMG.molCounter.get()-duplicate.get())+"\n");
		for (int i=0; i<atoms.length; i++) {
			for (int j=0; j<atoms.length; j++){
				writer.write(""+adjacency[i][j]);
			}
			writer.write("\n");
		}
		try {
			out.write(writer.toString());
		} catch (IOException e) {
			System.err.println("Could not output to Matrix.");
			e.printStackTrace();
		}
	}

	@Override
	public void run() {
		switch(method){
		case SEM_CAN:
		case MIN_CAN:
			generateOrderly(startLeft, startRight);
			PMG.availThreads.incrementAndGet();
			PMG.pendingTasks.decrementAndGet();
			return;
		}
		
		if (isComplete() && isConnectedDFS()) {
			long currentCount = PMG.molCounter.incrementAndGet();		
			writeToFile(currentCount);
		}	
		// get all possible ways to add one bond to the molecule
		ArrayList<MolProcessor> extMolList = addBond();

		for (MolProcessor molecule : extMolList) {
			if (method == BRT_FRC && molSet.add(molecule.canString) == false) continue;
			PMG.pendingTasks.incrementAndGet();
			PMG.startedTasks.incrementAndGet();
			if (PMG.taskQueue.size() > PMG.executorCount*2) {
				molecule.run();	// continue sequentially
			} else {
				PMG.executor.execute(molecule);
			}
		}
		PMG.pendingTasks.decrementAndGet();
	}
	
	
	/*
	 * This function should be called first with left=0 and right=1 as parameters.
	 */
	void generateOrderly (int left, int right) {
		if (left == atoms.length-1) return;
		if (right == atoms.length || isFull(left)) {
			submitIfAvailThread(left+1, left+2);
			return;
		}
		else if (incBondSemiCan(left, right)) {
			maxOpenings-=2;
			if (new SortCompare(left).isMinimal()) {
				checkMolecule();
				if (maxOpenings>0) submitIfAvailThread(left, right);	// if there are still open places for new bonds
			}
			decBond(left, right);
			maxOpenings+=2;
		}
		submitIfAvailThread(left, right+1);
	}

	private void submitIfAvailThread(int left, int right) {
		int avail = PMG.availThreads.get();
		while (avail>0) {
			if (PMG.availThreads.compareAndSet(avail, avail-1)) {
				PMG.pendingTasks.incrementAndGet();
				PMG.executor.execute(new MolProcessor(atoms, nH, maxOpenings, adjacency, graph, left, right, acontainer, method, hashMap, cdkCheck));
				return;
			}
			avail = PMG.availThreads.get();
		}
		generateOrderly(left, right);
	}
	
	private void checkMolecule() {
		if (isComplete() && isConnectedDFS()) {
			long currentCount = PMG.molCounter.incrementAndGet();
			MolProcessor newMol = this;
			if (hashMap) {
				// canonize 
				int[] perm1 = graph.canonize(this, true);	// ask for the automorphisms to be reported back
				newMol = new MolProcessor(atoms, nH, maxOpenings, adjacency, graph, perm1, acontainer, method, hashMap, cdkCheck);
			}
			if (hashMap ? !molSet.add(newMol.canString):!new SortCompare().isMinimal()) {
				duplicate.incrementAndGet();
			} else {
				writeToFile(currentCount);
			}
		}	
	}

	private void writeToFile(long currentCount) {
		BufferedWriter theOutFile = PMG.outFile;
		if (!acceptedByCDK()){
			PMG.rejectedByCDK.incrementAndGet();
			theOutFile = PMG.rejectedFile;
		}
		if(PMG.wFile){
			writeMol(theOutFile, currentCount);
//					outputMatrix(PMG.matrixFile);
		}
	}

	final static SaturationChecker satCheck = new SaturationChecker();
	private boolean acceptedByCDK() {
		try {
			IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();
			for (int r=0; r<atoms.length; r++)
				for (int c=r+1; c<atoms.length; c++) 
					if (adjacency[r][c] > 0) acprotonate.addBond(r, c, bondOrder(adjacency[r][c]));
			CDKAtomTypeMatcher typeMatcher = CDKAtomTypeMatcher.getInstance(acprotonate.getBuilder());
			for (IAtom atom : acprotonate.atoms()) {
				IAtomType type;
				type = typeMatcher.findMatchingAtomType(acprotonate, atom);
				if (type == null) return false;
				AtomTypeManipulator.configure(atom, type);	// TODO: What does this line mean? Is this method correct at all?!
			}
			CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(acprotonate.getBuilder());
			hAdder.addImplicitHydrogens(acprotonate);
			return (satCheck.isSaturated(acprotonate)&&(AtomContainerManipulator.getTotalHydrogenCount(acprotonate)==nH));
		} catch (CDKException e) {
			return false;
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return true;	// if things go wrong, we assume it is accepted! :P
	}


	private Order bondOrder(int i) {
		switch (i){
			case 1: return IBond.Order.SINGLE;
			case 2: return IBond.Order.DOUBLE;
			case 3: return IBond.Order.TRIPLE;
			case 4: return IBond.Order.QUADRUPLE;
		}
		return null;
	}

	private ArrayList<MolProcessor> addBondSemiCan() {
		ArrayList<MolProcessor> extMolList = new ArrayList<>();
		if (maxOpenings<=0) 
			return extMolList;	// the molecule is already saturated!

		for (int left = startLeft; left < atoms.length; left++){
			for (int right = left+1; right < atoms.length; right++){
				if (left == startLeft && right<startRight) continue;
				if (incBondSemiCan(left, right)) {	
					extMolList.add(new MolProcessor(atoms, nH, maxOpenings-2, adjacency, graph, left, right, acontainer, method, hashMap, cdkCheck)); 
					decBond(left,right);
				}
			}
		}
		return extMolList;	
	}

	/**
	 * Add one bond, with different orders, between left and right, and return the list.
	 * It then moves to the next place in the adjacency matrix for the next step.
	 * @return
	 */
	private ArrayList<MolProcessor> addBondSemiCan2() {
		ArrayList<MolProcessor> extMolList = new ArrayList<>();
		if (maxOpenings<=0) 
			return extMolList;	// the molecule is already saturated!

		int order = 0;
		int newLeft = startLeft;
		int newRight = startRight + 1;
		if (newRight == atoms.length) {
			newLeft ++;
			newRight = newLeft+1;
		}
		if (startLeft == atoms.length-1) return extMolList;
		do {
			extMolList.add(new MolProcessor(atoms, nH, maxOpenings-order, adjacency, graph, newLeft, newRight, acontainer, method, hashMap, cdkCheck));
			order += 2;
		}while(incBondSemiCan(startLeft, startRight));
		
		return extMolList;	
	}
		
	private ArrayList<MolProcessor> addBond() {
		ArrayList<MolProcessor> extMolList = new ArrayList<>();
		if (maxOpenings<=0) return extMolList;	// the molecule is already saturated!
		
    	Set<String> visited = new HashSet<>();

		// Note that the representative of an atom never has a bigger ID
		for (int left = 0; left < atoms.length; left++){
//			if (left>0 && rep[left] <= rep[left-1]) continue;	// make sure each orbit is considered only once
			for (int right = left+1; right < atoms.length; right++){
//				if (right>left+1 && rep[right] <= rep[right-1]) continue;	// make sure each orbit is considered only once
				// For the first iteration (in inner loop), we may consider the same orbit as "left"
				
					if (!incBond(left, right)) continue;	
					// canonize 
					int[] perm1 = graph.canonize(this, false);	// ask for the automorphisms to be reported back
					MolProcessor newMol = new MolProcessor(atoms, nH, maxOpenings-2, adjacency, graph, perm1, acontainer, method, hashMap, cdkCheck);
					if (visited.add(newMol.canString)) {	
						if (method == BRT_FRC || canString.equals("") || canString.equals(newMol.canDel())){ 
							extMolList.add(newMol); 
						}
					}	
					decBond(left,right);
			}
		}
		return extMolList;
	}
	
	private String canDel() {
		int left=0,right=0;
		for (int i = 0; i < atoms.length; i++){
			for (int j = i+1; j < atoms.length; j++){
				if (adjacency[i][j] != 0) {
					left = i;
					right = j;
					break;
				}
			}
		}
		decBond(left, right);
		int[] perm1 = graph.canonize(this, true);	// ask for the automorphisms to be reported back
		MolProcessor tempMol = new MolProcessor(atoms, nH, maxOpenings+2, adjacency, graph, perm1, acontainer, method, hashMap, cdkCheck);
		incBond(left, right);
		return tempMol.canString;
	}
}
