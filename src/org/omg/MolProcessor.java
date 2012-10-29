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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import fi.tkk.ics.jbliss.Graph;

public class MolProcessor implements Runnable{
	static final boolean semiCan = true;
	static final boolean canAug = !semiCan;
	static final boolean hashMap = false;	// otherwise, minimality check (with semiCan)
	
	public final Atom[] atoms;
	final int nH;
	final int maxOpenings;
	final String canString;
	final Graph graph;
	final int[][] adjacency;
//	int []perm;
	static AtomicLong duplicate = new AtomicLong(0);
	final int startLeft, startRight;
	
	private static final Set<String> molSet = Collections.synchronizedSet(new HashSet<String>());
	final LinkedList<LinkedList<int[]>> allPerm;
	
	
	public MolProcessor(final ArrayList<String> atomSymbols){
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
		allPerm = new LinkedList<LinkedList<int[]>>();
		//if (!hashMap && semiCan) populate(allPerm);
//		perm = new int[atomCount];
//		for (int i=0; i<atomCount; i++) perm[i] = i;
	}

	private int[] identity(int length) {
		int [] id = new int[length];
		for (int i=0; i<length; i++) id[i] = i;
		return id;
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
			            final int stL, final int stR, final LinkedList<LinkedList<int[]>> allPerm) {
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
		this.allPerm = allPerm;
	}
		 	
	public MolProcessor(final Atom[] atoms, final int nH, final int maxOpenings, 
			            final int[][] adjacency, final Graph gr, final int[] canPerm,
			            final LinkedList<LinkedList<int[]>> allPerm) {
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
		this.allPerm = allPerm;
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

	private int openings = 0;
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
		
	/***************************************************************************************/

	//---------------------------------------------------------------
	class RowCompare {
		boolean equal;
		private int base;
		private int[] row_i_sorted;
		
		RowCompare (int base){
			this.base = base;
		}
		
		boolean isSmallerThan (int i) {
			equal = false;
			row_i_sorted = sortInMatrix(i);		
			if (row_i_sorted == null) return false;
			for (int yCol=base+1; yCol<atoms.length; yCol++){	
				if (adjacency[base][yCol] < row_i_sorted[yCol]) return true;
				if (adjacency[base][yCol] > row_i_sorted[yCol]) return false;
			}
			equal = true;
			return false;
		}
		
		/**
		 * Sorts a given row (param i), considering all the rows and columns up to base (not including) to be fixed.
		 * It means that two numbers in row i can be swapped if the corresponding columns do not change anything in
		 * the rows above base.
		 *  
		 * @param i An integer showing the index of the row to be sorted.
		 * @param base An integer showing the base row, up to which (but not including base itself) all rows and columns are fixed.
		 * @return
		 */
		private int[] sortInMatrix(int i) {
			int []iRow = new int[atoms.length];
			int []tPerm = iPerm.clone();
			for (int tillBase=0; tillBase<atoms.length; tillBase++)
				if (iPerm[tillBase] != -1)
					iRow[iPerm[tillBase]] = adjacency[i][tillBase];
			for (int xCol=base+1; xCol<atoms.length; xCol++){
				int max=-1;
				for (int yCol=atoms.length-1; yCol>=0; yCol--) {
					if (iRow[xCol] <= adjacency[i][yCol] && tPerm[yCol] == -1 && atoms[xCol].symbol.equals(atoms[yCol].symbol) && areCompatible(yCol, xCol, base)) {
						iRow[xCol] = adjacency[i][yCol];
						max = yCol;
					}
				}
				if (max == -1) 
					return null;
				tPerm[max] = xCol;
			}
			return iRow;
		}

		private void swap(int[] iRow, int xCol, int yCol) {
			int t = iRow[xCol]; 
			iRow[xCol] = iRow[yCol]; 
			iRow[yCol] = t;
		}
	}
	
	LinkedList<Integer>[] trans;
	int[] iPerm;
	@SuppressWarnings("unchecked")
	boolean isMinimalSort(){
		trans = new LinkedList[atoms.length];
		iPerm = new int[atoms.length];
		Arrays.fill(iPerm, -1);
//		trans = new LinkedList[atoms.length];
		return sortCheckMin2(0); // && checkPerms(0, mPerm);
	}

	private boolean sortCheckMin2(int base) {
		if (base == atoms.length - 1) return true;
		trans[base] = new LinkedList<>();
		RowCompare baseRow = new RowCompare(base);
		for (int row=0; row<atoms.length; row++){	// even check row == base
			// check what happens if base <- row ?
			if (iPerm[row] != -1 || !atoms[base].symbol.equals(atoms[row].symbol)) continue;
//			if (row < base) {
//				if (!trans[row].contains(base)) continue;
//			} else if (row >= base) {
				if (!areCompatible(row, base, base)) continue;
				iPerm[row] = base;
				if (baseRow.isSmallerThan(row)) 
					return false;
				iPerm[row] = -1;
				if (! baseRow.equal) continue; 
//			}
			trans[base].add(row);
		}
		for (int row:trans[base]){
			iPerm[row] = base;
			if (!sortCheckMin2(base+1)) return false;
			iPerm[row] = -1;
		}
		return true;
	}

	private boolean canBeMapped(int row, int base) {
		int b = iPerm[base];
		while (b != -1 && b != base) {
			if (b == row) return true;
			b = iPerm[b];
		}
		return false;
	}

	private boolean areCompatible(int src, int dst, int base) {
		for (int above=0; above<atoms.length; above++) {
			if (iPerm[above] != -1 && iPerm[above]!=base && adjacency[above][src] != adjacency[iPerm[above]][dst]) return false;
		}
		return true;
	}

	//----------------------------------------------------------------------------------------
	/**
	 * This function works, but takes too long. Now we want to improve it by using "sorting" 
	 * instead of checking all permutations to find a bigger graph at each step.
	 * @return whether the current graph is minimal
	 */
	@SuppressWarnings("unchecked")
	boolean isMinimal0ToN(){
		iPerm = new int[atoms.length];
		Arrays.fill(iPerm, -1);
		return sortCheckMin2(0); // && checkPerms(0, mPerm);
	}
	private boolean checkMin2(int base) {
		if (base == atoms.length - 1) return true;
		RowPermute baseRow = new RowPermute(base);
		for (int row=0; row<atoms.length && atoms[base].symbol.equals(atoms[row].symbol); row++){	// even check row == base
			// check what happens if base <- row ?
			if (iPerm[row] != -1) continue;
			if (row < base) {
				if (!canBeMapped(row, base)) continue;
			} else if (row > base) {
				if (baseRow.existsBiggerPerm(row)) return false;
				if (! baseRow.equal) continue; 
			}
			iPerm[row] = base;
			if (!sortCheckMin2(base+1)) return false;
			iPerm[row] = -1;
		}
		return true;
	}
	
	class RowPermute{
		boolean equal;
		int base;
		int [] perm;
		
		RowPermute(int base){
			this.base = base;
		}
		
		boolean existsBiggerPerm(int row){
			equal = false;
			this.perm = iPerm.clone();
			this.perm[row] = base;
			if (this.existsBigger(base+1)) return true;
			return false;
		}

		private boolean existsBigger(int depth) {
			if (depth==atoms.length) return compare();
			for (int i=0; i<atoms.length; i++){
				if (perm[i] == -1) {
					perm[i] = depth;
					if (existsBigger(depth+1)) return true;
					perm[i] = -1;
				}
			}
			return false;
		}

		private boolean compare() {
			for (int row=0; row<atoms.length; row++)
				for (int col=row+1; col<atoms.length; col++) {
					if (adjacency[row][col] < adjacency[perm[row]][perm[col]]) return true;
					if (adjacency[row][col] > adjacency[perm[row]][perm[col]]) return false;
				}
			equal = true;
			return false;
		}
	}
	//----------------------------------------------------------------------------------------

	@Override
	public void run() {
		if (isComplete() && isConnectedDFS()) {
			//					if (!molSet.add(mol.canString)) System.err.println("Duplicate");
			long currentCount = PMG.molCounter.incrementAndGet();
			MolProcessor newMol = this;
			if (semiCan && hashMap) {
				// canonize 
				int[] perm1 = graph.canonize(this, true);	// ask for the automorphisms to be reported back
				newMol = new MolProcessor(atoms, nH, maxOpenings, adjacency, graph, perm1, allPerm);
			}
			if (semiCan && (hashMap ? !molSet.add(newMol.canString):!isMinimalSort())) {
				duplicate.incrementAndGet();
/*/			
			int[] perm1 = graph.canonize(this, true);	// ask for the automorphisms to be reported back
			int i;
			for (i=0; i<atoms.length; i++) if (perm1[i] != i) break;
			if (i<atoms.length){
				 duplicate.incrementAndGet();
*/
			} else {
				if(PMG.wFile){
					writeMol(PMG.outFile, currentCount);
					outputMatrix(PMG.matrixFile);
				}
			}
		}	

		// get all possible ways to add one bond to the molecule
		ArrayList<MolProcessor> extMolList;
		if (semiCan)
			extMolList = addBondSemiCan();
		else 
			extMolList = addBond();

		for (MolProcessor molecule : extMolList) {
			if (!canAug && !semiCan && molSet.add(molecule.canString) == false) continue;
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

	private ArrayList<MolProcessor> addBondSemiCan() {
		ArrayList<MolProcessor> extMolList = new ArrayList<>();
		if (maxOpenings<=0) 
			return extMolList;	// the molecule is already saturated!

		for (int left = startLeft; left < atoms.length; left++){
			for (int right = left+1; right < atoms.length; right++){
				if (left == startLeft && right<startRight) continue;
				if (incBondSemiCan(left, right)) {	
					extMolList.add(new MolProcessor(atoms, nH, maxOpenings-2, adjacency, graph, left, right, allPerm)); 
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
			extMolList.add(new MolProcessor(atoms, nH, maxOpenings-order, adjacency, graph, newLeft, newRight, allPerm));
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
					int[] perm1 = graph.canonize(this, true);	// ask for the automorphisms to be reported back
					MolProcessor newMol = new MolProcessor(atoms, nH, maxOpenings-2, adjacency, graph, perm1, allPerm);
					if (visited.add(newMol.canString)) {	
						if (!canAug || canString.equals("") || canString.equals(newMol.canDel())){ 
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
		MolProcessor tempMol = new MolProcessor(atoms, nH, maxOpenings+2, adjacency, graph, perm1, allPerm);
		incBond(left, right);
		return tempMol.canString;
	}
}
