package org.omg;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
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
	public static final boolean canAug = true;
	private static final boolean semiCan = true;
	
	public final Atom[] atoms;
	final int nH;
	final int maxOpenings;
	final String canString;
	final Graph graph;
	final int[][] adjacency;
//	int []perm;
	static AtomicLong duplicate = new AtomicLong(0);
	final int startLeft, startRight;
	final int[] mPerm;
	
	private static final Set<String> molSet = Collections.synchronizedSet(new HashSet<String>());
	final LinkedList<int[]>[] allPerm;
	
	
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
		mPerm = new int[atoms.length];
		allPerm = new LinkedList[atoms.length];
		populate(allPerm);
//		perm = new int[atomCount];
//		for (int i=0; i<atomCount; i++) perm[i] = i;
	}
	
	private void populate(LinkedList<int[]>[] allPerm2) {
		int n = atoms.length-1;
		int permutation[] = new int[atoms.length];
		for (int i=0; i<atoms.length; i++){
			allPerm[i] = new LinkedList<int[]>();
			permutation[i] = i;
		}
		allPerm[n].add(permutation);
		while (--n>=0){
			for (int k=n+1;k<atoms.length; k++){
				if (!atoms[k].symbol.equals(atoms[n].symbol)) break;	// consider only permutations between the same atom types
				for (int kk=n+1;kk<atoms.length; kk++)
				for (int[] p : allPerm[kk]){
					permutation = new int[atoms.length];
					for (int i=0;i<atoms.length; i++) permutation[i] = p[i];
					permutation[k] = p[n];
					permutation[n] = p[k];
					allPerm[n].add(permutation);
				}
			}
		}
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
			            final int[][] adjacency, final Graph gr, final int stL, final int stR, LinkedList<int[]>[] allPerm) {
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
		mPerm = new int[atoms.length];
		this.allPerm = allPerm;
	}
		 	
	public MolProcessor(final Atom[] atoms, final int nH, final int maxOpenings, 
			            final int[][] adjacency, final Graph gr, final int[] canPerm, LinkedList<int[]>[] allPerm) {
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
		mPerm = new int[atoms.length];
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
	private boolean genPerm(int i){
		if (i==0) 
			return checkMinimality();
		for (int p=0; p<atoms.length; p++) {
			if (mPerm[p] == 0 && atoms[p].symbol.equals(atoms[i].symbol)) {
				mPerm[p] = i;
				if (! genPerm(i-1))
					return false;
				mPerm[p] = 0;
			}
		}
		return true;
	}

	private boolean isMinimalOrderly(){
		nextLevel: for (int ip=atoms.length-2;ip>=0;ip--){
			nextPerm: for(int[]p : allPerm[ip]){
				for (int i=0; i<atoms.length; i++){
					for (int j=i+1; j<atoms.length; j++){
						if (adjacency[i][j] < adjacency[p[i]][p[j]]) return false;	// permuted is smaller
						if (adjacency[i][j] > adjacency[p[i]][p[j]]) continue nextPerm;	// original is smaller
					}
				}
				continue nextLevel;	// if (permuted == original) skip to next level
			}
		}
		return true;
	}

	private boolean checkMinimality() {
		int[][] pAdjacency = new int [atoms.length][atoms.length];
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++){
				pAdjacency[mPerm[i]][mPerm[j]] = adjacency[i][j];
			}
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++){
				if (adjacency[i][j] == pAdjacency[i][j]) continue;
				return (adjacency[i][j] > pAdjacency[i][j]);
			}
		return true;
	}

	private boolean isMinimal(){
		return genPerm(atoms.length-1);
	}

	@Override
	public void run() {
		if (isComplete() && isConnectedDFS()) {
			//					if (!molSet.add(mol.canString)) System.err.println("Duplicate");
//			MolProcessor newMol = this;
//			if (semiCan) {
//				// canonize 
//				int[] perm1 = graph.canonize(this, true);	// ask for the automorphisms to be reported back
//				newMol = new MolProcessor(atoms, nH, maxOpenings, adjacency, graph, perm1);
//			}
			long currentCount = PMG.molCounter.incrementAndGet();
			if (semiCan && !isMinimalOrderly()/*!molSet.add(newMol.canString)*/) duplicate.incrementAndGet();
			else {
				if(PMG.wFile){
					writeMol(PMG.outFile, currentCount);
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
			if (!canAug && molSet.add(molecule.canString) == false) continue;
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
	
	private int[] canForm(final int [] canPerm){
		int[] tempadj = new int [atoms.length];
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<atoms.length; j++)
				for (int k=0; k<adjacency[i][j]; k++){
					tempadj[canPerm[i]] <<= 8;
					tempadj[canPerm[i]] += canPerm[j];
				}
		return tempadj;
	}
}
