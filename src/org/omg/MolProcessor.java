package org.omg;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

import org.omg.tools.Atom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import fi.tkk.ics.jbliss.Graph;

public class MolProcessor implements Runnable{
	public static Atom[] atoms;
	static int nH;
	int maxOpenings = 0;
	String canString="";
	public int[] adjacency;
	int bondsCount=0;
	int []perm;
	Graph graph;
	
	public MolProcessor(int maxOpenings, String canString, int[] adjacency,	int bondsCount, int []perm) {
		super();
		this.maxOpenings = maxOpenings;
		this.canString = canString;
		this.adjacency = adjacency.clone();
		this.bondsCount = bondsCount;
		this.perm = perm;
		graph = new Graph(atoms.length);
	}

	public MolProcessor(ArrayList<String> atomSymbols){
		int nH=0;
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
		adjacency = new int [atomCount * atomCount];
		graph = new Graph(atomCount);
		perm = new int[atomCount];
		for (int i=0; i<atomCount; i++) perm[i] = i;
	}



	private int openings = 0;
	private boolean addUpOpenings(int atomNum) {
		if (atoms.length == atomNum) 
			return openings == nH;
		else {
			int bondSum = 0;
			int base = perm[atomNum]*atoms.length;
			for (int i=base;i<base+atoms.length;i++)
				bondSum += adjacency[i];
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
		return addUpOpenings (0);
	}

	public boolean isConnected() {
		boolean connectedToZero[] = new boolean [atoms.length];
		connectedToZero[0] = true;
		for (int i=0; i<atoms.length; i++){
			if (!connectedToZero[i]) continue;
			for (int j=0; j<atoms.length; j++) {
				if (adjacency[i*atoms.length+j]>0) connectedToZero[j] |= connectedToZero[i];
			}
		}
		boolean connected = true;
		for (boolean b:connectedToZero) connected &= b;
		return connected;
	}
	
	private String molString(int [] canPerm){
		StringBuffer buf = new StringBuffer();
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<i; j++){
				buf.append(adjacency[canPerm[perm[i]]*atoms.length+canPerm[perm[j]]]);
			}

		return buf.toString();
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
		writer.write(String.format("%n  PMG%n%n%3d%3d  0  0  0  0  0  0  0  0  0%n",atoms.length, bondsCount));
		for (Atom atom : atoms) {
			writer.write(String.format("    0.0000    0.0000    0.0000 %-3s 0  0  0  0  0  0%n", atom.symbol));
		}
		for (int i=0; i<atoms.length; i++)
			for (int j=0; j<i; j++) {
				int order = adjacency[perm[i]*atoms.length+perm[j]];
				if (order > 0)
					writer.write(String.format("%3d%3d%3d%3d%3d%3d%n",perm[i]+1, perm[j]+1, order, 0, 0, 0));
			}
		
		writer.write("M  END\n> <Id>\n"+mol_counter+"\n\n> <can_string>\n"+canString+"\n\n$$$$\n");
		try {
			outFile.write(writer.toString());
		} catch (IOException e) {
			System.err.println("Could not write molecule to output file.");
			e.printStackTrace();
		} 
	}
	

	private boolean incBond(int left, int right) {
		if (isFull(left))  return false;
		if (isFull(right)) return false;

		if (adjacency[left*atoms.length+right] > 2) return false;	// no more than Triple bonds
		adjacency[left*atoms.length+right]++;
		adjacency[right*atoms.length+left]++;
		return true;
	}

	private boolean decBond(int left, int right) {
		if (adjacency[left*atoms.length+right] <= 0) return false;	// no more than Triple bonds
		adjacency[left*atoms.length+right]--;
		adjacency[right*atoms.length+left]--;
		return true;
	}


	/**
	 * @param atom
	 */
	private boolean isFull(int atom) {
		int bondSum = 0;
		int base = atom*atoms.length;
		for (int i=base;i<base+atoms.length;i++)
			bondSum += adjacency[i];
		return (atoms[atom].maxValence <= bondSum);
	}
	

	@Override
	public void run() {
		if (isComplete() && isConnected()) {
			//					if (!molSet.add(mol.canString)) System.err.println("Duplicate");
			long currentCount = PMG.molCounter.incrementAndGet();
			if(PMG.wFile){
				writeMol(PMG.outFile, currentCount);
			}
		}	

		// get all possible ways to add one bond to the molecule
		ArrayList<MolProcessor> extMolList = addOneBond();

		for (MolProcessor molecule : extMolList) {
			if (PMG.taskQueue.size() > PMG.executorCount*2) {
				molecule.run();	// continue sequentially
			} else {
				PMG.executor.execute(molecule);
			}
		}
	}


	private ArrayList<MolProcessor> addBond(boolean canAug, boolean bigStep) {
		ArrayList<MolProcessor> extMolList = new ArrayList<>();
		if (maxOpenings<=0) return extMolList;	// the molecule is already saturated!
		
    	Set<String> visited = new HashSet<>();

		// Note that the representative of an atom never has a bigger ID
		for (int left = 0; left < atoms.length; left++){
//			if (left>0 && rep[left] <= rep[left-1]) continue;	// make sure each orbit is considered only once
			for (int right = left+1; right < atoms.length; right++){
//				if (right>left+1 && rep[right] <= rep[right-1]) continue;	// make sure each orbit is considered only once
				// For the first iteration (in inner loop), we may consider the same orbit as "left"
				
				if (bigStep && adjacency[left*atoms.length+right] > 0) continue;
				
				do {
					if (!incBond(left, right)) break;	

					// canonize 
					int[] perm1 = graph.canonize(this, true);	// ask for the automorphisms to be reported back
					String molString = molString(perm1);
					decBond(left,right);
					if (visited.add(molString) == false) break;	

					int[] orbit = graph.orbitRep;
					if (!canAug || bondsCount==0){ // no need to check canonical augmentation
						extMolList.add(new MolProcessor(maxOpenings-2, molString, adjacency, bondsCount, join(perm1, perm))); 
						continue;
					}

					assert (false);
/* for now let's forget about canonical augmentation
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
						extMolList.add(new MolHelper2(canExtMol, orbit, molString, maxOpenings-2)); 
					}
*/
				} while (bigStep);
			}
		}
		return extMolList;
	}

	private int[] join(int[] perm1, int[] perm2) {
		int perm[] = new int [atoms.length];
		for (int i=0; i<atoms.length; i++)
			perm[i] = perm1[perm2[i]];
		return perm;
	}

	ArrayList<MolProcessor> addOneBond() {
		return addBond(true, false);
	}

	ArrayList<MolProcessor> addOneBondNoCheck(){
		return addBond(false, false);
	}

	ArrayList<MolProcessor> addBigBond(){
		return addBond(false, true);
	}

}
