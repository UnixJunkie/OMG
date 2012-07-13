package org.omg;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;


public class MolManipulator {

	
	   public static IAtomContainer getcanonical(IAtomContainer ac) throws CloneNotSupportedException{

		   IAtomContainer ac2 = ac.getBuilder().newInstance(IAtomContainer.class);
			int vCount = ac.getAtomCount();
			int arr1[] = new int[vCount*vCount]; 
			int rr1[] = new int[vCount*vCount];
			int lab[] = new int[vCount];
			int lab1[] = new int[vCount];
			int ptn[] = new int[vCount];
			
			int c = 0;
			IBond bond = null;
			for (int i=0; i<vCount; i++){
				lab[i] = i;
				if(i == vCount-1){
					ptn[i]=0;
				}
				else 
				{
					if(ac.getAtom(i).getSymbol().equals(ac.getAtom(i+1).getSymbol())){
						ptn[i]=1;
					}
					else{
						ptn[i]=0;
					}
				}
				bond = null;
				for (int j=0; j<vCount; j++){
					bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
					if(bond == null){}
					else{
						arr1[c] = bond.getOrder().ordinal() + 1;
					}
					c++;
				}
			}
			int ret[] = OMGJNI.getcanmultig(vCount, arr1, lab, ptn);

			System.arraycopy(ret, 0, rr1, 0, vCount*vCount);
			for (int i = 0; i < vCount;i++){
				lab1[i] = ret[(vCount*vCount)+i];
				ac2.addAtom(ac.getAtom(lab1[i]));
			}

			for (int i=0; i<vCount; i++){
				for (int j=i+1; j<vCount; j++){
					if(rr1[i*vCount+j]==1){
						ac2.addBond(i, j, IBond.Order.SINGLE);
						if(ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("BondINfrag")!= null){
							ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("BondINfrag", ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("BondINfrag"));
						}
					}	
					else if(rr1[i*vCount+j]==2){
						ac2.addBond(i, j, IBond.Order.DOUBLE);
						if(ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("BondINfrag")!= null){
							ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("BondINfrag", ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("BondINfrag"));
						}
					}
					else if(rr1[i*vCount+j]==3){
						ac2.addBond(i, j, IBond.Order.TRIPLE);
						if(ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("BondINfrag")!= null){
							ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("BondINfrag", ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("BondINfrag"));
						}
					}					
				}
			}	

			return ac2;
		}
	
	   public static ArrayList<IAtomContainer> getcanonical_list(ArrayList<IAtomContainer> extMollist) throws CloneNotSupportedException{

		   int molcount = extMollist.size();
		
			int vCount = extMollist.get(0).getAtomCount();

			int arr1[] = new int[molcount*vCount*vCount]; 
			int rr1[] = new int[vCount*vCount];
			int lab[] = new int[molcount*vCount];
			int lab1[] = new int[vCount];
			int ptn[] = new int[molcount*vCount];
			int count = 0;
			
			int c = 0;
			for(IAtomContainer mol : extMollist){
				for (int i=0; i<vCount; i++){
					lab[(count*vCount)+i] = i;
					if(i == vCount-1){
						ptn[(count*vCount)+i]=0;
					}
					else 
					{
						if(mol.getAtom(i).getSymbol().equals(mol.getAtom(i+1).getSymbol())){
							ptn[(count*vCount)+i]=1;
						}
						else{
							ptn[(count*vCount)+i]=0;
						}
					}
					for (int j=0; j<vCount; j++){
						IBond bond = mol.getBond(mol.getAtom(i), mol.getAtom(j));
						if(bond == null){}
						else if(bond.getOrder() == IBond.Order.SINGLE){
							arr1[c] = 1;
						}
						else if(bond.getOrder() == IBond.Order.DOUBLE){
							arr1[c] = 2;
						}
						else if(bond.getOrder() == IBond.Order.TRIPLE){
							arr1[c] = 3;
						}
						else if(bond.getOrder() == IBond.Order.QUADRUPLE){
							arr1[c] = 4;
						}
						c++;
					}
					
				}
				count++;
			}
			int ret[] = OMGJNI.getcanmultig2(molcount,vCount, arr1, lab, ptn);

			ArrayList<IAtomContainer> listreturn = new ArrayList<IAtomContainer>();
			for(int mc = 0; mc < molcount;mc++){
				IAtomContainer ac2 = extMollist.get(mc).getBuilder().newInstance(IAtomContainer.class);

			
				System.arraycopy(ret, (mc*((vCount*vCount)+vCount)), rr1, 0, (vCount*vCount));
				for (int i = 0; i < vCount;i++){
					
					lab1[i] = ret[((((vCount*vCount)+vCount)*mc)+(vCount*vCount))+i];
					ac2.addAtom(extMollist.get(mc).getAtom(lab1[i]));
				}
	
				for (int i=0; i<vCount; i++){
					for (int j=i+1; j<vCount; j++){
						if(rr1[i*vCount+j]==1){
							ac2.addBond(i, j, IBond.Order.SINGLE);
							if(extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag")!= null){
								ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("BondINfrag", extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag"));
							}
						}	
						else if(rr1[i*vCount+j]==2){
							ac2.addBond(i, j, IBond.Order.DOUBLE);
							if(extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag")!= null){
								ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("BondINfrag", extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag"));
							}
						}
						else if(rr1[i*vCount+j]==3){
							ac2.addBond(i, j, IBond.Order.TRIPLE);
							if(extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag")!= null){
								ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("BondINfrag", extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag"));
							}
						}
						else if(rr1[i*vCount+j]==4){
							ac2.addBond(i, j, IBond.Order.QUADRUPLE);
							if(extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag")!= null){
								ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("BondINfrag", extMollist.get(mc).getBond(extMollist.get(mc).getAtom(lab1[i]), extMollist.get(mc).getAtom(lab1[j])).getProperty("BondINfrag"));
							}
						}
						
					}
				}
				listreturn.add(ac2);
			}	

			return listreturn;
		}


public static  ArrayList<int[]> extendMol(IAtomContainer ac) throws CloneNotSupportedException, CDKException {
	int vCount = ac.getAtomCount();
	int [] bondCounts = new int [vCount];
	for (IBond b:ac.bonds()) {
		bondCounts[ac.getAtomNumber(b.getAtom(0))] += b.getOrder().ordinal()+1;
		bondCounts[ac.getAtomNumber(b.getAtom(1))] += b.getOrder().ordinal()+1;
	}
	ArrayList<int[]> bondList = new ArrayList<int[]>();		
	
	for (int i = 0; i < vCount; i++){
		if (maxBondTable.get(ac.getAtom(i).getSymbol()) == bondCounts[i]) continue;
		for (int j = i+1; j < vCount; j++){
			if (maxBondTable.get(ac.getAtom(j).getSymbol()) == bondCounts[j]) continue;
					bondList.add(new int[] { i, j});	
			}		   			   
		}			  
	return bondList;
}

	public static boolean aresame(IAtomContainer ac1, IAtomContainer ac2) {
		return Arrays.equals(mol2array(ac1),mol2array(ac2));
	}
	public static String array2string(int[] array) {
		StringBuffer b = new StringBuffer();	   
		for(int i = 0;i<array.length;i++){
			b.append(array[i]);		   
		}
		return b.toString();
	}


	public static int[] mol2array(IAtomContainer ac) {
		int aCount = ac.getAtomCount();
		int[] ar = new int[aCount*aCount];
		for(IBond bond : ac.bonds()){
			//we read only the bonds and store the bond degree as ar[i*aCount+j]
			ar[ac.getAtomNumber(bond.getAtom(0)) * aCount + ac.getAtomNumber(bond.getAtom(1))] = bond.getOrder().ordinal()+1;		
		}	
		return ar;
	}


	public static IAtomContainer buildFromFragment(IAtomContainer acontainer,
			IAtomContainer frag) throws CloneNotSupportedException {
		List<IAtom> listcont = new ArrayList<IAtom>();
		for(IAtom atom: frag.atoms()){
			if(atom.getSymbol().equals("H")){
				listcont.add(atom);
			}
		}
		for(IAtom atom: listcont){
			frag.removeAtomAndConnectedElectronContainers(atom);
		}
		for(IAtom atomF : frag.atoms()){			
				atomF.setFlag(1, false);
		}
		for(IAtom atom : acontainer.atoms()){
			for(IAtom atomF : frag.atoms()){
				if(atomF.getSymbol().equals(atom.getSymbol())&&(!atomF.getFlag(1))&&(!atom.getFlag(1))){
					atomF.setID(atom.getID());
					atomF.setFlag(1, true);
					atom.setFlag(1, true);
					break;
				}
			}
		}
		List<IAtom> atomFinal = new ArrayList<IAtom>();
		for(IAtom atomF : frag.atoms()){
			for(IAtom atom : acontainer.atoms()){
				if(atomF.getID().equals(atom.getID())){
					List<IAtom> atomFs = frag.getConnectedAtomsList(atomF);
					for(IAtom atomF2: atomFs){
						if(!atomFinal.contains(atomF2)){
							for(IAtom atom2 : acontainer.atoms()){
								if(atomF2.getID().equals(atom2.getID())){
									atomFinal.add(atomF);
									acontainer.addBond(acontainer.getAtomNumber(atom),acontainer.getAtomNumber(atom2),frag.getBond(atomF, atomF2).getOrder());
									IBond newbond = acontainer.getBond(atom, atom2);
									newbond.setProperty("BondINfrag", ""+acontainer.getBond(atom, atom2).getOrder());
									newbond.setProperty("checked", "yes");
								}
							}
						
							}
						}
				}
			}
		}

		return acontainer;
	}

	private static Map<String, Double> maxBondTable; 
	static {

		// initialize the table
		maxBondTable = new HashMap<String, Double>();
		// TODO: read atom symbols from CDK
		maxBondTable.put("C", new Double(4));
		maxBondTable.put("N", new Double(5));
		maxBondTable.put("O", new Double(2));
		maxBondTable.put("S", new Double(6));
		maxBondTable.put("P", new Double(5));
	}

}
