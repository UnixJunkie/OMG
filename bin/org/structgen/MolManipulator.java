package org.structgen;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;


public class MolManipulator {


	   public static IAtomContainer getcanonical(IAtomContainer ac){

		   IAtomContainer ac2 = ac.getBuilder().newAtomContainer();



			int vCount = ac.getAtomCount();
			int arr1[] = new int[vCount*vCount]; 
			int rr1[] = new int[vCount*vCount];
			int lab[] = new int[vCount];
			int lab1[] = new int[vCount];
			int ptn[] = new int[vCount];
			
			int c = 0;
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
				for (int j=0; j<vCount; j++){
					IBond bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
					
					if(bond == null)
						arr1[c] = 0;
					else if(bond.getOrder() == IBond.Order.SINGLE){
						if(bond.getProperty("Bond2frag")!= null){
							arr1[c] = Integer.valueOf("1"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[c] = 1;
						}
						
					}
					else if(bond.getOrder() == IBond.Order.DOUBLE){
						
						if(bond.getProperty("Bond2frag")!= null){
							arr1[c] = Integer.valueOf("2"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[c] = 2;
						}
					}
					else if(bond.getOrder() == IBond.Order.TRIPLE){
						
						if(bond.getProperty("Bond2frag")!= null){
							arr1[c] = Integer.valueOf("3"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[c] = 3;
						}
					}
					else if(bond.getOrder() == IBond.Order.QUADRUPLE){
						
						if(bond.getProperty("Bond2frag")!= null){
							arr1[c] = Integer.valueOf("4"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[c] = 4;
						}
					}
					c++;
				}
			}
			int ret[] = StructGenJNI.getcanmultig(vCount, arr1, lab, ptn);

			for (int i = 0; i < vCount*vCount;i++){
				rr1[i] = ret[i];
			}
			for (int i = 0; i < vCount;i++){
				lab1[i] = ret[(vCount*vCount)+i];
			}

			
//			int rr1[] = StructGenJNI.getcanmultig(vCount, arr1, lab, ptn);
//			int lab1[] = StructGenJNI.getlabmultig(vCount, arr1, lab, ptn);
			
			for(int i = 0; i < lab1.length; i++){
				ac2.addAtom(ac.getAtom(lab1[i]));
			}
			
			for (int i=0; i<vCount; i++){
				for (int j=i+1; j<vCount; j++){
					String st = Integer.toString(rr1[i*vCount+j]);
//					System.out.println("rr1 " + rr1[i*vCount+j]+"="+st.substring(0, 1));
					if(st.substring(0, 1).equals("1")){
						
						
						ac2.addBond(i, j, IBond.Order.SINGLE);
						if(ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag")!= null){
							ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("Bond2frag", ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag"));
						}
					}	
					else if(st.substring(0, 1).equals("2")){
						ac2.addBond(i, j, IBond.Order.DOUBLE);
						if(ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag")!= null){
							ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("Bond2frag", ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag"));
						}
					}
					else if(st.substring(0, 1).equals("3")){
						ac2.addBond(i, j, IBond.Order.TRIPLE);
						if(ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag")!= null){
							ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("Bond2frag", ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag"));
						}
					}
					else if(st.substring(0, 1).equals("4")){
						ac2.addBond(i, j, IBond.Order.QUADRUPLE);
						if(ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag")!= null){
							ac2.getBond(ac2.getAtom(i), ac2.getAtom(j)).setProperty("Bond2frag", ac.getBond(ac.getAtom(lab1[i]), ac.getAtom(lab1[j])).getProperty("Bond2frag"));
						}
					}
					
				}
			}
			

			return ac2;
		}
//	   public static IAtomContainer getcanonical(IAtomContainer ac){
//		   int countFatoms = 0;
//		   
//		   for(IAtom afrag : ac.atoms()){
//			   if(afrag.getFlag(1)){
//				   countFatoms++;
//			   }
//		   }
//		   
//		   IAtomContainer ac2 = ac.getBuilder().newAtomContainer();
//
//		   
//
//			int vCount = ac.getAtomCount()-countFatoms+1;
//			int arr1[] = new int[vCount*vCount]; 
//			int rr1[] = new int[vCount*vCount];
//			int lab[] = new int[vCount];
//			int lab1[] = new int[vCount];
//			int ptn[] = new int[vCount];
//			
//			int c = 0;
//			for (int i=0; i<vCount; i++){
//				lab[i] = i;
//				if(i == vCount-1){
//					ptn[i]=0;
//				}
//				else 
//				{
//					if(ac.getAtom(i).getSymbol().equals(ac.getAtom(i+1).getSymbol())){
//						ptn[i]=1;
//					}
//					else{
//						ptn[i]=0;
//					}
//				}
//				for (int j=0; j<vCount; j++){
//					IBond bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
//					if(bond == null)
//						arr1[c] = 0;
//					else if(bond.getOrder() == IBond.Order.SINGLE){
//						arr1[c] = 1;
//					}
//					else if(bond.getOrder() == IBond.Order.DOUBLE){
//						arr1[c] = 2;
//					}
//					else if(bond.getOrder() == IBond.Order.TRIPLE){
//						arr1[c] = 3;
//					}
//					else if(bond.getOrder() == IBond.Order.QUADRUPLE){
//						arr1[c] = 4;
//					}
//					c++;
//				}
//			}
//			int ret[] = StructGenJNI.getcanmultig(vCount, arr1, lab, ptn);
//
//			for (int i = 0; i < vCount*vCount;i++){
//				rr1[i] = ret[i];
//			}
//			for (int i = 0; i < vCount;i++){
//				lab1[i] = ret[(vCount*vCount)+i];
//			}
//
//			
////			int rr1[] = StructGenJNI.getcanmultig(vCount, arr1, lab, ptn);
////			int lab1[] = StructGenJNI.getlabmultig(vCount, arr1, lab, ptn);
//			
//			for(int i = 0; i < lab1.length; i++){
//				ac2.addAtom(ac.getAtom(lab1[i]));
//			}
//			
//			for (int i=0; i<vCount; i++){
//				for (int j=i+1; j<vCount; j++){
//					if(rr1[i*vCount+j]==1){
////						System.out.println("simple " +i +" "+ j);
//						ac2.addBond(i, j, IBond.Order.SINGLE);
//					}	
//					else if(rr1[i*vCount+j]==2){
//						ac2.addBond(i, j, IBond.Order.DOUBLE);
//					}
//					else if(rr1[i*vCount+j]==3){
//						ac2.addBond(i, j, IBond.Order.TRIPLE);
//					}
//					else if(rr1[i*vCount+j]==4){
//						ac2.addBond(i, j, IBond.Order.QUADRUPLE);
//					}
//				}
//			}
//			
//
//			return ac2;
//		}
	public static ArrayList<IAtomContainer> extendMol(IAtomContainer ac) throws CloneNotSupportedException, CDKException {
		int vCount = ac.getAtomCount();

		ArrayList<IAtomContainer> extendedMol = new ArrayList<IAtomContainer>();		

		for (int i = 0; i < vCount; i++){
		   for (int j = i+1; j < vCount; j++){
			   IBond bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
			   IAtomContainer acCloned = (IAtomContainer) ac.clone();
			   
			   boolean accept = true;
			   if(ac.getAtom(i).getFlag(1) && ac.getAtom(j).getFlag(1))
				   accept = false;
			   
				if(bond == null){
					if(accept)
						acCloned.addBond(i, j, IBond.Order.SINGLE);
				}
				else if((bond.getOrder() == IBond.Order.SINGLE) && accept){
					acCloned.getBond(acCloned.getAtom(i), acCloned.getAtom(j)).setOrder(IBond.Order.DOUBLE);
				}
				else if(bond.getOrder() == IBond.Order.DOUBLE && accept){
					acCloned.getBond(acCloned.getAtom(i), acCloned.getAtom(j)).setOrder(IBond.Order.TRIPLE);
				}
				else if(bond.getOrder() == IBond.Order.TRIPLE && accept){
					acCloned.getBond(acCloned.getAtom(i), acCloned.getAtom(j)).setOrder(IBond.Order.QUADRUPLE);
				}
				else if(bond.getOrder() == IBond.Order.QUADRUPLE && accept){
					continue;
				}
				if(accept){
					CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(acCloned.getBuilder());
		            IAtomType type1 = matcher.findMatchingAtomType(acCloned,acCloned.getAtom(i));
		            IAtomType type2 = matcher.findMatchingAtomType(acCloned,acCloned.getAtom(j));
					if((type1 != null) && (type2 != null)){
						extendedMol.add(acCloned);
//						int countFlags = 0;
//						for(IAtom a1 : acCloned.atoms()){
//							if(a1.getFlag(1)){
//								countFlags++;
//							}
//						}
//						System.out.println("countFlags: "+countFlags);
	//						System.out.println("extended " + MolManipulator.matrixToString(MolManipulator.getMatrix(acCloned)));
					}
					else{
	//						System.out.println("adios" );
					}			
				}
	   		}			  
		}
		System.out.println("beforeRem(): "+extendedMol.size());
		ArrayList<IAtomContainer> noduplicateMol = removeduplicatecanon(extendedMol);
		System.out.println("afterRem(): "+noduplicateMol.size());

		return noduplicateMol;
	}

	public static ArrayList<IAtomContainer> removeduplicatecanon(ArrayList<IAtomContainer> extendedMol) {
		 HashMap<String,IAtomContainer> map = new HashMap<String,IAtomContainer>();
		 for( IAtomContainer ac : extendedMol){
			 map.put(MolManipulator.getCanonstr(ac), ac);
		 }
		 ArrayList<IAtomContainer> mollist = new ArrayList<IAtomContainer>();
		 for(String str : map.keySet()){
			 mollist.add(map.get(str));
		 }
		 return mollist;
	}

	public static String getCanonstr(IAtomContainer ac) {
		   IAtomContainer ac2 = null;
			try {
				ac2 = (IAtomContainer) ac.clone();
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			ac2.removeAllBonds();				
			int vCount = ac.getAtomCount();
			int arr1[] = new int[vCount*vCount];
			int rr1[] = new int[vCount*vCount];
			int lab[] = new int[vCount];
			int lab1[] = new int[vCount];
			int ptn[] = new int[vCount];
			int c = 0;
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
				for (int j=0; j<vCount; j++){
					IBond bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
					if(bond == null)
						arr1[c] = 0;
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
			//int rr1[] = StructGenJNI.getcanmultig(vCount, arr1, lab, ptn);
//			StructGenJNI.getcanmultig(vCount, arr1, lab, ptn);
			int ret[] = StructGenJNI.getcanmultig(vCount, arr1, lab, ptn);

			for (int i = 0; i < vCount*vCount;i++){
				rr1[i] = ret[i];
			}
			for (int i = 0; i < vCount;i++){
				lab1[i] = ret[(vCount*vCount)+i];
			}
//			int rr1[] = arr1;
			
			StringBuilder b = new StringBuilder();	   
			for(int i = 0;i<vCount*vCount;i++){
				b.append(rr1[i]);		   
			}
			return b.toString();
		}

	public static boolean aresame(IAtomContainer ac1, IAtomContainer ac2) {
		boolean aresame = false;
		int aCount1 = ac1.getAtomCount();
		int aCount2 = ac2.getAtomCount();
		int bCount1 = ac1.getBondCount();
		int bCount2 = ac2.getBondCount();
		
		if((aCount1==aCount2)&&(bCount1==bCount2)){
			String lab1 = MolManipulator.arrayToString(MolManipulator.getlab(ac1));
			String lab2 = MolManipulator.arrayToString(MolManipulator.getlab(ac2));
			if(lab1.equals(lab2)){
				String ptn1 = MolManipulator.arrayToString(MolManipulator.getptn(ac1));
				String ptn2 = MolManipulator.arrayToString(MolManipulator.getptn(ac2));
				if(ptn1.equals(ptn2)){
					String matrix1 = MolManipulator.matrixToString(MolManipulator.getMatrix(ac1));
					String matrix2 = MolManipulator.matrixToString(MolManipulator.getMatrix(ac2));
					if(matrix1.equals(matrix2)){
						aresame = true;
					}
				}
			}
		}	
		return aresame;
	}
	public static String arrayToString(int[] array) {
		StringBuilder b = new StringBuilder();	   
		for(int i = 0;i<array.length;i++){
			b.append(array[i]);		   
		}
		return b.toString();
	}
	public static String matrixToString(int[][] matrix) {
		StringBuilder b = new StringBuilder();	   
		for(int i = 0;i<matrix.length;i++){
			for(int j = 0;j<matrix.length;j++){
				b.append(matrix[i][j]);	
			}				   
		}		
		return b.toString();
	}
	public static int[] getlab(IAtomContainer ac) {
		int aCount = ac.getAtomCount();
		int lab[] = new int[aCount];
		for (int i=0; i<aCount; i++){		
			lab[i] = i;			
		}
		return lab;
	}
	public static int[] getptn(IAtomContainer ac) {
		int aCount = ac.getAtomCount();
		int ptn[] = new int[aCount];
		for (int i=0; i<aCount; i++){		
			if(i == aCount-1){
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
		}
		return ptn;
	}
	public static int[][] getMatrix(IAtomContainer ac) {
		int aCount = ac.getAtomCount();
		System.out.println("getBondCount: "+ac.getBondCount());

		int arr1[][] = new int[aCount][aCount];
		for (int i=0; i<aCount; i++){
			for (int j = i+1; j<aCount; j++){
				IBond bond = ac.getBond(ac.getAtom(i), ac.getAtom(j));
				if(StructGen.id_counter == 2){
					System.out.println(ac.getAtom(i).getID()+"["+i+"], "+ac.getAtom(j).getID()+"["+j+"], "+" getMatrix_bond: "+bond);

				}
				if(bond == null){
					arr1[i][j] = 0;
					arr1[j][i] = 0;
				}
				else{	
					if(bond.getOrder() == IBond.Order.SINGLE){
						if(bond.getProperty("Bond2frag")!= null){
							arr1[i][j] = Integer.valueOf("1"+(String)bond.getProperty("Bond2frag"));
							arr1[j][i] = Integer.valueOf("1"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[i][j] = 1;
							arr1[j][i] = 1;
						}
					}
					if(bond.getOrder() == IBond.Order.DOUBLE){
						if(bond.getProperty("Bond2frag")!= null){
							arr1[i][j] = Integer.valueOf("2"+(String)bond.getProperty("Bond2frag"));
							arr1[j][i] = Integer.valueOf("2"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[i][j] = 2;
							arr1[j][i] = 2;
						}
					}
					if(bond.getOrder() == IBond.Order.TRIPLE){
						if(bond.getProperty("Bond2frag")!= null){
							arr1[i][j] = Integer.valueOf("3"+(String)bond.getProperty("Bond2frag"));
							arr1[j][i] = Integer.valueOf("3"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[i][j] = 3;
							arr1[j][i] = 3;
						}
					}
					if(bond.getOrder() == IBond.Order.QUADRUPLE){
						if(bond.getProperty("Bond2frag")!= null){
							arr1[i][j] = Integer.valueOf("4"+(String)bond.getProperty("Bond2frag"));
							arr1[j][i] = Integer.valueOf("4"+(String)bond.getProperty("Bond2frag"));
						}
						else{
							arr1[i][j] = 4;
							arr1[j][i] = 4;
						}
					}
				}				
			}
		}	
		return arr1;
	}

	public static IAtomContainer builFromFragment(IAtomContainer acontainer,
			IMolecule frag) throws CloneNotSupportedException {
		List<IAtom> listcont = new ArrayList<IAtom>();
		System.out.println("before H " + MolManipulator.matrixToString(MolManipulator.getMatrix(frag)));
		for(IAtom atom: frag.atoms()){
			if(atom.getSymbol().equals("H")){
				listcont.add(atom);
			}
		}
		for(IAtom atom: listcont){
			frag.removeAtomAndConnectedElectronContainers(atom);
		}
		System.out.println("after H  " + MolManipulator.matrixToString(MolManipulator.getMatrix(frag)));
		for(IAtom atomF : frag.atoms()){			
				atomF.setFlag(1, false);
		}
		for(IAtom atom : acontainer.atoms()){
			for(IAtom atomF : frag.atoms()){
				if(atomF.getSymbol().equals(atom.getSymbol())&&(!atomF.getFlag(1))){
					atomF.setID(atom.getID());
					System.out.println(" atomF "+ atomF.getID());
					atomF.setFlag(1, true);
					atom.setFlag(1, true);
					break;
				}
			}
		}
		List<IAtom> atomFinal = new ArrayList<IAtom>();
		for(IAtom atomF : frag.atoms()){
			System.out.println("atomF : " + atomF.getID());
			for(IAtom atom : acontainer.atoms()){
				System.out.println("atom : " + atom.getID() + " is equal? " + atomF.getID().equals(atom.getID()));
				if(atomF.getID().equals(atom.getID())){
					List<IAtom> atomFs = frag.getConnectedAtomsList(atomF);
					for(IAtom atomF2: atomFs){
						System.out.println("atomF connected to : " + atomF2.getID());
						if(!atomFinal.contains(atomF2)){
							for(IAtom atom2 : acontainer.atoms()){
								if(atomF2.getID().equals(atom2.getID())){
									atomFinal.add(atomF);
									acontainer.addBond(acontainer.getAtomNumber(atom),acontainer.getAtomNumber(atom2),frag.getBond(atomF, atomF2).getOrder());
								}
							}
						
							}
						}
				}
			}
		}
//		System.out.println("atoms frag");
//		for(IAtom a : frag.atoms()){
//			System.out.println(a);
//		}
//		System.out.println("atoms acon");
//		for(IAtom a : acontainer.atoms()){
//			System.out.println(a);
//		}
		System.out.println("acontainer " + MolManipulator.matrixToString(MolManipulator.getMatrix(acontainer)));
		//begin new part
//		IAtomContainer canonacontainer = MolManipulator.getcanonical(acontainer);
//		String lastBondID[] = new String[2];
//		lastBondID[0] = canonacontainer.getBond(canonacontainer.getBondCount()-1).getAtom(0).getID();
//		lastBondID[1] = canonacontainer.getBond(canonacontainer.getBondCount()-1).getAtom(1).getID();
//		IAtomContainer m_ext_e = (IAtomContainer) acontainer.clone();
//		int[] lastBond = new int[2];
//		for(IAtom atom : m_ext_e.atoms()){
//			if(atom.getID().equals(lastBondID[0]))
//				lastBond[0]=m_ext_e.getAtomNumber(atom);
//			if(atom.getID().equals(lastBondID[1]))
//				lastBond[1] = m_ext_e.getAtomNumber(atom);
//		}
//		IBond bondremove = m_ext_e.getBond(m_ext_e.getAtom(lastBond[0]),m_ext_e.getAtom(lastBond[1]));
//		if(bondremove.getOrder() == IBond.Order.SINGLE){
//			m_ext_e.removeBond(bondremove);
//		}
//		else if(bondremove.getOrder() == IBond.Order.DOUBLE){
//			bondremove.setOrder(IBond.Order.SINGLE);
//		}
//		else if(bondremove.getOrder() == IBond.Order.TRIPLE){
//			bondremove.setOrder(IBond.Order.DOUBLE);
//		}
//		else if(bondremove.getOrder() == IBond.Order.QUADRUPLE){
//			bondremove.setOrder(IBond.Order.TRIPLE);
//		}
//		return m_ext_e;
		//end new part
		return acontainer;
	}
	public static IAtomContainer getreduction(IAtomContainer mExt) throws CloneNotSupportedException {
		IAtomContainer acred = (IAtomContainer) mExt.clone();
		System.out.println(" bondsBefo "+ acred.getBondCount());
		int countF = 0;
		for(IAtom a1 : acred.atoms()){

//			System.out.println(" a0 "+ a1.getID()+";"+acred.getConnectedAtomsList(a1).size());
			if(a1.getFlag(1))
				countF++;
		}
		System.out.println("coutF : "+countF);
		IAtom aR = new Atom("R");
		aR.setID("R");
		acred.addAtom(aR);
		List <IAtom> atomlist = new ArrayList<IAtom>();
		for(IAtom a1 : acred.atoms()){
			if(a1.getFlag(1)){
				atomlist.add(a1);
				System.out.println(" a1 "+ a1.getID());
				List<IAtom> latoms = acred.getConnectedAtomsList(a1);
				for(IAtom a2 : latoms){
//					System.out.println(" a2l "+ a2.getID());
					if(!a2.getFlag(1)){
//						System.out.println(" a2 "+ a2.getID());
						IBond bond = acred.getBond(a1, a2);
//						bond.setProperty("Bond2frag", ""+a1.getID());
//						IBond newbond = new Bond();
//						IAtom[] atoms = new Atom[2];
//						atoms[0] = aR;
//						atoms[1] = a2;
//						newbond.setAtoms(atoms);
//						newbond.setOrder(bond.getOrder());
//						newbond.setProperty("Bond2frag", ""+a1.getID());
//						acred.addBond(bond);
						acred.addBond(acred.getAtomNumber(aR), acred.getAtomNumber(a2), bond.getOrder());
						IBond newbond = acred.getBond(aR, a2);
						newbond.setProperty("Bond2frag", ""+a1.getID());

					}
				}
			}
			
		}
		for(IAtom a1 : atomlist){

			if(a1.getFlag(1)){
//				System.out.println(" a3 "+ a1.getID());
				acred.removeAtomAndConnectedElectronContainers(a1);
			}
		}

		System.out.println(" atoms "+ acred.getAtomCount());
		System.out.println(" bonds "+ acred.getBondCount());

		return acred;
	}
}
