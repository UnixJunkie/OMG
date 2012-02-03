package org.omg;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;
import junit.framework.TestCase;

import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.nonotify.NNAtom;
import org.openscience.cdk.nonotify.NNAtomType;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

public class AtomsPNTest extends  TestCase {
	public static Map<String, Integer> testedAtomTypes = new HashMap<String, Integer>();
	
	public void testP1() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);

//		IAtomContainer ac = atomP.getBuilder().newAtomContainer();
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
		}
	public void testP2() throws FileNotFoundException, CDKException, CloneNotSupportedException {
	IAtomContainer acontainer = MolecularFormulaManipulator.getAtomContainer(
			MolecularFormulaManipulator.getMolecularFormula("H3PO4", DefaultChemObjectBuilder.getInstance()));
//	IAtomContainer acprotonate = (IAtomContainer) acontainer.clone();
	IAtomType type = null;
	System.out.println(acontainer);
	int nH = 0;
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
	for (IAtom atom : acontainer.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(acontainer.getBuilder()).findMatchingAtomType(acontainer, atom);
		System.out.println(type);
		Assert.assertNotNull(type);
	}
	}
	public void test_Oxophosphanyl() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);
		
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		ac.addAtom(new Atom("O"));
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.DOUBLE);

		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
	}
	public void test_Methinophosphide() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);
		
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		ac.addAtom(new Atom("C"));
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.TRIPLE);

		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
	}
	public void test_P2double() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);
		
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		ac.addAtom(new Atom("C"));
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.DOUBLE);

		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 2, IBond.Order.DOUBLE);
		
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
	}
	public void test_Pphosphate() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);
		
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		ac.addAtom(new Atom("O"));
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.DOUBLE);

		ac.addAtom(new Atom("O"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 2, IBond.Order.SINGLE);
		
		ac.addAtom(new Atom("O"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 3, IBond.Order.SINGLE);
		
		ac.addAtom(new Atom("O"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 4, IBond.Order.SINGLE);
		
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
	}
	public void test_P2double_1single() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);
		
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		ac.addAtom(new Atom("C"));
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.DOUBLE);

		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 2, IBond.Order.DOUBLE);
		
		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 3, IBond.Order.SINGLE);
		
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
	}
	public void test_Ptriple_double() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);
		
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		ac.addAtom(new Atom("C"));
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.TRIPLE);

		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 2, IBond.Order.DOUBLE);
		
		
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
	}
	public void test_P5single() throws FileNotFoundException, CDKException {
		IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
		IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
		AtomTypeManipulator.configure(atomP, atomTypeP);
		
		IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomP);
		ac.addAtom(new Atom("C"));
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.SINGLE);

		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 2, IBond.Order.SINGLE);
		
		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 3, IBond.Order.SINGLE);
		
		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 4, IBond.Order.SINGLE);
		
		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(0);
		ac.addBond(0, 5, IBond.Order.SINGLE);
		
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
		type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
		Assert.assertNotNull(type);
		}
	}
	
	public void testN() throws CDKException, IOException, CloneNotSupportedException {
		IAtom atomN = new NNAtom(Elements.NITROGEN);
		IAtomType atomTypeN = new NNAtomType(Elements.NITROGEN);
		AtomTypeManipulator.configure(atomN, atomTypeN);
		
		IAtomContainer ac = (IAtomContainer) atomN.clone();
		ac.addAtom(atomN);
		
		ac.getAtom(0).setImplicitHydrogenCount(0);
		
		ac.addAtom(new Atom("C"));
		ac.getAtom(1).setImplicitHydrogenCount(2);
		ac.addBond(0, 1, IBond.Order.DOUBLE);
		
		ac.addAtom(new Atom("C"));
		ac.getAtom(2).setImplicitHydrogenCount(3);
		ac.addBond(0, 2, IBond.Order.SINGLE);
		
		ac.addAtom(new Atom("O"));
		ac.getAtom(3).setImplicitHydrogenCount(0);
		ac.addBond(0, 3, IBond.Order.DOUBLE);
		

		
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
			type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
			AtomTypeManipulator.configure(atom, type);
			Assert.assertNotNull(type);
		}
		SmilesGenerator sg = new SmilesGenerator();
		System.out.println(sg.createSMILES(new Molecule(ac)));
	}
	public void testC() throws CDKException, IOException, CloneNotSupportedException {
		IAtom atomN = new NNAtom(Elements.CARBON);
		IAtomType atomTypeN = new NNAtomType(Elements.CARBON);
		AtomTypeManipulator.configure(atomN, atomTypeN);
		
		IAtomContainer ac = atomN.getBuilder().newInstance(IAtomContainer.class);
		ac.addAtom(atomN);
		
		ac.getAtom(0).setImplicitHydrogenCount(0);
		
		ac.addAtom((IAtom)atomN.clone());
		ac.getAtom(1).setImplicitHydrogenCount(0);
		ac.addBond(0, 1, IBond.Order.TRIPLE);
		
	
		IAtomType type = null;
		for (IAtom atom : ac.atoms()) {
			type = CDKAtomTypeMatcher.getInstance(ac.getBuilder()).findMatchingAtomType(ac, atom);
			AtomTypeManipulator.configure(atom, type);
			Assert.assertNotNull(type);
		}
		SmilesGenerator sg = new SmilesGenerator();
		System.out.println(sg.createSMILES(new Molecule(ac)));
	}
    public void testPine() throws Exception {
    	IAtom atomP = new NNAtom(Elements.PHOSPHORUS);
    	IAtomType atomTypeP = new NNAtomType(Elements.PHOSPHORUS);
    	AtomTypeManipulator.configure(atomP, atomTypeP);

//    	IAtomContainer ac = atomP.getBuilder().newAtomContainer();
    	IAtomContainer ac = atomP.getBuilder().newInstance(IAtomContainer.class);
    	ac.addAtom(atomP);
    	IAtomType type = null;
    	for (IAtom atom : ac.atoms()) {
    		type = CDKAtomTypeMatcher.getInstance(
    				ac.getBuilder()
    			).findMatchingAtomType(ac, atom);
    		Assert.assertNotNull(type);
    	}
    }
//    public void testP() throws Exception {
//       	IAtom atomP = new NNAtom("P");
//       	IAtomContainer mol = new Molecule();
//       	mol.addAtom(atomP);
//        String[] expectedTypes = {"P.ine"};
//        assertAtomTypes(testedAtomTypes, expectedTypes, mol);
//    }
}
