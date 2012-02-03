package org.omg;

import java.io.FileNotFoundException;

import junit.framework.Assert;
import junit.framework.TestCase;

import org.openscience.cdk.exception.CDKException;

public class OMGTest extends TestCase {
	
	public void testC6H6() throws FileNotFoundException, CDKException, CloneNotSupportedException {
		int countTest = 217;
		OMG gen = new OMG();
		String formula = "C6H6";
		String fragments = null;
		String out = "testout_C6H6.sdf";
		gen.initializeMolecule(formula,fragments, out);
		
		int countValue = gen.getFinalCount();
		
		Assert.assertEquals(countTest, countValue);
	}
	public void testC4H6O5() throws FileNotFoundException, CDKException, CloneNotSupportedException {
		int countTest = 8070;
		OMG gen = new OMG();
		String formula = "C4H6O5";
		String fragments = null;
		String out = "testout_C4H6O5.sdf";
		gen.initializeMolecule(formula,fragments, out);
		
		int countValue = gen.getFinalCount();
		
		Assert.assertEquals(countTest, countValue);
	}
	public void testC9H12_benzene() throws FileNotFoundException, CDKException, CloneNotSupportedException {
		int countTest = 11;
		OMG gen = new OMG();
		String formula = "C9H12";
		String fragments = "fragment_benzene.sdf";
		String out = "testout_C9H12_benzene.sdf";
		gen.initializeMolecule(formula,fragments, out);
		
		int countValue = gen.getFinalCount();
		
		Assert.assertEquals(countTest, countValue);
	}
	public void testC10H14() throws FileNotFoundException, CDKException, CloneNotSupportedException {
		int countTest = 81909;
		OMG gen = new OMG();
		String formula = "C10H14";
		String fragments = null;
		String out = "testout_C10H14.sdf";
		gen.initializeMolecule(formula,fragments, out);
		
		int countValue = gen.getFinalCount();
		
		Assert.assertEquals(countTest, countValue);
	}
}
