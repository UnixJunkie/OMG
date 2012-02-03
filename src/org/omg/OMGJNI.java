package org.omg;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;



public class OMGJNI {
	static {    
		// Get input stream from jar resource
	     InputStream inputStream;
		try {
			inputStream = OMGJNI.class.getResource("libnautygetcan.so").openStream();
			File temporarySO = File.createTempFile("libnautygetcan", ".so");
		     FileOutputStream outputStream = new FileOutputStream(temporarySO);
		     byte[] array = new byte[8192];
		     int read = 0;
		     while ( (read = inputStream.read(array)) > 0)
		         outputStream.write(array, 0, read);
		     outputStream.close();  

		     // Delete on exit the dll
		     temporarySO.deleteOnExit();  
		     System.load(temporarySO.getPath());

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}  
	}  

	public static native int[] getcanmultig(int fw, int[] arr, int[] lab, int[] ptnn1);
	public static native int[] getcanmultig2(int mc, int fw, int[] arr, int[] lab, int[] ptnn1);
	
	public OMGJNI(){
		
	}


	
}