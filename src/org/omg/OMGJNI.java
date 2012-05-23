package org.omg;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;



public class OMGJNI {
	static {    
		// Get input stream from jar resource
	     InputStream inputStream = null;
	     File temporarySO = null;
		try {
			String OS = System.getProperty("os.name");
			
			String osName = System.getProperty("os.name").toLowerCase();
			String osBits = System.getProperty("sun.arch.data.model");
		    if (osName.contains("windows")) {

		    } else if (osName.contains("mac")) {
		    	if(osBits.equals("64")){
					inputStream = OMGJNI.class.getResource("libnautygetcanMac64.so").openStream();
					temporarySO = File.createTempFile("libnautygetcanMac64", ".so");
				}else if(osBits.equals("32")){
					inputStream = OMGJNI.class.getResource("libnautygetcanMac32.so").openStream();
					temporarySO = File.createTempFile("libnautygetcanMac32", ".so");
				}
		    }else if (osName.contains("linux")){
		    	if(osBits.equals("64")){
					inputStream = OMGJNI.class.getResource("libnautygetcanLinux64.so").openStream();
					temporarySO = File.createTempFile("libnautygetcanLinux64", ".so");
				}else if(osBits.equals("32")){
					inputStream = OMGJNI.class.getResource("libnautygetcanLinux32.so").openStream();
					temporarySO = File.createTempFile("libnautygetcanLinux32", ".so");
				}
		    }
			
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