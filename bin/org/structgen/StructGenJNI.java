package org.structgen;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;



public class StructGenJNI {
	static {    
		// Get input stream from jar resource
	     InputStream inputStream;
		try {
			inputStream = StructGenJNI.class.getResource("libnautygetcan.so").openStream();
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

	     // Copy resource to filesystem in a temp folder with a unique name
	     

		//System.out.println("Loading" + System.mapLibraryName("nautygetcan"));
		//System.loadLibrary("nautygetcan");
		
	}  
  
	public static native int[] getcan(int fw, int[] arr, int[] lab, int[] ptnn1);
	public static native int[] getlab(int fw, int[] arr, int[] lab, int[] ptnn1);
	public static native int[] getptn(int fw, int[] arr, int[] lab, int[] ptnn1);
	public static native int[] getcanmultig(int fw, int[] arr, int[] lab, int[] ptnn1);
	public static native int[] getlabmultig(int fw, int[] arr, int[] lab, int[] ptnn1);
	public static native int[] getcan2(int fw, int[] arr, int[] lab);
	public static native int[] getcan3(int fw, int[] arr, int[] lab, int[] ptnn1 ,int[] arr2, int[] lab2, int[] ptnn2);
	
	public StructGenJNI(){
		
	}
	public class graph{
		private int[] arr;
		private int[] lab;
	
		public graph(int[] arr, int lab[]){
			this.arr = arr;
			this.lab = lab;
		}
	}
}