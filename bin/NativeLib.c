#define MAXN 1000
#include "nauty.h" 
#include "JavaJNI.h"

JNIEXPORT jintArray JNICALL Java_JavaJNI_write(JNIEnv *env, jobject thisobject, jint fw, jintArray arr){
  jint *carr;
  graph g[MAXN*MAXN];
  int lab[MAXN],ptn[MAXN], orbits[MAXN];
  statsblk(stats);
  setword workspace[5*MAXN];

  int i,m,ll, n = fw;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  int k=0;
  static DEFAULTOPTIONS(options);
  options.writeautoms = TRUE;
  carr = (*env)->GetIntArrayElements(env, arr, NULL);
  for (k=0; k<fw*fw; ++k){
  	g[k] = carr[k];
  }
  nauty_check (WORDSIZE,m,n,NAUTYVERSIONID);
  nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,50*m,m,n,NULL);


   /*
   * Build a Java int[] with the Lab numbers
   */
   jintArray jArray = (*env)->NewIntArray(env, n);
   if ( n > 0) {
      jint* pArray = (*env)->GetIntArrayElements(env, jArray, 0);
	for ( i = 0; i < n; i++) {
		pArray[i] = lab[i];
	}
	(*env)->ReleaseIntArrayElements(env, jArray, pArray, 0);
  }

  /*
  * Return the LANA list
  */
   return jArray;   
}

