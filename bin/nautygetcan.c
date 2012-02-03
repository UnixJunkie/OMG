#define MAXN 1000
#include "naututil.h"
#include "Javanautygetcan.h"

JNIEXPORT jintArray JNICALL Java_Javanautygetcan_getcan(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label){
  jint *carr, *clab;
  graph g[MAXN*MAXN];
  int lab[MAXN],ptn[MAXN], orbits[MAXN];
  graph canong[MAXN*MAXN];
  statsblk(stats);
  setword workspace[160*MAXN];

  int i,m,ll, n = fw;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  int k=0, j=0;
  static DEFAULTOPTIONS(options);
  options.defaultptn = TRUE;
  options.digraph = FALSE;
  options.getcanon = TRUE;

  carr = (*env)->GetIntArrayElements(env, arr, NULL);
  for (k=0; k<fw*fw; ++k){
  	g[k] = carr[k];
  }
  clab = (*env)->GetIntArrayElements(env, label, NULL);
    for (j=0; j<fw; ++j){
    	lab[j] = clab[j];
    }
  nauty_check (WORDSIZE,m,n,NAUTYVERSIONID);
  nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,160*MAXN,m,n,canong);


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
JNIEXPORT jintArray JNICALL Java_Javanautygetcan_getcan3(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label, jintArray ptnn1, jintArray arr2, jintArray label2, jintArray ptnn2){
  jint *carr, *clab,*cptn, *carr2, *clab2,*cptn2;

  graph g[MAXN*MAXN], g2[MAXN*MAXN];
  int lab[MAXN],ptn[MAXN], orbits[MAXN];
  int lab2[MAXN],ptn2[MAXN], orbits2[MAXN];

  set  *gv;
  graph canong[MAXN*MAXN];
  graph canong2[MAXN*MAXN];
  statsblk(stats);
  setword workspace[5000*MAXN];

  int i,m,ll, u, v, n = fw;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  int k=0, j=0;
  static DEFAULTOPTIONS(options);
  options.defaultptn = FALSE;
    options.digraph = FALSE;
    options.getcanon = TRUE;
  carr = (*env)->GetIntArrayElements(env, arr, NULL);
  /*for (k=0; k<fw*fw; ++k){
  	g[k] = carr[k];
  }*/
  for(v=0;v<n;v++){
	  gv=GRAPHROW(g,v,m);
	  EMPTYSET(gv,m);
	  for(u=0;u<n;u++){
		  if(carr[v*n+u]==1)
			  ADDELEMENT(gv,u);
	  }
  }
  clab = (*env)->GetIntArrayElements(env, label, NULL);
    for (j=0; j<fw; ++j){
    	lab[j] = clab[j];
    }

    cptn = (*env)->GetIntArrayElements(env, ptnn1, NULL);
        for (j=0; j<fw; ++j){
        	ptn[j] = cptn[j];
        }
    carr2 = (*env)->GetIntArrayElements(env, arr2, NULL);
    for(v=0;v<n;v++){
  	  gv=GRAPHROW(g2,v,m);
  	  EMPTYSET(gv,m);
  	  for(u=0;u<n;u++){
  		  if(carr2[v*n+u]==1)
  			  ADDELEMENT(gv,u);

  	  }

    }

    /* for (k=0; k<fw*fw; ++k){
     	g2[k] = carr2[k];
     }*/
     clab2 = (*env)->GetIntArrayElements(env, label2, NULL);
       for (j=0; j<fw; ++j){
       	lab2[j] = clab2[j];
       }
     cptn2 = (*env)->GetIntArrayElements(env, ptnn2, NULL);
               for (j=0; j<fw; ++j){
               	ptn2[j] = cptn2[j];
     }
  nauty_check (WORDSIZE,m,n,NAUTYVERSIONID);
  printf("Before nauty\n");
  printf("Lab1[0]: %d\n", lab[0]);
    printf("Lab1[1]: %d\n", lab[1]);
    printf("Lab1[2]: %d\n", lab[2]);
    printf("Lab1[3]: %d\n", lab[3]);
    printf("Lab1[4]: %d\n", lab[4]);

    printf("\n");
    printf("Lab2[0]: %d\n", lab2[0]);
     printf("Lab2[1]: %d\n", lab2[1]);
     printf("Lab2[2]: %d\n", lab2[2]);
     printf("Lab2[3]: %d\n", lab2[3]);
     printf("Lab2[4]: %d\n", lab2[4]);

  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 5000*MAXN,
 m, n, canong);

   nauty(g2, lab2, ptn2, NULL, orbits, &options, &stats, workspace,
 5000*MAXN, m, n, canong2);

   FILE *outfile = fopen("out_graph1","w");
       fprintf(outfile,"Graph1\n");
       putgraph(outfile, g, 0, m, n);
       fprintf(outfile,"Canong1\n");
       putgraph(outfile, canong, 0, m, n);
       fprintf(outfile,"Graph2\n");
       putgraph(outfile, g2, 0, m, n);
       fprintf(outfile,"Canong2\n");
       putgraph(outfile, canong2, 0, m, n);
       fclose(outfile);

   int w;
   for (w = 0; w < m*n; ++w)
         if (canong[w] != canong2[w]) break;
     if (w < m*n) {printf(" NOT ISOMORPHIC\n"); } else { printf("ISOMORPHIC\n"); }
   if ( memcmp(canong,canong2,m*n*sizeof(graph)) ) {
     printf("Not identical\n");
   }
   else {
     printf("Identical\n");
   }
   int samerows = 4;
   printf("testcanlab1  %d\n", testcanlab(g,canong,clab,&samerows,m,n));
   printf("testcanlab2  %d\n", testcanlab(g2,canong2,clab2,&samerows,m,n));

   printf("Lab1[0]: %d\n", lab[0]);
   printf("Lab1[1]: %d\n", lab[1]);
   printf("Lab1[2]: %d\n", lab[2]);
   printf("Lab1[3]: %d\n", lab[3]);
   printf("Lab1[4]: %d\n", lab[4]);

   printf("\n");
   printf("Lab2[0]: %d\n", lab2[0]);
    printf("Lab2[1]: %d\n", lab2[1]);
    printf("Lab2[2]: %d\n", lab2[2]);
    printf("Lab2[3]: %d\n", lab2[3]);
    printf("Lab2[4]: %d\n", lab2[4]);

    printf("ptn2[0]: %d\n", ptn2[0]);
         printf("ptn2[1]: %d\n", ptn2[1]);
         printf("ptn2[2]: %d\n", ptn2[2]);
         printf("ptn2[3]: %d\n", ptn2[3]);
         printf("ptn2[4]: %d\n", ptn2[4]);
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
