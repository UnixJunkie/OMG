#include "naututil.h"
#include "org_structgen_StructGenJNI.h"


JNIEXPORT jintArray JNICALL Java_org_structgen_StructGenJNI_getcan(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label, jintArray ptnn1){
  jint *carr, *clab,*cptn;

  graph g[fw*fw];
  int lab[fw],ptn[fw], orbits[fw];
  set  *gv;
  graph canong[fw*fw];
  statsblk(stats);
  setword workspace[5000*fw];
  int m, n = fw;
  int i,ii,ll, u, v;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  int k=0, j=0;
  static DEFAULTOPTIONS(options);
  options.defaultptn = FALSE;
  options.digraph = FALSE;
  options.getcanon = TRUE;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  nauty_check (WORDSIZE,m,n,NAUTYVERSIONID);

 // printf("Before nauty5\n");
  carr = (*env)->GetIntArrayElements(env, arr, NULL);
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


      //  printf("Before nauty6\n");
  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 500000*fw,
 m, n, canong);

 // printf("After nauty\n");
   FILE *outfile = fopen("out_graph1","w");
       fprintf(outfile,"Graph1\n");
       putgraph(outfile, g, 0, m, n);
       fprintf(outfile,"Canong1\n");
       putgraph(outfile, canong, 0, m, n);
       fclose(outfile);


   /*
   * Build a Java int[] with the Lab numbers
   */
   jintArray jgraph = (*env)->NewIntArray(env, fw*fw);
  // if ( n > 0) {
	   jint* pArray = (*env)->GetIntArrayElements(env, jgraph, 0);
	   /*int w;
	   for (w = 0; w < m*n; ++w){
		   		pArray[w] = canong[w];
	   }*/
	   for(i=0;i<n;i++){
	   	   gv=GRAPHROW(canong,i,m);
	   	   for(ii=0;ii<n;ii++){
	   		   if(ISELEMENT(gv,ii)){
	   			   pArray[i*n+ii] = 1;
	   		   }
	   		   else{
	   			   pArray[i*n+ii] = 0;
	   		   }
	   	   }
	      }
	   (*env)->ReleaseIntArrayElements(env, jgraph, pArray, 0);
   //}


   //printf("After nauty2\n");
   return jgraph;
}

JNIEXPORT jintArray JNICALL Java_org_structgen_StructGenJNI_getlab(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label, jintArray ptnn1){
  jint *carr, *clab,*cptn;

  graph g[fw*fw];
  int lab[fw],ptn[fw], orbits[fw];
  set  *gv;
  graph canong[fw*fw];
  statsblk(stats);
  setword workspace[5000*fw];
  int m, n = fw;
  int i,ii,ll, u, v;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  int k=0, j=0;
  static DEFAULTOPTIONS(options);
  options.defaultptn = FALSE;
  options.digraph = FALSE;
  options.getcanon = TRUE;
  //printf("Before nautyyyy\n");
  carr = (*env)->GetIntArrayElements(env, arr, NULL);
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
       // printf("Before nautyyyy2\n");
  nauty_check (WORDSIZE,m,n,NAUTYVERSIONID);
 // printf("Before nautyyyy3\n");
  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 5000*fw,
 m, n, canong);
 // printf("After nautyyyy\n");
   jintArray jlab = (*env)->NewIntArray(env, fw);
 //  if ( n > 0) {
         jint* pArray = (*env)->GetIntArrayElements(env, jlab, 0);
   	for ( i = 0; i < n; i++) {
   		pArray[i] = lab[i];
   	}
	   (*env)->ReleaseIntArrayElements(env, jlab, pArray, 0);
 //  }
   return jlab;
}
JNIEXPORT jintArray JNICALL Java_org_structgen_StructGenJNI_getptn(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label, jintArray ptnn1){
  jint *carr, *clab,*cptn;

  graph g[fw*fw];
  int lab[fw],ptn[fw], orbits[fw];
  set  *gv;
  graph canong[fw*fw];
  statsblk(stats);
  setword workspace[5000*MAXN];
  int m, n = fw;
  int i,ii,ll, u, v;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  int k=0, j=0;
  static DEFAULTOPTIONS(options);
  options.defaultptn = FALSE;
  options.digraph = FALSE;
  options.getcanon = TRUE;
 // printf("Before nautyyyyoooo\n");
  carr = (*env)->GetIntArrayElements(env, arr, NULL);
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

  nauty_check (WORDSIZE,m,n,NAUTYVERSIONID);
  //printf("Before nautyyyyoooo2\n");
  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 5000*fw,
 m, n, canong);
 // printf("After nautyyyyoooo\n");
   jintArray jptn = (*env)->NewIntArray(env, fw);
//   if ( n > 0) {
         jint* pArray = (*env)->GetIntArrayElements(env, jptn, 0);
   	for ( i = 0; i < n; i++) {
   		pArray[i] = ptn[i];
   	}
	   (*env)->ReleaseIntArrayElements(env, jptn, pArray, 0);
  // }
 //  printf("After nautyyyyoooo2\n");
   return jptn;
}

JNIEXPORT jintArray JNICALL Java_org_structgen_StructGenJNI_getcanmultig(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label, jintArray ptnn1){
jint *carr, *clab,*cptn;
int i,j,k,m,new,edges,vertices = fw;
  
  char c;

  graph g[vertices*vertices*4*4]; 
  graph gcanon[vertices*vertices]; 
  nvector lab[vertices], ptn[vertices], orbits[vertices];
  static DEFAULTOPTIONS(options);
  statsblk(stats);
  setword workspace[5000*vertices];
  set *gv;

  options.writemarkers=FALSE;
  options.defaultptn=FALSE;
  options.getcanon=TRUE;
  options.writeautoms=FALSE;
m = (vertices*4 + WORDSIZE - 1) / WORDSIZE;
  nauty_check (WORDSIZE,m,vertices,NAUTYVERSIONID);
  
  for(i=0; i<vertices*4; i++) {
    gv=GRAPHROW(g,i,vertices*4);
    EMPTYSET(gv, vertices*4);
  }
  carr = (*env)->GetIntArrayElements(env, arr, NULL);
	/*car[0]=0; car[1]=1; car[2]=0; car[3]=0; 
	car[4]=1; car[5]=0; car[6]=0; car[7]=0; 
	car[8]=0; car[9]=0; car[10]=0; car[11]=2; 
	car[12]=0; car[13]=0; car[14]=2; car[15]=0; */
	
	clab = (*env)->GetIntArrayElements(env, label, NULL);
    for (j=0; j<vertices; ++j){
    	lab[j] = clab[j];
    }

    cptn = (*env)->GetIntArrayElements(env, ptnn1, NULL);
        for (j=0; j<vertices; ++j){
        	ptn[j] = cptn[j];
        }

  edges=0;
  for(i=0; i<vertices; i++) {
    for(j=0; j<i; j++) {
      for(k=0; k<carr[i*vertices+j]; k++) {
        /*    create new edge-vertex #vertices+edges and add edges {i,new} {j,new} */
        new = vertices+edges;
        gv = GRAPHROW(g,i,MAXM);
        ADDELEMENT(gv, new);
        gv = GRAPHROW(g,j,MAXM);
        ADDELEMENT(gv, new);
        gv = GRAPHROW(g,new,MAXM);
        ADDELEMENT(gv, i);
        ADDELEMENT(gv, j);
 
        edges++;
      }
      //s++;
    }
  }

 /* for(i=0; i<vertices+edges; i++) {
    lab[i]=i;
    ptn[i]=1;
  }
	lab[0]=0;
	lab[1]=1;
	lab[2]=2;
	lab[3]=3;
	ptn[0]=1;
	ptn[1]=0;
	ptn[2]=1;
	ptn[3]=0;*/
  ptn[vertices-1]=0;
  ptn[vertices+edges-1]=0;
 /* printf("Before nauty\n");
  printf("Lab[0]: %d\n", lab[0]);
  printf("Lab[1]: %d\n", lab[1]);
  printf("Lab[2]: %d\n", lab[2]);
  printf("Lab[3]: %d\n", lab[3]);

  printf("\n");
  printf("ptn[0]: %d\n", ptn[0]);
  printf("ptn[1]: %d\n", ptn[1]);
  printf("ptn[2]: %d\n", ptn[2]);
  printf("ptn[3]: %d\n", ptn[3]);*/

  nauty(g, lab, ptn, NILSET, orbits, &options, &stats, workspace, 50*MAXM, MAXM, vertices+edges, gcanon);
/*	printf("After nauty\n");
  printf("Lab[0]: %d\n", lab[0]);
  printf("Lab[1]: %d\n", lab[1]);
  printf("Lab[2]: %d\n", lab[2]);
  printf("Lab[3]: %d\n", lab[3]);

  printf("\n");
  printf("ptn[0]: %d\n", ptn[0]);
  printf("ptn[1]: %d\n", ptn[1]);
  printf("ptn[2]: %d\n", ptn[2]);
  printf("ptn[3]: %d\n", ptn[3]);*/
jintArray jgraph = (*env)->NewIntArray(env, fw*fw);
  // if ( n > 0) {
	   jint* pArray = (*env)->GetIntArrayElements(env, jgraph, 0);
  for(i=0; i<vertices; i++) {
    for(j=0; j<i; j++) {
      c='0';
      for(k=0; k<vertices+edges;k++) {
        if(i==k || j==k)
          continue;
        gv = GRAPHROW(gcanon,k,MAXM);
        if (ISELEMENT(gv, i) && ISELEMENT(gv, j)) {
          c++;
        }
      }
      //*(canon++)=c;
	  pArray[i*vertices+j]=c-'0';
	  pArray[i+(vertices*j)]=c-'0';
    }
  }
  /**canon=0;
  printf(" %d %d %d %d\n", res[0],res[1],res[2],res[3]);
  printf(" %d %d %d %d\n", res[4],res[5],res[6],res[7]);
  printf(" %d %d %d %d\n", res[8],res[9],res[10],res[11]);
  printf(" %d %d %d %d\n", res[12],res[13],res[14],res[15]);
	return canon;*/
	(*env)->ReleaseIntArrayElements(env, jgraph, pArray, 0);
	return jgraph;
}

JNIEXPORT jintArray JNICALL Java_org_structgen_StructGenJNI_getcan2(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label){
  jint *carr, *clab;
  graph g[fw*fw];
  int lab[fw],ptn[fw], orbits[fw];
  graph canong[fw*fw];
  statsblk(stats);
  setword workspace[5000*fw];

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
  nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,5000*fw,m,n,canong);


   /*
   * Build a Java int[] with the canonical graph
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
JNIEXPORT jintArray JNICALL Java_org_structgen_StructGenJNI_getcan3(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label, jintArray ptnn1, jintArray arr2, jintArray label2, jintArray ptnn2){
  jint *carr, *clab,*cptn, *carr2, *clab2,*cptn2;

  graph g[fw*fw], g2[fw*fw];
  int lab[fw],ptn[fw], orbits[fw];
  int lab2[fw],ptn2[fw], orbits2[fw];

  set  *gv;
  graph canong[fw*fw];
  graph canong2[fw*fw];
  statsblk(stats);
  setword workspace[5000*fw];

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

  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 5000*fw,
 m, n, canong);

   nauty(g2, lab2, ptn2, NULL, orbits, &options, &stats, workspace,
 5000*fw, m, n, canong2);

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
