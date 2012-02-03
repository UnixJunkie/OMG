#include "naututil.h"
#include "org_omg_OMGJNI.h"


#define MAXN 60
#define MAXNE 100
#define MAXNV 60



JNIEXPORT jintArray JNICALL Java_org_omg_OMGJNI_getcanmultig(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label, jintArray ptnn1){
	jint *carr, *clab,*cptn;
	int i,j,k,new,edges,vertices = fw;
	char c;
	int n = (vertices*4);
	int m = (n + WORDSIZE - 1) / WORDSIZE;
	graph g[n*n];
	graph gcanon[n*n];
	nvector lab[n], ptn[n], orbits[n];
	static DEFAULTOPTIONS(options);
	statsblk(stats);
	setword workspace[50*n];
	set *gv;
	options.writemarkers=FALSE;
	options.defaultptn=FALSE;
	options.getcanon=TRUE;
	options.writeautoms=FALSE;

	for(i=0; i<n; i++) {
		gv=GRAPHROW(g,i,m);
		EMPTYSET(gv, m);
	}
	carr = (*env)->GetIntArrayElements(env, arr, NULL);
	clab = (*env)->GetIntArrayElements(env, label, NULL);
	cptn = (*env)->GetIntArrayElements(env, ptnn1, NULL);

	edges=0;
	for(i=0; i<vertices; i++) {
		for(j=0; j<i; j++) {
			for(k=0; k<carr[i*vertices+j]; k++) {
				/*    create new edge-vertex #vertices+edges and add edges {i,new} {j,new} */
				new = vertices+edges;
				gv = GRAPHROW(g,i,m);

				ADDELEMENT(gv, new);
				gv = GRAPHROW(g,j,m);
				ADDELEMENT(gv, new);
				gv = GRAPHROW(g,new,m);
				ADDELEMENT(gv, i);
				ADDELEMENT(gv, j);

				edges++;
			}
		}
		lab[i] = clab[i];
		ptn[i] = cptn[i];
	}
	for(i=vertices; i<vertices+edges; i++) {
		lab[i]=i;
		ptn[i]=1;
	}
	ptn[vertices-1]=0;
	ptn[vertices+edges-1]=0;

	nauty(g, lab, ptn, NILSET, orbits, &options, &stats, workspace, 50*n, m, vertices+edges, gcanon);

	jintArray jgraph = (*env)->NewIntArray(env, (fw*fw)+fw);
	jint* pArray = (*env)->GetIntArrayElements(env, jgraph, 0);
	for(i=0; i<vertices; i++) {
		for(j=0; j<i; j++) {
			c='0';
			for(k=0; k<vertices+edges;k++) {
				if(i==k || j==k)
					continue;
				gv = GRAPHROW(gcanon,k,m);
				if (ISELEMENT(gv, i) && ISELEMENT(gv, j)) {
					c++;
				}
			}
			pArray[i*vertices+j]=c-'0';
			pArray[i+(vertices*j)]=c-'0';
		}
		pArray[(fw*fw)+i] = lab[i];
	}

	(*env)->ReleaseIntArrayElements(env, jgraph, pArray, 0);
	(*env)->ReleaseIntArrayElements(env, arr, carr, 0);
	(*env)->ReleaseIntArrayElements(env, label, clab, 0);
	(*env)->ReleaseIntArrayElements(env, ptnn1, cptn, 0);

	return jgraph;
}
JNIEXPORT jintArray JNICALL Java_org_omg_OMGJNI_getcanmultig2(JNIEnv *env, jobject thisobject,jint mc, jint fw, jintArray arr, jintArray label, jintArray ptnn1){
	jint *carr, *clab,*cptn;
	int i,j,k,new,edges,vertices = fw;
	char c;

	DYNALLSTAT(graph,g,g_sz);
	DYNALLSTAT(graph,gcanon,gcanon_sz);
	DYNALLSTAT(int,lab,lab_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);
	DYNALLSTAT(setword,workspace,workspace_sz);
	static DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;
	int n,m,v;
	set *gv;

	options.writemarkers=FALSE;
	options.defaultptn=FALSE;
	options.getcanon=TRUE;
	options.writeautoms=FALSE;


	n = vertices*vertices*4;
	m = (n + WORDSIZE - 1) / WORDSIZE;

	DYNALLOC2(graph,g,g_sz,m,n,"malloc");
	DYNALLOC2(graph,gcanon,gcanon_sz,m,n,"malloc");
	DYNALLOC1(setword,workspace,workspace_sz,50*m,"malloc");
	DYNALLOC1(int,lab,lab_sz,n,"malloc");
	DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
	DYNALLOC1(int,orbits,orbits_sz,n,"malloc");

	for(i=0; i<n; i++) {
		gv=GRAPHROW(g,i,m);
		EMPTYSET(gv, m);
	}
	carr = (*env)->GetIntArrayElements(env, arr, NULL);
	clab = (*env)->GetIntArrayElements(env, label, NULL);

	cptn = (*env)->GetIntArrayElements(env, ptnn1, NULL);

	edges=0;
	for(i=0; i<vertices; i++) {

		for(j=0; j<i; j++) {
			for(k=0; k<carr[i*vertices+j]; k++) {
				/*    create new edge-vertex #vertices+edges and add edges {i,new} {j,new} */
				new = vertices+edges;
				gv = GRAPHROW(g,i,m);

				ADDELEMENT(gv, (new+n-1)%n);
				gv = GRAPHROW(g,j,m);
				ADDELEMENT(gv, new);
				gv = GRAPHROW(g,new,m);
				ADDELEMENT(gv, (i+n-1)%n);
				ADDELEMENT(gv, (j+1)%n);

				edges++;
			}
		}
		lab[i] = clab[i];
		ptn[i] = cptn[i];
	}
	for(i=vertices; i<vertices+edges; i++) {
		lab[i]=i;
		ptn[i]=1;
	}
	ptn[vertices-1]=0;
	ptn[vertices+edges-1]=0;

	nauty(g, lab, ptn, NILSET, orbits, &options, &stats, workspace, 50*m, m, vertices+edges, gcanon);

	jintArray jgraph = (*env)->NewIntArray(env, (fw*fw)+fw);
	jint* pArray = (*env)->GetIntArrayElements(env, jgraph, 0);
	for(i=0; i<vertices; i++) {
		for(j=0; j<i; j++) {
			c='0';
			for(k=0; k<vertices+edges;k++) {
				if(i==k || j==k)
					continue;
				gv = GRAPHROW(gcanon,k,m);
				if (ISELEMENT(gv, i) && ISELEMENT(gv, j)) {
					c++;
				}
			}
			pArray[i*vertices+j]=c-'0';
			pArray[i+(vertices*j)]=c-'0';
		}
		pArray[(fw*fw)+i] = lab[i];
	}

	(*env)->ReleaseIntArrayElements(env, jgraph, pArray, 0);
	(*env)->ReleaseIntArrayElements(env, arr, carr, 0);
	(*env)->ReleaseIntArrayElements(env, label, clab, 0);
	(*env)->ReleaseIntArrayElements(env, ptnn1, cptn, 0);


	return jgraph;

}
