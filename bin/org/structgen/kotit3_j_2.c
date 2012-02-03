#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MAXN 60
#include "nauty.h"

void canon(int vertices, char *s, char *canon) {
  int i,j,k,new,edges;
  int car[vertices*vertices], res[vertices*vertices];
  char c;

  graph g[MAXN*MAXM]; 
  graph gcanon[MAXN*MAXM]; 
  nvector lab[MAXN], ptn[MAXN], orbits[MAXN];
  static DEFAULTOPTIONS(options);
  statsblk(stats);
  setword workspace[50*MAXM];
  set *gv;

  options.writemarkers=FALSE;
  options.defaultptn=FALSE;
  options.getcanon=TRUE;
  options.writeautoms=FALSE;

  for(i=0; i<MAXN; i++) {
    gv=GRAPHROW(g,i,MAXM);
    EMPTYSET(gv, MAXM);
  }
	car[0]=0; car[1]=1; car[2]=0; car[3]=0; 
	car[4]=1; car[5]=0; car[6]=0; car[7]=0; 
	car[8]=0; car[9]=0; car[10]=0; car[11]=2; 
	car[12]=0; car[13]=0; car[14]=2; car[15]=0; 
  edges=0;
  for(i=0; i<vertices; i++) {
    for(j=0; j<i; j++) {
      for(k=0; k<car[i*vertices+j]; k++) {
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

  for(i=0; i<vertices+edges; i++) {
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
	ptn[3]=0;
  ptn[vertices-1]=0;
  ptn[vertices+edges-1]=0;
  printf("Before nauty\n");
  printf("Lab[0]: %d\n", lab[0]);
  printf("Lab[1]: %d\n", lab[1]);
  printf("Lab[2]: %d\n", lab[2]);
  printf("Lab[3]: %d\n", lab[3]);

  printf("\n");
  printf("ptn[0]: %d\n", ptn[0]);
  printf("ptn[1]: %d\n", ptn[1]);
  printf("ptn[2]: %d\n", ptn[2]);
  printf("ptn[3]: %d\n", ptn[3]);

  nauty(g, lab, ptn, NILSET, orbits, &options, &stats, workspace, 50*MAXM, MAXM, vertices+edges, gcanon);
	printf("After nauty\n");
  printf("Lab[0]: %d\n", lab[0]);
  printf("Lab[1]: %d\n", lab[1]);
  printf("Lab[2]: %d\n", lab[2]);
  printf("Lab[3]: %d\n", lab[3]);

  printf("\n");
  printf("ptn[0]: %d\n", ptn[0]);
  printf("ptn[1]: %d\n", ptn[1]);
  printf("ptn[2]: %d\n", ptn[2]);
  printf("ptn[3]: %d\n", ptn[3]);

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
      *(canon++)=c;
	  res[i*vertices+j]=c-'0';
	  res[i+(vertices*j)]=c-'0';
    }
  }
  *canon=0;
  printf(" %d %d %d %d\n", res[0],res[1],res[2],res[3]);
  printf(" %d %d %d %d\n", res[4],res[5],res[6],res[7]);
  printf(" %d %d %d %d\n", res[8],res[9],res[10],res[11]);
  printf(" %d %d %d %d\n", res[12],res[13],res[14],res[15]);
}

int main(int argc, char**argv)
{
  char buffer[1024];
  canon(4,"130002",buffer);
  printf("%s\n",buffer);
  return 0;
}
