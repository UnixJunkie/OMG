#define MAXN 100
#include "nauty.h"
#include "Javanautygetcan.h"

JNIEXPORT jintArray JNICALL Java_Javanautygetcan_getcan(JNIEnv *env, jobject thisobject, jint fw, jintArray arr, jintArray label){
  jint *carr, *clab;
  static DEFAULTOPTIONS(options);
  options.defaultptn = FALSE;
  options.getcanon = TRUE;
  options.writeautoms = FALSE;
  statsblk(stats);
  setword workspace[50*MAXN];

  graph g[MAXN*MAXN];
  int lab[MAXN], ptn[MAXN], orbits[MAXN];

  int n,m;
  set *gv;


  n = 5;
  m = (n + WORDSIZE - 1) / WORDSIZE;
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

  gv = GRAPHROW(g, 0, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 1);


  gv = GRAPHROW(g, 1, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 0);
  ADDELEMENT(gv, 2);

  gv = GRAPHROW(g, 2, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 1);
  ADDELEMENT(gv, 3);

  gv = GRAPHROW(g, 3, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 2);
  ADDELEMENT(gv, 4);

  gv = GRAPHROW(g, 4, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 3);

  lab[0] = 0;
  lab[1] = 4;
  lab[2] = 1;
  lab[3] = 3;
  lab[4] = 2;

  ptn[0] = 1;
  ptn[1] = 0;
  ptn[2] = 1;
  ptn[3] = 0;
  ptn[4] = 0;




  graph g2[MAXN*MAXN];
  int lab2[MAXN], ptn2[MAXN];

  gv = GRAPHROW(g2, 0, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 1);
  ADDELEMENT(gv, 4);

  gv = GRAPHROW(g2, 1, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 0);

  gv = GRAPHROW(g2, 2, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 4);
  ADDELEMENT(gv, 3);

  gv = GRAPHROW(g2, 3, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 2);

  gv = GRAPHROW(g2, 4, m);
  EMPTYSET(gv, m);
  ADDELEMENT(gv, 0);
  ADDELEMENT(gv, 2);

  lab2[0] = 1;
  lab2[1] = 3;
  lab2[2] = 2;
  lab2[3] = 0;
  lab2[4] = 4;

  ptn2[0] = 1;
  ptn2[1] = 0;
  ptn2[2] = 1;
  ptn2[3] = 0;
  ptn2[4] = 0;



  graph canong[MAXN*MAXN];
  graph canong2[MAXN*MAXN];

  nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 50*MAXN,
m, n, canong);
  nauty(g2, lab2, ptn2, NULL, orbits, &options, &stats, workspace,
50*MAXN, m, n, canong2);

  if ( memcmp(canong,canong2,m*n*sizeof(graph)) ) {
    printf("Not identical\n");
  }
  else {
    printf("Identical\n");
  }
}
