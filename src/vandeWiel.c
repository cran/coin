
/**
    Exact Distribution of Two-Sample Permutation Tests
    van de Wiel split-up Algorithm

    Author: Mark van de Wiel (2001-2005) <m.a.v.d.wiel@TUE.nl>
            with modifications for R by 
            Torsten Hothorn <Torsten.Hothorn@R-project.org>
      
    *\file vandeWiel.c
    *\author $Author: thothorn $
    *\date $Date: 2012-03-15 10:32:24 +0100 (Thu, 15 Mar 2012) $
*/
                    
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
                    

/**
    The probability distribution for the independent two sample problem.
    
    REFERENCE
    
    M.A. van de Wiel (2001). The split-up algorithm: a fast symbolic method 
    for computing p-values of rank statistics, Computational Statistics, 
    16, 519-538.
    
*/
                                        
typedef struct {
    long length;
    double *c;
    double *x;
} celW;

double binomi(int m, int n) { 

    double bin = 1;
    double bin1 = 1;
    double bin2 = 1;
    int i, j;
        
    for (i = 1; i <= n; i++) bin1 = bin1 * (m + 1 -i);
    for (j = 1; j <= n; j++) bin2 = bin2 * j;
    bin = bin1/bin2;
    
    return(bin);
}

celW** reserveW(int a, int b)
{
    long res = 0;
    int i, j;
    celW** W;
    
    /* <FIXME>: need to free memory in case Calloc barfs
       Writing R Extension advertises .on.exit but
       I still need a pointer to the memory.
       </FIXME>
    */

    W = Calloc(a + 1, celW*);

    for (i = 0; i <= a; i++)
        W[i] = Calloc(b + 1, celW);
        
    for (i = 0; i <= a; i++) {
        for (j = i; j <= b; j++) {
            res = (long) binomi(j,i);
            /* the majority of memory is freed on exit and error
               thanks to S_alloc */
            W[i][j].c = (double *) S_alloc(res, sizeof(double));
            W[i][j].x = (double *) S_alloc(res, sizeof(double));
        }
        R_CheckUserInterrupt();
    }
    return(W);
}

void FreeW(int a, celW **W)
{
     int i;
 
     for (i = a; i >= 0; i--)
         Free(W[i]);
         
     Free(W);
}

void initW(int a, int b, celW **W) {

    int i, j;

    for (i = 1; i <= a; i++)
    for (j = 0; j <= b; j++) {
        W[i][j].length = 0;
    }
    for (j = 0; j <= b; j++) {
        W[0][j].length = 1;
        W[0][j].c[0] = 1;  
        W[0][j].x[0] = 0;  
    }
}

void mult(celW *tem, int a, int b, int rank, double *rs) {

    /*

    mult multiplies the polynomial c_i*x^(l_i) by x^(rs[rank]), 
    which means adding the exponents

    */

    int j;
    for (j = 0; j < tem[0].length; j++)
        tem[0].x[j] += rs[rank];
}

void plus(celW **W, celW *tempie, int a, int b, double tol) {

    /*

    plus adds terms with the same exponents after multiplication 
    with 1 + x^(rs[rank]), so c1*x^j + c2*x^j becomes (c1+c2)*x^j

    */

    int elep = 0;
    int k = 0;
    int test = 1;
    int i, j;
    
    for (i = 0; i < W[a][b-1].length; i++) {

        test = 1;
        
        for (j = elep; j < tempie[0].length && test==1; j++) {
        
            if (tempie[0].x[j] - tol <= W[a][b-1].x[i]
                && W[a][b-1].x[i] <= tempie[0].x[j] + tol) {

                tempie[0].c[j] += W[a][b-1].c[i];
                test = 0;
                elep = j;             
            }
        }
         
        if (test == 1) {
            tempie[0].c[tempie[0].length + k] = W[a][b-1].c[i];
            tempie[0].x[tempie[0].length + k] = W[a][b-1].x[i];
            k++;
        }
        R_CheckUserInterrupt();
    }
    tempie[0].length += k;
}

void mymergesort(celW temptw, long tijd)
{

    /*
    
    mymergesort composes one sorted list (increasing exponents of 
    the polynomial) from two separately sorted lists. c1*x^3 + c2*x^5 
    and  c3*x^4 + c4*x^7  becomes  c1*x^3 + c3*x^4 + c2*x^5 + c4*x^7.

    */

    celW copiep;
    int t1 = 0; 
    int t2 = 0;
    int i, j;

    copiep.c = Calloc(temptw.length, double);
    copiep.x = Calloc(temptw.length, double);
        
    for (i = 0; i < temptw.length; i++) {
        copiep.c[i] = temptw.c[i];
        copiep.x[i] = temptw.x[i];
    }
    
    for (j = 0; j < temptw.length; j++) {
        if (t1 <= tijd-1 && t2 <= temptw.length - tijd - 1) {
            if (copiep.x[t1] < copiep.x[tijd + t2]) {
                temptw.x[j] = copiep.x[t1];
                temptw.c[j] = copiep.c[t1];
                t1++;
            } else {
                temptw.x[j] = copiep.x[tijd + t2];
                temptw.c[j] = copiep.c[tijd + t2];
                t2++;
            }
        } else {
            if (t1 > tijd - 1) {
                temptw.x[j] = copiep.x[tijd + t2];
                temptw.c[j] = copiep.c[tijd + t2];
                t2++; 
            } else {   
                temptw.x[j] = copiep.x[t1];
                temptw.c[j] = copiep.c[t1];
                t1++;
            }
        }   
        R_CheckUserInterrupt();       
    } 
    Free(copiep.c);
    Free(copiep.x);
}

void fillcell(celW **W, int i1, int j1, int r, double *rs, double tol) {
    
    /*

    fillcell makes the new recursive polynomial W[i1][j1] from 
    W[i1 - 1][j1 - 1] and W[i1][j1 - 1]. j1 is the total number of 
    rank scores assigned so far to either of the two groups, 
    i1 is the number of rank scores assigned to the smallest sample.

    */

    long tijd;
    celW temp2;
    int j, j2;

    temp2.c = Calloc(W[i1 - 1][j1 - 1].length + 
                     W[i1][j1 - 1].length, double);
    temp2.x = Calloc(W[i1 - 1][j1 - 1].length + 
                     W[i1][j1 - 1].length, double);
    temp2.length = W[i1 - 1][j1 - 1].length;

    for (j = 0; j < temp2.length; j++) {
       temp2.c[j] = W[i1 - 1][j1 - 1].c[j];
       temp2.x[j] = W[i1 - 1][j1 - 1].x[j];
    }

    if (i1 == j1) {       
        mult(&temp2, i1 - 1, j1 - 1, r, rs); 
    } else {           
        mult(&temp2, i1 - 1, j1 - 1, r, rs);                        
        tijd = temp2.length;                                
        plus(W, &temp2, i1, j1, tol);                            
        mymergesort(temp2, tijd);                              
    }

    W[i1][j1].length = temp2.length;

    for (j2 = 0; j2 < temp2.length; j2++) {
        W[i1][j1].c[j2] = temp2.c[j2];
        W[i1][j1].x[j2] = temp2.x[j2];
    }          

    Free(temp2.c);
    Free(temp2.x);
}

void mirrorW(celW **W,int ce, int bep, int start, double *rs) {   

    /* 

    mirrorW contains a trick to speed op computations considerably. 
    By symmetry arguments it is easy to find W[i][tot] from W[tot-i][tot]

    */

    double totsum = 0;
    long len;
    int r, h;
    
    for (r = 0; r < bep; r++) totsum += rs[start + r];
    
    len = W[bep-ce][bep].length;
        
    for (h = 0; h < len; h++) {
        W[ce][bep].length = W[bep-ce][bep].length;
        W[ce][bep].c[len-1-h] = W[bep-ce][bep].c[h];
        W[ce][bep].x[len-1-h] = totsum - W[bep-ce][bep].x[h];
    }
}

void makeW(celW **W, int a, int b, int start, double *rs, double tol) {


    /* 

    makeW simply determines whether a new polynomial W[i][j] 
    can be found from mirrorW (if W[j-i][j] is available) or needs 
    to be constructed via multiplication etc.

    */

    long i,j;
    int rank;
    int hulp;

    for (j = 1; j <= b; j++) {  /* verander naar 0!! */

        if (j < a) {
            hulp = j; 
        } else {
            hulp = a;
        }

        for (i=1; i <= hulp; i++) {   
            if (i <= j/2 || j == 1) {
                rank = start+j;
                fillcell(W, i, j, rank - 1, rs, tol);
            } else {
                mirrorW(W, i, j, start, rs);
            }                               
            R_CheckUserInterrupt();
        }
    }
}

void cumulcoef(celW **W, int i1, int j1) {


    /*

    cumulcoef recursively adds the coefficients of the 
    sorted polynomial. So, 3*x^4 + 4*x^6 + 2*x^7  becomes  
    3*x^4 + 7*x^6 + 9*x^7.  

    */
    double coef = 0;
    int i;
     
    for(i = 0; i < W[i1][j1].length; i++) {
        W[i1][j1].c[i] += coef;
        coef = W[i1][j1].c[i];
    }
}

double numbersmall(int c, int b, double ob, celW **W1, celW **W2, double tol) {


    /*

    numbersmall is the core of the split-up algorithm. 
    
    It efficiently combines two polynomials which have used 
    complementary sets of rank scores and computes their contribution 
    to the tail-probability. 

    */

    double tot = 0;
    long le;
    int test = 1;
    int be = b/2;
    int bp = (b+1)/2;
    int tempel = 0;
    int i, j, h;
    double th;
    
    for (h = 0; h <= c; h++) {

        tempel = 0;
        le = W2[c-h][bp].length;
        
        for (i = 0; i < W1[h][be].length; i++) {
            test = 1;
            for (j = tempel; j < le && test == 1; j++) {
                th = W1[h][be].x[i] + W2[c-h][bp].x[le-j-1];
                if ((th < ob) | (th - ob < tol)) {
                    tot += W1[h][be].c[i] * W2[c - h][bp].c[le - j -1];
                    tempel = j;
                    test = 0;
                } 
            }
        }
    }
    return(tot);
}
       
SEXP R_split_up_2sample(SEXP scores, SEXP m, SEXP obs, SEXP tol) {
	
    /*

    R interface to the split-up algorithm. 

    `scores' is a REAL vector giving the scores of the total sample 
    and `m' is a scalar integer with the sample size of one group.
    `obs' is the scalar observed test statistic, namely the
    sum of the `m' scores measured in one group.

    */

    int b, c, u;
    double tot, bino, prob;
    double ob;  
    SEXP ans;

    celW **W1;
    celW **W2;
    double *rs;

    b = LENGTH(scores);
    rs = REAL(scores);
    c = INTEGER(m)[0];
    /* d = b - INTEGER(m)[0]; not used */
    ob = REAL(obs)[0];

    /* total number of possible permutations */
    bino = binomi(b, c);

    /* allocate and initialise memory */
    W1 = reserveW(c, (b+1)/2);
    initW(c, (b+1)/2, W1);
    W2 = reserveW(c, (b+1)/2);
    initW(c, (b+1)/2, W2);
    
    makeW(W1, c, b/2, 0, rs, REAL(tol)[0]);  
    makeW(W2, c, (b+1)/2, b/2, rs, REAL(tol)[0]);

    for (u = 0; u <= c; u++) cumulcoef(W2, u, (b+1)/2);

    /* number of permutations <= ob */
    tot = numbersmall(c, b, ob, W1, W2, REAL(tol)[0]);
    
    /* probability */
    prob = tot/bino; 
    
    /* free memory: this will _not_ take place 
       in case of an error */
    FreeW(c, W1);
    FreeW(c, W2);

    /* return to R */
    PROTECT(ans = allocVector(REALSXP, 1));
    REAL(ans)[0] = prob;
    UNPROTECT(1);
    return(ans);
}
