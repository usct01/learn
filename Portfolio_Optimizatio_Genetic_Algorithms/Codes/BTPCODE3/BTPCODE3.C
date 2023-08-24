/********************************************************************\
***                                                                ***
***                      GENETIC ALGORITHM                         ***
***                                                                ***
***               Developed By : Prof. Kalyanmoy Deb               ***
***               with assistance from his Students                ***
***                                                                ***
*** Last Edited : 15.11.2001                                       ***
***................................................................***
This is a GA implementation using binary and real coded
variables. Mixed variables can be used. Constraints can also be
handled. All constraints must be greater-than-equal-to type (g >= 0)
and normalized (see the sample problem in prob1 in objective()).

There are three sample input file (inp-r for real-coded variables only,
inp-b for binary-coded variables only, and inp-rb for a mixed real and binary
variables) which can be used to run this code. The template file for
each input data file is also included (input-real, input-binary, and
input-real+binary).

Code your objective function and constraints at the end of the code
(in objective())

Variable boundaries for real-coded variables can be fixed or flexible.

Following selection opeartor is coded:
Tournament selection: Set MINM=1 for minimization and -1 for maximization
in objective().

For binary strings, single-point crossover and for real parameters
 simulated binary crossover (SBX) are used.

Mutation: bit-wise for Binary coded GAs and polynomial mutation (with eta) for
          Real coded GAs

Constraints are handled using Deb's paramater-less
approach (see CMAME, 2000 paper)

Niching allows restricted tournament selection. Recommended for
complex and disconnected feasible regions. (Niching parameter of 0.1 is
recommended.)

The execution creates a file `result.out' which contains the input
data and best solution obtained by the GA. The feasiblilty of the best
solution and constraint values are also marked.
The report.out contains population record of each generation.
The file 'plot.out' contains a gnuplot-compatibale data file for
plotting best, avg, and worst population fitness versus generation number.

Send Comments to De. K. Deb (deb@iitk.ac.in)
**************************************************************************/

#include<stdio.h>
#include<math.h>
#define BITS_PER_BYTE 8
#define UINTSIZE (BITS_PER_BYTE*sizeof(unsigned))
#define INFINITY 1e7
#define EPSILON  1e-6
#define PI 3.1415927
#define MAXVECSIZE 30
#define MAXPOPSIZE 500
#define MAXCONSTR  10
#define TRUE 1
#define FALSE 0
#define square(x)  ((x)*(x))
/***** Current Objective Function ******/
#define yours  /* define your problem at the end in objfunc() */

/*=================
TYPE DEFINTIONS :
=================*/
struct indiv
{ double xreal[MAXVECSIZE];   /* real-coded variables */
  double  xbin[MAXVECSIZE];          /* binary-coded variables */
  double obj,penalty; /* objective fn. etc. */
  double cons[MAXCONSTR];
  unsigned *chrom;           /* chrosome string      */
  int parent1;
  int parent2;             /* s.no. of parents     */
  int cross_var;         /* cross over variable  */
};
typedef struct indiv INDIVIDUAL ;
typedef INDIVIDUAL *POPULATION ;        /* array of individuals */

/*====================
FUNCTION PROTOTYPES :
====================*/
void     objective();
double   randomperc();
double   get_beta();
double   get_delta();

/*==================
GLOBAL VARIABLES  :
==================*/
int     pop_size,               /* Population Size                      */
  gen_no,                 /* Current generation number            */
  max_gen,                /* Maximum no. of generations           */
  no_xover,               /* No. of cross overs done              */
  no_mutation,binmut,            /* No. of mutations done                */
  best_ever_gen,          /* Generation no. of best ever indiv.   */
  nvar_bin,                /* Number of total design variables     */
  nvar_real,
  lchrom,                 /* Length of chromosome                 */
  chromsize,              /* Number of bytes needed to store
			     lchrom strings          */
  maxrun,                 /* Maxm no. of GA runs for each set of
			     parameter values          */
  run,                    /* Actual run no.                       */
  SHARING,                /* Flag for Sharing ( True / False)     */
  REPORT,                 /* Flag for Full reports (True/False)   */
  RIGID,                  /* Flag for rigid boundaries (T/F)      */
  tourneylist[MAXPOPSIZE],/* List of indices of individuals for
			     tournament selection routine    */
  tourneypos,             /* Current position of tournament       */
  tourneysize,            /* Tournament size ( = 2 for binary )   */
  MINM,
  nc,                     /* Number of constraints */
  critical_size;          /* subpopulation size used in TS (0.25*N)  */

int chr_len[MAXVECSIZE];   /* chrom length for each variable */

double   seed,                   /* Random seed number                   */
  basic_seed,             /* Basic seed number                    */
  n_distribution_c, n_distribution_m,
  p_xover,                /* Cross over probability               */
  p_mutation_bin,p_mutation_real, /* Mutation probability                 */
  sum_obj,                /* Sum of objective fn. values          */
  avg_obj,                /* Average of objective fn. values      */
  max_obj,                /* Maximum objective fn. value          */
  min_obj,                /* Minimum objective fn. value          */
  minx_bin[MAXVECSIZE],       /* Minimum and maximum values of design */
  maxx_bin[MAXVECSIZE],       /*        variables in a population     */
  minx_real[MAXVECSIZE],       /* Minimum and maximum values of design */
  maxx_real[MAXVECSIZE],       /*        variables in a population     */
  xbin_lower[MAXVECSIZE],    /* Lower and Upper bounds on each       */
  xbin_upper[MAXVECSIZE],    /*        design variable               */
  xreal_lower[MAXVECSIZE],    /* Lower and Upper bounds on each       */
  xreal_upper[MAXVECSIZE],    /*        design variable               */
  sigma_share;            /* Sharing distance                     */

double V[8][8] = {{0.000653334,0.000151849,0.000280901,0.000219613,0.000335117,0.000289801,0.000214726,0.000267712},
{0.000151849,0.000426849,0.000148108,0.000139426,0.000126458,0.000116267,0.000113461,0.000117979},
{0.000280901,0.000148108,0.000783555,0.000267724,0.000249795,0.000291752,0.000221797,0.000216618},
{0.000219613,0.000139426,0.000267724,0.000870581,0.000252703,0.000257116,0.000188415,0.000234382},
{0.000335117,0.000126458,0.000249795,0.000252703,0.001665417,0.000243512,0.000259191,0.000575903},
{0.000289801,0.000116267,0.000291752,0.000257116,0.000243512,0.000526708,0.000228075,0.000204555},
{0.000214726,0.000113461,0.000221797,0.000188415,0.000259191,0.000228075,0.000480547,0.000254181},
{0.000267712,0.000117979,0.000216618,0.000234382,0.000575903,0.000204555,0.000254181,0.001031198}}; /*Variance-Covariance matrix */
double R[8] = {0.000741601,0.001586918,0.00154022,0.000909421,-0.000431511,0.001199846,0.000443096,0.000400247}; /*Returns matrix */
double sigma[8] = {0.000653334,0.000426849,0.000783555,0.000870581,0.001665417,0.000526708,0.000480547,0.001031198};

POPULATION oldpop, newpop;      /* Old and New populations              */
INDIVIDUAL best_ever,current_best;           /* Best fit individual till current gen.*/

/*====================================================================
SUBROUTINE FOR INPUTTING GLOBAL PARAMETERS :
====================================================================*/
input_parameters()
{
   int k;
   char ans;

   printf("       ");
   puts("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::");
   printf("       ");
   puts("::::::::       REAL-CODED GENETIC ALGORITHM        :::::::");
   printf("       ");
   puts("::::::::       ============================        :::::::");
   printf("       ");
   puts("::::::::  Kalyanmoy Deb and his students at KanGAL :::::::");
   printf("       ");
   puts("::::::::            All rights reserved.           :::::::");
   printf("       ");
   puts("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::");

   printf("\nHow many generations ? ------------- : ");
   scanf("%d",&max_gen);
   printf("\nPopulation Size ? ------------------ : ");
   scanf("%d", &pop_size );
   if (pop_size > MAXPOPSIZE)
   {
        printf("\n Increase the value of MAXPOPSIZE in program");
        printf("  and re-run the program");
        exit(-1);
   }
   printf("\nNumber of binary-coded variables (Maximum %d) ---- : ",MAXVECSIZE);
   scanf("%d",&nvar_bin);
   printf("\nNumber of real-coded variables (Maximum %d) ---- : ",MAXVECSIZE);
   scanf("%d",&nvar_real);
   if (nvar_bin > 0)
     for (k=0, lchrom = 0; k < nvar_bin; k++)
       {
	 printf("\nString length and Lower and Upper bounds of xbin[%d] ----- : ",k+1);
	 scanf("%d %lf %lf",&chr_len[k],&xbin_lower[k],&xbin_upper[k]);
	 lchrom += chr_len[k];
       }
   if (nvar_real > 0)
     for (k=0; k < nvar_real; k++)
       {
	 printf("\nLower and Upper bounds of xreal[%d] ----- : ",k+1);
	 scanf("%lf %lf",&xreal_lower[k],&xreal_upper[k]);
       }
   if (nvar_real > 0) {
     printf("\n Are the real-parameter bounds rigid ? (y/n) ");
     do { ans = getchar(); } while (ans!= 'y' && ans !='n');
     if (ans == 'y')      RIGID = TRUE;
     else  RIGID = FALSE;
   }
   printf("\nparameter-space niching to be done ? (y/n) --------- : ");
   do { ans = getchar(); } while (ans!= 'y' && ans !='n');
   if (ans == 'y')
   { SHARING = TRUE;
     printf("\nNiching parameter value  ? --------------- : ");
     scanf("%lf",&sigma_share);
   }
   else  SHARING = FALSE;
   printf ("\n Reports to be printed ? (y/n) ");
   do { ans = getchar(); } while (ans!= 'y' && ans !='n');
   if (ans == 'y') REPORT = TRUE;
   else            REPORT = FALSE;
   printf("\n How many runs ? ");
   scanf("%d",&maxrun);
   tourneysize=2;
   printf("\nCross Over Probability? (0 to 1): ");
   scanf("%lf",&p_xover);
   if (nvar_bin > 0)
     {
       printf("\nMutation Probability for binary strings? (0 to 1) : ");
       scanf("%lf",&p_mutation_bin);
     }
   if (nvar_real > 0)
     {
       printf("\nMutation Probability for real variables? (0 to 1) : ");
       scanf("%lf",&p_mutation_real);
     }
   if (nvar_real > 0)
     {
       printf("\n Give distr. index for SBX and mutation? ");
       scanf("%lf %lf",&n_distribution_c,&n_distribution_m);
     }
   printf("\n Give random seed (0 to 1.0) ");
   scanf("%lf",&basic_seed);

   critical_size = pop_size/4;
   input_app_parameters();
   
}

/*====================================================================
Initialses zero'th generation and global parameters
====================================================================*/
initialize()
{
   double u;
   int k,k1,i,j,j1,stop;
   double temp[MAXVECSIZE],coef;
   unsigned mask=1,nbytes;

   randomize();
   app_initialize();
   oldpop = (INDIVIDUAL *)malloc(pop_size*sizeof(INDIVIDUAL));
   newpop = (INDIVIDUAL *)malloc(pop_size*sizeof(INDIVIDUAL));

   if (oldpop == NULL) nomemory("oldpop in initialize()");
   if (newpop == NULL) nomemory("newpop in initialize()");

   chromsize = (lchrom/UINTSIZE);
   if(lchrom%UINTSIZE) chromsize++;
   nbytes = chromsize*sizeof(unsigned);
   if (nvar_bin > 0)
     for(j = 0; j < pop_size; j++)
       {
	 if((oldpop[j].chrom = (unsigned *) malloc(nbytes)) == NULL)
	   nomemory("oldpop chromosomes");

	 if((newpop[j].chrom = (unsigned *) malloc(nbytes)) == NULL)
	   nomemory("newpop chromosomes");
       }
   if (nvar_bin > 0)
     {
       if((best_ever.chrom = (unsigned *) malloc(nbytes)) == NULL)
	 nomemory("best_ever chromosomes");
       if((current_best.chrom = (unsigned *) malloc(nbytes)) == NULL)
	 nomemory("current_best chromosomes");
     }

   for (k=0; k < pop_size; k++)
   {
     oldpop[k].obj = 0.0;
     oldpop[k].parent1 = oldpop[k].parent2 = 0;
     oldpop[k].penalty = 0.0; oldpop[k].cross_var = 0;
     for (j=0; j < nc; j++)
       oldpop[k].cons[j] = 0.0;

     for (j=0; j < nvar_real; j++)
       {
         u = randomperc();
         oldpop[k].xreal[j] = xreal_lower[j] * (1-u) + xreal_upper[j] * u;
      }

      for(k1 = 0; k1 < chromsize; k1++)
      {
            oldpop[k].chrom[k1] = 0;
            if(k1 == (chromsize-1))
                stop = lchrom - (k1*UINTSIZE);
            else
                stop = UINTSIZE;
            /* A fair coin toss */

           for(j1 = 1; j1 <= stop; j1++)
            {
               if(flip(0.5))
                  oldpop[k].chrom[k1] = oldpop[k].chrom[k1]|mask;
               if (j1 != stop) oldpop[k].chrom[k1] = oldpop[k].chrom[k1]<<1;
            }
      }
   }
   no_xover = no_mutation = binmut = 0;

   copy_individual(&oldpop[0],&best_ever);
   decode_string(&best_ever);
   objective(&best_ever);
}

/*====================================================================
Decodes the string of the individual (if any) and puts the values in
the array of doubles.
====================================================================*/
decode_string(ptr_indiv)
INDIVIDUAL *ptr_indiv;
{
   double *temp,coef;
   int j;

   if (ptr_indiv == NULL) error_ptr_null("ptr_indiv in decode_string");
   if (nvar_bin > 0)
   {
     temp = (double *) malloc(nvar_bin * sizeof(double));
     for(j=0; j < nvar_bin; j++)
       temp[j] = 0.0;
     decodevalue(ptr_indiv->chrom,temp);
     for(j=0; j < nvar_bin; j++)
     {
       coef = pow(2.0,(double)(chr_len[j])) - 1.0;
       temp[j] = temp[j]/coef;
       ptr_indiv->xbin[j] = temp[j]*xbin_upper[j] + (1.0 - temp[j])*xbin_lower[j];
     }
     free(temp);
   }
}
/*====================================================================
Prints an error message and terminates the program
====================================================================*/
nomemory(string)
char *string;
{
   printf("\nmalloc: out of memory making %s!!\n",string);
   printf("\n Program is halting .....");
   exit(-1);
}
/*==============================================================
Gives error message of null pointer  and terminates the program.
==============================================================*/
error_ptr_null(string)
char *string;
{
   printf("\n Error !! Pointer %s found Null !",string);
   printf("\n Program is halting .....");
   exit(-1);
}
/*====================================================================
Copys contents of one individual into another.
====================================================================*/
copy_individual(indiv1,indiv2)
INDIVIDUAL *indiv1, *indiv2;
{
   int k;

   if (indiv1==NULL) error_ptr_null("indiv1 in copy_individual");
   if (indiv2==NULL) error_ptr_null("indiv2 in copy_individual");

   for (k=0; k < nvar_bin; k++)
      indiv2->xbin[k] = indiv1->xbin[k];
   for (k=0; k < nvar_real; k++)
      indiv2->xreal[k] = indiv1->xreal[k];
   for (k=0; k < nc; k++)
     indiv2->cons[k] = indiv1->cons[k];
   indiv2->obj = indiv1->obj;
   indiv2->penalty = indiv1->penalty;
   indiv2->parent1 = indiv1->parent1;
   indiv2->parent2 = indiv1->parent2;
   indiv2->cross_var = indiv1->cross_var;
   for (k=0; k < chromsize; k++)
     indiv2->chrom[k] = indiv1->chrom[k];
}

/*====================================================================
Calculates statistics of current generation :
====================================================================*/
statistics(gen)
int gen;
{
  // INDIVIDUAL current_best;
   int k,j, change_flag;
   double f;
   double pow(),bitpow,coef,temp[MAXVECSIZE];

   for (k=0; k < pop_size; k++)
   {
     decode_string(&oldpop[k]);
     objective(&(oldpop[k]));
   }
   copy_individual(&oldpop[0], &current_best);
   sum_obj = avg_obj = oldpop[0].obj;
   max_obj = min_obj = oldpop[0].obj;
   for (k=0;  k < nvar_bin; k++)
     maxx_bin[k] = minx_bin[k] = oldpop[0].xbin[k];
   for (k=0;  k < nvar_real; k++)
     maxx_real[k] = minx_real[k] = oldpop[0].xreal[k];
   for(k=1; k < pop_size; k++)
     {
       if (current_best.penalty > oldpop[k].penalty)
	 copy_individual(&oldpop[k], &current_best);
       else if ((current_best.penalty <= 0.0) && (oldpop[k].penalty <= 0.0))
	 if (MINM * current_best.obj  >  MINM * oldpop[k].obj)
	   copy_individual(&oldpop[k], &current_best);
     if(MINM * max_obj < MINM * oldpop[k].obj)
                   max_obj = oldpop[k].obj;
     if(MINM * min_obj > MINM * oldpop[k].obj)
                   min_obj = oldpop[k].obj;
     sum_obj += oldpop[k].obj;
     for (j=0; j < nvar_bin; j++)
     {
        if (MINM * maxx_bin[j] < MINM * oldpop[k].xbin[j]) maxx_bin[j] = oldpop[k].xbin[j];
        if (MINM * minx_bin[j] > MINM * oldpop[k].xbin[j]) minx_bin[j] = oldpop[k].xbin[j];
     }
     for (j=0; j < nvar_real; j++)
     {
        if (MINM * maxx_real[j] < MINM * oldpop[k].xreal[j]) maxx_real[j] = oldpop[k].xreal[j];
        if (MINM * minx_real[j] > MINM * oldpop[k].xreal[j]) minx_real[j] = oldpop[k].xreal[j];
     }
   }
   avg_obj = sum_obj/pop_size;
   change_flag = 0;
   if (best_ever.penalty > current_best.penalty) change_flag = 1;
   else if ((best_ever.penalty <= 0.0) && (current_best.penalty <= 0.0))
     if (MINM * best_ever.obj > MINM * current_best.obj)
       change_flag = 1;
   if (change_flag == 1)
     {
       copy_individual(&current_best, &best_ever);
       best_ever_gen = gen;
     }
   app_statistics();
}

/*====================================================================
Decodes the value of a group of binary strings and puts the decoded
values into an array 'value'.
====================================================================*/
decodevalue(chrom,value)
unsigned *chrom;
double value[];
{
    int k,j,stop,tp,bitpos,mask=1,position,bits_per_var,count;
    double pow(), bitpow;

    if (nvar_bin <= 0) return;
    if (chrom == NULL) error_ptr_null("chrom in decodevalue");

    position = 0; count = 0;
    for(k = 0; k < chromsize; k++)
    {
        if(k == (chromsize-1))
            stop = lchrom-(k*UINTSIZE);
        else
            stop = UINTSIZE;
        /* loop thru bits in current byte */
        tp = chrom[k];
        for(j = 0; j < stop; j++) {
	  bits_per_var = chr_len[position];
	  bitpos = j + UINTSIZE*k;
            /* test for current bit 0 or 1 */
            if((tp&mask) == 1) {
	      // position = bitpos / bits_per_var;
              //  bitpos -= position * bits_per_var;
                bitpow = pow(2.0,(double)(count));
                value[position] += bitpow;
            }
            tp = tp>>1; count++;
	    if (count >= chr_len[position])
	      {
		position += 1;
		count = 0;
	      }

	}
    }
}

/*====================================================================
GENERATION OF NEW POPULATION through SELECTION, XOVER & MUTATION :
====================================================================*/
generate_new_pop()
{
   int k,mate1,mate2;

   app_computation();

   preselect_tour();

   for (k=0; k < pop_size; k += 2)
   {
     // selection
     if (SHARING)
       {
	 mate1 = tour_select_constr();
	 mate2 = tour_select_constr();
       }
     else
       {
	 mate1 = tour_select();
	 mate2 = tour_select();
       }
     // crossover
     cross_over(mate1,mate2,k,k+1);
     // mutation
     mutation(&newpop[k]);
     mutation(&newpop[k+1]);
     newpop[k].parent1 = newpop[k+1].parent1 = mate1+1;
     newpop[k].parent2 = newpop[k+1].parent2 = mate2+1;
   }
}

/*====================================================================
Binary cross over routine.
====================================================================*/
binary_xover (parent1, parent2, child1, child2, x_site)
unsigned *parent1, *parent2, *child1, *child2;
int *x_site;
/* Cross 2 parent strings, place in 2 child strings */
{
    int j, jcross, k;
    unsigned mask, temp;

    if (nvar_bin <= 0) return;
    if (parent1 == NULL) error_ptr_null("parent1 in binary_xover");
    if (parent2 == NULL) error_ptr_null("parent2 in binary_xover");
    if (child1== NULL) error_ptr_null("child1 in binary_xover");
    if (child2== NULL) error_ptr_null("child2 in binary_xover");

    jcross = rnd(1 ,(lchrom - 1));/* Cross between 1 and l-1 */
    for(k = 1; k <= chromsize; k++)
    {
            if(jcross >= (k*UINTSIZE))
            {
                child1[k-1] = parent1[k-1];
                child2[k-1] = parent2[k-1];
            }
            else if((jcross < (k*UINTSIZE)) && (jcross > ((k-1)*UINTSIZE)))
            {
                mask = 1;
                for(j = 1; j <= (jcross-1-((k-1)*UINTSIZE)); j++)
                {
                    temp = 1;
                    mask = mask<<1;
                    mask = mask|temp;
                }
                child1[k-1] = (parent1[k-1]&mask)|(parent2[k-1]&(~mask));
                child2[k-1] = (parent1[k-1]&(~mask))|(parent2[k-1]&mask);
            }
            else
            {
                child1[k-1] = parent2[k-1];
                child2[k-1] = parent1[k-1];
            }
    }
    *x_site = jcross;
}
/*====================================================================
Creates two children from parents p1 and p2, stores them in addresses
pointed by c1 and c2.  low and high are the limits for x values and
rand_var is the random variable used to create children points.
====================================================================*/
create_children(p1,p2,c1,c2,low,high,rand_var)
double p1,p2,*c1,*c2,low,high,*rand_var;
{
    double difference,x_mean,beta,v2,v1;
    double u,distance,umax,temp,alpha;
    int flag;

    if (c1 == NULL) error_ptr_null("c1 in create_children");
    if (c2 == NULL) error_ptr_null("c2 in create_children");
    if (rand_var == NULL) error_ptr_null("rand_var in create_children");
    flag = 0;
    if ( p1 > p2) { temp = p1; p1 = p2; p2 = temp; flag = 1; }
    x_mean = (p1 + p2) * 0.5;
    difference = p2 - p1;
    if ( (p1-low) < (high-p2) ) distance = p1-low;
    else                        distance = high-p2;
    if (distance < 0.0) distance = 0.0;
    if (RIGID && (difference > EPSILON))
    {
      alpha = 1.0 + (2.0*distance/difference);
      umax = 1.0 - (0.5 / pow((double)alpha,(double)(n_distribution_c+1.0)));
      *rand_var = umax * randomperc();
    }
    else *rand_var = randomperc();
    beta = get_beta(*rand_var);
    if (fabs(difference*beta) > INFINITY) beta = INFINITY/difference;
    v2 = x_mean + beta * 0.5 * difference;
    v1 = x_mean - beta * 0.5 * difference;

    if (v2 < low) v2 = low;
    if (v2 > high) v2 = high;
    if (v1 < low) v2 = low;
    if (v1 > high) v2 = high;
    *c2 = v2; *c1 = v1;
    //  if (flag == 1) { temp = *c1; *c1 = *c2; *c2 = temp; }

}

/*====================================================================
CROSS - OVER  USING strategy of uniform 50% variables
  For one variable problem, it is crossed over as usual.
  For multivariables, each variable is crossed over with a probability
  of 50 % , each time generating a new random beta.
====================================================================*/
cross_over(first,second,childno1,childno2)
int first,second,childno1,childno2;
{
    double difference,x_mean,beta;
    double u = 0.0;
    int site,k,x_s;

    x_s = 0;
    if (flip(p_xover))   /* Cross over has to be done */
    {
     no_xover++;
     if (nvar_bin > 0)
       {
	 binary_xover(oldpop[first].chrom,oldpop[second].chrom,
		      newpop[childno1].chrom,newpop[childno2].chrom,&x_s);
	 newpop[childno1].cross_var = newpop[childno2].cross_var = x_s;
       }
     if (nvar_real > 0)
       {
	 for (site = 0; site < nvar_real; site++)
	   {
	     if(flip(0.5) || (nvar_real == 1))
	       {
		 create_children(oldpop[first].xreal[site],oldpop[second].xreal[site],
		    &(newpop[childno1].xreal[site]),&(newpop[childno2].xreal[site]),
		     xreal_lower[site],xreal_upper[site],&u);
	       }
	     else
	       {
		 newpop[childno1].xreal[site] = oldpop[first].xreal[site];
		 newpop[childno2].xreal[site] = oldpop[second].xreal[site];
	       }
	   }               /* for loop */
	 if (nvar_bin == 0)
	   newpop[childno1].cross_var = newpop[childno2].cross_var = 0;
       }                /* if REALGA      */
    }                 /* Cross over done */

    else              /* Passing x-values straight */
    {
      for (k=0; k < chromsize; k++)
	{
	  newpop[childno1].chrom[k] = oldpop[first].chrom[k];
	  newpop[childno2].chrom[k] = oldpop[second].chrom[k];
	}
      for (site=0; site < nvar_real; site++)
	{
	  newpop[childno1].xreal[site] = oldpop[first].xreal[site];
	  newpop[childno2].xreal[site] = oldpop[second].xreal[site];
	}
      for (site=0; site < nvar_bin; site++)
	{
	  newpop[childno1].xbin[site] = oldpop[first].xbin[site];
	  newpop[childno2].xbin[site] = oldpop[second].xbin[site];
	}
      newpop[childno1].cross_var = newpop[childno2].cross_var = 0;
    }
}

/*===================================================================
Calculates beta value for given random number u (from 0 to 1)
If input random numbers (u) are uniformly distributed for a set of
inputs, this results in uniform distribution of beta values in case
of BLX , and Binary Probability distribution simulation in case of
SBX.
====================================================================*/
double get_beta(u)
double u;
{
   double beta;

   if (1.0-u < EPSILON ) u = 1.0 - EPSILON;
   if ( u < 0.0) u = 0.0;
   if (u < 0.5) beta = pow(2.0*u,(1.0/(n_distribution_c+1.0)));
   else beta = pow( (0.5/(1.0-u)),(1.0/(n_distribution_c+1.0)));
   return beta;
}
/*==================================================================
For given u value such that   -1 <= u <= 1, this routine returns a
value of delta from -1 to 1. Exact value of delta depends on specified
n_distribution. This is called by mutation().
====================================================================*/
double get_delta(u, delta_l, delta_u)
     double u, delta_l, delta_u;
{
  double delta, aa;

  if (u >= 1.0-1.0e-9)      delta = delta_u;
  else if (u <= 0.0+1.0e-9) delta = delta_l;
  else
    {
      if (u <= 0.5)
	{
	  aa = 2.0*u + (1.0-2.0*u)*pow((1+delta_l),(n_distribution_m + 1.0));
	  delta = pow(aa, (1.0 / (n_distribution_m + 1.0))) - 1.0;
	}
      else
	{
	  aa = 2.0*(1-u) + 2.0*(u-0.5)*pow((1-delta_u),(n_distribution_m + 1.0));
	  delta = 1.0 - pow(aa, (1.0 / (n_distribution_m + 1.0)));
	}
    }
  if (delta < -1.0 || delta > 1.0)
    {
      printf("Error in mutation!! delta = %lf\n",delta);
      exit(-1);
    }
  return (delta);
}

/*==================================================================
Binary mutation routine ( borrowed from sga.c )
====================================================================*/
binmutation(child)
unsigned *child;
/* Mutate an allele w/ pmutation, count # of mutations */
{
    int j, k, stop;
    unsigned mask, temp = 1;

    if (nvar_bin <= 0) return;
    if (child== NULL) error_ptr_null(" child in binmutation");
    for(k = 0; k < chromsize; k++)
    {
        mask = 0;
        if(k == (chromsize-1))
	  stop = lchrom - (k*UINTSIZE);
        else
            stop = UINTSIZE;
        for(j = 0; j < stop; j++)
        {
            if(flip(p_mutation_bin))
            {
                mask = mask|(temp<<j);
		binmut++;
            }
        }
        child[k] = child[k]^mask;
    }
}

/*===================================================================
Mutation Using polynomial probability distribution. Picks up a random
site and generates a random number u between -1 to 1, ( or between
minu to maxu in case of rigid boudaries) and calls the routine
get_delta() to calculate the actual shift of the value.
====================================================================*/
mutation(indiv)
INDIVIDUAL  *indiv;
{
   double distance1,x,delta_l,delta_u,delta,u;
   int k, site;

   if (indiv == NULL) error_ptr_null("indiv in mutation");

   if (nvar_real > 0)
     for (site = 0; site < nvar_real; site++)
       {
	 if(flip(p_mutation_real))
	   {
	     no_mutation++;
	     if(RIGID)
	       {
		 x = indiv->xreal[site];
		 distance1 = xreal_lower[site] - x;
		 delta_l = distance1/(xreal_upper[site] - xreal_lower[site]);
		 if (delta_l < -1.0)  delta_l = -1.0;

		 distance1 = xreal_upper[site] - x;
		 delta_u = distance1/(xreal_upper[site] - xreal_lower[site]);
		 if (delta_u > 1.0)   delta_u = 1.0;

		 if (-1.0*delta_l < delta_u) delta_u = -1.0 * delta_l;
		 else delta_l = -1.0 * delta_u;
	       }
	     else
	       {
		 delta_l = -1.0;
		 delta_u =  1.0;
	       }
	     u = randomperc();
	     /* calculation of actual delta value */
	     delta = get_delta(u, delta_l, delta_u)
	       * (xreal_upper[site] - xreal_lower[site]);
	     indiv->xreal[site] += delta;

	     if (indiv->xreal[site] < xreal_lower[site])
	       indiv->xreal[site] = xreal_lower[site];
	     if (indiv->xreal[site] > xreal_upper[site])
	       indiv->xreal[site] = xreal_upper[site];

	   }    /* if flip() */
       }
   if (nvar_bin > 0)
     binmutation(indiv->chrom);
}

/*====================================================================
  Reporting the user-specified parameters :
  fp is the file pointer to output file.
====================================================================*/
initreport(fp)
FILE *fp;
{
   int k;

   if (fp == NULL) error_ptr_null(" File fp in initreport");
   fprintf(fp,"\n\n=============================================");
   fprintf(fp,"\n             INITIAL REPORT                  ");
   fprintf(fp,"\n=============================================");

   fprintf(fp,"\n Variable Boundaries : ");
   if (nvar_real > 0)
     if (RIGID) fprintf(fp," Rigid");
     else       fprintf(fp," Flexible");
   fprintf(fp,"\n Population size            : %d",pop_size);
   fprintf(fp,"\n Total no. of generations   : %d",max_gen);
   fprintf(fp,"\n Cross over probability     : %6.4f",p_xover);
   if (nvar_bin > 0)
     fprintf(fp,"\n Mutation probability (binary): %6.4f",p_mutation_bin);
   if (nvar_real > 0)
     fprintf(fp,"\n Mutation probability (real): %6.4f",p_mutation_real);
   if (SHARING)
   {
        fprintf(fp,"\n Niching to be done :");
        fprintf(fp,"\n Niching parameter value: %6.4f",sigma_share);
   }
   if (nvar_bin > 0)
     fprintf(fp,"\n Total String length              : %d",lchrom);
   if (nvar_bin > 0)
     fprintf(fp,"\n Number of binary-coded variables: %d",nvar_bin);
   if (nvar_real > 0)
     fprintf(fp,"\n Number of real-coded variables  : %d",nvar_real);
   fprintf(fp,"\n Total Runs to be performed : %d",maxrun);
   if (nvar_real > 0)
     {
       fprintf(fp,"\n Exponent (n for SBX)       : %7.2f",n_distribution_c);
       fprintf(fp,"\n Exponent (n for Mutation)  : %7.2f",n_distribution_m);
     }
   fprintf(fp,"\n Lower and Upper bounds     :");
   for (k=0; k < nvar_bin; k++)
     fprintf(fp,"\n   %8.4f   <=   x_bin[%d]   <= %8.4f, string length = %d",xbin_lower[k],k+1,xbin_upper[k],chr_len[k]);
   for (k=0; k < nvar_real; k++)
     fprintf(fp,"\n   %8.4f   <=   x_real[%d]   <= %8.4f",xreal_lower[k],k+1,xreal_upper[k]);

   fprintf(fp,"\n=================================================\n");

   app_initreport();
}

/*====================================================================
Writes a given string of 0's and 1's
puts a `-` between each substring (one substring for one variable)
Leftmost bit is most significant bit.
====================================================================*/
writechrom(chrom,fp)
unsigned *chrom;
FILE *fp;
{
    int j, k, stop,bits_per_var,count=0, bitcount, position;
    unsigned mask = 1, tmp;

   if (fp == NULL) error_ptr_null(" File fp in initreport");

   if (nvar_bin <= 0) return;
   if (chrom == NULL) error_ptr_null("chrom in writechrom");

   position = 0;
   bitcount = 0;
    for(k = 0; k < chromsize; k++)
    {
        tmp = chrom[k];
        if(k == (chromsize-1))
            stop = lchrom - (k*UINTSIZE);
        else
            stop = UINTSIZE;

        for(j = 0; j < stop; j++)
        {
	  bits_per_var = chr_len[position];
	  if(tmp&mask)
	    fprintf(fp,"1");
	  else
	    fprintf(fp,"0");
	  count++; bitcount++;
	  if (( (count % bits_per_var) == 0) && (count < lchrom))
	    fprintf(fp,"-");
	  tmp = tmp>>1;
	  if (bitcount >= chr_len[position])
	    {
	      bitcount = 0;
	      position += 1;
	    }
        }
    }
}

void writeindv(ind)
     INDIVIDUAL ind;
{
  int i;

  for (i=0; i<nvar_bin; i++)
    printf(" %lf",ind.xbin[i]);
  for (i=0; i<nvar_real; i++)
    printf(" %lf",ind.xreal[i]);
  printf("Obj=%lf penalty = %lf parent1 = %d parent2 = %d",ind.obj,ind.penalty,ind.parent1,ind.parent2);
  printf("cross_site = %lf ",ind.cross_var);
  for (i=0; i<chromsize; i++)
    printf(" %d",ind.chrom[i]);
}

/*====================================================================
Reporting the statistics of current population ( gen. no. 'num'):
  fp is file pointer to output file.
====================================================================*/
report(fr,fp,fl,num)
FILE *fr,*fp,*fl;
int num;
{
  int k,j;
  char string[30];


  if (fr == NULL) error_ptr_null(" file fr in report()");
  if (fp == NULL) error_ptr_null(" file fp in report()");
  if (num == 0)
    fprintf(fl,"\n# Generation Number  Best Fitness  Average Fitness  Worst Fitness");
  fprintf(fl,"\n %d  %lf  %lf  %lf",num,min_obj,avg_obj,max_obj);

  if (REPORT)
    {
      /* ----------------------------------------- */
      /* WRITING IN THE OUTPUT FILE FOR INSPECTION */
      /* ----------------------------------------- */
      fprintf(fr,"\n========== Generation No. : %3d ===================",num);
      fprintf(fr,"\n  No. Binary|Real x   Constr. violation    Fitness    Parents     Cross-site ");
      fprintf(fr,"\n===================================================");
      for (k=0; k < pop_size; k++)
	{
	  fprintf(fr,"\n%3d.  ",k+1);
	  if (nvar_bin > 0)
	    {
	      for (j=0; j < nvar_bin; j++)
		fprintf(fr," %8.5lf",oldpop[k].xbin[j]);
	    }
	  fprintf(fr,"|");

	  if (nvar_real > 0)
	    {
	      for (j=0; j<nvar_real; j++)
		fprintf(fr," %8.5lf",oldpop[k].xreal[j]);
	    }
	  for (j=0; j<nc; j++)
	    fprintf(fr," %8.5lf",oldpop[k].cons[j]);
	  fprintf(fr," %8.5lf",oldpop[k].penalty);
	  fprintf(fr,"  %8.5lf (%3d %3d)", oldpop[k].obj,
		  oldpop[k].parent1,oldpop[k].parent2);
	  if (nvar_bin > 0) fprintf(fr,"  %d",oldpop[k].cross_var);
	  if (nvar_bin > 0)
	    {
	      fprintf(fr,"\n String = ");
	      writechrom(oldpop[k].chrom,fr);
	    }
	}
    }
  if (num==max_gen)
    {
      fprintf(fp,"\n===================================================");
      fprintf(fp,"\nMax = %8.5lf  Min = %8.5lf   Avg = %8.5lf",
	      max_obj,min_obj,avg_obj);
      fprintf(fp,"\nMutations (real)= %d ; Mutations (binary) = %d ; Crossovers = %d",
	      no_mutation,binmut,no_xover);
      fprintf(fl,"\n");

      if (best_ever.penalty <= 0.0)
	{
	  fprintf(fp,"\nBest ever fitness: %lf (from generation : %d)\n",
		  best_ever.obj,best_ever_gen);
	  fprintf(fp,"Variable vector: Binary | Real -> ");
	  for (j=0; j < nvar_bin; j++)
	    fprintf(fp," %lf",best_ever.xbin[j]);
	  fprintf(fp,"|");
	  for (j=0; j < nvar_real; j++)
	    fprintf(fp," %lf",best_ever.xreal[j]);
	  if (nvar_bin > 0)
	    {
	      fprintf(fp,"\nBest_ever String = ");
	      writechrom(best_ever.chrom,fp);
	    }
	  fprintf(fp,"\nConstraint value:");
	  for (j=0; j < nc; j++)
	    fprintf(fp," %lf",best_ever.cons[j]);
	  fprintf(fp,"| Overall penalty: %lf",best_ever.penalty);
	}
      else
	fprintf(fp,"No feasible solution found!\n");
      fprintf(fp,"\n===================================================");
      fprintf(fp,"\n\n");
    }

  app_report();
}
/*====================================================================
Releases the memory for all mallocs
====================================================================*/
free_all()
{
  int i;

  if (nvar_bin > 0)
    {
      for(i = 0; i < pop_size; i++)
	{
	  free(oldpop[i].chrom);
	  free(newpop[i].chrom);
	}
      free(best_ever.chrom);
      free(current_best.chrom);
    }
  free(oldpop);
  free(newpop);

  app_free();
}

/*====================================================================
MAIN PROGRAM ;
====================================================================*/
main()
{
   FILE *fp_out, *fp_rep, *fp_plot; /* File pointer for output file         */
   int runno=0, k;
   POPULATION   temp;           /* A temporary pointer of population    */

/*---------------------------*/
/* Program starts here :     */
/*---------------------------*/
   input_parameters();
   fp_out = fopen("result.out","w+");
   fp_rep = fopen("report.out","w");
   fp_plot= fopen("plot.out","w");

   select_memory();
   initreport(fp_out);
   for (run = 1; run <= maxrun; run++)
     {
       printf("\nRun No. %d :  Wait Please .........",run);
       fprintf(fp_out,"\nRun No. %d ",run);
       seed = basic_seed + (1.0-basic_seed)*(double)(run-1)/(double)maxrun;
       if (seed > 1.0) printf("\n Warning !!! seed number exceeds 1.0");
       gen_no = 0;
       initialize();
       statistics(gen_no);
       report(fp_rep,fp_out,fp_plot,0);
       for(gen_no = 1; gen_no <= max_gen; gen_no++)
	 {
	   generate_new_pop();

	   temp = oldpop;
	   oldpop = newpop;
	   newpop = temp;

	   statistics(gen_no);
	   report(fp_rep,fp_out,fp_plot,gen_no);
	 };
     /* One GA run is over  */
       free_all();
     }                      /* for loop of run  */

   fclose(fp_out); fclose(fp_rep); fclose(fp_plot);
   app_closure();
   printf("\n Results are stored in file 'result.out', 'report.out', and 'plot.out'\n");
}

/**************** End of Main Program ***************************/

/*-------------------------------------------------------  */
/* random.c - contains random number generator and related */
/* utilities,                                              */
/* Source : sga.c  (c) E.Goldberg 1986
/*-------------------------------------------------------  */


/* variables are declared static so that they cannot       */
/* conflict with names of other global variables in other  */
/* files.  See K&R, p 80, for scope of static              */

static double oldrand[55];   /* Array of 55 random numbers */
static int jrand;                 /* current random number */

advance_random()
/* Create next batch of 55 random numbers */
{
    int j1;
    double new_random;

    for(j1 = 0; j1 < 24; j1++)
    {
        new_random = oldrand[j1] - oldrand[j1+31];
        if(new_random < 0.0) new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    }
    for(j1 = 24; j1 < 55; j1++)
    {
        new_random = oldrand [j1] - oldrand [j1-24];
        if(new_random < 0.0) new_random = new_random + 1.0;
        oldrand[j1] = new_random;
    }
}

int flip(prob)
/* Flip a biased coin - true if heads */
double prob;
{
    double randomperc();

    if(randomperc() <= prob)
        return(1);
    else
        return(0);
}

randomize()
/* Get seed number for random and start it up */
{
    int j1;

    for(j1=0; j1<=54; j1++) oldrand[j1] = 0.0;
    jrand=0;

    warmup_random(seed);
}

double randomperc()
/* Fetch a single random number between 0.0 and 1.0 -  */
/* Subtractive Method . See Knuth, D. (1969), v. 2 for */
/* details.Name changed from random() to avoid library */
/* conflicts on some machines                          */
{
    jrand++;
    if(jrand >= 55)
    {
        jrand = 1;
        advance_random();
    }
    return((double) oldrand[jrand]);
}

int rnd(low, high)
/* Pick a random integer between low and high */
int low,high;
{
    int i;
    double randomperc();

    if(low >= high)
        i = low;
    else
    {
        i = (randomperc() * (high - low + 1)) + low;
        if(i > high) i = high;
    }
    return(i);
}

double rndreal(lo ,hi)
/* real random number between specified limits */
double lo, hi;
{
    return((randomperc() * (hi - lo)) + lo);
}

warmup_random(random_seed)
/* Get random off and running */
double random_seed;
{
    int j1, ii;
    double new_random, prev_random;

    oldrand[54] = random_seed;
    new_random = 0.000000001;
    prev_random = random_seed;
    for(j1 = 1 ; j1 <= 54; j1++)
    {
        ii = (21*j1)%54;
        oldrand[ii] = new_random;
        new_random = prev_random-new_random;
        if(new_random<0.0) new_random = new_random + 1.0;
        prev_random = oldrand[ii];
    }

    advance_random();
    advance_random();
    advance_random();

    jrand = 0;
}
/*----------------------------------------------------------*/
/* Files for tournament selection :                         */
/* Source : sga.c (c) E.Goldberg                            */
/*----------------------------------------------------------*/

select_memory()
{
  unsigned nbytes;

  if(tourneysize > pop_size)
    {
      printf("FATAL: Tournament size (%d) > pop_size (%d)\n",
              tourneysize,pop_size);
      exit(-1);
    } ;
}


preselect_tour()
{
    reset1();
    tourneypos = 0;
}


int tour_select()
{
    int pick, winner, i;

    /* If remaining members not enough for a tournament, then reset list */
start_select :
    if((pop_size - tourneypos) < tourneysize)
    {
        reset1();
        tourneypos = 0;
    }

    /* Select tourneysize structures at random and conduct a tournament */
    winner=tourneylist[tourneypos];
/* Added by RBA */
    if( winner < 0 || winner > pop_size-1) {
                                             printf("\n Warning !! ERROR1");
                                             printf(" tourpos = %d",tourneypos);
                                             printf(" winner = %d",winner);
                                             preselect_tour();
                                             goto start_select; }
    for(i=1; i<tourneysize; i++)
    {
        pick=tourneylist[i+tourneypos];
/* Added by RBA */
        if (pick < 0 || pick > pop_size-1) { preselect_tour();
                                             printf("\n Warning !! ERROR2");
                                             goto start_select; }
	// case 1:
	if (oldpop[winner].penalty > oldpop[pick].penalty) winner = pick;
	else if ((oldpop[winner].penalty <= 0.0) && (oldpop[pick].penalty <= 0.0))
	  {
	      if(MINM * oldpop[pick].obj < MINM * oldpop[winner].obj) winner=pick;
	  }
    }

    /* Update tourneypos */
    tourneypos += tourneysize;
    return(winner);
}

double distanc(one,two)
     int one,two;
{
  int k;
  double sum;

  sum = 0.0;
  for (k=0; k<nvar_bin; k++)
    sum += square((oldpop[one].xbin[k]-oldpop[two].xbin[k])/(xbin_upper[k]-xbin_lower[k]));
  for (k=0; k<nvar_real; k++)
    sum += square((oldpop[one].xreal[k]-oldpop[two].xreal[k])/(xreal_upper[k]-xreal_lower[k]));

  return (sqrt(sum/(nvar_bin + nvar_real)));
};

int tour_select_constr()
{
  int pick, winner, i, minus=0, rand_pick, rand_indv, flag, indv;

  /* If remaining members not enough for a tournament, then reset list */
  start_select :
    if((pop_size - tourneypos) < tourneysize)
      {
        reset1();
        tourneypos = 0;
      }

  /* Select tourneysize structures at random and conduct a tournament */
  winner = tourneylist[tourneypos];
  /* Added by RBA */
  if( winner < 0 || winner > pop_size-1)
    {
      printf("\n Warning !! ERROR1");
      printf(" tourpos = %d",tourneypos);
      printf(" winner = %d",winner);
      preselect_tour();
      goto start_select;
    }
  for(i=1; i<tourneysize; i++)
    {
      pick = tourneylist[i+tourneypos];

      if((oldpop[winner].penalty>0.0) && (oldpop[pick].penalty<=0.0))
	winner = pick;
      else if((oldpop[winner].penalty>0.0)&&(oldpop[pick].penalty>0.0))
	{
	  if(oldpop[pick].penalty < oldpop[winner].penalty)
	    winner=pick;
	}
      else if((oldpop[winner].penalty<=0.0)&&(oldpop[pick].penalty<=0.0))
	{
	  if(distanc(winner,pick) < sigma_share)
	    {
	      if (MINM * oldpop[pick].obj < MINM * oldpop[winner].obj)
		winner=pick;
	    }
	  else
	    {
	      minus = -1;
	      for (indv = flag = 0; indv<critical_size && flag==0; indv++)
		{
		  rand_indv = rnd(0,pop_size-1);
		  rand_pick = tourneylist[rand_indv];
		  if(oldpop[rand_pick].penalty <= 0.0)
		    {
		      if(distanc(winner,rand_pick) < sigma_share)
			{
			  flag = 1;
			  if (MINM * oldpop[rand_pick].obj < MINM * oldpop[winner].obj)
			    winner=rand_pick;
			}
		    }
		}
	    }
	}
      if (pick < 0 || pick > pop_size-1)
	{
	  preselect_tour();
	  printf("\n Warning !! ERROR2");
	  goto start_select;
	}
    }
    /* Update tourneypos */
  tourneypos += tourneysize + minus;
  return(winner);
}

reset1()
/* Name changed from reset because of clash with lib. function - RBA */
/* Shuffles the tourneylist at random */
{
    int i, rand1, rand2, temp_site;

    for(i=0; i<pop_size; i++) tourneylist[i] = i;

    for(i=0; i < pop_size; i++)
    {
        rand1= rnd(0,pop_size-1);
        rand2=  rnd(0,pop_size-1);
        temp_site = tourneylist[rand1];
        tourneylist[rand1]=tourneylist[rand2];
        tourneylist[rand2]=temp_site;
    }
}
/******************* APPLICATION ORIENTED ROUTINES ***************/
/**** Change these routines for your particular application ******/

input_app_parameters()
/* Input your application dependent parameters here and put the
 output in global variables */
{
}
app_computation()
/* this routine should contain any application-dependent computations */
/* that should be performed before each GA cycle.
   called by generate_new_pop    */
{
}

app_free()
/* application dependent free() calls, called by free_all() */
{
}

app_initialize()
/* application dependent initialization routine called by intialize() */
{
}

app_initreport()
/* Application-dependent initial report called by initreport() */
{
}

app_report()
/* Application-dependent report, called by report() */
{
}


app_statistics()
/* Application-dependent statistics calculations called by statistics() */
{
}

app_closure()
/* Application-dependent work which should be done before closure of
   the main program. called by main() */
{
}

int constr, dim;
double mx[8], x_new[8], u[8], d[8];
double dv[8]; /* dv is just a copy of the vector x*/
double rdesired = 0.00115;
double riskdesired = 0.0180;

double compute_g(double uu[8])
{
  int i;
  double g;
  for (i=0; i<dim; i++)
    x_new[i] = (uu[i]) * sigma[i] + mx[i];
  
  g = 0;	
  for(i=0;i<dim;i++){
	g += x_new[i]*dv[i];
  }	  
  g -= rdesired;
  
  return(g);
}

void compute_gradient()
{
  int i;
  double sum;
  for (i=0; i<dim; i++)
    x_new[i] = u[i]*sigma[i] + mx[i];
  
  switch(constr){
  case 0: //Return constraint
    for (i=0; i<dim; i++)
      {
	d[i] = -1*sigma[i]*dv[i];
      }
    break;
  case 1: // Risk constraint
    for (i=0; i<dim; i++)
      {
        d[i] = -1*sigma[i]*dv[i];
      }
    break;
  }
  sum=0.0;
  for (i=0; i<dim; i++)
    {
      sum = sum + d[i]*d[i];
    }
  
  sum = sqrt(sum);
  
  for (i=0; i<dim; i++)
    {
      if(sum!=0)
	d[i] = d[i]/sum;
    }
}

double getb(double x[8], int num)
{
  int i,j,k,iter;
  double g, old_u[8], beta, diff, min, ans;
  double zero[8], sumu, oldu[8];

  constr = num;
  beta = 1.00;
  dim = 8;
  
  for(i=0;i<dim;i++){
    mx[i] = R[i];
    zero[i] = 0;
    old_u[i] = 0;
    u[i] = 0;
  }
  
  iter = 1;
  
  do
    {
      compute_gradient();
      for (i=0; i<dim; i++)
	u[i] = beta * d[i];
      g = compute_g(u);
      diff = 0.0;
      for (i=0; i<dim; i++)
	diff += pow((u[i]-old_u[i]),2);
      diff = sqrt(diff);
      iter += 1;
      for (i=0; i<dim; i++)
	old_u[i] = u[i];
    }while ((diff > 0.00001)&&(iter<10));
  
  //printf("mpp point=%f\t%f\n",a[0],a[1]);
  ans=g;

  //printf("ans= %f\n",ans);
  return ans;
}

/*====================================================================
OBJECTIVE FUNCTION  ( Supposed to be minimized) :
Change it for different applications
====================================================================*/
void objective(indv)
INDIVIDUAL *indv;
{
  int i,j;
  double term1, term2, term3, pi, your_func;
  double g[MAXCONSTR], gsum, x[2*MAXVECSIZE];
      
  for (i=0; i < nvar_bin; i++)
    x[i] = indv->xbin[i];
  for (i=nvar_bin; i < nvar_bin+nvar_real; i++)
    x[i] = indv->xreal[i-nvar_bin];
  
#ifdef yours // define `yours' in the beginning of the code

  term1=0; term2=0; term3=2;
  x[7] = 1-x[0]-x[1]-x[2]-x[3]-x[4]-x[5]-x[6];
  
  /*x[0] = 0.0662; x[1] = 0.3776;
    x[2] = 0.0481; x[3] = 0.0726;
    x[4] = 0.0072; x[5] = 0.1498;
    x[6] = 0.2227; x[7] = 0.0558;*/ // these correspond to the MVP from Matlab code
  
   MINM = -1; // use -1 for maximization
   // Put your function here
   your_func = 0;
   for(i=0;i<8;i++){
     for(j=0;j<8;j++){
       term1 += x[i]*x[j]*V[i][j];
     }	
   }
   term1 = sqrt(term1);

   for(i=0;i<8;i++){
     term2 = term2 + x[i]*R[i]; 
   }

   your_func = 1000*term2 - term1; 
    
   for(i=0;i<8;i++){
     dv[i] = x[i];
   }	
   
   // Constraints :
   nc = 4;
   g[0] = x[7]; //So that all the weights are positive
   g[1] = term2 - rdesired; //Deterministic constraint on return 
   g[2] = riskdesired - term1; //Deterministic constraint on risk
   g[3] = getb(x,0); //Reliability return
   //g[4] = getb(x,1); //Reliability risk
			
#endif
   
   indv->obj = your_func;
   for (i=0, gsum=0.0; i<nc; i++)
     {
       indv->cons[i] = g[i];
       if (g[i] < 0.0) gsum += -1.0 * g[i];
     }
     indv->penalty = gsum;
}

/********************** END  OF  FILE **************************/

