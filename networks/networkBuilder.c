#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define	N							4096
#define	DOUBLE_N					4096.0
#define	TAU_T						1.0
#define	RHO							0.72
#define	DIM							2

/* stuff that we don't touch usually */
#define	TIME_STEP						0.01
#define	HALF_TIME_STEP				0.5*TIME_STEP
#define	DELTA(a,b)					((a)==(b)?1.0:0.0)
#define	MAX_NEBZ						64
#define	MAX_CELLS					(N/2)
#define	X_COMP					0
#define 	Y_COMP 					1
#define	CUTOFF_SQRD					1.919383103666485
#define	C_0							-0.806140903539923
#define	C_2							0.7
#define	C_4							-0.156300219287607
#define	N_0							10.0
#define	DR							0.35

/*minimizer shit */
#define	MAX_ITERATIONS		500000
#define	TOL					1e-18	
#define	LINMIN_MAX_ITERATIONS	40
#define	LINMIN_G1			2.0	 // factors for growing and shrinking the interval- don't change
#define	LINMIN_G2			1.25
#define	LINMIN_G3			0.5
#define	LAST_X_DEFAULT		0.0001

/* constant constants... :P*******************************************************/
#define	LENGTH				1.0
#define	HALF_LENGTH			0.5

#define	TWO_PI					6.28318530717958
#define	UNI		((double)rand()/((double)RAND_MAX + 1.0))


void initializeSystem();
void updateNebzLists();
void calculateForces();
int min();

/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double T;			/*	temperature			*/
double t;			/*	time				*/
double u;			/*	potential energy 		*/
double kinetic;		/*	kinetic energy 		*/
double sxy;
double xi_T;
int contacts;
double coordination;

double rx[N];		/*	x component of position	*/
double ry[N];		/*	y component of position	*/
double fx[N];		/*	x component of position	*/
double fy[N];		/*	y component of position	*/
int type[N]; 		/*	for binary system	*/

double px[N];		/*	x component of momentum	*/
double py[N];		/*	y component of momentum	*/

int serial;

/******************************* for nebz list & cell subdivision *******************************/
int nebz[N][MAX_NEBZ]; //maximum we give here MAX_NEBZ nebz. 
int numOfNebz[N]; //number of elements in nebz[N][MAX_NEBZ]
double cellCutOff,maxD,listCutOffSqrd; //for updating nebz list.
int nebListCounter; //for counting how many steps have took place for nebz list updating.
int numInCell[MAX_CELLS];
int cells[MAX_CELLS][MAX_NEBZ];
const int shift[4][2] = {{-1,0},{-1,1},{0,1},{1,1}}; //{yShift,xShift};
/************************************************************************************/

/*******************************for minimizer*******************************************/
int restart;		/* whether to restart optimizer - fresh cg directions */
double lastx;		/* keeps track of typical step length */
double gtyp;		/* stores the rms gradient for linmin */
double p[DIM*N],q[DIM*N],g[DIM*N],h[DIM*N],mosesX[DIM*N],mosesY[DIM*N],scratch[DIM*N];
/************************************************************************************/

/*******************************for network*******************************************/
int networkNebz[N][MAX_NEBZ]; //maximum we give here MAX_NEBZ nebz. 
int numOfNetworkNebz[N]; //number of elements in nebz[N][MAX_NEBZ]
/************************************************************************************/

const double invSizeSqrd[2][2] = {{1.00,0.718184429761563},{0.718184429761563,0.510204081632653}}; // [0][0] = small-small; [0][1] = [1][0] = large-small; [1][1] = large-large;

int intSqrt(int number){
	int k=1;
	while ( k < number ){
		if ( k*k == number )
			return k;
		k++;
	}
	return 0;
}

void fixDrift(){
	int i;
	double Px,Py;
	Px = 0.0; Py = 0.0;
	for (i=0;i<N;i++){
		Px += px[i];
		Py += py[i];
	}
	Px = Px/DOUBLE_N; Py = Py/DOUBLE_N; 
	for (i=0;i<N;i++){
		px[i] -= Px;
		py[i] -= Py;
	}
	return;
}

//notice this initializes the strain increment to the default value... just so its something, but this can be changed.
void initializeSystem(){
	cellCutOff = 1.4*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	invL = 1.0/L;
	updateNebzLists();
	calculateForces();
	return;
}

void initializeGrid(){
	int i,sqrtN,j,k;
	double space,tempRx,tempRy;
	
	V = DOUBLE_N/RHO;
	sqrtN = intSqrt(N);
	L = sqrt(V);
	space = 1.0/(double)sqrtN;
	for (i=0; i<N; i++){
		rx[i] = space*(double)(i%sqrtN);
		ry[i] = space*(double)(i/sqrtN);
		px[i] = 0.0;
		py[i] = 0.0;
		type[i] = 0;
	}
	
	for (i=0; i<N/2; i++){
		j = (int)(UNI*DOUBLE_N);
		while ( type[j] )
			j = (int)(UNI*DOUBLE_N);
		type[j] = 1;
	}
	//permutate indices:
	for (i=0; i<N; i++){
		k = (int)(UNI*DOUBLE_N); //choose random index to switch with
		tempRx = rx[i];
		tempRy = ry[i];
		rx[i] = rx[k];
		ry[i] = ry[k];
		rx[k] = tempRx;
		ry[k] = tempRy;
	}
	initializeSystem();
	return;
}

void saveSnapShot(char *fileName){
	int i;
	FILE *file;
	file = fopen(fileName,"wb");
	fprintf(file,"%.15g\t%.15g\t%10g\n",L,u,sxy); //why not add strain in second place after L.
	for (i=0; i<N; i++)
		fprintf(file,"%.15g\t%.15g\t%d\n",rx[i],ry[i],type[i]);
	fclose(file);
	return;
}

//this is with strain
void readSnapShot(char *fileName){
	int i,k;
	double dummy;
	FILE *inFile;
	inFile = fopen(fileName, "rb");
	fscanf(inFile, "%lf",&L);
	V = L*L;
	fscanf(inFile, "%lf",&dummy);
	fscanf(inFile, "%lf",&dummy);
	for (i=0;i<N;i++)
		fscanf(inFile, "%lf %lf %d",&(rx[i]),&(ry[i]),&(type[i]));
	fclose(inFile);
	initializeSystem();
	return;
}

//this uses the cell subdivision.
void updateNebzLists(){
	int i,j,current,mSqrd,a,b,k,l,numHere,numThere,target,w,m;
	int ii,jj;
	double dx,dy,r2;
	double includeSqrd,invCellSize;
	
	m =(int)(L/cellCutOff);
	mSqrd = m*m; 		//the length of cells and numInCell...
	invCellSize = (double)m;
	
	//set number in each cell to zero
	for (i=0; i<mSqrd; i++)
		numInCell[i] = 0;
	
	/*********cells sorting********/
	for (i=0;i<N;i++){
		current = m*(int)(invCellSize*ry[i]) + (int)(invCellSize*rx[i]);
		cells[current][numInCell[current]] = i;
		numInCell[current]++;
	}

	//set all neb lists to zero
	for (i=0; i<N; i++)
		numOfNebz[i] = 0;
	
	//here's where it can get messy when sheering. maybe deal with that later.
	for (i=0;i<m;i++){
		for (j=0;j<m;j++){
			current = m*i + j;
			numHere = numInCell[current];
			//first check interactions within a cell
			for (k=0;k<numHere-1;k++){
				for (l=k+1;l<numHere;l++){
					ii = cells[current][k];
					jj = cells[current][l];
					dx = rx[jj] - rx[ii];
					dy = ry[jj] - ry[ii];
					r2 = V*( dx*dx + dy*dy );
					if (r2 < listCutOffSqrd){
						nebz[ii][numOfNebz[ii]] = jj;
						nebz[jj][numOfNebz[jj]] = ii;
						numOfNebz[ii]++;
						numOfNebz[jj]++;
					}
				}
			}
			
			//now there are 4 nebz that need to be checked.
			for (w=0;w<4;w++){ //w is an index... sorry.
				a = i + shift[w][0];
				if (a<0) a = m-1;
				else if (a==m) a = 0;
				b = j + shift[w][1];
				if (b==m) b = 0;
				target = m*a + b;
				numThere = numInCell[target];
				for (k=0;k<numHere;k++){
					for (l=0;l<numThere;l++){
						ii = cells[current][k];
						jj = cells[target][l];
						dx = rx[jj] - rx[ii];
						dy = ry[jj] - ry[ii];
						// mess due to periodic boundary conditions 
						if ( dx >= HALF_LENGTH )
							dx -= LENGTH;
						else if ( dx < -HALF_LENGTH )
							dx += LENGTH;
						if ( dy >= HALF_LENGTH )
							dy -= LENGTH;
						else if ( dy < -HALF_LENGTH )
							dy += LENGTH;
						// end of mess 
						r2 = V*( dx*dx + dy*dy );
						if (r2 < listCutOffSqrd){
							nebz[ii][numOfNebz[ii]] = jj;
							nebz[jj][numOfNebz[jj]] = ii;
							numOfNebz[ii]++;
							numOfNebz[jj]++;
						}
					}
				}
			}
		}
	}
	return;
}

void calculateForces(){
	int i,j,m,k,typeI;
	double g,dx,dy,rxi,ryi,r2OverSigma2,invSigmaSqrd;
	double sigmaOverR12,temp;
	
	//first set 'em all to zero
	u = 0.0; sxy = 0.0;
	for (i=0;i<N;i++){
		fx[i] = 0.0;
		fy[i] = 0.0;
	}
	
	for (i=0;i<N-1;i++){
		rxi = rx[i];
		ryi = ry[i];
		typeI = type[i];
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = rx[j] - rxi;
				dy = ry[j] - ryi;
				// mess due to periodic boundary conditions 
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				// end of mess 
				invSigmaSqrd = invSizeSqrd[typeI][type[j]];
				r2OverSigma2 = L*L*( dx*dx + dy*dy )*invSigmaSqrd;
				if (r2OverSigma2 < CUTOFF_SQRD){
					sigmaOverR12 = 1.0/(r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2);
					g = (N_0*sigmaOverR12 - r2OverSigma2*4.0*C_4  - 2.0*C_2)*invSigmaSqrd;
					temp = L*dx*g;
					fx[j] += temp;
					fx[i] -= temp;
					temp = L*dy*g;
					fy[j] += temp;
					fy[i] -= temp;
					u += r2OverSigma2*(sigmaOverR12 + C_2 + r2OverSigma2*C_4) + C_0;
					sxy += -temp*L*dx;
				}
			}
		}
	}
	sxy = sxy/V;
	return;
}


/****************************************************************minimizer functions***********************************************************/

void grad(double *point, double *res){
	int i,k,j,m,typeI;
	double rxi,ryi,dx,dy,invSigmaSqrd, halfL;
	double r2OverSigma2,sigmaOverR12,g,temp;
	
	halfL = 0.5*L;
	
	//set to zero
	for (i=0;i<DIM*N;i++)
		res[i] = 0.0;
	
	for (i=0;i<N-1;i++){
		rxi = point[DIM*i+X_COMP];
		ryi = point[DIM*i+Y_COMP];
		typeI = type[i];
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = point[DIM*j+X_COMP] - rxi;
				dy = point[DIM*j+Y_COMP] - ryi;
				// mess due to periodic boundary conditions 
				if ( dy >= halfL)
					dy -= L;
				else if ( dy < -halfL )
					dy += L;
				if ( dx >= halfL )
					dx -= L;
				else if ( dx < -halfL )
					dx += L;
				// end of mess 
				invSigmaSqrd = invSizeSqrd[typeI][type[j]];
				r2OverSigma2 = ( dx*dx + dy*dy )*invSigmaSqrd; //real space coordinates...
				if (r2OverSigma2 < CUTOFF_SQRD){
					sigmaOverR12 = 1.0/(r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2);
					g = (N_0*sigmaOverR12 - r2OverSigma2*4.0*C_4  - 2.0*C_2)*invSigmaSqrd;
					temp = dx*g;
					res[DIM*j+X_COMP] -= temp;
					res[DIM*i+X_COMP] += temp;
					temp = dy*g;
					res[DIM*j+Y_COMP] -= temp;
					res[DIM*i+Y_COMP] += temp;
				}
			}
		}
	}
	return;
}

//take care of periodic boundary conditions - move to model specific file 
//also check if neb lists need updating...
void fix(double *point){
	int i,xIndex,yIndex;
	double tmpd,dx,dy,length,halfLength;
	
	maxD = 0.0; length = L; halfLength = 0.5*L;
	
	for (i=0;i<N;i++){
		xIndex = DIM*i+X_COMP;
		yIndex = DIM*i+Y_COMP;
		
		if (point[yIndex] >= length)
			point[yIndex] -= length;
		else if (point[yIndex] < 0.0)
			point[yIndex] += length;
		if (point[xIndex] >= length)
			point[xIndex] -= length;
		else if (point[xIndex] < 0.0)
			point[xIndex] += length;
		
		//after fixing boundary conditions, I can now look for the maximal distance
		//traveled by some particle
		dx = point[xIndex] - rx[i]*L;
		dy = point[yIndex] - ry[i]*L;
		
		// mess due to periodic boundary conditions 
		if ( dy >= halfLength )
			dy -= length;
		else if ( dy < -halfLength )
			dy += length;
		if ( dx >= halfLength )
			dx -= length;
		else if ( dx < -halfLength )
			dx += length;
		
		// end of mess 
		tmpd = dx*dx + dy*dy;
		if ( tmpd > maxD )
			maxD = tmpd;
	}
	if ( 2.0*sqrt(maxD) > DR ){
		for (i=0;i<N;i++){
			rx[i] = point[DIM*i+X_COMP]*invL;
			ry[i] = point[DIM*i+Y_COMP]*invL;
		}
		updateNebzLists();
	}
	//that's it.
	return;
}

void initializeMinimizer(int start){ 
	int i;
	if ( start==0 ) 
		lastx = LAST_X_DEFAULT;
	for (i=0;i<DIM*N;i++){
		if ( restart !=2 )
			g[i] = -q[i];
		q[i] = h[i] = g[i];
	}
	restart = 0;
	return;
}

double prod(double x, double *gradScratch){
	int i;
	double s;
	for (i=0;i<DIM*N;i++)
		scratch[i] = p[i] + x*q[i];
	fix(scratch);
	grad(scratch, gradScratch);
	for (s = 0.0, i=0; i<DIM*N; i++)
		s += gradScratch[i]*q[i];
	return s;
}

double linmin(){
	double x,y,s,t,m,tmpd,step;
	double *gx,*gy,*gUnused;
	int it,i;
	
	gx = mosesX; gy = mosesY;
	
	x = lastx/gtyp;
	s = prod(x, gx); //result goes into gx
	it = 0;
	if ( s<0.0 ){
		do{
			y = x*LINMIN_G1;
			t = prod(y, gy);
			if ( t >= 0.0 ) break ;
			x = y; s = t; gUnused = gx; gx = gy; gy = gUnused; 
			it++;
		}while (it < LINMIN_MAX_ITERATIONS);
	}
	else if ( s>0 ){		
		do{
			y = x*LINMIN_G3;
			t = prod(y, gy);
			if ( t <= 0.0 ) break ;
			x = y ; s = t ; gUnused = gx ; gx = gy ; gy = gUnused;
			it++;
		}while (it < LINMIN_MAX_ITERATIONS);
	}
	else{		/* hole in one s = 0.0 */
		t = 1.0; y = x;
	}
	
	if ( it >= LINMIN_MAX_ITERATIONS){
		printf("Warning! linmin overran, ");
		tmpd = prod(0.0, gy);
		if ( tmpd > 0.0 )
			restart = 1;
	}
	
	if ( s < 0.0 ) s = - s;
	if ( t < 0.0 ) t = - t;
	m = ( s + t );
	s /= m; t /= m;
  
	m =  s * y + t * x;
	/* evaluate the step length, not that it necessarily means anything */
	for ( step=0.0, i=0; i<DIM*N; i++ ){
		tmpd = m*q[i];
		p[i] += tmpd; 			/* this is the point where the parameter vector steps */
		step += fabs(tmpd);
		q[i] = s*gy[i] + t*gx[i];	/* send back the estimated gradient in q (NB not like linmin) */
	}
	
	fix(p);			/* take care of periodic boundary conditions */
	lastx = m*LINMIN_G2*gtyp;
	return ( step / ((double)DIM*DOUBLE_N) ); 
}

int min(){
	int it,i,j,k;
	double gg, gam, dgg ,step, tmpd,maxGradSqr;
	
	lastx = LAST_X_DEFAULT;
	restart = 0;
	
	//initialize
	updateNebzLists();
	for (i=0;i<N;i++){
		p[DIM*i+X_COMP] = L*rx[i];
		p[DIM*i+Y_COMP] = L*ry[i];
	}
	
	grad(p, q);
	initializeMinimizer(1);
	
	//main loop
	for (it = 0; it < MAX_ITERATIONS; it++){
		maxGradSqr = 0.0; gg=0.0;
		for (i=0;i<DIM*N;i++){
			tmpd = g[i]*g[i];
			gg += tmpd;
			if (tmpd > maxGradSqr)
				maxGradSqr = tmpd;
		}
		
		gtyp = sqrt( gg / (DOUBLE_N*(double)DIM) );
		//fprintf(stderr, "\r");
		//fprintf(stderr, "iteration %d out of %d : gg = %6.3g tol = %6.3g: ", it , MAX_ITERATIONS , gg , TOL ) ;
		
		//printf("gtyp is %g\n",gtyp);
		if ( maxGradSqr < TOL ){
			for (i=0;i<N;i++){
				rx[i] = p[DIM*i+X_COMP]*invL;
				ry[i] = p[DIM*i+Y_COMP]*invL;
			}
			updateNebzLists();
			calculateForces();
			return 1;
		}
		
		step = linmin();
		//printf("here?\n");
		if (restart)
			grad(p, q);
		if (restart){
			printf("\nRestarting optimizer\n"); 
			initializeMinimizer(0);
		}
		else{
			for ( dgg=0.0, i=0; i<DIM*N; i++)
				dgg += ( q[i]+g[i] )*q[i];
			gam = dgg/gg;
			if (gam > 3.0){
				//printf("gamma was too big, restarting optimizer\n");
				restart = 1;
				initializeMinimizer(0);
			}
			else{//this is my addition, since <gam> can't be too large, it fucks everything up.
				for ( tmpd=0.0, i=0; i<DIM*N; i++){
					g[i] = -q[i];                /* g stores (-) the most recent gradient */
					q[i] = h[i] = g[i] + gam*h[i];
					/* h stores q, the current line direction */
					/* check that the inner product of gradient and line search is < 0 */
					tmpd -= q[i]*g[i];
				}
				if ( tmpd > 0.0 ){
					restart = 2; 	/* signifies that g[j] = -q[j] is already done */
					printf("Restarting optimizer (2)\n");
					initializeMinimizer(0);
				}
			}
		}
	}
	return 0;
}

/****************************************************************end of minimizer functions***********************************************************/


//runs, with thermostating.
void advanceTime(){
	int i;
	double temp,dx,dy,inst_T;
	
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		
		dy = TIME_STEP*py[i]; //real space displacement
		temp = ry[i]  + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;
		
		dx = TIME_STEP*px[i]; //real space displacement
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;
		
		temp = dx*dx + dy*dy;
		if (temp > maxD)
			maxD = temp;
	}
	
	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();
	
	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		kinetic += px[i]*px[i] + py[i]*py[i];
	}
	kinetic = 0.5*kinetic;
	
	/************************** thermostating part *************************************/
	inst_T = 2.0*kinetic/((double)DIM*DOUBLE_N);
	xi_T = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_T);
	for (i=0;i<N;i++){
		px[i] *= xi_T;
		py[i] *= xi_T;
	}
	/********************************************************************************/
	return;
}


//run the system for <duration>
void run(double duration){
	int k,steps;
	int t0,t1;	/*for timing purposes...*/
	//FILE *conFile = fopen("conserve.dat","wb");
	
	t0 = time(0);
	steps = (int)(duration/TIME_STEP);
	
	printf("Starting run at T = %g, u = %g\n\n",T,u);
	for (k=0; k<steps; k++){
		advanceTime();
		t = (double)k*TIME_STEP;
		//fprintf(conFile,"%.10g\t%.10g\t%.10g\n",t,kinetic,u);
		//~ if (!(k%256)){
			//~ fprintf(stderr,"\r");
			//~ fprintf(stderr,"%.10g\t%.10g\t%.10f",t,u,2.0*kinetic/((double)DIM*DOUBLE_N) );
		//~ }
		
		if ( !(k%65536) )
                        fixDrift(); 
	}
	
	t1 = time(0);
	printf("equilibration at T = %g - Time in seconds: %d\n\n\n",T, t1-t0);
	printf("Ended run T = %g, u = %g\n\n",T,u);
	//fclose(conFile);
	return;
}

void cool(double finalTemperature, double coolingRate){ //, char *coolingFileName ){
	int k,steps;
	
	steps = (int)((T - finalTemperature)/(coolingRate*TIME_STEP));
	
	k = 0;
	while ( T >= finalTemperature ){
		T -= coolingRate*TIME_STEP;
		advanceTime();
		k++;
	}
	return;
}


int countContacts(){
	int i,sum = 0;
	for (i=0;i<N;i++)
		sum += numOfNetworkNebz[i];
	return sum;
}

void getContactsStatistics(int *totalNumOfContacts, double *mean, double *flucs){
	int i,sum,sumSqrz;
	
	sum = 0; sumSqrz = 0;
	for (i=0;i<N;i++){
		sum += numOfNetworkNebz[i];
		sumSqrz += numOfNetworkNebz[i]*numOfNetworkNebz[i];
	}
	*totalNumOfContacts = sum;
	*mean = (double)sum/DOUBLE_N;
	*flucs = (double)sumSqrz/DOUBLE_N - (*mean)*(*mean);
	return;
}

void removeMaximallyConnectedBond(){
	int i,maxSum,minDifference;
	int m,k,j,iBond,jBond;
	int maxNebCoordination, nebWithMaxCoordination;
	int indexOfNebWithMaxCoordination;
	
	//first find the maxSum pair
	maxSum = 0; 
	for (i=0;i<N;i++){
		//find neighbor with largest coordination
		maxNebCoordination = 0; nebWithMaxCoordination = 0;
		m = numOfNetworkNebz[i];
		for (k=0;k<m;k++){
			j = networkNebz[i][k];
			if (j > i){ //go over every pair only once... 
				if ( numOfNetworkNebz[j] > maxNebCoordination ){
					maxNebCoordination = numOfNetworkNebz[j];
					nebWithMaxCoordination = k;
					indexOfNebWithMaxCoordination = j;
				}
			}
		}
		if ( m + maxNebCoordination >  maxSum )
			maxSum = m + maxNebCoordination;
	}
	
	//now find that pair that has maxSum, and the smallest difference.
	minDifference = 1000000;
	for (i=0;i<N;i++){
		//find neighbor with largest coordination
		m = numOfNetworkNebz[i];
		for (k=0;k<m;k++){
			j = networkNebz[i][k];
			if (j > i){ //go over every pair only once... 
				if ( numOfNetworkNebz[j] + m == maxSum ){
					if ( abs(numOfNetworkNebz[j] - m) < minDifference ){
						minDifference = abs(numOfNetworkNebz[j] - m);
						iBond = i;
						jBond = j;
					}
				}
			}
		}
	}
	
	//now remove bond
	m = numOfNetworkNebz[iBond];
	for (k=0;k<m;k++){
		j = networkNebz[iBond][k];
		if ( j == jBond ){
			if ( k < m-1 ) //not the last in the list, replace with last
				networkNebz[iBond][k] = networkNebz[iBond][m-1];
			numOfNetworkNebz[iBond]--; //removed last
			break;
		}
	}
	m = numOfNetworkNebz[jBond];
	for (k=0;k<m;k++){
		i = networkNebz[jBond][k];
		if ( i == iBond ){
			if ( k < m-1 ) //not the last in the list, replace with last
				networkNebz[jBond][k] = networkNebz[jBond][m-1];
			numOfNetworkNebz[jBond]--; //removed last
			break;
		}
	}
	return;
}
	

void makeNetwork(double z, int serialNumber){
	FILE *networkFile;
	char networkFileName[128];
	int i;
	double rxi,ryi,dx,dy,r2OverSigma2,invSigmaSqrd;
	int typeI,m,k,j;
	
	//set all neb lists to zero
	for (i=0; i<N; i++)
		numOfNetworkNebz[i] = 0;
	
	for (i=0;i<N-1;i++){
		rxi = rx[i];
		ryi = ry[i];
		typeI = type[i];
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = rx[j] - rxi;
				dy = ry[j] - ryi;
				// mess due to periodic boundary conditions 
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				// end of mess 
				invSigmaSqrd = invSizeSqrd[typeI][type[j]];
				r2OverSigma2 = L*L*( dx*dx + dy*dy )*invSigmaSqrd;
				if (r2OverSigma2 < CUTOFF_SQRD){
					networkNebz[i][numOfNetworkNebz[i]] = j;
					networkNebz[j][numOfNetworkNebz[j]] = i;
					numOfNetworkNebz[i]++;
					numOfNetworkNebz[j]++;
				}
			}
		}
	}
	
	while ( (double)countContacts()/DOUBLE_N > z )
		removeMaximallyConnectedBond();

	coordination = (double)countContacts()/DOUBLE_N;
	sprintf(networkFileName,"network_%d_%.3f_%.4d.dat",N,z,serialNumber);
	networkFile = fopen(networkFileName,"wb");
	fprintf(networkFile,"0.0");
	for (i=0;i<9;i++)
		fprintf(networkFile,"\t-1.1");
	fprintf(networkFile,"\n");
	
	for (i=0;i<N;i++){
		fprintf(networkFile,"%.15g\t%.15g",rx[i],ry[i]);
		for (k=0;k<numOfNetworkNebz[i];k++)
			fprintf(networkFile,"\t%d",networkNebz[i][k]);
		for (k=0;k<8-numOfNetworkNebz[i];k++)
			fprintf(networkFile,"\t-1");
		fprintf(networkFile,"\n");
	}
	fclose(networkFile);
	return;
}


int main(int n,char **inputStrings){
	int num,k,numOfNetworks;
	FILE *stuffFile;
	char networkFileName[128];
	double mean,flucs,targetCoordination;
	
	sscanf(inputStrings[1],"%d",&serial);
	sscanf(inputStrings[2],"%lf",&targetCoordination);
	sscanf(inputStrings[3],"%d",&numOfNetworks);
	
	//srand(time(0));
	srand(1000 + serial);
	
	//targetCoordination = 4.02;
	
	k=0;
	for (k=0;k<numOfNetworks;k++){
		initializeGrid();
		T = 1.00;
		run(200.0);
		T = 1e-3;
		run(50.0);
		min();
		makeNetwork(targetCoordination, serial*numOfNetworks + k);
		getContactsStatistics(&num, &mean, &flucs);
		printf("target coordination is %.15g, actual is %.15g\nflucs are %.10g\n\n",targetCoordination,coordination,flucs);
		printf("\nfinished number %d\n\n",k);
	}
	
	return 0;
}


