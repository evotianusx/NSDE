
#include "stdafx.h"

#include <iostream>
#include <ctime>
#include <complex>

#include "de.h"
#include <stdio.h>
using namespace std;
long  RND;
//double X  [MAXPOP][MAXDIM];			// Initial population
double Xo[MAXPOP][MAXDIM];				// Old population, for trial 
//double Xn[MAXPOP][MAXDIM];			// New population
//double Xs[MAXPOP][MAXDIM];			// Swap population

long double rnd_uni(long *idum);
double extern evaluate(int D, double Tmp[], long *FES, int Num); //Num is function number

clock_t Start = clock();

int _tmain(int argc, _TCHAR* argv[]) {
	int GenCount[23] = { 1500, 2000, 5000, 5000, 20000, 1500, 3000, 9000, 5000, 1500, 2000, 1500, 1500, 100, 4000, 100, 100, 100, 100, 200, 100, 100, 100 };
	double Low[23] = { -100, -10, -100, -100, -30, -100, -1.28, -500, -5.12, -32, -600, -50, -50, -65.536, -5, -5, -5, -2, 0, 0, 0, 0, 0 };
	double Upp[23] = { 100, 10, 100, 100, 30, 100, 1.28, 500, 5.12, 32, 600, 50, 50, 65.536, 5, 5, 15, 2, 1, 1, 10, 10, 10 };
	int DimCount[23] = { 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 2, 4, 2, 2, 2, 3, 6, 4, 4, 4 };
	int i, j, n, Trial, Func;	// counting variables    
	int r1, r2, r3;				// placeholders for random indexes    
	int D;						// Dimension of parameter vector      
	int NP = 100;					// number of population members    
	int MaxTrial = 1;

	double* result = new double[MaxTrial];	//store result in array 
	int Strategy = 7;			// choice parameter for screen output 
	int genCount, genMax;
	int Seed = 11;  // Seed for random number
	long  FES;					// number of function evaluations     
	double Tmp_Fit;				// Fitness of Trial Vector                  
	double UB;					// upper parameter bound              
	double LB;					// lower parameter bound              
	double Tmp[MAXDIM];			// Temporary Vector
	double Best[MAXDIM];		// Best Vector
	double TotBestfit = 0;		//Store the total best fitness from each trial
	double Fitness[MAXPOP];		// Fitness values                
	double F = 0.5, CR = 0.9;		// Parameters of DE            
	double Best_Fit;			// Fitness of the best   
	double stdev;				//store the stdev
	char ch;
	int stop;					//Stop Flag	
	int Best_index;
	FILE *pFileTXT;
	pFileTXT = fopen("Res.xls", "w"); //This will create a new file: delete the old data

	fprintf(pFileTXT, "Function\t");
	fprintf(pFileTXT, "Pop Size\t");
	fprintf(pFileTXT, "Dim\t");
	fprintf(pFileTXT, "Max Gen\t");
	fprintf(pFileTXT, "Max Trial\t");
	fprintf(pFileTXT, "Ave Time\t");
	fprintf(pFileTXT, "Crossover CR\t");
	fprintf(pFileTXT, "Scaling F\t");
	fprintf(pFileTXT, "Mean Best\t");
	fprintf(pFileTXT, "STD\n");
	fclose(pFileTXT);


	for (Func = 0; Func <2; Func++){

		genMax = GenCount[Func];//set variable of runs
		D = DimCount[Func];
		UB = Upp[Func];
		LB = Low[Func];

		printf("\n\nFunction %d", Func + 1);
	/*	printf("\n\Dimension %d", D = DimCount[Func]);

		printf("\n\Lower BOund %f , Upper Bound %f", LB = Low[Func],UB = Upp[Func]);*/
		stdev = 0;
		memset(Best, NULL, sizeof(Best));
		memset(Fitness, NULL, sizeof(Fitness));
		memset(result, NULL, sizeof(result));

		TotBestfit = 0;

		clock_t Start = clock();

		for (Trial = 0; Trial < MaxTrial; Trial++){

			RND = rand();
			//RND = -(long)Trial;		// initialization of rnd generator
			FES = 0;
			// Initialization
			for (i = 0; i < NP; i++) {
				for (j = 0; j < D; j++) {
					Xo[i][j] = LB + rnd_uni(&RND)*(UB - LB);
					//printf("\n\X %24.3E", Xo[i][j]);
				}
				Fitness[i] = evaluate(D, Xo[i], &FES, Func);
				//printf("\n\nFitness %20.3E", Fitness[i]);
			}

			Best_Fit = Fitness[0];
			//Best_Index = 0;
			for (i = 1; i < NP; i++) {
				if (Fitness[i] < Best_Fit) {
					Best_Fit = Fitness[i];
				}
			}

			//printf("\n\nFinish Populating, Now evaluating");

			// Iteration loop
			genCount = 0; // generation counter reset
			while ((genCount < genMax)) {
				genCount++;
				/*if (rnd_uni(&RND) < 0.1){ CR =rnd_uni(&RND); }
				if (rnd_uni(&RND) < 0.1){ F = 0.1 + rnd_uni(&RND)* 0.9; }*/
				for (i = 0; i < NP; i++) {
				/*	if (rnd_uni(&RND) < 0.1){ CR = rnd_uni(&RND); }
					if (rnd_uni(&RND) < 0.1){ F = 0.1 + rnd_uni(&RND)* 0.9; }*/
					// Pick a random population member 

					do {
						r1 = (int)(rnd_uni(&RND)*NP);
					} while (r1 == i);

					do {
						r2 = (int)(rnd_uni(&RND)*NP);
					} while ((r2 == i) || (r2 == r1));

					do {
						r3 = (int)(rnd_uni(&RND)*NP);
					} while ((r3 == i) || (r3 == r1) || (r3 == r2));

					if (Strategy == 7) {
						for (int k = 0; k < D; k++) {
							Tmp[k] = Xo[i][k];
						}

						n = (int)(rnd_uni(&RND)*D);
					
						// perform D binomial trials */
						for (int L = 0; L < D; L++) {
							//if (rnd_uni(&RND) < CR || n==(D-1)) {               // change at least one parameter
							if (rnd_uni(&RND) < CR) {
								Tmp[n] = Xo[r1][n] + F*(Xo[r2][n] - Xo[r3][n]);
							}
							n = (n + 1) % D;

						}
					}

					Tmp_Fit = evaluate(D, Tmp, &FES, Func);  // Evaluate new vector in Tmp[]
					//printf("%20.3E", Tmp_Fit);
					// improved objective function value?
					if (Tmp_Fit <= Fitness[i]) {
						Fitness[i] = Tmp_Fit;
						for (int k = 0; k < D; k++) { Xo[i][k] = Tmp[k]; }

						// Was this a new minimum?
						if (Tmp_Fit < Best_Fit) {
							Best_Fit = Tmp_Fit;
							Best_index = i;
							
							for (int k = 0; k < D; k++) { Best[k] = Tmp[k]; }
						}
					}

				}

			}
			//for (int k = 0; k < D; k++)
			//{
			//	printf("\nX %d: %20.15E", k,Xo[i][k]);
			//	/*printf("\nX1: %20.15E", Xo[i][k]);
			//	printf("\nX2: %20.15E", Xo[i][k]);
			//	printf("\nX3: %20.15E", Xo[i][k]);*/
			//}

			
			//end while Gen < GenMax
			printf("\nBest Fitness: %20.5e", Best_Fit);

			TotBestfit = TotBestfit + Best_Fit;
			result[Trial] = Best_Fit;

			//clear arrays for next trial
			memset(Fitness, NULL, sizeof(Fitness));
			memset(Best, NULL, sizeof(Best));
			memset(Tmp, NULL, sizeof(Tmp));

		} //Calculate avg and STDEV
		double avgbestfit = TotBestfit / MaxTrial;

		double sum2 = 0;

		if (MaxTrial > 1)
		{
			for (int i = 0; i < MaxTrial - 1; i++){
				sum2 += pow((result[i] - avgbestfit), 2);
			}

			stdev = sqrt(sum2 / (MaxTrial - 1));
		}
		else
		{
			stdev = 0;
		}


		pFileTXT = fopen("Res.xls", "a");
		fprintf(pFileTXT, "%d\t", Func + 1);	//Function 
		fprintf(pFileTXT, "%d\t", NP);			//Population
		fprintf(pFileTXT, "%d\t", D);			//Dimension 
		fprintf(pFileTXT, "%d\t", genMax);		//Generation
		fprintf(pFileTXT, "%d\t", MaxTrial);	//Max Trial 
		fprintf(pFileTXT, "%d\t", (clock() - Start) / MaxTrial / 1000); //AveTime
		fprintf(pFileTXT, "%20.5e\t", CR);		//CR
		fprintf(pFileTXT, "%20.5e\t", F);		//F
		fprintf(pFileTXT, "%20.5e\t", avgbestfit); //Mean Best
		fprintf(pFileTXT, "%20.5e\n", stdev);	//STD
		fclose(pFileTXT);
		sum2 = 0;
		avgbestfit = 0;
		TotBestfit = 0;
	}
	/*double X[29]; double mul = 1.0;
	for (i = 0; i < 30; i++) {
		X[i] = i + 1;
		mul = mul * X[i];

	}*/
	printf("\n\nSimulation is completed");
	//printf("%20.5e\n", mul);
	cin.get(ch);

	return(0);

}


//Copy from vector b to vector a
void CopyVector(double a[], double b[]) {
	for (int k = 0; k<MAXDIM; k++) {
		a[k] = b[k];
	}
}

//Copy from scr to dest
int CopyArray(double dest[MAXPOP][MAXDIM], double src[MAXPOP][MAXDIM]) {
	for (int j = 0; j<MAXPOP; j++) {
		for (int k = 0; k<MAXDIM; k++) {
			dest[j][k] = src[j][k];
		}
	}
	return 0;
}
