///////////////////////////////////////////////////////////////////////////////////////////
// EE4371: Data Structures and Algorithms
// Implementation Code for Question 1: Modified Knapsack Problem
// Author: Rahul Chakwate, AE16B005
////////////////////////////////////////////////////////////////////////////////////////////

// Please save the input1.txt file in the same directory as the code.
// To compile from Linux terminal, use "gcc AE16B005.c -lm" command.
// Total time to run the code is about 16 seconds for i7 processor.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


//global variables
int n=0; 			//total number of weights
int w[1000]; 		//weights array
double v[1000]; 	//value array, here v=log(w)
int W0 = 10000; 	//weight constraint
int W; 				//weight constraint times root of weights used
int weightsum=0; 	//sum of the weights used
int Woptimal; 		//optimal W for ordinary knapsack m(n,W)

#define WIDTH 350000 
#define LENGTH 1000 

double *value[LENGTH];  	//array used for dynamic programming to store obj. function value
int *numweight[LENGTH]; 	//array used to store the number of weights used upto that n,W

////////////////////////////////////////////////////////////////////////////////////////////////

//to read the string input1.txt into array
void read_ints (const char* file_name)
{
  FILE* file = fopen (file_name, "r");
  int i = 0;

  fscanf (file, "%d", &i);    
  while (!feof (file))
    {  
      n++;
      w[n]=i;
      fscanf (file, "%d", &i);      
    }
  fclose (file);        
}


void knapsack();  	//declaration
double m(int, int); //declaration


int main () 
{ 
	read_ints("input1.txt");
	for(int i=1;i<=n;i++) v[i] = log(w[i]) / log(100);
	printf("Input read from input1.txt file\n...\n");

	double time_spent = 0.0;    					// initialize time t=0
	clock_t begin = clock();
	
	W  = W0*sqrt(n); 			// W cannot exceed this value
	for(int i=0; i<LENGTH; i++) numweight[i] = (int *)malloc(WIDTH * sizeof(int)); //bigger arrays cannot be loaded directly into memory  

	for(int i=0; i<LENGTH; i++) value[i] = (double *)malloc(WIDTH * sizeof(double)); //hence array of pointers is used with malloc()


    knapsack();
    printf("\n");

	clock_t end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;   	// calculate time spent
	printf("Total Time elpased is %f seconds.\n", time_spent);
    return 0; 
} 


// This function recursively computes the entries of value and numweight arrays using dynamic programming
//An entry in the value array either comes from the entry above value[i-1][j] or from value[i-1][j-w[i]].
double m(int i, int j)
{
	if(i==0 || j<=0) return 0.0;  					//if out of bounds return 0
	if(value[i-1][j]==-1) value[i-1][j]=m(i-1,j); 	//if above entry unknown, first compute that entry
	if(w[i]>j)  									//if weight > W, the entry comes from the entry above
	{
		value[i][j]=value[i-1][j];
		numweight[i][j]=numweight[i-1][j];
	}
	else
	{
		if(value[i-1][j-w[i]]==-1) value[i-1][j-w[i]]=m(i-1,j-w[i]); //if j-w[i]th above entry unknown, first compute that entry.
		if(value[i-1][j]>value[i-1][j-w[i]]+v[i])  		//i,j th entry is the max one among the two options
		{
			value[i][j] = value[i-1][j];
			numweight[i][j] = numweight[i-1][j];
		}
		else if(value[i-1][j]<value[i-1][j-w[i]]+v[i])
		{
			value[i][j] = value[i-1][j-w[i]] + v[i];
			numweight[i][j] = numweight[i-1][j-w[i]] + 1;
		}
		else  											//the ties are broken by taking the entry using more number of weights 
														//so that the constraint is more relaxed
		{
			if(numweight[i-1][j]>numweight[i-1][j-w[i]]+1)
			{
				value[i][j] = value[i-1][j];
				numweight[i][j] = numweight[i-1][j];				
			}
			else
			{
				value[i][j] = value[i-1][j-w[i]] + v[i];
				numweight[i][j] = numweight[i-1][j-w[i]] + 1;				
			}
		}

	}
	return value[i][j];
}


void knapsack()
{
	for(int i=0;i<=n;i++) //initialization loop
	{
		for(int j=0;j<=W;j++) 
		{
			value[i][j]=-1.0;
			numweight[i][j]=0;
		}
	}
	for(int Wk=W;Wk>0;Wk--)  //the maximum W can go is 10000 * root(n). Hence start from that value and go back in the table.
	{
		weightsum=0;

		double res = m(n,Wk); //computes ordinary knapsack for Wk constraint and using n weights
		int j = Wk;

		for(int i=n;i>0 && res>0;i--)  //loop to backtrack the weights
		{
			if(w[i]>j) continue; 		// if weight greater than constraint, ignore it
			if(value[i-1][j]>value[i-1][j-w[i]]+v[i]) continue; //if data entry came from above entry, no weight is added.
			else if(value[i-1][j]<value[i-1][j-w[i]]+v[i]) //if data entry came from a new weight, print the weight.
			{
				res-=v[i];  		//subtract the weight from the residual till it is zero
				if(res<0) continue;
				else 
				{
					weightsum+=w[i];
				}
				j-=w[i]; 			//subtract weight from j pointer
			}

			else
			{
				if(numweight[i-1][j]>numweight[i-1][j-w[i]]+1) continue; //ties are broken using number of weights used to form that entry
				else
				{
					res-=v[i];
					if(res<0) continue;
					else 
					{
						weightsum+=w[i];
					}
					j-=w[i];
				}
			}	

		
		}
		if(weightsum<=(W0*sqrt(numweight[n][Wk])))  	//check if modified knapsack constraint is also satisfied.
		{
			printf("\nOptimum solution found at Wk=%d\n",Wk); //if satisfied, this is the optimal solution as all other larger solutions are checked before.
			Woptimal=Wk;
			break;
		}
	}
	
		int Wk=Woptimal;  		//set Wk to optimal and repeat the above to print the optimal weights and obj function value.
		weightsum=0;
		printf("\nFor W = %d\n",Wk);
		double res = m(n,Wk);
		int j = Wk;
		printf("Maximum objective function value = %lf\n",res);
		printf("\nThe weights which maximize are \n");

		for(int i=n;i>0 && res>0;i--)
		{
			if(w[i]>j) continue;
			if(value[i-1][j]>value[i-1][j-w[i]]+v[i]) continue;
			else if(value[i-1][j]<value[i-1][j-w[i]]+v[i])
			{
				res-=v[i];
				if(res<0) continue;
				else 
				{
					printf("%d, ",w[i]);
					weightsum+=w[i];
				}
				j-=w[i];
			}

			else
			{
				if(numweight[i-1][j]>numweight[i-1][j-w[i]]+1) continue;
				else
				{
					res-=v[i];
					if(res<0) continue;
					else 
					{
						printf("%d, ",w[i]);
						weightsum+=w[i];
					}
					j-=w[i];
				}
			}	

		
		}
		printf("\n\nNumber of Weights=%d, Sum of Weights=%d \n",numweight[n][Wk],weightsum);

		return;
	
}

