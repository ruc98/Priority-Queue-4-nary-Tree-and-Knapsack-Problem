///////////////////////////////////////////////////////////////////////////////////////////
// EE4371: Data Structures and Algorithms
// Implementation Code for Question 2
// Author: Rahul Chakwate, AE16B005
////////////////////////////////////////////////////////////////////////////////////////////


// To compile from Linux terminal, use "gcc AE16B005.c -lm" command.
// NOTE: 1000 iterations are run on a server with 220GB RAM.
// NOTE: 90-100 iterations requires 12-13GB of RAM. The algorithm is not memory optimised yet.

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 		//use -lm to compile
#include <time.h>


struct particle // struct designed to store the particle information 
{
	double x,y,vx,vy; // position: x,y velocity: vx,vy
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//global variables
#define PARTICLES 1000000 		// Total number of particles
#define BITS 7  				//1 for sign, 6 for fraction representation

struct particle p[PARTICLES]; 	//particles array
double dt = 0.0001; 			//timestamp
int num_of_times=50;			//number of timestamps, change if required
int nbidx[500]; 				//stores neighbor indices, change wisely according to PARTICLES and BITS.
int i,j,k;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//tree structures
struct leafnode{ 				//leaf of the 4nary tree, stores particle id
  int particleidx;
};

typedef struct leafnode leafnode;


struct treenode{ 				//non leaf treenode which stores pointer to 4 children and leaves if any
  int numleaves;
  struct treenode *child00;
  struct treenode *child01;
  struct treenode *child10;
  struct treenode *child11;
  struct leafnode *leaves[500];
};

typedef struct treenode treenode;

/////////////////////////////////////////////////////////////////////////////////
//tree utility functions

treenode *createtreenode(){ 			//create and initialize treenode
   treenode *new_node;
   new_node=malloc(sizeof(treenode));
   new_node->child00=NULL;
   new_node->child01=NULL;
   new_node->child10=NULL;
   new_node->child11=NULL;
   new_node->numleaves=0;
   new_node->leaves[0]=NULL;

   return new_node;
}

treenode *addleaf(treenode *current_node,int particleidx) 		//adds a leaf with particleidx to current_node 
{
  int leafidx = current_node->numleaves;
  current_node->leaves[leafidx] = malloc(sizeof(leafnode *));
  current_node->leaves[leafidx]->particleidx = particleidx;
  current_node->numleaves = leafidx+1;
}

treenode *deleteleaf(treenode *current_node) 					// free leaves memory
{
  int totalleaves = current_node->numleaves;
  if(totalleaves!=0){
  for(int i=totalleaves-1;i>=0;i--) 
  {
    free(current_node->leaves[i]);
  }
  current_node->numleaves=0;
  }
}


treenode *gotochild(treenode *current_node, int bx, int by) 	//directs pointer towards child by looking at binary bits
{
  switch (2*bx+by)
  {
  case 0:
      {
        if(current_node->child00) 				
        {
          current_node = current_node->child00;
        }
        else 													//if child deosn't exist, create childnode
        {
          current_node->child00 = malloc(sizeof(treenode *));
          current_node->child00=createtreenode();
          current_node = current_node->child00;
        }
      }
      break;
  case 1:
      {
        if(current_node->child01)
        {
          current_node = current_node->child01;
        }
        else
        {
          current_node->child01 = malloc(sizeof(treenode *));
          current_node->child01=createtreenode();
          current_node = current_node->child01;
        }
      }
      break;
  case 2:
      {
        if(current_node->child10)
        {
          current_node = current_node->child10;
        }
        else
        {
          current_node->child10 = malloc(sizeof(treenode *));
          current_node->child10=createtreenode();
          current_node = current_node->child10;
        }
      }
      break;
  case 3:
      {
        if(current_node->child11)
        {
          current_node = current_node->child11;
        }
        else
        {
          current_node->child11 = malloc(sizeof(treenode *));
          current_node->child11=createtreenode();
          current_node = current_node->child11;
        }
      }
      break;
  }
  return current_node;
} 

int num_of_children; 									//variables for getchildren()
int nbcount;
int getchildren(treenode *current_node,int pidx) 		//updates neighborindex and returns #children
{

  num_of_children = current_node->numleaves;
  nbcount=0;
  for(k=0;k<num_of_children;k++)
  {
    if(current_node->leaves[k]->particleidx==pidx) continue;
    nbidx[nbcount] = current_node->leaves[k]->particleidx;
    nbcount++;
  }
  return num_of_children;
}

void dec_to_bin(long double dFractional,  int *bits) 	//converts fractional decimal to binary
{
  int temp;

  if(dFractional<0)
  {
    dFractional+=1;
    bits[0]=1;
  }
  else bits[0]=0;

  for(int i=1;i<BITS;i++)
  {
    
    dFractional = dFractional * 2;
    temp =  dFractional;
    bits[i]=temp;
    if(temp ==1)dFractional = dFractional - temp;
  }
}

treenode *updatetree(treenode *root,int key) 		//creates the tree and adds particles as leaves
{

  treenode *current_node;
  int bx[BITS],by[BITS];

  for(k=0;k<PARTICLES;k++)
  {
    current_node=root;
    dec_to_bin(p[k].x,bx);
    dec_to_bin(p[k].y,by);

    for(i=0;i<BITS;i++)
    {
      current_node = gotochild(current_node,bx[i],by[i]);
    }

    if(key==1) addleaf(current_node,k);
    else deleteleaf(current_node);

  }  

  return root;
}

int getneighboridx(int k, treenode *root)  					// directs pointer to correct leaf
{
	treenode *current_node;
  int bx[BITS],by[BITS];
	current_node=root;
  dec_to_bin(p[k].x,bx);
  dec_to_bin(p[k].y,by);
	for(i=0;i<BITS;i++)
	{
	current_node = gotochild(current_node,bx[i],by[i]);
	}
	return getchildren(current_node,k) - 1; 				// also updates global nbidx array.

}
//////////////////////////////////////////////////////////////////////////////////////////////////

// general functions
double rand_double(double a, double b) 				//creates random double number within a range
{
    return a + (b-a)*rand()/(double)RAND_MAX;
}

double euclidean_distance(int i, int k)				//computes euclidean distance
{ 
	return sqrt(pow((p[i].x - p[k].x),2) + pow((p[i].y - p[k].y),2));
}
/////////////////////////////////////////////////////////////////////////////////////////////////

// Partial QuickSort utility functions and variables
#define NBINS 100 								//#of bins to segregate the velocities
double pivot; 									//pivot for quicksort
double dp; 										//increment in pivot value
int sortcount=0;  
int partitionidx[NBINS]; 						//stores partition indices of vmag[] array
int binstrength[NBINS];							//stores bin strengths of 100 bins
double binaverage[NBINS]; 						//stores bin averages of 100 bins

void swap(double* a, double* b)  				//swaps in place
{ 
    double t = *a; 
    *a = *b; 
    *b = t; 
} 

double maxarr(double arr[], int n) 				//compute max element in an array
{ 
    double max = arr[0];  
    for (int i = 1; i < n; i++) 
        if (arr[i] > max) 
            max = arr[i]; 
  
    return max; 
} 

double minarr(double arr[], int n)  			//computes min element in an array
{ 
    double min = arr[0];  
    for (int i = 1; i < n; i++) 
        if (arr[i] < min) 
            min = arr[i]; 
  
    return min; 
}

double avgarr(double arr[], int low, int high) 	//computes average element in an array
{
  double sum=0.0;
  for(int idx=low;idx<high;idx++) sum+=arr[idx];
  return sum / (double)(high-low);
}


void quickSort(double arr[], int low, int high)  //recursively calls quicksort NBINS times for partial sorting
{ 
    if (low >= high || sortcount>=NBINS) return; 
    
    pivot+=dp;
    int i = (low - 1); 							// Index of smaller element 

    for (int j = low; j <= high-1; j++) 
    { 
        										// If current element is smaller than the pivot 
        if (arr[j] < pivot) 
        { 
            i++; 								// increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    int pi = i + 1;  							// pi is partitioning index
    partitionidx[sortcount] = pi;
    sortcount++;
    
    

    quickSort(arr, pi, high); 

} 



void dist(double *vmag) 						//computes velocity distribution and updates binstrength[]
{ 												//and binavg[] arrays
    // pivot = minarr(vmag,PARTICLES);
    // dp = (maxarr(vmag,PARTICLES) - minarr(vmag,PARTICLES))/NBINS;
    pivot = 0;
    dp = 0.1; 									//denots bin size, min to max does not give good results
    sortcount=0;
    quickSort(vmag, 0, PARTICLES-1);

    binstrength[0]=partitionidx[0];
    binaverage[0] =avgarr(vmag,0,partitionidx[0]);
    for(int i=1;i<NBINS;i++)
    {
      binstrength[i] = partitionidx[i] - partitionidx[i-1];
      binaverage[i] = avgarr(vmag,partitionidx[i-1],partitionidx[i]);
    }


}


/////////////////////////////////////////////////////////////////////////////////////////////////

//main function
int main () 
{ 
	// initialization
	for(int i=0;i<PARTICLES;i++)  //loop to initialize the particles
	{
		p[i].x = rand_double(-1,1);
		p[i].y = rand_double(-1,1); 
	
		if(i%2)
		{
			double vdir = rand_double(0,2* M_PI);
			p[i].vx = cos(vdir);
			p[i].vy = sin(vdir);
		}
		else
		{
			p[i].vx=0;
			p[i].vy=0;
		}

	}
	printf("initialization complete\n");

	////////////////////////////////////////////////////////////////////////////////////////
	// simulation

	double time_spent = 0.0;    					// initialize time t=0
	clock_t begin = clock();
	struct particle pn;

	for(int t=0;t<num_of_times;t++) 				// time loop
	{
		printf("t = %d\n",t);
	    treenode *root; 							//define a root for the tree
	    root=createtreenode();
	    root = updatetree(root,1); 					// argument 1:addleaf, -1:deleteleaf
		for(int i=0;i<PARTICLES;i++) 				// loop to compute neighbors, forces and update velocities
		{
			int num_of_nb = getneighboridx(i,root); //updates neighbor index array globally
			for(int k=0;k<num_of_nb;k++) 			//loop to iterate through neighbors
			{
				pn=p[nbidx[k]];
				double denom= pow(euclidean_distance(i,nbidx[k]),3) * 10000;
		        if(denom==0)denom = 0.000000001;
				double Fx = (p[i].x - pn.x) / denom;
				double Fy = (p[i].y - pn.y) / denom;
		        p[i].vx = p[i].vx + Fx * dt;
		        p[i].vy = p[i].vy + Fy * dt;
			}
		}
    updatetree(root,-1);   							//deallocation of memory, for memory optimization 
    for(int i=0;i<PARTICLES;i++) 					// loop to update positions
		{
			p[i].x += p[i].vx * dt;
			p[i].y += p[i].vy * dt;
		    while(p[i].x<-1 || p[i].x>1 || p[i].y<-1 || p[i].y>1) //check boundary conditions
      {
  			if(p[i].x<-1) p[i].x = p[i].x - 2*(p[i].x+1);
  			if(p[i].x> 1) p[i].x = p[i].x - 2*(p[i].x-1);
  			if(p[i].y<-1) p[i].y = p[i].y - 2*(p[i].y+1);
  			if(p[i].y> 1) p[i].y = p[i].y - 2*(p[i].y-1);
      }

		}

    if(t%10==0)  //at regular intervals print time ellapsed and velocity histogram bins
    {
      printf("Simulation checkpoint\n");
      clock_t mid = clock();
      time_spent = (double)(mid - begin) / CLOCKS_PER_SEC;   	// calculate time spent
      printf("Time elpased is %f seconds.\n", time_spent);

      double vmag[PARTICLES];
      for(int i=0;i<PARTICLES;i++) vmag[i] = sqrt((p[i].vx*p[i].vx) + (p[i].vy*p[i].vy)); //compute magnitude

      dist(vmag); 												//updates binstrength[] and binaverage[] arrays
      printf("\nBinstrength\n");
      for(int i=0;i<NBINS;i++) printf("%d, ",binstrength[i]);
      printf("\nBinaverage\n");
      for(int i=0;i<NBINS;i++) printf("%lf, ",binaverage[i]);
	  printf("\n");
    }

	}

  printf("\nDone, Simulation Complete!\n");
}