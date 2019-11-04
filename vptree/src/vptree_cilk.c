#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include "quickselect.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>  
#include <string.h>

double *shared_points;
double *distances;
int threshold=8;


typedef struct vptree
{
	int idx;
	double *vantage_point;
	double median_value;
	struct vptree *bigger;
	struct vptree *smaller;
}vptree;


vptree* newnode(double *point ,int index ,int n , int d);
vptree* buildvp(double *X, int n, int d);
vptree * getInner(vptree * T);
vptree * getOuter(vptree * T);
double getMD(vptree * T);
double * getVP(vptree * T);
int getIDX(vptree * T);
void calculate_distances(vptree *node  , int start , int  end , int D );
vptree* calculate_median(vptree* ,int start , int end , int d  );
int check_equals( double *distances , int index , int  d);
void cilk_calculate_distances( vptree *node , int start , int  end , int D );

 
vptree* newnode(double *point ,int index ,int n , int d)
{
	vptree *node = (vptree*) malloc(sizeof(vptree));
	node->idx = index;
	node->vantage_point = (double *) malloc(sizeof(double)*d);
	for (int i = 0; i < d; i++)
	{
		node->vantage_point[i] = point[index*d+i];
	}
	node->median_value = n;
	node->bigger = node->smaller = NULL ;
	return node;
}

vptree* buildvp(double *X, int n, int d)
{

 
	vptree *vpt = (vptree *) malloc(sizeof(vptree));
	distances = (double *) malloc(n*sizeof(double));
	shared_points = (double *) malloc(n*d*sizeof(double));
	
	distances = (double *) malloc(n*sizeof(double));
	shared_points = X;

	/** typecasting to const char* **/
    char num_of_threads[3];
    sprintf(num_of_threads, "%d", threshold);
	    __cilkrts_end_cilk();
	__cilkrts_set_param("nworkers", num_of_threads ); 
	printf("Number of workers = %d\n", __cilkrts_get_nworkers());

	vpt = calculate_median(vpt,0,n-1,d);

	return vpt;
}

vptree * getInner(vptree * T)
{
	return T->smaller;
}

vptree * getOuter(vptree * T)
{
	return T->bigger;
}

double getMD(vptree * T)
{
	return T->median_value;
}

double * getVP(vptree * T)
{
	return T->vantage_point;
}

int getIDX(vptree * T)
{
	return T->idx;
}


//Recursevly calculate the vantage_point tree
vptree* calculate_median(vptree* node   ,int start , int end ,  int d  )
{

	int n = end -start + 1;

	//Initiate node with the selecetd(last point) vantage point
	node = newnode(shared_points,end, 0 , d);
	//End - start = 0 means that the is only one point int the array
	//so we are having a leaf
	if(end-start==0)
	{
		return node;
	}
	//In this occasion we have two points in the array. One of them will
	//be chosen as a vantage point,so it is going to be created a complete binary
	//which means there are gonna be nodes with only a right child node
	else if(end-start==1)
	{
		//Calculate the distances
		calculate_distances(node , start , end-1 , d);
		
		//Set the median as the distance of the first point
		double median = distances[start];
		node->median_value = median;
		
		//Set the last point as the vantage point of the inner-smaller part
		node->smaller = calculate_median(node->smaller,start , end-1  , d);

		return node;	
	}
	//In this a occasion we are having three points and the tree gonna be splitted in a perfect binary tree 
	else if(end-start==2)
	{

		calculate_distances(node , start, end-1 , d);
		
		if(distances[start]>distances[end-1])
		{
			//Set the median as the distance of the last point
			double median = distances[end-1];
			node->median_value = median;
			node->bigger = calculate_median(node->bigger,start , end-2  , d);
			node->smaller = calculate_median(node->smaller,start+1, end-1  , d);
			return node;
		}
		else if(distances[start]==distances[end-1])
		{
			double median = distances[start];
			node->median_value = median;
			node->smaller = calculate_median(node->smaller,start, end-1 , d);
			return node;
		}
		else
		{
			//Set the median as the distance of the first point
			double median = distances[start];
			node->median_value = median;

			node->bigger = calculate_median(node->bigger, start+1, end-1  , d);
			node->smaller = calculate_median(node->smaller  ,start , end-2  , d);
			return node;
		}
	}
	//General occasion
	else
	{
		//threshold-=1;
		//if(threshold>0 && ( (end-start)/threshold >= 12500000 || d*(end-start)/10000000 ))
		//	cilk_calculate_distances(node , start , end-1 , d);
		//else

		calculate_distances(node , start , end-1 , d);
		//Pick and and set the median to the node
		threshold+=2;
		double median = quickselect( distances , shared_points , start , end-1 ,start+(n-1)/2 , d);
		node->median_value = median;

		if(threshold>1)
		{
			threshold--;
			node->smaller =cilk_spawn calculate_median(node->smaller,start , start+(n-1)/2 ,  d);	
			
			node->bigger = calculate_median(node->bigger, start + (n-1)/2 +1 , end-1 , d);
			
			cilk_sync;
			threshold++;

		}
		else
		{
			node->smaller = calculate_median(node->smaller,start , start+(n-1)/2 ,  d);		
			node->bigger = calculate_median(node->bigger, start + (n-1)/2 +1 , end-1 , d);
		}
 		return node;
 	}	
}


void cilk_calculate_distances( vptree *node , int start , int  end , int D )
{

	cilk_for(int i=start ; i<=end; i++){
		double sum = 0;
		for(int j=0 ; j<D ;j++){
			sum+= pow(shared_points[i*D+j] - node->vantage_point[j],2);
			//printf("sum %lf \n", sum );
		}
		//printf("seg fault before here i %d \n" , i);
		distances[i] = sqrt(sum);
		//printf("Distance  : %d from  %d : %lf \n", i , getIDX(node) , distances[i] );
	}
}



//Calculating distances using for loops
void calculate_distances( vptree *node , int start , int  end , int D ){
	for(int i=start ; i<=end; i++){
		double sum = 0;
		for(int j=0 ; j<D ;j++){
			sum+= pow(shared_points[i*D+j] - node->vantage_point[j],2);
			//printf("sum %lf \n", sum );
		}
		//printf("seg fault before here i %d \n" , i);
		distances[i] = sqrt(sum);
		//printf("Distance  : %d from  %d : %lf \n", i , getIDX(node) , distances[i] );
	}
}