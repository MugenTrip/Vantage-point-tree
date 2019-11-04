#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include "quickselect.h"
#include <omp.h>

double *shared_points;
double *distances;
int nthreads=8;
volatile int threshold=8;

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
void delete_distances(double *distances);
void p_calculate_distances(vptree *node  , int start , int  end , int D );
vptree* p_calculate_median(vptree* node   ,int start , int end ,  int d  );


 
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

	omp_set_dynamic(0);
	omp_set_nested(1);
	omp_set_num_threads(threshold);
	#pragma omp parallel
	{
		#pragma omp master
		{
			//printf("Master id : %d \n" ,   omp_get_thread_num());
			vpt = p_calculate_median(vpt,0,n-1,d);
		}
	}
	
	free(distances);
	
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
vptree* p_calculate_median(vptree* node   ,int start , int end ,  int d )
{

	//printf("distances %lf \n" ,distances[0]);
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
		//printf("7 .median Input parameters, start %d , end %d \n", start+1 , end  );
		calculate_distances(node , start , end-1 , d);
		
		//Set the median as the distance of the last point
		double median = distances[start];
		node->median_value = median;
		
		//Set the last point as the vantage point of the inner-smaller part
		//printf("7 .median Input parameters, start %d , end %d \n", start+1 , end  );
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
			//printf("5 .median Input parameters, start %d , end %d \n", start+1 , end-1  );
			node->bigger = calculate_median(node->bigger,start , end-2  , d);
		    //printf("6 .median Input parameters, start %d , end %d \n", start+2 , end  );
			node->smaller = calculate_median(node->smaller,start+1, end-1  , d);
			return node;
		}
		else
		{
			//Set the median as the distance of the first point
			double median = distances[start];
			node->median_value = median;

			//printf("3 .median Input parameters, start %d , end %d \n", start+2 , end  );
			node->bigger = calculate_median(node->bigger, start+1, end-1  , d);
			//printf("4 .median Input parameters, start %d , end %d \n", start+1 , end-1  );
			node->smaller = calculate_median(node->smaller  ,start , end-2  , d);
			return node;
		}
	}
	//General occasion
	else
	{

		calculate_distances(node , start , end-1 , d);
			
		//Pick and and set the median to the node
		//printf("Qs arguments start: %d , end %d , half %d \n", start , end-1 , start+(n-1)/2 );
		double median = quickselect( distances , shared_points , start , end-1 ,start+(n-1)/2 , d);
		node->median_value = median;
		
	
		if (threshold>1)
		{	
			//omp_set_num_threads(threshold);
			//printf("max active level %d and threshold %d \n",omp_get_max_active_levels() , threshold );
			#pragma omp task
			{
				#pragma omp critical
				{
					threshold--;
				}	
						
				node->smaller = p_calculate_median(node->smaller,start , start+(n-1)/2 , d );

				#pragma omp critical
				{
					threshold++;
				}
			}

				
			node->bigger = calculate_median(node->bigger, start + (n-1)/2 +1  , end-1 , d);

				//printf("WTF 1");
			#pragma omp taskwait
			return node;
		}
		else
		{
			
			
			//calculate the distances in parallel
			//if((end-start)/nthreads>125000 || d*(end-start)/nthreads> 500000 ){
				//omp_set_num_threads(nthreads);
				//p_calculate_distances(node , start , end-1 , d);
			//}
			node->smaller = calculate_median(node->smaller,start , start+(n-1)/2 , d);
			node->bigger = calculate_median(node->bigger, start + (n-1)/2 +1  , end-1 , d);

 			return node;
 		}
 	}	
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
		
		//Set the median as the distance of the last point
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

		calculate_distances(node , start , end-1 , d);
			
		//Pick and and set the median to the node
		double median = quickselect( distances , shared_points , start , end-1 ,start+(n-1)/2 , d);
		node->median_value = median;
		
	
		node->smaller = calculate_median(node->smaller,start , start+(n-1)/2 , d);
		
		node->bigger = calculate_median(node->bigger, start + (n-1)/2 +1  , end-1 , d);

 		return node;
				
		//recurse for the inner and outer part
		node->smaller = calculate_median(node->smaller,start , start+(n-1)/2 , d);
		node->bigger = calculate_median(node->bigger, start + (n-1)/2 +1  , end-1 , d);

 		return node;
 	}	
}

//Calculating distances using for loops
void calculate_distances( vptree *node , int start , int  end , int D ){
	for(int i=start ; i<=end; i++){
		double sum = 0;
		for(int j=0 ; j<D ;j++){
			sum+= pow(shared_points[i*D+j] - node->vantage_point[j],2);
		}
		distances[i] = sqrt(sum);
	}
}




//Calculating distances using for loops
void p_calculate_distances( vptree *node , int start , int  end , int D )
{
		#pragma omp parallel for 
		for(int i=start ; i<=end; i++){
			double sum = 0;
			for(int j=0 ; j<D ;j++){
				sum+= pow(shared_points[i*D+j] - node->vantage_point[j],2);
			}
			distances[i] = sqrt(sum);
		}

}


void delete_distances(double *distances)
{
	free(distances);
}

