#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include "quickselect.h"


typedef struct vptree
{
	int idx;
	double *vantage_point;
	double median_value;
	struct vptree *bigger;
	struct vptree *smaller;
	//bool leaf;
}vptree;

typedef struct data 
{
	struct vptree* node;
	int start;
	int end;
	int previous_end;
	int d;

}data;



struct vptree* newnode(double *point ,int index ,int n , int d);
struct vptree* buildvp(double *X, int n, int d);
struct vptree * getInner(struct vptree * T);
struct vptree * getOuter(struct vptree * T);
double getMD(struct vptree * T);
double * getVP(struct vptree * T);
int getIDX(struct vptree * T);
void calculate_distances(struct vptree *node , int start , int  end , int D );
struct vptree* calculate_median(struct vptree* node   ,int start , int end ,  int d  );
void *parallel_calculate_median(void *arg);
void* p_calculate_distances(void *arg);


int nthreads=4;
double *shared_points;
pthread_mutex_t mute ,*wait;
int volatile threshold = 8;
double *distances;




struct data* Initiate_data(  struct vptree* vpnode , int new_start , int new_end , int new_d)
{
	struct data *new_data = (struct data*) malloc(sizeof(struct data));
	new_data->node= vpnode;
	new_data->start = new_start;
	new_data->end = new_end;
	new_data->d = new_d;
	return new_data;
}


struct vptree* newnode(double *point ,int index ,int n , int d)
{
	struct vptree *node = (struct vptree*) malloc(sizeof(struct vptree));
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

struct vptree* buildvp(double *X, int n, int d)
{

	pthread_mutex_init(&mute, NULL);

	struct vptree *vpt = (struct vptree*) malloc(sizeof(struct vptree));
	shared_points = (double *) malloc(n*d*sizeof(double));
	distances = (double *) malloc(n*sizeof(double));
	shared_points = X;


	vpt = calculate_median(vpt   ,0 , n-1, d  );
	

	pthread_mutex_destroy(&mute);

	return vpt;
}

struct vptree * getInner(struct vptree * T)
{
	return T->smaller;
}

struct vptree * getOuter(struct vptree * T)
{
	return T->bigger;
}

double getMD(struct vptree * T)
{
	return T->median_value;
}

double * getVP(struct vptree * T)
{
	return T->vantage_point;
}

int getIDX(struct vptree * T)
{
	return T->idx;
}




void *parallel_calculate_median(void *arg)
{
	struct data* args = (struct data *) arg;
	struct data *new_args1;
	//data *new_args2;
	pthread_t thread0;
	//pthread_t thread1;
	int start = args->start;
	int end = args->end;
	int d = args->d;
	int n=end-start+1;
	
	pthread_attr_t attr;       
    /* Initialize and set thread detached attribute */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Initiate node with the selecetd(first point) vantage point
	args->node = newnode(shared_points , args->end , 0 , args->d);
	
	if(end-start==0)
	{
	}
	//In this occasion we have two points in the array. One of them will
	//be chosen as a vantage point
	else if(end-start==1)
	{

		calculate_distances( args->node , start , end-1 , d);
		double median = distances[start];
		args->node->median_value = median;
		args->node->smaller = calculate_median(args->node->smaller   , start ,  end-1 ,   d );
	}
	//In this a occasion we are having three points and it is going to be created a pefect binary tree node 
	else if(end-start==2)
	{
		calculate_distances(args->node , start , end-1 , d);
			
		if(distances[start]>distances[end-1])
		{
			double median = distances[end-1];
			args->node->median_value = median;

			args->node->smaller = calculate_median(args->node->smaller   , start+1 ,  end-1 ,   d );
			args->node->bigger = calculate_median(args->node->bigger   , start ,  end-2 ,   d );

		}
		else
		{
			double median = distances[start];
			args->node->median_value = median;

			new_args1 = Initiate_data(args->node->smaller,start,end-2,d);
			
			args->node->smaller = calculate_median(args->node->smaller   , start ,  end-2 ,   d );
			args->node->bigger = calculate_median(args->node->bigger   , start+1,  end-1 ,   d );

		}
	
	}
		//General occasion
	else
	{
			calculate_distances(args->node , start , end-1 , d);
				
			double median = quickselect( distances , shared_points, start , end-1 , start+(n-1)/2 , d );
			args->node->median_value = median;

			new_args1 = Initiate_data(args->node->smaller,start, start+(n-1)/2 ,d);

		if(threshold>1){
			
			pthread_mutex_lock(&mute);
			threshold-=1;							//Locking code for counter increase 
			pthread_mutex_unlock(&mute);

			pthread_create( &thread0, &attr, parallel_calculate_median,(void *) new_args1);
			args->node->bigger = calculate_median(args->node->bigger   , start+ (n-1)/2 +1 ,  end-1 ,   d );
				
			pthread_join(thread0,(void **) &args->node->smaller );

			pthread_mutex_lock(&mute);
			threshold+=1;					//Locking code for counter decrease
			pthread_mutex_unlock(&mute);
		}
		else
		{
			args->node->smaller = calculate_median(args->node->smaller   ,start, start+(n-1)/2,   d );
			args->node->bigger = calculate_median(args->node->bigger   , start+(n-1)/2+1 ,  end-1 ,   d );
		}
	}
	pthread_exit((void *) args->node);
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
		//Calculate the distance
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
	
		//calculate the distances
		calculate_distances(node , start , end-1 , d);
				
		double median = quickselect( distances , shared_points, start , end-1 , start+(n-1)/2 , d );
		node->median_value = median;

		if(threshold>1){
			
			pthread_t thread0;
			pthread_attr_t attr;       
   			/* Initialize and set thread detached attribute */
   			pthread_attr_init(&attr);
   			
   			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			struct data *new_args1;

			new_args1 = Initiate_data(node->smaller,start, start+(n-1)/2 ,d);
			new_args1 = Initiate_data(node->smaller,start, start+(n-1)/2 ,d);
			pthread_mutex_lock(&mute);
			threshold-=1;							//Locking code for counter increase 
			pthread_mutex_unlock(&mute);
			pthread_create( &thread0, &attr, parallel_calculate_median,(void *) new_args1);
	
			node->bigger = calculate_median(node->bigger   , start+(n-1)/2+1 ,  end-1 ,   d );
				
				//printf("finish and waiting\n");
			pthread_join(thread0,(void **) &node->smaller );
		
			pthread_mutex_lock(&mute);
			threshold+=1;					//Locking code for counter decrease
			pthread_mutex_unlock(&mute);
			return node;
		}
		else
		{
			node->smaller = calculate_median(node->smaller   ,start, start+(n-1)/2,   d );
			node->bigger = calculate_median(node->bigger   , start+(n-1)/2+1 ,  end-1 ,   d );
			return node;
		}
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

