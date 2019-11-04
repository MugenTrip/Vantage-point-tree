#include <stdio.h>

#define SWAP(x, y) { double temp = x; x = y; y = temp; }
#define N (sizeof(A)/sizeof(A[0]))

// Partition using Lomuto partition scheme
int partition(double a[] , double *X , int left, int right, int pivotIndex , int d)
{
	// Pick pivotIndex as pivot from the array
	double pivot = a[pivotIndex];

	// Move pivot to end
	SWAP(a[pivotIndex], a[right]);
	for (int j = 0; j < d; j++)
	{
		SWAP(X[(pivotIndex)*d+j] , X[(right)*d+j]);
	}

	// elements less than pivot will be pushed to the left of pIndex
	// elements more than pivot will be pushed to the right of pIndex
	// equal elements can go either way
	int pIndex = left;
	int i;

	// each time we finds an element less than or equal to pivot, pIndex
	// is incremented and that element would be placed before the pivot.
	for (i = left; i < right; i++)
	{
		if (a[i] < pivot)
		{
			SWAP(a[i], a[pIndex]);
			for (int j = 0; j < d; j++)
			{
				SWAP(X[(i)*d+j] , X[(pIndex)*d+j]);	
			}
			pIndex++;
		}
	}

	// Move pivot to its final place
	for (int j = 0; j < d; j++)
	{
		SWAP(X[(pIndex)*d+j] , X[(right)*d+j]);
	}
	SWAP(a[pIndex], a[right]);
		
	// return pIndex (index of pivot element)
	return pIndex;
}

// Returns the k-th smallest element of list within left..right
// (i.e. left <= k <= right). The search space within the array is
// changing for each round - but the list is still the same size.
// Thus, k does not need to be updated with each round.
double quickselect(double A[], double *X , int left, int right, int k ,int d )
{
	// If the array contains only one element, return that element
	if (left == right)
		return A[left];

	// select a pivotIndex between left and right
	int pivotIndex = left + rand() % (right - left + 1);

	pivotIndex = partition(A, X , left, right, pivotIndex , d );

	// The pivot is in its final sorted position
	if (k == pivotIndex)
		return A[k];

	// if k is less than the pivot index
	else if (k < pivotIndex)
		return quickselect(A,X , left, pivotIndex - 1, k , d );

	// if k is more than the pivot index
	else
		return quickselect(A,X , pivotIndex + 1, right, k , d);
}
