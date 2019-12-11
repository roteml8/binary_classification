#include "StructsAndFuncsDecl.h"

// redefines the weights vector W according to the algorithms formula
void redefineWeights(double* W, double* w0, double alpha, int dimension, int sign, const double * point)
{
#pragma omp parallel for
	for (int i = 0; i < dimension; ++i)
	{
		W[i] += (alpha * sign) * point[i];
	}
	*w0 += (alpha * sign);
}

// calculates the discriminant function value of given point
double discFunction(const double * point,const double* W, double w0 ,int dimension)
{
	int i;
	double result = 0;
	for (i = 0; i < dimension; ++i)
	{
		result += W[i] * point[i];
	}
	result += w0;
	return result;
}


// gathers results from processes, checks for success and sends to master
int gatherAndCheck(double* myResults, int numprocs, const FileProperties* consts, int id, double* w0, double* W, int* minAlphaId,
	double* results)
{
	int i, dimensions = consts->dimension, finished = 0;;
	double* allResults = (double*)calloc(numprocs * 2, sizeof(double)); // alphas and q's from all processes

	MPI_Allgather(myResults, 2, MPI_DOUBLE, allResults, 2, MPI_DOUBLE, MPI_COMM_WORLD); 

	*minAlphaId = NONE; // initiate rank of process with minAlpha to nonexistent

	for (i = 0; i < numprocs; ++i) // check each process results
	{
		if (allResults[i * 2] < consts->qc || allResults[i * 2 + 1] > consts->alphaMax - consts->alpha0) // found required quality or reached the end
		{
			*minAlphaId = i;
			break;
		}
	}

	if (MASTER != *minAlphaId && *minAlphaId != NONE) {
		if (id == *minAlphaId) // rank with minAlpha sends to master
		{
			MPI_Send(w0, 1, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
			MPI_Send(W, dimensions, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
		}
		else if (id == MASTER) // master receives from rank with minAlpha
		{
			MPI_Recv(w0, 1, MPI_DOUBLE, *minAlphaId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(W, dimensions, MPI_DOUBLE, *minAlphaId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	results[0] = allResults[2 * (*minAlphaId)]; // q of process with minAlpha
	results[1] = allResults[2 * (*minAlphaId) + 1]; // alpha of process with minAlpha
	free(allResults);
	return *minAlphaId != NONE;
}

// get specific point by location 
const double *getPointAtLocation(const PointsData* points, int location)
{
	double *coords = points->coordinations;
	int dim = points->dimension;
	return  coords + (dim * location);
}

// receives value and returns its sign
int getSign(double value)
{
	if (value > 0)
		return 1;
	return -1;
}


// perform classification of the points and check accuracy
//returns 1 if classification failed, 0 else
int classification(int numOfPoints, const PointsData* points, double* W, double* w0, double alpha, double *classified)
{
	int i, k = points->dimension;
#pragma omp parallel for
	for (i = 0; i < numOfPoints; ++i)
	{
		const double* point = getPointAtLocation(points, i); // per each point
		classified[i] = discFunction(point, W, *w0, k); // calculate value of f

	}
	for (i = 0; i < numOfPoints; i++)
	{
		int sign = getSign(classified[i]); // get sign of result
		const double *point = getPointAtLocation(points, i);
		if (sign != points->tags[i]) // misclassified
		{
			redefineWeights(W, w0, alpha,k, points->tags[i], point); // update W
			return 1;
		}
	}
	return 0;
}

// calculate quality of classifier according to the classifications received and points actual sets
double getQuality(const PointsData *points, const double *classified)
{
	double count = 0; // counter for misclassified points
	int i;
#pragma omp parallel for reduction(+:count)
	for (i = 0; i < points->numberOfPoints; i++)
	{
		int sign = getSign(classified[i]);
		if (points->tags[i] != sign) // mismatch
			count++;
	}
	return count / points->numberOfPoints;
}

// main algorithm of the binary classifier
void mainAlgorithm(const FileProperties* consts, const PointsData* points, int numprocs, int id)
{
	int classifyFail; // 1 if classified failed, 0 else
	int dimension = consts->dimension; // K - num of coordinates
	int minAlphaId = NONE; // set rank of process with minAlpha to nonexistent
	double alpha = consts->alpha0 * (1 + id); // unique value of alpha based on proccess rank
	double quality; // quality of classifier
	double w0 = 0; // w0 sets to 0
	double *weights = (double*)calloc(dimension, sizeof(double)); // weights vector sets to zeroes
	double currentResults[2]; // current q and alpha
	double results[2]; // final results of q and alpha
	double *classified = (double*)calloc(consts->numberOfPoints, sizeof(double)); // results of current classification

	if (id == MASTER)
	{
		printf("\nqc = %lf alphaMax = %lf\n", consts->qc, consts->alphaMax);
		fflush(stdout);
	}

	do {
		int iteration = 0; // start counting iterations
		do {
			classifyFail = classification(points->numberOfPoints, points, weights, &w0, alpha, classified); // get current classification
			iteration++;
		} while (classifyFail && iteration < consts->limit); // classification fails and limit is not reached

		quality = getQuality(points, classified); // get quality of classification
		//set current results
		currentResults[0] = quality;
		currentResults[1] = alpha;
		alpha += consts->alpha0 * numprocs;

		gatherAndCheck(currentResults, numprocs, consts, id,
			&w0, weights, &minAlphaId, results); // synchronize and check results

	} while (quality > consts->qc &&
		currentResults[1] < consts->alphaMax
		&& minAlphaId == NONE); // q is too big & current alpha<alphaMax & rank of minAlpha is none

	alpha -= consts->alpha0 * numprocs;

	if (id == MASTER)
	{
		double minAlpha; // min value of alpha
		if (minAlphaId == NONE) // no alpha as required
			minAlpha = alpha;
		else // found a process with an alpha as required
		{
			minAlpha = results[1];
			quality = results[0];
		}
		fflush(stdout);
		writeOutput(minAlpha, consts->alphaMax, quality, weights, w0, dimension); // send results to the function that writes to output file
	}

	free(weights);
	free(classified);
}

