#define _CRT_SECURE_NO_WARNINGS
#include <mpi.h>
#include <cmath>
#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#define NONE -1 // represents nonexisting rank
#define MASTER 0 // represents master rank

// holds the constant values given in the input file 
struct FileProperties
{
	int numberOfPoints; // N - number of points
	int dimension; // K - number of coordinates of each point
	double alpha0; // increment value of alpha
	double alphaMax; // maximal value of alpha
	int limit; // max amount of iterations
	double qc; // quality of classifier to be reached
}
typedef FileProperties;

// holds the data of the points given in the input file
struct PointsData
{
	int numberOfPoints; //N - number of points
	int dimension;// K - number of coordinates of each point
	double* coordinations; // corrdinates of all points
	int* tags; // sets of points: 1 = A, -1 = B
}
typedef PointsData;

// all functions used in the program

PointsData readInput(FileProperties* consts);
void writeOutput(double alpha, double alphaMax, double q, double* W, double w0, int dimension);
void mainAlgorithm(const FileProperties* consts, const PointsData* points, int numprocs, int id);
void redefineWeights(double* W, double* w0, double alpha, int dimension, int sign, const double * point);
double discFunction(const double * point, const double* W, double w0, int dimension);
int gatherAndCheck(double* myResults, int numprocs, const FileProperties* consts, int id, double* w0, double* W, int* minAlphaId,
	double* results);
const double *getPointAtLocation(const PointsData* points, int location);
int getSign(double value);
int classification(int numOfPoints, const PointsData* points, double* W, double* w0, double alpha, double *classified);
double getQuality(const PointsData *points, const double *classified);
