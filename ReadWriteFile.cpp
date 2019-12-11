#include "StructsAndFuncsDecl.h"

// reads input file and loads data to fileProperties and pointsData
PointsData readInput(FileProperties* consts)
{
	const char* INPUT_PATH = "D:\\Rotem\\data1.txt";
	PointsData pointsData;
	FILE* f = fopen(INPUT_PATH, "r");

	if (f == NULL)
	{
		printf("error opening input file, make sure the .txt file is in the correct path.\n");
		exit(1);
	}

	// constants
	fscanf(f, "%d %d %lf %lf %d %lf\n", &consts->numberOfPoints, &consts->dimension, &consts->alpha0, &consts->alphaMax, &consts->limit, &consts->qc);
	// points coordinates array
	pointsData.coordinations = (double*)calloc(consts->numberOfPoints * consts->dimension, sizeof(double));
	//points sets array
	pointsData.tags = (int*)calloc(consts->numberOfPoints, sizeof(int));
	// K
	pointsData.dimension = consts->dimension;
	// N
	pointsData.numberOfPoints = consts->numberOfPoints;

	// read data to array
	for (int i = 0; i < consts->numberOfPoints; ++i)
	{
		for (int j = 0; j < consts->dimension; ++j)
		{
			fscanf(f, "%lf ", &pointsData.coordinations[(consts->dimension * i) + j]);
		}
		fscanf(f, "%d", &pointsData.tags[i]);
	}
	fclose(f);

	return pointsData;
}

// receives final results, writes to output file and prints to console
void writeOutput(double alpha, double alphaMax, double q, double* W, double w0, int dimension)
{
	const char* OUTPUT_PATH = "D:\\Rotem\\output.txt";
	FILE* f = fopen(OUTPUT_PATH, "w");

	if (f == NULL)
	{
		printf("error opening output file, make sure the .txt file is in the correct path.\n");
		exit(1);
	}
	if (alpha >= alphaMax) // failed to find minAlpha 
	{
		fprintf(f, "Alpha is not found.\n");
		printf("Alpha is not found.\n");
	}
	else //print found values 
	{
		
		fprintf(f, "Alpha minimum: %lf, q: %lf\n", alpha, q);
		printf("Alpha minimum: %lf, q: %lf\n", alpha, q);

		fprintf(f, "w0 = %lf\n", w0);
		printf("w0 = %lf\n", w0);

		for (int i = 0; i < dimension; ++i)
		{
			fprintf(f, "W%d = %lf\n", i + 1, W[i]);
			printf("W%d = %lf\n", i + 1, W[i]);
		}
	}


	fclose(f);
}