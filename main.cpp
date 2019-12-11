#include "mpi.h"
#include "StructsAndFuncsDecl.h"
// Rotem Levi
int main(int argc, char* argv[])
{
	double t0; // start time
	double t1; // end time

	int id; // rank
	int numOfProcesses; // number of processes 

	MPI_Datatype filePropertiesType; // new MPI_Type for fileProperties struct
	int blocklen[] = { 1, 1, 1, 1, 1, 1, 1 };
	MPI_Aint disp[6];
	MPI_Datatype dataType[] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_DOUBLE };
	disp[0] = offsetof(FileProperties, numberOfPoints);
	disp[1] = offsetof(FileProperties, dimension);
	disp[2] = offsetof(FileProperties, alpha0);
	disp[3] = offsetof(FileProperties, alphaMax);
	disp[4] = offsetof(FileProperties, limit);
	disp[5] = offsetof(FileProperties, qc);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	MPI_Type_create_struct(6, blocklen, disp, dataType, &filePropertiesType);
	MPI_Type_commit(&filePropertiesType);

	FileProperties consts; // constants from file
	PointsData points; // points data from file

	if (id == MASTER) // read from file
	{
		points = readInput(&consts);
	}

	MPI_Bcast(&consts, 1, filePropertiesType, 0, MPI_COMM_WORLD); // Broadcast constant values

	if (id != MASTER)
	{
		points.numberOfPoints = consts.numberOfPoints;
		points.dimension = consts.dimension;

		points.coordinations = (double*)calloc(consts.numberOfPoints * consts.dimension, sizeof(double));
		points.tags = (int*)calloc(consts.numberOfPoints, sizeof(int));
	}

	MPI_Bcast(points.coordinations, consts.numberOfPoints * consts.dimension, MPI_DOUBLE, 0, MPI_COMM_WORLD);// Broadcast points coordinates
	MPI_Bcast(points.tags, consts.numberOfPoints, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast points sets

	t0 = MPI_Wtime();
	{
		mainAlgorithm(&consts, &points, numOfProcesses, id);
	}
	t1 = MPI_Wtime();

	if (id == MASTER) // print total runtime
	{
		printf("Time is: %lf\n", t1 - t0);
		fflush(stdout);
	}

	MPI_Finalize();
	return 0;
}


