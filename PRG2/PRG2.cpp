#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include "mpi.h"
#include "Data.h"

using namespace std::chrono;

int main(int *argc, char **argv) {
	int rank, numtasks;
	high_resolution_clock::time_point start_time = high_resolution_clock::now();
	MPI_Init(argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	if (numtasks != PNum) {
		printf("Please set correct number of processors so it is the same with command line argument!\n");
		exit(1);
	}

	//Creating graph
	MPI_Comm graph_comm;
	int nnodes = PNum;
	int* index = new int[PNum];
	int* edges = new int[6 * P];

	index[0] = P;
	int accum = P;
	for (size_t i = 1; i < PNum; i++)
	{
		if (i % 3 == 1) {
			accum += 3;
			index[i] = accum;
		}
		else {
			accum += 1;
			index[i] = accum;
		}
	}
	
	for (size_t i = 0; i < P; i++)
	{
		edges[i] = i + 1;
	}

	int idx = 0;
	for (size_t i = 0; i < P; i++)
	{
		edges[P + idx] = 0;
		idx++;
		edges[P + idx] = 3 * i + 2;
		idx++;
		edges[P + idx] = 3 * i + 3;
		idx++;
		edges[P + idx] = 3 * i + 1;
		idx++;
		edges[P + idx] = 3 * i + 1;
		idx++;
	}

	MPI_Graph_create(MPI_COMM_WORLD, nnodes, index, edges, false, &graph_comm);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("T%d started\n", rank);
	
	int* buff = new int[N * N];
	int* rcvBuff = new int[N * N];
	int d, a;
	int** MD,** MCH,** MRH,** MUH;
	int* ZH = new int[H];

	//Input
	if (rank == 0)
	{
		int d_input = 1;
		buff[0] = d_input;
	}
	
	MPI_Bcast(&(buff[0]), 1, MPI_INT, 0, graph_comm);
	d = buff[0];

	if (rank == 0) {
		int** MD_input = new int* [N];
		createAndFillMatrix(MD_input, 1);
		matrixToArray(MD_input, buff, N, N);
	}

	MPI_Bcast(&(buff[0]), N * N, MPI_INT, 0, graph_comm);
	MD = arrayToMatrix(buff);	

	if (rank == 0) {
		int* Z = new int[N];
		fillVector(Z, 1);
		copyArray(Z, N, buff);
	}

	MPI_Scatter(&(buff[0]), H, MPI_INT, &(ZH[0]), H, MPI_INT, 0, graph_comm);
	
	if (rank == 0) {
		int** MC = new int* [N];
		createAndFillMatrix(MC, 1);
		copyMatrixColumnsIntoArr(MC, N, buff);
	}
	
	MPI_Scatter(&(buff[0]), H * N, MPI_INT, &(rcvBuff[0]), H * N, MPI_INT, 0, graph_comm);
	MCH = arrayColumnsIntoMatrix(rcvBuff, H);

	if (rank == 0) {
		int** MR = new int* [N];
		createAndFillMatrix(MR, 1);
		copyMatrixColumnsIntoArr(MR, N, buff);
	}
	MPI_Scatter(&(buff[0]), H * N, MPI_INT, &(rcvBuff[0]), H * N, MPI_INT, 0, graph_comm);
	MRH = arrayColumnsIntoMatrix(rcvBuff, H);

	//Calculation of local max(ZH)
	int a_local = INT32_MIN;
	for (size_t i = 0; i < H; i++)
	{
		if (ZH[i] > a_local) {
			a_local = ZH[i];
		}
	}

	//Gathering local max at one thread
	buff[0] = a_local;
	MPI_Gather(&(buff[0]), 1, MPI_INT, &(rcvBuff[0]), 1, MPI_INT, 0, graph_comm);

	//Calculation of max(Z)
	if (rank == 0) {
		a = rcvBuff[0];
		for (size_t i = 1; i < PNum; i++)
		{
			if (rcvBuff[i] > a) {
				a = rcvBuff[i];
			}
		}
		buff[0] = a;
	}
	
	//Sending a = max(Z)
	MPI_Bcast(&(buff[0]), 1, MPI_INT, 0, graph_comm);
	a = buff[0];
	
	//Main calculation
	MUH = calcMUH(MD, MCH, d, a, MRH, rank);
	
	//Gathering results from all threads
	matrixToArray(MUH, buff, H, N);
	MPI_Gather(&(buff[0]), H * N, MPI_INT, &(rcvBuff[0]), N * H, MPI_INT, 0, graph_comm);

	//Printing results
	if (rank == 0) {
		int** MU;
		MU = partitionedArrToMatr(rcvBuff);
		/*
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				printf("%d ", MU[i][j]);
			}
			printf("\n");
		}
		*/
		high_resolution_clock::time_point finish_time = high_resolution_clock::now();
		duration<double, std::milli> time_span = finish_time - start_time;
		std::cout << "Time: " << time_span.count() << "\n";
	}
	printf("T%d finished\n", rank);
	MPI_Finalize();
}