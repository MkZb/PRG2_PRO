#include "Data.h"

void createAndFillMatrix(int** MA, int val)
{
	for (size_t i = 0; i < N; i++)
	{
		MA[i] = new int[N];
		for (size_t j = 0; j < N; j++)
		{
			MA[i][j] = val;
		}
	}
}

void fillVector(int* A, int val)
{
	for (size_t i = 0; i < N; i++)
	{
		A[i] = val;
	}
}

void fillMatrix(int** MA, int val)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			MA[i][j] = val;
		}
	}
}

void matrixToArray(int** MA, int* arr, int col, int row)
{
	int idx = 0;
	for (size_t i = 0; i < row; i++)
	{
		for (size_t j = 0; j < col; j++)
		{
			arr[idx] = MA[i][j]; 
			idx++;
		}
	}
}

int** arrayToMatrix(int* A)
{
	int** MR = new int* [N];
	for (size_t i = 0; i < N; i++)
	{
		MR[i] = new int[N];
		for (size_t j = 0; j < N; j++)
		{
			MR[i][j] = A[i * N + j];
		}
	}
	return MR;
}

int** partitionedArrToMatr(int* A)
{
	int** MR = new int* [N];
	for (size_t i = 0; i < N; i++)
	{
		MR[i] = new int[N];
	}
	int idx = 0;
	for (size_t i = 0; i < PNum; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			for (size_t k = 0; k < H; k++)
			{
				MR[j][i * H + k] = A[idx];
				idx++;
			}
		}
	}
	return MR;
}

void copyArray(int* src, int srcLen, int* dest)
{
	for (size_t i = 0; i < srcLen; i++)
	{
		dest[i] = src[i];
	}
}

void copyMatrixColumnsIntoArr(int** src, int colNum, int* dest)
{
	for (size_t i = 0; i < colNum; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			dest[i * N + j] = src[j][i];
		}
	}
}

int** arrayColumnsIntoMatrix(int* src, int colNum)
{
	int** MR = new int* [N];
	for (size_t i = 0; i < N; i++)
	{
		MR[i] = new int[colNum];
	}

	for (size_t i = 0; i < colNum; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			MR[j][i] = src[i * N + j];
		}
	}

	return MR;
}

int** calcMUH(int** MD, int** MCH, int d, int a, int** MRH, int rank)
{
	int** MUH;
	MUH = sumMatrix(multMatrixByVal(multMatrByMatrPart(MD, MCH, rank), d, H, N), multMatrixByVal(MRH, a, H, N), H, N);
	return MUH;
}

int** multMatrByMatrPart(int** fullMA, int** partMB, int rank)
{
	int** MR = new int* [N];
	for (size_t i = 0; i < N; i++)
	{
		MR[i] = new int[H];
	}

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = rank * H; j < (rank + 1) * H; j++)
		{
			MR[i][j - rank * H] = 0;
			for (size_t k = 0; k < N; k++)
			{
				MR[i][j - rank * H] += fullMA[i][k] * partMB[k][j - rank * H];
			}
			
		}
	}
	
	return MR;
}

int** multMatrixByVal(int** MA, int val, int col, int row)
{
	int** MR = new int* [row];
	for (size_t i = 0; i < row; i++)
	{
		MR[i] = new int[col];
		for (size_t j = 0; j < col; j++)
		{
			MR[i][j] = MA[i][j] * val;
		}
	}
	return MR;
}

int** sumMatrix(int** MA, int** MB, int col, int row)
{
	int** MR = new int* [row];
	for (size_t i = 0; i < row; i++)
	{
		MR[i] = new int[col];
		for (size_t j = 0; j < col; j++)
		{
			MR[i][j] = MA[i][j] + MB[i][j];
		}
	}
	return MR;
}

