#pragma once

//Processors Number = 3 * P + 1
const int P = 4;
const int PNum = 3 * P + 1;

const int N = PNum * 185;
const int H = N / PNum;

void createAndFillMatrix(int** MA, int val);
void fillVector(int* A, int val);
void fillMatrix(int** MA, int val);
void matrixToArray(int** MA, int* arr, int col, int row);
int** arrayToMatrix(int* A);
int** partitionedArrToMatr(int* A);
void copyArray(int* src, int srcLen, int* dest);
void copyMatrixColumnsIntoArr(int** src, int colNum, int* dest);
int** arrayColumnsIntoMatrix(int* src, int colNum);
int** calcMUH(int** MD, int** MCH, int d, int a, int** MRH, int rank);
int** multMatrByMatrPart(int** fullMA, int** partMB, int rank);
int** multMatrixByVal(int** MA, int val, int col, int row);
int** sumMatrix(int** MA, int** MB, int col, int row);

