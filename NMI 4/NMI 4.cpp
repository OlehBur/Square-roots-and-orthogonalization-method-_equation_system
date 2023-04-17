#include <Windows.h>
#include <iostream>
#include <cstdio>

using namespace std;


void MulMatrix(float(&a)[4][4], float(&b)[4][4], float(&c)[4][4], int row, int col, int col1) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col1; j++) {
			c[i][j] = 0;
			for (int k = 0; k < col; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
}


constexpr int matrix_size = 4;

int main() {

	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	int a_matrix[matrix_size][matrix_size] = {
		{	20,	3,	5,	7	},
		{	3,	15,	7,	1	},
		{	5,	7,	17,	2	},
		{	7,	1,	2,	15	}
	},
		b_matrix[matrix_size] = {
		32,
		11,
		24,
		24
	}/*,
		buffA_matrix[matrix_size][matrix_size]*/;

	float lb_triangl_matrix[matrix_size][matrix_size],
		transp_lb_triangl[matrix_size][matrix_size],
		y_matrix[matrix_size], x_matrix[matrix_size],
		buffA_matrix[matrix_size][matrix_size];

	printf("\tМатриця:\n");
	for (int i = 0; i < matrix_size; i++){
		for (int j = 0; j < matrix_size; j++) {
			lb_triangl_matrix[i][j] = 0;
			buffA_matrix[i][j] = 0.0f;
		}

		printf("\t%d\t%d\t%d\t%d\t|\t%d\n",
			a_matrix[i][0], a_matrix[i][1],
			a_matrix[i][2], a_matrix[i][3],
			b_matrix[i]);
		}
		
		//невиродженість матриці / non-degeneracy  matrix


		//симетричність матриці / symmetry matrix 
	int cnt{ 0 };
	for (int i{ 0 }; i < matrix_size; i++)
		for (int j{ 0 }; j < matrix_size; j++)
			if (a_matrix[i][j] == a_matrix[j][i])
				cnt++;
	if (cnt == pow(matrix_size, 2))
		printf("\n\tМатриця є симетричною.\n\n");
	else {
		printf("\n\tМатриця не є симетричною, \n\tа отже не підходить для розв'язку даним методом.\n\n");
		return 0;
	}

		//нижньо-трикутна матриця / bottom triangl matrix
	lb_triangl_matrix[0][0] = sqrt(a_matrix[0][0]);
	for (int col = 1; col < matrix_size; col++)
		lb_triangl_matrix[0][col] = a_matrix[0][col] / lb_triangl_matrix[0][0];

	lb_triangl_matrix[1][1] = sqrt(a_matrix[1][1] - pow(lb_triangl_matrix[0][1], 2));
	for (int col = 2; col < matrix_size; col++)
		lb_triangl_matrix[1][col] = (a_matrix[1][col] - lb_triangl_matrix[0][1] * lb_triangl_matrix[0][col]) / a_matrix[1][1];

	lb_triangl_matrix[2][2] = sqrt(a_matrix[2][2] - pow(lb_triangl_matrix[0][2], 2) - pow(lb_triangl_matrix[1][2], 2));
	lb_triangl_matrix[2][3] = (a_matrix[2][3] - lb_triangl_matrix[0][2] * lb_triangl_matrix[0][3]
		- lb_triangl_matrix[1][2] * lb_triangl_matrix[1][3]) / lb_triangl_matrix[2][2];
	
	float Euk4{ 0.0f };
	for (int k = 0; k < matrix_size-1; k++)
		Euk4 += pow(lb_triangl_matrix[k][3], 2);

	lb_triangl_matrix[3][3] = sqrt(a_matrix[3][3] - Euk4);

	printf("\tНижньо-трикутна матриця (U):\n");
	for (int i = 0; i < matrix_size; i++)
		printf("	\t%g\t\t%g\t\t%g\t\t%g\n",
			lb_triangl_matrix[i][0], lb_triangl_matrix[i][1], 
			lb_triangl_matrix[i][2], lb_triangl_matrix[i][3]);


	//транспонована нижньо-трикутна матриця / transp bottom triangl matrix
	for (int i = 0; i < matrix_size; i++)	//transp to new matrix
		for (int j = 0; j < matrix_size; j++)
			transp_lb_triangl[j][i] = lb_triangl_matrix[i][j];

	printf("\n\tUt:\n");
	for (int i = 0; i < matrix_size; i++)
		printf("	\t%g\t\t%g\t\t%g\t\t%g\n",
			transp_lb_triangl[i][0], transp_lb_triangl[i][1],
			transp_lb_triangl[i][2], transp_lb_triangl[i][3]);
	
	//is Ut*U=A?
	printf("\nПеревірка на рівність (Ut*U=A), Ut*U:\n");

	MulMatrix(transp_lb_triangl, lb_triangl_matrix, buffA_matrix, 4, 4, 4);

	for (int i = 0; i < matrix_size; i++)
		printf("	\t%g\t%g\t%g\t%g\n",
			buffA_matrix[i][0], buffA_matrix[i][1],
			buffA_matrix[i][2], buffA_matrix[i][3]);


	//знаходження у значень / finding y`s
	y_matrix[0] = b_matrix[0] / lb_triangl_matrix[0][0];
	y_matrix[1] = (b_matrix[1] - lb_triangl_matrix[0][1] * y_matrix[0]) / lb_triangl_matrix[1][1];
	y_matrix[2] = (b_matrix[2] - lb_triangl_matrix[0][2] * y_matrix[0] - lb_triangl_matrix[1][2] 
		* y_matrix[1])/ lb_triangl_matrix[2][2];

	y_matrix[3] = (b_matrix[3] - lb_triangl_matrix[0][3] * y_matrix[0] - lb_triangl_matrix[1][3] 
		* y_matrix[1] - lb_triangl_matrix[2][3] * y_matrix[2]) / lb_triangl_matrix[3][3];


	//знаходження х значень / get X`s
	x_matrix[3] = y_matrix[3] / lb_triangl_matrix[3][3];
	x_matrix[2] = (y_matrix[2] - lb_triangl_matrix[2][3] * x_matrix[3]) / lb_triangl_matrix[2][2];
	x_matrix[1] = (y_matrix[1] - lb_triangl_matrix[1][2] * x_matrix[2] - lb_triangl_matrix[1][3]
		* x_matrix[3]) / lb_triangl_matrix[1][1];

	x_matrix[0] = (y_matrix[0] - lb_triangl_matrix[0][1] * x_matrix[1] - lb_triangl_matrix[0][2]
		* x_matrix[2] - lb_triangl_matrix[0][3] * x_matrix[3]) / lb_triangl_matrix[0][0];

	printf("\n\t%g\t%g\t%g\t%g\n\n", x_matrix[0], x_matrix[1], x_matrix[2], x_matrix[3]);

	system("pause");

}
//#include<iostream>
//using namespace std;
//int main()
//{
//    setlocale(LC_CTYPE, "ukr");
//    float A[10][10];
//    float U[10][10];
//    float b[10], x[10], y[10];
//    int n, k;
//    int i, j;
//    float temp;
//
//    std::cout << "введіть розмірність матриці" << endl;
//    cin >> n;
//label:
//    std::cout << "введіть елементи симетричної  матриці " << n << "x" << n << endl;
//    for (i = 0; i < n; i++)
//        for (j = 0; j < n; j++)
//            cin >> A[i][j];
//        
//    for (i = 0; i < n; i++)
//        for (j = 0; j < n; j++)
//            U[i][j] = 0;
//        
//    //перевірка на симетричність
//    for (i = 0; i < n; i++)
//        for (j = 0; j < n; j++)
//            if (A[i][j] != A[j][i])
//            {
//                std::cout << "матриця не симетрична" << endl;
//                goto label;
//            }
//
//    std::cout << "введіть елементи вектора b" << n << "x" << n << endl;
//    for (i = 0; i < n; i++)
//        cin >> b[i];
//
//    for (int i = 0; i < n; i++) {
//        temp = 0;
//        for (int k = 0; k < i; k++)
//            temp = temp + U[k][i] * U[k][i];
//
//        U[i][i] = sqrt(A[i][i] - temp);
//        for (j = i; j < n; j++) {
//            temp = 0;
//            for (k = 0; k < i; k++)
//                temp = temp + U[k][i] * U[k][j];
//
//            U[i][j] = (A[i][j] - temp) / U[i][i];
//        }
//    }
//
//    for (i = 0; i < n; i++) {
//        for (j = 0; j < n; j++) 
//            std::cout << U[i][j] << "\t ";
//
//        std::cout << endl;
//    }
//
//    for (i = 0; i < n; i++) {
//        temp = 0;
//        for (int k = 0; k < i; k++)
//            temp = temp + U[k][i] * y[k];
//
//        y[i] = (b[i] - temp) / U[i][i];
//    }
//
//    for (i = n - 1; i >= 0; i--) {
//        temp = 0;
//
//        for (int k = i + 1; k < n; k++)
//            temp = temp + U[i][k] * x[k];
//
//        x[i] = (y[i] - temp) / U[i][i];
//    }
//
//    for (i = 0; i < n; i++)
//        std::cout << "x" << i << "= " << x[i] << endl;
//
//    std::system("pause");
//}