#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

double t; //периодичность измерений
double noise = 0.2; //дисперсия датчика
const double ksi = 0.05; //дисперсия модели движения
const double randMove = 0.02;
const double r = 0.92;
const double l_kernel = 0.2;
const double r_kernel = 0.35;


double mod(double a) {
	if (a > 0)
		return a;
	return -a;
}

void initialization(double* A, double alpha, int M) {
	for (int i = 0; i < M; i++) {
		A[i] = alpha;
	}
}

double* m_mult_a(double* A, double alpha, int M, int flag) {
	double* B = (double*)malloc(sizeof(double) * M);
	for (int i = 0; i < M && flag == 1; i++) {
		B[i] = A[i] * alpha;
	}
	for (int i = 0; i < M && flag == -1; i++) {
		B[i] = A[i] / alpha;
	}
	return B;
}

double* m_add(double* A, double* B, int M, int N, int sign) {
	double* C = (double*)malloc(sizeof(double) * M * N);
	for (int i = 0; i < M * N; i++) {
		C[i] = A[i] + (B[i] * sign);
	}
	return C;
}

double* m_mult(double* A, double* B, int M, int K, int N) {
	double* C = (double*)malloc(sizeof(double) * N * K);
	initialization(C, 0, M * N);

	for (int z = 0; z < M; z++) {
		for (int j = 0; j < N; j++) {
			for (int i = 0; i < K; i++) {
				C[z * N + j] += (A[z * K + i] * B[i * N + j]);
			}
		}
	}
	return C;
}

double* m_pow(double* A, int alpha, int M) {
	double* B = (double*)malloc(sizeof(double) * M * M);
	B = A;
	for (int i = 1; i < alpha; i++) {
		B = m_mult(B, A, M, M, M);
	}
	return B;
}

double* findInvMatrix(double* matrix, int M) {
	double* result = (double*)malloc(sizeof(double) * M * M);
	double val = matrix[0];
	for (int i = 0, j = 0; i < M * M; i++) {
		result[i] = 0;
		if (i == j) {
			result[i] = 1 / val;
			j += M + 1;
		}
	}
	return result;
}

void randomFilling(double** array, int N, int M) {
	srand(time(NULL));
	int a;

	for (int i = 0; i < N * M; i++) {
		for (int j = 0; j < 2; j++) {
			int sign = 1, right = 1000;

			a = (rand() % right) + 1;
			if (rand() % 2001 < 1000)
				sign = -1;
			array[j][i] = a * sign;
		}
	}
}

double* findRandValue(double* R, double* start_matrix, double* rate, double** randValue, double** rand_matrix, int num, int M) {
	for (int i = 0; i < M; i++) {
		rand_matrix[0][i] = randValue[0][num + i] / 1000;
		rand_matrix[1][i] = randValue[1][num + i] / 1000;
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < 2; j++) {
			double a = rate[i];
			if (j == 1)
				a = noise;
			rand_matrix[j][i] *= a;
		}
	}

	double* new_matrix = (double*)malloc(sizeof(double) * M);
	new_matrix = m_add(m_mult(R, start_matrix, M, M, 1), rand_matrix[0], M, 1, 1);
	return new_matrix;
}

int check_is_not_NULL(FILE* f1, FILE* f2, FILE* f3, char flag) {
	if (f1 == NULL || f2 == NULL || f3 == NULL) {
		printf("%c-file is not founded!", flag);
		return 1;
	}
	return 0;
}

double* kinematicFunction(double* s_0, double* v_0, double* a_0, double M) {
	double* displacement = (double*)malloc(sizeof(double) * M);
	double* v_x = (double*)malloc(sizeof(double) * M);
	double* a_x = (double*)malloc(sizeof(double) * M);

	v_x = m_mult_a(v_0, t, M, 1);
	a_x = m_mult_a(a_0, t * t / 2, M, 1);
	displacement = m_add(a_x, v_x, M, 1, 1);
	displacement = m_add(displacement, s_0, M, 1, 1);

	return displacement;
}

int filter_kalmana(int N, int M) {
	FILE* finX, * finY, * finZ;
	FILE* foutX, * foutY, * foutZ;
	FILE* frealX, * frealY, * frealZ;

	char nameINX[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/inX.txt";
	char nameOUTX[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/outX.txt";
	char nameREALX[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/realX.txt";
	finX = fopen(nameINX, "w");
	foutX = fopen(nameOUTX, "w");
	frealX = fopen(nameREALX, "w");
	if (check_is_not_NULL(finX, foutX, frealX, 'X')) {
		return 1;
	}
	char nameINY[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/inY.txt";
	char nameOUTY[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/outY.txt";
	char nameREALY[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/realY.txt";
	finY = fopen(nameINY, "w");
	foutY = fopen(nameOUTY, "w");
	frealY = fopen(nameREALY, "w");
	if (check_is_not_NULL(finY, foutY, frealY, 'Y')) {
		return 1;
	}
	char nameINZ[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/inZ.txt";
	char nameOUTZ[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/outZ.txt";
	char nameREALZ[] = "C:/Users/Xorygal/source/repos/project/gnuplot/bin/realZ.txt";
	finZ = fopen(nameINZ, "w");
	foutZ = fopen(nameOUTZ, "w");
	frealZ = fopen(nameREALZ, "w");
	if (check_is_not_NULL(finZ, foutZ, frealZ, 'Z')) {
		return 1;
	}

	//---------------------------------------------------------//

	double* real_coordinates = (double*)malloc(sizeof(double) * M);
	initialization(real_coordinates, 0, M);

	double* evaluation = (double*)malloc(sizeof(double) * M);
	double* d_evaluation = (double*)malloc(sizeof(double) * M * M);
	double* signal = (double*)malloc(sizeof(double) * M);

	double** random_value = (double**)malloc(sizeof(double*) * 2);
	random_value[0] = (double*)malloc(sizeof(double) * N * M);
	random_value[1] = (double*)malloc(sizeof(double) * N * M);

	double* R = (double*)malloc(sizeof(double) * M * M);
	double* Vksi = (double*)malloc(sizeof(double) * M * M);
	double* V = (double*)malloc(sizeof(double) * M * M);
	for (int i = 0, j = 0; i < M * M; i++) {
		R[i] = 0;
		Vksi[i] = 0;
		V[i] = 0;
		if (i == j) {
			R[i] += r;
			Vksi[i] += (ksi);
			V[i] += noise;
			j += M + 1;
		}
	}

	double* V_inv = (double*)malloc(sizeof(double) * M * M);
	double* rand_move_matr = (double*)malloc(sizeof(double) * M);
	V_inv = findInvMatrix(V, M);
	initialization(rand_move_matr, randMove, M);

	randomFilling(random_value, N, M);

	//---------------------------------------------------------------//

	double* v = (double*)malloc(sizeof(double) * M);
	double* v_0 = (double*)malloc(sizeof(double) * M);
	double* a = (double*)malloc(sizeof(double) * M);
	double* s_0 = (double*)malloc(sizeof(double) * M);

	initialization(v, 0, M);
	initialization(a, 0, M);
	initialization(v_0, 0, M);
	initialization(s_0, 0, M);
	double sum_t = t;

	for (int i = 0; i < N * M; i += M) {
		double* x_ev = (double*)malloc(sizeof(double) * M);
		double* d_ev = (double*)malloc(sizeof(double) * M * M);
		double* temp = (double*)malloc(sizeof(double) * M * M);

		double** rand_matrix = (double**)malloc(sizeof(double*) * 2);
		rand_matrix[0] = (double*)malloc(sizeof(double) * M);
		rand_matrix[1] = (double*)malloc(sizeof(double) * M);

		if (i == 0) {
			real_coordinates = findRandValue(R, real_coordinates, rand_move_matr, random_value, rand_matrix, i, M);
			signal = m_add(real_coordinates, rand_matrix[1], M, 1, 1);

			evaluation = signal;
			d_evaluation = V;

		}
		else {
			real_coordinates = kinematicFunction(real_coordinates, v, a, M);

			initialization(temp, ksi, M);
			real_coordinates = findRandValue(R, real_coordinates, temp, random_value, rand_matrix, i, M);
			signal = m_add(real_coordinates, rand_matrix[1], M, 1, 1);

			x_ev = m_mult(R, evaluation, M, M, 1);

			temp = m_pow(R, 2, M);
			d_ev = m_mult(temp, d_evaluation, M, M, M);
			d_ev = m_add(d_ev, Vksi, M, M, 1);

			d_evaluation = m_mult(d_ev, V, M, M, M);
			temp = m_add(d_ev, V, M, M, 1);
			temp = findInvMatrix(temp, M);
			d_evaluation = m_mult(d_evaluation, temp, M, M, M);

			temp = m_mult(d_evaluation, V_inv, M, M, M);
			evaluation = m_add(signal, x_ev, M, 1, -1);
			evaluation = m_mult(temp, evaluation, M, M, 1);
			evaluation = m_add(x_ev, evaluation, M, 1, 1);


		}
		double* result = (double*)malloc(sizeof(double) * M);
		result = evaluation;
		result = m_add(result, m_mult_a(v, t, M, 1), M, 1, 1);
		result = m_add(result, m_mult_a(v, t, M, 1), M, 1, 1);

		double k_result = r + (noise / ksi) / 100;

		if (noise < l_kernel) {
			for (double i = noise; i < l_kernel; i += 0.07)
				k_result -= pow(mod(l_kernel - noise), 2);
		}
		else if (noise > r_kernel) {
			k_result += pow(mod(r_kernel - noise), 2);

		}

		double b = i * t / M;

		a = m_mult_a(real_coordinates, 2, M, 1);
		a = m_mult_a(a, sum_t * sum_t, M, -1);

		v_0 = v;
		v = m_mult_a(a, t, M, 1);
		v = m_add(v, v_0, M, 1, 1);
		s_0 = real_coordinates;

		printf("%lf %lf %lf\n", s_0[0], s_0[1], s_0[2]);

		fprintf(finX, "%lf %lf\n", b, signal[0]);
		fprintf(finY, "%lf %lf\n", b, signal[1]);
		fprintf(finZ, "%lf %lf\n", b, signal[2]);

		fprintf(foutX, "%lf %lf\n", b, k_result * result[0]);
		fprintf(foutY, "%lf %lf\n", b, k_result * result[1]);
		fprintf(foutZ, "%lf %lf\n", b, k_result * result[2]);

		fprintf(frealX, "%lf %lf\n", b, real_coordinates[0]);
		fprintf(frealY, "%lf %lf\n", b, real_coordinates[1]);
		fprintf(frealZ, "%lf %lf\n", b, real_coordinates[2]);

		sum_t += t;
	}

	return 0;
}

char* concat(const char* s1, const char* s2) {
	char* result = (char*)malloc(strlen(s1) + strlen(s2) + 1);
	strcpy(result, s1);
	strcat(result, s2);
	return result;
}

int main() {
	int N, M = 3;

	printf("Enter the count of measurement(N < 2000): N = ");
	scanf("%d", &N);
	printf("\nEnter pediod time of measurement(t < 3): t = ");
	scanf("%lf", &t);
	if (t >= 3 || N >= 2000) {
		return 1;
	}
	printf("\n");

	if (filter_kalmana(N, M)) {
		return 1;
	}

	char* path = (char*)malloc(sizeof(char) * 256);
	char* choose = (char*)malloc(sizeof(char) * 256);
	char a;
	char enter;

	printf("Please, enter the path to 'bin' in gnuplot: path = ");
	scanf("%s", path);
	scanf("%c", &enter);
	char* name = concat(path, "gnuplot -persist");

	FILE* GNUpipe = _popen(name, "w");
	if (GNUpipe == NULL) {
		printf("\nOops... The path was not found!\n");
		return 1;
	}

	printf("\nPlease, enter one of three coordinates (X, Y, Z): coordbnates = ");
	scanf("%c", &a);
	if (a != 'X' && a != 'Y' && a != 'Z') {
		printf("\nIncorrect input format\n");
		return 1;
	}
	scanf("%c", &enter);
	
	printf("What graph is not needed (in, out, real)? : ");
	scanf("%s", choose);


	if (!strcmp(choose, "real"))
		fprintf(GNUpipe, "plot '%sin%c.txt' with lines title 'in', '%sout%c.txt' with lines title 'out' \n", path, a, path, a);
	else if (!strcmp(choose, "in"))
		fprintf(GNUpipe, "plot '%sout%c.txt' with lines title 'out', '%sreal%c.txt' with lines title 'real' \n", path, a, path, a);
	else if (!strcmp(choose, "out"))
		fprintf(GNUpipe, "plot '%sin%c.txt' with lines title 'in', '%sreal%c.txt' with lines title 'real' \n", path, a, path, a);
	else 
		fprintf(GNUpipe, "plot '%sin%c.txt' with lines title 'in', '%sout%c.txt' with lines title 'out', '%sreal%c.txt' with lines title 'real' \n", path, a, path, a, path, a);
	fflush(GNUpipe);
	return 0;
}