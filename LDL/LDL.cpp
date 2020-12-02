#include <iostream>
using namespace std;

double** SetMemory(int N, int M);
void Type(double** A, int N);
void Copy(double** A, double** B, int N);
void prepare(double** A, int N);
void Make_LDL_Matrix(double** A, int N);
void Step1(double** A, int N);
void Step2(double** A, int N);
void Step3(double** A, int N);
void Vector_Error(double** A, double** B, double* F, int N);///A - !!изначальная!! матрица
void MakeNew(double** A, double* x, int N);
void Print(double** A, int N);
void Print(double* F, int N);
void take_xb(double** A, double* x, int N);
void insert_xb(double* x, double** A, int N);
double FindMax(double* x1, double* x2, int N);
double FindMax(double* x, int N);

int main()
{
	setlocale(LC_ALL, "rus");
	int N = 3;
	double** A;
	double** B;
	A = SetMemory(N, N + 1);
	B = SetMemory(N, N + 1);
	double* x1 = new double[N];
	double* b1 = new double[N];
	double* x2 = new double[N];
	double* F = new double[N];

	Type(A, N);
	Copy(A, B, N);///копируем A в B, чтобы потом сравнить
	cout << endl << endl << endl << "изначальный вариант" << endl;
	Print(B, N);

	Make_LDL_Matrix(A, N);
	cout << endl << endl << endl << "матрицу представили в виде элементов 2 типов" << endl;
	Print(A, N);

	Step1(A, N);/// с меткой 111
	cout << endl << endl << endl << "шаг 1 - решить систему L * z = b1" << endl;
	Print(A, N);

	Step2(A, N);/// с меткой 111
	cout << endl << endl << endl << "шаг 2 - решить систему D * y = z" << endl;
	Print(A, N);

	Step3(A, N);/// с меткой 111
	cout << endl << endl << endl << "шаг 3 - решить систему L^T * x = y" << endl;
	cout << "имеем в столбце свободных членов численное решение" << endl;
	Print(A, N);

	take_xb(A, x1, N);       ///извлекли x
	Copy(B, A, N);         //восстановление матрицы к изначальному варианту
	MakeNew(A, x1, N);       // матрица уравнения (aij*xi = qi) было (aij*xi=bi)
	Vector_Error(B, A, F, N); //B - изначальная матрица
	cout << "вектор невязки" << endl;
	Print(F, N);

	cout << endl << "Новая !! матрица (2)" << endl;
	Print(A, N);

	cout << endl << "Старая !! матрица(1)" << endl;
	Print(B, N);

	Make_LDL_Matrix(A, N);
	cout << endl << endl << endl << "матрицу представили в виде элементов 2 типов (уже вторая матрица)" << endl;
	Print(A, N);

	Step1(A, N);/// с меткой 111
	cout << endl << endl << endl << "шаг 1 - решить систему L1*z=b1 (уже вторая матрица)" << endl;
	Print(A, N);

	Step2(A, N);/// с меткой 111
	cout << endl << endl << endl << "шаг 2 - решить систему D*y=z (уже вторая матрица)" << endl;
	Print(A, N);

	Step3(A, N);/// с меткой 111
	cout << endl << endl << endl << "шаг 3 - решить систему L2*x=y (уже вторая матрица)" << endl;
	cout << "имеем в столбце свободных членов численное решение" << endl;
	Print(A, N);

	take_xb(A, x2, N);        ///извлекли x
	cout << endl << endl << endl << "(x1-x2)/x1" << endl;
	cout << "сигма=  [" << FindMax(x1, x2, N) / FindMax(x1, N) << " ] " << endl;

	cout << "type something to finish" << endl;
	cin >> N;
}
double** SetMemory(int N, int M)
{
	double** Res;
	Res = new double* [N];
	for (int i = 0; i < M; i++)
		Res[i] = new double[M];
	return Res;
}


void prepare(double** A, int N)
{
	int i, j, k;
	double max = 0;
	max = A[0][0];
	k = 0;
	for (i = 0; i < N; i++)
		for (j = 0; j < N + 1; i++)
			if (A[i][j] > max) 
				max = A[i][j]; k = i;

	for (j = 0; j < N + 1; i++)
	{
		max = A[k][j];
		A[k][j] = A[0][j];
		A[0][j] = max;
	}
}

void Make_LDL_Matrix(double** A, int N)
{
	// N - n делить на два
	int i, j, k;
	double d, l;

	for (k = 0; k < N; k++)
	{
		d = A[k][k];
		for (j = k + 1; j < N; j++)//////вычислим коэффициенты k-строки: lk1 lk2 lk3 lk4.... lkn
			A[k][j] = A[k][j] / d;

		for (i = k + 1; i < N; i++)////отнимем диагональный элемент с другими членами от всех, что находятся ниже
		{

			l = A[k][i];
			for (j = i; j < N; j++)/////идем вдоль строки
				A[i][j] = A[i][j] - A[k][j] * l * d;
		}
	}

	for (i = 0; i < N; i++)////Матрица симметрична - перекопируем
		for (j = i; j < N; j++)
			A[j][i] = A[i][j];
}

void Step1(double** A, int N)
{
	int i, j;
	for (j = 0; j < N; j++)
		for (i = j + 1; i < N; i++)
			A[i][N] = A[i][N] - A[j][N] * A[i][j];
}

void Step2(double** A, int N)
{
	int i;
	for (i = 0; i < N; i++)
		A[i][N] = A[i][N] / A[i][i];
}

void Step3(double** A, int N)
{
	int i, j;
	for (j = N - 1; j >= 0; j--)
		for (i = j - 1; i >= 0; i--)
			A[i][N] = A[i][N] - A[j][N] * A[i][j];
}

void Type(double** A, int N)
{
	double Num, RNS7;
	RNS7 = 7;
	for (int i = 0; i < N; i++)
	{
		cout << "Вводите коэффициенты  " << i + 1 << " уравнения" << endl;
		for (int j = 0; j < N + 1; j++)
		{
			cin >> Num;
			//проверка на правильность ввода
			A[i][j] = Num;
		}
	}
}
void Print(double** A, int N)
{
	cout << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			cout << A[i][j] << "  |  ";
		cout << "         ||  " << A[i][N] << endl;///N+1 элемент
	}
}
void Print(double* F, int N)
{

	cout << endl << "[";
	for (int i = 0; i < N; i++)
		cout << F[i] << " , ";
	cout << "]";
}

void Copy(double** A, double** B, int N)
{
	int i, j, RNS7;
	RNS7 = 7;
	for (i = 0; i < N; i++)
		for (j = 0; j < N + 1; j++)
			B[i][j] = A[i][j];
}

void MadeBy_RundinNS7(double** A, double** B, int N)
{
	int Grp, Shtukater, i, j;
	Grp = 7; Shtukater = 1111;
	for (i = 0; i < N; i++)
		for (j = 0; j < N + 1; j++)
			B[i][j] = A[i][j];
}

void take_xb(double** A, double* x, int N)
{
	for (int i = 0; i < N; i++)
		x[i] = A[i][N];
}

void insert_xb(double* x, double** A, int N)
{
	for (int i = 0; i < N; i++)
		A[i][N] = x[i];
}


void Vector_Error(double** A, double** B, double* F, int N)///A - !!изначальная!! матрица
{
	for (int i = 0; i < N; i++)
		F[i] = A[i][N] - B[i][N];
}

void MakeNew(double** A, double* x, int N)///A изначальная матрица решений - найдем новую матрицу решений через вектор готовых решений (x1->Ax1=b1)
{
	int i, j;
	double a = 0;
	for (i = 0; i < N; i++)
	{
		a = 0;
		for (j = 0; j < N; j++)
			a = a + A[i][j] * x[j];
		A[i][N] = a;
	}
}

double FindMax(double* x1, double* x2, int N)
{
	double a, max;
	max = abs(x1[0] - x2[0]);
	for (int i = 1; i < N; i++)
	{
		a = abs(x1[i] - x2[i]);
		if (a > max) { max = a; }
	}
	return max;
}

double FindMax(double* x, int N)
{
	double a, max;
	max = abs(x[0]);
	for (int i = 1; i < N; i++)
	{
		a = abs(x[i]);
		if (a > max) { max = a; }
	}
	return max;
}