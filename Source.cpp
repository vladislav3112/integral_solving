#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>
# define PI  3.14159265358979323846
using namespace std;
const double INTEGRAL_VALUE = 3.578861536040539915439859609644293194417;
const double EPS = 10e-6;
//используются: a = 0.1; b = 2.3; alpha = 0.2.
//LU + SLAE 
double sourse_function(double x) {
	return (2.5*cos(2 * x)*exp(2 * x / 3) + 4 * sin(3.5*x)*exp(-3 * x) + 3 * x);
}
void mult(vector <vector <double>> A, vector <vector <double>> B,
	vector <vector <double>> &R, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
}
void show(vector <vector <double>> A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << A[i][j] << "\t";
		}
		cout << endl;
	}
}
void LU(vector <vector <double>> A, vector <vector <double>> &L,
	vector <vector <double>> &U, int n, vector<int> &transp) // transp - vector перестановок
{
	U = A;
	vector<double> tmp(n);
	transp.clear();
	double max_elem;


	for (int i = 0; i < n; i++) {
		max_elem = U[i][i];
		for (int j = i; j < n; j++) {
			if (U[j][i] > max_elem) max_elem = U[j][i];
		}
		while (U[i][i] != max_elem) {
			transp.push_back(i);
			tmp = U[i];
			for (int idx = i; idx < n - 1; idx++)  //перестановка строки в конец
			{
				U[idx] = U[idx + 1];
			}
			U[n - 1] = tmp;
		}
	}
	while (U[0][0] * U[1][1] == U[0][1] * U[1][0]) {
		transp.push_back(0);
		tmp = U[0];
		for (int idx = 0; idx < n - 1; idx++)  //перестановка строки в конец
		{
			U[idx] = U[idx + 1];
		}
		U[n - 1] = tmp;
	}

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			if (abs(U[i][i]) > 0.01*EPS) {
				L[j][i] = U[j][i] / U[i][i];//Если U[i][i]=0,то U[j]- в конец.
			}
			else {
				transp.push_back(i);
				tmp = U[i];
				for (int idx = i; idx < n - 1; idx++)  //перестановка строки в конец
				{
					U[idx] = U[idx + 1];
				}
				U[n - 1] = tmp;
			}
		}
	}
	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			for (int j = i; j < n; j++) {
				L[j][i] = U[j][i] / U[i][i];
			}
		for (int i = k; i < n; i++)
			for (int j = k - 1; j < n; j++) {
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
			}
	}
}

void slae_solution(vector <vector <double>> L,
	vector <vector <double>> U, vector<double>b, vector<double>&x, vector<int> transp, int n ) {
	double sum = 0;
	vector <double> y(n);
	//������������ ������� b � ��������� L � U:
	double tmp;
	if (!transp.empty()) {
		for (int j = 0; j < transp.size(); j++) {
			tmp = b[transp[j]];

			for (int idx = transp[j]; idx < n - 1; idx++)  //������������ ������ � �����
			{
				b[idx] = b[idx + 1];
			}
			b[n - 1] = tmp;
		}
	}

	//1 step  Ly=b,  where Ux=y :
	y[0] = b[0];

	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			y[i] -= L[i][j] * y[j];
		}
		y[i] += b[i];
	}
	//2 step  Ux=y:
	for (int i = 0; i < n - 1; i++)x[i] = 0;
	x[n - 1] = y[n - 1] / U[n - 1][n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		for (int j = i + 1; j < n; j++)
		{
			x[i] -= U[i][j] * x[j];
		}
		x[i] = (x[i] + y[i]) / U[i][i];
	}
}


double zero_moment(double b) {
	return 1.25 * pow((-0.1 + b),(0.8));
}
double zero_moment(double a, double b) {
	return 1.25 * pow((-0.1 + b), (0.8)) -
		1.25 * pow((-0.1 + a), (0.8));
}
double first_moment(double b) {
	return (5.0 / 9) * pow((-0.1 + b), (0.8))*(0.125 + b);
}
double first_moment(double a, double b) {
	return (5.0 / 9) * pow((-0.1 + b), (0.8))*(0.125 + b) -
		(5.0 / 9) * pow((-0.1 + a), (0.8))*(0.125 + a);
}
double second_moment(double b) {
	return (5.0 / 14)*pow((-0.1 + b),(0.8))*(1.0 / 72 + 1.0/9 * b + b*b);
}
double second_moment(double a, double b) {
	return (5.0 / 14)*pow((-0.1 + b), (0.8))*(1.0 / 72 + 1.0 / 9 * b + b * b)-
		(5.0 / 14)*pow((-0.1 + a), (0.8))*(1.0 / 72 + 1.0 / 9 * a + a * a);
}

//
double third_moment(double b) {//6.310962062637622
	return 0.0726174 * pow((-0.5 + 5 * b), (0.8))*(0.0014881 + 0.0119048*b + 0.107143 *b*b + pow(b, 3));
}
double third_moment(double a, double b) {
	return (0.0726174 * pow((-0.5 + 5*b), (0.8))*(0.0014881 + 0.0119048*b + 0.107143 *b*b + pow(b, 3))) -
		(0.0726174 * pow((-0.5 + 5*a), (0.8))*(0.0014881 + 0.0119048*a + 0.107143 *a*a + pow(a,3)));
}
double fourth_moment(double b) {
	return 0.0574887 * pow((-0.5 + 5 * b), (0.8))*(0.000156642 + 0.00125313*b + 0.0112782*b*b + 0.105263*pow(b, 3) + pow(b, 4));
}
double fourth_moment(double a, double b) {
	return (0.0574887 * pow((-0.5 + 5 * b), (0.8))*(0.000156642 + 0.00125313*b + 0.0112782*b*b + 0.105263*pow(b, 3) + pow(b, 4))) -
		(0.0726174 * pow((-0.5 + 5 * a), (0.8))*(0.000156642 + 0.00125313*a + 0.0112782*a*a + 0.105263*pow(a, 3) + pow(a, 4)));
}
double fifth_moment(double b) {
	return 0.0475769 * pow((-0.5 + 5 * b), (0.8))*(0.0000163168 + 0.000130535*b + 0.00117481*b*b + 0.0109649*pow(b, 3) + 0.104167*pow(b, 4) + pow(b, 5));
}
double fifth_moment(double a, double b) {
	return (0.0475769 * pow((-0.5 + 5 * b), (0.8))*(0.0000163168 + 0.000130535*b + 0.00117481*b*b + 0.0109649*pow(b, 3) + 0.104167*pow(b, 4) + pow(b, 5)))-
		(0.0475769 * pow((-0.5 + 5 * a), (0.8))*(0.0000163168 + 0.000130535*a + 0.00117481*a*a + 0.0109649*pow(a, 3) + 0.104167*pow(a, 4) + pow(a, 5)));
}
double Integrate_Newton_Cotes_ikf(double a, double b, double alpha) {//исходный вариант
	const int nodes_number = 3;
	vector<double> x(nodes_number);
	vector<double> x_deg_k(nodes_number);
	for (int i = 0; i < nodes_number; i++) {
		x_deg_k[i] = 1;
		x[i] = a + i * (b-a) / (nodes_number - 1);
	}
	vector<double> moment(nodes_number);
	vector<double> A(nodes_number);
	vector<int> transporation(nodes_number);

	moment[0] = zero_moment(b);
	moment[1] = first_moment(b);
	moment[2] = second_moment(b);
	
	vector<vector<double>> matrix(nodes_number),L(nodes_number),U(nodes_number);
	for (int i = 0; i < nodes_number; i++)//инициализация L и U:
	{
		for (int j = 0; j < nodes_number; j++)
		{
			L[i].push_back(0);
			U[i].push_back(0);
			matrix[i].push_back(0);
		}
	}
	for (int i = 0; i < nodes_number; i++) {
		matrix[i][0] = x_deg_k[0];
		matrix[i][1] = x_deg_k[1];
		matrix[i][2] = x_deg_k[2];
		x_deg_k[0] *= x[0];
		x_deg_k[1] *= x[1];
		x_deg_k[2] *= x[2];
	}
	//решение СЛАУ:
	LU(matrix, L, U, nodes_number, transporation);
	slae_solution(L, U, moment, A, transporation, nodes_number);

	double res = 0;
	for (int i = 0; i < nodes_number; i++)
		res += sourse_function(x[i])*A[i];
	return res;
}
double Integrate_Newton_Cotes_ikf(double x1, double x2, double a, double alpha) {//если используется для скф
	const int nodes_number = 3;
	vector<double> x(nodes_number);
	vector<double> x_deg_k(nodes_number);
	for (int i = 0; i < nodes_number; i++) {
		x_deg_k[i] = 1;
		x[i] = x1 + i * (x2 - x1) / (nodes_number - 1);
	}
	vector<double> moment(nodes_number);
	vector<double> A(nodes_number);
	vector<int> transporation(nodes_number);

	moment[0] = zero_moment(x1, x2);
	moment[1] = first_moment(x1, x2);
	moment[2] = second_moment(x1, x2);

	vector<vector<double>> matrix(nodes_number), L(nodes_number), U(nodes_number);
	for (int i = 0; i < nodes_number; i++)//инициализация L и U:
	{
		for (int j = 0; j < nodes_number; j++)
		{
			L[i].push_back(0);
			U[i].push_back(0);
			matrix[i].push_back(0);
		}
	}
	for (int i = 0; i < 3; i++) {
		matrix[i][0] = x_deg_k[0];
		matrix[i][1] = x_deg_k[1];
		matrix[i][2] = x_deg_k[2];
		x_deg_k[0] *= x[0];
		x_deg_k[1] *= x[1];
		x_deg_k[2] *= x[2];
	}
	//решение СЛАУ:
	LU(matrix, L, U, nodes_number, transporation);
	slae_solution(L, U, moment, A, transporation, nodes_number);

	double res = 0;
	for (int i = 0; i < nodes_number; i++)
		res += sourse_function(x[i])*A[i];
	return res;
}
double Integrate_Newton_Cotes_skf(double a, double b, double alpha, double nodes_num) {
	vector<double> x(nodes_num);
	for (int i = 0; i < nodes_num; i++) {
		x[i] = a + i * (b - a) / (nodes_num-1);
	}
	double res = 0;
	for (int i = 0; i < nodes_num-1; i++)
		res += Integrate_Newton_Cotes_ikf(x[i], x[i + 1], a, alpha);
	return res;
}
double Integrate_Gaus_ikf(double a, double b) {
	const int nodes_number = 3;
	vector<double> x(nodes_number);
	vector<double> x_deg_k(nodes_number);
	for (int i = 0; i < nodes_number; i++) {
		x_deg_k[i] = 1;
		x[i] = a + i * (b - a) / (nodes_number - 1);
	}
	vector<double> moment(6);
	vector<double> A(nodes_number);
	vector<int> transporation(nodes_number);

	moment[0] = zero_moment(b);
	moment[1] = first_moment(b);
	moment[2] = second_moment(b);
	moment[3] = third_moment(b);
	moment[4] = fourth_moment(b);
	moment[5] = fifth_moment(b);

	vector<vector<double>> matrix(nodes_number), L(nodes_number), U(nodes_number), B(nodes_number);
	for (int i = 0; i < nodes_number; i++)//инициализация L и U:
	{
		for (int j = 0; j < nodes_number; j++)
		{
			L[i].push_back(0);
			U[i].push_back(0);
			B[i].push_back(0);
			matrix[i].push_back(0);
		}
	}
	for (int i = 0; i < nodes_number; i++) {
		for (int j = 0; j < nodes_number; j++) {
			matrix[i][j] = moment[i + j];	
		}
	}
	vector<double> _b(nodes_number);
	_b[0] = -moment[3];
	_b[1] = -moment[4];
	_b[2] = -moment[5];
	//решение СЛАУ:
	LU(matrix, L, U, nodes_number, transporation);
	slae_solution(L, U, _b, A, transporation, nodes_number);
	//узловой многочлен - кубичиский, его корни:
	double p, q, D; 
	p = A[1] - A[2] * A[2] / 3.0;
	q = A[0] + 2 * pow(A[2], 3) / 27.0 - A[1] * A[2] / 3.0;
	D = q*q/4.0 + pow(p, 3)/27.0;
	if (D > 0)cout << "NO!!!!!!!!!!!"<<endl;
	else {
		double fi = atan(sqrt(-D) / (-q / 2.0));
		if (q > 0) fi += PI / 2.0;
		if (abs(q) < EPS*0.001) fi = PI / 2.0;
		x[0] = 2 * sqrt(-p / 3)*cos(fi / 3.0) - A[2]/3.0;
		x[1] = 2 * sqrt(-p / 3)*cos(fi / 3.0 + 2 * PI / 3.0) - A[2] / 3.0;
		x[2] = 2 * sqrt(-p / 3)*cos(fi / 3.0 + 4 * PI / 3.0) - A[2] / 3.0;
	}
	//решение 2ого СЛАУ:

	for (int i = 0; i < nodes_number; i++) {
		A[i] = 0;
		matrix[i][0] = x_deg_k[0];
		matrix[i][1] = x_deg_k[1];
		matrix[i][2] = x_deg_k[2];
		x_deg_k[0] *= x[0];
		x_deg_k[1] *= x[1];
		x_deg_k[2] *= x[2];
	}
	LU(matrix, L, U, nodes_number, transporation);
	slae_solution(L, U, moment, A, transporation, nodes_number);
	double res = 0;
	for (int i = 0; i < nodes_number; i++)
		res += sourse_function(x[i])*A[i];
	return res;
}
double Integrate_Gaus_ikf(double x1, double x2, double a) {
	const int nodes_number = 3;
	vector<double> x(nodes_number);
	vector<double> x_deg_k(nodes_number);
	for (int i = 0; i < nodes_number; i++) {
		x_deg_k[i] = 1;
		x[i] = x1 + i * (x2 - x1) / (nodes_number - 1);
	}
	vector<double> moment(6);
	vector<double> A(nodes_number);
	vector<int> transporation(nodes_number);

	moment[0] = zero_moment(x1,x2);
	moment[1] = first_moment(x1,x2);
	moment[2] = second_moment(x1,x2);
	moment[3] = third_moment(x1,x2);
	moment[4] = fourth_moment(x1,x2);
	moment[5] = fifth_moment(x1,x2);

	vector<vector<double>> matrix(nodes_number), L(nodes_number), U(nodes_number), B(nodes_number);
	for (int i = 0; i < nodes_number; i++)//инициализация L и U:
	{
		for (int j = 0; j < nodes_number; j++)
		{
			L[i].push_back(0);
			U[i].push_back(0);
			B[i].push_back(0);
			matrix[i].push_back(0);
		}
	}
	for (int i = 0; i < nodes_number; i++) {
		for (int j = 0; j < nodes_number; j++) {
			matrix[i][j] = moment[i + j];	
		}
	}
	vector<double> _b(nodes_number);
	_b[0] = -moment[3];
	_b[1] = -moment[4];
	_b[2] = -moment[5];
	//решение СЛАУ:
	LU(matrix, L, U, nodes_number, transporation);
	slae_solution(L, U, _b, A, transporation, nodes_number);
	//узловой многочлен - кубичиский, его корни:
	double p, q, D; 
	p = A[1] - A[2] * A[2] / 3.0;
	q = A[0] + 2 * pow(A[2], 3) / 27.0 - A[1] * A[2] / 3.0;
	D = q*q/4.0 + pow(p, 3)/27.0;
	if (D > 0)cout << "NO!!!!!!!!!!!"<<endl;
	else {
		double fi = atan(sqrt(-D) / (-q / 2.0));
		if (q > 0) fi += PI / 2.0;
		if (abs(q) < EPS*0.001) fi = PI / 2.0;
		x[0] = 2 * sqrt(-p / 3)*cos(fi / 3.0) - A[2]/3.0;
		x[1] = 2 * sqrt(-p / 3)*cos(fi / 3.0 + 2 * PI / 3.0) - A[2] / 3.0;
		x[2] = 2 * sqrt(-p / 3)*cos(fi / 3.0 + 4 * PI / 3.0) - A[2] / 3.0;
	}
	//решение 2ого СЛАУ:

	for (int i = 0; i < nodes_number; i++) {
		A[i] = 0;
		matrix[i][0] = x_deg_k[0];
		matrix[i][1] = x_deg_k[1];
		matrix[i][2] = x_deg_k[2];
		x_deg_k[0] *= x[0];
		x_deg_k[1] *= x[1];
		x_deg_k[2] *= x[2];
	}
	LU(matrix, L, U, nodes_number, transporation);
	slae_solution(L, U, moment, A, transporation, nodes_number);
	double res = 0;
	for (int i = 0; i < nodes_number; i++)
		res += sourse_function(x[i])*A[i];
	return res;
}
double Integrate_Gaus_skf(double a, double b, double nodes_num) {
	vector<double> x(nodes_num);
	for (int i = 0; i < nodes_num; i++) {
		x[i] = a + i * (b - a) / (nodes_num - 1);
	}
	double res = 0;
	for (int i = 0; i < nodes_num - 1; i++)
		res += Integrate_Gaus_ikf(x[i], x[i + 1], a);
	return res;
}
int main() {
	setlocale(LC_ALL, "Ru");
	cout << "ИКФ Ньютона-Котса" << endl;
	double value = Integrate_Newton_Cotes_ikf(0.1, 2.3, 0.2);
	cout <<value<< endl;

	cout << "abs error= " << abs(INTEGRAL_VALUE - value)<<endl;

	cout << "CКФ Ньютона-Котса" << endl;
	double value1 = 1;
	double value2 = 1;
	double delta = 1;
	double h_opt = 0;
	double h_1, h_2, h_3;
	int n_opt = 1;
	int i = 3;
	double m;
	while (delta > EPS) {
		value = Integrate_Newton_Cotes_skf(0.1, 2.3, 0.2, i);
		delta = abs(INTEGRAL_VALUE - value);
		i++;
	}
	cout <<std::setprecision(7)<<value << endl;
	cout << "шагов: " << i << endl;
	
	//погрешность Ридчарсона:
	// "CКФ Ньютона-Котса, 2 грубые стеки:"
	//везде далее используется L = 2
	value = Integrate_Newton_Cotes_skf(0.1, 2.3, 0.2, 2);
	value1 = value;
	h_1 = 2.2 / 3;
	h_2 = h_1/ 2;
	h_3 = h_1 / 4;
	m = 4;//см далее, L=2, 2^4=16

	value = Integrate_Newton_Cotes_skf(0.1, 2.3, 0.2, 3);
	value2 = value;
	delta = (value2 - value1) / 15;
	cout << "погрешность: " <<delta<< endl;

	h_opt = (EPS / delta) * 0.95 * h_1;
	h_opt = (2.2) / ceil(2.2 / h_opt);
	cout << "Hopt= "<< h_opt << endl;
	cout <<"Nopt= "<< (2.2) / h_opt + 1;
	cout << " шагов: " << (int)(log10((2.2) / h_opt) / log10(3)) << endl;
	i = 8;
	while (delta > EPS) {
		value = Integrate_Newton_Cotes_skf(0.1, 2.3, 0.2, i);
		delta = abs(INTEGRAL_VALUE - value);
		i++;
	}
	cout << value <<" шагов= "<< i <<endl;
	value1 = Integrate_Newton_Cotes_skf(0.1, 2.3, 0.2, 5);// можно брать в качестве N любое число вида 2^k+1
	value2 = Integrate_Newton_Cotes_skf(0.1, 2.3, 0.2, 9);
	value = Integrate_Newton_Cotes_skf(0.1, 2.3, 0.2, 17);
	m = -log((value - value2) / (value2 - value1))/log(2);	

	cout << "скорость сходимости= " << m << endl;//m->4 при шаге->0

	cout << "КФ Гаусса:" << Integrate_Gaus_ikf(0.1, 2.3)<< endl;
	system("pause");
}