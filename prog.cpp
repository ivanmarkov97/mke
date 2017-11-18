/*
A(d2U/dx2) - BU + C = 0
U(5) = 0
dU/dx(50) = 1
Nl = 20
Nl = 40
L = (50 -5)/Nl
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
using namespace std;

#define _AREA_LEN 45.0
#define _CONSTR_TYPE_RIGHT 2
#define _CONSTR_TYPE_LEFT 1
#define _CONSTR_R_VALUE 1.0
#define _CONSTR_L_VALUE 0.0
#define _A 51.0
#define _B 8.0
#define C 20.0
//#define _LINER_FUNC 

class
Constrain{
	private:
		int type;
		double value;
	public:
		Constrain(): type(0), value(0.0){};
		Constrain(int _t, double _v): type(_t), value(_v){};
		int get_type() { return type; };
		double get_value() { return value; };
};

class
Element{
	private:
		double len;
		double ui;
		double uj;
		int size;
		Constrain l_constr;
		Constrain r_constr;
		double** rigit_matx;
		double* free_part;
	public:
		Element(): len(0.0), ui(0.0), uj(0.0), rigit_matx(NULL), free_part(NULL) { l_constr = Constrain(); r_constr = Constrain(); };
		Element(double _l, double _i, double _j, Constrain _lc, Constrain _rc){
			len = _l;
			ui = _i;
			uj = _j;
			l_constr = _lc;
			r_constr = _rc;
			#ifdef _LINER_FUNC
			size = 2;
			rigit_matx = (double**)calloc(2, sizeof(double*));
			free_part = (double*)calloc(2, sizeof(double));
			for(int i = 0; i < 2; i++){
				rigit_matx[i] = (double*)calloc(2, sizeof(double));
			}
			#else
			size = 4;
			rigit_matx = (double**)calloc(4, sizeof(double*));
			free_part = (double*)calloc(4, sizeof(double));
			for(int i = 0; i < 4; i++){
				rigit_matx[i] = (double*)calloc(4, sizeof(double));
			}
			#endif
		};
		double get_len() { return len; };
		double get_ui() { return ui; };
		double get_uj() { return uj; };
		int get_size() { return size; };
		double** get_rigit_matx() { return rigit_matx; };
		double* get_free_part() { return free_part; };
		void show_rigit_matx(){
			for(int i = 0; i < size; i++){
				for(int j = 0; j < size; j++){cout << rigit_matx[i][j] << " ";}
				cout << endl;
			}
		};
		void init_rigit_matx(){
			#ifdef _LINER_FUNC
			for(int i = 0; i < size; i++){
				rigit_matx[i][i] =  -(_A / len) - (_B * len / 3.0);
				rigit_matx[size - (i + 1)][i] = (_A / len) - (_B * len / 6.0);
			}
			#else
			rigit_matx[0][0] = rigit_matx[size - 1][size - 1] = -(188.7 / len) - (0.608 * len);
			rigit_matx[1][1] = rigit_matx[size - 2][size - 2] = -(550.8 / len) - (3.088 * len);

			rigit_matx[0][1] = rigit_matx[1][0] = (240.985 / len) - (0.472 * len);
			rigit_matx[0][2] = rigit_matx[2][0] = -(68.85 / len) + (0.168 * len);
			rigit_matx[0][3] = rigit_matx[3][0] = (16.575 / len) - (0.088 * len);

			rigit_matx[2][1] = rigit_matx[1][2] = (378.675 / len) + (0.384 * len);
			rigit_matx[3][1] = rigit_matx[1][3] = -(68.85 / len) + (0.168 * len);

			rigit_matx[2][3] = rigit_matx[3][2] = (240.975 / len)  - (0.472 * len);

			/*rigit_matx[0][0] = 4.0 - 1887.0/(10.0*len);
			rigit_matx[0][1] = 9639.0/(40.0*len) - 57.0/10.0;
			rigit_matx[0][2] = 12.0/5.0 - 1337.0/(20.0*len);
			rigit_matx[0][3] = 663.0/(40.0*len) - 0.7;
			
			rigit_matx[1][0] = 9639.0/(40.0*len) + 57.0/10.0;
			rigit_matx[1][1] = -2754.0/(5.0*len);
			rigit_matx[1][2] = 15147.0/(40.0*len) - 8.1;
			rigit_matx[1][3] = 12.0/5.0 - 1377.0/(20.0*len);

			rigit_matx[2][0] = -1377.0/(20.0*len) - 12.0/5.0;
			rigit_matx[2][1] = 15147.0/(40.0*len) + 8.1;
			rigit_matx[2][2] = -2754.0/(5.0*len);
			rigit_matx[2][3] = 9639.0/(40.0*len) - 5.7;
	
			rigit_matx[3][0] = 663.0/(40.0*len) + 0.7;
			rigit_matx[3][1] = -1377.0/(20*len) - 12.0/5.0;
			rigit_matx[3][2] = 9639.0/(40.0*len) + 5.7;
			rigit_matx[3][3] = -1887.0/(10.0*len) - 4.0;*/

			#endif
		};
		void init_free_part(){
			#ifdef _LINER_FUNC
			for(int i = 0; i < size; i++){ free_part[i] = -C * len / 2.0; }
			#else
			for(int i = 1; i < size - 1; i++){ free_part[i] = -7.5 * len; }
			free_part[0] = free_part[size - 1] = -2.5 * len;
			#endif
		};
		void show_free_part(){
			for(int i = 0; i < size; i++)
				cout << free_part[i] << endl;
		}
		Constrain get_l_constr() {return l_constr; };
		Constrain get_r_constr() {return r_constr; };
};

class 
EquationSystem{
	private:
		double** global_rigit_matx;
		double* global_free_part;
		int size;
	public: 
		EquationSystem(): global_rigit_matx(NULL), global_free_part(NULL), size(0) {};
		EquationSystem(int s){
			size = s;
			for(int i = 0; i < size; i++){
				global_rigit_matx = (double**)calloc(size, sizeof(double*));
				global_free_part = (double*)calloc(size, sizeof(double));
				for(int j = 0; j < size; j++){
					global_rigit_matx[j] = (double*)calloc(size, sizeof(double));	
				}
			}

		};

		double** get_global_rigit_matx(){ return global_rigit_matx; };
		double* get_global_free_part(){ return global_free_part; };
		int get_size(){ return size; };

		void show_global_rigit_matx(){
			for(int i = 0; i < size; i++){
				for(int j = 0; j < size; j++){
					cout << global_rigit_matx[i][j] << "\t";
				}
				cout << endl;
			}
		};
		void show_global_free_part(){
			for(int i = 0; i < size; i++){
				cout << global_free_part[i] << endl;
			}
		};
		void init_global_rigit_matx(Element elem){
			double** local_rigit_matx = elem.get_rigit_matx();
			int el_size = elem.get_size();
			for(int i = 0; i < size - el_size + 1; i++){
				for(int j = 0; j < el_size; j++){
					for(int k = 0; k < el_size; k++){
						global_rigit_matx[i + j][i + k] += local_rigit_matx[j][k];
					}
				}
			}
			#ifdef _LINER_FUNC
			global_rigit_matx[0][0] = 1.0;
			global_rigit_matx[1][0] = 1.0;
			#else
			global_rigit_matx[0][0] = 2.0;
			global_rigit_matx[1][0] = 2.0;
			global_rigit_matx[2][0] = 2.0;
			global_rigit_matx[3][0] = 2.0;
			#endif
		}
		void init_global_free_part(Element elem){
			double* local_free_part = elem.get_free_part();
			int el_size = elem.get_size();
			for(int i = 0; i < size - el_size + 1; i++){
				for(int j = 0; j < el_size; j++){
					global_free_part[i + j] += local_free_part[j];
				}
			}
			#ifdef _LINER_FUNC
			global_free_part[size - 1] += -_A;
			#else
			global_free_part[size - 1] += -_A;
			global_free_part[0] += -(3.7/el_size);
			global_free_part[1] += -(-4.725/el_size);
			global_free_part[2] += -(1.35/el_size);
			global_free_part[3] += -(-0.325/el_size);
			/*global_free_part[0] += -(4.0 - 1887.0/(10.0*el_size));
			global_free_part[1] += -(9639.0/(40.0*el_size) + 5.7);
			global_free_part[2] += -(-1377.0/(20.0*el_size) - 12.0/5.0);
			global_free_part[3] += -(663.0/(40.0*el_size) + 0.7);*/
			#endif
		};
};


double* gauss(double **a, double *y, int n) 
{
	double *x, max;
  	int k, index;
  	const double eps = 0.0000000000000000000000000001;  // точность
  	x = new double[n];
  	k = 0;
  	while (k < n) 
  	{
    	// Поиск строки с максимальным a[i][k]
    		max = fabs(a[k][k]);
    		index = k;
    		for (int i = k + 1; i < n; i++) {
      			if (fabs(a[i][k]) > max){
        			max = fabs(a[i][k]);
        			index = i;
      			}
    		}
    		// Перестановка строк
    		if (max < eps) {
      			// нет ненулевых диагональных элементов
      			cout << "Решение получить невозможно из-за нулевого столбца ";
      			cout << index << " матрицы A" << endl;
      			return 0;
    		}
    		for (int j = 0; j < n; j++) {
      			double temp = a[k][j];
      			a[k][j] = a[index][j];
      			a[index][j] = temp;
    		}
    		double temp = y[k];
    		y[k] = y[index];
    		y[index] = temp;
    		// Нормализация уравнений
    		for (int i = k; i < n; i++) {
      			double temp = a[i][k];
      			if (fabs(temp) < eps)  continue; // для нулевого коэффициента пропустить
      			for (int j = 0; j < n; j++) 
        			a[i][j] = a[i][j] / temp;
      			y[i] = y[i] / temp;
      			if (i == k)  continue; // уравнение не вычитать само из себя
      			for (int j = 0; j < n; j++)
        		a[i][j] = a[i][j] - a[k][j];
      			y[i] = y[i] - y[k];
    		}
    		k++;
  	}
  	// обратная подстановка
  	for (k = n - 1; k >= 0; k--){
    		x[k] = y[k];
    		for (int i = 0; i < k; i++)
      			y[i] = y[i] - a[i][k] * x[k];
  	}
	#ifndef _LINER_FUNC
	if(n > 2)
		x[n - 1] = 0.9384 * 45.0 / n + x[n - 2];
	#endif
	x[0] = 0.0;
  	return x;
}


int 
main(int argc, char** argv){

	int num_elem = 0;
	double len_elem = 0.0;

	cout << "num elements: Usage: 20/40" << endl;
	cin >> 	num_elem;
	len_elem = _AREA_LEN / num_elem;
	cout << "element len: " << len_elem << endl;

	Constrain l_constr = Constrain(_CONSTR_TYPE_LEFT, _CONSTR_L_VALUE);
	Constrain r_constr = Constrain(_CONSTR_TYPE_RIGHT, _CONSTR_R_VALUE);

	//cout << l_constr.get_type() << " " << l_constr.get_value() / 2 << endl;
	//cout << r_constr.get_type() << " " << r_constr.get_value() / 2 << endl;

	Element test_elem = Element(len_elem, 0.0, 0.0, l_constr, r_constr);
	//cout << test_elem.get_r_constr().get_value() / 2 << endl;
	test_elem.init_rigit_matx();
	test_elem.init_free_part();
	test_elem.show_rigit_matx();
	test_elem.show_free_part();

	cout << endl;

	#ifdef _LINER_FUNC
	EquationSystem eq_system = EquationSystem(num_elem + 1);
	#else
	EquationSystem eq_system = EquationSystem(num_elem + 3);
	#endif
	eq_system.show_global_rigit_matx();

	cout << endl;
	
	eq_system.init_global_rigit_matx(test_elem);
	eq_system.show_global_rigit_matx();

	cout << endl;
	
	eq_system.init_global_free_part(test_elem);
	eq_system.show_global_free_part();

	cout << endl;
	
	FILE* file;
	if((file = fopen("data.txt", "w"))== NULL){
		cout << "error open file" << endl;
	}
	#ifdef _LINER_FUNC
	fprintf(file, "linear\n");
	fprintf(file, "%d\n", num_elem);
	#else
	fprintf(file, "cube\n");
	fprintf(file, "%d\n", num_elem + 2);
	#endif
	double* result = gauss(eq_system.get_global_rigit_matx(), eq_system.get_global_free_part(), eq_system.get_size());
	for(int i = 0; i < eq_system.get_size(); i++){
		cout << result[i] << "," << endl;
		fprintf(file, "%lf,\n", result[i]);
	}
	fclose(file);
	return 0;
}