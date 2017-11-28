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
#define _LINER_FUNC 

class
Element{
	private:
		double len;
		int size;
		double** rigit_matx;
		double* free_part;
	public:
		Element(): len(0.0),rigit_matx(NULL), free_part(NULL) {};
		Element(double _l){
			len = _l;
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
			rigit_matx[0][0] = -8.0 * len / 3.0 - 51.0 / len;
			rigit_matx[0][1] = 51.0 / len - 4.0 * len / 3.0;
			rigit_matx[1][0] = 51.0 / len - 4.0 * len / 3.0;
			rigit_matx[1][1] = -8.0 * len / 3.0 - 51.0 / len;
			#else

			/*rigit_matx[0][0] = -64.0*len/105.0 - 1887.0/(10.0*len);
			rigit_matx[0][1] = 9639.0/(40.0*len) - 33.0*len/70.0;
			rigit_matx[0][2] = 6.0*len/35.0 - 1337.0/(20.0*len);
			rigit_matx[0][3] = 663.0/(40.0*len) - 19.0*len/210.0;
			
			rigit_matx[1][0] = 9639.0/(40.0*len) - 33.0*len/70.0;
			rigit_matx[1][1] = -108.0*len/35.0 - 2754.0/(5.0*len);
			rigit_matx[1][2] = 27.0*len/70.0 + 15147.0/(40*len);
			rigit_matx[1][3] = 6.0*len/35.0 - 1377.0/(20.0*len);

			rigit_matx[2][0] = 6.0*len/35.0 - 1377.0/(20.0*len);
			rigit_matx[2][1] = 27.0*len/70.0 + 15147.0/(40.0*len);
			rigit_matx[2][2] = -108.0*len/35.0 - 2754/(5.0*len);
			rigit_matx[2][3] = 9639.0/(40.0*len) - 33.0*len/70.0;
	
			rigit_matx[3][0] = 663.0/(40.0*len) - 19.0*len / 210.0;
			rigit_matx[3][1] = 6.0*len/35.0 - 1377/(20.0*len);
			rigit_matx[3][2] = 9639.0/(40.0*len) - 33.0*len/70.0;
			rigit_matx[3][3] = -64.0*len/105.0 - 1887.0/(10.0*len);
			*/
			rigit_matx[0][0] = rigit_matx[size - 1][size - 1] = -(188.7 / len) - (0.608 * len);
                        rigit_matx[1][1] = rigit_matx[size - 2][size - 2] = -(550.8 / len) - (3.088 * len);

                        rigit_matx[0][1] = rigit_matx[1][0] = (240.985 / len) - (0.472 * len);
                        rigit_matx[0][2] = rigit_matx[2][0] = -(68.85 / len) + (0.168 * len);
                        rigit_matx[0][3] = rigit_matx[3][0] = (16.575 / len) - (0.088 * len);

                        rigit_matx[2][1] = rigit_matx[1][2] = (378.675 / len) + (0.384 * len);
                        rigit_matx[3][1] = rigit_matx[1][3] = -(68.85 / len) + (0.168 * len);

                        rigit_matx[2][3] = rigit_matx[3][2] = (240.975 / len)  - (0.472 * len);

			#endif
		};
		void init_free_part(){
			#ifdef _LINER_FUNC
			for(int i = 0; i < size; i++){ free_part[i] = -20.0 * len / 2.0; }
			#else
			for(int i = 1; i < size - 1; i++){ free_part[i] = -7.5 * len; }
			free_part[0] = free_part[size - 1] = -2.5 * len;
			#endif
		};
		void show_free_part(){
			for(int i = 0; i < size; i++)
				cout << free_part[i] << endl;
		}
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
			global_free_part[size - 1] += -51.0;
			#else
			global_free_part[size - 1] += -51.0*1.5;
			#endif
		};
};


double* gauss(double **a, double *y, int n) 
{
	double *x, max;
  	int k, index;
  	const double eps = 0.00000000000000000000000000000000001;  // точность
  	x = new double[n];
  	k = 1;
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
	#endif
	x[0] = 0.0;
  	return x;
}

double* origin_func(int num_elem){
	double *x = (double*)calloc(num_elem, sizeof(double));
	double l = 0;
	int i = 0;
	double step = _AREA_LEN / num_elem;
	#ifndef _LINER_FUNC
	num_elem += 2;
	step = _AREA_LEN / num_elem;
	#endif
        printf("num == %d\n", num_elem);
	while(i <= num_elem){
		x[i] = 2.5;
		double a = exp(-2.0*sqrt(2.0/51.0)*l);
		double b = (exp(2.0*sqrt(2.0/51.0)*l) - 1.0);
		double c = -10.0*exp(2.0*sqrt(2.0/51.0)*l);
		double d = sqrt(102.0)*exp(2.0*sqrt(2.0/51.0)*(l+45.0));
		double e = sqrt(102.0)*exp(30.0*sqrt(6.0/17.0));
		double f = 10.0*exp(60.0*sqrt(6.0/17.0));
		double g = (4.0*(1.0 + exp(60.0*sqrt(6.0/17.0))));

		//printf ("%lf %lf %lf %lf %lf %lf %lf | %lf\n", a,b,c,d,e,f,g, (a*b*(c+d+e+f)/g));

		x[i] = ( exp(-2.0*sqrt(2.0/51.0)*l)  *
			 (exp(2.0*sqrt(2.0/51.0)*l) - 1.0) * 
			 (-10.0*exp(2.0*sqrt(2.0/51.0)*l) + 
							sqrt(102.0)*exp(2.0*sqrt(2.0/51.0)*(l+45.0)) + 
							sqrt(102.0)*exp(30.0*sqrt(6.0/17.0)) + 
							10.0*exp(60.0*sqrt(6.0/17.0))
			 )
			) / (4.0*(1.0 + exp(60.0*sqrt(6.0/17.0))));
		i++;
		l += step;
	}
	for(int i = 0; i < num_elem; i++)
		;//printf("f[i] == %lf\n",x[i]);
	return x;
}

double search_max_error(double* x, double *f, int size){
	double max = 0.0;
	for(int i = 0; i < size; i++){
		//printf("raznica %lf\n", fabs( (x[i] - f[i]) ));
		if( fabs( (x[i] - f[i]) ) >= max ){
			max = fabs( (x[i] - f[i]) );
			//printf("pishem max");
		}
	}
	return max;
}

int 
main(int argc, char** argv){

	int num_elem = 0;
	double len_elem = 0.0;

	cout << "num elements: Usage: 20/40" << endl;
	cin >> 	num_elem;
	len_elem = _AREA_LEN / num_elem;
	cout << "element len: " << len_elem << endl;

	Element test_elem = Element(len_elem);
	test_elem.init_rigit_matx();
	test_elem.init_free_part();
	test_elem.show_rigit_matx();
	test_elem.show_free_part();

	//cout << endl;

	#ifdef _LINER_FUNC
	EquationSystem eq_system = EquationSystem(num_elem + 1);
	#else
	EquationSystem eq_system = EquationSystem(num_elem + 3);
	#endif
	//eq_system.show_global_rigit_matx();

	//cout << endl;
	
	eq_system.init_global_rigit_matx(test_elem);
	//eq_system.show_global_rigit_matx();

	//cout << endl;
	
	eq_system.init_global_free_part(test_elem);
	//eq_system.show_global_free_part();

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
	double *f = origin_func(num_elem);
	#ifdef _LINER_FUNC
	printf("MAX error liner == %lf | %d\n", search_max_error(result, f, num_elem), num_elem);
	#else
	printf("MAX error cube == %lf | %d\n", search_max_error(result, f, num_elem + 2), num_elem + 2);
	#endif
	for(int i = 0; i < eq_system.get_size(); i++){
		cout << result[i] << "," << endl;
		fprintf(file, "%lf,\n", result[i]);
	}
	fclose(file);
	return 0;
}
