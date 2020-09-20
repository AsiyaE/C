#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 0.00001

//Metod Gaussa
using namespace std;

void chtenie(FILE *inp, int &n1, int &n2);                            //chtenie razmera matrix
double ** getmemory(double **arr, int n1, int n2);                   //videlenie pamiti pod massiv
void record(FILE *inp, double **arr, int n1, int n2);                 //zapis koef v matrix
double ** dano(double **arr, int &n1, int &n2);                      //matrix koef
void pechat(double **arr, int n1, int n2);                           //pechat matrix koef
void clean(double **arr, int n1);                                    //clean memory
int poisk(double **arr, int i, int j, int n1);                       //poisk nenulevoi stroki
void perestanovka(double **arr, int i, int no0);                     //perestanovka 2x strok
void subtract(double**arr, int n1, int n2, int k);                   //vichitanie tekyshei stroki ymnogennoi na koef
int proverka(double **arr, int n1, int n2);                          //proverka 
int oper_on_matrix(double **arr, int n1, int n2);                    //diagonal matrix
double * decision(double **arr, double *x, int n1, int n2, int flag); //poisk reshenia
double * memory(double *x, int n1);                                  //videlenie pamiti X
void pechat_x(double *x, int n1, int flag);                          //pechat x

/*chtenie razmera matrix
*/
void chtenie(FILE *inp, int &n1, int &n2) {
	fscanf_s(inp, "%d %d", &n1, &n2);
}

/*videlenie pamiti pod massiv
	int **arr - ykazatel
	int n1 - kol-vo strok
	int n2 - kol-vo stolbcov
*/
double ** getmemory(double **arr, int n1, int n2) {
	arr = (double**)calloc(n1, sizeof(double*));
	for (int i = 0; i < n1; i++)
	{
		arr[i] = (double*)calloc(n2, sizeof(double));
	}
	return arr;
}

/*zapis koef v matrix
*/
void  record(FILE *inp, double **arr, int n1, int n2) {
	float tmp;
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			fscanf_s(inp, "%f", &tmp);
			arr[i][j] = (double)tmp;
		}
	}
}

/*matrix koef
*/
double ** dano(double **arr, int &n1, int &n2) {
	FILE *inp;
	fopen_s(&inp, "input.txt", "r");
	if (inp == NULL) printf("%s \n", "file error");
	chtenie(inp, n1, n2);
	arr = getmemory(arr, n1, n2);
	record(inp, arr, n1, n2);
	fclose(inp);
	return arr;
}

/*pechat matrix koef
*/
void pechat(double **arr, int n1, int n2) {

	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			printf("%e\t", arr[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/*clean memory
*/
void clean(double **arr, int n1) {
	for (int i = 0; i < n1; i++)
		free(arr[i]);
	free(arr);
}

/*poisk nenulevoi stroki
*/
int poisk(double **arr, int i, int j, int n1) {
	for (int k = i + 1; k < n1; k++) {
		if (fabs(arr[k][j]) >= EPS) {
			return k;
		}
	}
	return(-1); //esli takoi stroki net
}

/*perestanovka 2x strok
*/
void perestanovka(double **arr, int i, int no0) {
	double *temp = arr[i];
	arr[i] = arr[no0];
	arr[no0] = temp;
}

/*vichitanie tekyshei stroki ymnogennoi na koef
*/
void subtract(double**arr, int n1, int n2, int k) {
	double mi;
	for (int i = k + 1; i < n1; i++) {
		mi = arr[i][k] / arr[k][k];
		for (int j = k; j < n2; j++) {
			arr[i][j] -= mi * arr[k][j];
			if (fabs(arr[i][j]) < EPS) {    //EPS=0.00001= 1.0*e-05
				arr[i][j] = 0;
			}
		}
	}
}

/*proverka 
*/
int proverka(double **arr, int n1, int n2) {
	if (n1 < n2 - 1) {                       //mnogo reshenii
		return 0;
	}else if (arr[n1 - 1][n2 - 2] == 0) {    //last diagonal element =0 
		return (-1);
	}else {
		return 1;
	}
}
/*diagonal matrix
return: -1 sistema virozgdena!
        0 sistema sovmestna i neopredelena (mnogo reshenii)
		1 est reshenie
*/
int oper_on_matrix(double **arr, int n1, int n2) {
	int f;
	for (int i = 0; i < n1 - 1; i++) {
		for (int j = 0; j < n2; j++) {
			if (i == j) {
				if (fabs(arr[i][j]) < EPS) {
					int no0 = poisk(arr, i, j, n1);  //poisk nenulevoi stroki
					if (no0 != -1) {
						perestanovka(arr, i, no0);
					}
					else {
						return (-1);  //nuli na diagonali
					}
				}
				subtract(arr, n1, n2, i);
			}
		}
	}
	f=proverka(arr, n1, n2); 
	return f;
}


/*poisk reshenia
*/
double * decision(double **arr, double *x, int n1, int n2, int flag) {
	if (flag == 1) {
		double tm;
		x = memory(x, n1);
		x[n1 - 1] = arr[n1 - 1][n2 - 1] / arr[n1 - 1][n2 - 2];
		for (int i = n1 - 2; i >= 0; i--) {
			tm = arr[i][n1];
			for (int j = i + 1; j < n1; j++) {
				tm = tm - arr[i][j] * x[j];
			}
			x[i] = tm / arr[i][i];
		}
	}
	return x;
}

/*videlenie pamiti X
*/
double * memory(double *x, int n1) {
	x = (double*)calloc(n1, sizeof(double));
	for (int i = 0; i < n1; i++) {
		x[i] = 0;
	}
	return x;
}

/*pechat x
*/
void pechat_x(double *x, int n1, int flag) {
	if (flag == 1) {
		printf("%s \n\n", "reshenie sistemi: ");
		for (int i = 0; i < n1; i++) {
			printf("%s%d %s %e \n", "X", i + 1, " = ", x[i]);
		}
	}else if (flag == 0) {
		printf("%s \n\n", "sistema sovmestna i neopredelena (mnogo reshenii)!");
	}else {
		printf("%s \n", "sistema virozgdena!");
	}

}

int main()
{
	int flag, n1 = 0, n2 = 0;
	double **arr = NULL;

	arr = dano(arr, n1, n2);
	printf("%s \n\n", "matrix koeff:");
	pechat(arr, n1, n2);

	flag = oper_on_matrix(arr, n1, n2);
	if (flag == 1) {
		printf("%s \n\n", "diagonal matrix:");
	}
	else { printf("%s \n\n", "diagonal matrix net"); }
	pechat(arr, n1, n2);

	double *x = NULL;
	x = decision(arr, x, n1, n2, flag);
	pechat_x(x, n1, flag);

	clean(arr, n1);
	free(x);
}