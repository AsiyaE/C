#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 0.00001
static const int COUNT = 4; //glubina proverki sxodimosti max

//Iterative method


double ** dano(double **arr, int &n1, int &n2);                         //matrix koef
void chtenie(FILE *inp, int &n1, int &n2);                              //chtenie razmera matrix
double ** getmemory(double **arr, int n1, int n2);                      //videlenie pamiti pod massiv
void record(FILE *inp, double **arr, int n1, int n2);                   //zapis koef v matrix
void pechat(double **arr, int n1, int n2);                              //pechat matrix koef
void clean(double **arr, int n1);                                       //clean memory
double * memory(double *x, int n1);                                     //videlenie pamiti X
void pechat_x(double *x, int n1);                                       //pechat x
bool nul_diag(double **arr, int n1, int n2);                            //proverka nulei na diagonali
int poisk_n0(double **arr, int i, int n1);                              //poisk nenulevoi stroki
double ** DUS(double **arr, int n1, int n2,bool &flagD);                                 //yslovie DUS
bool poisk(int *a, int n);                                               //poisk stroki diagonal element >= sum
void swap(int *a, int i, int j);
bool proverka(int *a, double **arr, int n1, int n2);
double ** change(double **arr, int *a, int n1);
double * decision(double **arr, double *x, int n1, int n2, int flagD);   //poisk reshenia sistem
void substitute(double **arr, double *x, double *lastx, int n1, int n2); //schitaem tekushee x
double max_diff(double *x, double *lasx, int n1, double *maxk10);        //max raznica x-x'
double sumstr(double **arr, int n2, int i);                             //sum elements stroki
bool check(double *maxk10);                                             //proverka monotonnosti dlya poslednix max
void perestanovka(double **arr, int i, int no0);                        //perestanovka 2x strok             

int main()
{
	int n1 = 0, n2 = 0;
	bool flagD;
	double **arr = NULL;
	double *x = NULL;

	arr = dano(arr, n1, n2);
	printf(" %s \n\n", "matrix koeff:");
	pechat(arr, n1, n2);

	arr = DUS(arr, n1, n2,flagD);
	if (flagD) printf(" %s \n\n", "DUS");
	else   	printf(" %s \n\n", "NE DUS");
	pechat(arr, n1, n2);

	x = decision(arr, x, n1, n2, flagD);
	if (x != NULL) {
		printf(" %s \n\n", "reshenie:");
		pechat_x(x, n1);
		free(x);
	}
	else
		printf(" %s \n\n", "ETIM METODOM RESHIT' NEL'ZYA");

	clean(arr, n1);
}

/*return -0 NE DUS
1- DUS
*/
double ** DUS(double **arr, int n1, int n2, bool &flagD) {
	int *a =(int*)calloc(n1, sizeof(int));     //dlya zapisi poryadka strok
	for (int i = 0; i < n1; i++) {
		a[i] = i;
	}
	while (!proverka(a, arr, n1, n2)) {  //vipolnyaetsya DUS?
		flagD = poisk(a, n1);            //poisk novoi perestanovki
		if (!flagD)return arr;           //bolshe net perestanovok
	}
	arr=change(arr, a, n1);     //pomenyat' stroki massiva arr
	return arr;
}

double ** change(double **arr, int *a, int n1) {
	double **b = NULL;
	b = (double**)calloc(n1, sizeof(double*));
	for (int i =0; i < n1; i++) {
		b[i] = arr[a[i]];
	}
	free(a);
	return b;
}


/*poisk set*/
bool poisk(int *a, int n) {    
	int j = n - 2;
	while ((j != -1) && (a[j] >= a[j + 1]))
		j--;
	if (j == -1)
		return false;  //bolshe variantov net
	int k = n - 1;
	while (a[j] >= a[k])
		k--;
	swap(a, j, k);
	int l = j + 1, r = n - 1;
	while (l < r)
		swap(a, l++, r--);
	return true;  //naidena novaya perestanovka
}

/*perestavit elementi*/
void swap(int *a, int i, int j) 
{
	int s = a[i];
	a[i] = a[j];
	a[j] = s;
}

/*proverka DUS 
*a-massiv poryadka strok
*/
bool proverka(int *a, double **arr, int n1, int n2) {
	bool  kb = 0;                 //flag
	double modul, sum; int id;
	for (int i = 0; i < n1; i++) {
		id = a[i];                 //nomer naidennoi stroki
		modul = fabs(arr[id][i]);           //diagonal element
		sum = sumstr(arr, n2, id) - modul;
		if ((modul < sum) || (modul < EPS)) 
			return false;
		if (modul > sum)       //dlya uslovia aii>sum
			kb = 1;
	}
	if (kb = 1) 
		return true;          //DUS 
	else 
		return false;         //NE DUS

}

/*sum elements stroki*/
double sumstr(double **arr, int n2, int i) {
	double sum = 0;
	for (int j = 0; j < n2 - 1; j++) {
		sum += fabs(arr[i][j]);
	}
	return sum;
}

/*poisk reshenia sistem
**arr-matrica koef
*x - reshenie*/
double * decision(double **arr, double *x, int n1, int n2, int flagD) {
	if (n1 < n2 - 1) {
		return NULL;
	}
	else {
		if (!flagD) {
			if (!nul_diag(arr, n1, n2)) {
				return NULL;
			}
		}
		double *lastx = NULL;        //dlya poiska x-x'
		lastx = memory(lastx, n1);	x = memory(x, n1);
		int k = 0;

		double *maxk10 = NULL;       //dlya proverki sxodimosti
		maxk10 = memory(maxk10, COUNT); // count = 4;                                       

		do {
			k++;
			substitute(arr, x, lastx, n1, n2); //schitaem tekusheie x
			if ((k == 10) && (!flagD)) {
				printf("\n %s %d %s\n\n", "poslednie", COUNT, " max:  ");
				pechat_x(maxk10, COUNT);
				if (!check(maxk10)) {
					printf(" %s \n\n", "ne sxoditsya  ");
					free(maxk10); free(lastx); free(x);
					return NULL;                         //nel'zya reshit dannim metodom
				}
			}
		} while (max_diff(x, lastx, n1, maxk10) >= EPS);
		free(maxk10);
		free(lastx);
		return x;
	}
}


/*schitaem tekushie x - podstanovka v uravnenie*/
void substitute(double **arr, double *x, double *lastx, int n1, int n2) {
	double tmp;

	for (int i = 0; i < n1; i++) {
		tmp = arr[i][n2 - 1];
		for (int j = 0; j < n2 - 1; j++) {
			if (i != j) {
				tmp -= arr[i][j] * x[j];
			}
		}
		lastx[i] = x[i];                //ubrat esli cherez ukazateli
		x[i] = tmp / arr[i][i];
	};

}

/*proverka nulei na diagonali*/
bool nul_diag(double **arr, int n1, int n2) {
	for (int i = 0; i < n1; i++) {
		if (fabs(arr[i][i]) < EPS) {
			int no0 = poisk_n0(arr, i, n1);  //poisk nenulevoi stroki
			if (no0 != -1) {
				perestanovka(arr, i, no0);
			}
			else {
				return 0; //NULL na diagonali
			}
		}
	}
	return 1;//Net nulei na diagonali
}

/*poisk nenulevoi stroki*/
int poisk_n0(double **arr, int i, int n1) {
	for (int k = i + 1; k < n1; k++) {
		if (fabs(arr[k][i]) >= EPS) {
			return k; //nomer stroki bez 0 na diag
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


/*max raznica x-x'*/
double max_diff(double *x, double *lastx, int n1, double *maxk10) {

	double max = 0;
	for (int i = 0; i < n1; i++) {
		if (fabs(x[i] - lastx[i]) > max) {
			max = fabs(x[i] - lastx[i]);
			for (int k = 0; k < COUNT - 1; k++)   //zapominaem poslednie 4 max znachenia
				maxk10[k] = maxk10[k + 1];
			maxk10[COUNT - 1] = max;
		}
	}
	return max;
}



/*proverka monotonnosti dlya poslednix max*/
bool check(double *maxk10) {
	for (int i = 0; i < COUNT - 1; i++) {
		if (maxk10[i] < maxk10[i + 1])
			return 0;         //raznica ne ymenshaetsa
	}
	return 1;              //schitaem dal'she
}

/*matrix koef*/
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

/*chtenie razmera matrix*/
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

/*zapis koef v matrix*/
void  record(FILE *inp, double **arr, int n1, int n2) {
	float tmp;
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			fscanf_s(inp, "%f", &tmp);
			arr[i][j] = (double)tmp;
		}
	}
}



/*pechat matrix koef*/
void pechat(double **arr, int n1, int n2) {

	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			printf("%16e", arr[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/*clean memory*/
void clean(double **arr, int n1) {
	for (int i = 0; i < n1; i++)
		free(arr[i]);
	free(arr);
}

/*videlenie pamiti pod reshenie x*/
double * memory(double *x, int n1) {
	x = (double*)calloc(n1, sizeof(double));
	return x;
}

/*pechat x*/
void pechat_x(double *x, int n1/*, int flag*/) {
	//printf("%s \n\n", "Pechat_x ");
	for (int i = 0; i < n1; i++) {
		printf("%s%d %s %e \n", "X", i + 1, " = ", x[i]);
	}
	printf(" \n");
}