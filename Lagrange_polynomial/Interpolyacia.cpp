#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#define EPS 0.0000001

struct polinom {
	double koef;
	int step;
	polinom * next;
};
struct setka {
	double x;
	double y;
};

polinom * adding(polinom *p1, polinom *p2);                          //p1+p2
polinom * umnog_K(polinom *p1, double k);                            //P1*chislo
polinom * delete_pol(polinom * h);                                   //delete polinom
polinom * umnogenie(polinom *p1, polinom *p2);                       //umnogenie polinomov p3=p1*p2
polinom * multi(int step, double koef);                              //rezultat umnogenia 2x chlenov
double y_val(polinom *h, double x);                                  //vichislenie znachenia v tochke
setka * zap_setka(int &N);                                           //zapolnenie setki znachenii x f(x)   
void print_setka(setka *tabl, polinom *lagran, polinom *newt, int N);//pecat setki znachenii
double fun_y(double x);                                              //vichislenie znachenia y(x)
polinom *lagrange(setka *tabl, int N);                               //polinom Lagrange
polinom * Newton(setka *tabl, int N);                                //polinom Newton
polinom * new_elem(double koef, int step, polinom *prev);            //add new element v spisok
void print_pol(polinom *h);                                          //pecat polinoma
polinom * delete_el(polinom *p1, polinom *prev);                     //delete luboi element spiska
polinom * add_pered(polinom * p1_2, polinom *p1_1, polinom *p2);     //dobavit element pered
polinom * addafter(polinom *prev, polinom *p2);                      //dobavit element posle
polinom * copy(polinom *h);                                          //copy polinom


int main() {
	int N;
	polinom *lagran = NULL;
	polinom *newt = NULL;
	setka *tabl = NULL;

	tabl = zap_setka(N);
	print_setka(tabl, lagran, newt, N);

	printf("\n%s\n", "Lagrange:");
	lagran = lagrange(tabl, N);
	print_pol(lagran);
	
    printf("\n%s\n", "Newton:");
	newt = Newton(tabl, N);
	print_pol(newt);

	print_setka(tabl,lagran,newt, N);
	
    free(tabl);
	delete_pol(lagran);
	delete_pol(newt);
	return 0;
}

/*vichislenie znachenia v tochke*/
double y_val(polinom *h, double x) {
	int k = 0;
	double rez = 0;
	polinom *prev = h;
	while (h != NULL) {
		if ((prev->step) - (h->step) - k <= 1) {  //<1 pri pervom zahode (prev==h)
			rez += h->koef;
			if (h->step != NULL)rez *= x;       //x^0 iskluchaetsya
			prev = h;
			h = h->next;
			k = 0;                     
		}
		else {                        //esli raznica stepen >1
			rez *= x;
			k++;                      //kol-vo umnogenii na x
		}
	}
	if (prev != NULL) {
		for (int i = 1; i < prev->step;i++ ) {         //end
			rez*=x;
		}
	}
	return rez;
}


/*function y(x)*/
double fun_y(double x) {
	return  sin(x);
}

/*zapolnenie setki znachenii x f(x)*/
setka * zap_setka(int &N) {
	float tmp;
	double Beg, End;
	int k = 0;
	double shag;
	printf("%s\n", "vvedite nachalo i konec intervala");
	scanf_s("%f", &tmp);

	Beg = (double)tmp;
	scanf_s("%f", &tmp);
	End = (double)tmp;
	printf("%s\n", "Skol'ko znachenii v intervale?");
	scanf_s("%d", &N);
	if (N == 1)shag = 0;
	else shag = (End - Beg) / (N - 1);
	struct setka *arr = NULL;;
	arr = (setka*)calloc(N, sizeof(setka));
	for (int i = 0; i < N; i++) {              //formirovanie setki
		arr[i].x = Beg + shag * k;
		arr[i].y = fun_y(arr[i].x);
		k++;
	}
	return arr;
}


/*pecat setki znachenii*/
void print_setka(setka *tabl, polinom *lagr, polinom *newt, int N) {
	if ((lagr != NULL) && (newt != NULL)) {
		double tmp = (tabl[0].x + tabl[1].x) / 2;
		printf("\n%16s%16s%16s%16s\n", "x", "f(x)", "Lagrange", "Newton");
		int M = N - 1;
		for (int i = 0; i < M; i++) {
			printf("%16e%s%16e%s%16e%s%16e ", tabl[i].x, " ", tabl[i].y, " ", y_val(lagr, tabl[i].x), " ", y_val(newt, tabl[i].x));
			printf("\n");
			double x = tabl[i].x + tmp;
			printf("%16e%s%16e%s%16e%s%16e ", x, " ", fun_y(x), " ", y_val(lagr, x), " ", y_val(newt, x));
			printf("\n");
		}
		printf("%16e%s%16e%s%16e%s%16e ", tabl[M].x, " ", tabl[M].y, " ", y_val(lagr, tabl[M].x), " ", y_val(newt, tabl[M].x));
		printf("\n");
	}
	else {
	printf("\n%s\n", "SETKA:");
	printf("%16s%16s\n", "x", "f(x)");
	for (int i = 0; i < N; i++) {
		printf("%16e%s%16e", tabl[i].x, " ", tabl[i].y);
		printf("\n");
	}
	printf("\n");
	}
}



/*polinom Newton
New- iscomii polinom
proiz-(x-x0)*(x-x1)*(x-x2)...
fun[0]- koefs polinoma
*/
polinom * Newton(setka *tabl, int N) {
	double *fun = (double*)calloc(N, sizeof(double));  //perepis tabl[i].y v fun[i]
	for (int i = 0; i < N; i++) {
		fun[i] = tabl[i].y;
	}
	polinom *New = NULL;
	if (fabs(fun[0]) > EPS) New = new_elem(fun[0], 0, NULL);    // P(x) New=fun[0]
	polinom *proiz = new_elem(1, 0, NULL);
	polinom *tmp = NULL;
	polinom *cop_proiz = NULL;   //copy of proiz

	for (int j = 1; j < N; j++) {
		for (int i = 0; i < N - j; i++) {
			fun[i] = (fun[i + 1] - fun[i]) / (tabl[i + j].x - tabl[i].x);     //razdelen raznosti
		}
		if (fabs(fun[0]) > EPS) {
			tmp = new_elem(1, 1, NULL);  //1x^1
			if (fabs(tabl[j - 1].x) > EPS) tmp = new_elem(-tabl[j - 1].x, 0, tmp); //esli xi!==0  1x^1-xi
			proiz = umnogenie(proiz, tmp);   // D(x)=(x-x0)*(x-x1)*(x-x2)*...
			tmp = delete_pol(tmp);         //delete (x-xi)
			cop_proiz = copy(proiz);
			tmp = umnog_K(cop_proiz, fun[0]);   //fun[0]*D(x)
			New = adding(New, tmp);        //P(x) =P(x)+fun[0]*D(x)
			tmp = delete_pol(tmp);
		}
	}
	delete_pol(proiz);      //D(x)
	free(fun);              //massiv yi
	return New;
}

/*copy polinom */
polinom * copy(polinom *p1) {
	polinom *nov;
	polinom *head = NULL;     //golova
	polinom *prev = NULL;
	bool k = 1;
	while (p1 != NULL) {
		nov = (polinom*)calloc(1, sizeof(polinom));
		if (k) {
			head = nov;
			k = 0;
		}
		else prev->next = nov;
		nov->koef = p1->koef;
		nov->step = p1->step;
		nov->next = p1->next;
		prev = nov;
		p1 = p1->next;
	}
	return head;
}


/*polinom Lagrange
*/
polinom *lagrange(setka *tabl, int N) {
	polinom *lag = NULL;     //iskomii polinom
	polinom *chisl = NULL;   //chislitel
	polinom *tmp = NULL;
	double znam;           //znamenatel
	double mnog;
	for (int i = 0; i < N; i++) {
		if (fabs(tabl[i].y) > EPS) {    //esli yi==0 slagaemoe=0
			znam = 1;
			chisl = new_elem(1, 0, NULL);   //chisl=1*x^0
			for (int j = 0; j < N; j++) {
				if (j != i) {
					tmp = new_elem(1, 1, NULL);    //1*x^1
					if (fabs(tabl[j].x) > EPS) tmp = new_elem(-tabl[j].x, 0, tmp);	  //esli xi ==0	mnogitel==x	
					chisl = umnogenie(chisl, tmp);    //chisl*(x-xi)
					znam *= (tabl[i].x - tabl[j].x);    //znam*(xi-xj)
					tmp = delete_pol(tmp);
				}
			}
			mnog = tabl[i].y / znam;
			if (fabs(mnog) > EPS) {
				chisl = umnog_K(chisl, mnog);   //  chisl*(yi/znam)
				lag = adding(lag, chisl);  //dobavlyaem novoe slagaemoe
			}
			chisl = delete_pol(chisl);
		}
	}
	return lag;
}


/*add new element v spisok*/
polinom * new_elem(double koef, int step, polinom *prev) {
	polinom *t = (polinom*)calloc(1, sizeof(polinom));
	t->koef = koef;
	t->step = step;
	t->next = NULL;
	if (prev != NULL) {
		prev->next = t;
		return prev;
	}
	return t;
}




/*a+b*/
polinom * adding(polinom *p1, polinom *p2) {
	polinom *t1 = p1; //safe ykazatel golovi
	polinom *prev = NULL;
	while ((p2 != NULL) && (p1 != NULL)) {
		if (p2->step > p1->step) {
			p1 = add_pered(p1, prev, p2);
			if (prev == NULL) t1 = p1;      //save golovy esli max stepen' p2 bol'she max stepeni p1
			prev = p1;
			p1 = p1->next;    //vozvrashaem tekushii element
			p2 = p2->next;    //perehodim k sled v p2
		}
		else if (p2->step == p1->step) {
			p1->koef += p2->koef;
			if (fabs(p1->koef) < EPS) {
				p1 = delete_el(p1, prev);
				if (prev == NULL) t1 = p1;  //save new golovy esli delete old golova
			}
			else {
				prev = p1;
				p1 = p1->next;
			}
			p2 = p2->next;
		}
		else {//if (p2->step < p1->step)
			prev = p1;
			p1 = p1->next;             //idem k sled elementy p1
		}
	}
	while (p2 != NULL) {
		prev = addafter(prev, p2);
		if (t1 == NULL) t1 = prev;  //esli golova 1 polinoma==NULL 
		p2 = p2->next;
	}
	return t1;
}



/*delete luboi element spiska*/
polinom * delete_el(polinom *p1, polinom *prev) {
	if (prev == NULL) {  //delete golova
		polinom *p = p1;
		p1 = p1->next;
		free(p);
		return p1;
	}
	else if (p1->next == NULL) { //delete end
		prev->next = NULL;
		free(p1);
		return NULL;
	}
	else {                      //delte serediny
		prev->next = p1->next;
		polinom *p = p1;
		p1 = p1->next;
		free(p);
		return(p1);
	}
}



/*umnogenie polinoma na chislo*/
polinom * umnog_K(polinom *p1, double k) {
	polinom *t = p1;
	for (; p1 != NULL; p1 = p1->next) {
		p1->koef *= k;
	}
	return t;
}



/*umnogenie polinomov*/
polinom * umnogenie(polinom *p1, polinom *p2) {
	polinom *t2 = p2;  //save golovy p2
	polinom *p3 = NULL; //xodit po p3
	polinom *h = NULL;   //save golova p3
	polinom *prev = NULL; //save pred
	polinom *x = NULL;  //rezultat umnogeniya
	double xk; int xs;
	bool fl = 0;
	for (; p1 != NULL; p1 = p1->next) {
		while (p2 != NULL) {
			//x = multi(p1, p2);
			xs = p1->step + p2->step;
			xk = p1->koef*p2->koef;
			if (!fl) {         //golova
				h = multi(xs, xk);          //soxranyaem golovy
				p3 = h;
				prev = p3;
				p2 = p2->next;
				fl = 1;        //chtobi bol'she ne zahodit v cikl
			}
			else {
				if (xs == p3->step) {
					p3->koef += xk;
					p2 = p2->next;
				}
				else if (xs > p3->step) {
					x = multi(xs, xk);
					x->next = p3;
					prev->next = x;
					p2 = p2->next;
				}
				else if (xs < p3->step) {
					if (p3->next == NULL) {  //esli conec
						x = multi(xs, xk);
						p3->next = x;
						p2 = p2->next;
					}
					prev = p3;
					p3 = p3->next;
				}

			}
		}
		p2 = t2;  //vostanavlivaem golovi polinomov
		p3 = h;
	}
	return h;
}



/*rezultat umnogenia 2x chlenov*/
polinom * multi(int xs, double xk) {
	polinom *x = NULL;
	x = (polinom*)calloc(1, sizeof(polinom));
	x->koef = xk;          //schitaem koef i stepen
	x->step = xs;
	return x;
}



/*delete polinom'*/
polinom * delete_pol(polinom * h) {
	polinom *prev = h;
	while (h != NULL) {
		h = h->next;
		free(prev);
		prev = h;
	}
	return h;
}

/*pecat polinoma*/
void print_pol(polinom *h) {
	printf("%s", "P(n)=");
	while (h != NULL) {
		printf("%e%s%d ", h->koef, "*x^", h->step);
		h = h->next;
	}
	printf("\n");
}


/*dobavit element posle*/
polinom * addafter(polinom *prev, polinom *p2) {
	polinom *t = (polinom*)calloc(1, sizeof(polinom));
	if (prev != NULL) prev->next = t;
	t->koef = p2->koef;
	t->step = p2->step;
	t->next = NULL;
	return t;
}

/*dobavit element pered
*/
polinom *add_pered(polinom * p1_2, polinom *p1_1, polinom *p2) {
	polinom *t = (polinom*)calloc(1, sizeof(polinom));
	t->koef = p2->koef;
	t->step = p2->step;
	t->next = p1_2;
	if (p1_1 != NULL) {    //dlya vsex krome pervogo elementa spiska
		p1_1->next = t;
	}
	return t;
};