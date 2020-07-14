#include <stdio.h>  /* per scrivere l'output eventualmente su un file*/
#include <stdlib.h> /* per la funzione rand */
// #include <time.h> /* per inizializzare la funzione rand */
#include <R.h> /* per comunicare con R */
#include <Rinternals.h>
#include <math.h>

/* la matrice zeta di un ordinamento lineare e' una matrice ''triangolare'' */
void linzeta(int *lin, int *n, int *result)
{
	int nval = *n;
	int pos;

	for (int i = 0; i < nval; i++)
	{
		for (int j = 0; j < nval; j++)
		{
			pos = nval * (lin[i] - 1) + (lin[j] - 1);
			if (i <= j)
				result[pos] = 1;
			else
				result[pos] = 0;
		}
	}
}

// funzione che genera nuove estensioni lineari... se necessario
// se lo fa restituisce 1 altrimenti 0
int new_linext(
	int n,
	int *linext,
	int *z)
{
	int c = (int)floor(unif_rand() * RAND_MAX) % 2; // scelgo se generare una nuova estensione lineare
	if (c == 0)
		return 0; // inutile fare altri conti e riciclare i precedenti

	int p = (int)floor(unif_rand() * RAND_MAX) % (n - 1); // posizione da scambiare con la successiva

	int pos = n * linext[p] + linext[p + 1]; // determino la posizione in z
	if (z[pos] == 0)
	{ // verifico se puo' avvenire lo scambio
		// faccio lo scambio
		int tmp = linext[p];
		linext[p] = linext[p + 1];
		linext[p + 1] = tmp;
		return 1; // lo scambio e' avvenuto
	}

	// sebbene si abbia provato a generare una nuova estensione lineare,
	// lo scambio ha generato un ordinamento completo non compatibile col poset
	return 0;
}

// funzione principale, che esegue il loop e lancia le altre funzioni necessarie
// a generare le estensioni lineari e ad eseguire i conti su di esse
void bd(
	int *linext,	   // !!! estensione lineare come ranghi da 0 a n-1
	int *n_scal,	   // !!! numero di elementi, _scal indica che e' uno
					   // scalare e permette di dare un nome piu' semplice e
					   // chiaro alla variabile che viene effettivamente usata
	int *nit_scal,	 // numero di iterazioni
	int *z,			   // !!! matrice di copertura traslata
	int *rankfreq,	 // matrice profili x ranghi che conta le rispettive
					   // frequenze assolute
	int *threshold,	// soglia indicata come vettore di 0 e 1 se incluso o no
	int *thrfreq,	  // frequenza con la quale un elemento e' threshold
					   // nell'estensione lineare
	int *loweqthr,	 // conta il numero di volte che un elemento e' sotto o
					   // uguale alla soglia
	double *weights,   // pesi da attribuire a ciascun profilo
	double *distances, // distanze tra i profili qualora siano confrontabili
					   // che rappresentano i pesi degli edges che li collegano
	// necessarie per valutare i gap pesati
	double *cumdist, // vettore delle distanze cumulate dal basso
					 // nell'estensione lineare, inizializzate a 0,
	// strumentale al calcolo del gap
	double *gapAP,			  // gap assoluto di poverta'
	double *gapRP,			  // gap relativo di poverta'
	double *gapAR,			  // gap assoluto di ricchezza
	double *gapRR			  // gap relativo di ricchezza
)
{
	// inizializzazione della random seed
	// srand(seed);
	GetRNGstate();

	// importo gli scalari
	int n = *n_scal;
	int nit = *nit_scal;

	// variabili di supporto
	int new = 0;
	int pos;
	int linthr = 0; // threshold nell'estensione lineare

	for (int w = 0; w < nit; w++)
	{ // LOOP principale

		new = new_linext(n, linext, z);

		// qua ci sono tutte le operazioni necessarie ogni qual volta che si
		// cambia l'estensione lineare
		if (w == 0 || new == 1)
		{

			linthr = 0;
			for (int i = 0; i < n; i++)
			{
				// calcolo le "distanze cumulate"
				if (i > 0)
				{
					pos = n * linext[i - 1] + linext[i];
					cumdist[linext[i]] = cumdist[linext[i - 1]] + distances[pos];
				}
				else
					cumdist[linext[0]] = 0;
				// trovo la soglia dell'estensione lineare
				if (threshold[linext[i]])
				{
					linthr = i;
				}
			}
		}

		// frequenza assoluta di soglia
		thrfreq[linext[linthr]]++;

		// conti che vanno a fatti a prescindere dal cambio di estensione
		// lineare
		for (int i = 0; i < n; i++)
		{
			// conto le frequenze dei ranghi
			pos = n * linext[i] + i;
			rankfreq[pos]++;
			if (i <= linthr)
			{
				// frequenze assolute delle volte che un elemento e' minore o
				// uguale alla soglia
				loweqthr[linext[i]]++;
				// calcolo il gap assoluto di poverta'
				// il +1 serve a spostarsi al primo profilo non di poverta'
				gapAP[linext[i]] +=
					cumdist[linext[linthr + 1]] - cumdist[linext[i]];
				// calcolo il gap relativo di poverta' dividendo per il massimo
				// valore possibile
				gapRP[linext[i]] +=
					1 - cumdist[linext[i]] / cumdist[linext[linthr + 1]];
			}
			if (i > linthr)
			{
				// calcolo il gap assoluto di ricchezza
				// il +1 non serve perche' il primo profilo di poverta' e' la
				// soglia stessa
				gapAR[linext[i]] +=
					cumdist[linext[i]] - cumdist[linext[linthr]];
				// calcolo il gap relativo di ricchezza dividendo per il massimo
				// valore possibile
				gapRR[linext[i]] +=
					(cumdist[linext[i]] - cumdist[linext[linthr]]) / (cumdist[linext[n - 1]] - cumdist[linext[linthr]]);
			}
		}
	}
	PutRNGstate();
}

// funzione principale, che esegue il loop e lancia le altre funzioni necessarie
// a generare le estensioni lineari e ad eseguire i conti su di esse
// (solo l'identification function)
void bd_simp(
	int *linext, // !!! estensione lineare come ranghi da 0 a n-1
	int *n_scal, // !!! numero di elementi, _scal indica che e' uno
	// scalare e permette di dare un nome piu' semplice e
	// chiaro alla variabile che viene effettivamente usata
	int *nit_scal, // numero di iterazioni
	int *z,		   // !!! matrice di copertura traslata
	int *rankfreq, // matrice profili x ranghi che conta le rispettive
	// frequenze assolute
	int *threshold, // soglia indicata come vettore di 0 e 1 se incluso o no
	int *thrfreq,   // frequenza con la quale un elemento e' threshold
	// nell'estensione lineare
	int *loweqthr, // conta il numero di volte che un elemento e' sotto o
	// uguale alla soglia
	double *weights //, // pesi da attribuire a ciascun profilo
					//	double *distances, // distanze tra i profili qualora siano confrontabili
					//	                   // che rappresentano i pesi degli edges che li collegano
					//					   // necessarie per valutare i gap pesati
					//  double *cumdist,   // vettore delle distanze cumulate dal basso
					//	                   // nell'estensione lineare, inizializzate a 0,
					//					   // strumentale al calcolo del gap
					//  double *gapAP,     // gap assoluto di poverta'
					//	double *gapRP,     // gap relativo di poverta'
					//	double *gapAR,     // gap assoluto di ricchezza
					//	double *gapRR,     // gap relativo di ricchezza
					//	double *polarization_scal // indice di polarizzazione
)
{
	// inizializzazione della random seed
	// srand(time(NULL));
	// srand(seed);
	GetRNGstate();

	// importo gli scalari
	int n = *n_scal;
	int nit = *nit_scal;

	// variabili di supporto
	int new = 0;
	int pos;
	int linthr = 0; // threshold nell'estensione lineare
					//	double tmppolar = 0;

	for (int w = 0; w < nit; w++)
	{ // LOOP principale

		new = new_linext(n, linext, z);

		// qua ci sono tutte le operazioni necessarie ogni qual volta che si
		// cambia l'estensione lineare
		if (w == 0 || new == 1)
		{

			linthr = 0;
			for (int i = 0; i < n; i++)
			{
				// calcolo le "distanze cumulate"
				//				if (i > 0) {
				//				    pos = n * linext[i-1] + linext[i];
				//				    cumdist[linext[i]] = cumdist[linext[i-1]] + distances[pos];
				//				}
				//				else
				//				    cumdist[linext[0]] = 0;
				// trovo la soglia dell'estensione lineare
				if (threshold[linext[i]])
				{
					linthr = i;
				}
			}

			// calcolo l'indice di polarizzazione
			//			tmppolar = 0;
			//			if (polarization_scal[0] >= 0)
			//			    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
			//			        tmppolar += weights[linext[j]]*weights[linext[i]]*(j - i);
			//			    }
		}

		// aggiorno la media della polarizzazione, media e non somma per non
		// avere numeri troppo elevati
		//		if (polarization_scal[0] >= 0)
		//		    polarization_scal[0] = (polarization_scal[0] * w + tmppolar)/(w+1);

		// frequenza assoluta di soglia
		thrfreq[linext[linthr]]++;

		// conti che vanno a fatti a prescindere dal cambio di estensione
		// lineare
		for (int i = 0; i < n; i++)
		{
			// conto le frequenze dei ranghi
			pos = n * linext[i] + i;
			rankfreq[pos]++;
			if (i <= linthr)
			{
				// frequenze assolute delle volte che un elemento e' minore o
				// uguale alla soglia
				loweqthr[linext[i]]++;
				// calcolo il gap assoluto di poverta'
				// il +1 serve a spostarsi al primo profilo non di poverta'
				//				gapAP[linext[i]] +=
				//			        cumdist[linext[linthr +1]] - cumdist[linext[i]];
				// calcolo il gap relativo di poverta' dividendo per il massimo
				// valore possibile
				//			    gapRP[linext[i]] +=
				//			        1 - cumdist[linext[i]]/cumdist[linext[linthr +1]];
			}
			//			if (i > linthr) {
			// calcolo il gap assoluto di ricchezza
			// il +1 non serve perche' il primo profilo di poverta' e' la
			// soglia stessa
			//				gapAR[linext[i]] +=
			//			        cumdist[linext[i]] - cumdist[linext[linthr]];
			// calcolo il gap relativo di ricchezza dividendo per il massimo
			// valore possibile
			//			    gapRR[linext[i]] +=
			//			        (cumdist[linext[i]] - cumdist[linext[linthr]])
			//					/(cumdist[linext[n-1]] - cumdist[linext[linthr]]);
			//			}
		}
	}
	PutRNGstate();
}

/*****************************/
/* definisco tutte le funzioni per il calcolo degli indici di polarizzazione disponibili */

int findmedian(
	int *linext,
	double *freqrel,
	int n)
{
	double cumfreqctr = 0;
	for (int i = 0; i < n; i++)
	{
		cumfreqctr += freqrel[linext[i]];
		if (cumfreqctr >= 0.5)
		{
			return i;
		}
	}

	return 0;
}

double abulnaga(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double tmppolar = 0;
	double cumfreqctr = 0;
	int median = findmedian(linext, freqrel, n);

	for (int i = 0; i < median; i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar += pow(cumfreqctr, *alpha);
	}

	for (int i = median; i < n; i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar -= pow(cumfreqctr, *beta);
	}

	tmppolar = (tmppolar + n - median) / (median * pow(0.5, *alpha) - 1 - (n - median - 1) * pow(0.5, *beta) + n - median);

	return tmppolar;
}

double apouey(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double cumfreqctr = 0;
	double tmppolar = 0;
	for (int i = 0; i < (n - 1); i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar += pow(fabs(cumfreqctr - 0.5), *alpha);
	}

	return tmppolar;
}

double blairLacyOne(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double cumfreqctr = 0;
	double tmppolar = 0;

	for (int i = 0; i < (n - 1); i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar += cumfreqctr * (1 - cumfreqctr);
	}

	return tmppolar;
}
double blairLacyTwo(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double tmppolar = 0;
	double cumfreqctr = 0;

	for (int i = 0; i < (n - 1); i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar += pow(cumfreqctr - 0.5, 2);
	}

	return sqrt(tmppolar);
}

double kobusMilos(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double tmppolar = 0;
	double tmppolar_a = 0;
	double tmppolar_b = 0;
	double cumfreqctr = 0;
	int median = findmedian(linext, freqrel, n);

	for (int i = 0; i < median; i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar_a += cumfreqctr;
	}

	for (int i = median; i < n; i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar_b += cumfreqctr;
	}

	tmppolar = (*alpha * tmppolar_a - *beta * tmppolar_b + *beta * (n - median)) / (*alpha * median + *beta * (n - median - 1));

	return tmppolar;
}

double kvalseth(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double tmppolar = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			tmppolar += freqrel[linext[i]] * freqrel[linext[j]] * abs(i - j);
		}
	}

	tmppolar *= 2.0 / (n - 1);
	return sqrt(1 - tmppolar);
}

double leik(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double tmppolar = 0;
	double cumfreqctr = 0;
	for (int i = 0; i < (n - 1); i++)
	{
		cumfreqctr += freqrel[linext[i]];
		if (cumfreqctr <= 0.5)
			tmppolar += cumfreqctr;
		else
			tmppolar += 1 - cumfreqctr;
	}

	return tmppolar;
}

double reardonOne(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double tmppolar = 0;
	double cumfreqctr = 0;

	for (int i = 0; i < (n - 1); i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar += cumfreqctr * log2(cumfreqctr) + (1 - cumfreqctr) * log2(1 - cumfreqctr);
	}
	return tmppolar;
}

double reardonThree(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double cumfreqctr = 0;
	double tmppolar = 0;
	for (int i = 0; i < (n - 1); i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar += sqrt(cumfreqctr * (1 - cumfreqctr));
	}
	return tmppolar;
}

double reardonFour(
	int *linext,
	double *freqrel,
	double *alpha,
	double *beta,
	int n)
{
	double cumfreqctr = 0;
	double tmppolar = 0;

	for (int i = 0; i < (n - 1); i++)
	{
		cumfreqctr += freqrel[linext[i]];
		tmppolar += fabs(2 * cumfreqctr - 1);
	}

	return tmppolar;
}

/* Functions to compute quantiles, Raatikainen (1987) */

int findK(
	double *q,
	double x)
{
	if (x < q[0])
	{
		q[0] = x;
		return 0;
	}
	for (int i = 1; i < 2 * 12; i++)
	{
		if (q[i] <= x && x < q[i + 1])
			return i;
	}
	q[24] = x;
	return 23;
}

void adjustQuantilePos(
	double *q,
	int *n,
	double *d)
{
	double offset;
	int dp;
	int dm;
	double qp;
	double qm;
	double qt;

	for (int i = 1; i < 2 * 12; i++)
	{
		offset = d[i] - n[i];
		dp = n[i + 1] - n[i];
		dm = n[i - 1] - n[i];
		qp = (q[i + 1] - q[i]) / dp;
		qm = (q[i - 1] - q[i]) / dm;

		if (offset >= 1 && dp > 1)
		{
			qt = q[i] + ((1 - dm) * qp + (dp - 1) * qm) / (dp - dm);
			if (q[i - 1] < qt && qt < q[i + 1])
				q[i] = qt;
			else
				q[i] += qp;
			n[i] += 1;
		}

		else if (offset <= -1 && dm < -1)
		{
			qt = q[i] - ((1 + dp) * qm - (dm + 1) * qp) / (dp - dm);
			if (q[i - 1] < qt && qt < q[i + 1])
				q[i] = qt;
			else
				q[i] -= qm;
			n[i] -= 1;
		}
	}
}

void computeQuantiles(
	double *q,
	int *n,
	double *f,
	double *d,
	double x)
{

	int k = findK(q, x);
	for (int i = 0; i < k + 1; i++)
	{
		d[i] += f[i];
	}
	for (int i = k + 1; i < 2 * 11 + 3; i++)
	{
		d[i] += f[i];
		n[i] += 1;
	}

	adjustQuantilePos(q, n, d);
}

/* standard double comparing function to use with qsort */
static int compare(const void *a, const void *b)
{
	if (*(double *)a > *(double *)b)
		return -1;
	else if (*(double *)a < *(double *)b)
		return 1;
	else
		return 0;
}

/*****************************/
/* funzione principale per il calcolo delle misure di polarizzazione */

void bd_polarization(
	int *measure,
	double *alpha,
	double *beta,
	int *linext,
	int *n_scal,
	int *nit_scal,
	int *z,
	double *relfreqs,
	double *polarization_mean,
	double *polarization_variance,
	double *quantiles_q,
	int *quantiles_n,
	double *quantiles_f,
	double *quantiles_d,
	double *polar_min,
	double *polar_max,
	int *firstit)
{
	GetRNGstate();
	int n = *n_scal;
	int nit = *nit_scal;
	int m = *measure;
	int new = 0;
	double tmppolar = 0;
	double prev_mean = 0;
	double delta = 0;
	int countelements = 0;

	double (*functionPtr[10])(int *le, double *f, double *a, double *b, int n) = {
		abulnaga, apouey, blairLacyOne, blairLacyTwo, kobusMilos, kvalseth, leik, reardonOne, reardonThree, reardonFour};

	for (int w = 0; w < nit; w++)
	{
		new = new_linext(n, linext, z);

		if (w == 0 || new == 1)
		{
			prev_mean = polarization_mean[0];
			tmppolar = (*functionPtr[m])(linext, relfreqs, alpha, beta, n);
			delta = tmppolar - prev_mean;
			polarization_mean[0] += delta / (w + 1);
			polarization_variance[0] += delta * (tmppolar - polarization_mean[0]);

			/* alla prima run devo tenere in memoria le prime 25 osservazioni generate e ordinarle per calcolare i quantili */
			if (firstit[0])
			{
				if (countelements <= 24)
					quantiles_q[countelements] = tmppolar;
					countelements += 1;

				if (countelements == 24)
				{
					qsort(quantiles_q, 25, sizeof(double), compare);
					firstit[0] = 0;
				}
			}

			else
				computeQuantiles(quantiles_q, quantiles_n, quantiles_f, quantiles_d, tmppolar);

			if (tmppolar < polar_min[0])
				polar_min[0] = tmppolar;
			else if (tmppolar > polar_max[0])
				polar_max[0] = tmppolar;
		}
	}

	PutRNGstate();
}
