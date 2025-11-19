/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/
const double M_PI = 3.141592653589793238;
#include"opt_alg.h"
#include<random>
#include<sstream>
#include<iomanip>
#include<vector>
#include<cmath>
using namespace std;

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

// Funkcja pomocnicza: konwertuje liczbê na string z przecinkiem zamiast kropki
string double_to_csv(double value) {
	ostringstream oss;
	oss << fixed << setprecision(6) << value;
	string result = oss.str();
	// Zamieñ kropkê na przecinek
	size_t pos = result.find('.');
	if (pos != string::npos) {
		result[pos] = ',';
	}
	return result;
}

// Funkcja pomocnicza: konwertuje MA£E liczby z wiêksz¹ precyzj¹ (10 cyfr)
string double_to_csv_high_precision(double value) {
	ostringstream oss;
	oss << fixed << setprecision(10) << value;
	string result = oss.str();
	// Zamieñ kropkê na przecinek
	size_t pos = result.find('.');
	if (pos != string::npos) {
		result[pos] = ',';
	}
	return result;
}

// Funkcja sprawdzaj¹ca typ znalezionego ekstremum
// Dla funkcji testowej f(x) = x^2 - 1: minimum globalne to x=0, f=-1
string check_minimum_type(double x_opt, double f_opt) {
	// Dla f(x) = x^2 - 1: minimum globalne f = -1

	double f_global = -1.0;
	double f_diff = abs(f_opt - f_global);

	// 3 kategorie:

	// 1. Minimum globalne: f bardzo blisko -1 (idealny wynik)
	if (f_diff < 1e-3) {
		return "Minimum globalne";
	}

	// 2. Minimum lokalne: f blisko -1 (dobry wynik)
	if (f_diff < 0.1) {
		return "Minimum lokalne";
	}

	// 3. Brak zbie¿noœci: f daleko od -1 (s³aby wynik)
	return "Brak zbieznosci";
}

int main()
{
	try
	{
		//lab2_wykres_pelny();
		lab2();  // Zamiast lab1()
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok³adnoœæ
	int Nmax = 10000;										// maksymalna liczba wywo³añ funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz górne ograniczenie
		a(2, 1);											// dok³adne rozwi¹zanie optymalne
	solution opt;											// rozwi¹zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Wahadlo
	Nmax = 1000;											// dok³adnoœæ
	epsilon = 1e-2;											// maksymalna liczba wywo³añ funkcji celu
	lb = 0, ub = 5;											// dolne oraz górne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad³a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumieñ do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumieñ
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
}

void lab1()
{



	// Parametry
	double epsilon = 1e-6;
	double gamma = 1e-6;
	int Nmax = 10000;
	double d = 1.0;  // poczatkowa odleglosc dla ekspansji

	// Trzy wspolczynniki ekspansji
	double alphas[3] = { 1.5, 2.0, 3.0 };
	int num_repeats = 10;

	// Generator losowych punktow startowych
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-100.0, 100.0);

	// Plik wynikowy - Tabela 1 (ROZSZERZONA)
	ofstream results_file("wyniki_lab1_tabela1.csv");
	results_file << "Nr;Alpha;x0;a_exp;b_exp;calls_exp;x_fib;f_fib;calls_fib;min_typ_fib;x_lag;f_lag;calls_lag;min_typ_lag" << endl;

	// Tablice do przechowywania srednich
	double sum_fib_x[3] = { 0, 0, 0 };
	double sum_fib_f[3] = { 0, 0, 0 };
	double sum_fib_calls[3] = { 0, 0, 0 };
	double sum_lag_x[3] = { 0, 0, 0 };
	double sum_lag_f[3] = { 0, 0, 0 };
	double sum_lag_calls[3] = { 0, 0, 0 };
	double sum_interval[3] = { 0, 0, 0 };
	double sum_exp_calls[3] = { 0, 0, 0 };

	// Liczniki dla 3 kategorii
	int count_fib_min_global[3] = { 0, 0, 0 };
	int count_fib_min_local[3] = { 0, 0, 0 };
	int count_fib_brak[3] = { 0, 0, 0 };
	int count_lag_min_global[3] = { 0, 0, 0 };
	int count_lag_min_local[3] = { 0, 0, 0 };
	int count_lag_brak[3] = { 0, 0, 0 };

	// NOWE: Szczegó³owe akumulatory dla KA¯DEJ kategorii osobno
	// [alpha_idx][0=min_global, 1=min_local, 2=brak]
	double sum_interval_cat[3][3] = { 0 };
	double sum_exp_calls_cat[3][3] = { 0 };
	double sum_fib_x_cat[3][3] = { 0 };
	double sum_fib_f_cat[3][3] = { 0 };
	double sum_fib_calls_cat[3][3] = { 0 };
	double sum_lag_x_cat[3][3] = { 0 };
	double sum_lag_f_cat[3][3] = { 0 };
	double sum_lag_calls_cat[3][3] = { 0 };

	// Petla po 100 powtorzeniach dla kazdego alpha
	for (int alpha_idx = 0; alpha_idx < 3; alpha_idx++) {
		double alpha = alphas[alpha_idx];
		cout << "Przetwarzanie alpha = " << alpha << endl;

		for (int rep = 0; rep < num_repeats; rep++) {
			if (rep % 10 == 0) {
				cout << "  Postep: " << rep << "/" << num_repeats << endl;
			}
			// Losowy punkt startowy
			double x0 = dis(gen);

			solution::clear_calls();

			// Metoda ekspansji
			double* interval = expansion(ff1T, x0, d, alpha, Nmax);
			double a = interval[0];
			double b = interval[1];
			int calls_exp = solution::f_calls;

			// WALIDACJA: sprawdz czy przedzial jest poprawny
			if (isnan(a) || isnan(b) || isinf(a) || isinf(b) || (b - a) < 1e-10 || a >= b) {
				// Przedzia³ niepoprawny - u¿yj szerokiego przedzia³u wokó³ minimum
				a = -10.0;
				b = 10.0;
			}

			// Metoda Fibonacciego
			solution::clear_calls();
			solution opt_fib = fib(ff1T, a, b, epsilon);
			int calls_fib = solution::f_calls;

			// WALIDACJA: sprawdŸ czy wynik jest poprawny
			if (isnan(m2d(opt_fib.x)) || isnan(m2d(opt_fib.y))) {
				opt_fib.x = 0.0;
				opt_fib.y = -1.0;
			}
			string min_typ_fib = check_minimum_type(m2d(opt_fib.x), m2d(opt_fib.y));

			// Metoda Lagrange'a
			solution::clear_calls();
			solution opt_lag = lag(ff1T, a, b, epsilon, gamma, Nmax);
			int calls_lag = solution::f_calls;

			// WALIDACJA: sprawdŸ czy wynik jest poprawny
			if (isnan(m2d(opt_lag.x)) || isnan(m2d(opt_lag.y))) {
				opt_lag.x = 0.0;
				opt_lag.y = -1.0;
			}
			string min_typ_lag = check_minimum_type(m2d(opt_lag.x), m2d(opt_lag.y));

			// Zapis do pliku Z PRZECINKAMI
			results_file << (alpha_idx * num_repeats + rep + 1) << ";"
				<< double_to_csv(alpha) << ";"
				<< double_to_csv(x0) << ";"
				<< double_to_csv(a) << ";"
				<< double_to_csv(b) << ";"
				<< calls_exp << ";"
				<< double_to_csv(m2d(opt_fib.x)) << ";"
				<< double_to_csv(m2d(opt_fib.y)) << ";"
				<< calls_fib << ";"
				<< min_typ_fib << ";"
				<< double_to_csv(m2d(opt_lag.x)) << ";"
				<< double_to_csv(m2d(opt_lag.y)) << ";"
				<< calls_lag << ";"
				<< min_typ_lag << endl;

			// Akumulacja do srednich
			sum_fib_x[alpha_idx] += m2d(opt_fib.x);
			sum_fib_f[alpha_idx] += m2d(opt_fib.y);
			sum_fib_calls[alpha_idx] += calls_fib;
			sum_lag_x[alpha_idx] += m2d(opt_lag.x);
			sum_lag_f[alpha_idx] += m2d(opt_lag.y);
			sum_lag_calls[alpha_idx] += calls_lag;
			sum_interval[alpha_idx] += (b - a);
			sum_exp_calls[alpha_idx] += calls_exp;

			// Zliczanie typow ekstremum dla Fibonacciego
			int fib_cat_idx = -1;
			if (min_typ_fib == "Minimum globalne") {
				count_fib_min_global[alpha_idx]++;
				fib_cat_idx = 0;
			}
			else if (min_typ_fib == "Minimum lokalne") {
				count_fib_min_local[alpha_idx]++;
				fib_cat_idx = 1;
			}
			else {
				count_fib_brak[alpha_idx]++;
				fib_cat_idx = 2;
			}

			// Zliczanie typow ekstremum dla Lagrange'a
			int lag_cat_idx = -1;
			if (min_typ_lag == "Minimum globalne") {
				count_lag_min_global[alpha_idx]++;
				lag_cat_idx = 0;
			}
			else if (min_typ_lag == "Minimum lokalne") {
				count_lag_min_local[alpha_idx]++;
				lag_cat_idx = 1;
			}
			else {
				count_lag_brak[alpha_idx]++;
				lag_cat_idx = 2;
			}

			// NOWE: Akumulacja szczegó³owa dla kategorii Fibonacciego
			sum_interval_cat[alpha_idx][fib_cat_idx] += (b - a);
			sum_exp_calls_cat[alpha_idx][fib_cat_idx] += calls_exp;
			sum_fib_x_cat[alpha_idx][fib_cat_idx] += m2d(opt_fib.x);
			sum_fib_f_cat[alpha_idx][fib_cat_idx] += m2d(opt_fib.y);
			sum_fib_calls_cat[alpha_idx][fib_cat_idx] += calls_fib;

			// NOWE: Akumulacja szczegó³owa dla kategorii Lagrange'a
			sum_lag_x_cat[alpha_idx][lag_cat_idx] += m2d(opt_lag.x);
			sum_lag_f_cat[alpha_idx][lag_cat_idx] += m2d(opt_lag.y);
			sum_lag_calls_cat[alpha_idx][lag_cat_idx] += calls_lag;

			delete[] interval;
		}
	}

	results_file.close();

	// Tabela 2: Srednie + statystyki minimum/brak zbie¿noœci
	ofstream avg_file("wyniki_lab1_tabela2.csv");
	avg_file << "Alpha;Sredni_przedzial;Srednie_calls_exp;Srednie_x_fib;Srednie_f_fib;Srednie_calls_fib;"
		<< "Fib_min_global;Fib_min_lokalne;Fib_brak;Srednie_x_lag;Srednie_f_lag;Srednie_calls_lag;"
		<< "Lag_min_global;Lag_min_lokalne;Lag_brak" << endl;

	for (int i = 0; i < 3; i++) {
		avg_file << double_to_csv(alphas[i]) << ";"
			<< double_to_csv(sum_interval[i] / num_repeats) << ";"
			<< (int)round(sum_exp_calls[i] / num_repeats) << ";"  // Zaokr¹glone do int
			<< double_to_csv(sum_fib_x[i] / num_repeats) << ";"
			<< double_to_csv(sum_fib_f[i] / num_repeats) << ";"
			<< (int)round(sum_fib_calls[i] / num_repeats) << ";"  // Zaokr¹glone do int
			<< count_fib_min_global[i] << ";"
			<< count_fib_min_local[i] << ";"
			<< count_fib_brak[i] << ";"
			<< double_to_csv(sum_lag_x[i] / num_repeats) << ";"
			<< double_to_csv(sum_lag_f[i] / num_repeats) << ";"
			<< (int)round(sum_lag_calls[i] / num_repeats) << ";"  // Zaokr¹glone do int
			<< count_lag_min_global[i] << ";"
			<< count_lag_min_local[i] << ";"
			<< count_lag_brak[i] << endl;
	}
	avg_file.close();

	// NOWA TABELA 2B: Szczegó³owe œrednie dla KA¯DEJ kategorii osobno (format pionowy)
	ofstream detailed_file("wyniki_lab1_tabela2b_szczegoly.csv");
	detailed_file << "Alpha;b-a;Liczba_wywolan_exp;Rodzaj_minimum;"
		<< "x_fib;y_fib;Liczba_wywolan_fib;Liczba_wystapien_fib;"
		<< "x_lag;y_lag;Liczba_wywolan_lag;Liczba_wystapien_lag" << endl;

	string kategorie[3] = { "Minimum globalne", "Minimum lokalne", "Brak zbieznosci" };

	for (int alpha_idx = 0; alpha_idx < 3; alpha_idx++) {
		for (int cat_idx = 0; cat_idx < 3; cat_idx++) {
			// Liczniki dla tej kategorii
			int count_fib_cat = 0;
			int count_lag_cat = 0;

			if (cat_idx == 0) {
				count_fib_cat = count_fib_min_global[alpha_idx];
				count_lag_cat = count_lag_min_global[alpha_idx];
			}
			else if (cat_idx == 1) {
				count_fib_cat = count_fib_min_local[alpha_idx];
				count_lag_cat = count_lag_min_local[alpha_idx];
			}
			else {
				count_fib_cat = count_fib_brak[alpha_idx];
				count_lag_cat = count_lag_brak[alpha_idx];
			}

			// Œrednie dla tej kategorii (jeœli count > 0)
			double avg_interval = (count_fib_cat > 0) ? sum_interval_cat[alpha_idx][cat_idx] / count_fib_cat : 0.0;
			double avg_exp_calls = (count_fib_cat > 0) ? sum_exp_calls_cat[alpha_idx][cat_idx] / count_fib_cat : 0.0;
			double avg_fib_x = (count_fib_cat > 0) ? sum_fib_x_cat[alpha_idx][cat_idx] / count_fib_cat : 0.0;
			double avg_fib_f = (count_fib_cat > 0) ? sum_fib_f_cat[alpha_idx][cat_idx] / count_fib_cat : 0.0;
			double avg_fib_calls = (count_fib_cat > 0) ? sum_fib_calls_cat[alpha_idx][cat_idx] / count_fib_cat : 0.0;
			double avg_lag_x = (count_lag_cat > 0) ? sum_lag_x_cat[alpha_idx][cat_idx] / count_lag_cat : 0.0;
			double avg_lag_f = (count_lag_cat > 0) ? sum_lag_f_cat[alpha_idx][cat_idx] / count_lag_cat : 0.0;
			double avg_lag_calls = (count_lag_cat > 0) ? sum_lag_calls_cat[alpha_idx][cat_idx] / count_lag_cat : 0.0;

			detailed_file << double_to_csv(alphas[alpha_idx]) << ";"
				<< double_to_csv(avg_interval) << ";"
				<< (int)round(avg_exp_calls) << ";"  // Zaokr¹glone do int
				<< kategorie[cat_idx] << ";"
				<< double_to_csv(avg_fib_x) << ";"
				<< double_to_csv(avg_fib_f) << ";"
				<< (int)round(avg_fib_calls) << ";"  // Zaokr¹glone do int
				<< count_fib_cat << ";"
				<< double_to_csv(avg_lag_x) << ";"
				<< double_to_csv(avg_lag_f) << ";"
				<< (int)round(avg_lag_calls) << ";"  // Zaokr¹glone do int
				<< count_lag_cat << endl;
		}
	}
	detailed_file.close();


	cout << "Optymalizacja bez zawezania przedzialu" << endl;

	double a_full = -100.0;
	double b_full = 100.0;

	// Wektory do zapisu historii dlugosci przedzialu
	vector<double> history_fib_ba;
	vector<double> history_lag_ba;

	//FIBONACCI Z HISTORIA 
	{
		// Generuj ciag Fibonacciego
		const int MAX_FIB = 100;
		long long fib_seq[MAX_FIB];
		fib_seq[0] = 1;
		fib_seq[1] = 1;
		for (int i = 2; i < MAX_FIB; i++) {
			fib_seq[i] = fib_seq[i - 1] + fib_seq[i - 2];
		}

		// Znajdz k
		int k = 0;
		double ratio = (b_full - a_full) / epsilon;
		for (int i = 0; i < MAX_FIB; i++) {
			if (fib_seq[i] > ratio) {
				k = i;
				break;
			}
		}

		// Inicjalizacja
		double a_i = a_full;
		double b_i = b_full;

		// ZAPISZ POCZATKOWY PRZEDZIAL
		history_fib_ba.push_back(b_i - a_i);

		double c_i = b_i - (double)fib_seq[k - 1] / (double)fib_seq[k] * (b_i - a_i);
		double d_i = a_i + b_i - c_i;

		solution C(c_i), D(d_i);
		C.fit_fun(ff1T);
		D.fit_fun(ff1T);

		// Petla glowna
		for (int i = 0; i <= k - 3; i++) {
			if (m2d(C.y) < m2d(D.y)) {
				b_i = d_i;
				d_i = c_i;
				D = C;
				c_i = b_i - (double)fib_seq[k - i - 2] / (double)fib_seq[k - i - 1] * (b_i - a_i);
				C.x = c_i;
				C.fit_fun(ff1T);
			}
			else {
				a_i = c_i;
				c_i = d_i;
				C = D;
				d_i = a_i + b_i - c_i;
				D.x = d_i;
				D.fit_fun(ff1T);
			}

			// ZAPISZ PRZEDZIAL PO KAZDEJ ITERACJI
			history_fib_ba.push_back(b_i - a_i);
		}

		cout << "Fibonacci (bez zawezania): x = " << C.x << ", f = " << C.y << endl;
		cout << "  Liczba iteracji: " << history_fib_ba.size() << endl;
	}

	//LAGRANGE Z HISTORIA
	{
		double a_i = a_full;
		double b_i = b_full;
		double c_i = (a_full + b_full) / 2.0;

		solution A(a_i), B(b_i), C(c_i);
		A.fit_fun(ff1T);
		B.fit_fun(ff1T);
		C.fit_fun(ff1T);

		// ZAPISZ POCZATKOWY
		history_lag_ba.push_back(b_i - a_i);

		double d_i = c_i;
		double d_prev = 0.0;

		int max_iter = 1000;
		for (int iter = 0; iter < max_iter; iter++) {
			// Interpolacja Lagrange'a
			double l = m2d(A.y) * (b_i * b_i - c_i * c_i) +
				m2d(B.y) * (c_i * c_i - a_i * a_i) +
				m2d(C.y) * (a_i * a_i - b_i * b_i);

			double m = m2d(A.y) * (b_i - c_i) +
				m2d(B.y) * (c_i - a_i) +
				m2d(C.y) * (a_i - b_i);

			if (m <= 0) break;

			d_prev = d_i;
			d_i = 0.5 * l / m;

			// Aktualizacja przedzialu
			if (a_i < d_i && d_i < c_i) {
				solution D(d_i);
				D.fit_fun(ff1T);

				if (m2d(D.y) < m2d(C.y)) {
					b_i = c_i;
					B = C;
					c_i = d_i;
					C = D;
				}
				else {
					a_i = d_i;
					A = D;
				}
			}
			else if (c_i < d_i && d_i < b_i) {
				solution D(d_i);
				D.fit_fun(ff1T);

				if (m2d(D.y) < m2d(C.y)) {
					a_i = c_i;
					A = C;
					c_i = d_i;
					C = D;
				}
				else {
					b_i = d_i;
					B = D;
				}
			}
			else {
				break;
			}

			// ZAPISZ PRZEDZIAL
			history_lag_ba.push_back(b_i - a_i);

			// Warunek zakonczenia
			if (b_i - a_i < epsilon || abs(d_i - d_prev) < gamma) {
				break;
			}
		}

		cout << "Lagrange (bez zawezania): x = " << C.x << ", f = " << C.y << endl;
		cout << "  Liczba iteracji: " << history_lag_ba.size() << endl << endl;
	}
	ofstream history_file("wyniki_lab1_wykres.csv");
	history_file << "Nr_iteracji;b-a_Fibonacci;b-a_Lagrange" << endl;

	size_t max_len = (history_fib_ba.size() > history_lag_ba.size()) ?
		history_fib_ba.size() : history_lag_ba.size();

	for (size_t i = 0; i < max_len; i++) {
		history_file << i << ";";

		if (i < history_fib_ba.size()) {
			history_file << double_to_csv(history_fib_ba[i]);
		}
		history_file << ";";

		if (i < history_lag_ba.size()) {
			history_file << double_to_csv(history_lag_ba[i]);
		}
		history_file << endl;
	}
	history_file.close();

	// Optymalizacja
	double a_real = 1.0;
	double b_real = 100.0;

	// Fibonacci
	cout << "Optymalizacja metoda Fibonacciego" << endl;
	solution::clear_calls();
	solution opt_fib_real = fib(ff1R, a_real, b_real, epsilon);
	int calls_fib_real = solution::f_calls;
	cout << "Wynik: DA = " << opt_fib_real.x << " cm^2" << endl;
	cout << "Funkcja celu = " << opt_fib_real.y << endl;
	cout << "Wywolan funkcji: " << calls_fib_real << endl << endl;

	// Lagrange
	cout << "Optymalizacja metoda Lagrange'a" << endl;
	solution::clear_calls();
	solution opt_lag_real = lag(ff1R, a_real, b_real, epsilon, gamma, Nmax);
	int calls_lag_real = solution::f_calls;
	cout << "Wynik: DA = " << opt_lag_real.x << " cm^2" << endl;
	cout << "Funkcja celu = " << opt_lag_real.y << endl;
	cout << "Wywolan funkcji: " << calls_lag_real << endl << endl;

	// Tabela 3: Wyniki dla problemu rzeczywistego
	ofstream real_file("wyniki_lab1_tabela3.csv");
	real_file << "Metoda;DA_opt[cm^2];F_celu;Wywolan_funkcji" << endl;
	real_file << "Fibonacci;" << double_to_csv(m2d(opt_fib_real.x)) << ";"
		<< double_to_csv_high_precision(m2d(opt_fib_real.y)) << ";" << calls_fib_real << endl;
	real_file << "Lagrange;" << double_to_csv(m2d(opt_lag_real.x)) << ";"
		<< double_to_csv_high_precision(m2d(opt_lag_real.y)) << ";" << calls_lag_real << endl;
	real_file.close();


	cout << "Symulacje dla optymalnych wartosci DA" << endl;

	// Warunki poczatkowe (te same dla obu symulacji)
	matrix Y0(4, 1);
	Y0(0) = 5.0;   // VA [m^3]
	Y0(1) = 1.0;   // VB [m^3]
	Y0(2) = 95.0;  // TA [st.C]
	Y0(3) = 20.0;  // TB [st.C]

	// SYMULACJA 1: DA z Fibonacciego
	double DA_fib = m2d(opt_fib_real.x);
	cout << "  Symulacja dla DA (Fibonacci) = " << DA_fib << " cm^2" << endl;
	matrix DA_param_fib(DA_fib);
	matrix* Y_fib = solve_ode(df1, 0, 1.0, 2000.0, Y0, DA_param_fib);

	//SYMULACJA 2: DA z Lagrange'a 
	double DA_lag = m2d(opt_lag_real.x);
	cout << "  Symulacja dla DA (Lagrange) = " << DA_lag << " cm^2" << endl;
	matrix DA_param_lag(DA_lag);
	matrix* Y_lag = solve_ode(df1, 0, 1.0, 2000.0, Y0, DA_param_lag);

	cout << "  Obie symulacje zakonczone." << endl << endl;

	//ZAPIS DO PLIKU CSV
	ofstream sim_file("wyniki_lab1_symulacja.csv");
	sim_file << "t[s];VA_Fib[m^3];VB_Fib[m^3];TA_Fib[C];TB_Fib[C];VA_Lag[m^3];VB_Lag[m^3];TA_Lag[C];TB_Lag[C]" << endl;

	int N_sim = get_size(Y_fib[0])[0];
	for (int i = 0; i < N_sim; i++) {
		sim_file << double_to_csv(Y_fib[0](i)) << ";"
			// Fibonacci
			<< double_to_csv(Y_fib[1](i, 0)) << ";"
			<< double_to_csv(Y_fib[1](i, 1)) << ";"
			<< double_to_csv(Y_fib[1](i, 2)) << ";"
			<< double_to_csv(Y_fib[1](i, 3)) << ";"
			// Lagrange
			<< double_to_csv(Y_lag[1](i, 0)) << ";"
			<< double_to_csv(Y_lag[1](i, 1)) << ";"
			<< double_to_csv(Y_lag[1](i, 2)) << ";"
			<< double_to_csv(Y_lag[1](i, 3)) << endl;
	}
	sim_file.close();

	// ZNAJDZ MAKSYMALNE TEMPERATURY TB
	double TB_max_fib = Y_fib[1](0, 3);
	double TB_max_lag = Y_lag[1](0, 3);

	for (int i = 0; i < N_sim; i++) {
		if (Y_fib[1](i, 3) > TB_max_fib) {
			TB_max_fib = Y_fib[1](i, 3);
		}
		if (Y_lag[1](i, 3) > TB_max_lag) {
			TB_max_lag = Y_lag[1](i, 3);
		}
	}

	cout << "Maksymalna temperatura TB (Fibonacci): " << TB_max_fib << " st.C" << endl;
	cout << "Maksymalna temperatura TB (Lagrange):  " << TB_max_lag << " st.C" << endl;

	// Zwolnij pamiec
	Y_fib[0].~matrix();
	Y_fib[1].~matrix();
	Y_lag[0].~matrix();
	Y_lag[1].~matrix();
}

void lab2()
{
	//PARAMETRY 
	double epsilon = 1e-3;
	int Nmax = 10000;
	double alpha_reduction = 0.5;  // wspó³czynnik zmniejszania kroku (Hooke-Jeeves)
	double alpha_expansion = 2.0;  // wspó³czynnik ekspansji (Rosenbrock)
	double beta_contraction = 0.5; // wspó³czynnik kontrakcji (Rosenbrock)

	// Trzy ró¿ne d³ugoœci kroku startowego
	double step_sizes[3] = { 0.5, 0.1, 0.01 };
	int num_repeats = 100;

	// Generator losowych punktów startowych
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(-1.0, 1.0);

	//CZÊŒÆ A: FUNKCJA TESTOWA
	cout << "LAB 2: funkcja testowa" << endl;

	// Plik wynikowy - Tabela 1
	ofstream results_file("wyniki_lab2_tabela1.csv");
	results_file << "Nr;Dlugosc_kroku;x1_0;x2_0;x1_HJ;x2_HJ;f_HJ;calls_HJ;Min_global_HJ;x1_Rosen;x2_Rosen;f_Rosen;calls_Rosen;Min_global_Rosen" << endl;

	// Tablice do przechowywania sum (dla œrednich)
	int count_HJ_global[3] = { 0, 0, 0 };
	int count_Rosen_global[3] = { 0, 0, 0 };

	double sum_HJ_x1[3] = { 0, 0, 0 };
	double sum_HJ_x2[3] = { 0, 0, 0 };
	double sum_HJ_f[3] = { 0, 0, 0 };
	double sum_HJ_calls[3] = { 0, 0, 0 };

	double sum_Rosen_x1[3] = { 0, 0, 0 };
	double sum_Rosen_x2[3] = { 0, 0, 0 };
	double sum_Rosen_f[3] = { 0, 0, 0 };
	double sum_Rosen_calls[3] = { 0, 0, 0 };

	// Zmienna do przechowania jednego przypadku dla wykresu
	//vector<matrix> history_HJ_x;  // Historia punktów bazowych dla Hooke-Jeeves
	bool saved_for_plot = false;

	// Pêtla po 100 powtórzeniach dla ka¿dego step_size
	for (int step_idx = 0; step_idx < 3; step_idx++) {
		double s = step_sizes[step_idx];
		cout << "Przetwarzanie dla dlugosci kroku = " << s << endl;

		for (int rep = 0; rep < num_repeats; rep++) {
			// Losowy punkt startowy z przedzia³u [-1, 1] x [-1, 1]
			matrix x0(2, 1);
			x0(0) = dis(gen);
			x0(1) = dis(gen);

			//METODA HOOKE'A-JEEVESA
			solution::clear_calls();

			// Dla pierwszego przypadku z pierwsz¹ d³ugoœci¹ kroku - zapisz historiê
			if (step_idx == 0 && rep == 0) {
				// Zmodyfikowana wersja HJ z zapisem historii
				// (u¿yjemy standardowej funkcji, a potem zapiszemy jeden przypadek osobno)
			}

			solution opt_HJ = HJ(ff2T, x0, s, alpha_reduction, epsilon, Nmax);
			int calls_HJ = solution::f_calls;

			// METODA ROSENBROCKA!!!!!!!!!!!
			solution::clear_calls();
			matrix s0(2, 1);
			s0(0) = s;
			s0(1) = s;
			solution opt_Rosen = Rosen(ff2T, x0, s0, alpha_expansion, beta_contraction, epsilon, Nmax);
			int calls_Rosen = solution::f_calls;

			// SprawdŸ czy znaleziono minimum globalne (f ~ 0)
			bool HJ_found_global = (m2d(opt_HJ.y) < 0.001);
			bool Rosen_found_global = (m2d(opt_Rosen.y) < 0.001);


			results_file << (step_idx * num_repeats + rep + 1) << ";"
				<< double_to_csv(s) << ";"
				<< double_to_csv(x0(0)) << ";"
				<< double_to_csv(x0(1)) << ";"
				<< double_to_csv(opt_HJ.x(0)) << ";"
				<< double_to_csv(opt_HJ.x(1)) << ";"
				<< double_to_csv(m2d(opt_HJ.y)) << ";"
				<< calls_HJ << ";"
				<< (HJ_found_global ? "TAK" : "NIE") << ";"
				<< double_to_csv(opt_Rosen.x(0)) << ";"
				<< double_to_csv(opt_Rosen.x(1)) << ";"
				<< double_to_csv(m2d(opt_Rosen.y)) << ";"
				<< calls_Rosen << ";"
				<< (Rosen_found_global ? "TAK" : "NIE") << endl;

			// Akumulacja do œrednich (tylko dla minimum globalnego)
			if (HJ_found_global) {
				count_HJ_global[step_idx]++;
				sum_HJ_x1[step_idx] += opt_HJ.x(0);
				sum_HJ_x2[step_idx] += opt_HJ.x(1);
				sum_HJ_f[step_idx] += m2d(opt_HJ.y);
				sum_HJ_calls[step_idx] += calls_HJ;
			}

			if (Rosen_found_global) {
				count_Rosen_global[step_idx]++;
				sum_Rosen_x1[step_idx] += opt_Rosen.x(0);
				sum_Rosen_x2[step_idx] += opt_Rosen.x(1);
				sum_Rosen_f[step_idx] += m2d(opt_Rosen.y);
				sum_Rosen_calls[step_idx] += calls_Rosen;
			}
		}
	}

	results_file.close();

	//TABELA 2: ŒREDNIE (tylko dla minimum globalnego)
	ofstream avg_file("wyniki_lab2_tabela2.csv");
	avg_file << "Dlugosc_kroku;Metoda;Liczba_min_globalnych;Srednie_x1;Srednie_x2;Srednie_f;Srednie_calls" << endl;

	for (int i = 0; i < 3; i++) {
		// Hooke-Jeeves
		avg_file << double_to_csv(step_sizes[i]) << ";Hooke-Jeeves;"
			<< count_HJ_global[i] << ";";

		if (count_HJ_global[i] > 0) {
			avg_file << double_to_csv(sum_HJ_x1[i] / count_HJ_global[i]) << ";"
				<< double_to_csv(sum_HJ_x2[i] / count_HJ_global[i]) << ";"
				<< double_to_csv(sum_HJ_f[i] / count_HJ_global[i]) << ";"
				<< (int)round(sum_HJ_calls[i] / count_HJ_global[i]);
		}
		else {
			avg_file << ";;;";
		}
		avg_file << endl;

		// Rosenbrock
		avg_file << double_to_csv(step_sizes[i]) << ";Rosenbrock;"
			<< count_Rosen_global[i] << ";";

		if (count_Rosen_global[i] > 0) {
			avg_file << double_to_csv(sum_Rosen_x1[i] / count_Rosen_global[i]) << ";"
				<< double_to_csv(sum_Rosen_x2[i] / count_Rosen_global[i]) << ";"
				<< double_to_csv(sum_Rosen_f[i] / count_Rosen_global[i]) << ";"
				<< (int)round(sum_Rosen_calls[i] / count_Rosen_global[i]);
		}
		else {
			avg_file << ";;;";
		}
		avg_file << endl;
	}
	avg_file.close();

	//WYKRES: Jeden wybrany przypadek dla Hooke-Jeevesa
	// Wykonaj jedn¹ optymalizacjê z zapisem historii punktów bazowych

//WYKRES: Jeden wybrany przypadek dla Hooke-Jeevesa
//WYKRES: Trajektorie dla OBU metod
	cout << endl << "Trajektorie dla obu metod:" << endl;

	double epsilon_wykres = 1e-5;
	double alpha_reduction_wykres = 0.5;
	double alpha_expansion_wykres = 2.0;
	double beta_contraction_wykres = 0.5;

	// Punkt startowy DALEKO od minimum
	matrix x0_plot(2, 1);
	x0_plot(0) = -0.95;
	x0_plot(1) = 0.95;
	double s_plot = 0.08;

	//HOOKE-JEEVES
	cout << "Hooke-Jeeves:" << endl;
	vector<matrix> history_HJ;

	solution XB_HJ(x0_plot);
	XB_HJ.fit_fun(ff2T);
	history_HJ.push_back(XB_HJ.x);

	double s_HJ = s_plot;
	int max_iter_HJ = 3000;
	int iter_HJ = 0;

	while (s_HJ >= epsilon_wykres && iter_HJ < max_iter_HJ) {
		iter_HJ++;

		if (iter_HJ % 200 == 0) {
			cout << "  HJ: " << iter_HJ << " iter, " << history_HJ.size() << " punktow" << endl;
		}

		solution X = HJ_trial(ff2T, XB_HJ, s_HJ);

		if (X.y < XB_HJ.y) {
			int pattern_moves = 0;

			while (pattern_moves < 100) {
				pattern_moves++;

				solution XB_prev = XB_HJ;
				XB_HJ = X;
				history_HJ.push_back(XB_HJ.x);

				solution X_pattern(2.0 * XB_HJ.x - XB_prev.x);
				X_pattern.fit_fun(ff2T);
				X = HJ_trial(ff2T, X_pattern, s_HJ);

				if (X.y >= XB_HJ.y) {
					X = XB_HJ;
					break;
				}
			}
		}
		else {
			s_HJ = alpha_reduction_wykres * s_HJ;
		}
	}

	cout << "HJ: " << history_HJ.size() << " punktow" << endl;

	// ROSENBROCK
	cout << "Rosenbrock:" << endl;
	vector<matrix> history_Rosen;

	matrix s_Rosen(2, 1);
	s_Rosen(0) = s_plot;
	s_Rosen(1) = s_plot;

	matrix D_wykres = ident_mat(2);

	solution XB_Rosen(x0_plot);
	XB_Rosen.fit_fun(ff2T);
	history_Rosen.push_back(XB_Rosen.x);

	int iter_Rosen = 0;

	while (iter_Rosen < 3000) {
		iter_Rosen++;

		if (iter_Rosen % 200 == 0) {
			cout << "  Rosen: " << iter_Rosen << " iter, " << history_Rosen.size() << " punktow" << endl;
		}

		for (int j = 0; j < 2; j++) {
			matrix d_j = D_wykres[j];
			solution X_trial(XB_Rosen.x + s_Rosen(j) * d_j);
			X_trial.fit_fun(ff2T);

			if (X_trial.y < XB_Rosen.y) {
				XB_Rosen = X_trial;
				history_Rosen.push_back(XB_Rosen.x);
				s_Rosen(j) = alpha_expansion_wykres * s_Rosen(j);
			}
			else {
				s_Rosen(j) = -beta_contraction_wykres * s_Rosen(j);
			}
		}

		double max_s_wykres = max(abs(s_Rosen(0)), abs(s_Rosen(1)));
		if (max_s_wykres < epsilon_wykres) break;
	}

	cout << "Rosenbrock: " << history_Rosen.size() << " punktow" << endl;
	ofstream plot_file("wyniki_lab2_wykres_pelny.csv");
	plot_file << "Iteracja;x1_HJ;x2_HJ;x1_Rosen;x2_Rosen" << endl;

	size_t max_len = max(history_HJ.size(), history_Rosen.size());

	for (size_t i = 0; i < max_len; i++) {
		plot_file << i << ";";

		if (i < history_HJ.size()) {
			plot_file << double_to_csv(history_HJ[i](0)) << ";"
				<< double_to_csv(history_HJ[i](1));
		}
		else {
			plot_file << ";";
		}
		plot_file << ";";

		if (i < history_Rosen.size()) {
			plot_file << double_to_csv(history_Rosen[i](0)) << ";"
				<< double_to_csv(history_Rosen[i](1));
		}
		plot_file << endl;
	}

	plot_file.close();

	//CZÊŒÆ B: PROBLEM RZECZYWISTY (ROBOT)
	cout << endl << "LAB 2: problem rzeczywisty" << endl;

	// Punkt startowy dla problemu robota: k1, k2 z  [0, 20]
	random_device rd_real;
	mt19937 gen_real(rd_real());
	uniform_real_distribution<> dis_real(0.0, 20.0);

	matrix x0_real(2, 1);
	x0_real(0) = dis_real(gen_real);  // k1
	x0_real(1) = dis_real(gen_real);  // k2

	double s_real = 1.0;  // d³ugoœæ kroku dla problemu rzeczywistego

	cout << "Punkt startowy: k1 = " << x0_real(0) << ", k2 = " << x0_real(1) << endl;

	// HOOKE-JEEVES
	cout << endl << "Optymalizacja metoda Hooke'a-Jeevesa:" << endl;
	solution::clear_calls();
	solution opt_HJ_real = HJ(ff2R, x0_real, s_real, alpha_reduction, epsilon, Nmax);
	int calls_HJ_real = solution::f_calls;


	cout << "Wynik HJ: k1 = " << opt_HJ_real.x(0) << ", k2 = " << opt_HJ_real.x(1) << endl;
	cout << "Q = " << opt_HJ_real.y << endl;
	cout << "Wywolan funkcji: " << calls_HJ_real << endl;

	//ROSENBROCK 
	cout << endl << "Optymalizacja metoda Rosenbrocka:" << endl;
	solution::clear_calls();
	matrix s0_real(2, 1);
	s0_real(0) = s_real;
	s0_real(1) = s_real;
	solution opt_Rosen_real = Rosen(ff2R, x0_real, s0_real, alpha_expansion, beta_contraction, epsilon, Nmax);
	int calls_Rosen_real = solution::f_calls;

	cout << "Wynik Rosenbrock: k1 = " << opt_Rosen_real.x(0) << ", k2 = " << opt_Rosen_real.x(1) << endl;
	cout << "Q = " << opt_Rosen_real.y << endl;
	cout << "Wywolan funkcji: " << calls_Rosen_real << endl << endl;

	// Weryfikacja dla k1=5, k2=5
	cout << "Weryfikacja dla k1=5, k2=5:" << endl;
	matrix x_verify(2, 1);
	x_verify(0) = 5.0;
	x_verify(1) = 5.0;
	solution::clear_calls();
	solution verify_sol(x_verify);
	verify_sol.fit_fun(ff2R);
	cout << "Q(5, 5) = " << verify_sol.y << " (oczekiwane: ~775.229)" << endl << endl;

	//TABELA 3: Wyniki dla problemu rzeczywistego
	ofstream real_file("wyniki_lab2_tabela3.csv");
	real_file << "Metoda;k1_0;k2_0;k1_opt;k2_opt;Q;Wywolan_funkcji" << endl;
	real_file << "Hooke-Jeeves;"
		<< double_to_csv(x0_real(0)) << ";"
		<< double_to_csv(x0_real(1)) << ";"
		<< double_to_csv(opt_HJ_real.x(0)) << ";"
		<< double_to_csv(opt_HJ_real.x(1)) << ";"
		<< double_to_csv(m2d(opt_HJ_real.y)) << ";"
		<< calls_HJ_real << endl;
	real_file << "Rosenbrock;"
		<< double_to_csv(x0_real(0)) << ";"
		<< double_to_csv(x0_real(1)) << ";"
		<< double_to_csv(opt_Rosen_real.x(0)) << ";"
		<< double_to_csv(opt_Rosen_real.x(1)) << ";"
		<< double_to_csv(m2d(opt_Rosen_real.y)) << ";"
		<< calls_Rosen_real << endl;
	real_file.close();

	//SYMULACJE
	cout << "start" << endl;

	// Warunki pocz¹tkowe: alpha(0) = 0, omega(0) = 0
	matrix Y0(2, 1);
	Y0(0) = 0.0;  // alpha
	Y0(1) = 0.0;  // omega

	//SYMULACJA 1: Wyniki z Hooke-Jeeves
	matrix k_HJ(2, 1);
	k_HJ(0) = opt_HJ_real.x(0);
	k_HJ(1) = opt_HJ_real.x(1);
	matrix* Y_HJ = solve_ode(df2, 0, 0.1, 100.0, Y0, k_HJ);

	// SYMULACJA 2: Wyniki z Rosenbrock
	matrix k_Rosen(2, 1);
	k_Rosen(0) = opt_Rosen_real.x(0);
	k_Rosen(1) = opt_Rosen_real.x(1);
	matrix* Y_Rosen = solve_ode(df2, 0, 0.1, 100.0, Y0, k_Rosen);

	cout << "stop" << endl << endl;

	ofstream sim_file("wyniki_lab2_symulacja.csv");
	sim_file << "t[s];alpha_HJ[rad];omega_HJ[rad/s];alpha_Rosen[rad];omega_Rosen[rad/s]" << endl;

	int N_sim = get_size(Y_HJ[0])[0];
	for (int i = 0; i < N_sim; i++) {
		sim_file << double_to_csv(Y_HJ[0](i)) << ";"
			<< double_to_csv(Y_HJ[1](i, 0)) << ";"
			<< double_to_csv(Y_HJ[1](i, 1)) << ";"
			<< double_to_csv(Y_Rosen[1](i, 0)) << ";"
			<< double_to_csv(Y_Rosen[1](i, 1)) << endl;
	}
	sim_file.close();

	//ANALIZA WYNIKÓW SYMULACJI 
	double alpha_final_HJ = Y_HJ[1](N_sim - 1, 0);
	double omega_final_HJ = Y_HJ[1](N_sim - 1, 1);
	double alpha_final_Rosen = Y_Rosen[1](N_sim - 1, 0);
	double omega_final_Rosen = Y_Rosen[1](N_sim - 1, 1);

	const double M_PI = 3.141592653589793238;
	cout << "Wyniki koncowe symulacji:" << endl;
	cout << "Hooke-Jeeves:" << endl;
	cout << "  alpha(100s) = " << alpha_final_HJ << " rad (cel: " << M_PI << " rad)" << endl;
	cout << "  omega(100s) = " << omega_final_HJ << " rad/s (cel: 0 rad/s)" << endl;
	cout << "  Blad polozenia: " << abs(alpha_final_HJ - M_PI) << " rad" << endl;
	cout << "  Blad predkosci: " << abs(omega_final_HJ) << " rad/s" << endl << endl;

	cout << "Rosenbrock:" << endl;
	cout << "  alpha(100s) = " << alpha_final_Rosen << " rad (cel: " << M_PI << " rad)" << endl;
	cout << "  omega(100s) = " << omega_final_Rosen << " rad/s (cel: 0 rad/s)" << endl;
	cout << "  Blad polozenia: " << abs(alpha_final_Rosen - M_PI) << " rad" << endl;
	cout << "  Blad predkosci: " << abs(omega_final_Rosen) << " rad/s" << endl << endl;

	// Zwolnij pamiêæ
	Y_HJ[0].~matrix();
	Y_HJ[1].~matrix();
	Y_Rosen[0].~matrix();
	Y_Rosen[1].~matrix();
}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}