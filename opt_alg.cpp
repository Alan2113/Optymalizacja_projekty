#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejściowe:
	// ff - wskaźnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zakłądana dokładność rozwiązania
	// Nmax - maksymalna liczba wywołań funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosując rozkład jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwiązanie do przedziału [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartość funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwiązanie z zadaną dokładnością
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywołań funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

// ========== METODA EKSPANSJI ==========
double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		solution X0(x0);
		solution X1(x0 + d);
		int i = 0;
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);

		if (X1.y == X0.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}
		if (X1.y > X0.y) {
			d = -d;
			X1.x = x0 + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y) {
				p[0] = min(m2d(X1.x), x0 - d);
				p[1] = max(m2d(X1.x), x0 - d);
				return p;
			}
		}
		solution X_prev = X0;
		solution X_curr = X1;

		while (true) {
			if (solution::f_calls > Nmax) {
				throw string("Expansion method: exceeded max function calls");
			}

			i++;
			solution X_next(x0 + pow(alpha, i) * d);
			X_next.fit_fun(ff, ud1, ud2);
			if (X_curr.y <= X_next.y) {
				break;
			}

			X_prev = X_curr;
			X_curr = X_next;
		}

		if (d > 0) {
			p[0] = m2d(X_prev.x);
			p[1] = m2d(X_curr.x) + abs(d); // x_next
		}
		else {
			p[0] = m2d(X_curr.x) + abs(d); // x_next
			p[1] = m2d(X_prev.x);
		}
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

// ========== METODA FIBONACCIEGO ==========
solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		const int MAX_FIB = 100;
		long long* fib_seq = new long long[MAX_FIB];
		fib_seq[0] = 1;
		fib_seq[1] = 1;
		for (int i = 2; i < MAX_FIB; i++) {
			fib_seq[i] = fib_seq[i - 1] + fib_seq[i - 2];
		}

		int k = 0;
		double ratio = (b - a) / epsilon;
		for (int i = 0; i < MAX_FIB; i++) {
			if (fib_seq[i] > ratio) {
				k = i;
				break;
			}
		}
		double a_i = a;
		double b_i = b;
		double c_i = b_i - (double)fib_seq[k - 1] / (double)fib_seq[k] * (b_i - a_i);
		double d_i = a_i + b_i - c_i;

		solution C(c_i), D(d_i);
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);
		for (int i = 0; i <= k - 3; i++) {
			if (C.y < D.y) {
				b_i = d_i;
				d_i = c_i;
				D = C;
				c_i = b_i - (double)fib_seq[k - i - 2] / (double)fib_seq[k - i - 1] * (b_i - a_i);
				C.x = c_i;
				C.fit_fun(ff, ud1, ud2);
			}
			else {
				a_i = c_i;
				c_i = d_i;
				C = D;
				d_i = a_i + b_i - c_i;
				D.x = d_i;
				D.fit_fun(ff, ud1, ud2);
			}
		}
		Xopt = C;
		Xopt.flag = 1;

		delete[] fib_seq;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}
}
solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double a_i = a;
		double b_i = b;
		double c_i = (a + b) / 2.0; // punkt środkowy
		solution A(a_i), B(b_i), C(c_i);
		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);
		double d_i = c_i;
		double d_prev = 0.0;
		while (true) {
			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				break;
			}

			double l = m2d(A.y) * (b_i * b_i - c_i * c_i) +
				m2d(B.y) * (c_i * c_i - a_i * a_i) +
				m2d(C.y) * (a_i * a_i - b_i * b_i);

			double m = m2d(A.y) * (b_i - c_i) +
				m2d(B.y) * (c_i - a_i) +
				m2d(C.y) * (a_i - b_i);
			if (m <= 0) {
				Xopt.flag = 0;
				break;
			}
			d_prev = d_i;
			d_i = 0.5 * l / m;
			if (a_i < d_i && d_i < c_i) {
				solution D(d_i);
				D.fit_fun(ff, ud1, ud2);

				if (D.y < C.y) {
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
				D.fit_fun(ff, ud1, ud2);

				if (D.y < C.y) {
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
				Xopt.flag = 0;
				break;
			}
			if (b_i - a_i < epsilon || abs(d_i - d_prev) < gamma) {
				Xopt.flag = 1;
				break;
			}
		}
		Xopt.x = d_i;
		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(XB.x);  // liczba zmiennych
		solution X = XB;        // kopia punktu bazowego

		// Dla każdej zmiennej
		for (int j = 0; j < n; j++) {
			// Próba: x + s*ej
			X.x(j) = XB.x(j) + s;
			X.fit_fun(ff, ud1, ud2);

			if (X.y < XB.y) {
				XB = X;  // Akceptuj ruch
			}
			else {
				// Próba: x - s*ej
				X.x(j) = XB.x(j) - s;
				X.fit_fun(ff, ud1, ud2);

				if (X.y < XB.y) {
					XB = X;  // Akceptuj ruch
				}
				else {
					X.x(j) = XB.x(j);  // Przywróć oryginalną wartość
				}
			}
		}

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution XB(x0);  // Punkt bazowy
		XB.fit_fun(ff, ud1, ud2);

		while (true) {
			// Faza eksploracji
			solution X = HJ_trial(ff, XB, s, ud1, ud2);

			if (solution::f_calls > Nmax) {
				Xopt = XB;
				Xopt.flag = 0;
				return Xopt;
			}

			// Jeśli znaleziono lepszy punkt
			if (X.y < XB.y) {
				// Faza analizy wzorca (pattern move)
				while (true) {
					solution XB_prev = XB;
					XB = X;

					// Ruch wzorcowy: 2*XB - XB_prev
					solution X_pattern(2.0 * XB.x - XB_prev.x);
					X_pattern.fit_fun(ff, ud1, ud2);  // ← DODAJ TĘ LINIĘ!

					// Próba eksploracji z nowego punktu
					X = HJ_trial(ff, X_pattern, s, ud1, ud2);

					if (solution::f_calls > Nmax) {
						Xopt = XB;
						Xopt.flag = 0;
						return Xopt;
					}

					// Jeśli nie ma poprawy, zakończ fazę wzorcową
					if (X.y >= XB.y) {
						X = XB;
						break;
					}
				}
			}
			else {
				// Nie znaleziono poprawy - zmniejsz krok
				s = alpha * s;
			}

			// Warunek zakończenia
			if (s < epsilon) {
				Xopt = XB;
				Xopt.flag = 1;
				break;
			}
		}

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int n = get_len(x0);

		matrix s = s0;
		matrix D = ident_mat(n);  // Stała baza (osie współrzędnych)

		solution XB(x0);
		XB.fit_fun(ff, ud1, ud2);

		int max_iterations = 1000;
		int iter_count = 0;

		while (iter_count < max_iterations) {
			iter_count++;

			bool any_improvement = false;

			for (int j = 0; j < n; j++) {
				matrix d_j = D[j];
				solution X_trial(XB.x + s(j) * d_j);
				X_trial.fit_fun(ff, ud1, ud2);

				if (solution::f_calls > Nmax) {
					Xopt = XB;
					Xopt.flag = 0;
					return Xopt;
				}

				if (X_trial.y < XB.y) {
					XB = X_trial;
					s(j) = alpha * s(j);
					any_improvement = true;
				}
				else {
					s(j) = -beta * s(j);
				}
			}

			// Warunek stopu
			double max_s = abs(s(0));
			for (int j = 1; j < n; j++) {
				if (abs(s(j)) > max_s) max_s = abs(s(j));
			}

			if (max_s < epsilon) {
				Xopt = XB;
				Xopt.flag = 1;
				return Xopt;
			}
		}

		Xopt = XB;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
