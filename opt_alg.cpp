#include"opt_alg.h"
#include <cmath>

// ====== GLOBALNE DO FUNKCJI KARY ======
static matrix(*__ff_base)(matrix, matrix, matrix);
static double __c_pen;
static matrix __ud1_pen;
static matrix __ud2_pen;

// ====== CACHE FF(x) oraz GG(x) ======
#include <unordered_map>
#include <sstream>
#include <iomanip>

static std::unordered_map<std::string, matrix> __cache_ff;
static std::unordered_map<std::string, matrix> __cache_gg;

// HASH punktu
static std::string __hash_x(const matrix& x)
{
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(8);
	for (int i = 0; i < get_len(x); i++)
		oss << x(i) << ",";
	return oss.str();
}

// pobranie ff(x) z cache
static matrix __cached_ff(matrix(*ff)(matrix, matrix, matrix), const matrix& x)
{
	std::string key = __hash_x(x);
	auto it = __cache_ff.find(key);
	if (it != __cache_ff.end()) return it->second;

	matrix xx = x;
	matrix val = ff(xx, __ud1_pen, __ud2_pen);
	__cache_ff[key] = val;
	return val;
}

// pobranie gg(x) z cache
static matrix __cached_gg(const matrix& x)
{
	std::string key = __hash_x(x);
	auto it = __cache_gg.find(key);
	if (it != __cache_gg.end()) return it->second;

	matrix xx = x;
	matrix g = gg3T(xx, __ud1_pen, __ud2_pen);
	__cache_gg[key] = g;
	return g;
}


solution MC(matrix(*ff)(matrix, matrix, matrix), int Ndim, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	solution best;
	best.x = lb;
	best.fit_fun(ff, ud1, ud2);
	best.flag = 0;

	std::random_device rd;
	std::mt19937 gen(rd());
	vector<std::uniform_real_distribution<double>> dist;
	for (int i = 0; i < Ndim; ++i) {
		double a = lb(i);
		double b = ub(i);
		dist.emplace_back(a, b);
	}

	int iter = 0;
	while (solution::f_calls < Nmax) {
		matrix x(Ndim, 1);
		for (int i = 0; i < Ndim; i++) x(i) = dist[i](gen);
		solution s(x);
		s.fit_fun(ff, ud1, ud2);
		if (m2d(s.y) < m2d(best.y)) { best = s; best.flag = 1; }
		iter++;
		if (m2d(best.y) < epsilon) break;
	}
	return best;
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

matrix ff_penalty_wrapper(matrix x, matrix, matrix)
{
	if (!__ff_base) throw string("ff_penalty_wrapper: base function null");

	// użyj cache zamiast surowego wywołania
	matrix f = __cached_ff(__ff_base, x);
	matrix g = __cached_gg(x);

	double S = 0.0;
	int ng = get_len(g);
	for (int i = 0; i < ng; i++) {
		double gi = g(i);
		if (gi > 0.0) S += gi * gi;   // L2 (bezpieczna wersja)
		// jeśli chcesz karę L1 zamiast L2 --> zamień na: S += gi;
	}

	matrix y(1, 1);
	y(0) = m2d(f) + __c_pen * S;
	return y;
}



solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution X = x0;
		__ff_base = ff;
		__ud1_pen = ud1;
		__ud2_pen = ud2;
		// wyczyszczenie cache dla nowej optymalizacji
		__cache_ff.clear();
		__cache_gg.clear();
		int iter = 0;
		while (true) {
			__c_pen = c;
			// run NM on penalized function
			solution Xm = sym_NM(ff_penalty_wrapper, X.x, 0.2, 1.0, 0.5, 2.0, 0.5, epsilon, Nmax / 10, ud1, ud2);
			X = Xm;
			if (m2d(X.y) < epsilon) { X.flag = 1; return X; }
			c *= dc;
			iter++;
			if (solution::f_calls > Nmax || iter > 1) { X.flag = 0; return X; }
		}
	}
	catch (string ex) {
		throw string("solution pen(...):\n") + ex;
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix),
	matrix x0,
	double s,
	double alpha,
	double beta,
	double gamma,
	double delta,
	double epsilon,
	int Nmax,
	matrix ud1,
	matrix ud2)
{
	try {
		int n = get_len(x0);
		// simplex P: n+1 points
		vector<solution> P;
		P.reserve(n + 1);
		// p0 = x0
		P.emplace_back(x0);
		for (int i = 0; i < n; i++) {
			matrix xi = x0;
			xi(i) = xi(i) + s;
			P.emplace_back(xi);
		}
		// evaluate
		for (int i = 0; i <= n; i++) P[i].fit_fun(ff, ud1, ud2);

		while (true) {
			if (solution::f_calls > Nmax) { solution X = P[0]; X.flag = 0; return X; }
			// find best and worst
			int idx_min = 0, idx_max = 0;
			double fmin = m2d(P[0].y);
			double fmax = fmin;
			for (int i = 1; i <= n; i++) {
				double fi = m2d(P[i].y);
				if (fi < fmin) { fmin = fi; idx_min = i; }
				if (fi > fmax) { fmax = fi; idx_max = i; }
			}
			// compute centroid p (n x 1) of all points except idx_max
			matrix p(n, 1);
			for (int i = 0; i < n; i++) p(i) = 0.0;
			for (int i = 0; i <= n; i++) {
				if (i == idx_max) continue;
				for (int j = 0; j < n; j++) {
					p(j) = p(j) + P[i].x(j);
				}
			}
			for (int j = 0; j < n; j++) p(j) = p(j) / static_cast<double>(n);

			// reflection: x_ref = p + alpha*(p - x_max)
			matrix x_ref(n, 1);
			for (int j = 0; j < n; j++) {
				double val = p(j) + alpha * (p(j) - P[idx_max].x(j));
				x_ref(j) = val;
			}
			solution Sref(x_ref);
			Sref.fit_fun(ff, ud1, ud2);
			double f_ref = m2d(Sref.y);

			if (f_ref < m2d(P[idx_min].y)) {
				// expansion
				matrix x_exp(n, 1);
				for (int j = 0; j < n; j++) x_exp(j) = p(j) + gamma * (x_ref(j) - p(j));
				solution Sexp(x_exp); Sexp.fit_fun(ff, ud1, ud2);
				if (m2d(Sexp.y) < f_ref) P[idx_max] = Sexp; else P[idx_max] = Sref;
			}
			else {
				if (f_ref < m2d(P[idx_max].y) && f_ref >= m2d(P[idx_min].y)) {
					P[idx_max] = Sref;
				}
				else {
					// contraction
					matrix x_con(n, 1);
					for (int j = 0; j < n; j++) x_con(j) = p(j) + beta * (P[idx_max].x(j) - p(j));
					solution Scon(x_con); Scon.fit_fun(ff, ud1, ud2);
					if (m2d(Scon.y) < m2d(P[idx_max].y)) {
						P[idx_max] = Scon;
					}
					else {
						// shrink towards best
						matrix x_min = P[idx_min].x;
						for (int i = 0; i <= n; i++) {
							if (i == idx_min) continue;
							matrix xi(n, 1);
							for (int j = 0; j < n; j++) xi(j) = x_min(j) + delta * (P[i].x(j) - x_min(j));
							P[i].x = xi;
							P[i].fit_fun(ff, ud1, ud2);
						}
					}
				}
			}

			// convergence test: max distance to best
			double maxd = 0.0;
			for (int i = 0; i <= n; i++) {
				if (i == idx_min) continue;
				double s2 = 0.0;
				for (int j = 0; j < n; j++) {
					double d = P[idx_min].x(j) - P[i].x(j);
					s2 += d * d;
				}
				double dn = sqrt(s2);
				if (dn > maxd) maxd = dn;
			}
			if (maxd < epsilon) { solution X = P[idx_min]; X.flag = 1; return X; }
		}
	}
	catch (string ex) {
		throw string("solution sym_NM(...):\n") + ex;
	}
}

// =============================================================
// POMOCNICZA FUNKCJA DO SZUKANIA KROKU (Line Search)
// =============================================================
static double h_golden(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix d, int Nmax, matrix ud1, matrix ud2)
{
	double a = 0.0;
	double b = 1.0;
	double epsilon = 1e-7;
	double alpha = 0.6180339887; // (sqrt(5)-1)/2

	double c = b - alpha * (b - a);
	double d_pt = a + alpha * (b - a);

	matrix xc = x0 + d * c;
	matrix xd = x0 + d * d_pt;
	double fc = m2d(ff(xc, ud1, ud2));
	double fd = m2d(ff(xd, ud1, ud2));

	for (int i = 0; i < Nmax; ++i)
	{
		if (fc < fd)
		{
			b = d_pt;
			d_pt = c;
			fd = fc;
			c = b - alpha * (b - a);
			xc = x0 + d * c;
			fc = m2d(ff(xc, ud1, ud2));
		}
		else
		{
			a = c;
			c = d_pt;
			fc = fd;
			d_pt = a + alpha * (b - a);
			xd = x0 + d * d_pt;
			fd = m2d(ff(xd, ud1, ud2));
		}

		if ((b - a) < epsilon) break;
	}

	return (a + b) / 2.0;
}

// =============================================================
// METODA NAJSZYBSZEGO SPADKU (Steepest Descent)
// =============================================================
// =============================================================
// METODA NAJSZYBSZEGO SPADKU (Steepest Descent)
// =============================================================
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution x(x0);
		matrix d;
		double h;

		x.fit_fun(ff, ud1, ud2);

		while (solution::f_calls < Nmax)
		{
			// ============================================================
			// WYPISYWANIE DANYCH DO EXCELA (Wykres ścieżki optymalizacji)
			// Wypisuje: x1, x2
			// ============================================================
			std::cout << x.x(0) << "," << x.x(1) << std::endl;
			// ============================================================

			matrix grad = gf(x.x, ud1, ud2);
			d = -grad;

			if (h0 < 0) h = h_golden(ff, x.x, d, 100, ud1, ud2);
			else h = h0;

			matrix x_new = x.x + d * h;
			solution next_x(x_new);
			next_x.fit_fun(ff, ud1, ud2);

			matrix diff = next_x.x - x.x;
			double dist = sqrt(m2d(trans(diff) * diff));

			if (dist < epsilon) {
				Xopt = next_x;
				Xopt.flag = 1;
				return Xopt;
			}
			x = next_x;
		}

		Xopt = x;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

// =============================================================
// METODA GRADIENTÓW SPRZĘŻONYCH (Conjugate Gradient)
// =============================================================
solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution x(x0);
		x.fit_fun(ff, ud1, ud2);

		matrix grad = gf(x.x, ud1, ud2);
		matrix d = -grad;
		matrix grad_old = grad;

		while (solution::f_calls < Nmax)
		{
			double h;
			if (h0 < 0) h = h_golden(ff, x.x, d, 100, ud1, ud2);
			else h = h0;

			matrix x_new = x.x + d * h;
			solution next_x(x_new);
			next_x.fit_fun(ff, ud1, ud2);

			matrix diff = next_x.x - x.x;
			double dist = sqrt(m2d(trans(diff) * diff));
			if (dist < epsilon) {
				Xopt = next_x;
				Xopt.flag = 1;
				return Xopt;
			}

			x = next_x;

			matrix grad_new = gf(x.x, ud1, ud2);
			double num = m2d(trans(grad_new) * grad_new);
			double den = m2d(trans(grad_old) * grad_old);
			double beta = num / den;

			d = -grad_new + d * beta;
			grad_old = grad_new;
		}

		Xopt = x;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

// =============================================================
// METODA NEWTONA
// =============================================================
solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution x(x0);
		x.fit_fun(ff, ud1, ud2);
		matrix d;

		while (solution::f_calls < Nmax)
		{
			matrix grad = gf(x.x, ud1, ud2);
			matrix H = Hf(x.x, ud1, ud2);
			d = -inv(H) * grad;

			double h;
			if (h0 < 0) h = h_golden(ff, x.x, d, 100, ud1, ud2);
			else h = h0;

			matrix x_new = x.x + d * h;
			solution next_x(x_new);
			next_x.fit_fun(ff, ud1, ud2);

			matrix diff = next_x.x - x.x;
			double dist = sqrt(m2d(trans(diff) * diff));
			if (dist < epsilon) {
				Xopt = next_x;
				Xopt.flag = 1;
				return Xopt;
			}
			x = next_x;
		}

		Xopt = x;
		Xopt.flag = 0;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

// =============================================================
// METODA ZŁOTEGO PODZIAŁU 1D
// =============================================================
solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double alpha = 0.6180339887;
		double a_i = a;
		double b_i = b;
		double c_i = b_i - alpha * (b_i - a_i);
		double d_i = a_i + alpha * (b_i - a_i);

		solution C(c_i); C.fit_fun(ff, ud1, ud2);
		solution D(d_i); D.fit_fun(ff, ud1, ud2);

		while (solution::f_calls < Nmax)
		{
			if (C.y < D.y) {
				b_i = d_i;
				d_i = c_i;
				D = C;
				c_i = b_i - alpha * (b_i - a_i);
				C.x = c_i;
				C.fit_fun(ff, ud1, ud2);
			}
			else {
				a_i = c_i;
				c_i = d_i;
				C = D;
				d_i = a_i + alpha * (b_i - a_i);
				D.x = d_i;
				D.fit_fun(ff, ud1, ud2);
			}

			if ((b_i - a_i) < epsilon) {
				Xopt.x = (a_i + b_i) / 2.0;
				Xopt.fit_fun(ff, ud1, ud2);
				Xopt.flag = 1;
				return Xopt;
			}
		}

		Xopt.x = (a_i + b_i) / 2.0;
		Xopt.fit_fun(ff, ud1, ud2);
		Xopt.flag = 0;
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
		// Tu wpisz kod funkcji (na razie puste)
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
		// Tu wpisz kod funkcji (na razie puste)
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}


