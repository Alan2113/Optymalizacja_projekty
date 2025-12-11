#include "user_funs.h"
#include "ode_solver.h"
#include <cmath>
#include <limits>

// Stała PI (MSVC nie posiada M_PI)
static const double PI = 3.14159265358979323846;

// Stałe fizyczne
static const double rho_air = 1.2;
static const double m_ball = 0.6;
static const double r_ball = 0.12;
static const double S_ball = PI * r_ball * r_ball;
static const double C_d = 0.47;
static const double g_acc = 9.81;

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera wartość funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera współrzędne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartość funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki początkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment siły działający na wahadło oraz czas działania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwiązujemy równanie różniczkowe
	int n = get_len(Y[0]);									// długość rozwiązania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahadła
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// wartość funkcji celu (ud1 to założone maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pamięci rozwiązanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z położenia to prędkość
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z prędkości to przyspieszenie
	return dY;
}


// ========== FUNKCJA TESTOWA (LAB 1) ==========
// f(x) = -cos(0.1x) * e^(-(0.1x-2π)²) + 0.002 * (0.1x)²
matrix ff1T(matrix x, matrix ud1, matrix ud2)
{

	//double M_PI = 3.141592653589793238;
	matrix y;
	double x_val = m2d(x);
	double x_scaled = 0.1 * x_val;
	double term1 = -cos(x_scaled) * exp(-pow(x_scaled - 2 * M_PI, 2));
	double term2 = 0.002 * pow(x_scaled, 2);
	y = term1 + term2;
	return y;
}

// ========== PROBLEM RZECZYWISTY (LAB 1) - ZBIORNIKI ==========
// Funkcja celu: minimalizuje (T_B_max - 50)^2
// x to pole przekroju DA w cm^2
matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	// x = DA [cm^2]
	double DA = m2d(x);

	// Parametry zbiorników
	const double PA = 2.0;        // [m^2]
	const double PB = 1.0;        // [m^2]
	const double DB = 36.5665;    // [cm^2]
	const double FBin = 0.01;     // [m^3/s] (10 l/s)
	const double TBin = 20.0;     // [°C]

	// Warunki początkowe
	double VA = 5.0;              // [m^3]
	double VB = 1.0;              // [m^3]
	double TA = 95.0;             // [°C]
	double TB = 20.0;             // [°C]

	// Parametry symulacji
	const double t0 = 0.0;
	const double dt = 1.0;
	const double tend = 2000.0;
	const int N = (int)((tend - t0) / dt) + 1;

	// Stałe fizyczne
	const double a_coeff = 0.98;
	const double b_coeff = 0.63;
	const double g = 9.81;        // [m/s^2]

	// Konwersja DA i DB z cm^2 na m^2
	double DA_m2 = DA * 1e-4;
	double DB_m2 = DB * 1e-4;

	double T_B_max = TB;

	// Symulacja Eulera
	for (int i = 0; i < N; i++) {
		// Przepływ z A do B
		double FAout = 0.0;
		if (VA > 0) {
			FAout = a_coeff * b_coeff * DA_m2 * sqrt(2.0 * g * VA / PA);
		}

		// Przepływ z B na zewnątrz
		double FBout = 0.0;
		if (VB > 0) {
			FBout = a_coeff * b_coeff * DB_m2 * sqrt(2.0 * g * VB / PB);
		}

		// Aktualizacja objętości
		double dVA_dt = -FAout;
		double dVB_dt = FAout + FBin - FBout;

		VA = VA + dVA_dt * dt;
		VB = VB + dVB_dt * dt;

		// Zabezpieczenie przed ujemnymi objętościami
		if (VA < 0) VA = 0;
		if (VB < 0) VB = 0;

		// Aktualizacja temperatury zbiornika A
		if (VA > 1e-6) {
			// Brak dopływu do A, tylko wypływ
			// TA pozostaje stałe (lub można modelować chłodzenie)
		}

		// Aktualizacja temperatury zbiornika B
		if (VB > 1e-6) {
			double Vin_total = FAout + FBin;
			double Tin_avg = 0.0;
			if (Vin_total > 1e-9) {
				Tin_avg = (FAout * TA + FBin * TBin) / Vin_total;
			}
			double dTB_dt = (Vin_total / VB) * (Tin_avg - TB);
			TB = TB + dTB_dt * dt;
		}

		// Śledź maksymalną temperaturę w B
		if (TB > T_B_max) {
			T_B_max = TB;
		}
	}

	// Funkcja celu: minimalizuj odchylenie od 50°C
	matrix y;
	y = pow(T_B_max - 50.0, 2);
	return y;
}
// ========== FUNKCJA TESTOWA (LAB 2) ==========
// f(x1, x2) = x1^2 + x2^2 - cos(2.5π*x1) - cos(2.5π*x2) + 2
matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	//const double M_PI = 3.141592653589793238;
	matrix y;
	double x1 = x(0);
	double x2 = x(1);

	y = pow(x1, 2) + pow(x2, 2) - cos(2.5 * M_PI * x1) - cos(2.5 * M_PI * x2) + 2.0;
	return y;
}

// ========== PROBLEM RZECZYWISTY (LAB 2) - ROBOT ==========
// x = [k1, k2]^T - współczynniki wzmocnienia regulatora
matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	double k1 = x(0);
	double k2 = x(1);

	// Parametry ramienia
	const double l = 2.0;           // długość ramienia [m]
	const double mr = 1.0;          // masa ramienia [kg]
	const double mc = 5.0;          // masa ciężarka [kg]
	const double b = 0.25;          // współczynnik tarcia [Nms]

	// Moment bezwładności
	double I = (1.0 / 3.0) * mr * pow(l, 2) + mc * pow(l, 2);

	// Referencje
	//float M_PI = 3.141592653589793238;
	const double alpha_ref = M_PI;  // π rad
	const double omega_ref = 0.0;   // 0 rad/s

	// Parametry symulacji
	const double t0 = 0.0;
	const double dt = 0.1;
	const double tend = 100.0;
	const int N = (int)((tend - t0) / dt);

	// Warunki początkowe: alpha(0) = 0, omega(0) = 0
	double alpha = 0.0;
	double omega = 0.0;

	// Całka dla funkcjonału jakości (metoda prostokątów)
	double Q = 0.0;

	for (int i = 0; i < N; i++) {
		// Moment siły
		double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

		// Składniki funkcjonału
		double term1 = 10.0 * pow(alpha_ref - alpha, 2);
		double term2 = pow(omega_ref - omega, 2);
		double term3 = pow(M, 2);

		Q += (term1 + term2 + term3) * dt;

		// Równanie ruchu: I * d²α/dt² + b * dα/dt = M(t)
		// omega = dα/dt
		// d_omega/dt = (M - b*omega) / I

		double d_omega = (M - b * omega) / I;

		// Aktualizacja metodą Eulera
		omega = omega + d_omega * dt;
		alpha = alpha + omega * dt;
	}

	matrix y;
	y = Q;
	return y;
}

// ========== RÓWNANIA RÓŻNICZKOWE DLA ROBOTA (LAB 2) ==========
// Y = [alpha, omega]^T
// ud1 = [k1, k2]^T
matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	// Parametry
	const double l = 2.0;
	const double mr = 1.0;
	const double mc = 5.0;
	const double b = 0.25;
	double I = (1.0 / 3.0) * mr * pow(l, 2) + mc * pow(l, 2);

	// Referencje
	//const double M_PI = 3.141592653589793238;
	const double alpha_ref = M_PI;
	const double omega_ref = 0.0;

	// Współczynniki regulatora
	double k1 = ud1(0);
	double k2 = ud1(1);

	// Stan
	double alpha = Y(0);
	double omega = Y(1);

	// Moment siły
	double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

	// Pochodne
	matrix dY(2, 1);
	dY(0) = omega;                    // dα/dt = ω
	dY(1) = (M - b * omega) / I;      // dω/dt = (M - b*ω) / I

	return dY;
}






// ========== RÓWNANIA RÓŻNICZKOWE DLA ZBIORNIKÓW (LAB 1) ==========
// Używane do generowania wykresów po znalezieniu optymalnego DA
// Y = [VA, VB, TA, TB]^T
// ud1 = DA [cm^2]
matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	// Parametry
	const double PA = 2.0;        // [m^2]
	const double PB = 1.0;        // [m^2]
	const double DB = 36.5665;    // [cm^2]
	const double FBin = 0.01;     // [m^3/s]
	const double TBin = 20.0;     // [°C]
	const double a_coeff = 0.98;
	const double b_coeff = 0.63;
	const double g = 9.81;

	double DA = m2d(ud1);         // [cm^2]
	double DA_m2 = DA * 1e-4;     // [m^2]
	double DB_m2 = DB * 1e-4;     // [m^2]

	// Stan aktualny
	double VA = Y(0);
	double VB = Y(1);
	double TA = Y(2);
	double TB = Y(3);

	// Przepływy
	double FAout = 0.0;
	if (VA > 0) {
		FAout = a_coeff * b_coeff * DA_m2 * sqrt(2.0 * g * VA / PA);
	}

	double FBout = 0.0;
	if (VB > 0) {
		FBout = a_coeff * b_coeff * DB_m2 * sqrt(2.0 * g * VB / PB);
	}

	// Pochodne
	matrix dY(4, 1);

	// dVA/dt
	dY(0) = -FAout;

	// dVB/dt
	dY(1) = FAout + FBin - FBout;

	// dTA/dt (temperatura w A nie zmienia się, bo nie ma dopływu)
	dY(2) = 0.0;

	// dTB/dt
	if (VB > 1e-6) {
		double Vin_total = FAout + FBin;
		double Tin_avg = 0.0;
		if (Vin_total > 1e-9) {
			Tin_avg = (FAout * TA + FBin * TBin) / Vin_total;
		}
		dY(3) = (Vin_total / VB) * (Tin_avg - TB);
	}
	else {
		dY(3) = 0.0;
	}

	return dY;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
	matrix y(1, 1);
	double x1 = x(0);
	double x2 = x(1);

	// domyślne a
	double a = 4.0;
	// bezpieczne pobranie ud1
	try { if (get_len(ud1) == 1) a = ud1(0); }
	catch (...) {}

	y(0) = (x1 - a) * (x1 - a) + (x2 - 1.0) * (x2 - 1.0);
	return y;
}

matrix gg3T(matrix x, matrix ud1, matrix ud2)
{
	matrix g(3, 1);
	double x1 = x(0);
	double x2 = x(1);
	g(0) = x1 + x2 - 2.0;
	g(1) = x1 * x1 + x2 - 2.0;
	g(2) = -x1;
	return g;
}

// ================= df3 =================
// Y = [vx; vy; x; y]  (kolumnowy wektor stanu)
// ud1(0) = omega (rad/s)
// ================= df3 =================
// Y = [vx; vy; x; y]  (kolumnowy wektor stanu)
// ud1(0) = omega (rad/s)
matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(4, 1);
	// stan
	double vx = Y(0);
	double vy = Y(1);
	// parametry z instrukcji
	const double rho = 1.2;
	const double m = 0.6;
	const double r = 0.12;
	const double S = 3.14159265358979323846 * r * r;
	const double C = 0.47;
	const double g = 9.81;

	double omega = 0.0;
	if (get_len(ud1) == 1) {
		try { omega = ud1(0); }
		catch (...) { omega = 0.0; }
	}

	// --- Siły wg instrukcji ---
	double Dx = 0.5 * C * rho * S * vx * fabs(vx);   // Dx = 1/2 C rho S vx |vx|
	double Dy = 0.5 * C * rho * S * vy * fabs(vy);   // Dy = 1/2 C rho S vy |vy|

	double FMx = rho * vy * omega * 3.14159265358979323846 * r * r * r; // FMx = rho * vy * omega * pi * r^3
	double FMy = rho * vx * omega * 3.14159265358979323846 * r * r * r; // FMy = rho * vx * omega * pi * r^3

	// równania ruchu zgodne z zadaniem:
	// m d2x/dt2 + Dx + FMx = 0  -> ax = -(Dx + FMx)/m
	// m d2y/dt2 + Dy + FMy = - m g  -> ay = -g - (Dy + FMy)/m
	double ax = -(Dx + FMx) / m;
	double ay = -g - (Dy + FMy) / m;

	dY(0) = ax;
	dY(1) = ay;
	dY(2) = vx;
	dY(3) = vy;
	return dY;
}


// ================= ff3R =================
// x = [v0x; omega]
// Cel: maksymalizacja x_end (pierwsze przecięcie y<=0) z ograniczeniem x przy y=50 w [3,7].
// Zwracamy wartość do minimalizacji: yout = -x_end + PEN_COEFF*penalty
matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	// decyzje
	double v0x = x(0);
	double omega = x(1);

	// jeżeli wartości wyraźnie poza dopuszczalnym przedziałem, dodaj karę
	const double LIM = 10.0;
	bool out_of_bounds = (v0x < -LIM || v0x > LIM || omega < -LIM || omega > LIM);

	// warunki poczatkowe
	matrix Y0(4, 1);
	Y0(0) = v0x;    // vx(0)
	Y0(1) = 0.0;    // vy(0) = 0
	Y0(2) = 0.0;    // x(0)
	Y0(3) = 100.0;  // y(0)

	// ud dla ODE (omega)
	matrix udw(1, 1); udw(0) = omega;

	// integracja
	double t0 = 0.0, dt = 0.01, tend = 7.0;
	matrix* S = nullptr;
	try {
		S = solve_ode(df3, t0, dt, tend, Y0, udw, matrix());
	}
	catch (...) {
		// solver nie wykonał się -> ogromny koszt
		matrix yout(1, 1); yout(0) = 1e9; return yout;
	}
	// rozmiary: S[1] jest macierzą (N x 4)
	int* sz = get_size(S[1]);
	int N = sz[0];
	delete[] sz;
	if (N < 2) { delete[] S; matrix yout(1, 1); yout(0) = 1e9; return yout; }

	// --- znajdz pierwsze przecięcie z ziemią (y_prev > 0, y_curr <= 0) i interpoluj x_end
	double x_end = NAN;
	bool hit_ground = false;
	for (int i = 1; i < N; ++i)
	{
		double y_prev = S[1](i - 1, 3);
		double y_curr = S[1](i, 3);
		if (y_prev > 0.0 && y_curr <= 0.0) {
			double x_prev = S[1](i - 1, 2);
			double x_curr = S[1](i, 2);
			if (y_curr == y_prev) x_end = x_curr;
			else {
				double tau = (0.0 - y_prev) / (y_curr - y_prev);
				x_end = x_prev + tau * (x_curr - x_prev);
			}
			hit_ground = true;
			break;
		}
	}
	// jeśli nie uderzyła w ziemię w czasie tend, weź pozycję końcową (karane)
	if (!hit_ground) x_end = S[1](N - 1, 2);

	// --- znajdz x przy y=50 (pierwsze przecięcie przy spadku: y_prev >=50 && y_curr <=50)
	double x_at_y50 = NAN;
	bool reached50 = false;
	for (int i = 1; i < N; ++i)
	{
		double y_prev = S[1](i - 1, 3);
		double y_curr = S[1](i, 3);
		if (y_prev >= 50.0 && y_curr <= 50.0) {
			double x_prev = S[1](i - 1, 2);
			double x_curr = S[1](i, 2);
			if (y_curr == y_prev) x_at_y50 = x_curr;
			else {
				double tau = (50.0 - y_prev) / (y_curr - y_prev);
				x_at_y50 = x_prev + tau * (x_curr - x_prev);
			}
			reached50 = true;
			break;
		}
	}

	// free memory from solve_ode
	delete[] S;

	// --- oblicz cel i kary
	double obj = -x_end; // minimalizator: chcemy maksymalizowac x_end

	const double PEN_COEFF = 1e5;
	double penalty = 0.0;

	if (out_of_bounds) {
		penalty += 1.0; // naruszenie prostych ograniczen
	}

	// jeśli nie osiągnął 50m podczas spadku -> kara (może oznaczać zbyt wolny spadek)
	if (!reached50) {
		penalty += 1.0;
	}
	else {
		// sprawdz warunek x_at_y50 in [3,7]
		if (!(x_at_y50 >= 3.0 && x_at_y50 <= 7.0)) {
			double dist = 0.0;
			if (x_at_y50 < 3.0) dist = 3.0 - x_at_y50;
			else if (x_at_y50 > 7.0) dist = x_at_y50 - 7.0;
			penalty += dist * dist;
		}
	}

	matrix yout(1, 1);
	yout(0) = obj + PEN_COEFF * penalty;
	// zabezpieczenie numeryczne: nie zwracamy NaN/inf
	if (!(yout(0) == yout(0)) || isinf(yout(0))) { yout(0) = 1e9; }
	return yout;
}

// ==========================================================
//              LAB 4 - METODY GRADIENTOWE
// ==========================================================

// --- Funkcja testowa (Test Function) ---
// f(x) = 1/6*x1^6 - 1.05*x1^4 + 2*x1^2 + x2^2 + x1*x2
matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    double x1 = x(0);
    double x2 = x(1);

    y = (1.0 / 6.0) * pow(x1, 6) - 1.05 * pow(x1, 4) + 2.0 * pow(x1, 2) + pow(x2, 2) + x1 * x2;
    return y;
}

// Gradient funkcji testowej
matrix gf4T(matrix x, matrix ud1, matrix ud2)
{
    matrix g(2, 1);
    double x1 = x(0);
    double x2 = x(1);

    g(0) = pow(x1, 5) - 4.2 * pow(x1, 3) + 4.0 * x1 + x2; // pochodna po x1
    g(1) = 2.0 * x2 + x1;                                 // pochodna po x2
    return g;
}

// Hesjan funkcji testowej (do metody Newtona)
matrix Hf4T(matrix x, matrix ud1, matrix ud2)
{
    matrix H(2, 2);
    double x1 = x(0);
    // double x2 = x(1); // niewykorzystane w drugich pochodnych

    H(0, 0) = 5.0 * pow(x1, 4) - 12.6 * pow(x1, 2) + 4.0; // d2f / dx1^2
    H(0, 1) = 1.0;                                        // d2f / dx1 dx2
    H(1, 0) = 1.0;                                        // d2f / dx2 dx1
    H(1, 1) = 2.0;                                        // d2f / dx2^2
    return H;
}


// --- Problem Rzeczywisty (Real Problem) - Regresja Logistyczna ---
// ud1 = X (macierz danych wejściowych, 3 x m)
// ud2 = Y (wektor etykiet, 1 x m)
// theta = x (szukane parametry, 3 x 1)

// Wklej to do user_funs.cpp zamiast starego ff4R i gf4R

// --- Problem Rzeczywisty (Real Problem) - Regresja Logistyczna ---
matrix ff4R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix X = ud1;
	matrix Y = ud2;

	// --- POPRAWKA TUTAJ ---
	// Pobieramy rozmiar macierzy Y (1x100), żeby znać liczbę m
	int* s = get_size(Y);
	int m = s[1]; // liczba kolumn to liczba próbek
	delete[] s;   // trzeba zwolnić pamięć po get_size
	// ----------------------

	double J = 0.0;

	for (int i = 0; i < m; ++i)
	{
		matrix xi = X[i];
		double z = (trans(x) * xi)(0,0);
		double h = 1.0 / (1.0 + exp(-z));

		// Zabezpieczenie logarytmu
		if (h < 1e-15) h = 1e-15;
		if (h > 1.0 - 1e-15) h = 1.0 - 1e-15;

		double yi = Y(0, i);

		J += yi * log(h) + (1.0 - yi) * log(1.0 - h);
	}

	y = -1.0 / m * J;
	return y;
}

matrix gf4R(matrix x, matrix ud1, matrix ud2)
{
	matrix X = ud1;
	matrix Y = ud2;

	// --- POPRAWKA TUTAJ ---
	int* s = get_size(Y);
	int m = s[1];
	delete[] s;
	// ----------------------

	matrix g(3, 1);
	g(0) = 0; g(1) = 0; g(2) = 0;

	for (int i = 0; i < m; ++i)
	{
		matrix xi = X[i];
		double z = (trans(x) * xi)(0,0);
		double h = 1.0 / (1.0 + exp(-z));
		double yi = Y(0, i);

		double error = h - yi;

		g = g + xi * error;
	}

	g = g * (1.0 / m);
	return g;
}

