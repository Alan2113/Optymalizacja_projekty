#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);  // Funkcja testowa lab1
matrix ff1R(matrix, matrix = NAN, matrix = NAN);  // Problem rzeczywisty lab1
matrix df1(double, matrix, matrix = NAN, matrix = NAN);  // Równania ró¿niczkowe lab1
matrix ff2T(matrix, matrix = NAN, matrix = NAN);  // Funkcja testowa lab2
matrix ff2R(matrix, matrix = NAN, matrix = NAN);  // Problem rzeczywisty lab2
matrix df2(double, matrix, matrix = NAN, matrix = NAN);  // Równania ró¿niczkowe lab2