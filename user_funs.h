#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix ff1T(matrix, matrix = NAN, matrix = NAN);  // Funkcja testowa lab1
matrix ff1R(matrix, matrix = NAN, matrix = NAN);  // Problem rzeczywisty lab1
matrix df1(double, matrix, matrix = NAN, matrix = NAN);  // R�wnania r�niczkowe lab1
matrix ff2T(matrix, matrix = NAN, matrix = NAN);  // Funkcja testowa lab2
matrix ff2R(matrix, matrix = NAN, matrix = NAN);  // Problem rzeczywisty lab2
matrix df2(double, matrix, matrix = NAN, matrix = NAN);  // R�wnania r�niczkowe lab2
matrix ff3T(matrix x, matrix ud1, matrix ud2);
matrix gg3T(matrix x, matrix ud1, matrix ud2);
matrix ff3R(matrix x, matrix ud1, matrix ud2);
matrix df3(double t, matrix Y, matrix ud1, matrix ud2);
matrix ff4T(matrix x, matrix ud1, matrix ud2);
matrix gf4T(matrix x, matrix ud1, matrix ud2);
matrix Hf4T(matrix x, matrix ud1, matrix ud2);

matrix ff4R(matrix x, matrix ud1, matrix ud2);
matrix gf4R(matrix x, matrix ud1, matrix ud2);
