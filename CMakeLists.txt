cmake_minimum_required(VERSION 3.20)
project(Optymalizacja)

set(CMAKE_CXX_STANDARD 17)

# Tu wymieniamy wszystkie pliki .cpp, które widzę w projekcie
add_executable(Optymalizacja
        main.cpp
        matrix.cpp
        ode_solver.cpp
        opt_alg.cpp
        solution.cpp
        user_funs.cpp
)