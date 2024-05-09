#ifndef DATAFILE_H
#define DATAFILE_H

#include <iostream>
#include <fstream> //library for file handling
#include <string>
#include <chrono>

void write_solverData(std::fstream* file, double alpha, const int xdiv, const int nEqs, int nEqsCond){
    file->open("data.txt", std::ios::app); //opening file in write mode

    if(file->is_open()){
        *file << "xdiv = " << xdiv << std::endl;
        *file << "alpha = " << alpha << std::endl;
        *file << "\n";
        *file << nEqs << std::endl;
        *file << nEqsCond << std::endl;
        file->close();
    }
}

void write_directData(std::fstream* file, const int xdiv, const int nEqs, int nEqsCond){
    file->open("dataDirect.txt", std::ios::app); //opening file in write mode

    if(file->is_open()){
        *file << "xdiv = " << xdiv << std::endl;
        *file << "\n";
        *file << nEqs << std::endl;
        *file << nEqsCond << std::endl;
        file->close();
    }
}

void write_solverResults(std::fstream* file, const int iterations, std::chrono::steady_clock::time_point beginS, 
                std::chrono::steady_clock::time_point endS, std::chrono::steady_clock::time_point beginI, std::chrono::steady_clock::time_point endI){ 
    file->open("data.txt", std::ios::app); //opening file in write mode

    if(file->is_open()){
        if (iterations > 200) {
            *file << "diverged" << std::endl;
        } else {
            *file << iterations << std::endl;
        }
        *file << std::chrono::duration_cast<std::chrono::milliseconds>(endS - beginS).count() << std::endl;
        *file << std::chrono::duration_cast<std::chrono::milliseconds>(endI - beginI).count() << std::endl;
        *file << "\n";
        *file << "\n";
        file->close();
    }
}

void write_directResults(std::fstream* file, std::chrono::steady_clock::time_point beginS, 
                std::chrono::steady_clock::time_point endS){ 
    file->open("dataDirect.txt", std::ios::app); //opening file in write mode

    if(file->is_open()){
        *file << std::chrono::duration_cast<std::chrono::milliseconds>(endS - beginS).count() << std::endl;
        *file << "\n";
        *file << "\n";
        file->close();
    }
}

void write_Norms(std::fstream* file, double norm_dsol, double norm_rhs){ 
    file->open("norms.txt", std::ios::app); //opening file in write mode

    if(file->is_open()){
        *file << norm_dsol << " " << norm_rhs << std::endl;
        file->close();
    }
}


#endif

//write_file(&my_file, alpha, xdiv, nit, nEquationsFull, nEquationsCondensed, begin, end, begin1, end1);