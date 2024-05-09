#ifndef CREATEFILE_H
#define CREATEFILE_H

#include <iostream>
#include <fstream> //library for file handling
#include <string>

void create_file(std::fstream* file){
    file->open("data.txt", std::ios::out); //opening file in write mode

    if(file->is_open()){
        *file << "xdiv = " << std::endl;
        *file << "alpha = " << std::endl;
        *file << "\n";
        *file << "Number of equations = " <<std::endl;
        *file << "Number of equations condensed = " <<std::endl;
        *file << "Number of iterations = " <<std::endl;
        *file << "Time for Initial Solver = " <<std::endl;
        *file << "Time for Iterative Process = " <<std::endl;
        *file << "\n";
        *file << "\n";
        file->close();
    }
    else{
        std::cout << "File not created!";
    }
}

void create_file2(std::fstream* file){
    file->open("dataDirect.txt", std::ios::out); //opening file in write mode

    if(file->is_open()){
        *file << "xdiv = " << std::endl;
        *file << "\n";
        *file << "Number of equations = " <<std::endl;
        *file << "Number of equations condensed = " <<std::endl;
        *file << "Time for Solver = " <<std::endl;
        *file << "\n";
        *file << "\n";
        file->close();
    }
    else{
        std::cout << "File not created!";
    }
}

void creat_file3(std::fstream* file, double xdiv, double alpha){
    file->open("norms.txt", std::ios::app); //opening file in write mode

    if(file->is_open()){
        *file << "\n" << std::endl;
        *file << xdiv << " " << alpha << std::endl;
        file->close();
    }
    else{
        std::cout << "File not created!";
    }
}


//alpha REAL
//xdiv const int
//iterations int
//nEqs const int
//nEqsCond const int
//timeSolver

#endif 