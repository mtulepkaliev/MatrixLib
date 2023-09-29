//
//  Helper functions for matrices
//
//  Created by Moses T on 9/9/22.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include "Matrix.hpp"
#include "Vector.hpp"

using std::cout,std::vector,std::ofstream,std::endl,std::string,std::ifstream,std::getline,std::stringstream;
#ifndef PA1_helper_h
#define PA1_helper_h

//function prototypes
vector<double> cramersRuleSolution2x2(vector<vector<double>> matA, vector<double> vecB);
double get2x2Determinant(vector<vector<double>> matrix);
double get2x2ColDeterminant(vector<vector<double>> matrix, vector<double> colVec,int col);
vector<double> getRoots(double a, double b,double c);
vector<double> getColVec(vector<vector<double>> matrix,unsigned int col);
double vecDotProduct(vector<double> vec1,vector<double> vec2);
void printMatrix(vector<vector<double>> matrix,char colDelimiter,char rowDelimiter);
void writeMatrixToFile(vector<vector<double>> matrix,string outFile,bool append,char colDelimiter,char rowDelimiter);
void writeVecToFile(vector<double> vec,string outFile,char rowDelimiter);
vector<vector<double>> readMatrixFromFile(string filename,char colDelimeter,char rowDelimiter);
vector<vector<double>> addMatrices(vector<vector<double>> matA,vector<vector<double>> matB);
vector<vector<double>> multMatrices(vector<vector<double>> matA,vector<vector<double>> matB);



//solves a 2x2 system of equations using cramer's rule
Vec cramersRuleSolution2x2(Matrix matA, Vec vecB)
{
    Vec vecX(2);
    
    //get our determinants
    double matADeter = matA.determinant();
    
    Matrix b0Mat(matA);
    Matrix b1Mat(matA);
    b0Mat.insertColVec(vecB, 0);
    b1Mat.insertColVec(vecB, 1);
    double b0Deter = b0Mat.determinant();
    double b1Deter = b1Mat.determinant();
    
    //if the determinarte is not 0 then find the 2 solutions
    if(matADeter != 0)
    {
        vecX[0] = b0Deter/matADeter;
        vecX[1] = b1Deter/matADeter;
        return vecX;
    }
    //if both the numertor determiantes are not zero then it is inconsistent
    else if (b0Deter != 0 || b1Deter != 0)
    {
        throw std::logic_error("System Inconsistent");
    }
    //otherwise it is unsolveable
    else
    {
        throw std::logic_error("System Underdetermined");
    }
}

//uses quadratic formula to get the roots
vector<double> getRoots(double a, double b,double c)
{
    //vector for the 2 roots
    vector<double> roots(2,0);
    
    //caclualte the discriminant (value inside the square root)
    double discrim = pow(b,2)-4*a*c;
    
    //no real eigenvalkues if the discriminant is less than 0
    if(discrim < 0)
    {
        throw std::logic_error("No real eigenvalues");
    }
    
    //calculate the 2 roots
    roots[0] = (-1 * b + sqrt(discrim)) / (2 * a);
    roots[1] = (-1 * b - sqrt(discrim)) / (2 * a);
    
    //swap if needed so the dominant eigenvalue is the first one
    if(abs(roots[0]) < abs(roots[1]))
    {
        std::swap(roots[0], roots[1]);
    }
    
    return roots;
}

//calculates 2d area of a triangle using area of triangle in determinant form
double areaTriangle2D(Matrix pointMat)
{
    //create a copy of the input matrix
    Matrix areaMat(pointMat);
    
    //insert a row vector of 1's as the first row of the matrix
    areaMat.insertRowVec(Vec(vector<double>({1,1,1})), -1);
    
    //area = 1/2 * determinant
    return .5 * abs(areaMat.determinant());
}

//calculates area of a triangle using half cross product formula
double areaTriangle3D(Matrix pointMat)
{
    //extract points from the matrix
    Vec pointA(pointMat.inRep,0);
    Vec pointB(pointMat.inRep,1);
    Vec pointC(pointMat.inRep,2);
    
    //declare 2 vectors between points b and a and c and a
    Vec pointBA(pointB - pointA);
    Vec pointCA(pointC - pointA);
    
    //calculate the cross product of the vectors
    Vec cross(pointBA ^ pointCA);
    
    //return 1/2 multiplied by the magnitude of the cross product
    return .5 * cross.magnitude();
}

#endif /* PA1_helper_h */
