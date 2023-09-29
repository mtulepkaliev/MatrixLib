//
//  Matrix.hpp
//  Program3
//
//  Created by Moses T on 10/28/22.
//
#include <vector>
#include "Vector.hpp"
#include <math.h>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

using std::cout,std::vector,std::ofstream,std::endl,std::string,std::ifstream,std::getline,std::stringstream;
#ifndef Matrix_h
#define Matrix_h

class Matrix;
class Vec;
Matrix operator*(double scalar, Matrix MatA);
Matrix operator+(const Matrix& A,const Matrix& B);
vector<vector<double>> readMatrixFromFile(string filename,char colDelimeter,char rowDelimiter);
bool doubleEquals(double d1, double d2);


class Matrix
{
public:
    double inf_norm();
    int numRows;
    int numCols;
    vector<vector<double>> inRep;
    
    void swapRows(unsigned int row1, unsigned int row2)
    {
        std::swap(inRep[row1],inRep[row2]);
    }
    //gets the eigenvector of an input matrix given the eigenvalue
    Vec getEigenvector(double lambda)
    {
        Vec eigenvector(2);
        //calculates A - lambda I
        
        //gets the matrix of -lambda I
        Matrix negLamdaI((-1 * lambda * identity(2)));
        
        //add A to it
        Matrix aMinusLamdaI((*this) + negLamdaI);
        
        //create a shear matrix by declaring and dientity matrix and set (1,0) to -v2/v1
        Matrix gaussMatrix(identity(2));
        gaussMatrix[1][0] = -1 * (inRep[1][0]/inRep[0][0]);
        
        //multiply aMinusLambdaI by the shear matrix to perform gauss elimination
        gaussMatrix = gaussMatrix * aMinusLamdaI;
        
        //let the system of 2 equatons be defined by
        // ar1 + br2 = 0 and
        // cr1 + dr2 = 0
        // since we performed gauss elimination c and d are
        // now equal to 0 so we are left with
        // ar1 + br2 = 0
        //if we let r1 = 1 then we get
        // a + br2 = 0
        //therefor r2 = -a/b
        
        //ser r1 to 1
        eigenvector[0] = 1;
        
        //define a and b for the previous equation
        double a = gaussMatrix[0][0];
        double b = gaussMatrix[0][1];
        
        //caculate r2 as a-b
        eigenvector[1] = -1 * a/b;
        
        //normalize the eigenvector
        eigenvector.normalize();
        //eigenvector = normalizeVector(eigenvector);
        return eigenvector;
    }
    
    //determinant of a matrix
    double determinant()
    {
        double deter = 0;
        //check the size to see if we are calculating the determinant of a 2x2 or something larger
        if(this->numRows == 2 && this->numCols == 2)
        {
            //return the difference of the diagonals
            deter = inRep[0][0] * inRep[1][1] - inRep[0][1] * inRep[1][0];
        }
        //calculate the determinant for a larger square matrix
        else if(this->numRows == this->numCols && this->numCols >= 3)
        {
            double sum = 0;
            for(int col = 0;col <numCols;col++)
            {
                //create a copy of the matrix and remove the row and column
                Matrix cofactor(*this);
                cofactor.removeRow(0);
                cofactor.removeCol(col);
                
                //compute the cofactor expansion
                sum += pow(-1,col) * inRep[0][col] * cofactor.determinant();
            }
            deter = sum;
        }
        else
        {
            throw std::invalid_argument("Cannot calculate the determinant of a non square matrix");
        }
        return deter;
    }
    
    
    //returns the transpose of a matrix
    Matrix transpose()
    {
        //decalre an output matrix of size nxm where the input matrix is size mxn
        Matrix outMat(inRep[0].size(),inRep.size());
        
        //loop through the rows of the input matrix
        for(int i = 0;i <inRep.size();i++)
        {
            //set the cooresponding column vector of the output to the row vector of the input
            outMat.insertColVec(Vec(inRep[i]), i);
        }
        return outMat;
    }
    
/*-----EVERYTHING BELOW THIS IS MATRIX IMPLEMENTATION AND NOT SPECIFICALLY RELEVANT TO THE PROJECT----*/

    
    //copy constructor
    Matrix(const Matrix &in):Matrix(in.inRep)
    {}
    //constructor converts a 2d vector into a matrix object
    Matrix(vector<vector<double>> inMat)
    {
        numRows = inMat.size();
        for(int i = 1;i < inRep.size();i++)
        {
            if(inMat[i].size() != inMat[i-1].size())
            {
                throw std::invalid_argument("Matrix has missing elements");
            }
        }
        numCols = inMat[0].size();
        
        inRep.assign(inMat.begin(), inMat.end());
    }
    
    //constructs a mxn zero matrix
    Matrix(int m, int n)
    {
        vector<vector<double>> newMat(m,vector<double>(n,0));
        inRep.assign(newMat.begin(), newMat.end());
        
        numRows = m;
        numCols = n;
    }
    
    //constructs matrix from input file
    Matrix(string inFile, char colDelimiter = ' ', char rowDelimiter = '\n') : Matrix(readMatrixFromFile(inFile, colDelimiter, rowDelimiter))
    {
    }
    
    //overloads the subscript operator
    vector<double>& operator[](int index)
    {
        return inRep[index];
    }

    
    //matrix multiplication
    Matrix operator*(Matrix &B)
    {
        //check if columns of matrix a are equal to rows of matrix b
        if(this->numCols != B.numRows)
        {
            throw std::runtime_error("Unable to multiply matrices: Invalid Dimensions\n");
        }
        
        //declare output matrix with size of Rows in A x Columns in B
        Matrix outMat(this->numRows,B.numCols);
        
        //iterate through matrix
        for(int row=0;row<outMat.numRows;row++)
        {
            for(int col=0;col<outMat.numCols;col++)
            {
                Vec rowVec((*this)[row]);
                Vec colVec(B.inRep,col);
                
                
                outMat[row][col] = rowVec*colVec;
            }
        }
        //return the output matrix
        return outMat;
    }

    //matrix equality
    bool operator==(Matrix &B)
    {
        //different sizes means false
        if(numRows != B.numRows || numCols != B.numCols)
        {
            return false;
        }
        for(int row =0;row<numRows;row++)
        {
            for(int col =0;col<numCols;col++)
            {
                //compare the elements
                if(!doubleEquals(inRep[row][col],B[row][col]))
                {
                    return false;
                }
            }
        }
    return true;
    }
    
    //return's a matrix's column vector
    Vec getColVec(unsigned int col)
    {
        
        return Vec(inRep,col);
    }
    
    //returns a matrix's row vector
    Vec getRowVec(unsigned int row)
    {
        return Vec(inRep[row]);
    }
    
    //insetrs a column vector in the given position
    void insertColVec(Vec colVec,unsigned int col)
    {
        if(col < 0)
        {
            for(int row=0;row<numRows;row++)
            {
                inRep[row].insert(inRep[row].begin(), colVec[row]);
            }
            numCols++;
        }
        if(col >=0 && col < numCols)
        {
            for(int i = 0;i<numRows;i++)
            {
                (*this)[i][col] = colVec[i];
            }
        }
        if(col >= numCols)
        {
            for(int row=0;row<numRows;row++)
            {
                inRep[row].insert(inRep[row].end(), colVec[row]);
            }
            numCols++;
        }
    }
    
    //insetrs a row vector in the given position
    void insertRowVec(Vec rowVec,int row)
    {
        if(row < 0)
        {
            inRep.insert(inRep.begin(),rowVec.inRep);
            numRows++;
        }
        if(row >=0 && row < numRows)
        {
            inRep[row] = rowVec.inRep;
        }
        if(row >= numRows)
        {
            inRep.insert(inRep.end(), rowVec.inRep);
            numRows++;
        }
    }
    
    //removes a column from the matrix
    void removeCol(unsigned int col)
    {
        
        for(int i = 0;i<numRows;i++)
        {
            (*this)[i].erase((*this)[i].begin() + col);
        }
        numCols--;
    }
    
    //removes a row from the matrix
    void removeRow(unsigned int row)
    {
        inRep.erase(inRep.begin() + row);
        numRows--;
    }
    
    //reuturns the matrix as a 2d vector
    vector<vector<double>> getInRep()
    {
        return inRep;
    }
    
    //function to print a matrix with given delimiters to cout, defaults to a space between columns and a newline between rows
    void print(char colDelimiter = ' ',char rowDelimiter = '\n')
    {
        for(auto row:inRep)
        {
            for(auto col:row)
            {
                cout << col << colDelimiter;
            }
            cout << rowDelimiter;
        }
    }
    
    //function to write a matrix to a given output file with given delimiters
    //defaults to space between columns and a line between rows
    void writeToFile(string outFile,bool append = false,char colDelimiter = ' ',char rowDelimiter = '\n')
    {
        ofstream fileStream;
        
        //determine whether to open in append more or not
        std::ios::openmode mode = std::ios::out;
        if(append)
        {
            mode = mode | std::ios::app;
        }
        fileStream.open(outFile,mode);
        fileStream << std::setprecision(4);
        if(fileStream.is_open())
        {
            for(auto row :inRep)
            {
                for(auto col:row)
                {
                    fileStream << col << colDelimiter;
                }
                fileStream << rowDelimiter;
            }
            fileStream.close();
        }
        else
        {
            throw std::runtime_error("Could not open file:" + outFile);
        }
    }

    //static method to return an identity of the given size
    static Matrix identity(int size)
    {
        Matrix outMat(size,size);
        
        for(int i = 0;i<size;i++)
        {
            outMat[i][i] = 1;
        }
        return outMat;
    }
};


//compare 2 doubles with epsilon floating point comaprison 
bool doubleEquals(double d1, double d2)
{
    return abs(d1 - d2) <= 0.00001;
    //return abs(d1 - d2) <= __FLT_EPSILON__;
}

//define scalar multilication by matrix
Matrix operator*(double scalar, Matrix MatA)
{
    for(int row =0;row<MatA.numRows;row++)
    {
        for(int col =0; col<MatA.numCols;col++)
        {
            MatA[row][col] *= scalar;
        }
    }
    return MatA;
}

//define scalar multilication by matrix
Matrix operator/(Matrix MatA,double scalar)
{
    return (1/scalar) * MatA;
}

Matrix operator+(const Matrix& A,const Matrix& B)
{
    if(A.numRows != B.numRows || A.numCols != B.numCols)
    {
        throw std::invalid_argument("Cannot add matrices of different sizes");
    }
    Matrix outputMat(A.numRows,A.numCols);
    
    for(int row = 0;row<A.numRows;row++)
    {
        for(int col = 0;col<A.numCols;col++)
        {
            outputMat[row][col] = A.inRep[row][col] + B.inRep[row][col];
        }
    }
    return outputMat;
}

Matrix operator-(Matrix& A,Matrix& B)
{
    Matrix outMat = A + (-1.0 * B);
    return outMat;
}

vector<vector<double>> readMatrixFromFile(string filename,char colDelimeter,char rowDelimiter)
{
    //vecotr for output
    vector<vector<double>> output;
    
    //open the file in a filestream
    ifstream fileStream;
    fileStream.open(filename);
    if(fileStream.is_open())
    {
        //string to input a row
        string rowString = "";
        //while loop that reads rows into rowString
        while(getline(fileStream,rowString,rowDelimiter))
        {
            //convert rowString to a stream
            stringstream rowStream(rowString);
            
            //declare vector for a roe
            vector<double> row;
            
            //decalre string to hold the column
            string colString = "";
            //read a column from the row
            while(getline(rowStream,colString,colDelimeter))
            {
                //convert the column to a stringstream
                //getline returns a string so we need to convert it to a stringstream in order to read in as double from the string
                stringstream colStream(colString);
                
                //declare a variable for the column value
                double col;
                
                //read in the value from the column stringstream
                colStream >> col;
                
                //push it onto the row
                row.push_back(col);
            }
            //push the row onto the output matrix
            output.push_back(row);
        }
        fileStream.close();
    }
    else
    {
        throw std::runtime_error("Could not open file:" + filename);
    }
    for(int i = 0;i<output.size();)
    {
        if(output[i].size() == 0)
        {
            output.erase(output.begin() + i);
        }
        else
        {
            i++;
        }
    }
    return output;
}

#endif /* Matrix_h */
