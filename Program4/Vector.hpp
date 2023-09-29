//
//  Vector.hpp
//  Program3
//
//  Created by Moses T on 10/28/22.
//
#include <vector>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>

using std::cout,std::vector,std::ofstream,std::endl,std::string,std::ifstream,std::getline,std::stringstream;
#ifndef Vector_h
#define Vector_h

class Vec
{
public:
    vector<double> inRep;
    
    int size;
    //get the magnitude of a vector
    double magnitude()
    {
        double sum = 0;
        for(double element : inRep)
        {
            //sum of squares of elements
            sum+= pow(element,2);
        }
        //squareroot of the sum
        return sqrt(sum);
    }
    
    //normalize the vector
    void normalize()
    {
        double mag = magnitude();
        
        //divide each element by the magnitude
        for(int i=0;i<size;i++)
        {
            inRep[i] = inRep[i]/mag;
        }
    }
    
    //dot product of 2 vectors using the * operator
    double operator*(Vec &vecB)
    {
        if(vecB.size != this->size)
        {
            throw std::invalid_argument("Cannot take dot product of vectors with different sizes");
        }
        //add the prodcucts of the components to the results
        double result = 0;
        for(int i = 0;i<this->size;i++)
        {
            result += vecB[i] * (*this)[i];
        }
        return result;
    }
    
    //vector cross product as overload of the ^ operator
    Vec operator^(Vec &vec2)
    {
        if(vec2.size != 3 || this->size != 3)
        {
            throw std::invalid_argument("Cannot take cross product of 2 non 3d vectors");
        }
        
        //implementation of the cross product formula
        Vec result(size);
        result[0] = (*this)[1] * vec2[2] - (*this)[2] * vec2[1];
        result[1] = -1*((*this)[0] * vec2[2] - (*this)[2] * vec2[0]);
        result[2] = (*this)[0] * vec2[1] - (*this)[1] * vec2[0];
        
        return result;
    }
    
    //vecor subtraction
    Vec operator-(Vec &vecB)
    {
        if(vecB.size != this->size)
        {
            throw std::invalid_argument("Cannot subtract vectors with different sizes");
        }
        Vec outVec(size);
        for(int i=0;i<size;i++)
        {
            outVec[i] = (*this)[i] - vecB[i];
        }
        return outVec;
    }
    
    
    //declare a vector with the given size
    Vec(int inSize)
    {
        size = inSize;
        inRep.resize(inSize, 0);
    }
    //declares a vector from a vector object, can be used to declare row vectors
    Vec(vector<double> inVec)
    {
        inRep.assign(inVec.begin(),inVec.end());
        size = inVec.size();
    }
    
    //declares a vector from a column vector
    Vec(vector<vector<double>> inMat,unsigned int col)
    {
        for(int row = 0;row<inMat.size();row++)
        {
            inRep.push_back(inMat[row][col]);
        }
        size = inMat.size();
    }
    
    double& operator[](int index)
    {
        return inRep[index];
    }
    
    //define scalar multilication by vector
    Vec operator*(double scalar)
    {
        Vec VecA((*this));
        for(int comp =0; comp<VecA.size;comp++)
        {
            VecA[comp] *= scalar;
        }
        return VecA;
    }
    
    //vector addition
    Vec operator+(Vec &vec2)
    {
        //check sizes
        if(vec2.size != this->size)
        {
            throw std::invalid_argument("Cannot add 2 vectors of different size");
        }
        
        //declare an output vector
        Vec result(size);
        
        //set the output vector
        for(int i =0;i<size;i++)
        {
            result[i] = (*this)[i] + vec2[i];
        }
        return result;
    }
    
    void writeToFile(string outFile,char rowDelimiter = '\n',bool append = false)
    {
        ofstream fileStream;
        //determine whether to open in append more or not
        std::ios::openmode mode = std::ios::out;
        if(append)
        {
            mode = mode | std::ios::app;
        }
        fileStream.open(outFile,mode);
        if(fileStream.is_open())
        {
            fileStream << std::setprecision(4);
            for(double num :inRep)
            {
                fileStream << num;
                
                fileStream << rowDelimiter;
            }
            fileStream << endl;
            fileStream.close();
        }
        else
        {
            throw std::runtime_error("Could not open file:" + outFile);
        }
    }

    void print(char rowDelimiter = '\n')
    {
        cout<<std::setprecision(4);
        for(double num :inRep)
        {
            cout << num;
            
            cout << rowDelimiter;
        }
        cout << endl;
    }
};

#endif /* Vector_h */
