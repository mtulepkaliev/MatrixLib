/*
 Moses Tulepkaliev
 Programming Assignment 3
 CS 2300

 */


#include <iostream>
#include "Matrix.hpp"
#include <fstream>
#include <map>
using std::string, std::pair;

// function prototypes
void part1(string inFile, string outFileParalell, string outFilePersepctive);
Vec parallelProject(Vec point, Vec planePoint, Vec planeNormal, Vec projDir);
Vec perspectiveProjection(Vec point, Vec planePoint, Vec planeNormal);

void part2_1(string inFile, string outFile);
double distanceToPlane(Vec point, Vec planePoint, Vec planeNormal);

void part2_2(string inFile, string outFile);
Vec lineTriangleIntersection(Vec linePoint1, Vec linePoint2, Vec vertex1, Vec vertex2, Vec vertex3);
Vec gaussElim(Matrix A, Vec B);

void part3(string inFile, string outFile);
Vec powerMethod(Matrix A, Vec x, double tol, int maxIter);
void writeStringToFile(string outString, string outFile, bool append = false);

int main(int argc, const char *argv[])
{
    part1("inputHW4_Part1_2.txt","mtulepkaliev_output_1_A1.txt","mtulepkaliev_output_1_A2.txt");
    part2_1("inputHW4_Part1_2.txt", "mtulepkaliev_output_1_B1.txt");
    part2_2("inputHW4_Part1_2.txt","mtulepkaliev_output_1_B2.txt");
    part3("InputHW4_Part3.txt","mtulepkaliev_output_1_C.txt");
    return 0;
}

// part 1 function pass in input file, output file for parralell projction, and outputfile for persepctive projection
void part1(string inFile, string outFileParalell, string outFilePersepctive)
{
    // read in the file as a matrix
    Matrix points(inFile);

    // extract 3 vectors for the points and vectors
    Vec planePoint(vector<double>{points[0][0], points[0][1], points[0][2]});
    Vec normalVec(vector<double>{points[0][3], points[0][4], points[0][5]});
    Vec projDir(vector<double>{points[0][6], points[0][7], points[0][8]});

    // remove the top row
    points.removeRow(0);

    // define 2 matries to store the ouputs
    Matrix parProjPoints(points.numRows, 9);
    Matrix perspProjPoints(points.numRows, 9);

    // loop through the rows
    for (unsigned int k = 0; k < points.numRows; k++)
    {
        vector<double> nextRowParallel(0);
        vector<double> nextRowPerspective(0);
        // get 3 points from each line
        for (unsigned int index = 0; index < 9; index += 3)
        {
            // get one point at a time
            Vec point(vector<double>{points[k][index], points[k][index + 1], points[k][index + 2]});

            Vec projectedPoint = parallelProject(point, planePoint, normalVec, projDir);

            // store the projected point in nextRowParallel
            nextRowParallel.insert(nextRowParallel.end(), projectedPoint.inRep.begin(), projectedPoint.inRep.end());

            projectedPoint = perspectiveProjection(point, planePoint, normalVec);

            // store the projected point in nextRowPerspective
            nextRowPerspective.insert(nextRowPerspective.end(), projectedPoint.inRep.begin(), projectedPoint.inRep.end());
        }
        // store the next row in the output matrices
        parProjPoints.insertRowVec(Vec(nextRowParallel), k);
        perspProjPoints.insertRowVec(Vec(nextRowPerspective), k);
    }

    parProjPoints.writeToFile(outFileParalell);
    perspProjPoints.writeToFile(outFilePersepctive);
}

// paralell projects the point onto the plane in the direction of projDir
Vec parallelProject(Vec point, Vec planePoint, Vec planeNormal, Vec projDir)
{
    // normalize the plane normal
    planeNormal.normalize();

    // convert projection direction to 3x1 matrix
    Matrix V(3, 1);
    V.insertColVec(projDir, 0);

    // get the normal direction transpose as a matrix
    Matrix Ntranspose(1, 3);
    Ntranspose.insertRowVec(planeNormal, 0);

    Matrix A(3, 3);
    // A = I_3 - (V * Ntranspose)/(V*N)
    A = Matrix::identity(3) + ((-1 / (projDir * planeNormal)) * (V * Ntranspose));

    Matrix P(3, 1);
    // P = ((Q*N)/(V*N))V
    P.insertColVec(projDir * ((planePoint * planeNormal) / (projDir * planeNormal)), 0);

    // convert the point to a 3x1 matrix named X
    Matrix X(3, 1);
    X.insertColVec(point, 0);

    // X' = AX + P
    Matrix Xprime((A * X + P));

    // convert X' back to a vector
    Vec XprimeVec(Xprime.getColVec(0));

    return XprimeVec;
}

// perspective projects the point onto the plane
Vec perspectiveProjection(Vec point, Vec planePoint, Vec planeNormal)
{
    // define the vector going from the point to the origin
    Vec pointToOrigin = point * -1;

    // run a parallel projection with the vector from the point to the origin as the projection direction.
    return parallelProject(point, planePoint, planeNormal, pointToOrigin);
}

// part 2 subpart 1 function, finds the distance from a given point to a pane in 3D space
void part2_1(string inFile, string outFile)
{
    // read in the file as a matrix
    Matrix inMatrix(inFile);

    Vec outVec(inMatrix.numRows);

    // loop through the rows
    for (int row = 0; row < inMatrix.numRows; row++)
    {
        // get the plane and point
        Vec planePoint(vector<double>{inMatrix[row][0], inMatrix[row][1], inMatrix[row][2]});
        Vec planeNormal(vector<double>{inMatrix[row][3], inMatrix[row][4], inMatrix[row][5]});
        Vec point(vector<double>{inMatrix[row][6], inMatrix[row][7], inMatrix[row][8]});

        // calculate the distance
        outVec[row] = distanceToPlane(point, planePoint, planeNormal);
    }

    // write the output to a file
    outVec.writeToFile(outFile);
}

// calculates the distance from the point to the plane
double distanceToPlane(Vec point, Vec planePoint, Vec planeNormal)
{
    planeNormal.normalize();
    // calculate C coefficient
    // yes i'm taking the dot pruduct of a vector and a point
    // but the point is defined as a vector and the formula is the same
    double c = -1 * (planeNormal * planePoint);

    // t = c + n.p(n is the normal vector and p is the point we are finding the distance of
    // since n is normalized then ||p-q|| = t
    double distance = c + planeNormal * point;

    return abs(distance);
}

// part 2 subpart 2 function intersects a line with triangles and finds intersection point if it exists
void part2_2(string inFile, string outFile)
{
    // read in the file as a matrix
    Matrix triangleMatrix(inFile);

    // get the points that define the line
    Vec linePoint1 = Vec(vector<double>{triangleMatrix[0][0], triangleMatrix[0][1], triangleMatrix[0][2]});
    Vec linePoint2 = Vec(vector<double>{triangleMatrix[0][3], triangleMatrix[0][4], triangleMatrix[0][5]});

    // remove the row that defines the line
    triangleMatrix.removeRow(0);
    writeStringToFile("", outFile);
    // loop through the rows the contain the triangles
    for (int row = 0; row < triangleMatrix.numRows; row++)
    {
        // get the points that define the triangle
        Vec vertex1(vector<double>{triangleMatrix[row][0], triangleMatrix[row][1], triangleMatrix[row][2]});
        Vec vertex2(vector<double>{triangleMatrix[row][3], triangleMatrix[row][4], triangleMatrix[row][5]});
        Vec vertex3(vector<double>{triangleMatrix[row][6], triangleMatrix[row][7], triangleMatrix[row][8]});

        try
        {
            // get the intersection point and write it to the file
            Vec intersectionPoint = lineTriangleIntersection(linePoint1, linePoint2, vertex1, vertex2, vertex3);
            intersectionPoint.writeToFile(outFile, ' ', true);
        }
        // if the triangle does not intersect catch the exception and print it to the file
        catch (const std::exception &e)
        {
            writeStringToFile(e.what(),outFile, true);
        }
    }
}

// calculates the intersection point of a line and a triangle
Vec lineTriangleIntersection(Vec linePoint1, Vec linePoint2, Vec vertex1, Vec vertex2, Vec vertex3)
{
    //define line vector v
    Vec lineVec(linePoint2 - linePoint1);

    Matrix A(3, 3);

    //b = p-p_1
    Vec b(linePoint1 - vertex1);

    //our coefficients are p_2-p_1, p_3-p_1, and -v
    A.insertColVec(Vec(vertex2 - vertex1), 0);
    A.insertColVec(Vec(vertex3 - vertex1), 1);
    A.insertColVec(lineVec * -1, 2);

    // outvec = (u1,u2,t)
    Vec outVec = gaussElim(A, b);

    //make sure u1 and u2 mett conditions for being in triangle, throw exception if not
    if (outVec[0] <= 0 || outVec[1] <= 0 || outVec[0] >= 1 || outVec[1] >= 1 || outVec[0] + outVec[1] > 1)
    {
        throw std::runtime_error("Does not intersect");
    }

    //calculate and return the intersection point
    double t = outVec[2];
    Vec outPoint = lineVec * t;
    outPoint = linePoint1 + outPoint;
    return (outPoint);
}

// implement gaussian elimination
Vec gaussElim(Matrix A, Vec B)
{
    // create augmented matrix
    Matrix aug(A);
    aug.insertColVec(B, A.numCols);

    // implementation of gaussian elimination from notes
    for (int col = 0; col < A.numCols - 1; col++)
    {
        // find the row which has the largest value in the column
        int maxRow = 0;
        for (int row = col; row < A.numRows; row++)
        {
            if (abs(aug[row][col]) > abs(aug[maxRow][col]))
            {
                maxRow = row;
            }
        }
        // swap the rows if needed
        aug.swapRows(col, maxRow);

        // forward elimination loop
        for (int row = col + 1; row < A.numRows; row++)
        {
            double factor = aug[row][col] / aug[col][col];
            aug[row][col] = 0;
            for (int k = col + 1; k < aug.numCols; k++)
            {
                aug[row][k] -= factor * aug[col][k];
            }
        }
    }

    // back substitution
    int n = A.numRows;
    Vec X(n);

    // set X_n
    X[n - 1] = aug[n - 1][n] / aug[n - 1][n - 1];
    for (int j = n - 2; j >= 0; j--)
    {
        // start with b_j
        double sum = aug[j][n];

        // a_j,j+1*u_j+1 + a_j,j+2*u_j+2 + ... + a_j,n*u_n
        for (int col = j + 1; col < n; col++)
        {
            sum -= aug[j][col] * X[col];
        }

        // divide by a_j,j
        X[j] = sum / aug[j][j];
    }

    return X;
}


//part 3 function, calcualtes page ranking given stoachastic google matrix
void part3(string inFile, string outFile)
{
    Matrix stochasticMatrix(inFile);

    // test for invalid input
    try
    {
        for (int col = 0; col < stochasticMatrix.numCols; col++)
        {
            double colSum = 0;
            for (int row = 0; row < stochasticMatrix.numRows; row++)
            {
                // check for negative values
                if (stochasticMatrix[row][col] < 0)
                {
                    throw std::runtime_error("Invalid input");
                }
                colSum += stochasticMatrix[row][col];
            }
            // check for a column sum of 1
            if (doubleEquals(colSum, 1) == false)
            {
                throw std::runtime_error("Invalid input");
            }
        }

        //use the sum of each row as the initial guess for it's eigenvector
        vector<double> inLinkSum(stochasticMatrix.numRows, 0);
        for (int row = 0; row < stochasticMatrix.numRows; row++)
        {
            for (int col = 0; col < stochasticMatrix.numCols; col++)
            {
                inLinkSum[row] += stochasticMatrix[row][col];
            }
        }

        //caculate the eigenvector using the power method
        Vec estimate(inLinkSum);
        Vec eigenvec = powerMethod(stochasticMatrix, estimate, 0.00001, 1000);

        //write the eigenvector to the file
        eigenvec.writeToFile(outFile, ' ');

        //pair together the eigenvector and the indeces 
        vector<pair<double, int>> eigenvecPairs;
        for (int i = 0; i < eigenvec.size; i++)
        {
            eigenvecPairs.push_back(pair<double, int>(eigenvec[i], i));
        }

        //sort by the eigenvalues in reverse order
        std::sort(eigenvecPairs.rbegin(), eigenvecPairs.rend());

        //write the sorted indeces to the file
        Vec indexVec(eigenvec.size);
        for (int i = 0; i < eigenvecPairs.size(); i++)
        {
            indexVec[i] = eigenvecPairs[i].second;
        }
        indexVec.writeToFile(outFile, ' ', true);
    }
    catch (const std::exception &e)
    {
        writeStringToFile(e.what(), outFile);
    }
}

//power method for finding the eigenvector of a matrix
Vec powerMethod(Matrix A, Vec estimate, double tol, int maxIter)
{
    // convert to matrix because it's easier to do matrix operations between 2 matrices, i didn't implment matrix-vector multiplication
    Matrix Rprev(estimate.size, 1);
    Rprev.insertColVec(estimate, 0);

    //calculate r^2
    Matrix R = A * Rprev;

    //divide r by it's infinite norm
    double lambda = R.inf_norm();
    R = R / lambda;
    int iter = 0;

    //keep going until we reached max iteration or we are withing convergence tolerance
    while ((R - Rprev).inf_norm() > tol && iter < maxIter)
    {
        //cauculate new R
        Rprev = R;
        R = A * Rprev;

        //normalize R
        lambda = R.inf_norm();
        R = R / lambda;
        iter++;
    }
    return Vec(R.getColVec(0));
}

//find infinite norm of a Vector(passed in as a 1 column matrix)
inline double Matrix::inf_norm()
{
    //make sure it is a 1 column matrix
    if (numCols != 1)
    {
        throw std::runtime_error("Not a vector as a matrix");
    }
    else
    {
        //find and return the maximum absolute value
        double max = 0;
        for (int row = 0; row < numRows; row++)
        {
            if (abs((inRep[row][0])) > max)
            {
                max = abs(inRep[row][0]);
            }
        }
        return max;
    }
}

//helper function to write a string to a file
void writeStringToFile(string outString, string outFile, bool append)
{
    ofstream fileStream;
    // determine whether to open in append more or not
    std::ios::openmode mode = std::ios::out;
    if (append)
    {
        mode = mode | std::ios::app;
    }
    fileStream.open(outFile, mode);
    if (fileStream.is_open())
    {
        fileStream << outString << endl;
        fileStream.close();
    }
    else
    {
        throw std::runtime_error("Could not open file:" + outFile);
    }
}
