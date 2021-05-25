#include<iostream>
#include <itp/getopt>
#include <itp/core>
#include <fstream>

void cvtBinaryToTextFile(std::string inputFileName, std::string outputFileName)
{
    std::ifstream inputFile(inputFileName, std::ios::binary);
    auto begin = inputFile.tellg();
    inputFile.seekg(0, std::ios::end);
    auto length = (inputFile.tellg() - begin) / sizeof(float);
    inputFile.seekg(std::ios::beg);
    std::vector<float> data(length);
    inputFile.read((char*)(data.data()), length * sizeof(float));
    inputFile.close();

    int column = int(std::sqrt(length));
    fmt::print("column: {}\n", column);
    std::ofstream outputFile(outputFileName);
    for (int i = 0; i < column; i++) {
        for (int j = 0; j < column; j++) {
            fmt::print(outputFile, "{:15.8e} ", data[i * column + j]);
        }
        fmt::print(outputFile, "\n");
    }
    outputFile.close();
}

void cvtTextToBinaryFile(std::string inputFileName, std::string outputFileName)
{
    std::string line;
    int column = 0;
    float tmp;
    std::ifstream inputFile(inputFileName);
    std::getline(inputFile, line);
    std::stringstream ss;
    ss.str(line);
    while (ss >> tmp) {
        column++;
    }
    ss.clear();
    std::vector<float> data(column * column);
    int ndx = 0;
    do {
        ss.str(line);
        for (int i = 0; i < column; i++) {
            ss >> data[ndx++];
        }
        ss.clear();
    } while (std::getline(inputFile, line));

    std::ofstream outputFile(outputFileName, std::ios::binary);
    outputFile.write((char*)(data.data()), sizeof(float) * column * column);
    outputFile.close();
}

void main()
{
    // cvtTextToBinaryFile("rMatrix_wall.gro.dat", "rMatrix_wall.gro.bin");
    cvtBinaryToTextFile("rMatrix_wall.gro.bin", "rMatrix_wall.gro.dat1");

}

