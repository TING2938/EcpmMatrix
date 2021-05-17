﻿#ifndef __CONP_H__
#define __CONP_H__

#include <itp/core>
#include <itp/getopt>
#include <thread>
#include <mutex>
#include <unordered_map>

class Conp
{
public:
    Conp() = default;

    Conp(int gc, char** gv);

    void readGro(const std::string& fileName);

    void calc_rReal();

    void calc_rSelf();

    void calc_rSlab();

    void calc_rKspace();

    void calc_Matrix();

    void get_EcpmMatrix();

    void get_getPotFile();

    void get_CpmControlFile();

    void cvtBinaryToTextFile();

    void cvtTextToBinaryFile();

public:
    void kspaceThreadFunc(int rank);

    void calcKspaceFcCache();

    void saveToTextFile(std::string fnm, const Eigen::ArrayXXd& data, std::string fmt = "{:15.8e} ");

    void saveToBinaryFile(std::string fnm, const Eigen::ArrayXXd& data);

    double calc_hash(double dx[3]);

    void load_cache();

    void save_cache();

    void printProgressInfo();

public:
    int argc;
    char** argv;

    std::vector<std::string> fnm = { "wall.gro" };
    double eta = 1.979;
    double g_ewald = 0.2602844;
    double cutoff = 1.2; // nm;
    int kxmax = 1;
    int kymax = 1;
    int kzmax = 1;
    int natoms = 0;
    int nthread = 16;
    double volume;
    double box[3];
    bool b3dc = true;
    bool binary = false;
    bool calcInv = true;

    Eigen::ArrayX3d x;
    Eigen::ArrayXXd rReal, rSelf, rSlab, rKspace, rMatrix, rInvMatrix;
    std::vector<std::thread> thread;

    double unitk[3], squnitk[3];
    int kmax;
    double gsqmx;
    double g_ewald_sq_inv;
    std::mutex mut;
    int cacheHitTimes = 0;
    size_t progressCount;

    static std::unordered_map<double, double> cache;
    // Fourier coefficient of the Gaussian function used in the Ewald sum
    std::vector<std::array<double, 3>> fcCache1;
    std::vector<std::vector<std::vector<double>>> fcCache2;  
};


#endif // !__CONP_H__