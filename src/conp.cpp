#include "conp.h"
#include <fstream>
#include <chrono>
#include <iomanip>
#include <itp/utility>
#include <itp/timer>
#include <filesystem>
#include "fast_math.h"

std::unordered_map<double, double> Conp::cache = {};

Conp::Conp(int gc, char** gv) : argc(gc), argv(gv)
{}

void Conp::readGro(const std::string& fileName)
{
    // load .gro file;
    std::ifstream file(fileName);
    std::stringstream ss;
    std::string line;
    int allatoms;

    std::getline(file, line);
    std::getline(file, line);

    ss.str(line);
    ss >> allatoms;
    ss.clear();

    x.resize(natoms, 3);
    for (int i = 0; i != natoms; ++i) {
        std::getline(file, line);
        ss.str(line.substr(20, 24));
        ss >> x(i, 0) >> x(i, 1) >> x(i, 2);
        ss.clear();
    }

    for (int i = 0; i != allatoms - natoms; ++i) {
        std::getline(file, line);
    }

    std::getline(file, line);
    ss.str(line);
    ss >> box[0] >> box[1] >> box[2];
    ss.clear();
    for (auto&& i : box) {
        i *= 10;
    }

    x *= 10;
    for (int i = 0; i != natoms; ++i) {
        x(i, 2) -= box[2] / 2;
    }
    volume = box[0] * box[1] * box[2];
}

void Conp::calc_rReal()
{
    Eigen::ArrayX3d xall(natoms * 9, 3);
    Eigen::ArrayXXd ratio(9, 2);
    ratio << 0, 0, 0, 1, 0, -1, 1, 0, 1, 1, 1, -1, -1, 0, -1, 1, -1, -1;
    for (int i = 0; i != natoms; ++i) {
        for (int j = 0; j != ratio.rows(); ++j) {
            xall(i + natoms * j, 0) = x(i, 0) + box[0] * ratio(j, 0);
            xall(i + natoms * j, 1) = x(i, 1) + box[1] * ratio(j, 1);
            xall(i + natoms * j, 2) = x(i, 2);
        }
    }

    double r;
    for (int i = 0; i != natoms; ++i) {
        for (int j = 0; j != natoms * 9; ++j) {
            r = std::pow(x(i, 0) - xall(j, 0), 2) +
                std::pow(x(i, 1) - xall(j, 1), 2) +
                std::pow(x(i, 2) - xall(j, 2), 2);
            if ((r < cutoff * cutoff) && r > 1e-7) {
                r = std::sqrt(r);
                rReal(i, j % natoms) = (fast_erfc(g_ewald * r) - fast_erfc(eta * r / sqrt(2))) / r;
            }
        }
    }
}

void Conp::calc_rSelf()
{
    for (int i = 0; i != natoms; ++i) {
        rSelf(i, i) = (std::sqrt(2) * eta - 2 * g_ewald) / std::sqrt(itp::pi);
    }
}

void Conp::calc_rSlab()
{
    for (int i = 0; i != natoms; ++i) {
        for (int j = 0; j != natoms; ++j) {
            if (b3dc) {
                //if system is not slab, no need to add slab correction.
                rSlab(i, j) = 4 * itp::pi * x(i, 2) * x(j, 2) / volume;
            }
        }
    }
}

void Conp::calc_rKspace()
{
    for (int i = 0; i < 3; i++) {
        unitk[i] = 2.0 * itp::pi / box[i];
        squnitk[i] = unitk[i] * unitk[i];
    }

    kmax = std::max({ kxmax, kymax, kzmax});

    double gsqxmx = squnitk[0] * kxmax * kxmax;
    double gsqymx = squnitk[1] * kymax * kymax;
    double gsqzmx = squnitk[2] * kzmax * kzmax;
    gsqmx = std::max({ gsqxmx, gsqymx, gsqzmx }) * 1.00001;

    g_ewald_sq_inv = -1.0 / (g_ewald * g_ewald) / 4;
    
    fmt::print("cache size: {}\n", cache.size());

    itp::Timer timer;
    timer.start();
    for (int i = 0; i != nthread; ++i) {
        thread.emplace_back(&Conp::kspaceThreadFunc, this, i);
    }
    for (auto&& i : thread) {
        if (i.joinable())
            i.join();
    }
    timer.stop();
    fmt::print("calc_kspace: {} s\n", timer.span());

    for (int ii = 0; ii < natoms; ii++) {
        for (int jj = 0; jj < ii; jj++) {
            rKspace(jj, ii) = rKspace(ii, jj);
        }
    }
    rKspace *= (32.0 * itp::pi / volume);
}

void Conp::kspaceThreadFunc(int rank)
{
    int k, l, m;
    double sqk;
    double fc; // Fourier coefficient of the Gaussian function used in the Ewald sum
    double dx[3];
    double hash;

    for (int ii = rank; ii < natoms; ii += nthread) 
    {
        for (int jj = 0; jj <= ii; ++jj) 
        {
            for (m = 0; m < 3; m++) {
                dx[m] = x(ii, m) - x(jj, m);
            }
            hash = calc_hash(dx);
            mut.lock();
            auto it = cache.find(hash);
            if (it != cache.end()) {
                rKspace(ii, jj) = it->second;
                cacheHit++;
                mut.unlock();
                continue;
            }
            mut.unlock();

            // (k,0,0), (0,l,0), (0,0,m)
            for (m = 1; m <= kmax; m++) {
                sqk = m * m * squnitk[0];
                if (sqk <= gsqmx) {
                    fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                    rKspace(ii, jj) += fc / 4 * fast_cos(m * unitk[0] * dx[0]);
                }

                sqk = m * m * squnitk[1];
                if (sqk <= gsqmx) {
                    fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                    rKspace(ii, jj) += fc / 4 * fast_cos(m * unitk[1] * dx[1]);
                }

                sqk = m * m * squnitk[2];
                if (sqk <= gsqmx) {
                    fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                    rKspace(ii, jj) += fc / 4 * fast_cos(m * unitk[2] * dx[2]);
                }
            }

            // 1 = (k,l,0), 2 = (k,-l,0)

            for (k = 1; k <= kxmax; k++) {
                for (l = 1; l <= kymax; l++) {
                    sqk = squnitk[0] * k * k + squnitk[1] * l * l;
                    if (sqk <= gsqmx) {
                        fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                        rKspace(ii, jj) += fc / 2 *
                            fast_cos(k * unitk[0] * dx[0]) *
                            fast_cos(l * unitk[1] * dx[1]);
                    }
                }
            }

            // 1 = (0,l,m), 2 = (0,l,-m)
            for (m = 1; m <= kzmax; m++) {
                for (l = 1; l <= kymax; l++) {               
                    sqk = squnitk[1] * l * l + squnitk[2] * m * m;
                    if (sqk <= gsqmx) {
                        fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                        rKspace(ii, jj) += fc / 2 *
                            fast_cos(l * unitk[1] * dx[1]) *
                            fast_cos(m * unitk[2] * dx[2]);
                    }
                }
            }

            // 1 = (k,0,m), 2 = (k,0,-m)
            for (m = 1; m <= kzmax; m++) {
                for (k = 1; k <= kxmax; k++) {                
                    sqk = squnitk[0] * k * k + squnitk[2] * m * m;
                    if (sqk <= gsqmx) {
                        fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                        rKspace(ii, jj) += fc / 2 *
                            fast_cos(k * unitk[0] * dx[0]) *
                            fast_cos(m * unitk[2] * dx[2]);
                    }
                }
            }

            // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
            for (m = 1; m <= kzmax; m++) {
                for (k = 1; k <= kxmax; k++) {
                    for (l = 1; l <= kymax; l++) {                   
                        sqk = squnitk[0] * k * k + squnitk[1] * l * l + squnitk[2] * m * m;
                        if (sqk <= gsqmx) {
                            fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                            rKspace(ii, jj) += fc *
                                fast_cos(k * unitk[0] * dx[0]) *
                                fast_cos(l * unitk[1] * dx[1]) *
                                fast_cos(m * unitk[2] * dx[2]);
                        }
                    }
                }
            }
            mut.lock();
            cache[hash] = rKspace(ii, jj);
            mut.unlock();
        }
    }
}

void Conp::calc_Matrix()
{
    rMatrix.resize(natoms, natoms);
    rReal.resize(natoms, natoms);
    rSelf.resize(natoms, natoms);
    rSlab.resize(natoms, natoms);
    rKspace.resize(natoms, natoms);

    rMatrix.fill(0);
    rReal.fill(0);
    rSelf.fill(0);
    rSlab.fill(0);
    rKspace.fill(0);

    calc_rKspace();
    calc_rReal();
    calc_rSelf();
    calc_rSlab();

    rMatrix = rReal + rSlab + rSelf + rKspace;
}

void Conp::get_EcpmMatrix()
{
    // get comment line args;
    itp::Getopt getopt(argc, argv, "calculate Matrix for ECPM\n");

    getopt.getArray(fnm, "-f", true, "gro file name");
    getopt(natoms, "-n", true, "number of total electrode atoms");
    getopt(kxmax, "-kx", true, "kxmax");
    getopt(kymax, "-ky", true, "kymax");
    getopt(kzmax, "-kz", true, "kzmax");
    getopt(b3dc, "-3dc", true, "3d(0) or 3dc(1)");

    // optional 
    nthread = std::thread::hardware_concurrency();
    getopt(nthread, "-thread", false, "number of thread");
    getopt(cutoff, "-cutoff", false, "cut off(nm)");
    getopt(binary, "-binary", false, "output binary file ?");
    getopt(calcInv, "-calcInv", false, "calc inverse ?");
    getopt(eta, "-eta", false, "eta");
    getopt(g_ewald, "-ewald", false, "g_ewald");
    getopt.finish();

    // print input args
    fmt::print("\ninput args:\n");
    // fmt::print(".gro file name: {}\n", fnm);
    fmt::print("number of total electrode atoms: {}\n", natoms);
    fmt::print("kxmax: {}\n", kxmax);
    fmt::print("kymax: {}\n", kymax);
    fmt::print("kzmax: {}\n", kzmax);
    fmt::print("3d(0) or 3dc(1): {}\n", b3dc);
    fmt::print("output binary file: {}\n", binary);
    fmt::print("number of thread: {}\n", nthread);
    fmt::print("cutoff: {}\n", cutoff);

    itp::Timer timer;
    timer.start();
    load_cache();
    timer.stop();
    fmt::print("load_cache: {} s\n\n", timer.span());

    cutoff *= 10; // transfer nm to A;

    for (auto&& fileName : fnm) {
        readGro(fileName);

        cacheHit = 0;
        timer.start();
        calc_Matrix();
        timer.stop();
        fmt::print("{} *** calculate matrix spend time: {} s\n", fileName, timer.span());
        fmt::print("cache hit count: {}\n", cacheHit);
        fmt::print("cache hit rate: {:4.2f} %\n\n", double(cacheHit) / (natoms * (natoms+1)) * 200);

        if (calcInv)
            rInvMatrix = rMatrix.matrix().inverse().array();

        if (binary) {
            saveToBinaryFile("rMatrix_{}.bin"_format(fileName), rMatrix);
            if (calcInv) saveToBinaryFile("rInvMatrix_{}.bin"_format(fileName), rInvMatrix);
        } else {
            saveToTextFile("rMatrix_{}.dat"_format(fileName), rMatrix, "{:15.8e} ");
            if (calcInv) saveToTextFile("rInvMatrix_{}.dat"_format(fileName), rInvMatrix, "{:15.8e} ");
        }
    }
    timer.start();
    save_cache();
    timer.stop();
    fmt::print("save_cache: {} s\n", timer.span());
}

void Conp::get_getPotFile()
{
    // get comment line args;
    itp::Getopt getopt(argc, argv, "get \"getPotFile\" for ECPM\n");

    int allatoms = 0;
    std::vector<double> voltage = { 0, 1, 2 };
    double lowPos = 0.0;
    double upPos = 10.0;

    getopt.getArray(fnm, "-f", true, "gro file name");
    getopt(natoms, "-n", true, "number of total electrode atoms");
    getopt(allatoms, "-N", true, "number of total atoms of system");
    getopt.getArray(voltage, "-v", true, "voltage");
    getopt(lowPos, "-low", true, "bulk low position(z)(nm)");
    getopt(upPos, "-up", true, "bulk up position(z)(nm)");

    getopt(cutoff, "-cutoff", false, "cut off(nm)");

    getopt.finish();

    readGro(fnm[0]);

    for (auto&& v : voltage) {
        std::ofstream file("getPot_parameters_{}V.dat"_format(v));
        fmt::print(file, "NAllatom: {:8d}\n", allatoms);
        fmt::print(file, "rcutoff: {:8.3f}\n", cutoff);
        fmt::print(file, "freqCal: {:8d}\n", 1);
        fmt::print(file, "freqOut: {:8d}\n", 2000);
        fmt::print(file, "boundraylow1: {:10.4f}\n", lowPos);
        fmt::print(file, "boundrayup1: {:10.4f}\n", upPos);
        fmt::print(file, "boundraylow2: {:10.4f}\n", lowPos);
        fmt::print(file, "boundrayup2: {:10.4f}\n", upPos);
        fmt::print(file, "NUserSelectGrid: {:8d}\n", natoms);
        for (int i = 0; i != natoms / 2; ++i) {
            fmt::print(file, "{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:8d}\n", x(i, 0) / 10, x(i, 1) / 10,
                (x(i, 2) + box[2] / 2) / 10, v / 2, i);
        }
        for (int i = natoms / 2; i != natoms; ++i) {
            fmt::print(file, "{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:8d}\n", x(i, 0) / 10, x(i, 1) / 10,
                (x(i, 2) + box[2] / 2) / 10, -v / 2, i);
        }
        file.close();
    }
}

void Conp::saveToTextFile(std::string fnm, const Eigen::ArrayXXd& data, std::string fmt)
{
    std::ofstream ofile(fnm);
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            fmt::print(ofile, fmt, data(i, j));
        }
        fmt::print(ofile, "\n");
    }
    ofile.close();
}

void Conp::saveToBinaryFile(std::string fnm, const Eigen::ArrayXXd& data)
{
    Eigen::ArrayXXf outputData = data.cast<float>();
    std::ofstream ofile(fnm, std::ios::binary);
    ofile.write((char*)(outputData.data()), sizeof(float) * outputData.size());
    ofile.close();
}

double Conp::calc_hash(double dx[3])
{   
    return std::abs(dx[0]) * 800043.022364 + std::abs(dx[1]) * 4400.3234231 + std::abs(dx[2]) * 103.9262402;
}

void Conp::load_cache()
{
    std::ifstream ifile("EcpmMatrixCache/kspaceCache.bin", std::ios::binary);
    if (ifile) {
        fmt::print("Success of load cache file.\n");
        int length = 0;
        double tmp[2];
        ifile.read((char*)(&length), sizeof(int));
        for (int i = 0; i < length; i++) {
            ifile.read((char*)(tmp), 2 * sizeof(double));
            cache[tmp[0]] = tmp[1];
        }
    }
}

void Conp::save_cache()
{
    if (!std::filesystem::is_directory("EcpmMatrixCache")) {
        std::filesystem::create_directory("EcpmMatrixCache");
    }
    std::ofstream ofile("EcpmMatrixCache/kspaceCache.bin", std::ios::binary);
    int length = cache.size();
    ofile.write((char*)(&length), sizeof(int));
    for (auto it = cache.begin(); it != cache.end(); it++) {
        ofile.write((char*)(&(*it)), sizeof(double) * 2);
    }
}
