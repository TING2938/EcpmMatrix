#include <filesystem>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <itp/utility>
#include <itp/timer>
#include "fast_math.h"
#include "conp.h"


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

    for (int i = 0; i < natoms; i++) {
        for (int j = 0; j != ratio.rows(); ++j) {
            xall(i + natoms * j, 0) = x(i, 0) + box[0] * ratio(j, 0);
            xall(i + natoms * j, 1) = x(i, 1) + box[1] * ratio(j, 1);
            xall(i + natoms * j, 2) = x(i, 2);
        }
    }

#pragma omp parallel for num_threads(nthread)
    for (int i = 0; i < natoms; i++) {
        for (int j = 0; j != natoms * 9; ++j) {
            double r = std::pow(x(i, 0) - xall(j, 0), 2) +
                std::pow(x(i, 1) - xall(j, 1), 2) +
                std::pow(x(i, 2) - xall(j, 2), 2);
            if ((r < cutoff * cutoff) && r > 1e-7) {
                r = std::sqrt(r);
                rMatrix(i, j % natoms) += (fast_erfc(g_ewald * r) - fast_erfc(eta * r / sqrt(2))) / r;
            }
        }
    }
}

void Conp::calc_rSelf()
{
    for (int i = 0; i != natoms; ++i) {
        rMatrix(i, i) += (std::sqrt(2) * eta - 2 * g_ewald) / std::sqrt(itp::pi);
    }
}

void Conp::calc_rSlab()
{
    for (int i = 0; i != natoms; ++i) {
        for (int j = 0; j != natoms; ++j) {
            if (b3dc) {
                //if system is not slab, no need to add slab correction.
                rMatrix(i, j) += 4 * itp::pi * x(i, 2) * x(j, 2) / volume;
            }
        }
    }
}

void Conp::calc_rKspace()
{
    omp_lock_t lock;
    omp_init_lock(&lock);

#pragma omp parallel for schedule(dynamic,2) num_threads(nthread)
    for (int ii = 0; ii < natoms; ii++) {
        for (int jj = 0; jj <= ii; ++jj) { // loop for each atom
            int k, l, m;
            double sqk;
            double fc; // Fourier coefficient of the Gaussian function used in the Ewald sum
            double dx[3];
            double hash;

            for (m = 0; m < 3; m++) {
                dx[m] = x(ii, m) - x(jj, m);
            }
            hash = calc_hash(dx);
            omp_set_lock(&lock);
            auto it = cache.find(hash);
            if (it != cache.end()) {
                rMatrix(ii, jj) = it->second;
                cacheHitTimes++;
                if (showProgressBar) progressBar.update();
                omp_unset_lock(&lock);
                continue;
            }
            omp_unset_lock(&lock);

            // (k,0,0), (0,l,0), (0,0,m)
            for (m = 1; m <= kmax; m++) {
                sqk = m * m * squnitk[0];
                if (sqk <= gsqmx) {
                    rMatrix(ii, jj) += fcCache1[m - 1][0] * fast_cos(m * unitk[0] * dx[0]);
                }

                sqk = m * m * squnitk[1];
                if (sqk <= gsqmx) {
                    rMatrix(ii, jj) += fcCache1[m - 1][1] * fast_cos(m * unitk[1] * dx[1]);
                }

                sqk = m * m * squnitk[2];
                if (sqk <= gsqmx) {
                    fc = std::exp(sqk * g_ewald_sq_inv) / sqk;
                    rMatrix(ii, jj) += fcCache1[m - 1][2] * fast_cos(m * unitk[2] * dx[2]);
                }
            }

            // 1 = (k,l,0), 2 = (k,-l,0)

            for (k = 1; k <= kxmax; k++) {
                for (l = 1; l <= kymax; l++) {
                    sqk = squnitk[0] * k * k + squnitk[1] * l * l;
                    if (sqk <= gsqmx) {
                        rMatrix(ii, jj) += fcCache2[k][l][0] *
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
                        rMatrix(ii, jj) += fcCache2[0][l][m] *
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
                        rMatrix(ii, jj) += fcCache2[k][0][m] *
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
                            rMatrix(ii, jj) += fcCache2[k][l][m] *
                                fast_cos(k * unitk[0] * dx[0]) *
                                fast_cos(l * unitk[1] * dx[1]) *
                                fast_cos(m * unitk[2] * dx[2]);
                        }
                    }
                }
            }
            omp_set_lock(&lock);
            cache[hash] = rMatrix(ii, jj);
            if (showProgressBar) progressBar.update();
            omp_unset_lock(&lock);
        }
    }

    omp_destroy_lock(&lock);

#pragma omp parallel for num_threads(nthread)
    for (int ii = 0; ii < natoms; ii++) {
        for (int jj = 0; jj < ii; jj++) {
            rMatrix(jj, ii) = rMatrix(ii, jj);
        }
    }
    rMatrix *= (32.0 * itp::pi / volume);
}

void Conp::calcKspaceFcCache()
{
    int k, l, m;
    double sqk;

    // (k,0,0), (0,l,0), (0,0,m)
    for (m = 1; m <= kmax; m++) {
        sqk = m * m * squnitk[0];
        if (sqk <= gsqmx) {
            fcCache1[m - 1][0] = std::exp(sqk * g_ewald_sq_inv) / sqk / 4;
        }

        sqk = m * m * squnitk[1];
        if (sqk <= gsqmx) {
            fcCache1[m - 1][1] = std::exp(sqk * g_ewald_sq_inv) / sqk / 4;
        }

        sqk = m * m * squnitk[2];
        if (sqk <= gsqmx) {
            fcCache1[m - 1][2] = std::exp(sqk * g_ewald_sq_inv) / sqk / 4;
        }
    }

    // 1 = (k,l,0), 2 = (k,-l,0)

    for (k = 1; k <= kxmax; k++) {
        for (l = 1; l <= kymax; l++) {
            sqk = squnitk[0] * k * k + squnitk[1] * l * l;
            if (sqk <= gsqmx) {
                fcCache2[k][l][0] = std::exp(sqk * g_ewald_sq_inv) / sqk / 2;
            }
        }
    }

    // 1 = (0,l,m), 2 = (0,l,-m)
    for (m = 1; m <= kzmax; m++) {
        for (l = 1; l <= kymax; l++) {
            sqk = squnitk[1] * l * l + squnitk[2] * m * m;
            if (sqk <= gsqmx) {
                fcCache2[0][l][m] = std::exp(sqk * g_ewald_sq_inv) / sqk / 2;
            }
        }
    }

    // 1 = (k,0,m), 2 = (k,0,-m)
    for (m = 1; m <= kzmax; m++) {
        for (k = 1; k <= kxmax; k++) {
            sqk = squnitk[0] * k * k + squnitk[2] * m * m;
            if (sqk <= gsqmx) {
                fcCache2[k][0][m] = std::exp(sqk * g_ewald_sq_inv) / sqk / 2;
            }
        }
    }

    // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)
    for (m = 1; m <= kzmax; m++) {
        for (k = 1; k <= kxmax; k++) {
            for (l = 1; l <= kymax; l++) {
                sqk = squnitk[0] * k * k + squnitk[1] * l * l + squnitk[2] * m * m;
                if (sqk <= gsqmx) {
                    fcCache2[k][l][m] = std::exp(sqk * g_ewald_sq_inv) / sqk;
                }
            }
        }
    }
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
    nthread = omp_get_max_threads();
    getopt(nthread, "-thread", false, "number of thread");
    getopt(cutoff, "-cutoff", false, "cut off(nm)");
    getopt(binary, "-binary", false, "output binary file ?");
    getopt(calcInv, "-calcInv", false, "calc inverse ?");
    getopt(eta, "-eta", false, "eta");
    getopt(g_ewald, "-ewald", false, "g_ewald");
    getopt(progressBarStyle, "-style", false, "progress bar style");
    getopt(showProgressBar, "-showBar", false, "show progress bar ?");
    getopt.finish();

    // print input args
    fmt::print("\ninput args:\n");
    // fmt::print(".gro file name: {}\n", fnm);
    fmt::print("number of total electrode atoms: {}\n", natoms);
    fmt::print("kxmax: {}\n", kxmax);
    fmt::print("kymax: {}\n", kymax);
    fmt::print("kzmax: {}\n", kzmax);
    fmt::print("b3dc: {}\n", b3dc);
    fmt::print("output binary file: {}\n", binary);
    fmt::print("number of thread: {}\n", nthread);
    fmt::print("cutoff: {}\n", cutoff);

    readGro(fnm[0]);
    itp::Timeit timer;
    timer.start();
    load_cache();
    timer.stop();
    timer.printSpan("load cache: ", "s\n");

    cutoff *= 10; // transfer nm to A;
    g_ewald_sq_inv = -1.0 / (g_ewald * g_ewald) / 4;
    kmax = std::max({ kxmax, kymax, kzmax });
    for (int i = 0; i < 3; i++) {
        unitk[i] = 2.0 * itp::pi / box[i];
        squnitk[i] = unitk[i] * unitk[i];
    }
    double gsqxmx = squnitk[0] * kxmax * kxmax;
    double gsqymx = squnitk[1] * kymax * kymax;
    double gsqzmx = squnitk[2] * kzmax * kzmax;
    gsqmx = std::max({ gsqxmx, gsqymx, gsqzmx }) * 1.00001;

    fcCache1.resize(kmax);
    fcCache2.resize(kxmax + 1);
    for (auto&& y : fcCache2) {
        y.resize(kymax + 1);
        for (auto&& z : y) {
            z.resize(kzmax + 1, 0);
        }
    }
    timer.start();
    calcKspaceFcCache();
    timer.stop();
    timer.printSpan("calc kspace fc cache: ", "s\n\n");

    for (auto&& fileName : fnm) {
        readGro(fileName);
        fmt::print("cache size: {}\n", cache.size());
        cacheHitTimes = 0;

        timer.start();
        rMatrix.resize(natoms, natoms);
        rMatrix.fill(0);

        // calculate matrix
        {
            if (showProgressBar) {
                progressBar.ncols = 30;
                progressBar.totalNum = natoms * (natoms + 1) / 2;
                progressBar.style = progressBarStyle;
                progressBar.start();
            }
            itp::Timeit timer;
            timer.start();
            calc_rKspace();
            timer.stop();
            if (showProgressBar) progressBar.finish();
            timer.printSpan("calc_kspace: ", "s\n");

            timer.start();
            calc_rReal();
            timer.stop();
            timer.printSpan("calc_real: ", "s\n");

            timer.start();
            calc_rSelf();
            timer.stop();
            timer.printSpan("calc_self: ", "s\n");

            timer.start();
            calc_rSlab();
            timer.stop();
            timer.printSpan("calc_slab: ", "s\n");
        }
        timer.stop();
        fmt::print("[{}] *** calculate matrix spend time: {} s\n", fileName, timer.span());
        fmt::print("cache hit times: {}\n", cacheHitTimes);
        fmt::print("cache hit rate: {:4.2f} %\n", double(cacheHitTimes) / (natoms * (natoms + 1)) * 200);

        if (calcInv) {
            timer.start();
            rInvMatrix = rMatrix.matrix().inverse().array();
            timer.stop();
            timer.printSpan("calculate inv: ", "s\n");
        }

        timer.start();
        if (binary) {
            saveToBinaryFile("rMatrix_" + fileName + ".bin", rMatrix);
            if (calcInv) saveToBinaryFile("rInvMatrix_" + fileName + ".bin", rInvMatrix);
        } else {
            saveToTextFile("rMatrix_" + fileName + ".dat", rMatrix, "{:15.8e} ");
            if (calcInv) saveToTextFile("rInvMatrix_" + fileName + ".dat", rInvMatrix, "{:15.8e} ");
        }
        timer.stop();
        timer.printSpan("save matrix spend time: ", "s\n\n");
    }
    timer.start();
    save_cache();
    timer.stop();
    timer.printSpan("save cache: ", "s\n");
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

void Conp::get_CpmControlFile()
{
    int NelectrodeAtom = 10800;
    int NelectrodeAtomOuter = 2700;
    int NelectrodeAtomInner = 2700;
    int CalcQfreq = 1;
    int OutQfreq = 2000;

    itp::Getopt getopt(argc, argv);
    getopt(NelectrodeAtom, "-n", true, "number of total electrode atoms");
    getopt(NelectrodeAtom, "-outer", true, "number of outer electrode atoms");
    getopt(NelectrodeAtom, "-inner", true, "number of inner electrode atoms");
    getopt(CalcQfreq, "-calcQ", false, "CalcQ freq");
    getopt(OutQfreq, "-outQ", false, "OutQ freq");
    getopt.finish();

    std::ofstream ofile("CPM_control.dat");
    fmt::print(ofile, "NelectrodeAtom  {}\n", NelectrodeAtom);
    fmt::print(ofile, "NelectrodeAtomOuter  {}\n", NelectrodeAtomOuter);
    fmt::print(ofile, "NelectrodeAtomInner  {}\n", NelectrodeAtomInner);
    fmt::print(ofile, "CalcQfreq  {}\n", CalcQfreq);
    fmt::print(ofile, "OutQfreq   {}\n", OutQfreq);
    ofile.close();
}

void Conp::cvtBinaryToTextFile()
{
    std::string inputFileName = "rMatrix.bin";
    std::string outputFileName = "rMatrix.dat";
    itp::Getopt getopt(argc, argv);
    getopt(inputFileName, "-f", true, "input binary file name");
    getopt(outputFileName, "-o", true, "output text file name");
    getopt.finish();

    std::ifstream inputFile(inputFileName, std::ios::binary);
    auto begin = inputFile.tellg();
    inputFile.seekg(std::ios::end);
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

void Conp::cvtTextToBinaryFile()
{
    std::string inputFileName = "rMatrix.dat";
    std::string outputFileName = "rMatrix.bin";
    itp::Getopt getopt(argc, argv);
    getopt(inputFileName, "-f", true, "input text file name");
    getopt(outputFileName, "-o", true, "output binary file name");
    getopt.finish();

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
    std::ifstream ifile("EcpmMatrixCache/kspaceCache_{:.3f}.bin"_format(calc_hash(box)), std::ios::binary);
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
    std::ofstream ofile("EcpmMatrixCache/kspaceCache_{:.3f}.bin"_format(calc_hash(box)), std::ios::binary);
    int length = int(cache.size());
    ofile.write((char*)(&length), sizeof(int));
    for (auto it = cache.begin(); it != cache.end(); it++) {
        ofile.write((char*)(&(*it)), sizeof(double) * 2);
    }
}

