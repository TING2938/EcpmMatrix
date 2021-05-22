#include <map>
#include <itp/core>
#include <itp/color>
#include <itp/timer>
using std::map;


bool findDict(const map<double, map<double, map<double, double>>>& dict, double x, double y, double z)
{
    auto dx = dict.find(x);
    if (dx != dict.end()) {
        auto dy = dx->second.find(y);
        if (dy != dx->second.end()) {
            auto dz = dy->second.find(z);
            if (dz != dy->second.end()) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
}


int main()
{
    map<double, map<double, map<double, double>>> dict;


    dict[2.3][4.5][6.9] = 3.56;

    bool a = findDict(dict, 2.3, 4.5, 6.93);
    bool b = findDict(dict, 2.3, 4.5, 6.9);
    // fmt::print("{} {} {}\n", findDict(dict, 2.3, 4.5, 6.93), a, b);
    fmt::print("{} {}\n", findDict(dict, 2.3, 4.5, 6.93), dict[2.3][4.5][6.93]);
}
