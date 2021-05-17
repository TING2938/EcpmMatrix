#include "src/conp.h"
#include <itp/timer>

int main(int argc, char** argv)
{
    itp::Getopt subprogram(argc, argv);

    subprogram.addSubProgram("matrix", "Calculate Matrix for ECPM",
        [argc, argv] {
            Conp conp(argc, argv);
            conp.get_EcpmMatrix();
        });

    subprogram.addSubProgram("potFile", "Get \"getPot_parameters.dat\" for different voltage",
        [argc, argv] {
            Conp conp(argc, argv);
            conp.get_getPotFile();
        });

    subprogram.finish();
}
