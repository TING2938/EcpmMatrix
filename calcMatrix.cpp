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

    subprogram.addSubProgram("potFile", "Get \"getPot_parameters.dat\" file for different voltage",
        [argc, argv] {
            Conp conp(argc, argv);
            conp.get_getPotFile();
        });

    subprogram.addSubProgram("cpmControlFile", "Get \"CPM_control.dat\" file",
        [argc, argv] {
            Conp conp(argc, argv);
            conp.get_CpmControlFile();
        });

    subprogram.addSubProgram("cvtBinaryToText", "convert binary file to text file",
        [argc, argv] {
            Conp conp(argc, argv);
            conp.cvtBinaryToTextFile();
        });

    subprogram.addSubProgram("cvtTextToBinary", "convert text file to binary file",
        [argc, argv] {
            Conp conp(argc, argv);
            conp.cvtTextToBinaryFile();
        });

    subprogram.finish();
}
