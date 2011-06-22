#include "System.h"
#include "BBQTable.h"

using namespace MSL;

int main() {
    System sys;
    // Read in a pdb file that only includes C-alpha atoms.
    sys.readPdb("../exampleFiles/example0009_caOnly.pdb");
    BBQTable bbq("../tables/PiscesBBQTable.txt");

    // Now fill in the missing bacbone atoms for each chain
    for(int chainNum = 0; chainNum < sys.chainSize(); ++ chainNum) {
        bbq.fillInMissingBBAtoms(sys.getChain(chainNum));
    }

    // Now output a pdb with all of the backbone atoms.
    // Note: Due to the way the BBQ algorithm works, no backbone
    // atoms will be generated for the first and last resiude in a chain.
    sys.writePdb("output.pdb");
}
