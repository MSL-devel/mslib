#include "HelixGenerator.h"
#include "PDBWriter.h"

int main(int argc, char **argv) {
    HelixGenerator hg;
    AtomVector av;
    int numCAlphas;
    PDBWriter pdbw;
    stringstream ss;

    if(argc < 2) {
        cout << "Usage: testHelixGenerator <num C-alpha atoms> <bbq filename>\n";
        exit(1);
    }

    numCAlphas = atoi(argv[1]);

    if(argc > 2)
        hg.setBBQTableFileName(string(argv[2]));

    hg.generateHelix(av, numCAlphas);

    pdbw.open(ss);
    pdbw.write(av);

    cout << ss.str();

    return 0;
}
