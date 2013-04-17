#include "pdb2crd.h"

int main(int argc, char *argv[]) {

	if (argc == 1) {
		usage();
		return 0;
	}

	Options opt = parseOptions(argc, argv);
	if(opt.errorFlag) {
		cout << opt.OPerrors ;
		cout << opt.errorMessages ;
		exit(0);
	}

	if (opt.help) {
		help();
		return 0;
	}

	/********************************************************
	 * Checkt the options
	 ********************************************************/
	bool writeCrd_flag = true;
	
	if (MslTools::pathExtension(opt.input) != "pdb" && MslTools::pathExtension(opt.input) != "crd") {
		cerr << "ERROR: input file format not recognized (not a pdb or crd)" << endl;
		exit(1);
	}

	if (opt.output == "") {
		// create output name from input
		if (MslTools::pathExtension(opt.input) == "pdb") {
			opt.output = MslTools::pathRoot(opt.input) + ".crd";
		} else if (MslTools::pathExtension(opt.input) == "crd") {
			opt.output = MslTools::pathRoot(opt.input) + "-crd.crd";
		}
	} else {
		if (MslTools::pathExtension(opt.output) == "pdb") {
			// write in pdb format
			writeCrd_flag = false;
		} else if (MslTools::pathExtension(opt.output) == "crd") {
			// write in crd format
			writeCrd_flag = true;
		} else {
			// error
			cerr << "ERROR: output file format not recognized (not a pdb or crd)" << endl;
			exit(1);
		}
	}


	System structure;
	if (!structure.readPdb(opt.input)) {
		cerr << "ERROR: cannot read input file " << opt.input << endl;
	}

	cout << "Converting " << opt.input << " to " << opt.output;
	if (!opt.noNameTranslation) {
		string inputFormat = "PDB3";
		if (opt.pdb2) {
			inputFormat = "PDB2";
		}
		cout << ", with translation from " << inputFormat;
		string outputFormat = "CHARMM22";
		if (opt.charmm19) {
			outputFormat = "CHARMM19";
		}
		cout << " to " << outputFormat << endl;
		FormatConverter fc(inputFormat, outputFormat);
		fc.convert(structure.getAtomPointers());
	} else {
		cout << ", without name translation)" << endl;
	}
	if (writeCrd_flag) {
		//structure.writeCrd(opt.output);
		CRDWriter writer;
		writer.open(opt.output);
		writer.write(structure.getAtomPointers());
		writer.close();
	} else {
		structure.writePdb(opt.output);
	}

	return 0;
}
