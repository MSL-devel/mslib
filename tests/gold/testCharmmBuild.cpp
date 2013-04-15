/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries) 
 Copyright (C) 2008-2013 The MSL Developer Group (see README.TXT)
 MSL Libraries: http://msl-libraries.org

If used in a scientific publication, please cite: 
 Kulp DW, Subramaniam S, Donald JE, Hannigan BT, Mueller BK, Grigoryan G and 
 Senes A "Structural informatics, modeling and design with a open source 
 Molecular Software Library (MSL)" (2012) J. Comput. Chem, 33, 1645-61 
 DOI: 10.1002/jcc.22968

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/




/******************************************************************
 *  This test builds a system from polymer sequence, writes it to
 *  PDB, then it builds a second system using the sequence from
 *  the PBD and reading the coordinates from the PDB.
 *  The test checks both for correct coordinates and correct
 *  energies with the two systems
 ******************************************************************/

#include <iostream>

#include "System.h"
#include "CharmmSystemBuilder.h"
#include "PDBWriter.h"
#include "AtomSelection.h"
#include "Transforms.h"

using namespace std;

using namespace MSL;

#include "SysEnv.h"
static SysEnv SYSENV;

bool checkCoor(AtomPointerVector & _atoms, unsigned int _set, double _epsilon);

int main() {

	bool result = true;
	double epsilon = 1e-8;

	System sys;

	PolymerSequence seq("\
A: ALA ILE VAL ILE\n\
B: ARG HSD PHE GLY");

	cout << seq << endl;

	CharmmSystemBuilder CSB(sys, SYSENV.getEnv("MSL_CHARMM_TOP"),SYSENV.getEnv("MSL_CHARMM_PAR"));
	CSB.buildSystem(seq);
	sys.printIcTable();

	if (!sys.seed("A 1 C", "A 1 CA", "A 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 A" << endl;
		result = false;
	}
	if (!sys.seed("B 1 C", "B 1 CA", "B 1 N")) {
		cerr << "Cannot seed atoms C, CA, N on residue 1 B" << endl;
		result = false;
	}
	sys.buildAtoms();
	AtomSelection sel(sys.getAllAtomPointers());
	sel.select("chainB, chain B");

	Transforms tr;
	tr.translate(sel.getSelection("chainB"), CartesianPoint(10,5,5));
	string filename = "/tmp/buildFromCharmmTopology.pdb";

	AtomPointerVector atoms = sys.getAtomPointers();
	cout << atoms;

	double E1 = sys.calcEnergy();
	cout << setprecision(15) << E1 << endl;
	cout << sys.getEnergySummary();

//	for (unsigned int i=0; i<atoms.size(); i++) {
//		cout << i << " " << setprecision(15) << atoms[i]->getX() << " " << atoms[i]->getY() << " " << atoms[i]->getZ() << endl;
//	}

	if (!sys.writePdb(filename)) {
		cerr << "Cannot write output file " << filename << endl;
		result = false;
	}
	cout << "Written pdb file " << filename << endl;
	cout << endl;

	cout << "=====================================================" << endl;
	cout << "Create a new System from the PDB we previously saved with the buildSystemFromPDB function" << endl;
	
	System sys2;
	CharmmSystemBuilder CSB2(sys2, SYSENV.getEnv("MSL_CHARMM_TOP"),SYSENV.getEnv("MSL_CHARMM_PAR"));

	CSB2.buildSystemFromPDB("/tmp/buildFromCharmmTopology.pdb");

	AtomPointerVector atoms2 = sys2.getAtomPointers();
	//for (unsigned int i=0; i<atoms2.size(); i++) {
	//	cout << i << " " << setprecision(15) << atoms2[i]->getX() << " " << atoms2[i]->getY() << " " << atoms2[i]->getZ() << endl;
	//}

	AtomSelection sel2(atoms2);
	AtomPointerVector noCrd = sel2.select("noCrd, HASCOOR 0");

	if (noCrd.size() != 0) {
		cerr << "Error building from PDB, there are " << noCrd.size() << " atoms without coordinates out of " << sys2.atomSize() << endl;
		result = false;
	} else { 
		cout << "System build OK from PDB " << filename << endl;
	}

	// NOTE the energies will be slightly different because writing the PDB rounds the coordinates to 3 decimal digits
	cout << "Calculate the energies (a small difference will occur due to rounding to 3 digits when the PDB file was written)" << endl;

	double E2 = sys2.calcEnergy();
	cout << setprecision(15) << E2 << endl;
	cout << sys2.getEnergySummary();

	cout << "=====================================================" << endl;
	cout << "    RUN THE FINAL TESTS: " << endl;

	cout << endl;
	cout << " - check the atoms' coordinates of the system built using a polymer sequence:";
	if (checkCoor(atoms, 0, epsilon)) {
		cout << " OK" << endl;
	} else {
		cout << " NOT OK" << endl;
		result = false;
	}

	cout << endl;
	cout << " - check the energy of the system built using a polymer sequence:";
	if (fabs(E1 - 3695.01125126651) < epsilon) {
		cout << " OK" << endl;
	} else {
		cout << " NOT OK" << endl;
		result = false;
	}

	cout << endl;
	cout << " - check the atoms' coordinates of the system built from the PDB:";
	if (checkCoor(atoms2, 1, epsilon)) {
		cout << " OK" << endl;
	} else {
		cout << " NOT OK" << endl;
		result = false;
	}

	cout << endl;
	cout << " - check the energy of the system built from the PDB:";
	if (fabs(E2 - 3690.06144455409) < epsilon) {
		cout << " OK" << endl;
	} else {
		cout << " NOT OK" << endl;
		result = false;
	}
	cout << endl;
	cout << endl;
	cout << "    DONE WITH THE FINAL TESTS: " << endl;
	cout << "=====================================================" << endl;


	if (result) {
		cout << "GOLD" << endl;
	} else {
		cout << "LEAD" << endl;
	}
	
	return 0;
}

bool checkCoor(AtomPointerVector & _atoms, unsigned int _set, double _epsilon) {

	bool out = true;

	vector<CartesianPoint> expectedCoor1(138);
	expectedCoor1[0].setCoor(2.14272956301791, 1.32844843887078, 0);
	expectedCoor1[1].setCoor(1.539, 0, 0);
	expectedCoor1[2].setCoor(1.87414722514083, -0.48503732226756, -0.909058514886646);
	expectedCoor1[3].setCoor(3.17886640099377, 1.23889145706922, 2.40107869948561e-16);
	expectedCoor1[4].setCoor(1.84011162925603, 1.84730485850798, 0.849005537064952);
	expectedCoor1[5].setCoor(1.84011162925603, 1.84730485850798, -0.849005537064952);
	expectedCoor1[6].setCoor(2.06222622934674, -0.78348494098001, 1.21844232534049);
	expectedCoor1[7].setCoor(1.68131865491636, -1.82697261546193, 1.19936634709838);
	expectedCoor1[8].setCoor(3.17201314176048, -0.81458383055927, 1.20532904573132);
	expectedCoor1[9].setCoor(1.7306885350093, -0.300884031791201, 2.161458295345);
	expectedCoor1[10].setCoor(0, 0, 0);
	expectedCoor1[11].setCoor(-0.66107931076106, 1.03688776387981, 5.07911737539159e-16);
	expectedCoor1[12].setCoor(-0.612144477169287, -1.20974079003361, -2.96291251590652e-16);
	expectedCoor1[13].setCoor(-0.0524505108006792, -2.035784072199, -4.89342133394197e-16);
	expectedCoor1[14].setCoor(-2.05155298920284, -1.46173449131471, -7.1313130149769e-17);
	expectedCoor1[15].setCoor(-2.47501240194844, -1.09314124240874, 0.925892532974584);
	expectedCoor1[16].setCoor(-2.80025710052581, -0.813679468706866, -1.17446893588584);
	expectedCoor1[17].setCoor(-2.33582181331223, -1.15395985474854, -2.12709726470363);
	expectedCoor1[18].setCoor(-4.28486963237707, -1.23713223309194, -1.18998011265759);
	expectedCoor1[19].setCoor(-4.81871916700626, -0.759674397118511, -2.03656023245947);
	expectedCoor1[20].setCoor(-4.38907348541936, -2.33285154282886, -1.31552528891219);
	expectedCoor1[21].setCoor(-4.78017931733158, -0.936677659825328, -0.241973590848046);
	expectedCoor1[22].setCoor(-2.67167741750261, 0.721450905789588, -1.09036282188673);
	expectedCoor1[23].setCoor(-3.24016229550911, 1.07950093431709, -0.202634253826989);
	expectedCoor1[24].setCoor(-1.60771411829703, 0.998262393354835, -0.932914867643168);
	expectedCoor1[25].setCoor(-3.17141310761921, 1.45461395467687, -2.33091561148095);
	expectedCoor1[26].setCoor(-3.03889714733243, 2.55117127292878, -2.20402665192398);
	expectedCoor1[27].setCoor(-2.6012385226534, 1.13047183690783, -3.22617471930799);
	expectedCoor1[28].setCoor(-4.24905020321679, 1.25460831963044, -2.50127017156574);
	expectedCoor1[29].setCoor(-2.22140121348204, -2.97120876426844, -9.27808641787354e-16);
	expectedCoor1[30].setCoor(-1.23920874602462, -3.71161267894604, -2.14443262326882e-15);
	expectedCoor1[31].setCoor(-3.47377733711965, -3.46579575858029, -6.05106911099473e-16);
	expectedCoor1[32].setCoor(-4.26896413672448, -2.86505495446071, 1.56855806228405e-16);
	expectedCoor1[33].setCoor(-3.79084348619396, -4.87732328119938, -1.61927156526465e-15);
	expectedCoor1[34].setCoor(-3.46042507423338, -5.32887896513626, -0.925857644209815);
	expectedCoor1[35].setCoor(-3.15665663560327, -5.65618696551472, 1.15839571813461);
	expectedCoor1[36].setCoor(-3.4551395640119, -5.17080462334018, 2.11540270000766);
	expectedCoor1[37].setCoor(-3.64218766719152, -7.11944846490544, 1.19300513543242);
	expectedCoor1[38].setCoor(-3.14638929473014, -7.66656113616527, 2.02196454867989);
	expectedCoor1[39].setCoor(-4.73800931684174, -7.1686787513641, 1.35956305644807);
	expectedCoor1[40].setCoor(-3.40432400329822, -7.6315050800708, 0.237560244156815);
	expectedCoor1[41].setCoor(-1.6218699427734, -5.60931579744845, 1.0592463806863);
	expectedCoor1[42].setCoor(-1.17192207233184, -6.1833841015889, 1.89821562092483);
	expectedCoor1[43].setCoor(-1.28332953305476, -6.06030151608424, 0.103159421997868);
	expectedCoor1[44].setCoor(-1.25186233893336, -4.5654624563886, 1.11311673857218);
	expectedCoor1[45].setCoor(-5.30692618632189, -4.95359415619627, -4.25316226658603e-16);
	expectedCoor1[46].setCoor(-5.98707420463395, -3.92911342763854, 1.58703653524312e-15);
	expectedCoor1[47].setCoor(-5.87516960498206, -6.17497772704918, -1.34980252093925e-15);
	expectedCoor1[48].setCoor(-5.32706356061084, -7.00875503792046, -2.75713715338086e-15);
	expectedCoor1[49].setCoor(-7.30393609253405, -6.40459639799809, -1.12494125122506e-16);
	expectedCoor1[50].setCoor(-7.72220786274085, -6.03012666139716, 0.925892532974585);
	expectedCoor1[51].setCoor(-8.04351896208669, -5.74615100337265, -1.17446893588584);
	expectedCoor1[52].setCoor(-7.58388000227863, -6.09288274948251, -2.12709726470363);
	expectedCoor1[53].setCoor(-9.53389910372206, -6.14883406337654, -1.18998011265759);
	expectedCoor1[54].setCoor(-10.0610302484288, -5.66396906460222, -2.03656023245946);
	expectedCoor1[55].setCoor(-9.65339142999417, -7.24299165380001, -1.31552528891219);
	expectedCoor1[56].setCoor(-10.0249655069387, -5.84149317407206, -0.241973590848044);
	expectedCoor1[57].setCoor(-7.89351804536266, -4.21296551995262, -1.09036282188673);
	expectedCoor1[58].setCoor(-8.45694835067373, -3.84701310431957, -0.202634253826986);
	expectedCoor1[59].setCoor(-6.8257935651244, -3.9510362623191, -0.932914867643166);
	expectedCoor1[60].setCoor(-8.38296846854121, -3.47289653703255, -2.33091561148095);
	expectedCoor1[61].setCoor(-8.2351550943493, -2.37829631852693, -2.20402665192398);
	expectedCoor1[62].setCoor(-7.81737519261174, -3.80496793910762, -3.22617471930799);
	expectedCoor1[63].setCoor(-9.46329303514065, -3.65783651286267, -2.50127017156573);
	expectedCoor1[64].setCoor(-7.4948433127056, -7.91155208191397, -2.15473871842149e-15);
	expectedCoor1[65].setCoor(-6.46549433169933, -8.63821615512893, -4.92511228730838e-15);
	expectedCoor1[66].setCoor(-8.67287977740024, -8.35857566949046, -1.24743921629057e-15);
	expectedCoor1[67].setCoor(12.0167931637329, 6.36790032734568, 5);
	expectedCoor1[68].setCoor(11.5227, 5, 5);
	expectedCoor1[69].setCoor(11.8467047936869, 4.52718978841486, 4.08158062003857);
	expectedCoor1[70].setCoor(13.0567726350233, 6.3613658576215, 5);
	expectedCoor1[71].setCoor(11.6737103532535, 6.86093560314915, 5.84900553706495);
	expectedCoor1[72].setCoor(11.6737103532535, 6.86093560314915, 4.15099446293505);
	expectedCoor1[73].setCoor(12.0428937095205, 4.17699499446226, 6.19884637274138);
	expectedCoor1[74].setCoor(11.8308275635746, 4.73349873118282, 7.13709965984445);
	expectedCoor1[75].setCoor(11.4883482081183, 3.21084627540568, 6.24639571069609);
	expectedCoor1[76].setCoor(13.5424289015865, 3.87162542368778, 6.10643834819603);
	expectedCoor1[77].setCoor(13.7158575209919, 3.22285053362766, 5.21773881179544);
	expectedCoor1[78].setCoor(14.0945577808191, 4.82314063207849, 5.94697967355416);
	expectedCoor1[79].setCoor(14.0846939594489, 3.18218632189438, 7.35344576274488);
	expectedCoor1[80].setCoor(13.954663195506, 3.84026203169544, 8.24108708704104);
	expectedCoor1[81].setCoor(13.5432681371555, 2.22325187299011, 7.51696851852693);
	expectedCoor1[82].setCoor(15.53055184705, 2.90825134158581, 7.1587717906142);
	expectedCoor1[83].setCoor(15.9861398247332, 3.25186438414399, 6.34611031646001);
	expectedCoor1[84].setCoor(16.2578676996444, 2.20337110131501, 8.0400948059954);
	expectedCoor1[85].setCoor(15.6999907758145, 1.70845050378761, 9.15062091104704);
	expectedCoor1[86].setCoor(16.2495923130756, 1.18799892537464, 9.79259018972607);
	expectedCoor1[87].setCoor(14.7363166652198, 1.8726006065728, 9.31452451785658);
	expectedCoor1[88].setCoor(17.5560335427169, 1.99086433145408, 7.79986542260673);
	expectedCoor1[89].setCoor(18.0955325229571, 1.46886786639365, 8.44880288552917);
	expectedCoor1[90].setCoor(17.9694322685379, 2.35505841215018, 6.97396227751315);
	expectedCoor1[91].setCoor(10, 5, 5);
	expectedCoor1[92].setCoor(9.36066908076843, 6.04739218333655, 5);
	expectedCoor1[93].setCoor(9.38409346862468, 3.79744839835905, 5);
	expectedCoor1[94].setCoor(9.924565375249, 2.95751294955964, 5);
	expectedCoor1[95].setCoor(7.94022048135497, 3.60607690980119, 5);
	expectedCoor1[96].setCoor(7.52853211410902, 4.01607203191243, 5.9132395349082);
	expectedCoor1[97].setCoor(6.97052322304357, 6.58357155755545, 4.69899090850055);
	expectedCoor1[98].setCoor(6.47694637984147, 6.31863920776392, 5.52675476364801);
	expectedCoor1[99].setCoor(7.41717474252918, 5.72867836893231, 3.71260547866748);
	expectedCoor1[100].setCoor(7.23566890024863, 4.2442641927652, 3.78245240280982);
	expectedCoor1[101].setCoor(7.63021196284581, 3.79209720294777, 2.84889103628061);
	expectedCoor1[102].setCoor(6.14529070148807, 4.02927814705623, 3.83030784991963);
	expectedCoor1[103].setCoor(7.95294358746725, 7.83196398281636, 3.18487202274138);
	expectedCoor1[104].setCoor(8.01945146292611, 6.51031806548459, 2.7819481320347);
	expectedCoor1[105].setCoor(8.49824772881789, 6.22371783258494, 1.85319267860793);
	expectedCoor1[106].setCoor(7.31624456700255, 7.83379205555088, 4.33850630588447);
	expectedCoor1[107].setCoor(7.09308019730544, 8.72192700070761, 4.93209679880794);
	expectedCoor1[108].setCoor(7.67295643260847, 2.11321205046273, 5);
	expectedCoor1[109].setCoor(8.6077312132555, 1.31793763275127, 5);
	expectedCoor1[110].setCoor(6.38401426494826, 1.70878789388033, 5);
	expectedCoor1[111].setCoor(5.62766172105791, 2.36096311523621, 5);
	expectedCoor1[112].setCoor(5.97578453682297, 0.312751096702417, 5);
	expectedCoor1[113].setCoor(6.30946305402684, -0.168755141704573, 4.08900274449821);
	expectedCoor1[114].setCoor(6.45919661631349, -0.476252555242374, 6.23975602382906);
	expectedCoor1[115].setCoor(6.08077160971131, 0.00905990148437219, 7.16497587539083);
	expectedCoor1[116].setCoor(6.05758513080508, -1.51243797144926, 6.20618165938394);
	expectedCoor1[117].setCoor(7.95798450057635, -0.558935186206423, 6.32525871690785);
	expectedCoor1[118].setCoor(8.65768076125708, -1.57334317429885, 5.65625466872602);
	expectedCoor1[119].setCoor(8.1120875802738, -2.30465409530323, 5.07560445759692);
	expectedCoor1[120].setCoor(8.6819988796628, 0.377499473167671, 7.07919828185917);
	expectedCoor1[121].setCoor(8.15526599895976, 1.15748980307658, 7.60908842862063);
	expectedCoor1[122].setCoor(10.0546622125637, -1.64414347039235, 5.72774244934801);
	expectedCoor1[123].setCoor(10.5832488819972, -2.42826216608633, 5.20481763261006);
	expectedCoor1[124].setCoor(10.0788909499777, 0.315655377861483, 7.15287250915182);
	expectedCoor1[125].setCoor(10.6244810541682, 1.04595586498205, 7.73349396648941);
	expectedCoor1[126].setCoor(10.7667849148979, -0.696596721532331, 6.47668892643769);
	expectedCoor1[127].setCoor(11.8450530254107, -0.748302691819968, 6.53154496063379);
	expectedCoor1[128].setCoor(4.45288602130871, 0.310624727859039, 5);
	expectedCoor1[129].setCoor(3.82798067426889, 1.36854580111594, 5);
	expectedCoor1[130].setCoor(3.82884969432682, -0.884570466494985, 5);
	expectedCoor1[131].setCoor(4.36056161463987, -1.73055100634032, 5);
	expectedCoor1[132].setCoor(2.38597156839324, -1.04973100567552, 5);
	expectedCoor1[133].setCoor(1.98859394841277, -0.619572061140091, 4.09089039738059);
	expectedCoor1[134].setCoor(1.98393498961672, -0.619922659501641, 5.90758264906561);
	expectedCoor1[135].setCoor(2.06423611712287, -2.5118509424098, 5);
	expectedCoor1[136].setCoor(3.02363184099489, -3.32864947633116, 5);
	expectedCoor1[137].setCoor(0.850592162884773, -2.85047810071421, 5);

	vector<CartesianPoint> expectedCoor2(138);
	expectedCoor2[0].setCoor(2.143, 1.328, 0);
	expectedCoor2[1].setCoor(1.539, 0, 0);
	expectedCoor2[2].setCoor(1.874, -0.485, -0.909);
	expectedCoor2[3].setCoor(3.179, 1.239, 0);
	expectedCoor2[4].setCoor(1.84, 1.847, 0.849);
	expectedCoor2[5].setCoor(1.84, 1.847, -0.849);
	expectedCoor2[6].setCoor(2.062, -0.783, 1.218);
	expectedCoor2[7].setCoor(1.681, -1.827, 1.199);
	expectedCoor2[8].setCoor(3.172, -0.815, 1.205);
	expectedCoor2[9].setCoor(1.731, -0.301, 2.161);
	expectedCoor2[10].setCoor(0, 0, 0);
	expectedCoor2[11].setCoor(-0.661, 1.037, 0);
	expectedCoor2[12].setCoor(-0.612, -1.21, -0);
	expectedCoor2[13].setCoor(-0.052, -2.036, -0);
	expectedCoor2[14].setCoor(-2.052, -1.462, -0);
	expectedCoor2[15].setCoor(-2.475, -1.093, 0.926);
	expectedCoor2[16].setCoor(-2.8, -0.814, -1.174);
	expectedCoor2[17].setCoor(-2.336, -1.154, -2.127);
	expectedCoor2[18].setCoor(-4.285, -1.237, -1.19);
	expectedCoor2[19].setCoor(-4.819, -0.76, -2.037);
	expectedCoor2[20].setCoor(-4.389, -2.333, -1.316);
	expectedCoor2[21].setCoor(-4.78, -0.937, -0.242);
	expectedCoor2[22].setCoor(-2.672, 0.721, -1.09);
	expectedCoor2[23].setCoor(-3.24, 1.08, -0.203);
	expectedCoor2[24].setCoor(-1.608, 0.998, -0.933);
	expectedCoor2[25].setCoor(-3.171, 1.455, -2.331);
	expectedCoor2[26].setCoor(-3.039, 2.551, -2.204);
	expectedCoor2[27].setCoor(-2.601, 1.13, -3.226);
	expectedCoor2[28].setCoor(-4.249, 1.255, -2.501);
	expectedCoor2[29].setCoor(-2.221, -2.971, -0);
	expectedCoor2[30].setCoor(-1.239, -3.712, -0);
	expectedCoor2[31].setCoor(-3.474, -3.466, -0);
	expectedCoor2[32].setCoor(-4.269, -2.865, 0);
	expectedCoor2[33].setCoor(-3.791, -4.877, -0);
	expectedCoor2[34].setCoor(-3.46, -5.329, -0.926);
	expectedCoor2[35].setCoor(-3.157, -5.656, 1.158);
	expectedCoor2[36].setCoor(-3.455, -5.171, 2.115);
	expectedCoor2[37].setCoor(-3.642, -7.119, 1.193);
	expectedCoor2[38].setCoor(-3.146, -7.667, 2.022);
	expectedCoor2[39].setCoor(-4.738, -7.169, 1.36);
	expectedCoor2[40].setCoor(-3.404, -7.632, 0.238);
	expectedCoor2[41].setCoor(-1.622, -5.609, 1.059);
	expectedCoor2[42].setCoor(-1.172, -6.183, 1.898);
	expectedCoor2[43].setCoor(-1.283, -6.06, 0.103);
	expectedCoor2[44].setCoor(-1.252, -4.565, 1.113);
	expectedCoor2[45].setCoor(-5.307, -4.954, -0);
	expectedCoor2[46].setCoor(-5.987, -3.929, 0);
	expectedCoor2[47].setCoor(-5.875, -6.175, -0);
	expectedCoor2[48].setCoor(-5.327, -7.009, -0);
	expectedCoor2[49].setCoor(-7.304, -6.405, -0);
	expectedCoor2[50].setCoor(-7.722, -6.03, 0.926);
	expectedCoor2[51].setCoor(-8.044, -5.746, -1.174);
	expectedCoor2[52].setCoor(-7.584, -6.093, -2.127);
	expectedCoor2[53].setCoor(-9.534, -6.149, -1.19);
	expectedCoor2[54].setCoor(-10.061, -5.664, -2.037);
	expectedCoor2[55].setCoor(-9.653, -7.243, -1.316);
	expectedCoor2[56].setCoor(-10.025, -5.841, -0.242);
	expectedCoor2[57].setCoor(-7.894, -4.213, -1.09);
	expectedCoor2[58].setCoor(-8.457, -3.847, -0.203);
	expectedCoor2[59].setCoor(-6.826, -3.951, -0.933);
	expectedCoor2[60].setCoor(-8.383, -3.473, -2.331);
	expectedCoor2[61].setCoor(-8.235, -2.378, -2.204);
	expectedCoor2[62].setCoor(-7.817, -3.805, -3.226);
	expectedCoor2[63].setCoor(-9.463, -3.658, -2.501);
	expectedCoor2[64].setCoor(-7.495, -7.912, -0);
	expectedCoor2[65].setCoor(-6.465, -8.638, -0);
	expectedCoor2[66].setCoor(-8.673, -8.359, -0);
	expectedCoor2[67].setCoor(12.017, 6.368, 5);
	expectedCoor2[68].setCoor(11.523, 5, 5);
	expectedCoor2[69].setCoor(11.847, 4.527, 4.082);
	expectedCoor2[70].setCoor(13.057, 6.361, 5);
	expectedCoor2[71].setCoor(11.674, 6.861, 5.849);
	expectedCoor2[72].setCoor(11.674, 6.861, 4.151);
	expectedCoor2[73].setCoor(12.043, 4.177, 6.199);
	expectedCoor2[74].setCoor(11.831, 4.733, 7.137);
	expectedCoor2[75].setCoor(11.488, 3.211, 6.246);
	expectedCoor2[76].setCoor(13.542, 3.872, 6.106);
	expectedCoor2[77].setCoor(13.716, 3.223, 5.218);
	expectedCoor2[78].setCoor(14.095, 4.823, 5.947);
	expectedCoor2[79].setCoor(14.085, 3.182, 7.353);
	expectedCoor2[80].setCoor(13.955, 3.84, 8.241);
	expectedCoor2[81].setCoor(13.543, 2.223, 7.517);
	expectedCoor2[82].setCoor(15.531, 2.908, 7.159);
	expectedCoor2[83].setCoor(15.986, 3.252, 6.346);
	expectedCoor2[84].setCoor(16.258, 2.203, 8.04);
	expectedCoor2[85].setCoor(15.7, 1.708, 9.151);
	expectedCoor2[86].setCoor(16.25, 1.188, 9.793);
	expectedCoor2[87].setCoor(14.736, 1.873, 9.315);
	expectedCoor2[88].setCoor(17.556, 1.991, 7.8);
	expectedCoor2[89].setCoor(18.096, 1.469, 8.449);
	expectedCoor2[90].setCoor(17.969, 2.355, 6.974);
	expectedCoor2[91].setCoor(10, 5, 5);
	expectedCoor2[92].setCoor(9.361, 6.047, 5);
	expectedCoor2[93].setCoor(9.384, 3.797, 5);
	expectedCoor2[94].setCoor(9.925, 2.958, 5);
	expectedCoor2[95].setCoor(7.94, 3.606, 5);
	expectedCoor2[96].setCoor(7.529, 4.016, 5.913);
	expectedCoor2[97].setCoor(6.971, 6.584, 4.699);
	expectedCoor2[98].setCoor(6.477, 6.319, 5.527);
	expectedCoor2[99].setCoor(7.417, 5.729, 3.713);
	expectedCoor2[100].setCoor(7.236, 4.244, 3.782);
	expectedCoor2[101].setCoor(7.63, 3.792, 2.849);
	expectedCoor2[102].setCoor(6.145, 4.029, 3.83);
	expectedCoor2[103].setCoor(7.953, 7.832, 3.185);
	expectedCoor2[104].setCoor(8.019, 6.51, 2.782);
	expectedCoor2[105].setCoor(8.498, 6.224, 1.853);
	expectedCoor2[106].setCoor(7.316, 7.834, 4.339);
	expectedCoor2[107].setCoor(7.093, 8.722, 4.932);
	expectedCoor2[108].setCoor(7.673, 2.113, 5);
	expectedCoor2[109].setCoor(8.608, 1.318, 5);
	expectedCoor2[110].setCoor(6.384, 1.709, 5);
	expectedCoor2[111].setCoor(5.628, 2.361, 5);
	expectedCoor2[112].setCoor(5.976, 0.313, 5);
	expectedCoor2[113].setCoor(6.309, -0.169, 4.089);
	expectedCoor2[114].setCoor(6.459, -0.476, 6.24);
	expectedCoor2[115].setCoor(6.081, 0.009, 7.165);
	expectedCoor2[116].setCoor(6.058, -1.512, 6.206);
	expectedCoor2[117].setCoor(7.958, -0.559, 6.325);
	expectedCoor2[118].setCoor(8.658, -1.573, 5.656);
	expectedCoor2[119].setCoor(8.112, -2.305, 5.076);
	expectedCoor2[120].setCoor(8.682, 0.377, 7.079);
	expectedCoor2[121].setCoor(8.155, 1.157, 7.609);
	expectedCoor2[122].setCoor(10.055, -1.644, 5.728);
	expectedCoor2[123].setCoor(10.583, -2.428, 5.205);
	expectedCoor2[124].setCoor(10.079, 0.316, 7.153);
	expectedCoor2[125].setCoor(10.624, 1.046, 7.733);
	expectedCoor2[126].setCoor(10.767, -0.697, 6.477);
	expectedCoor2[127].setCoor(11.845, -0.748, 6.532);
	expectedCoor2[128].setCoor(4.453, 0.311, 5);
	expectedCoor2[129].setCoor(3.828, 1.369, 5);
	expectedCoor2[130].setCoor(3.829, -0.885, 5);
	expectedCoor2[131].setCoor(4.361, -1.731, 5);
	expectedCoor2[132].setCoor(2.386, -1.05, 5);
	expectedCoor2[133].setCoor(1.989, -0.62, 4.091);
	expectedCoor2[134].setCoor(1.984, -0.62, 5.908);
	expectedCoor2[135].setCoor(2.064, -2.512, 5);
	expectedCoor2[136].setCoor(3.024, -3.329, 5);
	expectedCoor2[137].setCoor(0.851, -2.85, 5);

	if (_set == 0) {
		if (_atoms.size() != expectedCoor1.size()) {
			return false;
		}
		for (unsigned int i=0; i<_atoms.size(); i++) {
			if (_atoms[i]->getCoor().distance(expectedCoor1[i]) > _epsilon) {
				out = false;
			}
		}
	} else {
		if (_atoms.size() != expectedCoor2.size()) {
			return false;
		}
		for (unsigned int i=0; i<_atoms.size(); i++) {
			if (_atoms[i]->getCoor().distance(expectedCoor2[i]) > _epsilon) {
				out = false;
			}
		}
	}
	return out;
}


