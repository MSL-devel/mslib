from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain
from glob import glob

# Set background color
cmd.bg_color('grey50')

# Set Helices to be fancy
cmd.do("set cartoon_fancy_helices,1")
cmd.do("set cartoon_highlight_color,deepteal")

# Remove shadows when ray tracing
cmd.do("set ray_shadows,0")

cmd.do("set label_color,black")
cmd.do("set label_size,20")

# Add a global axis object
obj = [
        CYLINDER, 0., 0., -25.0, 0, 0., 25, 0.2, 1.0, 0.0, 0.0, 1.0, 0.75, 0.75,
        #CYLINDER, -26.135, .595, 2.623, -26.223, .250, 12.617, 0.2, 1.0, 0.0, 0.0, 1.0, 0.75, 0.75,
        #CYLINDER, -23.936, 9.207, -247.173, -26.223, .250, 12.617, 0.2, 1.0, 0.0, 0.0, 1.0, 0.75, 0.75,
        CYLINDER, 0,-10,0, 0, 10, 0, 0.2, 0.0, 0.0, 1.0, 0.75, 0.75, 1.0,
        #CYLINDER, -25.643,-6.979,2.776, -25.673, -7.618, 12.756, 0.2, 0.0, 0.0, 1.0, 0.75, 0.75, 1.0,
        #CYLINDER, -24.890, 9.020, -247.188, -25.673, -7.618, 12.756, 0.2, 0.0, 0.0, 1.0, 0.75, 0.75, 1.0,
        CYLINDER, -10,0,0,10,0,0, 0.2, 0.0, 1.0, 0.0, 0.75, 1.0, 0.75,
]

cmd.load_cgo(obj,"Axis")


# A function to "load *.pdb" from inside pymol
def loadAll(str):
    file_list = glob(str)
    for file in file_list:
        cmd.load(file)


def createSelectionDictionary(selection):
    # Object that will be returned.
    selectionDictionary = {}
    
    # Get atoms from this selection.
    atoms = cmd.get_model("byres (" + selection + " )")
    
    # Loop through atoms in the current selection
    for currAtom in atoms.atom:
        resNum = currAtom.resi
        chain = currAtom.chain
        key = (chain, resNum)

        # If the key is not in the dictionary yet, add it.
        if(key not in selectionDictionary):
            selectionDictionary[key] = currAtom.resn

    # Return the dictionary.
    return selectionDictionary

def diffSelections(selection1, selection2, selection3):
    chainsAndResidues = ''
    conjunction = ''
    selection1Dictionary = createSelectionDictionary(selection1)
    selection2Dictionary = createSelectionDictionary(selection2)

    for key in selection1Dictionary:
        if key in selection2Dictionary:
            if selection1Dictionary[key] != selection2Dictionary[key]:
                chainsAndResidues +=  conjunction + ' ( chain ' + key[0] + ' and resi ' + key[1] + ' )'
                conjunction = ' or'
                
    if chainsAndResidues != '':
        cmd.do( 'select %s, (%s or %s) and (%s)' % (selection3, selection1, selection2, chainsAndResidues) )
    
cmd.extend("diffSelections", diffSelections)

#####################################################
# MSL Stuff
import PythonMSL

def localSamplingCCD(system="all",fragSel="sele and name CA", bbqTable="/export/home/brettth/projectsS/code/tree/mslib/trunk/tables/PiscesBBQTable.txt",angle=10,numResults=25):
    cmd.read_pdbstr(PythonMSL.localSamplingCCD(cmd.get_pdbstr(system),cmd.get_pdbstr(fragSel),numResults,angle,bbqTable),"fragmentsCCD")

def localSamplingBR(system="all",fragSel="sele and name CA",bbqTable="/home/dwkulp/software/msl/tables/PiscesBBQTable.txt",numResults=25):
    cmd.read_pdbstr(PythonMSL.localSamplingBR(cmd.get_pdbstr(system),cmd.get_pdbstr(fragSel),numResults),"fragmentsBR")

def localSamplingPDB(system="all",fragSel="sele and name CA",fragDB="/var/lib/python-support/python2.6/pmg_tk/startup/lbbs/nr1000.fragdb",numRes=-1, rmsdTol=0.0,bbqTable="/var/lib/python-support/python2.6/pmg_tk/startup/lbbs/PiscesBBQTable.txt",numResults=25):
    cmd.read_pdbstr(PythonMSL.localSamplingPDB(cmd.get_pdbstr(system),cmd.get_pdbstr(fragSel),fragDB,numRes,rmsdTol,numResults,bbqTable),"fragmentsPDB")

def quench(sel="sele",rounds=10):
    cmd.read_pdbstr(PythonMSL.quickQuench(cmd.get_pdbstr(sel),rounds),"quench")

def getSasa(sel="sele",normalized=0,probeSize=1.4):
    name="SASA %3.1f" % probeSize
    if (normalized):
        name="normSASA %3.1f" % probeSize

    cmd.read_pdbstr(PythonMSL.getSasa(cmd.get_pdbstr(sel),normalized,probeSize),name)

cmd.extend("localSamplingBR",localSamplingBR)
cmd.extend("localSamplingCCD",localSamplingCCD)
cmd.extend("localSamplingPDB",localSamplingPDB)
cmd.extend("quench",quench)
cmd.extend("getSasa",getSasa)
#####################################################

print ".PyMOLRC.py loaded"
