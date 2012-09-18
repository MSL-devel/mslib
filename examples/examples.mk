
EXAMPLES = example_multipleAtomsCoordinates example_multipleResidueIdentities example_SasaCalculator_usage example_read_write_PDBs_with_the_AtomContainer \
           example_read_write_PDBs_with_the_System example_looping_over_Chain_Residues_Atoms example_AtomPointerVector \
           example_selecting_atoms example_measurements \
	   example_mutation_rotamers example_add_atoms_to_System_and_AtomContainer example_multiple_coordinates_from_NMR_multiModel_PDB example_add_identity_to_position \
	   example_ccd example_backrub example_bbq example_coiled_coils_and_symmetric_bundles 

# EXAMPLES_THAT_DO_NOT_COMPLILE = 

ifeq ($(MSL_GSL),T)
	EXAMPLES +=  example_molecular_alignment  example_pdbfrag
endif


ifeq ($(MSL_BOOST),T)
	EXAMPLES +=  example_regular_expressions 
endif