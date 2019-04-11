Author: Eric Witt
Version: v01

##########################################################################################################################################################

    SBMLComp is a comparison tool for SBML files at reaction level. This can be used for comparative analyses
    of different organisms in relation to their metabolic similarities and potentials. 
    For the comparison first the Species of both SBML files are matched, for this the database IDs, names and Inchis of the Species existing in the SBML 
    file are used. Since often different naming conventions are used for the species, the existing database IDs are used to extract synonyms from the 
    databases to increase the matching success.
    After matching the species, the reactions are matched. Therefore the reactions where the reactant and product numbers are identical are compared. 
    The matched species are now used to compare the reactions and determine the identity of the reactions. 
    Furthermore, the reactions are classified in a match. Therefore the reversibility and stoichiometry is checked. A perfect match would mean that the 
    reactants, products, reversibility and stoichiometry of the compared reactions match. If the reactions do not match, it will also be checked whether 
    they are simply upside down. This means that the reactants of one reaction are the products of the other reaction and the products are the reactants. 
    After Matching the reactions of the SBML files, the Jaccard similarity coefficient is calculated based on the matched reactions as a measure for the 
    similarity of the SBML files.

##########################################################################################################################################################

"SBMLComp" is written in Python version 2.7.14

It works under windows.

INPUT:
    - 2 SBML files (.xml) which you want compare

OUTPUT:
    - 01_matching_report.txt
    - 02_species_matches.txt
    - 03_species_id_assignment.txt
    - 04_reaction_matches.txt
    - 05_wrong_reactants_and_products.txt
    - 06_wrong_reactants.txt
    - 07_wrong_products.txt
    - 08_reactants_products_changed.txt
    - 09_reaction_id_assignment.txt
    - 10_SBML1_species_no_match.txt
    - 11_SBML2_species_no_match.txt
    - 12_SBML1_reactions_no_match.txt
    - 13_SBML2_reactions_no_match.txt

##########################################################################################################################################################

The scripts uses the biopython package from Cock et al. (2009) for accessing the KEGG database (Kanehisa and Goto, 2000), libSBML Python libary from 
Bornstein et al. (2008) for parsing the SBML files, libChEBIpy form Swainston et al. (2016) for accessing the ChEBI database (Hastings et al., 2016), 
PubChempy from Swain (2014) for accessing PubChem database (Kim et al., 2016).

Bornstein, B. J., Keating, S. M., Jouraku, A., and Hucka, M. (2008) LibSBML: An API Library for SBML. Bioinformatics, 24(6), 880-881

Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) 
Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25, 1422-1423

Hastings, J., Owen, G., Dekker, A., Ennis, M., Kale, N., Muthukrishnan, V., Turner, S., Swainston, N., Mendes, P., Steinbeck, C. (2016) ChEBI in 2016: 
Improved services and an expanding collection of metabolites. Nucleic acids research 44 (D1), D1214-9.

Kanehisa, M. and Goto, S. (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 29-34 

Kim, S., Thiessen, P.A., Bolton, E.E., Chen, J., Fu, G., Gindulyte, A., Han, L., He, J., He, S., Shoemaker, B.A., Wang, J., Yu, B., Zhang, J., Bryant, 
S.H. (2016) PubChem Substance and Compound databases. Nucleic acids research 44 (D1), D1202-13

Swain, M., 2014, PubChemPy. https://pubchempy.readthedocs.io/en/latest/index.html. Accessed 25 October 2018

Swainston, N., Hastings, J., Dekker, A., Muthukrishnan, V., May, J., Steinbeck, C., Mendes, P. (2016) libChEBI: an API for accessing the ChEBI database. 
Journal of Cheminformatics, 8 


##########################################################################################################################################################
###Executables############################################################################################################################################

This folder contains the applications of the GUI version of "SBMLComp" for Windows. 
These applications are based on the "GUI_sbmlcomp_v01.py" script which is stored in the scripts directory.

##########################################################################################################################################################
###Scripts################################################################################################################################################

"CLI_sbmlcomp_v01.py"
    This is the command line version of the script that could be included in pipelines. To run the script furthe spectifications are necessary otherwise 
    the script doesn't work.
		-SBML1-file-path		Enter the path of SBML file 1 (e.g.: C:\Input\Sfu2014_1_PWYTools.xml)
		-SBML2-file-path		Enter the path of SBML file 2 (e.g.: C:\Input\Sfu2017_1_PWYTools.xml).
						For running the Synonym-matching functions, this file should contain Kegg-, ChEBI, or PubChem-ids.
		-Report_storage_path:		Enter the path of the location for storing the output files (e.g.: C:\Output).
		-Species-id-Kegg-id-check	y or n
						y means the species-ids will be checked for Kegg-ids, n means the opposite	

	You can run it by using python 2.7: 
		python CLI_sbmlcomp_v01.py SBML1-file-path SBML2-file-path Report_storage_path Species-id-Kegg-id-check

	e.g.:	python CLI_sbmlcomp_v01.py C:\Input\Sfu2014_1_PWYTools.xml C:\Input\Sfu2017_1_PWYTools.xml C:\Output y

"GUI_sbmlcomp_v01.py"
	This is the graphical user interface guided version of the script "SBMLComp" controlled by the user, but takes much more processing time than the 
	CLI version of SBMLComp. You can use the Executable (Executables directory) if you do not know how to deal with the Python file.

	You can run it by using python 2.7: 
		python GUI_sbmcomp_v01.py

##########################################################################################################################################################
###Sample_Input###########################################################################################################################################

Sfu2014_1_PWYTools.xml
	This SBML file is exported from the draft reconstruction generated by Pathway Tools (Pathway-Prediction-Score-Cutoff = 0.99) based on the Syntrophomonas 
    	fumaroxidans annotation from 2014. The used annotation can be found under the following link: https://www.ncbi.nlm.nih.gov/nuccore/CP000478

Sfu2017_1_PWYTools.xml	
    	This SBML file is exported from the draft reconstruction generated by Pathway Tools (Pathway-Prediction-Score-Cutoff = 0.99) based on the Syntrophomonas 
    	fumaroxidans annotation from 2017. The used annotation can be found under the following link: https://www.ncbi.nlm.nih.gov/nuccore/NC_008554.1

##########################################################################################################################################################
###Sample_Output##########################################################################################################################################

01_matching_report.txt
    This file is the main report. All process steps, the matches made during the step, the process times, and the summary statistics are listed.
    At the end of the comparison the similarity is determined by the jaccard simmilarity coeficient.
    
02_species_matches.txt
    This file lists all species matches that occurred during the matching process. Each line in this file, except the first, is a match (tab-separated). 
    The specifications in the line are explained below.
    
    The first line of the file contains the following information:
    
    SBML1-file-name: Number-Of-Species	SBML2-file-name: Number-Of-Species
    
    The following lines of the file contains the following information:
    
    SBML1-Species-ID    SBML1-Species-Name  SBML1-Species-Formula  ||    SBML2-Species-ID    SBML2-Species-Name  SBML2-Species-Formula  || Match-type   (SBML1-ID:   Database-ID SBML2-ID:   Database-ID)*
    
    ()* - This part is only given for ID- or INCHI-Match-types. The ID of the SBML-Species in the respective database is given as "Database-ID"
    
03_species_id_assignment.txt
    This file lists for each SBML1-Species-ID the matched SBML2-Species-ID. The specifications in the line are explained below.
    
    - Single-Match:
        SBML1-Species-ID    ['SBML2-Species-ID']
    
    - Multiple-Matches:
        SBML1-Species-ID    ['SBML2-Species-ID', 'SBML2-Species-ID', ...]

04_reaction_matches.txt
    This file lists all reaction matches that occurred during the matching process. The file has the following structure:
    
    First line:
        SBML1-file-name: Number-Of-Reactions	SBML2-file-name: Number-Of-Reactions
        
    For each reaction-match:
    
    Reaction-Qualification: Number  (Number can be 123, 12, 13, 1, r123, r12, r13, r1)
        SBML1-Reaction-ID   SBML1-Reaction-Name (reversible: SBML1-Reaction-Reversebility) :
            Reactants (Stoichiometry, Species-ID, (Species-Name))   |   Productss (Stoichiometry, Species-ID, (Species-Name))
    SBML1-Reaction-Notes
        SBML2-Reaction-ID   SBML2-Reaction-Name (reversible: SBML2-Reaction-Reversebility) :
            Reactants (Stoichiometry, Species-ID, (Species-Name))   |   Productss (Stoichiometry, Species-ID, (Species-Name))
    SBML2-Reaction-Notes
        
05_wrong_reactants_and_products.txt
    This file lists the reactions of the matching that have differences in reactants and products.

06_wrong_reactants.txt
    This file lists the reactions of the matching that have the same products but differences in the reactants.

07_wrong_products.txt
    This file lists the reactions of the matching that have the same reactants but differences in the products.

08_reactants_products_changed.txt
    This file lists the reactions of the matching where the reactants and products of SBML1-Reaction are the products and reactants of SBML2-Reaction.

09_reaction_id_assignment.txt
    This file lists for each SBML1-Reaction-ID the matched SBML2-Reaction-ID. The specifications in the line are explained below.
    
    - Single-Match:
        SBML1-Reaction-ID    ['SBML2-Reaction-ID']
    
    - Multiple-Matches:
        SBML1-Reaction-ID    ['SBML2-Reaction-ID', 'SBML2-Reaction-ID', ...]

10_SBML1_species_no_match.txt
    This file lists the unmatched species ids of the SBML1 file.

11_SBML2_species_no_match.txt
    This file lists the unmatched species ids of the SBML2 file.

12_SBML1_reactions_no_match.txt
    This file lists the unmatched reaction ids of the SBML1 file.

13_SBML2_reactions_no_match.txt
    This file lists the unmatched reaction ids of the SBML2 file.
