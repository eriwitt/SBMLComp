#!/usr/bin/python
""" SBMLComp is a comparison tool for SBML files at reaction level. This can be used for comparative analyses
    of different organisms in relation to their metabolic similarities and potentials. 
    For the comparison first the Species of both SBML files are matched, for this the database IDs, names and Inchis of the Species existing in the SBML file are used. Since often different naming conventions are used for the species, the existing database IDs are used to extract synonyms from the databases to increase the matching success.
    After matching the species, the reactions are matched. Therefore the reactions where the reactant and product numbers are identical are compared. The matched species are now used to compare the reactions and determine the identity of the reactions. 
    Furthermore, the reactions are classified in a match. Therefore the reversibility and stoichiometry is checked. A perfect match would mean that the reactants, products, reversibility and stoichiometry of the compared reactions match. If the reactions do not match, it will also be checked whether they are simply upside down. This means that the reactants of one reaction are the products of the other reaction and the products are the reactants. 
    After Matching the reactions of the SBML files, the Jaccard similarity coefficient is calculated as a measure for the similarity of the SBML files.

    Author: Eric Witt
"""
import sys
sys.path.append("C:\Python27\Lib\site-packages\libsbml")
sys.path.append("C:\Python27\Lib\site-packages\Bio\KEGG")
sys.path.append("C:\Python27\Lib\site-packages\libchebipy")
sys.path.append("C:\Python27\Lib\site-packages")
from Bio.KEGG import Compound
from Bio.KEGG import REST
from libsbml import *
from tkinter import *
from tkinter import messagebox
import datetime
import libchebipy as lbc
import os
import platform
import pubchempy as pcp
import re
import socket
import tkFileDialog
import ttk
import urllib2

operating_system = platform.system()
if operating_system == "Windows":
    from os import startfile
else:
    import subprocess

time1 = datetime.datetime.now()
print time1
#####################################################################################################
# SBMLComp Functions

def species_id_kegg_check(species_id):
    #Check KEGG-Identifier in Species-ID
    d = ""
    res = re.search("[^a-z]c[0-9][0-9][0-9][0-9][0-9]", species_id)
    if res != None:
        d = check_compoundId(res.group()[1:])
    
    return d
        
    
def check_compoundId(ident):    
    #Function check the KEGG-ID, if it is existing
    aaaa = ""
    try:
        #print "Check KEGG-ID: "+id
        r = REST.kegg_get("compound:" + ident)
        aaaa = ident
        return aaaa
    except urllib2.HTTPError:
        pass 
    except urllib2.URLError:
        check_compoundId(ident)
    except socket.error:
        aaaa = check_compoundId(ident)
    return aaaa


def check_possible_identifier_to_comp_from_notes_sbml(sbml_Model, report, window1, progbar, progress_label):
    # Function checks the Species note section for ids from databases (KEGG, ChEBI, PubChem, BioCyc, HMDB, CAS, Chemspider & Metabolights) and Inchis. An Array with the numbers of the identifiers is returned.
    timea = datetime.datetime.now()
    output = []
    kegg = {}
    chebi = {}
    pubchem = {}
    inchi = {}
    biocyc = {}
    hmdb = {}
    cas = {}
    chemspider = {}
    metabolights = {}
    
    current_i1_nr = 0
    for i in sbml_Model.getListOfSpecies():
        current_i1_nr = current_i1_nr+1
        progress_nr = (100 / float(len(sbml_Model.getListOfSpecies()))) * current_i1_nr
        progress_label.config(text = "Species " + i.getId() + " is checked for identifiers!")
        progbar["value"] = progress_nr
        window1.update()
        if "<p>KEGG: " in i.getNotesString():
            kegg.update({i:i.getNotesString().split("<p>KEGG: ")[1].split("</p>")[0]})
        if "<p>CHEBI: " in i.getNotesString():
            chebi.update({i:i.getNotesString().split("<p>CHEBI: ")[1].split("</p>")[0]})
        if "<p>PUBCHEM: " in i.getNotesString():
            pubchem.update({i:i.getNotesString().split("<p>PUBCHEM: ")[1].split("</p>")[0]})
        if "<p>INCHI: " in i.getNotesString():
            inchi.update({i:i.getNotesString().split("<p>INCHI: ")[1].split("</p>")[0]})
        if "<p>BIOCYC: " in i.getNotesString():
            biocyc.update({i:i.getNotesString().split("<p>BIOCYC: ")[1].split("</p>")[0]})
        if "<p>HMDB: " in i.getNotesString():
            hmdb.update({i:i.getNotesString().split("<p>HMDB: ")[1].split("</p>")[0]})
        if "<p>CAS: " in i.getNotesString():
            cas.update({i:i.getNotesString().split("<p>CAS: ")[1].split("</p>")[0]})
        if "<p>CHEMSPIDER: " in i.getNotesString():
            chemspider.update({i:i.getNotesString().split("<p>CHEMSPIDER: ")[1].split("</p>")[0]})
        if "<p>METABOLIGHTS: " in i.getNotesString():
            metabolights.update({i:i.getNotesString().split("<p>METABOLIGHTS: ")[1].split("</p>")[0]})
        
    output.append(kegg)  #[0]
    output.append(chebi) #[1]
    output.append(pubchem) #[2]
    output.append(inchi) #[3]
    output.append(biocyc) #[4]
    output.append(hmdb) #[5]
    output.append(cas) #[6]
    output.append(chemspider) #[7]
    output.append(metabolights) #[8]
    timeaa = datetime.datetime.now() - timea
    report.write("\t\t- process-time: " + str(timeaa) + "\n")
    window1.destroy()
    return output


def db_id_Matching(identifier_sbml1, identifier_sbml2, i_s_match, j_s_match,species_ids, identifier_index, report, window1, progbar, progress_label, species_matches):
    #Function compares the entries of the database identifier lists and searches for identical IDs. If a hit was found, the Species-ID of the SBML file 2 is assigned to the Species-ID of the SBML file 1 (if not yet available in the list species_ids). Inputs are the lists with the respective database identifiers of both SBML files, list of sbml1-species-matches, list of sbml2-species-matches, species assignment list, and the identifier index which is needed for writing the report entries. The Function returns list of sbml1-species-matches, list of sbml2-species-matches, and species assignment list.
    timea = datetime.datetime.now()
    a = []
    if identifier_index == 0:
        print "KEGG-ID-Matching: " + "\n"
        identifier_index = 7
    if identifier_index == 1:
        print "ChEBI-ID-Matching: " + "\n"
        identifier_index = 8
    if identifier_index == 2:
        print "PubChem-ID-Matching: " + "\n"
        identifier_index = 9
    if identifier_index == 3:
        print "Inchi-Matching: " + "\n"
        identifier_index = 10
    if identifier_index == 4:
        print "BioCyc-ID-Matching: " + "\n"
        identifier_index = 9
    if identifier_index == 11:
        print "HMDB-ID-Matching: " + "\n"
        identifier_index = 9
    if identifier_index == 12:
        print "CAS-ID-Matching: " + "\n"
        identifier_index = 9
    if identifier_index == 13:
        print "Chemspider-ID-Matching: " + "\n"
        identifier_index = 9
    if identifier_index == 14:
        print "Metabolights-ID-Matching: " + "\n"
        identifier_index = 15
    
    current_i1_nr = 0
    for i in identifier_sbml1:
        current_i1_nr = current_i1_nr + 1
        progress_nr = (100 / float(len(identifier_sbml1))) * current_i1_nr
        #if i.getId() not in i_s_match:
        for j in identifier_sbml2:
            progress_label.config(text = "SBML1-Species-ID\n\n" + i.getId() + "\n\nis tried to match with SBML2-Species-ID\n\n" + j.getId())
            progbar["value"] = progress_nr
            window1.update()
            if identifier_sbml1[i].lower() == identifier_sbml2[j].lower():
                if j.getId() not in species_ids[i.getId()]:
                    species_ids[i.getId()].append(j.getId())
                e = "SBML1-ID:\t" +str(identifier_sbml1[i]) + "\tSBML2-ID:\t" + str(identifier_sbml2[j])
                write_species(i,j,e,identifier_index, species_matches)
                if i.getId() not in i_s_match:
                    i_s_match.append(i.getId())
                if j.getId() not in j_s_match:
                    j_s_match.append(j.getId())
                
    #window1.quit()
    window1.destroy()
    a.append(i_s_match)
    a.append(j_s_match)
    a.append(species_ids)
    timeaa = datetime.datetime.now() - timea
    report.write("\t\t- process-time: " + str(timeaa) + "\n")
    return a


def Name_Matching(sbml1_Model, sbml2_Model, i_s_match, j_s_match,species_ids, report, window1, progbar, progress_label, species_matches):
    #The function compares the species names, in a match the assignment of the SBML file 2 species ID to the SBML file 1 species ID (species_ids) takes place. If the names do not match, they are adapted using the adapt_name function. Inputs are the SBML models, list of sbml1-species-matches, list of sbml2-species-matches, and the species assignment list. The Function returns list of sbml1-species-matches, list of sbml2-species-matches, and species assignment list.
    timea = datetime.datetime.now()
    l = []
    print "Name-Matching: " + "\n"
    current_i1_nr = 0
    for i in sbml1_Model.getListOfSpecies():
        current_i1_nr = current_i1_nr + 1
        progress_nr = (100 / float(len(sbml1_Model.getListOfSpecies()))) * current_i1_nr
        for j in sbml2_Model.getListOfSpecies():
            progress_label.config(text = "SBML1-Species-Name\n\n" + i.getName() + "\n\nis tried to match with SBML2-Species-Name\n\n" + j.getName())
            progbar["value"] = progress_nr
            window1.update()
            if i.getName() == j.getName():
                if j.getId() not in species_ids[i.getId()]:
                    species_ids[i.getId()].append(j.getId())
                write_species(i,j,"",3,species_matches)
                if i.getId() not in i_s_match:
                    i_s_match.append(i.getId())
                if j.getId() not in j_s_match:
                    j_s_match.append(j.getId())
            else:
                cm3 = adapt_name(i,j,0)
                if cm3 == 1:
                    if j.getId() not in species_ids[i.getId()]:
                        species_ids[i.getId()].append(j.getId())
                    write_species(i,j,"",3,species_matches)
                    if i.getId() not in i_s_match:
                        i_s_match.append(i.getId())
                    if j.getId() not in j_s_match:
                        j_s_match.append(j.getId())
    #window1.quit()
    window1.destroy()
    l.append(i_s_match)
    l.append(j_s_match)
    l.append(species_ids)
    timeaa = datetime.datetime.now()-timea
    report.write("\t\t- process-time: " +str(timeaa)+"\n")
    return l


def write_species(i,j,e, match_index, species_matches):
    #This function generates the entry for the matched species. The entry is written in the species report file (02_species.txt). The entry contains the following information. 
    #SBML1-Species-ID    SBML1-Species-Name    SBML1-Species-Formula    ||    SBML2-Species-ID    SBML2-Species-Name    SBML2-Species-Formula    ||    MATCH-TYPE    SBML1-ID:    identifier    SBML2-ID:    identifier
    a = ""
    if match_index == 0:
        a = "KEGG-SYN-NAME-MATCH"
    if match_index == 1:
        a = "ChEBI-SYN-NAME-MATCH"
    if match_index == 2:
        a = "PubChem-SYN-NAME-MATCH"
    if match_index == 3:
        a = "NAME-MATCH"
    if match_index == 4:
        a = "KEGG-SYN2-NAME-MATCH"
    if match_index == 5:
        a = "ChEBI-SYN2-NAME-MATCH"
    if match_index == 6:
        a = "PubChem-SYN2-NAME-MATCH"
    if match_index == 7:
        a = "KEGG-ID-MATCH"
    if match_index == 8:
        a = "ChEBI-ID-MATCH"
    if match_index == 9:
        a = "PubChem-ID-MATCH"
    if match_index == 10:
        a = "Inchi-MATCH"
    if match_index == 11:
        a = "BioCyc-ID-MATCH"
    if match_index == 12:
        a = "HMDB-ID-MATCH"
    if match_index == 13:
        a = "CAS-ID-MATCH"
    if match_index == 14:
        a = "Chemspider-ID-MATCH"
    if match_index == 15:
        a = "Metabolights-ID-MATCH"
    if e == "":
        if "<p>FORMULA: " in j.getNotesString() and  "<p>FORMULA: " in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+ i.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+ j.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+a+"\n")
        if "<p>FORMULA: " not in j.getNotesString()and  "<p>FORMULA: " in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+ i.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+"-"+"\t"+"||"+"\t"+a+"\n")
        if "<p>FORMULA: "  in j.getNotesString()and  "<p>FORMULA: " not in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+"-"+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+ j.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+a+"\n")
        if "<p>FORMULA: " not in j.getNotesString()and  "<p>FORMULA: " not in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+ "-"+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+ "-"+"\t"+"||"+"\t"+a+"\n")
    else:
        if "<p>FORMULA: " in j.getNotesString() and  "<p>FORMULA: " in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+ i.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+ j.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+a+"\t" +e+"\n")
        if "<p>FORMULA: " not in j.getNotesString()and  "<p>FORMULA: " in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+ i.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+"-"+"\t"+"||"+"\t"+a+"\t" +e+"\n")
        if "<p>FORMULA: "  in j.getNotesString()and  "<p>FORMULA: " not in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+"-"+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+ j.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]+"\t"+"||"+"\t"+a+"\t" +e+"\n")
        if "<p>FORMULA: " not in j.getNotesString()and  "<p>FORMULA: " not in i.getNotesString():
            species_matches.write(i.getId()+"\t"+i.getName() + "\t"+ "-"+"\t"+"||"+"\t"+j.getId()+"\t"+j.getName() + "\t"+ "-"+"\t"+"||"+"\t"+a+"\t" +e+"\n")


def adapt_name(i,j,index):
    #This function adapts the species names to ensure matches for different name conventions. The function returns an integer which is set to 0 at the beginning. in the course of the adaption process it is tried again and again to match, this happens with the check_match function. When a match occurs, the integer is set to 1.
    cm2 = 0
    if index == 0:
        ii = i.getName()
        jj = j.getName()
# Name = fructofuranose --> fructose
        if "fructofuranose" in ii:
            ii = ii.replace("fructofuranose","fructose")
        if "fructofuranose" in jj:
            jj = jj.replace("fructofuranose","fructose")
# Name = Glucopyranose --> Glucose
        if "glucopyranose" in ii:
            ii = ii.replace("glucopyranose","glucose")
        if "glucopyranose" in jj:
            jj = jj.replace("glucopyranose","glucose")
# Name== Orthophospahte changed to phosphate
        if "orthophosphate" in ii:
            ii = ii.replace("orthophosphate","phosphate")
        if "orthophosphate" in jj:
            jj = jj.replace("orthophosphate","phosphate")
# Name== NH3 but FORMULA == H4N and CHARGE == 1 --> ammonium is meant
        if ii == "NH3" and i.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]=="H4N" and i.getNotesString().split("<p>CHARGE: ")[1].split("</p>")[0]=="1":
            ii="ammonium"
        if jj == "NH3" and j.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0] == "H4N" and j.getNotesString().split("<p>CHARGE: ")[1].split("</p>")[0] == "1":
            jj = "ammonium"
        cm2 = check_match(ii,jj)
    if index  == 1:
        ii = i.getName()
        jj = j
# some species name starts with "a ", it will be removed
    if cm2 == 0:
        if ii.startswith("a "):
            ii = ii[2:]
        if jj.startswith("a "):
            jj = jj[2:]
        cm2 = check_match(ii,jj)
# some species name starts with "an ", it will be removed
    if cm2 == 0:
        if ii.startswith("an "):
            ii = ii[3:]
        if jj.startswith("an "):
            jj = jj[3:]
        cm2 = check_match(ii,jj)
    if cm2 == 0:
        if "(2E)" in ii.lower():
            ii = ii.lower().replace("(2E)","trans")
        if "(2E)" in jj.lower():
            jj = jj.lower().replace("(2E)","trans")
        cm2 = check_match(ii,jj)
    #check it
    #if cm2==0:
        #if "quinol" in ii:
        #    ii=ii.replace("quinol", "quinone")
        #if "quinol" in jj:
        #    jj=jj.replace("quinol", "quinone")
        #cm2=check_match(ii,jj)
# sometimes acyl-carrier-protein is named like ACP, [acp], [acyl-carier-protein] or [acyl-carier protein]
    if cm2 == 0:
        if "ACP" in ii:
            ii = ii.replace("ACP", "[acyl-carrier protein]")
        if "[acp]" in ii:
            ii = ii.replace("[acp]", "[acyl-carrier protein]")
        if "ACP" in jj:
            jj = jj.replace("ACP", "[acyl-carrier protein]")
        if "[acp]" in jj:
            jj = jj.replace("[acp]", "[acyl-carrier protein]")
        cm2 = check_match(ii,jj)
# remove " " and "-", e.g. ii = "...[acyl-carrier protein]" and jj="...[acyl-carier-protein]", they don't match, after remove they match [acylcarrierprotein]
    if cm2 == 0:
        if " " in ii:
            ii = ii.replace(" ","")
        if "-" in ii:
            ii = ii.replace("-","")
        if " " in jj:
            jj = jj.replace(" ","")
        if "-" in jj:
            jj = jj.replace("-","")
        if "(" in ii:
            ii = ii.replace("(","")
        if ")" in ii:
            ii = ii.replace(")","")
        if "(" in jj:
            jj = jj.replace("(","")
        if ")" in jj:
            jj = jj.replace(")","")
        cm2 = check_match(ii,jj)
    if cm2 == 0:
        #if "alpha" in ii:
        #    ii=ii.replace("alpha","")
        #if "alpha" in jj:
        #    jj=jj.replace("alpha","")
        #if "beta" in ii:
        #    ii=ii.replace("beta","")
        #if "beta" in jj:
        #    jj=jj.replace("beta","")
        #if "R" in ii:
        #    ii=ii.replace("R","")
        #if "R" in jj:
        #    jj=jj.replace("R","")
        if "N" in ii:
            ii = ii.replace("N","")
        if "N" in jj:
            jj = jj.replace("N","")
        if "D" in ii:
            ii = ii.replace("D","")
        if "D" in jj:
            jj = jj.replace("D","")
        if "L" in ii:
            ii = ii.replace("L","")
        if "L" in jj:
            jj = jj.replace("L","")
        #if "S" in ii:
        #    ii=ii.replace("S","")
        #if "S" in jj:
        #    jj=jj.replace("S","")
        if "Z" in ii:
            ii = ii.replace("Z","")
        if "Z" in jj:
            jj = jj.replace("Z","")
        cm2 = check_match(ii,jj)
    if cm2 == 0:
        if "hexadecanoyl" in ii.lower():
            ii=ii.lower().replace("hexadecanoyl","palmitoyl")
        if "hexadecanoyl" in jj.lower():
            jj=jj.lower().replace("hexadecanoyl","palmitoyl")
        if "hexadecenoyl" in ii.lower():
            ii=ii.lower().replace("hexadecenoyl","palmitoleoyl")
        if "hexadecenoyl" in jj.lower():
            jj=jj.lower().replace("hexadecenoyl","palmitoleoyl")
        if "butyryl" in ii.lower():
            ii=ii.lower().replace("butyryl","butanoyl")
        if "butyryl" in jj.lower():
            jj=jj.lower().replace("butyryl","butanoyl")
        if "tetradecanoyl" in ii.lower():
            ii=ii.lower().replace("tetradecanoyl","myristoyl")
        if "tetradecanoyl" in jj.lower():
            jj=jj.lower().replace("tetradecanoyl","myristoyl")
        cm2 = check_match(ii,jj)
    return cm2


def check_match( ii, jj):
    #This function checks if the given species names match. If a match is made, the integer 1 is returned, otherwise 0.
    cm = 0
    if ii.lower() == jj.lower():
        cm = 1
    return cm


def get_syn_from_chebi (chebiId):
    #This function extracts the synonyms stored in the respective ChEBI-Entry by using libchebipy package. The input is the ChEBI-ID of the species. The synonyms are returned.
    resi = []
    try:
        for i in lbc.ChebiEntity(chebiId).get_names():
            if "'_Name__name': u'" in str(i):
                resi.append(str(i).split("'_Name__name': u'")[1].split("',")[0])
            if "'_Name__name': u\"" in str(i):
                resi.append(str(i).split("'_Name__name': u\"")[1].split("',")[0])
        return resi
    except lbc._chebi_entity.ChebiException:
        pass
    return resi


def get_syn_from_KEGGID (KEGGID):
    #This function extracts the synonyms stored in the respective KEGG-Entry by using biopython package. The input is the KEGG-ID of the species. The synonyms are returned.
    synkegg = []
    try:
        ccc = REST.kegg_get("compound:" + KEGGID).read()
        
        for line in ccc.rstrip().split("\n"):
            section = line[:12].strip()
            if not section == "":
                current_section = section
            if current_section == "NAME":
                name = line[12:]
                if name[-1] == ";":
                    synkegg.append(name[:-1])
                else:
                    synkegg.append(name)
        return synkegg
    except urllib2.HTTPError:
        pass 
    return synkegg


def get_syn_from_cid(cid):
    #This function extracts the synonyms stored in the respective PubChem-Entry by using pubchempy package. The input is the PubChem-ID of the species. The synonyms are returned.
    synpcp=[]
    try:
        synpcp = pcp.Compound.from_cid(cid).synonyms
        #print synpcp
        return synpcp
    except urllib2.HTTPError:
        pass
    except urllib2.URLError:
        get_syn_from_cid(cid)
    except pcp.NotFoundError:
        pass
    return synpcp


def synonym_matching(restrict_can, i_s_match, j_s_match,species_ids, identifier_index,report, window1, progbar, progress_label, species_matches):
    #This function match synonyms of the restricted candidates. The synonyms are extracted from KEGG, ChEBI and PubChem, if an ID of this databases exists in the candidate species note block. Otherwise no synonyms can be determined. 
    #If the synonym not matches than the synonym will be adapted by the adapt_name function. 
    timea = datetime.datetime.now()
    l=[]
    if identifier_index == 0:
        print "KEGG-SYN-Matching: " + "\n"
    if identifier_index == 1:
        print "ChEBI-SYN-Matching: " + "\n"
    if identifier_index == 2:
        print "PubChem-SYN-Matching: " + "\n"
    current_i1_nr = 0
    for i in restrict_can:
        current_i1_nr = current_i1_nr+1
        progress_nr = (100 / float(len(restrict_can))) * current_i1_nr
        for j in restrict_can[i]:
            progress_label.config(text ="SBML1-Species-Name\n\n"+i.getName() + "\n\nis tried to match with SBML2-Species-Synonym\n\n"+j.getName())
            progbar["value"] = progress_nr
            window1.update()
            syn=[]
            if identifier_index == 0:
                if "<p>KEGG: " in j.getNotesString():
                    if i.getId() not in i_s_match:
                        syn=  get_syn_from_KEGGID(j.getNotesString().split("<p>KEGG: ")[1].split("</p>")[0])
            if identifier_index == 1:
                if "<p>CHEBI: " in j.getNotesString():
                    if i.getId() not in i_s_match:
                        syn = get_syn_from_chebi(j.getNotesString().split("<p>CHEBI: ")[1].split("</p>")[0])
            if identifier_index == 2:
                if "<p>PUBCHEM: " in j.getNotesString():
                    if i.getId() not in i_s_match:
                        syn = get_syn_from_cid(j.getNotesString().split("<p>PUBCHEM: ")[1].split("</p>")[0])
            for e in syn:
                if i.getName() == e:
                    if j.getId() not in species_ids[i.getId()]:
                        species_ids[i.getId()].append(j.getId())
                    write_species(i,j,e, identifier_index, species_matches)
                    if i.getId() not in i_s_match:
                        i_s_match.append(i.getId())
                    if j.getId() not in j_s_match:
                        j_s_match.append(j.getId())
                    break
                else:
                    cm = adapt_name(i,e,1)
                    if cm == 1:
                        if j.getId() not in species_ids[i.getId()]:
                            species_ids[i.getId()].append(j.getId())
                        write_species(i,j,e,identifier_index, species_matches)
                        if i.getId() not in i_s_match:
                            i_s_match.append(i.getId())
                        if j.getId() not in j_s_match:
                            j_s_match.append(j.getId())
    #window1.quit()
    window1.destroy()
    l.append(i_s_match)
    l.append(j_s_match)
    l.append(species_ids)
    timeaa = datetime.datetime.now()-timea
    report.write("\t\t- process-time: " + str(timeaa)+"\n")
    return l


def restrict_candidates_for_syn_matching(sbml1_Model,sbml2_Model,i_s_match, report, window1, progbar, progress_label):
    #This function determines candidates for the synonym matching. Candidates are those SBML2 species that have the same formula to the SBML1 species. For each SBML1 species a set of candidats is determined, if it has not yet been assigned. The function returns a list with the SBML1 species as key and the SBML2 species candidates as value.
    timea = datetime.datetime.now()
    restrict_can = {}
    current_i1_nr = 0
    for i in sbml1_Model.getListOfSpecies():
        current_i1_nr = current_i1_nr+1
        #print current_i1_nr 
        progress_nr = (100 / float(len(sbml1_Model.getListOfSpecies()))) * current_i1_nr
        if i not in i_s_match:
            progress_label.config(text = "SBML1-Species \n\n" + i.getName() +  "\n\ncandidates will be determined!")
            progbar["value"] = progress_nr
            window1.update()
            i_cand = []
            for j in sbml2_Model.getListOfSpecies():
                #Restriction to the same formula
                if "<p>FORMULA: " in i.getNotesString() and "<p>FORMULA: " in j.getNotesString():
                    if i.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0] == j.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0]:
                        if j not in i_cand:
                            i_cand.append(j)
            restrict_can.update({i:i_cand})
    #window1.quit()
    window1.destroy()
    print len(restrict_can)
    timeaa = datetime.datetime.now() - timea
    report.write("\t\t- process-time: " + str(timeaa) + "\n")
    return restrict_can
    

def id_match_counter(ids):
    a = 0
    for identifier in ids:
        if ids[identifier] != []:
            a = a + 1
    return a
    

def reaction_comparison(sbml1_Model, sbml2_Model, i_r_match, j_r_match, species_ids,wrplist,wplist,wrlist,crplist,reaction_ids, report, progress_label, progbar, window1, reaction_matches, wrp, wp, wr):
    #This function identifies the reactions that will be studied in more detail. Identified reactions are reactions that match in number of reactants and products. The components of these are checked with the function check_reactants_and_products.
    timea = datetime.datetime.now()
    count_list = [0,0,0,0] # count_list=[count123,count12,count13,count1]
    count_list_no_match = [0,0,0] #count_list_no_match=[count_products_reactants_missing, count_products_missing, count_reactants_missing]
    o = []
    r = 0
    current_i1_nr = 0
    for i in sbml1_Model.getListOfReactions():
        current_i1_nr=current_i1_nr+1
        #print current_i1_nr 
        progress_nr = (100 / float(len(sbml1_Model.getListOfReactions()))) * current_i1_nr
        for j in sbml2_Model.getListOfReactions():
            progress_label.config(text = "SBML1-Reaction\n\n" + i.getName() + "\n\nis tried to match with SBML2-Reaction\n\n" + j.getName())
            progbar["value"] = progress_nr
            window1.update()
            if len(i.getListOfReactants()) == len(j.getListOfReactants()) and len(i.getListOfProducts()) == len(j.getListOfProducts()):
                #if i.getId not in i_r_match and j.getId not in j_r_match:
                x = check_reactants_and_products(i, j, i_r_match, j_r_match, species_ids,count_list ,count_list_no_match,wrplist,wplist,wrlist,crplist ,reaction_ids, reaction_matches, sbml1_Model,sbml2_Model, wrp, wp, wr)
                i_r_match = x[0]
                j_r_match = x[1]
                reaction_ids = x[2]
                count_list = x[3]
                count_list_no_match = x[4]
                r = r + x[5]
                wrplist = x[6]
                wplist = x[7]
                wrlist = x[8]
    print "reactions_matched: " + str(r)
    window1.destroy()
    o.append(i_r_match)
    o.append(j_r_match)
    o.append(reaction_ids)
    o.append(count_list)
    o.append(count_list_no_match)
    o.append(wrplist)
    o.append(wplist)
    o.append(wrlist)
    o.append(crplist)
    timeaa = datetime.datetime.now() - timea
    report.write("\t\t- process-time: " + str(timeaa) + "\n")
    return o
    

def check_reactants_and_products_reversed(i,j, species_ids):
    #This function checks whether the reactants and products of the SBML1 reaction match those of the SBML2 reaction in reversed order. So reactants of SBML1 reaction could be the products of SBML2 reaction . For this the assigned species are used.
    r = len(i.getListOfReactants())
    p = len(i.getListOfProducts())
    stoich = []
    for l in i.getListOfReactants():
        for k in j.getListOfProducts():
            for m in species_ids[l.getSpecies()]:
                if k.getSpecies() == m:
                    r = r - 1
                    if l.getStoichiometry() == k.getStoichiometry():
                        stoich.append(1)

    for u in i.getListOfProducts():
        for v in j.getListOfReactants():
            for n in species_ids[u.getSpecies()]:
                if v.getSpecies() == n:
                    p = p - 1
                    if u.getStoichiometry() == v.getStoichiometry():
                        stoich.append(1)
    return r,p,stoich


def check_reactants_and_products(i, j, i_r_match, j_r_match, species_ids,count_list, count_list_no_match, wrplist,wplist,wrlist,crplist,reaction_ids, reaction_matches, sbml1_Model,sbml2_Model, wrp, wp, wr):
    #This function checks whether the reactants and products of the SBML1 reaction match those of the SBML2 reaction. For this the assigned species are used. Reversibility and stoichiometry are also checked. When the reactions match, a classification of the reactions takes place.
    x = []
    r = len(i.getListOfReactants())
    p = len(i.getListOfProducts())
    revers = 0
    reversedrp = 0
    stoich = []
    match = 0
    if i.getReversible() == j.getReversible():
        revers = 1
    for l in i.getListOfReactants():
        for k in j.getListOfReactants():
            for m in species_ids[l.getSpecies()]:
                if k.getSpecies() == m:
                    r = r - 1
                    if l.getStoichiometry() == k.getStoichiometry():
                        stoich.append(1)

    for u in i.getListOfProducts():
        for v in j.getListOfProducts():
            for n in species_ids[u.getSpecies()]:
                if v.getSpecies() == n:
                    p = p - 1
                    if u.getStoichiometry() == v.getStoichiometry():
                        stoich.append(1)
    if p != 0 and r != 0:
        r,p,stoich = check_reactants_and_products_reversed(i,j,species_ids)
        if p == 0 and r == 0:
            reversedrp = 1
#all reactants and products available
    if p == 0 and r == 0:
        match = 1
        #if j.getId() not in reaction_ids[i.getId()]:
        reaction_ids[i.getId()].append(j.getId())
        #perfect match (stoichiometry and reversibility are the same) reversed
        if revers == 1 and reversedrp == 1 and sum(stoich) == (len(i.getListOfReactants()) + len(i.getListOfProducts())):
            reaction_matches.write("Reaction-Qualification: " + "\t" + "r123\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
            crplist.append(get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,5))
        # match (stoichiometry  are the same and reversibility not) reversed
        if revers != 1 and reversedrp == 1 and sum(stoich) == (len(i.getListOfReactants()) + len(i.getListOfProducts())):
            reaction_matches.write("Reaction-Qualification: " + "\t" + "r13\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
            crplist.append(get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,3))
        # match (reversibility are the same and stoichiometry not) reversed
        if revers == 1 and reversedrp == 1 and sum(stoich) != (len(i.getListOfReactants()) +len(i.getListOfProducts())):
            reaction_matches.write("Reaction-Qualification: " + "\t" + "r12\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
            crplist.append(get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,6))
        #only reactants and products matches but stoichiometry and reversibility differs and reverserd
        if revers != 1 and reversedrp == 1 and sum(stoich) != (len(i.getListOfReactants()) + len(i.getListOfProducts())):
            reaction_matches.write("Reaction-Qualification: " + "\t" + "r1\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
            crplist.append(get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,4))
        #perfect match
        if revers == 1 and sum(stoich) == (len(i.getListOfReactants()) + len(i.getListOfProducts())) and reversedrp != 1:
            count_list[0] = count_list[0] + 1
            reaction_matches.write("Reaction-Qualification: " + "\t" + "123\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
        #stoichiometry differs
        if revers == 1 and sum(stoich) != (len(i.getListOfReactants()) + len(i.getListOfProducts())) and reversedrp != 1 :
            count_list[1] = count_list[1] + 1
            reaction_matches.write("Reaction-Qualification: " + "\t" + "12\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
        #reversibility differs
        if revers == 0 and sum(stoich) == (len(i.getListOfReactants()) + len(i.getListOfProducts())) and reversedrp != 1:
            count_list[2] = count_list[2] + 1
            reaction_matches.write("Reaction-Qualification: " + "\t" + "13\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
        #only reactants and products matches but stoichiometry and reversibility differs
        if revers == 0 and sum(stoich) != (len(i.getListOfReactants()) + len(i.getListOfProducts())) and reversedrp != 1:
            count_list[3] = count_list[3] + 1
            reaction_matches.write("Reaction-Qualification: " + "\t" + "1\n")
            write_reaction_in_file(sbml1_Model,sbml2_Model,i,j,reaction_matches)
        if i.getId() not in i_r_match:
            i_r_match.append(i.getId())
        if j.getId() not in j_r_match:
            j_r_match.append(j.getId())
#some(i) products and some (j) reactants missing.    i between 0<i<len(productlist),    j beetween 0<j<len(reactantlist)
    elif p > 0 and p < len(i.getListOfProducts()) and r > 0 and r < len(i.getListOfReactants()):
        if i.getId() not in wrplist:
            wrplist[i.getId()] = [get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,0)]
        else:
            wrplist[i.getId()].append(get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,0))
#some(i) products missing    zwischen 0<i<len(productlist)
    elif p > 0 and p < len(i.getListOfProducts()) and r == 0:
        if i.getId() not in wplist:
            wplist[i.getId()] = [get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,1)]
        else:
            wplist[i.getId()].append(get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,1))
#some(j) reactants missing    zwischen 0<j<len(reactantlist)
    elif r > 0 and r < len(i.getListOfReactants()) and p == 0:
        if i.getId() not in wrlist:
            wrlist[i.getId()] = [get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,2)]
        else:
            wrlist[i.getId()].append(get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,2))
    x.append(i_r_match)
    x.append(j_r_match)
    x.append(reaction_ids)
    x.append(count_list)
    x.append(count_list_no_match)
    x.append(match)
    x.append(wrplist)
    x.append(wplist)
    x.append(wrlist)
    x.append(crplist)
    return x


def write_reaction_in_file(sbml1_Model,sbml2_Model,i,j, file1):
    #This function generates the reaction report entries. 
    file1.write("\t" + i.getId() + "\t" + i.getName() + "\t" + "(reversible: " + str(i.getReversible()) + ") :" + "\n\t\t")
    for l in i.getListOfReactants():
        file1.write(str(l.getStoichiometry()) + "\t" + l.getSpecies() + "\t")
        for r in sbml1_Model.getListOfSpecies():
            if r.getId() == l.getSpecies():
                file1.write("(" + r.getName() + ")" + "\t")
    file1.write("|" + "\t")
    for u in i.getListOfProducts():
        file1.write(str(u.getStoichiometry()) + "\t" + u.getSpecies() + "\t")
        for r in sbml1_Model.getListOfSpecies():
            if r.getId() == u.getSpecies():
                file1.write("(" + r.getName() + ")" + "\t")
    file1.write("\n")
    file1.write(i.getNotesString() + "\n")
    file1.write("\t" + j.getId() + "\t" + j.getName() + "\t" + "(reversible: " + str(j.getReversible()) + ") :" + "\n\t\t")
    for k in j.getListOfReactants():
        file1.write(str(k.getStoichiometry()) + "\t" + k.getSpecies() + "\t")
        for r in sbml2_Model.getListOfSpecies():
            if r.getId() == k.getSpecies():
                file1.write("(" + r.getName() + ")" + "\t")
    file1.write("|" + "\t")
    for v in j.getListOfProducts():
        file1.write(str(v.getStoichiometry()) + "\t" + v.getSpecies() + "\t")
        for r in sbml2_Model.getListOfSpecies():
            if r.getId() == v.getSpecies():
                file1.write("(" + r.getName() + ")" + "\t")
    file1.write("\n")
    file1.write(j.getNotesString() + "\n")
    file1.write("\n")


def get_reaction_string(sbml1_Model,sbml2_Model,i,j,r,p,index):
    # This function generates reaction strings 
    s = ""
    if index == 0:
        s = str(r) + "/" + str(len(i.getListOfReactants())) + " reactants missing ||" + str(p) + "/" + str(len(i.getListOfProducts())) + " products missing\n"
    if index == 1:
        s = str(p) + "/" + str(len(i.getListOfProducts())) + " products missing\n"
    if index == 2:
        s = str(r) + "/" + str(len(i.getListOfReactants())) + " reactants missing\n"
    if index == 3:
        s = "reactants and products changed - Stoichiometry match\n"
    if index == 4:
        s = "reactants and products changed - Stoichiometry match\n"
    if index == 5:
        s = "reactants and products changed - Stoichiometry match\n"
    if index == 6:
        s = "reactants and products changed - Stoichiometry match\n"
    s = "\t" + i.getId() + "\t" + i.getName() + "\t" + "(reversible: " + str(i.getReversible())+") :" + "\n\t\t"
    for l in i.getListOfReactants():
        s = s + str(l.getStoichiometry()) + "\t" + l.getSpecies() + "\t"
        for r in sbml1_Model.getListOfSpecies():
            if r.getId() == l.getSpecies():
                s = s + "(" + r.getName() + ")" + "\t"
    s = s + "|" + "\t"
    for u in i.getListOfProducts():
        s = s + str(u.getStoichiometry()) + "\t" + u.getSpecies() + "\t"
        for r in sbml1_Model.getListOfSpecies():
            if r.getId() == u.getSpecies():
                s = s + "(" + r.getName() + ")" + "\t"
    s = s + "\n"
    s = s + i.getNotesString() + "\n"
    s = s + "\t" + j.getId() + "\t" + j.getName() + "\t" + "(reversible: " + str(j.getReversible()) + ") :" + "\n\t\t"
    for k in j.getListOfReactants():
        s = s + str(k.getStoichiometry()) + "\t" + k.getSpecies() + "\t"
        for r in sbml2_Model.getListOfSpecies():
            if r.getId() == k.getSpecies():
                s = s + "(" + r.getName() + ")" + "\t"
    s = s + "|" + "\t"
    for v in j.getListOfProducts():
        s = s + str(v.getStoichiometry()) + "\t" + v.getSpecies() + "\t"
        for r in sbml2_Model.getListOfSpecies():
            if r.getId() == v.getSpecies():
                s = s + "(" + r.getName() + ")" + "\t"
    s = s + "\n"
    s = s + j.getNotesString() + "\n"
    s = s + "\n"
    return s


def get_pwytools_syn_data (file_path):
    file3 = open(file_path,"r")
    syn_data = []
    for line in file:
        syn_data.append(line)
    file3.close()
    return syn_data


def ptoos_syn_matching(sbml1_Model,sbml2_Model,syn_data,i_s_match, j_s_match,species_ids, species_matches, report):
    format_list = ["&",";","<i>","</i>","<I>","</I>", "<sub>","</sub>","<SUB>","</SUB>","<sup>","</sup>","<SUP>","</SUP>"]
    timea = datetime.datetime.now()
    l = []
    for i in sbml1_Model.getListOfSpecies():
        if i.getId() not in i_s_match:
            ii = i.getName()
            if ii.startswith("a "):
                ii = ii[2:]
            if ii.startswith("an "):
                ii = ii[3:]
            if " " in ii:
                ii = ii.replace(" ","")
            if "-" in ii:
                ii = ii.replace("-","")

            for  j in sbml2_Model.getListOfSpecies():
                if "<p>BIOCYC: " in j.getNotesString():
                    #print "PWYTools-Synonym-Matching"
                    #print j.getNotesString().split("<p>BIOCYC: ")[1].split("</p>")[0]
                    if "|" in j.getNotesString().split("<p>BIOCYC: ")[1].split("</p>")[0]:
                        frameid = j.getNotesString().split("<p>BIOCYC: ")[1].split("</p>")[0].split("|")[1]
                    else:
                        frameid = j.getNotesString().split("<p>BIOCYC: ")[1].split("</p>")[0]
                    for line in syn_data:
                        if line.split("\t")[0] == frameid:
                            if line.split("\t")[1] != '':
                                s = line.split("\t")[1].split("$")
                                for u in range(len(s)):
                                    for z in format_list:
                                        if z in s[u]:
                                            s[u] = s[u].replace(z,"")
                                for e in s:
                                    ee = e
                                    if ee.startswith("a "):
                                        ee = ee[2:]
                                    if ee.startswith("an "):
                                        ee = ee[3:]
                                    if " " in ee:
                                        ee = ee.replace(" ","")
                                    if "-" in ee:
                                        ee = ee.replace("-","")

                                    if ee.lower() == ii.lower():
                                        species_ids[i.getId()].append(j.getId())
                                        species_matches.write(i.getId() + "\t" + i.getName() + "\t" + i.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0] + "\t"+"||" + "\t" + j.getId() + "\t" + j.getName() + "\t" + j.getNotesString().split("<p>FORMULA: ")[1].split("</p>")[0] + "\t"+"||" + "\t" + "ptools-SYN-NAME-MATCH" + "\t" + e + "\n")
                                        if i.getId() not in i_s_match:
                                            i_s_match.append(i.getId())
                                        if j.getId() not in j_s_match:
                                            j_s_match.append(j.getId())
                                        break
                                break
    l.append(i_s_match)
    l.append(j_s_match)
    l.append(species_ids)
    timeaa=datetime.datetime.now()-timea
    report.write("\t\t- process-time: " + str(timeaa) + "\n")
    return l


def jaccardc(m01, m10, m11):
    # returns the jaccard similarity coefficient
    #m01 = reactions only in SBML2
    #m10 = reactions only in SBML1
    #m11 = reactions in both
    return m11 * (m01 + m10 + m11) ** -1


def jaccardd(m01,m10,m11):
    #return the jaccard distanz
    #m01 = reactions only in SBML2
    #m10 = reactions only in SBML1
    #m11 = reactions in both
    
    return (m01 + m10) * (m01 + m10 + m11) ** -1            


#####################################################################################################
# Placeholder for file paths
#Path of SBML file 1
sbml1 = ""
#Path of SBML file 2
sbml2 = ""
#Directory for the storage of report files
reportdir = ""
#####################################################################################################
# GUI myButtonSettings
button_font = "Arial 12 bold"
button_borderwidth = 4
button_bg = "light grey"
button_cursor = "hand2"
button_overrelief = "sunken"

# GUI myLabelSettings
label_font = "Arial 10"
label_wraplength = 800

# GUI myFrameLabelSettings
framelabel_font = "Arial 12 bold"
framelabel_borderwidth = 4 
framelabel_relief = "groove"

#####################################################################################################
# GUI Main Window
main_window = Tk()
main_window.title("SBMLComp")
# GUI Fixed window width and height, centered on the screen
window_width = 900
window_height = 600
print main_window.winfo_screenwidth()
print main_window.winfo_screenheight()
window_pos_x = (main_window.winfo_screenwidth() /2) - (window_width / 2)
window_pos_y = (main_window.winfo_screenheight() /2) - (window_height / 2)
main_window.geometry('%dx%d+%d+%d' % (window_width, window_height, window_pos_x, window_pos_y))
main_window.resizable(0,0)

# GUI Program start information
main_frame = Frame(main_window)
main_frame.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
prog_label = Label(main_frame, text = "##########################################################################################################################\n\n SBMLComp\n\n ###########################################################################################################################\n\nAuthor: Eric Witt\n\n###########################################################################################################################",  font="Arial 10 bold", relief = "groove", borderwidth = 4)
prog_label.pack(fill = X, padx = 5, pady = 10)

info_label = Label(main_frame, text = "SBMLComp is a comparison tool for SBML files at reaction level. This can be used for comparative analyses of different organisms in relation to their metabolic similarities and potentials.\n\nFor the comparison first the Species of both SBML files are matched, for this the database IDs, names and Inchis of the Species existing in the SBML file are used. Since often different naming conventions are used for the species, the existing database IDs are used to extract synonyms from the databases to increase the matching success.\n\nAfter matching the species, the reactions are matched. Therefore the reactions where the reactant and product numbers are identical are compared. The matched species are now used to compare the reactions and determine the identity of the reactions. Furthermore, the reactions are classified in a match. Therefore the reversibility and stoichiometry is checked. A perfect match would mean that the reactants, products, reversibility and stoichiometry of the compared reactions match. If the reactions do not match, it will also be checked whether they are simply upside down. This means that the reactants of one reaction are the products of the other reaction and the products are the reactants. \n\nAfter Matching the reactions of the SBML files, the Jaccard similarity coefficient is calculated as a measure for the similarity of the SBML files.", font="Arial 10 bold",  wraplength = 700, relief = "groove", borderwidth = 4)
info_label.pack(fill = X, padx = 5, pady = 10, ipadx = 10, ipady = 10)

def choose_files():
    main_frame.pack_forget()
   
    # GUI Generate Tabs
    noteb = ttk.Notebook(main_window)
    noteb.pack(fill=BOTH, expand = 1, padx = 5, pady = 5)
    noteb.pressed_index = None
    tab1 = Frame(noteb)
    tab1.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
    noteb.add(tab1, text = "Choose Files")
    
    # GUI Information for loading sbml files
    info_SBML2_framelabel = LabelFrame(tab1, text = "General Informations", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
    info_SBML2_framelabel.pack(fill = X, padx = 5, pady = 10)
    info_SBML2_label = Label(info_SBML2_framelabel, text = "SBMLComp contains some Synonym-Matching functions which uses the KEGG-, ChEBI-, and Pubchem-identifier to extract the synonyms of the respective compound from the associoated database. The Synonym-Matching functions only works if the SBML file which contains the identifiers is loaded as SBML file 2.",  font = label_font, wraplength = label_wraplength) 
    info_SBML2_label.pack(fill = X, padx = 5, pady = 5)
    
    # GUI LOAD SBML1
    def load_SBML1():
        # Function for loading SBML file 1 by user in the GUI. The selected file path is displayed in the GUI.
        sbml1 = str(tkFileDialog.askopenfilename(initialdir = os.getcwd, title = "Select the SBML file 1 (.xml)!", filetypes = (("SBML files","*.xml"),("all files","*.*file"))))
        load_SBML1_label.config(text = sbml1)
        load_SBML2_button.config(state = "normal")

    
    load_SBML1_framelabel = LabelFrame(tab1, text = "Load SBML file 1", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
    load_SBML1_framelabel.pack(fill = X, padx = 5, pady = 10)
    load_SBML1_label = Label(load_SBML1_framelabel, text = sbml1,  font = label_font, wraplength = label_wraplength) 
    load_SBML1_label.pack(fill = X, padx = 5, pady = 5)
    load_SBML1_button = Button(load_SBML1_framelabel, text = "Load", command = load_SBML1,  font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
    load_SBML1_button.pack(fill = X, padx = 5, pady = 5)
    
    # GUI LOAD SBML2
    def load_SBML2():
        # Function for loading SBML file 2 by user in the GUI. The selected file path is displayed in the GUI.
        sbml2 = str(tkFileDialog.askopenfilename(initialdir = os.getcwd, title = "Select the SBML file 2 (.xml)!", filetypes = (("SBML files","*.xml"),("all files","*.*file"))))
        load_SBML2_label.config(text = sbml2)
        select_directory_button.config(state = "normal")

    
    load_SBML2_framelabel = LabelFrame(tab1, text = "Load SBML file 2", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
    load_SBML2_framelabel.pack(fill = X, padx = 5, pady = 10)
    load_SBML2_label = Label(load_SBML2_framelabel, text = sbml2,  font = label_font, wraplength = label_wraplength) 
    load_SBML2_label.pack(fill = X, padx = 5, pady = 5)
    load_SBML2_button = Button(load_SBML2_framelabel, text = "Load", command = load_SBML2,  font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth, state = DISABLED)
    load_SBML2_button.pack(fill = X, padx = 5, pady = 5)
    
    # GUI SELECT REPORT LOCATION
    def select_report():
        # Function for selecting the directory for storage of the report files. The selected directory path is displayed in the GUI.
        reportdir = str(tkFileDialog.askdirectory())
        select_directory_label.config(text = reportdir)
        analyze_sbml_files_button.config(state = "normal")


    select_directory_framelabel = LabelFrame(tab1, text = "Select Location For Report File Storing", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
    select_directory_framelabel.pack(fill = X, padx = 5, pady = 10)
    select_directory_label = Label(select_directory_framelabel, text = reportdir, font = label_font, wraplength = label_wraplength) 
    select_directory_label.pack(fill = X, padx = 5, pady = 5)
    select_directory_button = Button(select_directory_framelabel, text = "Select directory", state = DISABLED, command = select_report, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
    select_directory_button.pack(fill = X, padx = 5, pady = 5)

    def analyze_SBMl_files():
        #Function for analyzing the SBML files. The Species note section will be checked for possible identifiers (database-IDs & Inchis) for the comparison of the SBML files. The identifier type and the number of identifiers is displayed in the GUI. The user select the mode for comparison which are possible with the founded informations.
        
        # Disable "Choose Files" Buttons
        load_SBML1_button.config(state = DISABLED)
        load_SBML2_button.config(state = DISABLED)
        select_directory_button.config(state = DISABLED)
        analyze_sbml_files_button.config(state = DISABLED)

        # GUI New Tab "General Informations"
        tab2 = Frame(noteb)
        tab2.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
        noteb.add(tab2, text = "Analysis Results & Mode Selection")
        noteb.select(tab2)
        
        # GUI add objects for a scrollable tab
        mycanvas = Canvas(tab2)
        mycanvas.pack(side = "left", fill = BOTH, expand = 1, padx = 60)

        unacc_frame = Frame(mycanvas)
        unacc_frame.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
        mycanvas.create_window(0, 0, window = unacc_frame, anchor = "nw")
        
        #Load SBML-file. libSBML package is used for reading the SBML files
        read = SBMLReader()
        sbml1_SBML = read.readSBML(load_SBML1_label["text"])
        sbml2_SBML = read.readSBML(load_SBML2_label["text"])
        #Get model-data from SBML-file
        sbml1_Model = sbml1_SBML.getModel()
        sbml2_Model = sbml2_SBML.getModel()
        #report
        report_path = select_directory_label["text"]+"/01_matching_report.txt"
        report = open(report_path,"w")
        report.write("SBML1:\t" + load_SBML1_label["text"].split("/")[-1]+"\n")
        report.write("SBML2:\t" + load_SBML2_label["text"].split("/")[-1]+"\n")
        report.write("Analysing SBML files:\n")
        tt = datetime.datetime.now()
        # species-id-check for kegg-id
        report.write("\t1. Check SBML1-Species-Ids for KEGG-Id:\n")
        # dict for species_id_kegg_id
        kegg_id1 = {}
        kegg_id2 = {}
        
        answer1 = messagebox.askyesno("Question","Sometimes database identifiers are stored in the Species ID.\nTherefore it is recommended to search the Species-IDs for Kegg-IDs.\nDo you want to allow this?")
        
        #KEGG-ID-CHECK in SPECIES-ID SBML1
        if answer1 == True:
            window1 = Tk()
            window1.title("SBML1 Species-ID Check Progress!")
            progress_nr = 0
            progress_label  =Label(window1, text = 'SBML1 Species-ID check starts!', width = 75)
            progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate",maximum = 100, value = progress_nr)
            progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
            progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
            def id_check1():
                current_i1_nr = 0
                for i in sbml1_Model.getListOfSpecies():
                    current_i1_nr = current_i1_nr+1
                    progress_nr = (100/float(len(sbml2_Model.getListOfSpecies())))*current_i1_nr
                    progress_label.config(text = "Species-ID " + i.getId() + " is checked!")
                    progbar["value"] = progress_nr
                    window1.update()
                    res = species_id_kegg_check(i.getId())
                    #print res
                    if res != "":
                        kegg_id1.update({i: res})
                report.write("\t\t- KEGG-Ids found:\t" + str(len(kegg_id1)) + "\n")
                window1.destroy()
            
            id_check1()
        else:
            report.write("\t\t\t- KEGG-ID-CHECK disabled\n")
        #KEGG-ID-CHECK in SPECIES-ID SBML2
        report.write("\t2. Check SBML2-Species-Ids for KEGG-Id:\n")
        if answer1 == True:
            window1 = Tk()
            window1.title("SBML2 Species-ID Check Progress!")
            progress_nr = 0
            progress_label = Label(window1, text = 'SBML2 Species-ID check starts!', width = 75)
            progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate",maximum = 100, value = progress_nr)
            progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
            progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
            def id_check2():
                current_i2_nr = 0
                for i in sbml2_Model.getListOfSpecies():
                    current_i2_nr = current_i2_nr+1
                    progress_nr = (100/float(len(sbml2_Model.getListOfSpecies())))*current_i2_nr
                    progress_label.config(text = "Species-ID " + i.getId() + " is checked!")
                    progbar["value"] = progress_nr
                    window1.update()
                    res = species_id_kegg_check(i.getId())
                    #print res
                    if res != "":
                        kegg_id2.update({i: res})
                report.write("\t\t- KEGG-Ids found:\t" + str(len(kegg_id2)) + "\n")
                window1.destroy()
            
            id_check2()
        else:
            report.write("\t\t\t- KEGG-ID-CHECK disabled\n")
        
        # check notes of DB-Identifiers
        report.write("\t3. Check SBML1-Model for Ids (KEGG,ChEBI,PubChem, Inchi, BioCyc, HMDB, CAS, Chemspider, Metabolights):\n")
        #Get Identifiers from SBML file 1
        window1 = Tk()
        window1.title("Check SBML file 1 for identifiers!")
        progress_nr = 0
        progress_label = Label(window1, text = 'Check starts!', width = 75)
        progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate",   maximum = 100, value = progress_nr)
        progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
        progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
        identifier_sbml1 = check_possible_identifier_to_comp_from_notes_sbml(sbml1_Model, report, window1, progbar, progress_label)
        
        #write KEGG-IDs found in Species-ID SBML1 to identifier_sbml1
        if len(kegg_id1) != 0 :
            for i in kegg_id1:
                identifier_sbml1[0].update({i:kegg_id1[i]})
        #write report for SBML1
        report.write("\t\t- KEGG:\t" + str(len(identifier_sbml1[0])) + "\n")
        report.write("\t\t- ChEBI:\t" + str(len(identifier_sbml1[1])) + "\n")
        report.write("\t\t- PubChem:\t" + str(len(identifier_sbml1[2])) + "\n")
        report.write("\t\t- Inchi:\t" + str(len(identifier_sbml1[3])) + "\n")
        report.write("\t\t- BioCyc:\t" + str(len(identifier_sbml1[4])) + "\n")
        report.write("\t\t- HMDB:\t" + str(len(identifier_sbml1[5])) + "\n")
        report.write("\t\t- CAS:\t" + str(len(identifier_sbml1[6])) + "\n")
        report.write("\t\t- Chemspider:\t" + str(len(identifier_sbml1[7])) + "\n")
        report.write("\t\t- Metabolights:\t" + str(len(identifier_sbml1[8])) + "\n")
        
        report.write("\t4. Check SBML2-Model for Ids (KEGG,ChEBI,PubChem, Inchi, BioCyc, HMDB, CAS, Chemspider, Metabolights):\n")
        #Get Identifiers from SBML file 2
        window1 = Tk()
        window1.title("Check SBML file 2 for identifiers!")
        progress_nr = 0
        progress_label = Label(window1, text = 'Check starts!', width = 75)
        progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate", maximum = 100, value = progress_nr)
        progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
        progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
        identifier_sbml2 = check_possible_identifier_to_comp_from_notes_sbml(sbml2_Model, report, window1, progbar, progress_label)
        #write KEGG-IDs found in Species-ID SBML2 to identifier_sbml2
        if len(kegg_id2) != 0 :
            for i in kegg_id2:
                identifier_sbml2[0].update({i:kegg_id2[i]})
        
        #write report for SBML2
        report.write("\t\t- KEGG:\t" + str(len(identifier_sbml2[0])) + "\n")
        report.write("\t\t- ChEBI:\t" + str(len(identifier_sbml2[1])) + "\n")
        report.write("\t\t- PubChem:\t" + str(len(identifier_sbml2[2])) + "\n")
        report.write("\t\t- Inchi:\t" + str(len(identifier_sbml2[3])) + "\n")
        report.write("\t\t- BioCyc:\t" + str(len(identifier_sbml2[4])) + "\n")
        report.write("\t\t- HMDB:\t" + str(len(identifier_sbml2[5])) + "\n")
        report.write("\t\t- CAS:\t" + str(len(identifier_sbml2[6])) + "\n")
        report.write("\t\t- Chemspider:\t" + str(len(identifier_sbml2[7])) + "\n")
        report.write("\t\t- Metabolights:\t" + str(len(identifier_sbml2[8])) + "\n")
        
        structure_frame1 = Frame(unacc_frame)
        structure_frame1.grid(row = 0, column = 0, sticky = W+E+N+S)
        structure_frame2 = Frame(unacc_frame)
        structure_frame2.grid(row = 0, column = 1, sticky = W+E+N+S)
        # GUI information of selected SBML files
        sbml_framelabel = LabelFrame(structure_frame1, text = "Analyzed SBML files:", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
        sbml_framelabel.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
        sbml1_label = Label(sbml_framelabel, text = "SBML file 1:\n\n\t\t" + load_SBML1_label["text"].split("/")[-1] , justify = LEFT, font = label_font, wraplength = label_wraplength)
        sbml1_label.pack(fill = BOTH, padx = 5, pady = 5)
        sbml2_label = Label(sbml_framelabel, text = "SBML file 2:\n\n\t\t" + load_SBML2_label["text"].split("/")[-1] , justify = LEFT, font = label_font, wraplength = label_wraplength)
        sbml2_label.pack(fill = BOTH, padx = 5, pady = 5)
        
        # Analysing results - number of species, database IDs and Inchis found in the files is displayed.
        analysis_framelabel = LabelFrame(structure_frame1, text = "Analysis results", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
        analysis_framelabel.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
        # Headlines
        oversbml1_label = Label(analysis_framelabel, text = "# SBML 1", font = label_font, wraplength = label_wraplength)
        oversbml1_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W)
        
        oversbml2_label = Label(analysis_framelabel, text = "# SBML 2", font = label_font, wraplength = label_wraplength)
        oversbml2_label.grid(row = 0, column = 2, padx = 5, pady = 5, sticky = W)
        
        analysis_result_label = ["Species", "KEGG-IDs", "ChEBI-IDs", "PubChem-IDs", "Inchis", "BioCyc-IDs", "HMDB-IDs", "CAS-IDs", "Chemspider-IDs", "Metabolights-IDs"]
        index_criteria_label = 0
        criteria_row_index = 1
        for entry_label in analysis_result_label:
            criteria_label = Label(analysis_framelabel, text = entry_label+":", font = label_font, wraplength = label_wraplength)
            criteria_label.grid(row = criteria_row_index, column = 0, padx = 5, pady = 5, sticky = W)
            if entry_label == "Species":
                sbml1_label = Label(analysis_framelabel, text = len(sbml1_Model.getListOfSpecies()), font = label_font, wraplength = label_wraplength)
                sbml1_label.grid(row = criteria_row_index, column = 1, padx = 5 , pady = 5, sticky = W+E+N+S)
                sbml2_label = Label(analysis_framelabel, text = len(sbml2_Model.getListOfSpecies()), font = label_font, wraplength = label_wraplength)
                sbml2_label.grid(row = criteria_row_index, column = 2, padx = 5, pady = 5, sticky = W+E+N+S)
            else:
                sbml1_label = Label(analysis_framelabel, text = str(len(identifier_sbml1[index_criteria_label])), font = label_font, wraplength = label_wraplength)
                sbml1_label.grid(row = criteria_row_index, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
                sbml2_label = Label(analysis_framelabel, text = str(len(identifier_sbml2[index_criteria_label])), font = label_font, wraplength = label_wraplength)
                sbml2_label.grid(row = criteria_row_index, column = 2, padx = 5, pady = 5, sticky = W+E+N+S)
                index_criteria_label = index_criteria_label + 1
            criteria_row_index = criteria_row_index + 1
        # Modes can be selected by user
        mode_framelabel = LabelFrame(structure_frame2, text = "Mode-Selection", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
        mode_framelabel.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
        all_mode_yes_checkbuttons = []
        all_intvar_states = []
        
        mode_label = Label(mode_framelabel, text ="Mode", font = label_font, wraplength = label_wraplength)
        mode_label.grid(row = 0, column = 0, padx = 5, pady = 5, sticky = W)
        overyes_label = Label(mode_framelabel, text = "ENABLE?", font = label_font, wraplength = label_wraplength)
        overyes_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W)
        
        mode_select_labels = ["Name-Matching", "KEGG-ID-Matching", "ChEBI-ID-Matching", "PubChem-ID-Matching", "Inchi-Matching", "BioCyc-ID-Matching", "HMDB-ID-Matching", "CAS-ID-Matching", "Chemspider-ID-Matching", "Metabolights-ID-Matching", "KEGG-Synonym-Matching", "ChEBI-Synonym-Matching", "PubChem-Synonym-Matching", "Pathway Tools-Synonym-Matching"]
        index_mode_label = 0
        mode_row_index = 1
        for mode_select_entry in mode_select_labels:
            mode_selection_label = Label(mode_framelabel, text = mode_select_entry, font = label_font, wraplength = label_wraplength)
            mode_selection_label.grid(row = mode_row_index, column = 0, padx = 5, pady = 5, sticky = W)
            var1 = IntVar()
            all_intvar_states.append(var1)
            yes1_button = Checkbutton(mode_framelabel, variable = var1)
            yes1_button.grid(row = mode_row_index, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
            if mode_select_entry == "Name-Matching":
                print "yes"
                if len(sbml1_Model.getListOfSpecies()) == 0 or len(sbml2_Model.getListOfSpecies()) == 0:
                    mode_selection_label.config(fg = "grey")
                    yes1_button.config(state = DISABLED)
            elif mode_select_entry == "KEGG-Synonym-Matching":
                if len(identifier_sbml2[0]) == 0:
                    mode_selection_label.config(fg = "grey")
                    yes1_button.config(state = DISABLED)
            elif mode_select_entry == "ChEBI-Synonym-Matching":
                if len(identifier_sbml2[1]) == 0:
                    mode_selection_label.config(fg = "grey")
                    yes1_button.config(state = DISABLED)
            elif mode_select_entry == "PubChem-Synonym-Matching":
                if len(identifier_sbml2[2]) == 0:
                    mode_selection_label.config(fg = "grey")
                    yes1_button.config(state = DISABLED)
            elif mode_select_entry == "Pathway Tools-Synonym-Matching":
                if len(identifier_sbml2[4]) == 0:
                    mode_selection_label.config(fg = "grey")
                    yes1_button.config(state = DISABLED)
                else:
                    # IF the Pathway-Tools-Synonym-MAtching function is transposed to web query, remove the "else:" 
                    mode_selection_label.config(fg = "grey")
                    yes1_button.config(state = DISABLED)
            else:
                print mode_select_entry
                print index_mode_label
                if len(identifier_sbml1[index_mode_label]) == 0 or len(identifier_sbml2[index_mode_label]) == 0:
                    mode_selection_label.config(fg = "grey")
                    yes1_button.config(state = DISABLED)
                index_mode_label = index_mode_label + 1              
            all_mode_yes_checkbuttons.append(yes1_button)
            mode_row_index = mode_row_index + 1
            
        #GUI De/select all
        def select_all():
            #Function selects all modes
            for but in all_mode_yes_checkbuttons:
                if but["state"] != DISABLED:
                    but.select()
                else:
                    pass


        def deselect_all():
            #Function deselects all modes
            for but in all_mode_yes_checkbuttons:
                if but["state"]!=DISABLED:
                    but.deselect()
                else:
                    pass


        select_all_button = Button(mode_framelabel, text = "Select all", command = select_all, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
        deselect_all_button = Button(mode_framelabel, text = "Deselect all", command = deselect_all, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
        select_all_button.grid(row = 15, column = 0, padx = 5, pady = 5, sticky = W+E+N+S)
        deselect_all_button.grid(row = 15, column=1, padx = 5, pady = 5, sticky = W+E+N+S)
        
        #GUI confirm mode
        def confirm_mode_selection():
            #Function checks input, if no mode selected an error message is shown. If a mode was selected the mode selection is deactivated and the state of the checkbuttons is frozen.No more changes of the mode after confirmation possible.
            a = 0
            for var in all_intvar_states:
                if var.get() == 1:
                    a = a + 1
            
            if a == 0:
                messagebox.showerror("Error", "Your input is empty!")
                start_comparison_button.config(state = DISABLED)
            else:
                select_all_button.config(state = DISABLED)
                deselect_all_button.config(state = DISABLED)
                for butt in all_mode_yes_checkbuttons:
                    butt.config(state = DISABLED)
                
                start_comparison_button.config(state = "normal")
                confirm_mode_selection_button.config(state = DISABLED)

        confirm_mode_selection_button = Button(mode_framelabel, text = "Confirm mode selection!" , command = confirm_mode_selection, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
        confirm_mode_selection_button.grid(row = 16, column = 0, columnspan = 2, padx = 5, pady = 5, sticky = W+E+N+S)
        
        def start_comparison():
            print "yo"
            start_comparison_button.config(state=DISABLED)
                                  
            
            #Get model-data from SBML-file
            sbml1_Model = sbml1_SBML.getModel()
            sbml2_Model = sbml2_SBML.getModel()
            
            # list of sbml1-species-matches
            i_s_match =[]
            # list of sbml2-species-matches
            j_s_match=[]
            # list of sbml1-reaction-matches
            i_r_match =[]
            # list of sbml2-reaction-matches
            j_r_match = []
            #wrong reaction and product matches in reaction comparison
            wrplist = {}
            #wrong reaction matches in reaction comparison
            wrlist = {}
            #wrong product matches in reaction comparison
            wplist = {}
            #reactants and products changed
            crplist = []
            #Dictionary-Entry {sbml1 species-id: sbml2 species ids): Assignment SBML1-Species to SBML2-Species}
            species_ids = {}
            for i in sbml1_Model.getListOfSpecies():
                species_ids.update({i.getId():[]})
                print 1
            
            #Dictionary-Entry {sbml1 reaction-id: sbml2 reaction ids): Assignment SBML1-Reactions to SBML2-Reactions}
            reaction_ids = {}
            for i in sbml1_Model.getListOfReactions():
                reaction_ids.update({i.getId():[]})
                print 2
            
            print "1a"
                        
            
            #Report for species matches
            s_matchlist_path = select_directory_label["text"] + "/02_species_matches.txt"
            species_matches = open(s_matchlist_path,"w")
            species_matches.write(load_SBML1_label["text"].split("/")[-1] + ": " + str(len(sbml1_Model.getListOfSpecies())) + " Species" + "\t" + load_SBML2_label["text"].split("/")[-1] + ": " + str(len(sbml2_Model.getListOfSpecies())) + " Species" + "\n")
            print "1b"
            report.write("Species-Matching:\n")
            # KEGG-ID-Matching
            report.write("\t5. ID-Matching:\n")
            
            job_titles = ["KEGG-ID-Matching!", "ChEBI-ID-Matching!", "PubChem-ID-Matching!", "INCHI-Matching!", "BioCyc-ID-Matching!", "HMDB-ID-Matching!", "CAS-ID-Matching!", "Chemspider-ID-Matching!", "Metabolights-ID-Matching!"]
            identifier_index_for_matching = 0
            for checkbutton_variable_index in range(len(all_intvar_states)):
                if checkbutton_variable_index > 0 and checkbutton_variable_index < 10:
                    report.write("\t\t5." + str(checkbutton_variable_index) + " " + job_titles[checkbutton_variable_index - 1][:-1] + ":\n")
                    if all_intvar_states[checkbutton_variable_index].get() == 1:
                        window1 = Tk()
                        window1.title(job_titles[checkbutton_variable_index - 1])  
                        progress_nr = 0
                        progress_label = Label(window1, text = job_titles[checkbutton_variable_index - 1][:-1] + " starts!", width = 75)
                        progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate", maximum = 100, value = progress_nr)
                        progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
                        progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
                        a = len(i_s_match)
                        b = len(j_s_match)
                        t = db_id_Matching(identifier_sbml1[identifier_index_for_matching], identifier_sbml2[identifier_index_for_matching], i_s_match, j_s_match,species_ids,identifier_index_for_matching, report, window1, progbar, progress_label, species_matches)
                        i_s_match = t[0]
                        j_s_match = t[1]
                        c = len(i_s_match) - a
                        d = len(j_s_match) - b
                        report.write("\t\t\t- SBML1-Matches:\t" + str(c) + "\n")
                        report.write("\t\t\t- SBML2-Matches:\t" + str(d) + "\n")
                        species_ids = t[2]
                        species_count = id_match_counter(species_ids)
                        report.write("\t\t\t- total SBML1-Species with assignment to SBML2-Species: " + str(species_count) + "\n")
                    if all_intvar_states[checkbutton_variable_index].get() == 0 and (len(identifier_sbml1[identifier_index_for_matching]) == 0 or len(identifier_sbml2[identifier_index_for_matching]) == 0):
                        report.write("\t\t\t- no comparison possible!\t" + str(len(identifier_sbml1[identifier_index_for_matching])) + "\tSBML1- KEGG-IDs\t|\t" + str(len(identifier_sbml2[identifier_index_for_matching])) + "\tSBML2- KEGG-IDs\n")
                    if all_intvar_states[checkbutton_variable_index].get() == 0 and (len(identifier_sbml1[identifier_index_for_matching]) != 0 or len(identifier_sbml2[identifier_index_for_matching]) != 0):
                        report.write("\t\t\t- " + job_titles[checkbutton_variable_index - 1][:-2].upper() + " disabled\n")
                    identifier_index_for_matching = identifier_index_for_matching + 1

            # Name-Matching
            report.write("\t6. Name-Matching:\n")
            if all_intvar_states[0].get() == 1:
                window1 = Tk()
                window1.title("Name-Matching!")
                progress_nr = 0
                progress_label = Label(window1, text = 'Name-Matching starts!', width = 75)
                progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate", maximum = 100, value = progress_nr)
                progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
                progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
                a = len(i_s_match)
                b = len(j_s_match)
                t = Name_Matching(sbml1_Model,sbml2_Model,i_s_match,j_s_match,species_ids, report, window1, progbar, progress_label,species_matches)
                i_s_match = t[0]
                j_s_match = t[1]
                c = len(i_s_match) - a
                d = len(j_s_match) - b
                report.write("\t\t- SBML1-Matches:\t" + str(c) + "\n")
                report.write("\t\t- SBML2-Matches:\t" + str(d) + "\n")
                species_ids = t[2]
                species_count = id_match_counter(species_ids)
                report.write("\t\t- total SBML1-Species with assignment to SBML2-Species: " + str(species_count) + "\n")
            if all_intvar_states[0].get() == 0 and (len(sbml1_Model.getListOfSpecies()) == 0 or len(sbml2_Model.getListOfSpecies()) == 0):
                report.write("\t\t\t- no comparison possible!\t" + str(len(identifier_sbml1[8])) + "\tSBML1- Metabolights-IDs\t|\t" + str(len(identifier_sbml2[8])) + "\tSBML2- Metabolights-IDs\n")
            if all_intvar_states[0].get() == 0 and (len(sbml1_Model.getListOfSpecies()) != 0 and len(sbml2_Model.getListOfSpecies()) != 0):
                report.write("\t\t- NAME-MATCHING disabled\n")
            
            # restrict candidates for synonym-matching
            report.write("\t7. Restrict candidates for synonym-matching by formula:\n")
            if all_intvar_states[10].get() == 1 or all_intvar_states[11].get() == 1 or all_intvar_states[12].get() == 1:
                window1 = Tk()
                window1.title("Candidate restriction for Synonym Matching!")
                progress_nr = 0
                progress_label = Label(window1, text = 'Candidate restriction starts!', width = 75)
                progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate", maximum = 100, value = progress_nr)
                progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
                progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
                res_cand = restrict_candidates_for_syn_matching(sbml1_Model,sbml2_Model,i_s_match, report, window1, progbar, progress_label)
                report.write("\t\t- total candidates: " + str(len(res_cand)) + "\n")
            else:
                report.write("\t\tCandidate restriction disabled!\n")
                
            # Synonym-Matching
            syn_match_list = ["KEGG-Symonym-Matching!", "ChEBI-Synonym-Matching!", "PubChem-Synonym-Matching!"]
            report_num = 8
            syn_match_index = 0
            for checkbutton_variable_index in range(len(all_intvar_states)):
                if checkbutton_variable_index > 9 and checkbutton_variable_index < 13:
                    report.write("\t" + str(report_num) + ". " + syn_match_list[syn_match_index][:-1] + "\n")
                    if all_intvar_states[checkbutton_variable_index].get() == 1:
                        window1 = Tk()
                        window1.title(syn_match_list[syn_match_index])
                        progress_nr = 0
                        progress_label = Label(window1, text = syn_match_list[syn_match_index][:-1] + ' starts!', width = 75)
                        progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate",maximum = 100, value = progress_nr)
                        progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
                        progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
                        a = len(i_s_match)
                        b = len(j_s_match)
                        t = synonym_matching(res_cand, i_s_match, j_s_match, species_ids, syn_match_index, report, window1, progbar, progress_label, species_matches)
                        i_s_match = t[0]
                        j_s_match = t[1]
                        c = len(i_s_match) - a
                        d = len(j_s_match) - b
                        report.write("\t\t- SBML1-Matches:\t" + str(c) + "\n")
                        report.write("\t\t- SBML2-Matches:\t" + str(d) + "\n")
                        species_ids = t[2]
                        species_count = id_match_counter(species_ids)
                        report.write("\t\t- total SBML1-Species with assignment to SBML2-Species: " + str(species_count) + "\n")
                    else:
                        report.write("\t\t- " + syn_match_list[syn_match_index][:-2] + " disabled\n")
                     
                    report_num = report_num + 1
                    syn_match_index = syn_match_index + 1
                    
#===============================================================================
#             report.write("\t11. PathwayTools-Synonym-Matching:\n")
#             if all_intvar_states[13].get() == 1:
#                 #######
#                 syn_data = s_func.get_pwytools_syn_data("")
#                 a = len(i_s_match)
#                 b = len(j_s_match)
# ####FUNCTION auf web mdode umstellen
#                 t = s_func.ptoos_syn_matching(sbml1_Model,sbml2_Model,syn_data,i_s_match, j_s_match,species_ids, species_matches, report)
#                 i_s_match = t[0]
#                 j_s_match = t[1]
#                 c = len(i_s_match) - a
#                 d = len(j_s_match) - b
#                 report.write("\t\t- SBML1-Matches:\t" + str(c) + "\n")
#                 report.write("\t\t- SBML2-Matches:\t" + str(d) + "\n")
#                 species_ids = t[2]
#                 species_count = s_func.id_match_counter(species_ids)
#                 report.write("\t\t- total SBML1-Species with assignment to SBML2-Species: " + str(species_count) + "\n")
#             else:
#                 report.write("\t\t- ptools-SYN-MATCHING disabled\n")
#===============================================================================
            report.write("Species-Matching completed!\n")
            ft = datetime.datetime.now() - tt
            report.write("\t- process-time: " + str(ft) + "\n")
            report.write("\t- SBML1-Matches:\t" + str(len(i_s_match)) + "\t/\t" + str(len(sbml1_Model.getListOfSpecies())) + "\tSpecies\n")
            report.write("\t- SBML2-Matches:\t" + str(len(j_s_match)) + "\t/\t" + str(len(sbml2_Model.getListOfSpecies())) + "\tSpecies\n")
            species_count = id_match_counter(species_ids)
            report.write("\t- total SBML1-Species with assignment to SBML2-Species: " + str(species_count) + "\n")
            species_matches.close()
            
            #save id-assignment list (SBML1-Species-Id:[SBML2-Species-Ids])
            species_id_assignment = open(select_directory_label["text"] + "/03_species_id_assignment.txt","w")
            for i in species_ids:
                species_id_assignment.write(i + "\t" + str(species_ids[i]) + "\n")
            
            #reaction match list
            r_matchlist_path = select_directory_label["text"] + "/04_reaction_matches.txt"
            reaction_matches = open(r_matchlist_path,"w")
            reaction_matches.write(load_SBML1_label["text"].split("/")[-1] + ": " + str(len(sbml1_Model.getListOfReactions())) + " Reactions" + "\t" + load_SBML2_label["text"].split("/")[-1] + ": " + str(len(sbml2_Model.getListOfReactions())) + " Reactions" + "\n")
            
            report.write("Reaction-Matching:\n")
            tt2 = datetime.datetime.now()
            
            #reaction some reactants and products wrong
            wrp = open(select_directory_label["text"] + "/05_wrong_reactants_and_products.txt","w")
            #reaction some reactants wrong
            wr = open(select_directory_label["text"] + "/06_wrong_reactants.txt","w")
            #reaction some products wrong
            wp = open(select_directory_label["text"] + "/07_wrong_products.txt","w")
            #reactants and products changed
            crp = open(select_directory_label["text"] + "/08_reactants_products_changed.txt","w")
            
            
            
            #Reaction-Matching
            report.write("\t1.Reaction-Comparison:\n")
            window1 = Tk()
            window1.title("Reaction-Comparison!")
            progress_nr = 0
            progress_label = Label(window1, text = 'Reaction-Comparison starts!', width = 75)
            progbar = ttk.Progressbar(window1, orient = "horizontal", length = 300, mode = "determinate",maximum = 100, value = progress_nr)
            progress_label.grid(row = 0, column = 1, padx = 5, pady = 5, sticky = W+E+N+S) 
            progbar.grid(row = 1, column = 1, padx = 5, pady = 5, sticky = W+E+N+S)
            p = reaction_comparison(sbml1_Model, sbml2_Model, i_r_match, j_r_match, species_ids,wrplist,wplist,wrlist,crplist,reaction_ids, report, progress_label, progbar, window1, reaction_matches, wrp, wp, wr)
            i_r_match = p[0]
            j_r_match = p[1]
            report.write("\t\t- SBML1-Matches:\t" + str(len(i_r_match)) + "\n")
            report.write("\t\t- SBML2-Matches:\t" + str(len(j_r_match)) + "\n")
            reaction_ids = p[2]
            reaction_count = id_match_counter(reaction_ids)
            report.write("\t\t- total SBML1-Reactions with assignment to SBML2-Reactions: " + str(reaction_count) + "\n")
            count_list = p[3]
            report.write("\t\t- Reaction-Matches: " + "\n")
            report.write("\t\t\t- 1 = Reaction_SBML1(Reactants&Products) == Reaction_SBML2(Reactants&Products)" + "\n")
            report.write("\t\t\t- 2 = Reaction_SBML1(Reversebility) == Reaction_SBML2(Reversebility)" + "\n")
            report.write("\t\t\t- 3 = Reaction_SBML1(Stoichiometry) == Reaction_SBML2(Stoichiometry)" + "\n")
            report.write("\t\t- total 123_Reaction_Match:\t" + str(count_list[0]) + "\n")
            report.write("\t\t- total 12_Reaction_Match:\t" + str(count_list[1]) + "\n")
            report.write("\t\t- total 13_Reaction_Match:\t" + str(count_list[2]) + "\n")
            report.write("\t\t- total 1_Reaction_Match:\t" + str(count_list[3]) + "\n")
            count_list_no_match = p[4]
            wrplist = p[5]
            wplist = p[6]
            wrlist = p[7]
            crplist = p[8]
            print str(len(wrplist))
            print str(len(wplist))
            print str(len(wrlist))
            for ir in i_r_match:
                if ir in wrplist:
                    wrplist.pop(ir)
                if ir in wrlist:
                    wrlist.pop(ir)
                if ir in wplist:
                    wplist.pop(ir)

            for iwrp in wrplist:
                for iwrp1 in wrplist[iwrp]:
                    wrp.write(iwrp1)
            
            for iwr in wrlist:
                for iwr1 in wrlist[iwr]:
                    wr.write(iwr1)
            
            for iwp in wplist:
                for iwp1 in wplist[iwp]:
                    wp.write(iwp1)
            
            for icrp in crplist:
                crp.write(icrp)
            
            report.write("\t\t- Reactions reactants and products changed:\n")
            report.write("\t\t- total:\t" + str(len(crplist)) + "\n")
            report.write("\t\t- Reactions with some(i) products and (j) reactants missing, i between 0 < i < length(productlist),    j beetween 0 < j < length(reactantlist):\n")
            report.write("\t\t- total:\t" + str(len(wrplist)) + "\n")
            report.write("\t\t- Reactions with some(i) products missing, i between 0 < i < length(productlist):\n")
            report.write("\t\t- total:\t" + str(len(wplist)) + "\n")
            report.write("\t\t- Reactions with some(j) reactants missing,j between 0 < j < length(reactantlist):\n")
            report.write("\t\t- total:\t" + str(len(wrlist)) + "\n")

            report.write("Reaction-Matching completed!\n")
            ft2 = datetime.datetime.now() - tt2
            report.write("\t- process-time: " + str(ft2) + "\n")
            reaction_matches.close()

            wrp.close()
            wr.close()
            wp.close()
            
            #for i in reaction_ids:
                #if len(reaction_ids[i])>1:
                    #print reaction_ids[i]
            #save id-assignment list (SBML1-Species-Id:[SBML2-Species-Ids])
            reaction_id_assignment = open(select_directory_label["text"] + "/09_reaction_id_assignment.txt","w")
            for i in reaction_ids:
                reaction_id_assignment.write(i + "\t" + str(reaction_ids[i]) + "\n")
            
            report.write("No-Matches:\n")
            #species sbml1 no match
            s_sbml1_nomatch_path = select_directory_label["text"] + "/10_SBML1_species_no_match.txt"
            species_sbml1_no_match = open(s_sbml1_nomatch_path,"w")
            species_sbml1_no_match.write(load_SBML1_label["text"].split("/")[-1] + ": " + str(len(sbml1_Model.getListOfSpecies())) + " Species" + "\t" + load_SBML2_label["text"].split("/")[-1] + ": " + str(len(sbml2_Model.getListOfSpecies())) + " Species" + "\n\n")
            
            print "SBML1 Species No-Match: " + "\n"
            count_sbml1_species_no_match=0
            for i in sbml1_Model.getListOfSpecies():
                if i.getId() not in i_s_match:
                    #print i.getName()
                    count_sbml1_species_no_match=count_sbml1_species_no_match+1
                    species_sbml1_no_match.write(i.getId() + "\t" + i.getName() + "\n")
            
            species_sbml1_no_match.close()
            report.write("\t- SBML1-Species:\t" + str(count_sbml1_species_no_match) + "\tof\t" + str(len(sbml1_Model.getListOfSpecies())) + "\tSpecies" + "\n")
            #sbml2 species no match
            s_sbml2_nomatch_path = select_directory_label["text"] + "/11_SBML2_species_no_match.txt"
            species_sbml2_no_match = open(s_sbml2_nomatch_path,"w")
            species_sbml2_no_match.write(load_SBML1_label["text"].split("/")[-1] + ": " + str(len(sbml1_Model.getListOfSpecies())) + " Species" + "\t" + load_SBML2_label["text"].split("/")[-1] + ": " + str(len(sbml2_Model.getListOfSpecies())) + " Species" + "\n\n")

            print "SBML2 Species No-Match: " + "\n"
            count_sbml2_species_no_match = 0
            for j in sbml2_Model.getListOfSpecies():
                if j.getId() not in j_s_match:
                    #print j.getName()
                    count_sbml2_species_no_match = count_sbml2_species_no_match + 1
                    species_sbml2_no_match.write(j.getId() + "\t" + j.getName() + "\n")
            
            species_sbml2_no_match.close()
            report.write("\t- SBML2-Species:\t" + str(count_sbml2_species_no_match) + "\tof\t" + str(len(sbml2_Model.getListOfSpecies())) + "\tSpecies" + "\n")
            #reaction sbml1 no match list
            r_sbml1_nomatch_path = select_directory_label["text"] + "/12_SBML1_reactions_no_match.txt"
            reactions_sbml1_no_match = open(r_sbml1_nomatch_path,"w")
            reactions_sbml1_no_match.write(load_SBML1_label["text"].split("/")[-1] + ": " + str(len(sbml1_Model.getListOfReactions())) + " Reactions" + "\t" + load_SBML2_label["text"].split("/")[-1] + ": " + str(len(sbml2_Model.getListOfReactions())) + " Reactions"+"\n\n")
            
            #sbml1 reactions no match
            print "SBML1 Reactions No-Match: " + "\n"
            count_sbml1_reactions_no_match=0
            for i in sbml1_Model.getListOfReactions():
                if i.getId() not in i_r_match:
                    #print i.getName()
                    count_sbml1_reactions_no_match = count_sbml1_reactions_no_match+1
                    reactions_sbml1_no_match.write(i.getId() + "\t" + i.getName() + "\n")
            
            reactions_sbml1_no_match.close()
            report.write("\t- SBML1-Reactions:\t" + str(count_sbml1_reactions_no_match) + "\tof\t" + str(len(sbml1_Model.getListOfReactions())) + "\tReactions" + "\n")
            #reaction sbml2 no match list
            r_sbml2_nomatch_path = select_directory_label["text"] + "/13_SBML2_reactions_no_match.txt"
            reactions_sbml2_no_match = open(r_sbml2_nomatch_path,"w")
            reactions_sbml2_no_match.write(load_SBML1_label["text"].split("/")[-1] + ": " + str(len(sbml1_Model.getListOfReactions())) + " Reactions" + "\t" + load_SBML2_label["text"].split("/")[-1] + ": " + str(len(sbml2_Model.getListOfReactions())) + " Reactions" + "\n\n")

            #sbml2 reactions no match
            print "SBML2 Reactions No-Match: " + "\n"
            count_sbml2_reactions_no_match = 0
            for j in sbml2_Model.getListOfReactions():
                if j.getId() not in j_r_match:
                    #print j.getName()
                    count_sbml2_reactions_no_match = count_sbml2_reactions_no_match + 1
                    reactions_sbml2_no_match.write(j.getId() + "\t" + j.getName() + "\n")
            
            reactions_sbml2_no_match.close()
            report.write("\t- SBML2-Reactions:\t" + str(count_sbml2_reactions_no_match) + "\tof\t" + str(len(sbml2_Model.getListOfReactions())) + "\tReactions" + "\n")
            report.write("\nModel-Similarity:\n")
            jaccardcoef = jaccardc(count_sbml2_reactions_no_match,count_sbml1_reactions_no_match,reaction_count)
            jaccarddist = jaccardd(count_sbml2_reactions_no_match,count_sbml1_reactions_no_match,reaction_count)
            report.write("Jaccard similarity coefficient:\t" + str(jaccardcoef) + "\n")
            report.write("Jaccard distanz:\t" + str(jaccarddist) + "\n")
            ft3 = datetime.datetime.now() - time1
            report.write("-total process-time:\t" + str(ft3) + "\n")
            report.close() 
            #print report
            #reporting = ""
            #print reporting
            reporting = open(report_path,"r").read()
            # GUI New Tab "General Informations"
            tab3 = Frame(noteb)
            tab3.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
            noteb.add(tab3, text = "Comparison Results")
            noteb.select(tab3)
            
            # GUI add objects for a scrollable tab
            mycanvas1 = Canvas(tab3)
            mycanvas1.pack(side = "left", fill = BOTH, expand = 1)
 
            unacc_frame1 = Frame(mycanvas1)
            unacc_frame1.pack(fill = BOTH, expand = 1, padx = 5, pady = 5)
            mycanvas1.create_window(0, 0, window = unacc_frame1, anchor = "nw")
             
            matching_report_framelabel = LabelFrame(unacc_frame1, text = "Matching-Report", relief = framelabel_relief, font = framelabel_font, borderwidth = framelabel_borderwidth)
            matching_report_framelabel.pack(fill = X, padx = 5, pady = 5)
             
            matching_report_label = Label(matching_report_framelabel, text = reporting, justify = LEFT, font = label_font, wraplength = label_wraplength) 
            matching_report_label.pack(fill = X, padx = 5, pady = 5)
             
            def open_reportdir():
                if operating_system == "Windows":
                    startfile(select_directory_label["text"])
                else:
                    subprocess.call(['xdg-open', select_directory_label["text"]])
                 
            open_report_directory_button = Button(unacc_frame1, text = "Open Report directory!", command = open_reportdir, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
            open_report_directory_button.pack(fill = X, padx = 5, pady = 5)
            
            quit_button2 = Button(unacc_frame1, text = "Quit", command = close_program, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
            quit_button2.pack(fill = X, padx = 5, pady = 5)
            mycanvas1.update()
             
            if unacc_frame1.winfo_height() > (main_window.winfo_height()-50):
                myscrollbar1 = Scrollbar(tab3, command = mycanvas1.yview)
                mycanvas1.config(yscrollcommand = myscrollbar1.set, scrollregion = mycanvas1.bbox("all"))
                myscrollbar1.pack(side=RIGHT, fill=Y)
                     
                def mousewheel1(event):
                    if operating_system == "Windows":
                        mycanvas1.yview_scroll(int(-1*(event.delta/120)), "units")
                    elif operating_system == "Linux":
                        if event.num == 4:
                            mycanvas1.yview_scroll(-1, "units")
                        if event.num == 5:
                            mycanvas1.yview_scroll(1, "units")
                    elif operating_system == "Darwin":
                        mycanvas1.yview_scroll((event.delta/120), "units")
 
 
                if operating_system == "Linux":
                    mycanvas1.bind_all("<Button-4>", mousewheel1)
                    mycanvas1.bind_all("<Button-5>", mousewheel1)
                else:    
                    mycanvas1.bind_all("<MouseWheel>", mousewheel1)

        
        start_comparison_button = Button(unacc_frame, text = "Start SBML file comparison!",  state = DISABLED, command = start_comparison, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
        start_comparison_button.grid(row = 1, column = 0, columnspan = 2, padx = 5, pady = 5, sticky = W+E+N+S)
        
        # If the content of the unacc_frame is to big for the window, than add a scrollbar
        mycanvas.update()
        main_window.update()
        if unacc_frame.winfo_height() > (main_window.winfo_height()-50):
            myscrollbar = Scrollbar(tab2, command = mycanvas.yview)
            mycanvas.config(yscrollcommand = myscrollbar.set, scrollregion = mycanvas.bbox("all"))
            myscrollbar.pack(side=RIGHT, fill=Y, padx = 0)
                
            def mousewheel(event):
                if operating_system == "Windows":
                    mycanvas.yview_scroll(int(-1*(event.delta/120)), "units")
                elif operating_system == "Linux":
                    if event.num == 4:
                        mycanvas.yview_scroll(-1, "units")
                    if event.num == 5:
                        mycanvas.yview_scroll(1, "units")
                elif operating_system == "Darwin":
                    mycanvas.yview_scroll((event.delta/120), "units")


            if operating_system == "Linux":
                mycanvas.bind_all("<Button-4>", mousewheel)
                mycanvas.bind_all("<Button-5>", mousewheel)
            else:    
                mycanvas.bind_all("<MouseWheel>", mousewheel)
                
        pass

    analyze_sbml_files_button = Button(tab1, text = "Analyze SBML files!", state = DISABLED, command = analyze_SBMl_files, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
    analyze_sbml_files_button.pack(fill = BOTH, expand = 1, padx = 5, pady = 10)

start_button = Button(main_frame, text = "Start", command = choose_files, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
start_button.pack(side = LEFT, expand = 1, fill = BOTH, padx =  5, pady = 10)

def close_program():
    main_window.destroy()


quit_button = Button(main_frame, text = "Quit", command = close_program, font = button_font, cursor = button_cursor, overrelief = button_overrelief, bg = button_bg, borderwidth = button_borderwidth)
quit_button.pack(side = RIGHT, expand = 1, fill = BOTH, padx =  5, pady = 10)
   
main_window.mainloop()
