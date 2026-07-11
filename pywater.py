#!/usr/bin/env python

"""
PyWATER finds conserved water molecules in X-ray protein structure (pdb).
version 1.0

Important or conserved waters are the waters molecules which are present in most or all available pdb structures when superimposed.

Copyright 2013 Hitesh Patel and B. Gruening

Run it by following command:

GUI of the tool can be run in PyMOL as a plugin.
Steps to install in PyMOL:
    Run PyMOL as administrator.
    Install the PyWATER plugin in pymol by following the path: Plugins -> Manage Plugins -> install
    Restart the PyMol

After installation as plugin. It can be run from command in pymol 
    pywater [PDB id , Chain id [, sequence identity cutoff [, resolution cutoff [, refinement assessing method [, user defined proteins list [, linkage method [, inconsistency coefficient threshold [, degree of conservation]]]]]]]] 

    API extension:

    from pymol import cmd
    
    cmd.pywater(PDB id , Chain id [, sequence identity cutoff [, resolution cutoff [, refinement assessing method [, user defined proteins list [, linkage method [, inconsistency coefficient threshold [, degree of conservation]]]]]]])

    PDB id
            string: The PDB id of the protein for which you like to find conserved waters. {default: None}

    Chain id
            string: The chain identifier of the protein for which you like to find conserved waters in above mentioned PDB. {default: None}

    sequence identity cutoff
            string: '30', '40', '50', '70', '90', '95'or '100'. All the protein structures, clustered by BlastClust, having sequence identity more than given cutoff will be superimposed to find the conserved water molecules in query protein chain. {default: '95'} 

    resolution cutoff
            float: All the protein structures to be superimposed will be filtered first according to the structure resolution cutoff. Only structures with better resolution than given cutoff will be used further. {default: 2.0}

    refinement assessing method
            string: Choose either 'Mobility' or 'Normalized B-factor' or 'No refinement' as criteria to assess the refinement quality of crystal structure. Program will filter out the water molecules with bad refinement quality. {default: 'Mobility'}

    user defined proteins list
            string: Give a custom list of protein structures to superimpose. Specifying this list will disable 'sequence identity' and 'resolution cutoff' parameters. {default: disabled}

    linkage method
            string: Linkage method for hierarchical clustering. Choose one from single, complete, average. {default: complete}

    inconsistency coefficient threshold
            float: Any two clusters of water molecules will not be closer than given inconsistency coefficient threshold. Value ranges from 0 to 2.8. {default: 2.4} 

    degree of conservation
            float: Water molecules will be considered CONSERVED if their probability of being conserved is above given cutoff. Value ranges from 0 to 1. {default: 0.7} 

"""

import os
import glob
import shutil
import re
import collections
import tempfile
import json
import threading
import urllib.request as urllib
from concurrent.futures import ThreadPoolExecutor

from pymol.Qt import QtWidgets, QtCore

import pymol.cmd as cmd
import logging

try:
    import numpy as np
except ImportError:
    raise ImportError('Numpy not found. Install NumPy in PyMOL\'s Python environment.')

try:
    import scipy.cluster.hierarchy as hcluster
except ImportError:
    raise ImportError('Scipy not found. Install SciPy in PyMOL\'s Python environment.')


# setup output directory

home_dir = os.path.expanduser("~")

outdir = os.path.join( home_dir, 'PyWATER_outdir' )
if not os.path.exists(outdir):
    os.mkdir(outdir)


# setup logging

logger = logging.getLogger('PyWATER')
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler( os.path.join( outdir, 'pywater.log' ) )
fh.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)


RCSB_DATA_URL = "https://data.rcsb.org"
RCSB_GRAPHQL_URL = "%s/graphql" % RCSB_DATA_URL
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/%s.pdb"
HTTP_TIMEOUT = 20
JSON_HEADERS = {"Content-Type": "application/json"}
DOWNLOAD_HEADERS = {"User-Agent": "Mozilla/5.0"}


def _decode_response(response):
    data = response.read()
    if isinstance(data, bytes):
        return data.decode("utf-8")
    return data


def _get_json(url, timeout=HTTP_TIMEOUT):
    with urllib.urlopen(url, timeout=timeout) as response:
        return json.loads(_decode_response(response))


def _post_json(url, payload, timeout=HTTP_TIMEOUT):
    request = urllib.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers=JSON_HEADERS,
    )
    with urllib.urlopen(request, timeout=timeout) as response:
        return json.loads(_decode_response(response))


def _graphql(query, timeout=HTTP_TIMEOUT):
    return _post_json(RCSB_GRAPHQL_URL, {"query": query}, timeout=timeout)


# PyMOL plugin registration (modern Qt plugin API).
def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PyWATER', main)


# ---------------------------------------------------------------------------
# User messaging
#
# PyWATER historically used tkinter.messagebox for help and error dialogs.
# On macOS, Tk cannot run inside PyMOL's Qt GUI (it aborts the process during
# Cocoa color allocation), so messaging now goes through Qt and the console.
# ---------------------------------------------------------------------------

def _qt_info(title, message):
    """Native Qt information dialog, used by the GUI help buttons."""
    logger.info('%s: %s' % (title, message))
    parent = QtWidgets.QApplication.activeWindow()
    QtWidgets.QMessageBox.information(parent, title, message)


class _PopupFreeMessageBox(object):
    """
    Drop-in stand-in for tkinter.messagebox that never opens a modal dialog.
    Validation warnings are written to the PyMOL console / log instead, so the
    command line stays popup-free and error paths cannot crash PyMOL.
    """
    def showinfo(self, title='', message=''):
        logger.info('%s: %s' % (title, message))


tkMessageBox = _PopupFreeMessageBox()


# Display help messages
def pdb_id_help():
    _qt_info('PDB Identifier',
        "The PDB id of the protein for which you like to find conserved waters, e.g. 4lyw")

def chain_help():
    _qt_info('Chain Identifier',
        "The chain identifier of the protein for which you like to find conserved waters in above mentioned PDB, e.g. A.")

def seq_id_help():
    _qt_info('Sequence identity cutoff',
        """All the protein structures, clustered by sequence identity, having sequence identity more than given cutoff will be superimposed to find the conserved water molecules in query protein chain.
Minimum suggested Sequence identity cutoff is 95.""")

def resolution_help():
    _qt_info('Structure resolution cutoff',
        """All the protein structures to be superimposed will be filtered first according to the structure resolution cutoff. Only structures with better resolution than given cutoff will be used further.
Maximum suggested structure resolution cutoff is 2.5 A.""")

def refinement_quality_help():
    _qt_info('Filter by refinement quality',
        "Choose either Mobility or Normalized B-factor as criteria to assess the refinement quality of crystal structure. Program will filter out the water molecules with bad refinement quality.")

def user_defined_lists_help():
    _qt_info('User defined pdb-chains lists',
        """The user defined list of pdbchains to superimpose to find conserved waters in query protein structure.
Enter the pdb chains list in format: xxxx_x,yyyy_y,zzzz_z""")

def clustering_method_help():
    _qt_info('Clustering linkage method',
        """Choose any of the linkage method for hierarchical clustering. Default method is 'complete'.""")

def inconsistency_coefficient_help():
    _qt_info('Inconsistency coefficient threshold',
        """Any two clusters of water molecules will not be closer than given inconsistency coefficient threshold. The less threshold ensures each water molecule in a cluster from different structures.
Maximum suggested inconsistency coefficient threshold is 2.4 A.""")

def prob_help():
    _qt_info('Degree of conservation',
        """Water molecules will be considered CONSERVED if their probability of being conserved is above given cutoff.
Value ranges from 0 to 1.
Minimum suggested value is 0.5""")

def save_sup_files_help():
    _qt_info('Save superimposed files', """Save superimposed intermediate files.""")

def local_files_help():
    _qt_info('Use local files (skip RCSB)',
        """Find conserved waters across your own local .pdb files instead of RCSB homologs (e.g. unpublished structures on your workstation).
Point 'Local files folder' at a folder containing at least two .pdb files (all are compared), set 'Chain id' to the chain to analyse (applied to every file), and enter the 'Reference file' whose conserved waters are reported.
Sequence identity, resolution and maximum-structures settings are ignored in this mode. Currently supports .pdb files with single-character chain ids, up to 100 files.""")

# Display imput parameters

def displayInputs( selectedStruturePDB, selectedStrutureChain,
        seq_id,resolution, refinement, user_def_list,
        clustering_method, inconsistency_coefficient, prob):
    logger.info( 'Input values' )
    logger.info( 'PDB id : %s' % selectedStruturePDB)
    logger.info( 'Chain id : %s' % selectedStrutureChain)
    logger.info( 'Seqence identity cutoff : %s' % seq_id)
    logger.info( 'Structure resolution cutoff : %s' % resolution)
    logger.info( 'X-ray structure refinement assessing method : %s' % refinement)
    logger.info( 'User defined protein-chains list : %s' % user_def_list)
    logger.info( 'Hierarchical clustering linkage method : %s' % clustering_method)
    logger.info( 'Inconsistency coefficient threshold : %s' % inconsistency_coefficient)
    logger.info( 'probability cutoff : %s' % prob)

# display PyMOL session with identified conserved waters, showing H-bonds with other conserved waters, ligands or protein.
def displayInPyMOL(outdir, selectedPDBChain, atomNumbersProbDic):
    pdbCWMs = os.path.join(outdir, '%s_withConservedWaters.pdb' % selectedPDBChain)
    pdb = os.path.join(outdir, '%s.pdb' % selectedPDBChain)
    h_bond_dist = 4.0

    queryProteinCWMs = '%s_withConservedWaters' % selectedPDBChain

    cmd.load(pdbCWMs)
    cmd.orient(queryProteinCWMs)
    cmd.h_add(queryProteinCWMs)

    cmd.select('cwm_protein','polymer and %s' % queryProteinCWMs)
    cmd.select('cwm_waters','resn hoh and %s' % queryProteinCWMs)
    cmd.select('cwm_ligand','organic and %s' % queryProteinCWMs)

    cmd.select('don', '(elem n,o and (neighbor hydro)) and %s' % queryProteinCWMs)
    cmd.select('acc', '(elem o or (elem n and not (neighbor hydro))) and %s' % queryProteinCWMs)
    # h bonds between protein and conserved waters
    cmd.distance ('PW_HBA', '(cwm_protein and acc)','(cwm_waters and don)', h_bond_dist)
    cmd.distance ('PW_HBD', '(cwm_protein and don)','(cwm_waters and acc)', h_bond_dist)
    # h bonds between ligands and conserved waters
    cmd.distance ('LW_HBA', '(cwm_ligand and acc)','(cwm_waters and don)', h_bond_dist)
    cmd.distance ('LW_HBD', '(cwm_ligand and don)','(cwm_waters and acc)', h_bond_dist)
    # h bonds in between conserved waters
    cmd.distance ('HW_HBA', '(cwm_waters and acc)','(cwm_waters and don)', h_bond_dist)
    cmd.distance ('HW_HBD', '(cwm_waters and don)','(cwm_waters and acc)', h_bond_dist)

    cmd.delete('don')
    cmd.delete('acc')

    cmd.set('dash_color','yellow')
    cmd.set('dash_gap',0.3)
    cmd.set('dash_length',0.2)
    cmd.set('dash_round_ends','on')
    cmd.set('dash_width',3)

    # color and display
    cmd.util.cbam('cwm_ligand')
    cmd.show_as('sticks','cwm_ligand')

    MinDoc = min(atomNumbersProbDic.values())
    MaxDoc = max(atomNumbersProbDic.values())
    cmd.create ('conserved_waters','cwm_waters')
    for key, value in atomNumbersProbDic.items():
        cmd.alter('/conserved_waters//A/HOH`%s/O' % key, 'b=%s' % value)

    cmd.spectrum('b', 'red_blue', 'conserved_waters',minimum=MinDoc, maximum=MaxDoc)
    cmd.ramp_new('DOC', 'conserved_waters', range = [MinDoc,MaxDoc], color = '[red,blue]')
    cmd.set('sphere_scale',0.40,'conserved_waters')
    cmd.show_as('spheres','conserved_waters')

    cmd.hide('labels','*_HB*')
    cmd.remove('(hydro) and conserved_waters')
    cmd.remove('(hydro) and %s' % queryProteinCWMs)

    cmd.util.cbac('cwm_protein')
    cmd.set('transparency', 0.2)
    cmd.show('surface', 'cwm_protein')
    cmd.set('surface_color', 'gray', 'cwm_protein')

    cmd.load(pdb)
    cmd.create('all_waters', 'resn hoh and %s' % selectedPDBChain)
    cmd.color('red','ss h and %s' % selectedPDBChain)
    cmd.color('yellow','ss s and %s' % selectedPDBChain)
    cmd.color('green','ss l+ and %s' % selectedPDBChain)
    cmd.show_as('cartoon', selectedPDBChain)
    cmd.set('ray_shadows', 0)


def okMobility( pdbFile, mobilityCutoff = 2.0 ):
    """
    Check if the mobility of water molecules is acceptable. 
    Water oxygen atoms with mobility >= mobilityCutoff(default=2.0) are 
    removed from the PDB. If more than 50 % of water oxygen atoms 
    are removed than whole PDB is discarded 
    """
    normBfactors = []
    OccupancyAndBfactor = []
    with open(pdbFile) as f:
        for line in f:
            line = line.strip()
            if line.startswith('HETATM'):
                OccupancyAndBfactor.append([float(line[54:60]),float(line[60:66])])
    occupancy = [a[0] for a in OccupancyAndBfactor]
    Bfactors = [a[1] for a in OccupancyAndBfactor]
    avgB = np.mean(Bfactors)
    avgO = np.mean(occupancy)
    with open(pdbFile) as f:
        pdbFileLines = f.readlines()
    nWaters = len(pdbFileLines)-1
    logger.debug( 'Number of water molecules: %s' %nWaters )
    count = 0
    for line in reversed(pdbFileLines):
        if line.startswith('HETATM'):
            m = ((float(line[60:66])/avgB)/(float(line[54:60])/avgO))
            if m >= mobilityCutoff:
                count+=1
                pdbFileLines.remove(line)
    logger.debug( 'Water oxygen atoms with a higher mobility than %s: %s: ' % (mobilityCutoff, count))
    if count > (nWaters/2):
        considerPDB = False
    elif count > 0:
            with open(pdbFile,'w') as outfile:
                outfile.write("".join(pdbFileLines))
            considerPDB = True
    else:
        considerPDB = True
    if considerPDB:
        logger.debug( '%s is included in the prediction.' % (pdbFile))
    else:
        logger.debug( '%s is excluded from the prediction.' % (pdbFile))
    return considerPDB

def okBfactor( pdbFile, normBCutoff = 1.0 ):
    """
    Check if the normalized B-factor of water molecules is acceptable.
    Water oxygen atoms with normalized B-factor >= normBCutoff(default=1.0) 
    are removed from the PDB. If more than 50 % of water 
    oxygen atoms are removed than whole PDB is discarded 
    """
    Bfactors = []
    normBfactors = []
    with open(pdbFile) as f:
        for line in f:
            line = line.strip()
            if line.startswith('HETATM'):
                Bfactors.append( float(line[60:66]) )
    avg = np.mean( Bfactors )
    stddev = np.sqrt( np.var(Bfactors) )
    with open(pdbFile) as f:
        pdbFileLines = f.readlines()
    nWaters = len(pdbFileLines) - 1
    logger.debug( 'Number of water molecules is : %s' %nWaters )
    count = 0
    for line in reversed(pdbFileLines):
        if line.startswith('HETATM'):
            normB = ( (float(line[60:66]) - avg) / stddev )
            if normB >= normBCutoff:
                count+=1
                pdbFileLines.remove(line)
    logger.debug( 'water oxygen atoms having higher normalized B factor are %s: '% count)
    if count > (nWaters/2):
        considerPDB = False
    elif count > 0:
            with open(pdbFile,'w') as outfile:
                outfile.write("".join(pdbFileLines))
            considerPDB = True
    else:
        considerPDB = True
    logger.debug( '%s is considered : %s ' % (pdbFile, considerPDB))
    return considerPDB


class Protein():
    def __init__(self, pdb_id, chain=False):
        self.pdb_id = pdb_id.lower()
        self.chain = chain.upper()
        self.pdb_id_folder = self.pdb_id[1:3]
        self.pdb_filename = self.pdb_id + '.pdb'
        self.water_coordinates = list()
        self.water_ids = list()
        self.waterIDCoordinates = {}

    def __repr__(self):
        return "%s_%s" % (self.pdb_id, self.chain)

    def calculate_water_coordinates(self, tmp_dir = False):
        path = os.path.join( tmp_dir, 'cwm_%s_Water.pdb' % self.__repr__() )
        logger.debug( 'Creating water coordinates of cwm_%s_Water.pdb.' % self.__repr__())
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('HETATM'):
                    # water oxygen atom number
                    key = "%s_%s" % (self.__repr__(), int(line[22:30]))
                    coordinates = [ line[30:38], line[38:46], line[46:54] ]
                    self.waterIDCoordinates[key] = coordinates
                    self.water_coordinates.append( coordinates )
                    self.water_ids.append( key )
        return self.waterIDCoordinates


class ProteinsList():
    def __init__(self, ProteinName):
        self.ProteinName = ProteinName
        self.proteins = list()
        self.selectedPDBChain = ''
        self.probability = 0.7
        self.inconsistency_coefficient = 2.0
        self.refinement = ''

    def add_protein(self, protein):
        self.proteins.add(protein)

    def add_protein_from_string(self, s):
        """
            Accepts a '1234:B' formatted string and creates a protein object out of it
            if no double point is found, it is assumed that the whole string is the pdb id
        """
        s = s.strip()
        if s.find(':') > -1:
            pdb_id, chain = s.split(':')
        else:
            pdb_id = s
            chain = False
        if len(pdb_id) != 4:
            raise ValueError('%s is not of length 4' % pdb_id)
        p = Protein(pdb_id, chain)
        self.proteins.append(p)

    def pop(self, index):
        return self.proteins.pop(index)

    def remove(self, protein):
        return self.proteins.remove(protein)

    def __iter__(self):
        for protein in self.proteins:
            yield protein

    def __len__(self):
        return len(self.proteins)

    def __getitem__( self, key ) :
        if isinstance( key, slice ) :
            # Get the start, stop, and step from the slice
            return [self.proteins[ii] for ii in range(*key.indices(len(self)))]
        elif isinstance( key, int ):
            # Handle negative indices
            if key < 0 :
                key += len( self.proteins )
            if key >= len( self.proteins ) :
                raise IndexError("The index (%d) is out of range." % key)
            return self.proteins[key] #Get the data from elsewhere
        else:
            raise TypeError("Invalid argument type.")


def makePDBwithConservedWaters(ProteinsList, temp_dir, outdir, save_sup_files):
    logger.info( 'Minimum desired degree of conservation is : %s' % ProteinsList.probability )
    cmd.delete('cwm_*')
    logger.info( 'Loading all pdb chains ...' )
    for protein in ProteinsList:
        cmd.load(os.path.join(temp_dir, protein.pdb_filename),'cwm_%s' % protein.pdb_id)
        cmd.remove('(hydro) and cwm_%s' % protein.pdb_id)
        cmd.select('dods','resn dod')
        cmd.alter('dods', 'resn="HOH"')
        cmd.create('cwm_%s' % protein, 'cwm_%s & chain %s' % (protein.pdb_id, protein.chain))
        cmd.delete( 'cwm_%s' % protein.pdb_id )

    logger.info( 'Superimposing all pdb chains ...' )
    for protein in ProteinsList[1:]:
        logger.info( 'Superimposing %s' % protein )
        cmd.super('cwm_%s////CA' % protein, 'cwm_%s////CA' % ProteinsList[0])
        cmd.orient( 'cwm_%s' % ProteinsList[0] )

    logger.info( 'Creating new, water only, pymol objects for each pdb chain ...' )
    for protein in ProteinsList:
        cmd.create('cwm_%s_Water' % protein, 'cwm_%s & resname hoh' % protein)

    logger.info( 'Storing water molecules and proteins in separate pdb files for each pdb chain ...' )
    for protein in ProteinsList:
        cmd.save(os.path.join(temp_dir, 'cwm_%s.pdb' % protein), 'cwm_%s' % protein)
        cmd.save(os.path.join(temp_dir, 'cwm_%s_Water.pdb' % protein), 'cwm_%s_Water' % protein)

    cmd.delete('cwm_*')

    ### filter ProteinsList by mobility or normalized B factor cutoff
    logger.debug( 'Protein chains list is %s proteins long.' % len(ProteinsList.proteins) )
    if ProteinsList.refinement != 'No refinement':
        length = len(ProteinsList.proteins)
        if ProteinsList.refinement == 'Mobility':
            logger.info( 'Filtering water oxygen atoms by mobility ...' )
            for protein in reversed(ProteinsList.proteins):
                if str(protein) != str(ProteinsList.selectedPDBChain):
                    if not okMobility(os.path.join(temp_dir, 'cwm_%s_Water.pdb' % protein)):
                        ProteinsList.proteins.remove(protein)

        if ProteinsList.refinement == 'Normalized B-factor': 
            logger.info( 'Filtering water oxygen atoms by Normalized B-factor' )
            for protein in reversed(ProteinsList.proteins):
                if str(protein) != str(ProteinsList.selectedPDBChain):
                    if not okBfactor(os.path.join(temp_dir, 'cwm_%s_Water.pdb' % protein)):
                        ProteinsList.proteins.remove(protein)

        retained = len(ProteinsList.proteins)
        logger.info( '%s refinement filtering: %i of %i structures retained (%i discarded for poor water refinement).' % (ProteinsList.refinement, retained, length, length - retained) )
        logger.debug( 'filtered proteins chains list is %s proteins long :' % len(ProteinsList.proteins) )

    """ 
        Filtered ProteinsList
    """

    selectedPDBChain = str(ProteinsList.selectedPDBChain)
    if not os.path.exists(os.path.join(outdir,selectedPDBChain)):
        os.mkdir(os.path.join(outdir,selectedPDBChain))

    if save_sup_files:
        for file in glob.glob(os.path.join(temp_dir, 'cwm_????_?.pdb')):
            shutil.copy(file, os.path.join(outdir,selectedPDBChain))

    # Only if ProteinsList has more than one protein
    if len(ProteinsList.proteins) > 1:
        water_coordinates = list()
        waterIDCoordinates = {}
        water_ids = list()
        for protein in ProteinsList:
            protein.calculate_water_coordinates( temp_dir )
            logger.debug( 'Protein %s has %i coordinates.' % (protein, len(protein.water_coordinates)))
            water_coordinates += protein.water_coordinates
            water_ids += protein.water_ids

        if water_coordinates:
            # Only if there are any water molecules list of similar protein structures.
            logger.info( 'Number of water molecules to cluster: %i' % len(water_coordinates) )
            if len(water_coordinates) != 1:
                # Only if the total number of water molecules to cluster is less than 50000.
                if len(water_coordinates) < 50000:
                    cwm_count = 0
                    logger.info( 'Clustering the water coordinates ...' )
                    # The clustering returns a list of clusternumbers
                    # Available optoins are: single, complete, average
                    FD = hcluster.fclusterdata(water_coordinates, 
                            t = ProteinsList.inconsistency_coefficient,
                            criterion='distance', 
                            metric='euclidean',
                            depth=2,
                            method= ProteinsList.clustering_method
                        )
                    FDlist = list(FD)
                    fcDic = {}
                    for a,b in zip(water_ids,FDlist):
                        if b in fcDic:
                            fcDic[b].append(a)
                        else:
                            fcDic[b]=[a]

                    conservedWaterDic = {}
                    clusterPresenceDic = {}
                    clusterPresencePath = os.path.join( outdir, selectedPDBChain, '%s_clusterPresence.txt' % selectedPDBChain )
                    clusterPresenceLines = []
                    clusterPresenceLines.append('Water Conservation Score'+'\t')
                    proteins_numbers = {}
                    l=0
                    for i in ProteinsList.proteins:
                        proteins_numbers[(str(i))]=l
                        l += 1
                        clusterPresenceLines.append('%s' % i +'\t')
                    clusterPresenceLines.append('\n')
                    logger.info( 'Extracting conserved waters from clusters ...' )
                    logger.debug( 'Start iterating over all clusters ...')
                    for clusterNumber, waterMols in fcDic.items():
                        logger.debug('Cluster No.: %s - Included Water Molecules: %s' % ( clusterNumber, ', '.join( waterMols ) ))
                        waterMolsNumber = len(waterMols)
                        uniquePDBs = set([a[:6] for a in waterMols])
                        uniquePDBslen = len(uniquePDBs)
                        # Update waterMols if there are two waterMolecules from same PDB is present in waterMols
                        if uniquePDBslen < waterMolsNumber:
                            logger.debug('Removed water molecules from the same protein in one cluster.')
                            waterMols.sort()
                            PDBs = [i[:6] for i in waterMols]
                            PDBsCounter = collections.Counter(PDBs)
                            for i in PDBsCounter.keys():
                                if PDBsCounter.get(i) > 1:
                                    for j in reversed(waterMols):
                                        if j[:6] == i:
                                            waterMols.remove(j)
                        c = []
                        for i in range(len(ProteinsList.proteins)):
                            for prot in proteins_numbers.keys():
                                if proteins_numbers[prot] == i:
                                    for d in waterMols:
                                        if prot == d[:6]:
                                            c.append(d)
                        waterMols = c
                        waterMolsNumber = len(waterMols)
                        uniquePDBs = set([a[:6] for a in waterMols])
                        uniquePDBslen = len(uniquePDBs)
                        probability = float(uniquePDBslen) / len(ProteinsList)
                        logger.info( 'Degree of conservation is: %s' % probability )
                        if probability >= ProteinsList.probability:
                            clusterPresenceLines.append(str(probability)+'\t')
                            k=0
                            for waterMol in waterMols:
                                sr_no = proteins_numbers[waterMol[:6]]
                                if sr_no == k:
                                    clusterPresenceLines.append(str(waterMol[7:])+'\t')
                                    k += 1
                                elif k < sr_no:
                                    for j in range(sr_no-k):
                                        clusterPresenceLines.append('NoWater'+'\t')
                                        k += 1
                                    clusterPresenceLines.append(str(waterMol[7:])+'\t')
                                    k += 1
                            for j in range(len(ProteinsList.proteins)-k):
                                clusterPresenceLines.append('NoWater'+'\t')
                            clusterPresenceLines.append('\n')

                            if selectedPDBChain in uniquePDBs:
                                cwm_count += 1
                                for waterMol in waterMols:
                                    if waterMol[:6] in conservedWaterDic:
                                        conservedWaterDic[waterMol[:6]].append('_'.join([waterMol[7:], str(probability)]))
                                    else:
                                        conservedWaterDic[waterMol[:6]] = ['_'.join([waterMol[7:], str(probability)])]
                    with open(clusterPresencePath, 'w') as clusterPresenceOut:
                        clusterPresenceOut.write(''.join(clusterPresenceLines))
                    logger.debug( 'conservedWaterDic is: ')
                    for pdb_id, atoms in conservedWaterDic.items():
                        logger.debug('Oxygen atom numbers for %s: %s' % ( pdb_id, ', '.join( atoms ) ))
                    if selectedPDBChain in conservedWaterDic.keys():
                        # save pdb file of only conserved waters for selected pdb
                        atomNumbersProbDic = {}
                        atomNumbers_Prob = conservedWaterDic[selectedPDBChain]
                        logger.info( """Degree of conservation for each conserved water molecule is stored in cwm_%s_withConservedWaters.pdb with the format 'atomNumber'_'DegreeOfConservation'""" % ( selectedPDBChain ) )
                        for probability in atomNumbers_Prob:
                            atom, prob = probability.split('_')
                            atomNumbersProbDic[ atom ] = float( prob )
                        atomNumbers = list(atomNumbersProbDic.keys())
                        conservedWatersPath = os.path.join(temp_dir, 'cwm_'+selectedPDBChain+'_ConservedWatersOnly.pdb')
                        waterPath = os.path.join(temp_dir, 'cwm_'+selectedPDBChain+'_Water.pdb')
                        with open(conservedWatersPath,'w+') as selectedPDBChainConservedWatersOut, open(waterPath) as selectedPDBChainIn:
                            for line in selectedPDBChainIn:
                                if line.startswith('HETATM'):
                                    if str(int(line[22:30])) in atomNumbers:
                                        selectedPDBChainConservedWatersOut.write( line )
                            selectedPDBChainConservedWatersOut.write('END')

                        # add conserved waters to pdb file
                        cmd.delete('cwm_*')
                        cmd.load( os.path.join(temp_dir, 'cwm_%s.pdb' % selectedPDBChain) )
                        cmd.load( os.path.join(temp_dir, 'cwm_%s_ConservedWatersOnly.pdb' % selectedPDBChain) )
                        cmd.remove( 'resname hoh and '+'cwm_%s' % selectedPDBChain )
                        cmd.save( os.path.join(temp_dir, 'cwm_%s_withConservedWaters.pdb' % selectedPDBChain), 'cwm_*')
                        cmd.delete('cwm_*')
                        shutil.copy( os.path.join(temp_dir, 'cwm_%s_withConservedWaters.pdb' % selectedPDBChain), os.path.join(outdir,selectedPDBChain))
                        shutil.copy(os.path.join(temp_dir, 'cwm_%s.pdb' % selectedPDBChain),os.path.join(outdir,selectedPDBChain))
                        if os.path.exists(os.path.join(outdir, selectedPDBChain, 'cwm_%s_withConservedWaters.pdb' % selectedPDBChain)):
                            logger.info( "%s structure has %s conserved water molecules." % (selectedPDBChain,cwm_count))
                            displayInPyMOL(os.path.join(outdir, selectedPDBChain), 'cwm_%s' % selectedPDBChain, atomNumbersProbDic)
                        logger.info("""PDB file of query protein with conserved waters "cwm_%s_withConservedWaters.pdb" and logfile (pywater.log) is saved in %s""" % ( selectedPDBChain, os.path.abspath(outdir)))
                    else:
                        logger.info( "%s has no conserved waters" % selectedPDBChain )
                else:
                    logger.error( "%s has too many waters to cluster. Memory is not enough..." % selectedPDBChain )
            else:
                logger.info( "%s has only one water molecule..." % selectedPDBChain )
        else:
            logger.info( "%s and other structures from the same cluster do not have any water molecules." % selectedPDBChain )
    else:
        logger.error( "%s has only one PDB structure. We need atleast 2 structures to superimpose." % selectedPDBChain )


def pdbIdFormat( pdbId ):
    """
        Check whether the given PDB ID is valid or not.
    """
    if not re.compile('^[a-z0-9]{4}$').match(pdbId):
        logger.error( 'The entered PDB id %s is not valid.' % pdbId)
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The entered PDB id is not valid.""")
        return False
    else:
        return True


def chainIdFormat(chainId):
    """
        Check whether the given PDB Chain ID is valid or not.
    """
    if not re.compile('^[A-Z0-9]{1}$').match(chainId):
        logger.error( 'The entered PDB chain id %s is not valid.' % chainId)
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The entered PDB chain id is not valid.""")
        return False
    else:
        return True


def isXray( pdb ):
    """
        Check whether the PDB structure is determined by X-ray or not using the modern RCSB Data API.
    """
    url = "%s/rest/v1/core/entry/%s" % (RCSB_DATA_URL, pdb.lower())
    try:
        data = _get_json(url)
        for exptl in data.get("exptl", []):
            if exptl.get("method") == "X-RAY DIFFRACTION":
                return True
        return False
    except Exception as e:
        logger.error( 'Error checking isXray for %s: %s' % (pdb, e) )
        return False

def chainPresent(pdb,chain):
    """
        Check whether the given chain id is valid for a given PDB ID using the modern RCSB Data API.
    """
    query = """
    {
      entry(entry_id: "%s") {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            auth_asym_ids
          }
        }
      }
    }
    """ % pdb.upper()
    try:
        gql_data = _graphql(query)
        entities = gql_data.get("data", {}).get("entry", {}).get("polymer_entities", [])
        for entity in entities:
            identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
            auth_ids = identifiers.get("auth_asym_ids", [])
            if chain in auth_ids:
                return True
        return False
    except Exception as e:
        logger.error( 'Error checking chainPresent for %s chain %s: %s' % (pdb, chain, e) )
        return False

def fetchpdbChainsList( selectedStruture, seq_id ):
    """
        Fetch sequence cluster data from RCSB PDB for a given query protein using the modern Data/Search API.
    """
    pdb_id, chain_id = selectedStruture.split('.')
    pdbChainsList = []
    query = """
    {
      entry(entry_id: "%s") {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            auth_asym_ids
          }
          rcsb_polymer_entity_group_membership {
            group_id
            similarity_cutoff
          }
        }
      }
    }
    """ % pdb_id.upper()
    cluster_id = None
    try:
        gql_data = _graphql(query)
        entities = gql_data.get("data", {}).get("entry", {}).get("polymer_entities", [])
        
        for entity in entities:
            identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {})
            auth_ids = identifiers.get("auth_asym_ids", [])
            if chain_id in auth_ids:
                memberships = entity.get("rcsb_polymer_entity_group_membership", [])
                if memberships:
                    for membership in memberships:
                        if membership.get("similarity_cutoff") == float(seq_id):
                            cluster_id = membership.get("group_id")
                            break
                break
    except Exception as e:
        logger.error( 'Error fetching cluster ID: %s' % e )

    if not cluster_id:
        return pdbChainsList

    search_payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_group_membership.group_id",
                "operator": "exact_match",
                "value": cluster_id,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True},
    }
    try:
        search_data = _post_json(RCSB_SEARCH_URL, search_payload)
        members = search_data.get("result_set", [])
        
        entity_ids = [m.get("identifier") for m in members if m.get("identifier")]
        
        if entity_ids:
            for i in range(0, len(entity_ids), 100):
                chunk = entity_ids[i:i+100]
                chunk_str = ", ".join(['"%s"' % e for e in chunk])
                raw_query = "{ polymer_entities(entity_ids: [%s]) { rcsb_polymer_entity_container_identifiers { entry_id auth_asym_ids } entry { exptl { method } } } }" % chunk_str
                bulk_data = _graphql(raw_query)
                bulk_entities = bulk_data.get("data", {}).get("polymer_entities", [])
                
                if not bulk_entities:
                    continue
                    
                for poly in bulk_entities:
                    identifiers = poly.get("rcsb_polymer_entity_container_identifiers", {})
                    entry_info = poly.get("entry", {})
                    if not identifiers or not entry_info:
                        continue
                        
                    ent_pdb = identifiers.get("entry_id")
                    auth_ids = identifiers.get("auth_asym_ids", [])
                    
                    is_xray = False
                    exptl_list = entry_info.get("exptl") or []
                    for exptl in exptl_list:
                        if exptl.get("method") == "X-RAY DIFFRACTION":
                            is_xray = True
                            break
                            
                    if is_xray and ent_pdb:
                        for auth_id in auth_ids:
                            pdbChainsList.append("%s:%s" % (ent_pdb.lower(), auth_id))
    except Exception as e:
        logger.error( 'Error fetching members of cluster %s: %s' % (cluster_id, e) )
    
    return pdbChainsList


def filterbyResolution( pdbChainsList, resolutionCutoff ):
    """
        Filter the list of PDB structures by given resolution cutoff using modern Data API.
    """
    pdbsList = sorted(set([pdbChain.split(':')[0] for pdbChain in pdbChainsList]))

    pdbsResolution = {}
    
    for i in range(0, len(pdbsList), 100):
        chunk = pdbsList[i:i+100]
        chunk_str = ", ".join(['"%s"' % p.upper() for p in chunk])
        raw_query = "{ entries(entry_ids: [%s]) { rcsb_id rcsb_entry_info { resolution_combined } } }" % chunk_str
        
        try:
            data = _graphql(raw_query)
            entries = data.get("data", {}).get("entries", [])
            
            if not entries:
                continue
                
            for entry in entries:
                pdb = entry.get("rcsb_id", "").lower()
                res_info = entry.get("rcsb_entry_info", {})
                if res_info and "resolution_combined" in res_info and res_info["resolution_combined"]:
                    pdbsResolution[pdb] = res_info["resolution_combined"][0]
                    logger.debug('Resolution of %s is: %s' % (pdb, pdbsResolution[pdb]))
                else:
                    pdbsResolution[pdb] = 'null'
        except Exception as e:
            logger.error( 'Error fetching resolution data: %s' % e )

    for pdb in pdbsList:
        if pdb not in pdbsResolution:
            pdbsResolution[pdb] = 'null'

    scored = []
    for pdbChain in pdbChainsList:
        resolution = pdbsResolution[pdbChain.split(':')[0]]
        if resolution != 'null' and resolution is not None:
            if float(resolution) <= resolutionCutoff:
                scored.append((float(resolution), pdbChain))

    # Return best (lowest) resolution first, so a later structure cap keeps the
    # highest-quality structures.
    scored.sort(key=lambda rc: rc[0])
    return [pdbChain for resolution, pdbChain in scored]


# Cap on how many structures a sequence-cluster query may use. Popular protein
# families now return hundreds of homologs, whose pooled waters exceed what the
# clustering step can hold in memory (see the 50000-water limit in
# makePDBwithConservedWaters). Keeping the best-resolution structures up to this
# cap lets broad queries degrade gracefully instead of failing outright.
DEFAULT_MAX_STRUCTURES = 200


def _prioritize_and_cap(pdbChainsList, queryStr, max_structures):
    """
        Force the query chain to the front of the (resolution-sorted) list and
        cap the list to the best `max_structures` entries so broad sequence
        clusters stay within the clustering memory limit. Returns
        (selected_list, capped_from) where capped_from is the pre-cap length if
        a cap was applied, else None.
    """
    pdbChainsList = list(pdbChainsList)
    if queryStr in pdbChainsList:
        pdbChainsList.remove(queryStr)
    pdbChainsList.insert(0, queryStr)
    capped_from = None
    if max_structures and len(pdbChainsList) > max_structures:
        capped_from = len(pdbChainsList)
        pdbChainsList = pdbChainsList[:max_structures]
    return pdbChainsList, capped_from


# Maximum number of local files supported, bounded by the synthetic 4-char id
# scheme (lc00..lc99). Water identity strings depend on a 4-char id + 1-char
# chain (see the [:6]/[7:] slicing and the cwm_????_?.pdb glob).
MAX_LOCAL_STRUCTURES = 100


def collectLocalStructures(local_dir, reference_filename):
    """
        Gather local .pdb files for a conserved-water search over unpublished
        structures. Returns an ordered list of (synthetic_id, source_path) with
        the reference first (synthetic_id 'lc00'), plus the query id 'lc00'.
        Each file is assigned a synthetic 4-char id (lc00..lc99) so the
        downstream water-identity slicing keeps working. Raises ValueError on a
        missing folder, fewer than 2 .pdb files, a reference not in the folder,
        or more than MAX_LOCAL_STRUCTURES files.
    """
    if not os.path.isdir(local_dir):
        raise ValueError('Local files folder not found: %s' % local_dir)
    pdb_files = sorted(
        f for f in os.listdir(local_dir)
        if f.lower().endswith('.pdb') and os.path.isfile(os.path.join(local_dir, f))
    )
    if len(pdb_files) < 2:
        raise ValueError('Need at least 2 .pdb files in %s (found %i).' % (local_dir, len(pdb_files)))
    if len(pdb_files) > MAX_LOCAL_STRUCTURES:
        raise ValueError('Too many local files (%i); the current limit is %i.' % (len(pdb_files), MAX_LOCAL_STRUCTURES))
    ref_base = os.path.basename(reference_filename)
    if ref_base not in pdb_files:
        raise ValueError('Reference file "%s" is not among the .pdb files in %s.' % (ref_base, local_dir))
    ordered = [ref_base] + [f for f in pdb_files if f != ref_base]
    structures = [('lc%02d' % idx, os.path.join(local_dir, fname)) for idx, fname in enumerate(ordered)]
    return structures, 'lc00'


def _chainHasCA(pdb_path, chain):
    """True if the PDB file has at least one CA atom in the given chain."""
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')) and len(line) > 21:
                    if line[12:16].strip() == 'CA' and line[21] == chain:
                        return True
    except (IOError, OSError):
        return False
    return False


def _writeLocalMap(outdir, selectedPDBChain, structures):
    """Write the synthetic-id -> original-filename map so results stay legible."""
    map_dir = os.path.join(outdir, selectedPDBChain)
    if not os.path.exists(map_dir):
        os.makedirs(map_dir)
    lines = ['synthetic_id\toriginal_file']
    lines += ['%s\t%s' % (pdb_id, os.path.basename(src)) for pdb_id, src in structures]
    with open(os.path.join(map_dir, 'local_files_map.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    logger.info('Local file id mapping:\n%s' % '\n'.join(lines))


def FindConservedWatersLocal(local_files_dir, selectedStrutureChain, local_reference, refinement='Mobility', clustering_method='complete', inconsistency_coefficient=2.4, prob=0.7, save_sup_files=True):
    """
        Conserved-water search over a set of LOCAL PDB files (e.g. unpublished
        structures on the user's workstation) instead of RCSB homologs. Every
        .pdb in local_files_dir is compared; local_reference is the query whose
        conserved waters are reported. Reuses the same superimpose/cluster/score
        pipeline as FindConservedWaters (via makePDBwithConservedWaters).
    """
    selectedStrutureChain = str(selectedStrutureChain).upper()
    if not chainIdFormat(selectedStrutureChain):
        return None
    if inconsistency_coefficient > 2.8:
        logger.info( 'The maximum allowed inconsistency coefficient threshold is 2.8 A' )
        tkMessageBox.showinfo(title = 'Error message',
            message = """The maximum allowed inconsistency coefficient threshold is 2.8 A.""")
        return None
    if prob > 1.0 or prob < 0.4:
        logger.info( 'The degree of conservation is allowed from 0.4 A to 1.0 A.' )
        tkMessageBox.showinfo(title = 'Error message',
            message = """The degree of conservation is allowed from 0.4 A to 1.0 A.""")
        return None
    try:
        structures, query_id = collectLocalStructures(local_files_dir, local_reference)
    except ValueError as e:
        logger.error('Local files error: %s' % e)
        tkMessageBox.showinfo(title = 'Error message', message = str(e))
        return None

    displayInputs(query_id, selectedStrutureChain, 'n/a (local)', 'n/a (local)', refinement,
        ', '.join(os.path.basename(src) for _, src in structures),
        clustering_method, inconsistency_coefficient, prob)

    tmp_dir = tempfile.mkdtemp()
    try:
        up = ProteinsList(ProteinName = '.'.join([query_id, selectedStrutureChain]))
        up.refinement = refinement
        up.probability = prob
        up.clustering_method = clustering_method
        up.inconsistency_coefficient = inconsistency_coefficient
        up.selectedPDBChain = Protein(query_id, selectedStrutureChain)
        selectedPDBChain = str(up.selectedPDBChain)
        logger.info( 'Local mode: %i structures from %s (query %s = %s, chain %s).' % (
            len(structures), local_files_dir, query_id, os.path.basename(local_reference), selectedStrutureChain) )
        for pdb_id, src in structures:
            up.add_protein_from_string('%s:%s' % (pdb_id, selectedStrutureChain))
            shutil.copyfile(src, os.path.join(tmp_dir, pdb_id + '.pdb'))

        # Keep only structures that actually contain the requested chain.
        usable = []
        skipped = []
        names = dict((pid, os.path.basename(src)) for pid, src in structures)
        for protein in up.proteins:
            pdb_path = os.path.join(tmp_dir, protein.pdb_filename)
            if _chainHasCA(pdb_path, protein.chain):
                usable.append(protein)
            else:
                skipped.append(protein.pdb_id)
                logger.warning('Excluding %s (%s): chain %s not found.' % (protein.pdb_id, names.get(protein.pdb_id), protein.chain))
        up.proteins = usable

        if query_id in skipped:
            logger.error('The reference structure %s does not contain chain %s. Choose a different reference or chain.' % (os.path.basename(local_reference), selectedStrutureChain))
            return None
        if len(up.proteins) > 1:
            _writeLocalMap(outdir, selectedPDBChain, structures)
            logger.info( 'Save PDB file with conserved water molecules ...' )
            makePDBwithConservedWaters(up, tmp_dir, outdir, save_sup_files)
            # Also expose the result under the reference file's name.
            produced = os.path.join(outdir, selectedPDBChain, 'cwm_%s_withConservedWaters.pdb' % selectedPDBChain)
            if os.path.exists(produced):
                stem = os.path.splitext(os.path.basename(local_reference))[0]
                shutil.copy(produced, os.path.join(outdir, selectedPDBChain, '%s_withConservedWaters.pdb' % stem))
        else:
            logger.error('Not enough local structures with chain %s to superimpose (need at least 2).' % selectedStrutureChain)
    finally:
        shutil.rmtree(tmp_dir)


def FindConservedWaters(selectedStruturePDB,selectedStrutureChain,seq_id,resolution,refinement,user_def_list,clustering_method,inconsistency_coefficient,prob,save_sup_files=True,max_structures=DEFAULT_MAX_STRUCTURES):# e.g: selectedStruturePDB='3qkl',selectedStrutureChain='A'
    """
        The main function: Identification of conserved water molecules from a given protein structure.
    """
    try:
        with urllib.urlopen(RCSB_DATA_URL, timeout=HTTP_TIMEOUT):
            pass
    except Exception:
        logger.error('The PDB webserver is not reachable.')
        return None
    if not pdbIdFormat(selectedStruturePDB):
        return None
    if not isXray(selectedStruturePDB):
        logger.error( 'The entered PDB structure is not determined by X-ray crystallography.' )
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The entered PDB structure is not determined by X-ray crystallography.""")
        return None
    if not chainIdFormat(selectedStrutureChain):
        return None
    if not chainPresent(selectedStruturePDB,selectedStrutureChain):
        logger.error( 'The entered PDB chain id is not valid for given PDB.' )
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The entered PDB chain id is not valid for given PDB.""")
        return None
    if seq_id not in ['30', '40', '50', '70', '90', '95', '100']:
        logger.error( 'The entered sequence identity value is not valid. Please enter a value from list 30, 40, 50, 70, 90, 95 or 100' )
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The entered sequence identity value is not valid. Please enter a value from list 30, 40, 50, 70, 90, 95 or 100.""")
        return None
    if resolution > 3.0:
        logger.info( 'The maximum allowed resolution cutoff is 3.0 A' )
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The maximum allowed resolution cutoff is 3.0 A.""")
        return None
    UD_pdbChainsList = []
    if user_def_list != '':
        if len(user_def_list.split(',')) != 1:
            for i in user_def_list.split(','):
                if len(i.split('_')) == 2:
                    pdbid,chainid = i.split('_')[0].lower(),i.split('_')[1].upper()
                    j = ':'.join([pdbid,chainid])
                    if not pdbIdFormat(pdbid):
                        return None
                    if not chainIdFormat(chainid):
                        return None
                    else:
                        UD_pdbChainsList.append(j)
                else:
                    logger.info( 'Please enter atleast two pdb chains identifier in the format: xxxx_x,yyyy_y,zzzz_z' )
                    tkMessageBox.showinfo(title = 'Error message',
                        message = """Please enter atleast two pdb chains identifier in the format: xxxx_x,yyyy_y,zzzz_z""")
                    return None
        else:
            logger.info( 'Please enter atleast two pdb chains identifier in the format: xxxx_x,yyyy_y,zzzz_z' )
            tkMessageBox.showinfo(title = 'Error message', 
                message = """Please enter atleast two pdb chains identifier in the format: xxxx_x,yyyy_y,zzzz_z""")
            return None
    if inconsistency_coefficient > 2.8:
        logger.info( 'The maximum allowed inconsistency coefficient threshold is 2.8 A' )
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The maximum allowed inconsistency coefficient threshold is 2.8 A.""")
        return None
    if prob > 1.0 or prob < 0.4:
        logger.info( 'The degree of conservation is allowed from 0.4 A to 1.0 A.' )
        tkMessageBox.showinfo(title = 'Error message', 
            message = """The degree of conservation is allowed from 0.4 A to 1.0 A.""")
        return None
    displayInputs(selectedStruturePDB,selectedStrutureChain,seq_id,resolution,refinement,user_def_list,clustering_method,inconsistency_coefficient,prob)

    tmp_dir = tempfile.mkdtemp()
    try:
        selectedStruture = ".".join([selectedStruturePDB.lower(),selectedStrutureChain.upper()]) # 3qkl.A
        up = ProteinsList(ProteinName = selectedStruture) # ProteinsList class instance up
        up.refinement = refinement
        up.probability = prob
        up.clustering_method = clustering_method
        up.inconsistency_coefficient = inconsistency_coefficient
        logger.info( 'selectedStruture is : %s' % selectedStruture )
        up.selectedPDBChain = Protein(selectedStruturePDB, selectedStrutureChain) # up.selectedPDBChain = 3qkl_a
        logger.info( 'up selectedPDBChain is : %s' % up.selectedPDBChain )
        selectedPDBChain = str(up.selectedPDBChain)
        logger.info( 'selectedPDBChain name is : %s' % selectedPDBChain )
        if UD_pdbChainsList == []:
            logger.info( 'Fetching protein chains list from PDB clusters ...' )
            pdbChainsList = fetchpdbChainsList(selectedStruture,seq_id) # ['3QKL:A', '4EKL:A', '3QKM:A', '3QKK:A', '3OW4:A', '3OW4:B', '3OCB:A', '3OCB:B', '4EKK:A', '4EKK:B']
            candidate_count = len(pdbChainsList)
            logger.info( 'Protein chains list contains %i pdb chains.' % candidate_count)
            logger.debug( 'Protein chains list: "%s"' % ', '.join(pdbChainsList))
            logger.info( 'Filtering by resolution ...')
            pdbChainsList = filterbyResolution(pdbChainsList,resolution)
            resolution_count = len(pdbChainsList)
            # make sure query structure is not filtered out, and cap broad
            # clusters to the best-resolution structures.
            queryStr = "%s:%s" % (selectedStruturePDB.lower(), selectedStrutureChain.upper())
            pdbChainsList, capped_from = _prioritize_and_cap(pdbChainsList, queryStr, max_structures)
            if capped_from is not None:
                logger.info( 'Limiting to the %i best-resolution structures (out of %i) to keep clustering within memory limits.' % (max_structures, capped_from) )
            logger.info( 'Structure selection summary: %i candidate chains, %i passed the %.2f A resolution cutoff, %i selected for clustering.' % (candidate_count, resolution_count, resolution, len(pdbChainsList)) )
            logger.debug( 'Filtered protein chains list: "%s"' % ', '.join(pdbChainsList) )
        else:
            pdbChainsList = UD_pdbChainsList
            logger.info( 'Structure selection summary: using %i user-defined pdb chains.' % len(pdbChainsList) )
        for pdbChain in pdbChainsList:
            up.add_protein_from_string(pdbChain)

        if len(up.proteins)>1:
            def download_pdb(protein):
                pdb_path = os.path.join(tmp_dir, protein.pdb_id+'.pdb')
                if not os.path.exists(pdb_path):
                    logger.debug('Retrieving structure: %s' % protein.pdb_id)
                    url = RCSB_DOWNLOAD_URL % protein.pdb_id.upper()
                    try:
                        req = urllib.Request(url, headers=DOWNLOAD_HEADERS)
                        with urllib.urlopen(req, timeout=HTTP_TIMEOUT) as response, open(pdb_path, 'wb') as out_file:
                            out_file.write(response.read())
                    except Exception as e:
                        logger.debug("Failed to download %s: %s" % (protein.pdb_id, e))

            with ThreadPoolExecutor(max_workers=min(30, len(up.proteins))) as executor:
                list(executor.map(download_pdb, up.proteins))

            # Filter out any proteins whose downloads failed. Missing entries
            # (e.g. RCSB 404s) are common for broad clusters, so each miss is
            # logged at DEBUG and summarised once at WARNING instead of raising
            # an alarming per-structure ERROR on otherwise-normal runs.
            successful_proteins = []
            failed_ids = []
            for protein in up.proteins:
                pdb_path = os.path.join(tmp_dir, protein.pdb_filename)
                if os.path.exists(pdb_path):
                    successful_proteins.append(protein)
                else:
                    failed_ids.append(protein.pdb_id)
                    logger.debug('Excluding %s from superimposition because its PDB file could not be downloaded.' % protein.pdb_id)

            up.proteins = successful_proteins

            if failed_ids:
                logger.warning('%i of %i structures could not be downloaded and were skipped (set logging to DEBUG for the list).' % (len(failed_ids), len(failed_ids) + len(successful_proteins)))

            if len(up.proteins) > 1:
                logger.info( 'Save PDB file with conserved water molecules ...' )
                makePDBwithConservedWaters(up, tmp_dir, outdir, save_sup_files)
            else:
                logger.error('Not enough valid PDB structures remaining to superimpose.')
        else:
            logger.info( "%s has only one PDB structure. We need atleast 2 structures to superimpose." % selectedPDBChain)
    finally:
        shutil.rmtree(tmp_dir)


class _LogCollector(logging.Handler):
    """
        Collects PyWATER log records emitted during one search so the GUI can
        report a meaningful outcome (success, "too many waters", etc.) instead
        of leaving the result only in the console / log file.
    """
    def __init__(self):
        logging.Handler.__init__(self)
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())

    def _find(self, needle):
        for message in reversed(self.messages):
            if needle in message:
                return message
        return None

    def found_waters(self):
        return self._find('conserved water molecules.') is not None

    def summary(self):
        selection_summary = self._find('Structure selection summary:')
        for needle in ('conserved water molecules.',
                       'too many waters to cluster',
                       'has no conserved waters',
                       'do not have any water molecules',
                       'only one water molecule',
                       'only one PDB structure',
                       'not reachable',
                       'not determined by X-ray',
                       'is not valid'):
            message = self._find(needle)
            if message:
                if selection_summary and message != selection_summary:
                    return '%s\n%s' % (selection_summary, message)
                return message
        if selection_summary:
            return selection_summary
        return 'Finished. See ~/PyWATER_outdir/pywater.log for details.'


class ConservedWaters(QtWidgets.QDialog):
    """
        PyWATER input dialog, implemented with Qt (pymol.Qt) so it runs
        natively inside PyMOL's Qt GUI without depending on Tk.
    """

    SEQ_ID_CHOICES = ['30', '40', '50', '70', '90', '95', '100']
    REFINEMENT_CHOICES = ['Mobility', 'Normalized B-factor', 'No refinement']
    CLUSTERING_CHOICES = ['complete', 'average', 'single']

    # Emitted from the worker thread when a search finishes: (message, success).
    _finished = QtCore.Signal(str, bool)

    def __init__(self, parent=None):
        super(ConservedWaters, self).__init__(parent)
        self.setWindowTitle("PyWATER - Find Conserved Waters")
        self._worker_thread = None
        self.makeWindow()
        self._finished.connect(self._on_finished)

    def _help_button(self, callback):
        button = QtWidgets.QPushButton("Help")
        button.clicked.connect(callback)
        return button

    def makeWindow(self):
        grid = QtWidgets.QGridLayout(self)
        row = 0

        # PDB id / Chain id
        grid.addWidget(QtWidgets.QLabel("PDB id"), row, 0)
        self.pdb_id = QtWidgets.QLineEdit()
        grid.addWidget(self.pdb_id, row, 1)
        grid.addWidget(self._help_button(pdb_id_help), row, 2)

        grid.addWidget(QtWidgets.QLabel("Chain id"), row, 3)
        self.chain_id = QtWidgets.QLineEdit()
        grid.addWidget(self.chain_id, row, 4)
        grid.addWidget(self._help_button(chain_help), row, 5)
        row += 1

        # Sequence identity cutoff + user-defined-list checkbox
        grid.addWidget(QtWidgets.QLabel("Sequence identity cutoff"), row, 0)
        self.seq_id = QtWidgets.QComboBox()
        self.seq_id.addItems(self.SEQ_ID_CHOICES)
        self.seq_id.setCurrentText("95")
        grid.addWidget(self.seq_id, row, 1)
        grid.addWidget(self._help_button(seq_id_help), row, 2)

        self.use_user_list = QtWidgets.QCheckBox("User defined pdb-chains list")
        self.use_user_list.stateChanged.connect(self.varcheck)
        grid.addWidget(self.use_user_list, row, 3)
        grid.addWidget(self._help_button(user_defined_lists_help), row, 5)
        row += 1

        # Structure resolution cutoff + user-defined list entry
        grid.addWidget(QtWidgets.QLabel("Structure resolution cutoff"), row, 0)
        self.resolution = QtWidgets.QLineEdit("2.0")
        grid.addWidget(self.resolution, row, 1)
        grid.addWidget(self._help_button(resolution_help), row, 2)

        self.user_list = QtWidgets.QLineEdit()
        self.user_list.setPlaceholderText("xxxx_x,yyyy_y,zzzz_z")
        self.user_list.setEnabled(False)
        grid.addWidget(self.user_list, row, 3, 1, 2)
        row += 1

        # Refinement assessing method
        grid.addWidget(QtWidgets.QLabel("Refinement assessing method"), row, 0)
        self.refinement = QtWidgets.QComboBox()
        self.refinement.addItems(self.REFINEMENT_CHOICES)
        grid.addWidget(self.refinement, row, 1)
        grid.addWidget(self._help_button(refinement_quality_help), row, 2)
        row += 1

        # Clustering method
        grid.addWidget(QtWidgets.QLabel("Clustering method"), row, 0)
        self.clustering = QtWidgets.QComboBox()
        self.clustering.addItems(self.CLUSTERING_CHOICES)
        grid.addWidget(self.clustering, row, 1)
        grid.addWidget(self._help_button(clustering_method_help), row, 2)
        row += 1

        # Inconsistency coefficient threshold
        grid.addWidget(QtWidgets.QLabel("Inconsistency coefficient threshold"), row, 0)
        self.inconsistency = QtWidgets.QLineEdit("2.4")
        grid.addWidget(self.inconsistency, row, 1)
        grid.addWidget(self._help_button(inconsistency_coefficient_help), row, 2)
        row += 1

        # Degree of conservation
        grid.addWidget(QtWidgets.QLabel("Degree of conservation"), row, 0)
        self.prob = QtWidgets.QLineEdit("0.7")
        grid.addWidget(self.prob, row, 1)
        grid.addWidget(self._help_button(prob_help), row, 2)
        row += 1

        # Maximum structures
        grid.addWidget(QtWidgets.QLabel("Maximum structures"), row, 0)
        self.max_structures = QtWidgets.QLineEdit(str(DEFAULT_MAX_STRUCTURES))
        grid.addWidget(self.max_structures, row, 1)
        grid.addWidget(self._help_button(lambda: _qt_info('Maximum structures',
            "Upper limit on how many structures a sequence-identity query uses. "
            "Broad clusters are trimmed to this many best-resolution structures "
            "so clustering stays within memory. Ignored for a user-defined list.")),
            row, 2)
        row += 1

        # Local files mode (skip RCSB; use the user's own structures)
        self.use_local_files = QtWidgets.QCheckBox("Use local files (skip RCSB)")
        self.use_local_files.stateChanged.connect(self.varcheck)
        grid.addWidget(self.use_local_files, row, 0)
        grid.addWidget(self._help_button(local_files_help), row, 2)
        row += 1

        grid.addWidget(QtWidgets.QLabel("Local files folder"), row, 0)
        self.local_dir = QtWidgets.QLineEdit()
        self.local_dir.setPlaceholderText("/path/to/folder containing .pdb files")
        self.local_dir.setEnabled(False)
        grid.addWidget(self.local_dir, row, 1, 1, 3)
        self.browse_button = QtWidgets.QPushButton("Browse...")
        self.browse_button.clicked.connect(self._browse_local_dir)
        self.browse_button.setEnabled(False)
        grid.addWidget(self.browse_button, row, 4)
        row += 1

        grid.addWidget(QtWidgets.QLabel("Reference file"), row, 0)
        self.local_reference = QtWidgets.QLineEdit()
        self.local_reference.setPlaceholderText("my_apo.pdb (query; waters reported on this)")
        self.local_reference.setEnabled(False)
        grid.addWidget(self.local_reference, row, 1, 1, 3)
        row += 1

        # Save superimposed files
        self.save_sup = QtWidgets.QCheckBox("Save superimposed pdb files")
        grid.addWidget(self.save_sup, row, 1)
        grid.addWidget(self._help_button(save_sup_files_help), row, 2)
        row += 1

        # Run
        self.run_button = QtWidgets.QPushButton("Find Conserved Water Molecules")
        self.run_button.clicked.connect(self.run)
        grid.addWidget(self.run_button, row, 1)
        row += 1

        self.status_label = QtWidgets.QLabel("")
        self.status_label.setWordWrap(True)
        grid.addWidget(self.status_label, row, 0, 1, 6)

    def varcheck(self, *args):
        # A user-defined list disables the sequence-identity and resolution
        # filters. Local-files mode disables all the RCSB-only inputs and
        # enables the folder / reference fields instead.
        local = self.use_local_files.isChecked()
        use = self.use_user_list.isChecked()
        # Local-files widgets.
        self.local_dir.setEnabled(local)
        self.browse_button.setEnabled(local)
        self.local_reference.setEnabled(local)
        # RCSB-only widgets (chain id is used in both modes, so stays enabled).
        self.pdb_id.setEnabled(not local)
        self.use_user_list.setEnabled(not local)
        self.max_structures.setEnabled(not local)
        self.user_list.setEnabled(use and not local)
        self.resolution.setEnabled(not use and not local)
        self.seq_id.setEnabled(not use and not local)

    def _browse_local_dir(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select folder with .pdb files")
        if folder:
            self.local_dir.setText(folder)

    def run(self):
        if self._worker_thread is not None and self._worker_thread.is_alive():
            return
        try:
            inconsistency = float(self.inconsistency.text())
            prob = float(self.prob.text())
        except ValueError:
            _qt_info('Invalid input',
                'Inconsistency coefficient and degree of conservation must be numbers.')
            return

        if self.use_local_files.isChecked():
            func = FindConservedWatersLocal
            args = (
                str(self.local_dir.text()),
                str(self.chain_id.text()).upper(),
                str(self.local_reference.text()),
                str(self.refinement.currentText()),
                str(self.clustering.currentText()),
                inconsistency,
                prob,
                bool(self.save_sup.isChecked()),
            )
        else:
            try:
                resolution = float(self.resolution.text())
                max_structures = int(self.max_structures.text())
            except ValueError:
                _qt_info('Invalid input',
                    'Resolution and maximum structures must be numbers.')
                return
            func = FindConservedWaters
            args = (
                str(self.pdb_id.text()).lower(),
                str(self.chain_id.text()).upper(),
                str(self.seq_id.currentText()),
                resolution,
                str(self.refinement.currentText()),
                str(self.user_list.text()),
                str(self.clustering.currentText()),
                inconsistency,
                prob,
                bool(self.save_sup.isChecked()),
                max_structures,
            )

        self.run_button.setEnabled(False)
        self.status_label.setText(
            "Searching... fetching, superimposing and clustering structures. "
            "This can take several minutes for large protein families; the "
            "window stays responsive. Progress is logged to the PyMOL console "
            "and ~/PyWATER_outdir/pywater.log.")
        self._worker_thread = threading.Thread(target=self._worker, args=(func, args))
        self._worker_thread.daemon = True
        self._worker_thread.start()

    def _worker(self, func, args):
        # Runs off the main thread so the Qt GUI does not freeze. PyMOL's cmd
        # API is thread-safe; UI updates are marshalled back to the main thread
        # through the _finished signal.
        records = _LogCollector()
        logger.addHandler(records)
        try:
            func(*args)
            message = records.summary()
            success = records.found_waters()
        except Exception as e:
            logger.error('PyWATER failed: %s' % e)
            message = 'PyWATER failed: %s' % e
            success = False
        finally:
            logger.removeHandler(records)
        self._finished.emit(message, success)

    def _on_finished(self, message, success):
        self.run_button.setEnabled(True)
        self.status_label.setText(message)
        _qt_info('PyWATER', message)


def toPyWATER( v1, v2, v3 = '95', v4 = 2.0, v5 = 'Mobility', v6 = '', v7 = 'complete', v8 = 2.0, v9 = 0.7, v10 = DEFAULT_MAX_STRUCTURES):
    """
        Convert data types of input parameters given by command line.
    """
    selectedStruturePDB = str(v1).lower()
    selectedStrutureChain = str(v2).upper()
    seq_id = str(v3)
    resolution = float(v4)
    refinement = str(v5)
    user_def_list = str(v6)
    clustering_method = str(v7)
    inconsistency_coefficient = float(v8)
    prob = float(v9)
    max_structures = int(v10)
    FindConservedWaters(selectedStruturePDB,selectedStrutureChain,seq_id,resolution,refinement,user_def_list,clustering_method,inconsistency_coefficient,prob,max_structures=max_structures)


def toPyWATERLocal( v1, v2, v3, v4 = 'Mobility', v5 = 'complete', v6 = 2.4, v7 = 0.7, v8 = False ):
    """
        Convert data types for the local-files command line entry point.
        v1 = folder path, v2 = chain id, v3 = reference filename.
    """
    local_files_dir = str(v1)
    selectedStrutureChain = str(v2).upper()
    local_reference = str(v3)
    refinement = str(v4)
    clustering_method = str(v5)
    inconsistency_coefficient = float(v6)
    prob = float(v7)
    save_sup_files = str(v8).lower() in ('1', 'true', 'yes', 'on') if not isinstance(v8, bool) else v8
    FindConservedWatersLocal(local_files_dir, selectedStrutureChain, local_reference, refinement, clustering_method, inconsistency_coefficient, prob, save_sup_files)


# Keep a reference so the dialog is not garbage-collected while open.
_dialog = None


def main(parent=None):
    """Open the PyWATER input dialog (native Qt)."""
    global _dialog
    if parent is None:
        parent = QtWidgets.QApplication.activeWindow()
    _dialog = ConservedWaters(parent)
    _dialog.show()
    _dialog.raise_()
    return _dialog


#Extends PyMOL API to use this tool from command line.
cmd.extend('pywater', toPyWATER)
cmd.extend('pywater_local', toPyWATERLocal)
