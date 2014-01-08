
"""
Run it by following command:
Install in pymol by following the path: Plugins -> Manage Plugins -> install
Restart the PyMol

Script to find important or conserved waters in protein structure (pdb).
Important or conserved waters are the waters molecules which are present in most or all available pdb structures when superimposed.

Copyright 2013 Hitesh Patel and B. Gruening

"""

import os, urllib, shutil
import scipy.cluster.hierarchy as hcluster
from xml.dom.minidom import parseString
from Tkinter import *
import tkMessageBox

#import __main__
#__main__.pymol_argv = [ 'pymol','-qc'] # Quiet and no GUI -qc

import pymol

#pymol.finish_launching()

"""
    TODO:
    - To make running in Windows and Mac
    - Include in galaxy
    - Write a log file with enough outputs..
    - Give an input parameter for cluster size cutoff
    - Ask for Superimposition options.. like CA or other atoms..
    - make it write out pymol session files with colored and styled protein, waters
    - write out images
"""

def __init__(self):
   self.menuBar.addmenuitem('Plugin', 'command',
                            'Find Conserved Waters',
                            label = 'Find Conserved Waters',
                            command = main)


class ConservedWaters(Frame):
    def __init__(self,parent):
        Frame.__init__(self, parent, background="white")
        self.parent=parent
        self.parent.title("Find Conserved Waters")
        self.grid()
        self.makeWindow()

    def makeWindow(self):
        frame1 = Frame(self.parent)
        frame1.grid()

        Label(frame1, text="PDB id").grid(row=0, column=0, sticky=W)
        Button(frame1,text=" Help  ",command=pdb_id_help).grid(row=0, column=2, sticky=W)
        v1 = StringVar(master=frame1)
        v1.set('3qkl')
        Entry(frame1,textvariable=v1).grid(row=0, column=1, sticky=W)

        Label(frame1, text="Chain").grid(row=1, column=0, sticky=W)
        Button(frame1,text=" Help  ",command=chain_help).grid(row=1, column=2, sticky=W)
        v2 = StringVar(master=frame1)
        v2.set('A')
        Entry(frame1,textvariable=v2).grid(row=1, column=1, sticky=W)

        Label(frame1, text="Sequence identity cutoff").grid(row=2, column=0, sticky=W)
        Button(frame1,text=" Help  ",command=seq_id_help).grid(row=2, column=2, sticky=W)
        v3 = StringVar(master=frame1)
        v3.set("95")
        OptionMenu(frame1, v3, '40','50','60','70','80','90','95','100').grid(row=2, column=1, sticky=W)

        Label(frame1, text="Probability cutoff").grid(row=3, column=0, sticky=W)
        Button(frame1,text=" Help  ",command=prob_help).grid(row=3, column=2, sticky=W)
        v4 = StringVar(master=frame1)
        v4.set("0.7")
        Entry(frame1,textvariable=v4).grid(row=3, column=1, sticky=W)

        frame2 = Frame(self.parent)
        frame2.grid()

        Button(frame2,text=" Find Conserved Water Molecules ",command= lambda: FindConservedWaters(str(v1.get()).lower(),str(v2.get()).upper(),int(v3.get()),float(v4.get()))).grid(row=0, column=1, sticky=W)

def main():
    root = Tk()
    app = ConservedWaters(root)
    root.mainloop()  

if __name__ == '__main__':
    main()

def pdb_id_help():
   tkMessageBox.showinfo(title = 'PDB Identifier', message = "The PDB id of the protein for which you like to find conserved waters.")

def chain_help():
   tkMessageBox.showinfo(title = 'Chain Identifier', message = "The chain identifier of the protein for which you like to find conserved waters in above mentioned PDB.")

def seq_id_help():
   tkMessageBox.showinfo(title = 'Sequence identity cutoff', message = """All the protein structures, clustered by BlastClust, having sequence identity more than given cutoff will be superimposed to find the conserved water molecules in query protein chain.
Minimum suggested Sequence identity cutoff is 95.""")

def prob_help():
   tkMessageBox.showinfo(title = 'Probability cutoff', message = """Water molecules will be considered CONSERVED if their probability of being conserved is above given cutoff.
Value ranges from 0 to 1.
Minimum suggested value is 0.5
""")

def displayInputs(selectedStruturePDB,selectedStrutureChain,seq_id,prob):
    print 'Input values are : '
    print 'PDB id : ' + selectedStruturePDB
    print 'Chain id : ' + selectedStrutureChain
    print 'Seqence identity cutoff : ' + str(seq_id)
    print 'probability cutoff : ' + str(prob)
#   tkMessageBox.showinfo(title = 'Inputs', message = v1+' '+v2+' '+v3+' '+v4 )


class Protein():
    def __init__(self, pdb_id, chain=False):
        self.pdb_id = pdb_id.lower()
        self.chain = chain.lower()
        self.pdb_id_folder = self.pdb_id[1:3]
        self.pdb_filename = self.pdb_id + '.pdb'
        self.water_coordinates = list()
        self.water_ids = list()
        self.waterIDCoordinates = {}

    def __repr__(self):
        return "%s_%s" % (self.pdb_id, self.chain)

    def calculate_water_coordinates(self, tmp_dir = False):
        path = os.path.join( tmp_dir, '%s_Water.pdb' % self.__repr__() )
        print 'making water coordinates dictcionary...'
        i = 1
        for line in open(path):
            line = line.strip()
            if line.startswith('HETATM'):
                key = self.__repr__()+'_'+str(i)
                self.waterIDCoordinates[key] = [line[30:38],line[38:46],line[46:54]]
                self.water_coordinates.append([line[30:38],line[38:46],line[46:54]])
                self.water_ids.append( key )
                i += 1
        return self.waterIDCoordinates


class ProteinsList():
    def __init__(self, ProteinName):
        self.ProteinName = ProteinName
        self.proteins = list()
        self.selectedPDBChain = ''
        self.probability = 0.5

    def add_protein(self, protein):
        self.proteins.add(protein)

    def add_protein_from_string(self, s):
        """
            accepts a '1234:B' formatted string and creates a protein object out of it
            if no double point is found, it is assumed that the whole string is the pdb id
        """
        s = s.strip()
        if s.find(':') > -1:
            pdb_id, chain = s.split(':')
        else:
            pdb_id = s
            chain = False
        assert len(pdb_id) == 4, '%s is not of length 4' % pdb_id
        p = Protein(pdb_id, chain)
        self.proteins.append(p)

    def pop(self, index):
        return self.proteins.pop(index)

    def __iter__(self):
        for protein in self.proteins:
            yield protein

    def __len__(self):
        return len(self.proteins)

    def __getitem__( self, key ) :
        if isinstance( key, slice ) :
            #Get the start, stop, and step from the slice
            return [self.proteins[ii] for ii in xrange(*key.indices(len(self)))]
        elif isinstance( key, int ) :
            if key < 0 : #Handle negative indices
                key += len( self.proteins )
            if key >= len( self.proteins ) :
                raise IndexError, "The index (%d) is out of range."%key
            return self.proteins[key] #Get the data from elsewhere
        else:
            raise TypeError, "Invalid argument type."


def makePDBwithConservedWaters(ProteinsList, temp_dir, outdir):
    print 'prob is....' + str(ProteinsList.probability)
    pymol.cmd.delete('*')
    print 'loading all pdb chains...'
    for protein in ProteinsList:
        pymol.cmd.load(os.path.join(temp_dir, protein.pdb_filename))
        pymol.cmd.create(str(protein), '%s & chain %s' % (protein.pdb_id, protein.chain))
        pymol.cmd.delete( protein.pdb_id )

    print 'Superimposing all pdb chains...'
    for protein in ProteinsList[1:]:
        pymol.cmd.super(str(protein), ProteinsList[0])
        pymol.cmd.orient( ProteinsList[0] )

    print 'Creating new pymol objects of only water molecules for each pdb chains...'
    for protein in ProteinsList:
        pymol.cmd.create('%s_Water' % protein, '%s & resname hoh' % protein)

    print 'saving water molecules and proteins in separate pdb files for each pdb chains...'
    for protein in ProteinsList:
        pymol.cmd.save(os.path.join(temp_dir, '%s.pdb' % protein), str(protein))
        pymol.cmd.save(os.path.join(temp_dir, '%s_Water.pdb' % protein), '%s_Water' % protein)

    pymol.cmd.delete('*')

    water_coordinates = list()
    waterIDCoordinates = {}
    water_ids = list()
    for protein in ProteinsList:
        protein.calculate_water_coordinates( temp_dir )
        print 'protein '+str(protein)+ ' coordinates length is :'+str(len(protein.water_coordinates))
        water_coordinates += protein.water_coordinates
        water_ids += protein.water_ids

    selectedPDBChain = str(ProteinsList.selectedPDBChain)

    if water_coordinates:
        print 'number of water molecules to cluster : '+ str(len(water_coordinates))
        if len(water_coordinates)<20000 and len(water_coordinates) != 1:
            print 'clustering the water coordinates...'
            # returns a list of clusternumbers
            FD = hcluster.fclusterdata(water_coordinates, t=1.0, criterion='distance', metric='euclidean', depth=2, method='single')
            FDlist = list(FD)
            fcDic = {}
            print 'making flat cluster dictionary...'
            for a,b in zip(water_ids,FDlist):
                if fcDic.has_key(b):
                    fcDic[b].append(a)
                else:
                    fcDic[b]=[a]

            conservedWaterDic = {}
            print 'extracting conserved waters from clusters'
            for clusterNumber, waterMols in fcDic.items():
                waterMolsNumber = len(waterMols)
                #print '--',waterMols
                uniquePDBs = len(set([a[:6] for a in waterMols]))
                #print uniquePDBs, set([a[:6] for a in waterMols])
                if uniquePDBs == waterMolsNumber:
                    probability = float(uniquePDBs) / len(ProteinsList)
                    if probability > ProteinsList.probability:
                        print 'probability is : '+str(probability) 
                        #print str(clusterNumber)
                        for waterMol in waterMols:
                            if conservedWaterDic.has_key(waterMol[:6]):
                                conservedWaterDic[waterMol[:6]].append(waterMol[7:])
                            else:
                                conservedWaterDic[waterMol[:6]] = [waterMol[7:]]
                elif uniquePDBs < waterMolsNumber:
                    print 'Warning : cutoff distance is too large...'

            print 'conservedWaterDic keys are: '
            #print conservedWaterDic.keys()
            #selectedPDBChain = str(ProteinsList[0])
            if selectedPDBChain in conservedWaterDic.keys():
                #####  save pdb file of conserved waters only for selected pdb
                atomNumbers = conservedWaterDic[selectedPDBChain]
                selectedPDBChainConservedWatersOut = open(os.path.join(temp_dir, selectedPDBChain+'_ConservedWatersOnly.pdb'),'w+')
                selectedPDBChainIn = open(os.path.join(temp_dir, selectedPDBChain+'_Water.pdb'))
                for line in selectedPDBChainIn:
                    if line.startswith('HETATM'):
                        if line.split()[1] in atomNumbers:
                            selectedPDBChainConservedWatersOut.write( line )
                selectedPDBChainConservedWatersOut.write('END')
                selectedPDBChainConservedWatersOut.close()

                #####  add conserved waters to pdb file
                pymol.cmd.delete('*')
                pymol.cmd.load( os.path.join(temp_dir, selectedPDBChain+'.pdb') )
                pymol.cmd.load( os.path.join(temp_dir, selectedPDBChain+'_ConservedWatersOnly.pdb') )
                pymol.cmd.remove( 'resname hoh and '+selectedPDBChain )
                pymol.cmd.save( os.path.join(temp_dir, selectedPDBChain+'_withConservedWaters.pdb') )
                pymol.cmd.delete('*')
                shutil.copy( os.path.join(temp_dir, selectedPDBChain+'_withConservedWaters.pdb'), outdir )
            else:
                pymol.cmd.delete('*')
                print selectedPDBChain+" has no conserved waters"
                ### save pdb without any water
                pymol.cmd.load(os.path.join(temp_dir, selectedPDBChain+'.pdb'))
                pymol.cmd.remove('resname hoh and '+selectedPDBChain)
                pymol.cmd.save(os.path.join(temp_dir, selectedPDBChain+'_withConservedWaters.pdb'))
                pymol.cmd.delete('*')
                shutil.copy(os.path.join(temp_dir, selectedPDBChain+'_withConservedWaters.pdb'),outdir)
                shutil.copy(os.path.join(temp_dir, selectedPDBChain+'.pdb'),outdir)
        else:
            print selectedPDBChain+" has too many waters to cluster. Memory is not enough... OR only one water molecule..."
    else:
        pymol.cmd.delete('*')
        print selectedPDBChain+" and other structures from the same CD-HIT cluster do not have any water molecules."
        ### save pdb without any water
        pymol.cmd.load(os.path.join(temp_dir, selectedPDBChain+'.pdb'))
        pymol.cmd.remove('resname hoh and '+selectedPDBChain)
        pymol.cmd.save(os.path.join(temp_dir, selectedPDBChain+'_withConservedWaters.pdb'))
        pymol.cmd.delete('*')
        shutil.copy(os.path.join(temp_dir, selectedPDBChain+'_withConservedWaters.pdb'),outdir)
    if os.path.exists(os.path.join(outdir, selectedPDBChain+'_withConservedWaters.pdb')):
        pymol.cmd.load(os.path.join(outdir, selectedPDBChain+'.pdb'))
        pymol.cmd.load(os.path.join(outdir, selectedPDBChain+'_withConservedWaters.pdb'))
        pymol.cmd.show_as('cartoon', selectedPDBChain+'_withConservedWaters')
        pymol.cmd.show(representation = 'dots', selection ='resname hoh and '+selectedPDBChain+'_withConservedWaters')
    shutil.rmtree(temp_dir)


def fetchpdbChainsList(selectedStruture,seq_id):
    pdbChainsList = []
    seqClustAddress = 'http://pdb.org/pdb/rest/sequenceCluster?cluster=%d&structureId=%s' % (seq_id, selectedStruture) #http://pdb.org/pdb/rest/sequenceCluster?cluster=95&structureId=3qkl.A
    seqClustURL = urllib.urlopen(seqClustAddress)
    toursurl_string= seqClustURL.read()
    if toursurl_string.startswith('An error has occurred'):
        return pdbChainsList
    else:
        seqCluster = parseString( toursurl_string )
        for xmlTag in seqCluster.getElementsByTagName('pdbChain'):
            pdbChainsList.append(str(xmlTag.getAttribute('name')).replace('.',':'))
        return pdbChainsList


def FindConservedWaters(selectedStruturePDB='3qkl',selectedStrutureChain='A',seq_id=95,prob=0.7):# e.g: selectedStruturePDB='3qkl',selectedStrutureChain='A'
    online_pdb_db = 'http://www.pdb.org/pdb/files/%s.pdb'
    displayInputs(selectedStruturePDB,selectedStrutureChain,seq_id,prob)
    tmp_dir = './temp/'
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    outdir = './ConservedWaters_plugin_outdir/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    selectedStruture = ".".join([selectedStruturePDB.lower(),selectedStrutureChain.upper()]) # 3qkl:A
    up = ProteinsList(ProteinName = selectedStruture) # ProteinsList class instance up
    up.probability = prob
    print 'selectedStruture is : '
    print selectedStruture
    up.selectedPDBChain = Protein(selectedStruturePDB, selectedStrutureChain) # up.selectedPDBChain = 3qkl_a
    print 'up selectedPDBChain is : '
    print up.selectedPDBChain
    selectedPDBChain = str(up.selectedPDBChain)
    print 'selectedPDBChain name is : '+ selectedPDBChain
    pdbChainsList = fetchpdbChainsList(selectedStruture,seq_id) # ['3QKL:A', '4EKL:A', '3QKM:A', '3QKK:A', '3OW4:A', '3OW4:B', '3OCB:A', '3OCB:B', '4EKK:A', '4EKK:B']
    print 'pdbChainsList is : '
    print pdbChainsList
    for pdbChain in pdbChainsList:
        up.add_protein_from_string(pdbChain)
    print 'up is : '
    print up.proteins #[3qkl_a, 4ekl_a, 3qkm_a, 3qkk_a, 3ow4_a, 3ow4_b, 3ocb_a, 3ocb_b, 4ekk_a, 4ekk_b]
    if len(up.proteins)>1:
        print 'length of ProteinsList is : '+str(len(up.proteins))
        for protein in up:
            print protein
            print protein.pdb_id
            if not os.path.exists(os.path.join(tmp_dir, protein.pdb_id+'.pdb')):
                print 'retrieving pdb from website : '+ protein.pdb_id
                urllib.urlretrieve(online_pdb_db % protein.pdb_id.upper(), os.path.join(tmp_dir, protein.pdb_id+'.pdb'))
        print 'making pdb with conserved waters...'
        makePDBwithConservedWaters(up, tmp_dir,outdir)
    else:
        print selectedPDBChain+" has only one PDB structure. We need atleast 2 structures to superimpose."


################################################################################
################################################################################

