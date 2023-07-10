# Copyrights 2023: Nuttachat Wisittipanit
"""
Tool Usage (Linux style)
RECSALMO: Rapid Typing and Characterization Tool for Whole Genome Sequencing Data of Salmonella Enterica
>recsalmo -in X1 -out X2 (X1 is a folder, X2 is also a folder created if NOT YET existed)
"""

from const import File
from lib.seqsero2_cls import SeqSero2_CLS
from lib.sistr_cls import SISTR_CLS
from lib.amrfinder_cls import AMRFinder_CLS
from lib.mlst_cls import MLST_CLS
from lib.spifinder_cls import SPIFinder_CLS
from lib.crispr_cls import CrisprCas_CLS
import os
from openpyxl import Workbook
from openpyxl.utils import get_column_letter
from sistr.sistr_cmd import sistr_predict


class WGSSal:

    def __init__(self, fileCls, file_input):
        self.fileCls = fileCls
        self.file_input = file_input
        self.inputFilePath = os.path.join(self.fileCls.inputFolder, self.file_input)
        self.outFolderName = None
        self.outFolderPath = None
        self.mlst_obj = None
        self.mlst_res = None
        self.seqsero2_obj = None
        self.seqsero2_res = None
        self.sistr_obj = None
        self.sistr_res = None
        self.amrfinder_obj = None
        self.amrfinder_res = None
        self.spifinder_obj = None
        self.spifinder_res = None
        self.crisprfinder_obj = None
        self.crisprfinder_res = None

    def setupIO(self):
        """
        Setup IO of input files and main names
        There are many type of file names
        EX.
        (1) 'BK_SAL1.scaffolds.fasta' -> 'BK_SAL1'
        (2) 'SAL_HC4765AA_AS_genomic.fna' -> 'SAL_HC4765AA_AS'
        (3) SAL_IC4008AA_AS.result.fasta -> 'SAL_IC4008AA_AS'

        This means the first text before '.' for each genome file must be unique.
        :return:
        """
        # self.file_input -> 'BK_SAL1.scaffolds.fasta'
        # u -> ['BK_SAL1','scaffolds','fasta']
        #      ['SAL_HC4765AA_AS_genomic','fna']
        #      ['SAL_IC4008AA_AS','result','fasta']
        u = self.file_input.split('.')
        # v -> ['BK','SAL1'] or ['SAL','HC4765AA','AS','genomic']
        v = u[0].split('_')
        if 'AS' not in v:
            # this case, the file name does not contain '_AS_' meaning it is not downloaded from Enterobase
            # thus, it can just use the name before '.'
            name = u[0]
        else:
            # this case, file name contains '_AS_' meaning it is from Enterobase
            # file name -> 'SAL_HC4765AA_AS_genomic.fna'
            name = '_'.join(v[0:2]) + '_AS'
        self.outFolderName = name
        self.outFolderPath = os.path.join(self.fileCls.outputFolder, self.outFolderName)
        # make folder right here
        # util.createFolderAndClear(self.fileCls.DBSC_Folder)

    def runSeqSero2(self):
        """
        Run seqsero2
        :return:
        """
        self.seqsero2_obj = SeqSero2_CLS(self.fileCls, self.inputFilePath, self.outFolderPath, self.outFolderName)
        self.seqsero2_obj.perform()
        # self.seqsero2_res -> RES_SeqSero2 ins
        self.seqsero2_res = self.seqsero2_obj.collectResult()
        # print('identification: ', self.seqsero2_res.identification)
        # print('serotype: ', self.seqsero2_res.serotype)

    def runMLST(self):
        """
        Run amr genes
        :return:
        """
        self.mlst_obj = MLST_CLS(self.fileCls, self.inputFilePath, self.outFolderPath, self.outFolderName)
        self.mlst_obj.perform()
        # .mlst_res -> RES_MLST ins
        self.mlst_res = self.mlst_obj.collectResult()

    def runSISTR(self):
        """
        Run SISTR
        :return:
        """
        self.sistr_obj = SISTR_CLS(self.fileCls, self.inputFilePath, self.outFolderPath, self.outFolderName)
        self.sistr_obj.perform()
        # .sistr_res -> RES_SISTR ins
        self.sistr_res = self.sistr_obj.collectResult()

    def runAMRFinderPlus(self):
        """
        Run AMR finder
        :return:
        """
        # Use blastp
        self.amrfinder_obj = AMRFinder_CLS(self.fileCls, self.inputFilePath, self.outFolderPath, self.outFolderName)
        self.amrfinder_obj.perform()
        # .amrfinder_res -> RES_AMRFinder ins
        self.amrfinder_res = self.amrfinder_obj.collectResult()
        # print('amrClass: ', self.amrfinder_res.amrClass)

    def runSPIFinder(self):
        """
        Run SPI finder
        :return:
        """
        self.spifinder_obj = SPIFinder_CLS(self.fileCls, self.inputFilePath, self.outFolderPath, self.outFolderName)
        self.spifinder_obj.perform()
        # self.spifinder_res -> RES_SPIFinder ins
        self.spifinder_res = self.spifinder_obj.collectResult()

    def runCrisprFinder(self):
        """
        Run crispr cas finder
        :return:
        """
        self.crisprfinder_obj = CrisprCas_CLS(self.fileCls, self.inputFilePath, self.outFolderPath, self.outFolderName)
        self.crisprfinder_obj.perform()
        self.crisprfinder_res = self.crisprfinder_obj.collectResult()



if __name__ == '__main__':
    import sys
    print('pass')
    print(f"Arguments count: {len(sys.argv)}")
    for i, arg in enumerate(sys.argv):
        print(f"Argument {i:>6}: {arg}")