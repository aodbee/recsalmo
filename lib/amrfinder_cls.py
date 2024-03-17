# Copyrights 2022: Nuttachat Wisittipanit

from const import (File, AMRFinder_CONST, )
from structobj import RES_AMRFinder
from util import read_file_normal
import os

"""
amrfinder --nucleotide fasta_file --organism Salmonella -o output_file 

"""


class AMRFinder_CLS:

    def __init__(self, fileCls, inputFilePath, outFolderPath, outFolderName):
        self.fileCls = fileCls
        self.inputFilePath = inputFilePath
        self.outFolderPath = outFolderPath
        self.outFolderName = outFolderName  # BK_SAL1
        self.sistr_results = None
        self.allele_results = None
        # self.resSISTR = RES_SISTR(outFolderName)

    def perform(self):
        """
        Perform SeqSero2 run
        :return:
        """
        # self.file_input -> 'BK_CAL1.scaffolds.fasta'
        # SeqSero2_package.py -m k -t 4 -i assembly.fasta -d outFolder
        # run SISTR serovar prediction
        # genome_fasta_path = self.inputFilePath
        # genome_name = self.outFolderName
        # os.system('amrfinder -u')
        # os.system('amrfinder --force_update')
        outFilePath = os.path.join(self.outFolderPath, AMRFinder_CONST.mainResultFile)
        cmdLine = 'amrfinder --nucleotide %s --organism Salmonella -o %s' % (self.inputFilePath,
                                                                             outFilePath,
                                                                            )
        os.system(cmdLine)

    def collectResult(self):
        """
        Collect result
        :return: res -> RES_SeqSero2 obj
        """
        mainResultFilePath = os.path.join(self.outFolderPath, AMRFinder_CONST.mainResultFile)
        # read file
        contentLines = read_file_normal(mainResultFilePath)
        res = RES_AMRFinder(self.outFolderName)
        amrTotal = []
        amrTotal_Gene = []
        for i, line in enumerate(contentLines):
            if i == 0:
                continue
            a = line.split('\t')
            amrGene = a[5]  # tet(B)
            amrCls = a[10]  # AMINOGLYCOSIDE
            amrSubCls = a[11]  # STREPTOMYCIN
            # amrOutput -> tet(B):TETRACYCLINE—TETRACYCLINE
            amrOutput = '%s:%s--%s' % (amrGene, amrCls, amrSubCls,)
            amrOutput_Gene = amrGene
            amrTotal.append(amrOutput)
            amrTotal_Gene.append(amrOutput_Gene)
        # res.amrClass -> tet(B):TETRACYCLINE—TETRACYCLINE|
        res.amrClass = '|'.join(amrTotal)
        res.amrClass_Gene = '|'.join(amrTotal_Gene)
        return res


if __name__ == '__main__':
    obj = AMRFinder_CLS(None, None,None,None)
    obj.perform()
