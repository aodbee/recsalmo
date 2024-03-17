# Copyrights 2022: Nuttachat Wisittipanit

from const import (File, MLST_CONST)
from structobj import RES_MLST
from util import read_file_normal

import os

"""
>fastmlst --scheme senterica -to output_file inputfile 

"""


class MLST_CLS:

    def __init__(self, fileCls, inputFilePath, outFolderPath, outFolderName):
        self.fileCls = fileCls
        self.inputFilePath = inputFilePath
        self.outFolderPath = outFolderPath
        self.outFolderName = outFolderName  # BK_SAL1
        # self.sistr_results = None
        # self.allele_results = None
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
        outFilePath = os.path.join(self.outFolderPath, MLST_CONST.mainResultFile)
        cmdLine = 'fastmlst --scheme senterica -to %s %s ' % (outFilePath,
                                                              self.inputFilePath,
                                                              )
        os.system(cmdLine)

    def collectResult(self):
        """
        Collect result
        :return: res -> RES_SeqSero2 obj
        """
        # /out_file/BK_SAL1/SeQSero_result.txt
        mainResultFilePath = os.path.join(self.outFolderPath, MLST_CONST.mainResultFile)
        # print('mainResultFilePath: ', mainResultFilePath)
        # read file
        contentLines = read_file_normal(mainResultFilePath)
        res = RES_MLST(self.outFolderName)
        for i, line in enumerate(contentLines):
            # line -> 'Predicted serotype:	I 4,[5],12:i:-'
            if i != 1:
                continue
            a = line.split(',')  # [BK_SAL1.scaffolds.fasta,senterica,34,10,19,12,9,5,9,2]
            res.call_ST = a[2]
            res.call_aroC = a[3]
            res.call_dnaN = a[4]
            res.call_hemD = a[5]
            res.call_hisD = a[6]
            res.call_purE = a[7]
            res.call_sucA = a[8]
            res.call_thrA = a[9]
        return res


if __name__ == '__main__':
    pass
