# Copyrights 2023: Nuttachat Wisittipanit


from const import SeqSero2_CONST
from structobj import RES_SeqSero2
from util import read_file_normal
import os


class SeqSero2_CLS:

    def __init__(self, fileCls, inputFilePath, outFolderPath, outFolderName):
        self.fileCls = fileCls
        self.inputFilePath = inputFilePath
        self.outFolderPath = outFolderPath
        self.outFolderName = outFolderName  # BK_SAL1

    def perform(self):
        """
        Perform SeqSero2 run
        :return:
        """
        # self.file_input -> 'BK_CAL1.scaffolds.fasta'
        # SeqSero2_package.py -m k -t 4 -i assembly.fasta -d outFolder
        cmdLine = 'SeqSero2_package.py -m k -t 4 -i %s -d %s' % (self.inputFilePath,
                                                                 self.outFolderPath,
                                                                 )
        os.system(cmdLine)

    def collectResult(self):
        """
        Collect result
        :return: res -> RES_SeqSero2 obj
        """
        # /out_file/BK_SAL1/SeQSero_result.txt
        mainResultFilePath = os.path.join(self.outFolderPath, SeqSero2_CONST.mainResultFile)
        # read file
        contentLines = read_file_normal(mainResultFilePath)
        res = RES_SeqSero2(self.outFolderName)
        for line in contentLines:
            # line -> 'Predicted serotype:	I 4,[5],12:i:-'
            if line.find(SeqSero2_CONST.identificationKey) != -1:
                # find identification key
                b = line.split(':')  # ['Predicted serotype',' I 4,[5],12:i:-']
                identificationValue = b[1].strip()
                res.identification = identificationValue
            if line.find(SeqSero2_CONST.serotypeKey) != -1:
                # find identification key
                b = line.split(':')  # ['Predicted serotype',' I 4,[5],12:i:-']
                serotypeValue = b[1].strip()
                res.serotype = serotypeValue
        return res
