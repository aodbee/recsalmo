# Copyrights 2023: Nuttachat Wisittipanit

import os
import itertools
import util
import time
from openpyxl import (Workbook, load_workbook, )
from openpyxl.utils import get_column_letter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import (DistanceMatrix, DistanceTreeConstructor, )
from const import (File, ResultKey, Record, )
from structobj import (RES_ALL, RecordRES, )
from main import WGSSal
from binfutil import PerformListAlignment


class BatchRun:

    def __init__(self, input_path, output_path, mainName=None):
        # self.useInputFolder = File.inputFolder
        # self.useOutputFolder = File.outputFolder
        self.useInputFolder = input_path
        self.useOutputFolder = output_path
        self.mainName = mainName
        self.fileCls = None  # File ins
        self.wgsFileList = []
        self.numGenomes = 0
        self.resultHash = {}  # {'BK_SAL1': RES_ALL ins,'BK_SAL4': RES_ALL ins, ..}
        self.genomeNameList = []  # list of genome names run in order
        self.headLabels = Record.headLabels
        self.recordTitle = Record.recordTitle
        self.recordRESList = []  # RecordRES ins
        self.spProfileList = []  # list of SpacerProfile ins

    def setupFolderPath(self):
        """
        Setup all folder paths for the program
        :return:
        """
        self.fileCls = File(self.useInputFolder, self.useOutputFolder, self.mainName)
        self.fileCls.setupAllPaths()

    def databaseInitialize(self):
        """
        Initialize some databases
        :return:
        """
        # copy 'sptemp.fasta' from 'spacerall.fasta' if not already exists
        util.checkCopySpacerFasta(self.fileCls)
        self.makeSpacerDB()
        self.makeSpacerTempDB()

    def collectInfo(self):
        """
        Collect information from all the input files
        :return:
        """
        self.wgsFileList = [f for f in os.listdir(self.useInputFolder) if
                            os.path.isfile(os.path.join(self.useInputFolder, f))]
        self.wgsFileList.sort()
        self.numGenomes = len(self.wgsFileList)

    def run(self):
        """
        Main run method
        :return:
        """
        # clear DBSC folder first
        util.createFolderAndClear(self.fileCls.DBSC_Folder)
        for i, wgsFile in enumerate(self.wgsFileList):
            wSal = WGSSal(self.fileCls, wgsFile)
            wSal.setupIO()
            # create folder if not already existing
            util.createFolderNotExist(wSal.outFolderPath)
            # if wSal.outFolderName not in ('SAL_YC8093AA_AS','BK_SAL1', ):
            #     continue
            self.makeBlastDB(wSal.inputFilePath, wSal.outFolderName)
            wSal.runSeqSero2()
            wSal.runMLST()
            wSal.runSISTR()
            wSal.runAMRFinderPlus()
            wSal.runSPIFinder()
            wSal.runCrisprFinder()
            self.genomeNameList.append((wSal.outFolderName, wSal.outFolderPath))
            self.collectResult(wSal)
            # if i == 2:
            #     break
        # Perform Spacer Profile Alignment
        finalMatrix = self.performSpacerProfileAlignment()
        if len(finalMatrix) != 0:
            self.constructPhylogeneticTreeSpacerProfile(finalMatrix)
        self.removeTempSpacerDB()
        # Record all results in the excel file
        self.recordResult()
        # Read from the record
        self.readFromRecord()
        # Construct pie charts
        self.constructPieChartEach('serotype')
        # self.constructPieChartEach('identification')
        self.constructPieChartEach('call_ST')
        # Construct ParSNP tree
        self.constructSNPTree()
        # Remove all folders inside /DBSC/
        util.createFolderAndClear(self.fileCls.DBSC_Folder)

    def run_CRISPR(self):
        """
        Main run method
        :return:
        """
        # clear DBSC folder first
        util.createFolderAndClear(self.fileCls.DBSC_Folder)
        for i, wgsFile in enumerate(self.wgsFileList):
            wSal = WGSSal(self.fileCls, wgsFile)
            wSal.setupIO()
            # Test only
            # if wSal.outFolderName not in ('SAL_BA3995AA_AS', ):
            #     continue
            # if wSal.outFolderName not in ('SAL_HC6463AA_AS', ):
            #     continue
            # if wSal.outFolderName not in ('SAL_EA2914AA_AS', ):
            #     continue
            # if wSal.outFolderName not in ('SAL_KA9071AA_AS', ):
            #     continue
            # if wSal.outFolderName not in ('SAL_HC4765AA_AS', 'SAL_HC4766AA_AS', 'SAL_FA1054AA_AS', ):
            #     continue
            # if wSal.outFolderName not in ('SAL_VC9662AA_AS', ):
            #     continue
            # if wSal.outFolderName not in ('SAL_KA9071AA_AS', 'SAL_LA0087AA_AS','SAL_QC3063AA_AS'):
            #     continue
            # if wSal.outFolderName not in ('SAL_BA3995AA_AS', 'SAL_BA
            # 804AA_AS', 'SAL_HC6481AA_AS','SAL_IB8216AA_AS'):
            #     continue
            # if wSal.outFolderName not in (
            # 'SAL_BA3995AA_AS', 'SAL_BA5804AA_AS', 'SAL_HC6463AA_AS', 'SAL_HC6481AA_AS', 'SAL_IB8216AA_AS',
            # 'SAL_EA2914AA_AS','SAL_KA9071AA_AS', 'SAL_LA0087AA_AS','SAL_QC3063AA_AS'):
            #     continue
            # if wSal.outFolderName not in ('SAL_BA3995AA_AS', 'SAL_BA5804AA_AS', 'SAL_EA2914AA_AS'):
            #     continue
            # if wSal.outFolderName not in ('SAL_BA3995AA_AS','SAL_BA5804AA_AS','SAL_HC6463AA_AS','SAL_HC6481AA_AS'):
            #     continue
            # make folder right here
            util.createFolderAndClear(wSal.outFolderPath)
            # make DB in the "/maindb/DBSC/BK_SAL1_DB"
            self.makeBlastDB(wSal.inputFilePath, wSal.outFolderName)
            # call runCrisprFinder method
            wSal.runCrisprFinder()
            self.genomeNameList.append((wSal.outFolderName, wSal.outFolderPath))
            self.collectResult_CRISPR(wSal)
            # if i == 2:
            #     break
        # # Perform Spacer Profile Alignment
        # finalMatrix = self.performSpacerProfileAlignment()
        # if len(finalMatrix) != 0:
        #     self.constructPhylogeneticTreeSpacerProfile(finalMatrix)
        self.removeTempSpacerDB()
        # Record all results in the excel file
        self.recordResult_CRISPR()
        # Remove all folders inside /DBSC/
        util.createFolderAndClear(self.fileCls.DBSC_Folder)

    def run_AMRandSPI(self):
        """
        Main run method
        :return:
        """
        # clear DBSC folder first
        util.createFolderAndClear(self.fileCls.DBSC_Folder)
        for i, wgsFile in enumerate(self.wgsFileList):
            wSal = WGSSal(self.fileCls, wgsFile)
            wSal.setupIO()
            # Test only
            # if wSal.outFolderName not in ('SAL_BA3995AA_AS', 'SAL_BA5799AA_AS',):
            #     continue
            # make folder right here
            util.createFolderAndClear(wSal.outFolderPath)
            # make DB in the "/maindb/DBSC/BK_SAL1_DB"
            self.makeBlastDB(wSal.inputFilePath, wSal.outFolderName)
            # call runCrisprFinder method
            wSal.runAMRFinderPlus()
            wSal.runSPIFinder()
            self.genomeNameList.append((wSal.outFolderName, wSal.outFolderPath))
            self.collectResult_AMRandSPI(wSal)
        # Record all results in the excel file
        self.recordResult_AMRandSPI()
        # Remove all folders inside /DBSC/
        util.createFolderAndClear(self.fileCls.DBSC_Folder)

    def constructSNPTree(self):
        """
        Construct SNP-based tree
        :return:
        """
        cmdLine = 'parsnp -c -r %s -d %s -o %s' % (self.fileCls.refGenome_LT2_FilePath,
                                                   self.useInputFolder,
                                                   self.fileCls.treeOutputFolder,
                                                   )
        os.system(cmdLine)

    def collectResult(self, wSal):
        """
        Collect all results
        :return:
        """
        self.resultHash[wSal.outFolderName] = RES_ALL(wSal.outFolderName)
        # res -> RES_ALL ins
        resAll = self.resultHash[wSal.outFolderName]
        # resAll.res_SeqSero2 -> RES_SeqSero2 ins
        resAll.res_SeqSero2 = wSal.seqsero2_res
        resAll.res_MLST = wSal.mlst_res
        resAll.res_SISTR = wSal.sistr_res
        resAll.res_AMRFinder = wSal.amrfinder_res
        resAll.res_SPI = wSal.spifinder_res
        resAll.res_CRISPR = wSal.crisprfinder_res
        # after all results are collected, delete the temp folder
        tempDir_SISTR = os.path.join(self.fileCls.outputFolder, '%s_%s' % (wSal.outFolderName, 'tmpSistr',))
        util.removeFolderAndContent(tempDir_SISTR)

    def collectResult_CRISPR(self, wSal):
        """
        Collect all results
        :return:
        """
        # self.resultHash -> {'BK_SAL1: RES_ALL ins}
        self.resultHash[wSal.outFolderName] = RES_ALL(wSal.outFolderName)
        # res -> RES_ALL ins
        resAll = self.resultHash[wSal.outFolderName]
        # print('resAll: ', resAll)
        resAll.res_CRISPR = wSal.crisprfinder_res
        # print('resAll.res_CRISPR: ', resAll.res_CRISPR)

    def collectResult_AMRandSPI(self, wSal):
        """
        Collect all results
        :return:
        """
        # self.resultHash -> {'BK_SAL1: RES_ALL ins}
        self.resultHash[wSal.outFolderName] = RES_ALL(wSal.outFolderName)
        # res -> RES_ALL ins
        resAll = self.resultHash[wSal.outFolderName]
        print('resAll: ', resAll)
        resAll.res_AMRFinder = wSal.amrfinder_res
        resAll.res_SPI = wSal.spifinder_res
        print('resAll.res_AMRFinder: ', resAll.res_AMRFinder)
        print('resAll.res_SPI: ', resAll.res_SPI)

    def makeBlastDB(self, scaffoldFilePath, genomeName):
        """
        Need to transform a scaffolds file into database for the search to begin
        For instance, 'BK_SAL1' is associated with the file 'BK_SAL1.scaffolds.fasta'
        This one file is split into 8 sub-files which are all put into "/maindb/DBSC/BK_SAL1" folder.
        :param scaffoldFilePath: e.g. 'in_file/BK_SAL1.scaffolds.fasta'
        :param genomeName: e.g. 'BK_SAL1' or 'SAL_IC4008AA_AS'
        :return:
        """
        srList = []
        for seq_record in SeqIO.parse(scaffoldFilePath, "fasta"):
            try:
                newID = seq_record.id[0:50]
            except:
                newID = seq_record.id
            record = SeqRecord(seq_record.seq, id=newID, description='')
            srList.append(record)
        # Write into the file
        # fastaFilePath = os.path.join(scaffoldFilePath, 'NodeFound.fasta')
        with open(scaffoldFilePath, "w") as output_handle:
            SeqIO.write(srList, output_handle, "fasta")
        # create main db folder e.g. "/maindb/DBSC/BK_SAL1_DB"
        genomeName_DB = '%s_DB' % (genomeName,)  # 'BK_SAL1_DB'
        # dbFolder -> "/maindb/DBSC/BK_SAL1_DB"
        dbFolder = os.path.join(self.fileCls.DBSC_Folder, genomeName_DB)
        util.createFolderAndClear(dbFolder)
        # dbFolder -> "/maindb/DBSC/BK_SAL1_DB/BK_SAL1"
        dbFolder_X = os.path.join(dbFolder, genomeName)
        cmdLine = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out %s' % (scaffoldFilePath,
                                                                             dbFolder_X,
                                                                             )
        os.system(cmdLine)

    def makeRefLT2_DB(self):
        """
        Need to transform a scaffolds file into database for the search to begin
        For instance, 'BK_SAL1' is associated with the file 'BK_SAL1.scaffolds.fasta'
        This one file is split into 8 sub-files which are all put into "/maindb/DBSC/BK_SAL1" folder.
        :param scaffoldFilePath: e.g. 'in_file/BK_SAL1.scaffolds.fasta'
        :param genomeName: e.g. 'BK_SAL1' or 'SAL_IC4008AA_AS'
        :return:
        """
        RefGenome_DB = 'LT2DB'
        # dbFolder -> "/maindb/REF/LT2_DB"
        dbFolder = os.path.join(self.fileCls.REFDBFolder, RefGenome_DB)
        util.createFolderAndClear(dbFolder)
        # dbFolder -> "/maindb/REF/LT2_DB/BK_SAL1"
        dbFolder_X = os.path.join(dbFolder, 'LT2')
        cmdLine = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out %s' % (self.fileCls.refGenome_LT2_FilePath,
                                                                             dbFolder_X,
                                                                             )
        os.system(cmdLine)

    def makeSpacerDB(self):
        """
        Transform a 'spacerall.fasta' into the database used in Blast operation
        :param genomeName: genome name e.g. 'BK_SAL1'
        :return:
        """
        # get the list of scaffolds files in "in_file" folder
        spacerFile = 'spacerall.fasta'  # Spacer file
        spacerFilePath = os.path.join(self.fileCls.SPACER_Folder, spacerFile)
        # create main db folder -> "/SPACER/SP_DB"
        spacerDB_name = 'SP_DB'
        # dbFolder -> "/SPACER/SP_DB"
        dbFolder = os.path.join(self.fileCls.SPACER_Folder, spacerDB_name)
        if os.path.exists(dbFolder) is False:
            # DB not exist yet, create new one
            util.createFolderAndClear(dbFolder)
            # dbFolder -> "/SPACER/SP_DB/SPALL"
            dbFolder_X = os.path.join(dbFolder, 'SPALL')
            cmdLine = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out %s' % (spacerFilePath,
                                                                                 dbFolder_X,
                                                                                 )
            os.system(cmdLine)

    def makeSpacerTempDB(self):
        """
        Transform a 'spacerall.fasta' into the database used in Blast operation
        :param genomeName: genome name e.g. 'BK_SAL1'
        :return:
        """
        # get the list of scaffolds files in "in_file" folder
        spacerFile = 'spacerall.fasta'  # Spacer file
        spacerFilePath = os.path.join(self.fileCls.SPACER_Folder, spacerFile)
        # create main db folder -> "/SPACER/SP_DB"
        spacerDB_name = 'SPTEMP_DB'
        # dbFolder -> "/SPACER/SPTEMP_DB"
        dbFolder = os.path.join(self.fileCls.SPACER_Folder, spacerDB_name)
        util.createFolderAndClear(dbFolder)
        # dbFolder -> "/SPACER/SP_DB/SPALL"
        dbFolder_X = os.path.join(dbFolder, 'SPALL')
        cmdLine = 'makeblastdb -in %s -dbtype nucl -parse_seqids -out %s' % (spacerFilePath,
                                                                             dbFolder_X,
                                                                             )
        os.system(cmdLine)

    def performSpacerProfileAlignment(self):
        """
        Perform spacer profile alignment among all genomes
        :return:
        """
        # print(self.genomeNameList)
        for genomeName, folderPath in self.genomeNameList:
            spProfileFilePath = os.path.join(folderPath, 'spProfile.txt')
            spProfile = util.readSpacerProfile(spProfileFilePath)
            self.spProfileList.append((genomeName, spProfile))
        # find all combinations in pairs
        matrix = []
        finalMatrix = []
        # comb -> [(a,b),(c,d),..]
        combList = list(itertools.combinations(self.spProfileList, 2))
        if len(combList) == 0:
            return finalMatrix
        n = 0
        j = 0
        k = 0
        startNumber = len(self.genomeNameList)  # if len = 5, n -> 4
        u = []
        while True:
            # a -> ('BK_SAL1',spProfile)
            # b -> ('BK_SAL10',spProfile)
            a, b = combList[n]
            # align C1 first
            pair_C1 = PerformListAlignment(a[1].C1, b[1].C1)
            pair_C1.run()
            pair_C2 = PerformListAlignment(a[1].C2, b[1].C2)
            pair_C2.run()
            totalScore = pair_C1.totalScore + pair_C2.totalScore
            u.append(totalScore)
            n += 1
            j += 1
            if j == startNumber - 1:
                startNumber -= 1
                j = 0
                u = [0] * (k + 1) + u
                k += 1
                matrix.append(u)
                u = []
            if n == len(combList):
                matrix.append([0] * (len(self.genomeNameList)))
                break
        # transform to numpy array
        X = np.array(matrix)
        # copy to lower triangle
        X = np.triu(X)
        X = X + X.T - np.diag(np.diag(X))
        # get only the lower triangle
        for i in range(X.shape[0]):
            ab = X[i][0:i + 1]
            ab = ab.tolist()
            finalMatrix.append(list(ab))
        return finalMatrix

    def constructPhylogeneticTreeSpacerProfile(self, distanceMatrix):
        """
        Construct the tree by UPGMA method
        :param distanceMatrix:
        :return:
        """

        def get_label(leaf):
            a = str(leaf.name)
            if a.find('Inner') != -1:
                # find 'Inner'
                return None
            return leaf.name

        pureList = [u[0] for u in self.genomeNameList]
        # print(pureList)
        b = DistanceMatrix(names=pureList, matrix=distanceMatrix)
        # Creating a DistanceTreeConstructor object
        constructor = DistanceTreeConstructor()
        # Construct the phlyogenetic tree using UPGMA algorithm
        UPGMATree = constructor.upgma(b)
        # Draw the phlyogenetic tree
        # Phylo.draw(UPGMATree,axes=fig_axes, do_show=False)
        plt.clf()
        # fig1, ax1 = plt.subplots(figsize=(11, 11))
        fig = plt.figure(figsize=(15, 20), dpi=150)
        fig_axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(UPGMATree, label_func=get_label, axes=fig_axes, do_show=False)
        # save tree to PNG file
        plt.savefig(self.fileCls.CRISPR_TreeFilePath)

    def recordResult(self):
        """
        Record result to excel/csv
        :return:
        """
        wb = Workbook()
        ws_sum = wb.active
        ws_sum.title = self.recordTitle
        ws_sum.column_dimensions[get_column_letter(1)].width = 15
        ws_sum.column_dimensions[get_column_letter(2)].width = 30
        ws_sum.column_dimensions[get_column_letter(3)].width = 25
        currentRow = 1
        for i, label in enumerate(self.headLabels):
            ws_sum.cell(row=currentRow, column=i + 1, value=label)
        currentRow = 2
        for name in self.resultHash:
            resAll = self.resultHash[name]  # RES_ALL ins
            resList = [name,
                       resAll.res_SeqSero2.identification,
                       resAll.res_SeqSero2.serotype,
                       resAll.res_MLST.call_ST,  # MLST
                       resAll.res_SISTR.cgmlst_ST,  # cgMLST
                       resAll.res_AMRFinder.amrClass,  # AMR gene/classes/subclasses
                       resAll.res_SPI.spiList,  # SPI info
                       resAll.res_CRISPR.spacerC1,  # spacers of Crispr locus 1
                       resAll.res_CRISPR.spacerC2,  # spacers of Crispr locus 2
                       ]
            for j, resValue in enumerate(resList):
                ws_sum.cell(row=currentRow, column=j + 1, value=resValue)
            currentCol = len(resList) + 1
            spC1 = resAll.res_CRISPR.spacerC1
            spC2 = resAll.res_CRISPR.spacerC2
            C1_len = len(spC1.split('|')) if spC1 != '' else 0
            C2_len = len(spC2.split('|')) if spC2 != '' else 0
            ws_sum.cell(row=currentRow, column=currentCol, value=C1_len)
            ws_sum.cell(row=currentRow, column=currentCol + 1, value=C2_len)
            rc = resAll.res_CRISPR
            ws_sum.cell(row=currentRow, column=currentCol + 2, value="%s-%s" % (rc.posStartC1, rc.posEndC1,))
            ws_sum.cell(row=currentRow, column=currentCol + 3, value="%s-%s" % (rc.posStartC2, rc.posEndC2,))
            ws_sum.cell(row=currentRow, column=currentCol + 4, value=rc.DR_OneC1)
            ws_sum.cell(row=currentRow, column=currentCol + 5, value=rc.DR_OneC2)
            # ws_sum.cell(row=currentRow, column=currentCol + 6, value=rc.NodeC1)
            # ws_sum.cell(row=currentRow, column=currentCol + 7, value=rc.NodeC2)
            # ws_sum.cell(row=currentRow, column=currentCol + 8, value=rc.spacerC1)
            # ws_sum.cell(row=currentRow, column=currentCol + 9, value=rc.spacerC1)
            currentRow += 1
        wb.save(self.fileCls.recordFilePath)

    def recordResult_CRISPR(self):
        """
        Record result to excel/csv
        :return:
        """
        wb = Workbook()
        # summary sheet
        ws_sum = wb.active
        ws_sum.title = self.recordTitle
        ws_sum.column_dimensions[get_column_letter(1)].width = 15
        ws_sum.column_dimensions[get_column_letter(2)].width = 30
        ws_sum.column_dimensions[get_column_letter(3)].width = 25
        currentRow = 1
        for i, label in enumerate(Record.headLabels_CRISPR):
            ws_sum.cell(row=currentRow, column=i + 1, value=label)
        currentRow = 2
        for name in self.resultHash:
            resAll = self.resultHash[name]  # RES_ALL ins
            resList = [name,
                       resAll.res_CRISPR.spacerC1,  # spacers of Crispr locus 1
                       resAll.res_CRISPR.spacerC2,  # spacers of Crispr locus 2
                       ]
            for j, resValue in enumerate(resList):
                ws_sum.cell(row=currentRow, column=j + 1, value=resValue)
            currentCol = len(resList) + 1
            spC1 = resAll.res_CRISPR.spacerC1
            spC2 = resAll.res_CRISPR.spacerC2
            C1_len = len(spC1.split('|')) if spC1 != '' else 0
            C2_len = len(spC2.split('|')) if spC2 != '' else 0
            ws_sum.cell(row=currentRow, column=currentCol, value=C1_len)
            ws_sum.cell(row=currentRow, column=currentCol + 1, value=C2_len)
            rc = resAll.res_CRISPR
            ws_sum.cell(row=currentRow, column=currentCol + 2, value="%s-%s" % (rc.posStartC1, rc.posEndC1,))
            ws_sum.cell(row=currentRow, column=currentCol + 3, value="%s-%s" % (rc.posStartC2, rc.posEndC2,))
            ws_sum.cell(row=currentRow, column=currentCol + 4, value=rc.DR_OneC1)
            ws_sum.cell(row=currentRow, column=currentCol + 5, value=rc.DR_OneC2)
            ws_sum.cell(row=currentRow, column=currentCol + 6, value=rc.NodeC1)
            ws_sum.cell(row=currentRow, column=currentCol + 7, value=rc.NodeC2)
            currentRow += 1
        # spacer sheet
        currentRow = 1
        ws_sp = wb.create_sheet('spacer')
        for name in self.resultHash:
            resAll = self.resultHash[name]  # RES_ALL ins
            spC1 = resAll.res_CRISPR.spacerC1
            spC2 = resAll.res_CRISPR.spacerC2
            spC1_s = spC1.split('|')
            spC2_s = spC2.split('|')
            ws_sp.cell(row=currentRow, column=1, value=name)
            currentCol = 2
            for k, sp in enumerate(spC1_s):
                ws_sp.cell(row=currentRow, column=k + currentCol, value=sp)
            for k, sp in enumerate(spC2_s):
                ws_sp.cell(row=currentRow + 1, column=k + currentCol, value=sp)
            currentRow += 2
        wb.save(self.fileCls.recordFilePath)

    def recordResult_AMRandSPI(self):
        """
        Record result to excel/csv
        :return:
        """
        wb = Workbook()
        # summary sheet
        ws_sum = wb.active
        ws_sum.title = self.recordTitle
        ws_sum.column_dimensions[get_column_letter(1)].width = 15
        ws_sum.column_dimensions[get_column_letter(2)].width = 50
        ws_sum.column_dimensions[get_column_letter(3)].width = 50
        currentRow = 1
        for i, label in enumerate(Record.headLabels_AMRandSPI):
            ws_sum.cell(row=currentRow, column=i + 1, value=label)
        currentRow = 2
        for name in self.resultHash:
            resAll = self.resultHash[name]  # RES_ALL ins
            resList = [name,
                       resAll.res_AMRFinder.amrClass_Gene,  # tet(B)|...
                       resAll.res_SPI.spiList,  # SPI1|SPI2|SPI3
                       ]
            for j, resValue in enumerate(resList):
                ws_sum.cell(row=currentRow, column=j + 1, value=resValue)
            currentRow += 1
        wb.save(self.fileCls.recordFilePath_AMRandSPI)

    def readFromRecord(self):
        """
        Read summary result from record
        :return:
        """
        wb = load_workbook(self.fileCls.recordFilePath)
        ws = wb[self.recordTitle]
        rows = list(ws.rows)
        # data_x = []
        # data_y = []
        for i in range(len(rows)):
            if i == 0:
                # skip row 0
                continue
            row = rows[i]
            genomeName = row[0].value
            recordRES = RecordRES(genomeName)
            recordRES.identification = row[1].value
            recordRES.serotype = row[2].value
            recordRES.call_ST = row[3].value
            recordRES.cgmlst_ST = row[4].value
            recordRES.serovar_cgmlst = row[5].value
            recordRES.serogroup = row[6].value
            recordRES.cgmlst_subspecies = row[7].value
            recordRES.amrClass = row[8].value
            self.recordRESList.append(recordRES)

    def constructPieChartEach(self, mainName):
        """
        Construct pie chart for each name
        :param mainName: 'serotype','identification'
        :return:
        """
        # create pie chart for serotype
        # define labels
        # count
        countDict = {}  # {'I 4,[5],12':1, 'Agona':2, 'Weltevreden':2}
        # labelSerotype = []  # ['I 4,[5],12','Agona','Weltevreden']
        for record in self.recordRESList:
            # record -> RecordRES ins
            # keyUse = record.serotype
            keyUse = getattr(record, mainName)
            if countDict.get(keyUse) is None:
                # no key yet, initialize with 1
                countDict[keyUse] = 1
            else:
                countDict[keyUse] += 1
        labelType = list(countDict.keys())
        lavelValues = list(countDict.values())
        sumA = sum(lavelValues)
        prList = []
        prSumSoFar = 0.0
        # calculate percentage of each serotype
        for i, name in enumerate(labelType):
            # name -> 'Agona'
            # get the value
            rawValue = countDict[name]  # 2
            if i != len(labelType):
                # NOT last item
                prValue = round(rawValue / sumA, 2)
                prSumSoFar += prValue
            else:
                # last item
                prValue = 100.0 - prSumSoFar
            prList.append(prValue)
        # define Seaborn color palette to use
        colors = sns.color_palette('pastel')[0:len(labelType)]
        plt.clf()
        fig1, ax1 = plt.subplots(figsize=(9, 9))
        # fig1.subplots_adjust(0.3, 0, 1, 1)
        plt.pie(prList, labels=labelType, colors=colors, autopct='%.0f%%', textprops={'fontsize': 11})
        # plt.show()
        saveFile = 'pieChart_%s.png' % (mainName,)
        # saveFile = 'pieChart_%s_test.png' % (mainName,)
        saveFilePath = os.path.join(self.fileCls.outputFolder, saveFile)
        plt.savefig(saveFilePath)

    def removeTempSpacerDB(self):
        """
        Delete the folder /SPACER/SPTEMP_DB and file /SPACER/sptemp.fasta
        :return:
        """
        dbFolder = os.path.join(self.fileCls.SPACER_Folder, 'SPTEMP_DB')
        util.removeFolderAndContent(dbFolder)
        spFastaTempFilePath = os.path.join(self.fileCls.SPACER_Folder, 'sptemp.fasta')
        util.removeFile(spFastaTempFilePath)


if __name__ == '__main__':
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        prog='RECSALMO',
        description='Specialized Python Program for Typing and Characterization of Salmonella Genomes',
        epilog='Thank you for using RECSALMO. Have a superb day!')
    requirement = parser.add_argument_group("requirement")
    requirement.add_argument("input_dir", help='Directory containing Salmonella genome assemblies')
    requirement.add_argument("output_dir", help='Main output directory')
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    if not input_dir.exists():
        # if the input_dir does not exist, raise error
        print("The input directory does not exist")
        raise SystemExit(1)
    if not output_dir.exists():
        # if the output_dir does not exist, create one
        print("The output directory does not exist, creating output folder path %s .." % (output_dir,))
        print('output_dir: ', output_dir)
        os.mkdir(output_dir)

    start_time = time.time()
    obj = BatchRun(input_dir, output_dir)
    obj.setupFolderPath()
    obj.databaseInitialize()
    obj.collectInfo()
    obj.run()
    time_elapse = time.time() - start_time
    print("RECSALMO has finished the analysis of %s Salmonella genomes in %s seconds." % (obj.numGenomes, time_elapse))
    avgTime_genome = round(time_elapse / obj.numGenomes, 2)
    print("Average analysis time per genome equals %s seconds." % (avgTime_genome,))

    # Test Only
    # input_dir = '/home/nuttachat/Python/Project/WGS/in_UK'
    # output_dir = '/home/nuttachat/Python/Project/WGS/out_UK'
    # name = 'UK'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/China'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/China_OUTRS'
    # name = 'China'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/Thai'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/Thai_OUTRS'
    # name = 'Thai'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/USA'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/USA_OUTRS'
    # name = 'USA'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/UK'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/UK_OUTRS'
    # name = 'UK'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/Thai'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/Thai_OUTRSTEST'
    # name = 'Thai'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/Thai'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/Thai_OUTV3'
    # name = 'Thai'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/Thai'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/Thai_OUTV4'
    # name = 'Thai'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/UK'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/UK_OUTV3'
    # name = 'UK'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/UK'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/UK_OUTV4'
    # name = 'UK'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/USA'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/USA_OUTV3'
    # name = 'USA'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/USA'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/USA_OUTV4'
    # name = 'USA'
    # input_dir = '/home/nuttachat/Downloads/SalStrains/China'
    # output_dir = '/home/nuttachat/Downloads/SalStrains/China_OUTV3'
    # name = 'China'

    # obj = BatchRun(input_dir, output_dir, mainName=name)
    # obj.setupFolderPath()
    # obj.databaseInitialize()
    # obj.collectInfo()
    # # obj.makeRefLT2_DB()
    # # obj.run_CRISPR()
    # obj.run_AMRandSPI()
