import os


class File:

    def __init__(self, inputFolder, outputFolder, mainName):
        self.mainFolder = None  # folder where recsalmo.py resides
        self.mainDBFolder = None  # main database folder path
        self.inputFolder = inputFolder
        self.outputFolder = outputFolder
        self.mainName = mainName
        # DATABASE
        self.DBSC_Folder = None
        self.REFDBFolder = None
        self.SPACER_Folder = None
        self.SPIFS_Folder = None
        self.SpacerFasta_File = None
        # REFERENCE GENOMES
        self.refGenome_LT2_FilePath = None
        self.refGenome_CT18_FilePath = None
        self.refGenome_28791_FilePath = None
        # OVERALL OUTPUT
        self.recordFilePath = None
        self.treeOutputFolder = None
        self.CRISPR_TreeFilePath = None

    def setupAllPaths(self):
        self.mainFolder = os.path.dirname(__file__)
        self.mainDBFolder = os.path.join(self.mainFolder, 'maindb')
        self.DBSC_Folder = os.path.join(self.mainDBFolder, 'DBSC')
        self.REFDBFolder = os.path.join(self.mainDBFolder, 'REF')
        self.SPACER_Folder = os.path.join(self.mainDBFolder, 'SPACER')
        self.SPIFS_Folder = os.path.join(self.mainDBFolder, 'SPI_FS')
        self.SpacerFasta_File = os.path.join(self.SPACER_Folder, 'spacerall.fasta')
        # self.refGenome_LT2_FilePath = os.path.join(self.REFDBFolder, 'refgenome_so46998_09.fasta')
        self.refGenome_LT2_FilePath = os.path.join(self.REFDBFolder, 'REF_LT2.fasta')
        self.refGenome_CT18_FilePath = os.path.join(self.REFDBFolder, 'refgenome_CT18.fasta')
        self.refGenome_28791_FilePath = os.path.join(self.REFDBFolder, 'refgenome_gallinarum28791.fasta')
        if self.mainName is None:
            self.recordFilePath = os.path.join(self.outputFolder, 'summary.xlsx')
        else:
            # self.mainName -> 'UK'
            self.recordFilePath = os.path.join(self.outputFolder, 'summary_%s.xlsx' % (self.mainName,))
            self.recordFilePath_AMRandSPI = os.path.join(self.outputFolder, 'summary_amrspi_%s.xlsx' % (self.mainName,))
        self.treeOutputFolder = os.path.join(self.outputFolder, 'treeoutput1')
        self.CRISPR_TreeFilePath = os.path.join(self.outputFolder, 'cripr_tree.png')
        # self.CRISPR_TreeFilePath = os.path.join(self.outputFolder, 'cripr_tree_test.png')


class Parameter:
    class SPI:
        E_VALUE_THRESH = 1e-4

    class SPC:
        E_VALUE_THRESH = 1e-4


class Record:
    genomeName = 'genomeName'
    identification = 'identification'
    serotype = 'serotype'
    call_ST = 'call_ST'
    cgmlst_ST = 'cgmlst_ST'
    serovar_cgmlst = 'serovar_cgmlst'
    serogroup = 'serogroup'
    cgmlst_subspecies = 'cgmlst_subspecies'
    amrClass = 'amrClass'
    recordTitle = 'summary'
    # headLabels = ('Assembly Name', 'Identification', 'Serotype', 'MLST',
    #               'cgMLST', 'Serovar cgMLST', 'Serogroup', 'cgMLST subspecies',
    #               'AMR Class',)
    headLabels = ('Assembly', 'Identification', 'Serotype', 'MLST', 'cgMLST',
                  'AMR:Gene/Antibiotic', 'SPI', 'Spacer-C1', 'Spacer-C2', 'NumSP-C1', 'NumSP-C2',
                  'Pos-C1', 'Pos-C2', 'DR-C1', 'DR-C2',)
    # headLabels = ('Assembly', 'Identification', 'Serotype', 'MLST', 'cgMLST',
    #               'AMR:Gene/Antibiotic', 'SPI', 'Spacer-C1', 'Spacer-C2', 'NumSP-C1', 'NumSP-C2',
    #               'Pos-C1', 'Pos-C2', 'DR-C1', 'DR-C2', 'Node-C1', 'Node-C2')
    headLabels_CRISPR = ('Assembly', 'Spacer-C1', 'Spacer-C2', 'NumSP-C1', 'NumSP-C2',
                         'Pos-C1', 'Pos-C2', 'DR-C1', 'DR-C2', 'Node-C1', 'Node-C2')
    headLabels_AMRandSPI = ('Assembly', 'AMR Gene', 'SPI', )


class ResultKey:
    SeqSero2 = 'SeqSero2'
    SISTR = 'SISTR'
    LIST = (SeqSero2, SISTR,)


class SeqSero2_CONST:
    mainResultFile = 'SeqSero_result.txt'
    identificationKey = 'Predicted identification'
    serotypeKey = 'Predicted serotype'


class MLST_CONST:
    mainResultFile = 'mlst_result.txt'


class AMRFinder_CONST:
    mainResultFile = 'amrfinder_result.txt'


class SPIFinder_CONST:
    mainResultFile = 'spifinder_result.txt'


class CRISPRFinder_CONST:
    mainResultFile = 'crisprfinder_result.txt'


class CRName:
    C1 = 'C1'
    C2 = 'C2'
    C3 = 'C3'
    C4 = 'C4'
    C1X = 'C1X'
    CList = (C1, C2, C3, C4)
    CListUse = (C1, C2,)


if __name__ == '__main__':
    File.inputFolder = '/home/nuttachat/Python/Project/WGS/in_file'
    File.outputFolder = '/home/nuttachat/Python/Project/WGS/out_file'
    File.mainFolder = os.path.dirname(__file__)
    File.mainDBFolder = os.path.join(File.mainFolder, 'maindb')
    # print(File.mainFolder)
    # print(File.mainDBFolder)
    # print(File.SpacerFasta_File)
