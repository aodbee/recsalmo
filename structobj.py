# Copyrights 2023: Nuttachat Wisittipanit

class RES_ALL:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.res_SeqSero2 = None  # RES_SeqSero2 ins
        self.res_MLST = None  # RES_MLST ins
        self.res_SISTR = None  # RES_SISTR ins
        self.res_AMRFinder = None  # RES_AMRFinder ins
        self.res_SPI = None  # RES_SPI ins
        self.res_CRISPR = None  # RES_CRISPR ins



class RES_SeqSero2:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.identification = None
        self.serotype = None


class RES_SISTR:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.serovar_cgmlst = None
        self.cgmlst_matching_alleles = None
        self.h1 = None
        self.h2 = None
        self.serovar_antigen = None
        self.cgmlst_distance = None
        self.cgmlst_genome_match = None
        self.cgmlst_ST = None
        self.serovar = None
        self.serogroup = None
        self.cgmlst_subspecies = None


class RES_MLST:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.call_ST = None
        self.call_aroC = None
        self.call_dnaN = None
        self.call_hemD = None
        self.call_hisD = None
        self.call_purE = None
        self.call_sucA = None
        self.call_thrA = None


class RES_AMRFinder:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.amrClass = None  # 'class1:subclass1|class2:subclass2'
        self.amrClass_Gene = None


class RES_SPI:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.spiList = None


class RES_CRISPR:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.crisprRes = None
        self.spacerC1 = 'None'
        self.spacerC2 = 'None'
        self.posStartC1 = 'None'
        self.posEndC1 = 'None'
        self.posStartC2 = 'None'
        self.posEndC2 = 'None'
        self.DR_OneC1 = 'None'  # representative
        self.DR_OneC2 = 'None'  # representative
        self.NodeC1 = 'None'  # node
        self.NodeC2 = 'None'  # node
        self.numSpC1 = None
        self.numSpC2 = None





class RES_BLAST:

    def __init__(self, genomeName='X'):
        """
        qacc sallseqid evalue bitscore pident qcovs qcovhsp qstart qend sstart send
        :param genomeName: BK_SAL1
        """
        self.genomeName = genomeName  # 'BK_SAL1'
        # self.SPIindex = SPIindex  # 'SPI1'
        self.qacc = None  # C1_0_SAL_BA3995AA_AS
        self.sallseqid = None  # NODE_1_length_580998_cov_10.4184_ID_1
        self.evalue = None
        self.bitscore = None
        self.pident = None
        self.qcovs = None
        self.qcovhsp = None
        self.qstart = None
        self.qend = None
        self.sstart = None
        self.send = None
        self.sstrand = None  # 'plus' or 'minus'


class RecordRES:

    def __init__(self, genomeName):
        self.genomeName = genomeName
        self.identification = None
        self.serotype = None
        self.call_ST = None
        self.cgmlst_ST = None
        self.serovar_cgmlst = None
        self.serogroup = None
        self.cgmlst_subspecies = None
        self.amrClass = None


class CLocus:

    def __init__(self):
        self.name = None
        self.spList = []  # list of SpacerOne obj


class SpacerOne:

    def __init__(self):
        self.spID = None  # original spacer ID
        self.desc = None  # original desc (Genome Name)
        self.spID_adj = None  # new spacer ID after adjusting locus (reverse to only C1 and C2)
        self.spID_new = None  # new spacer ID after assignment (STM1)
        self.desc_new = None  # new desc
        self.spLength = None
        self.spSequence = None
        self.locusName = None
        self.locusName_new = None
        self.foundMatch = False  # found match in the existing spacer database or not


class LocusOne:

    def __init__(self, locus):
        self.locus = locus  # 'C1','C2'
        self.sstrand = None  # sstrand of the node this locus resides 'plus' or 'minus'
        self.seqStart = None  # start pos on the node
        self.seqEnd = None  # end pos on the node
        self.refStart = None  # start pos on the ref seq
        self.refEnd = None  # end pos on the req seq


class LocusANT:

    def __init__(self):
        self.name = None  # C1,C2 or C3
        self.posStart = None
        self.posEnd = None
        self.DR = []  # list of DR sequences
        self.DR_One = None  # representative of DR
        self.SP = []  # lisf of SpacerOne ins
        self.Node = None  # node name
        # self.sstrand = None  # 'plus' or 'minus'
        self.NodeStrand = None  # strand from cLocusBlastResult 'plus' or 'minus'
        self.NodeStart = None  # start position on node strain
        self.NodeEnd = None  # end position on node strain
        self.locusStrand = None  # strand of locus got from spBlastResult 'plus' or 'minus'
        self.NodeStrand_use = None  # strand of node got from NodeInfo (node strand on LT2)


class LocusUSE:

    def __init__(self):
        self.name = None  # C1,C2
        self.SP = []  # lisf of SpacerOne ins
        self.DR = []  # list of DR sequences
        self.DR_One = None  # representative of DR
        self.NodeList = []  # list of Node
        self.refStart = None  # pos start on LT2
        self.refEnd = None  # pos end on LT2




class NodeInfo:

    def __init__(self):
        self.name = None  # SAL_HC6463AA_AS_NODE_12_length_59827_cov_11.211156
        self.nodeStart = None  # pos start on the node itself
        self.nodeEnd = None  # pos end on the node itself
        self.posStart = None  # pos start on LT2 e.g. 3077012
        self.posEnd = None  # pos end on LT2 e.g. 3027175
        self.sstrand = None  # minus
        self.locusOrder_gtemp = []  # temp original locus order on the node [('C2',4050),('C1',4151)]
        self.locusOrder_g = []  # original locus order on the node ['C2','C1']
        self.locusOrder_use = []  # newly assigned locus on the node ['C1']


class SpacerProfile:

    def __init__(self):
        self.genomeName = None  # genome name
        self.C1 = []  # original spacer ID
        self.C2 = []  # original desc (Genome Name)
        self.C1_pos = []  # [start, end]
        self.C2_pos = []  # [start, end]


class StrainAnn:
    def __init__(self, name):
        self.mainName = name  # main name used for identification: same as AssemblyBarcode
        self.assemblyFile = None  # assembly file name
        self.assemblyFilePath = None  # absolute path of assembly file
        # properties from annotation file
        self.AssemblyBarcode = None
        self.Uberstrain = None  # Uberstrain
        # self.Name = None  # Name
        self.DataSource = None  # Data source
        self.Barcode = None  # Barcode
        self.CollectionYear = None  # Collection Year
        self.Country = None  # Country
        # Result properties
        self.resultOutFolderPath = None  # main output folder path
        self.crisResFile = None  # CRISPR result file
        self.crisResFilePath = None
        self.spResANT = {}
        # self.spRes ->
        # {'NODE_1': {'C1':['TTA','TCG','AAC',..],'C2':['TTA','TCG','AAC',..]},
        #          'NODE_2': {'C1':['TTA','TCG','AAC',..],'C2':['TTA','TCG','AAC',..]},
        # }
        self.spRes = {}
        self.spResObj = {}
        self.resCris = None  # RES_CRISPR ins
        # self.spacerOne_list = []  # list of SpacerOne ins
        # self.spacerRecord = []  # list of SecRecord ins


class spVariant:
    def __init__(self, name, num):
        self.name = name  # '29-34' -> 'numSp1-numSp2'
        self.num = num # number of variant