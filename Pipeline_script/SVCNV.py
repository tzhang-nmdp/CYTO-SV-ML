"""
SVCNV.py
========
SV/CNV data class for parsing VCF records from Lumpy, Delly, Manta, CNVnator.

Classes:
  SVCNV        - Parses a single SV record from VCF format
  SVCNV_merged - Container for merged SV records

Functions:
  merge_by_overlap    - Merge SVs with >= percent reciprocal overlap
  merge_by_breakpoint - Merge SVs with breakpoints within distance bp

Attributes parsed: chr, start_pos, end_pos, length, svcnv_type,
                   ci_start, ci_end, pe, sr, qual, caller_name
"""

#chr: chromsome

class SVCNV:
    def __init__(self,cnv_line):
        self.chr = ""
        self.start_pos = -1
        self.end_pos = -1
        self.length = -1
        self.svcnv_type = ""
        self.ci_start = [0,0]
        self.ci_end = [0,0]
        self.pe = 0
        self.sr = 0
        self.qual = -1.0
        cols = cnv_line.strip().split('\t')
        self.caller_name = cols[-1]
        if cols[-1] in ("lumpy","delly","manta","cnvnator"):
            self.init_lumpy_delly_manta_cnvnator(cols)
        elif cols[-1] in ("ALU","LINE1","SVA"):
            self.init_melt(cols)
            
    def init_lumpy_delly_manta_cnvnator(self,cols):
        self.chr = cols[0]
        self.start_pos = int(cols[1])
        if cols[5] == '.':
            self.qual = -1.0
        else:
            self.qual = float(cols[5])
        info_list = cols[7].split(';')
        info_dict = {}
        for info in info_list:
            tmp = info.split('=')
            if len(tmp) == 2:
                info_dict[tmp[0]] = tmp[1]
        if "END" in info_dict:
            self.end_pos = int(info_dict["END"])
        self.length = self.end_pos - self.start_pos
        if "SVTYPE" in info_dict:
            self.svcnv_type = info_dict["SVTYPE"]
        if "SVLEN" in info_dict:
            self.length = abs(int(info_dict["SVLEN"]))
        if "CIPOS" in info_dict:
            self.ci_start = [int(x) for x in info_dict["CIPOS"].split(",")]
        if "CIEND" in info_dict:
            self.ci_end = [int(x) for x in info_dict["CIEND"].split(",")]
        if "SR" in info_dict:
            self.sr = int(info_dict["SR"])
            
        #different callers use different "PE" header 
        if "PE" in info_dict:
            self.pe = int(info_dict["PE"])
        if "PAIR_COUNT" in info_dict:
            self.pe = int(info_dict["PAIR_COUNT"])
        
    def init_melt(self,cols):
        self.chr = cols[0]
        self.start_pos = self.end_pos = int(cols[1])
        if cols[5] == '.':
            self.qual = -1.0
        else:
            self.qual = float(cols[5])
        self.svcnv_type = "INS"
        info_list = cols[7].split(';')
        info_dict = {}
        for info in info_list:
            tmp = info.split('=')
            if len(tmp) == 2:
                info_dict[tmp[0]] = tmp[1]
        if "SVLEN" in info_dict:
            self.length = abs(int(info_dict["SVLEN"]))
        if "SR" in info_dict:
            self.sr = int(info_dict["SR"])
    
    def __lt__(self, other):
        if self.chr < other.chr:
            return True
        elif self.chr == other.chr:
            if self.start_pos < other.start_pos:
                return True
            else:
                return False
        else:
            return False

class SVCNV_merged:
    def __init__(self,svcnv):
        self.chr = svcnv.chr
        self.start_pos = svcnv.start_pos
        self.end_pos = svcnv.end_pos
        self.svcnv_type = svcnv.svcnv_type
        self.merged_list = [svcnv]

def merge_by_overlap(svcnv_list,percent=0.5):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sm_list = []
    svcnv_merged_current = None
    for idx,svcnv_current in enumerate(svcnv_list_sorted):
        if idx == 0:
            svcnv_merged_current = SVCNV_merged(svcnv_current)
            continue
        if svcnv_merged_current.chr == svcnv_current.chr and svcnv_current.start_pos <= svcnv_merged_current.end_pos:
            overlap_ratio = float(min(svcnv_current.end_pos,svcnv_merged_current.end_pos) - svcnv_current.start_pos) / float(max(svcnv_current.end_pos,svcnv_merged_current.end_pos) - svcnv_merged_current.start_pos)
            if overlap_ratio >= percent:
                svcnv_merged_current.end_pos = max(svcnv_current.end_pos,svcnv_merged_current.end_pos)
                svcnv_merged_current.merged_list.append(svcnv_current)
            else:
                sm_list.append(svcnv_merged_current)
                svcnv_merged_current = SVCNV_merged(svcnv_current)
        else:
            sm_list.append(svcnv_merged_current)
            svcnv_merged_current = SVCNV_merged(svcnv_current)
    sm_list.append(svcnv_merged_current)
    return sm_list

def merge_by_breakpoint(svcnv_list,distance=1000):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sm_list = []
    svcnv_merged_current = None
    for idx,svcnv_current in enumerate(svcnv_list_sorted):
        if idx == 0:
            svcnv_merged_current = SVCNV_merged(svcnv_current)
            continue
        if svcnv_merged_current.chr == svcnv_current.chr and (abs(svcnv_current.start_pos - svcnv_merged_current.start_pos) <= distance or abs(svcnv_current.end_pos - svcnv_merged_current.end_pos) <= distance):
            svcnv_merged_current.end_pos = max(svcnv_current.end_pos,svcnv_merged_current.end_pos)
            svcnv_merged_current.merged_list.append(svcnv_current)
        else:
            sm_list.append(svcnv_merged_current)
            svcnv_merged_current = SVCNV_merged(svcnv_current)
    sm_list.append(svcnv_merged_current)
    return sm_list
    
