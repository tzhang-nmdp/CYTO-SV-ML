"""
SVCNV_sim.py
=============
Simplified SV/CNV class for merging SVs across samples in database annotation.

Classes:
  SVCNV        - Parses SV from simplified format with sample ID
  SVCNV_merged - Container for merged SV records with sample tracking

Functions:
  merge_by_overlap    - Merge SVs with reciprocal overlap >= percent
  merge_by_breakpoint - Merge SVs with breakpoints within distance bp

Helper functions:
  getchr, getstart, getend - Sort key extractors for multi-level sorting
"""

#chr: chromsome

class SVCNV:
    def __init__(self,cnv_line):
        cols = cnv_line.strip().split('\t')
        self.sample=cols[-1]
        self.chr=cols[0]
        self.start_pos=int(cols[1])
        self.end_pos=int(cols[2])
        self.svcnv_type=cols[3]
#        self.info1=":".join(cols[4:])
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
        self.sample=svcnv.sample
#        self.info1=svcnv.info1
        self.merged_list = [svcnv]
def getchr(svcnv):
    return svcnv.chr
def getstart(svcnv):
    return svcnv.start_pos
def getend(svcnv):
    return svcnv.end_pos    
    
def merge_by_overlap(svcnv_list,percent):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(sorted(sorted(svcnv_list,key=getend),key=getstart),key=getchr)
    sm_list = []
    SVCNV_merged_current = None
    for idx,svcnv_current in enumerate(svcnv_list_sorted):
        if idx == 0:
            SVCNV_merged_current = SVCNV_merged(svcnv_current)
            continue
        if SVCNV_merged_current.chr == svcnv_current.chr and svcnv_current.start_pos <= SVCNV_merged_current.end_pos:
            overlap_ratio = float(min(svcnv_current.end_pos,SVCNV_merged_current.end_pos) - svcnv_current.start_pos) / float(max(svcnv_current.end_pos,SVCNV_merged_current.end_pos) - SVCNV_merged_current.start_pos+1)
            if overlap_ratio >= percent: #or svcnv_current.end_pos <= SVCNV_merged_current.end_pos:
                SVCNV_merged_current.end_pos = max(svcnv_current.end_pos,SVCNV_merged_current.end_pos)
                SVCNV_merged_current.merged_list.append(svcnv_current)
            else:
                sm_list.append(SVCNV_merged_current)
                SVCNV_merged_current = SVCNV_merged(svcnv_current)
        else:
            sm_list.append(SVCNV_merged_current)
            SVCNV_merged_current = SVCNV_merged(svcnv_current)
    sm_list.append(SVCNV_merged_current)
    return sm_list
 
    
def merge_by_breakpoint(svcnv_list,distance):
    if len(svcnv_list) == 0:
        return []      
    svcnv_list_sorted = sorted(sorted(sorted(svcnv_list,key=getend),key=getstart),key=getchr)
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
