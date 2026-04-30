"""
sv_database_mapping.py
======================
Annotates nonTRS (DEL/DUP/INV) SVs against a tabix-indexed SV database
using region overlap-based matching.

For each query SV, queries the database via tabix for overlapping SVs,
merges database hits, then computes the overlap ratio between the query
SV and the merged database region. Reports PASS/FAIL based on the
overlap threshold.

Note: This script uses Python 2.7 (commands module) and requires:
  - tabix (from htslib) in PATH
  - SVCNV, SVCNV_sim, SVCNV_set modules

Usage:
  python sv_database_mapping.py -i <query_sv> -t <database.gz> -d <distance> -p <percent> -o <output>

Arguments:
  -i : Input nonTRS SV file (bed-like format)
  -t : Tabix-indexed database file (.gz + .tbi)
  -d : Max distance for breakpoint-based merging (default: 1000)
  -p : Min overlap ratio threshold (default: 0.7)
  -o : Output file path

Output columns:
  chr, start, end, chr, svtype, sv_id, PASS/FAIL, info, overlap_ratio, db_matches
"""

import sys,getopt,os,commands,copy,subprocess,SVCNV,SVCNV_sim,SVCNV_set
#parameter setting
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"i:d:p:t:o:")
inFile = ""
percent = 0.7
distance = 1000

for op, value in opts:
	if op == "-i":
	    inFile = value
	if op == "-d":
	    distance = int(value)
	if op == "-p":
	    percent = float(value)
	if op == "-t":
	    template_database = value
	if op == "-o":
	    filter_pass = value

if inFile == "":
	print("-i invalid")
	sys.exit()
      
fp=open(filter_pass,'w')

# read sv_list from analysis file
with open(inFile,'r') as f:
    for line in f:
        if not line.startswith('#') and line[:3]=="chr":
# read basic sv info
            item = line.strip().split("\t")
            chr = item[0]    #chrom
            start = min(int(item[1]),int(item[2]))     #start_pos        
            end = max(int(item[1]),int(item[2]))     #end_pos
            svtype=item[3]  
            sv_id=item[4]            
            sub_len=int(start-end)

# prepare mapping sv range
            start1 = max(0,int(start-1000)) # min(int(start-(end-start)*percent),int(start-1000))
            if start1<0:
                start1=0
            end1 = int(end+1000) # max(int(end+(end-start)*percent),int(end+1000))

# initiate sv 
            svcnv1=SVCNV_set.SVCNV(line)
            svcnv_list=[]
            sv_dict_ori=[]
#            print(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos))
# read database sv_list in the range into dictionary
            svcnv_list = commands.getoutput("tabix -f " + template_database +" "+ chr + ":" + str(start1) + "-" + str(end1) +"|grep "+svtype).split("\n")
#             #svcnv_list = subprocess.check_output("tabix -f " + template_database +" "+ chr + ":" + str(start1) + "-" + str(end1) +"|grep "+svtype).split("\n")
#             process = subprocess.Popen(["tabix -f ",template_database," ",chr, ":",str(start1),"-",str(end1),"|grep ",svtype], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             out, err = process.communicate()
#             svcnv_list = out.split("\n")
            if len(svcnv_list)==1 and svcnv_list[0]=="":
                fp.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+chr+"\t"+svtype+"\t"+sv_id+"\tNot_in_database\n")
                continue
            else:
                sv_dict=[]
                sv_dict_ori=[]    
#                print(svcnv_list)
                for svcnvs in svcnv_list:
#                    result = svcnvs.strip().split("\t") 
                    s = SVCNV_set.SVCNV(svcnvs)
                    sv_dict.append(s)

# save sv info from class instance as string into a list
                for smo in sv_dict:
                    sv_dict_ori.append(smo.chr + ":" + str(smo.start_pos) + ":" + str(smo.end_pos) + ":" + str(smo.length) + ":" + smo.svcnv_type + ":" + smo.info) 

# simplify sv_dict from database
                if svtype=="DEL" or svtype=="DUP": 
                    sm1=SVCNV_set.merge_by_overlap(sv_dict,percent)
                else:
                    sm1=SVCNV_set.merge_by_breakpoint(sv_dict,distance)

# calculate svcnv exclude ratio
                sm2=[]
                sv_dict_sim=[]
                subtract_list_len_list=[]
                svcnv_ratio=0                
                svcnv2=copy.copy(svcnv1)
#                 for sv in sm1:
#                     print("sim:"+sv.chr + "\t" + str(sv.start_pos) + "\t" + str(sv.end_pos))
# check the difference of original sv and total databse svs by mapping out the gap sv          
                if svtype=="DEL" or svtype=="DUP" or svtype=="INV":
                    sm2= SVCNV_set.subtract_by_overlap(sm1,svcnv2,percent)
                    for sm in sm2:
#                         print("sub:"+sm.chr + "\t" + str(sm.start_pos) + "\t" + str(sm.end_pos))
                        subtract_list_len_list.append(sm.length) 
                else:
                    sm2 = SVCNV_set.subtract_by_breakpoint(sm1,svcnv2,distance)
                    for sm in sm2:  
#                         print("sub:"+sm.chr + "\t" + str(sm.start_pos) + "\t" + str(sm.end_pos))
                        subtract_list_len_list.append(sm.length)             
# calculate the sv consensus by 1 - the difference between original sv and total databse svs by summing up the gap sv length 
            if len(sm2)==0 or svcnv1.length==0:
                svcnv_ratio=1-float(len(filter(None,sm2)))
            else:                        
                svcnv_ratio=float(1-float(min(subtract_list_len_list))/float(svcnv1.length)) 
#             print(subtract_list_len_list)
#             print(svcnv_ratio)            
#  determine whether the sv diference is passsing the cutoff and save the data file
            if svcnv_ratio<percent:
                fp.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + svcnv1.chr + "\t"+ svcnv1.svcnv_type + "\t" + sv_id + "\tPASS\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio) + "\t" + "|".join(sv_dict_ori) + "\n" )
            if svcnv_ratio>=percent:
                fp.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + svcnv1.chr + "\t" + svcnv1.svcnv_type  + "\t" + sv_id + "\tFAIL\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio) + "\t" + "|".join(sv_dict_ori) + "\n" )
