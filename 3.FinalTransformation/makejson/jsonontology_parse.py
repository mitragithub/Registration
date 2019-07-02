# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 19:37:38 2016

@author: balamurali
"""

import json 
#from itertools import islice
#def rec_func(df,id_no,count_n):
#    global df_len
#    global c
#    global val
#    if df['id']==id_no:
#        c=1
#        val = df['parent_structure_id']
#        return 
#    le=len(df['children'])
#    if le ==0 or c==1:
#        return  
#    else:
#        for i in range(le):
#            rec_func(df['children'][i],id_no,count_n)
#        if df['id'] == val:
#            count_n.append([df['id'],df['name']])
#            #print df['id'],df['name']
#            val = df['parent_structure_id']


def details(node):
    #print node
    return {'id':node['id'],'name':node['name'],'acronym':node['acronym'], 'isleaf': len(node['children'])==0 }
    
def dfs_find(node, idval):
    if node['id']==idval:
        return (details(node),None)
    else:
        for ci in node['children']:
            cn = dfs_find(ci,idval)
            #print cn , node['id']
            if cn is not None:
                #print node['id']
                return cn,details(node)
        return None


def regiondetails(elt,level):
    region_name = ' '*level
    if elt['isleaf']:
        region_name = region_name + '*'
    region_name = region_name + "%s (%s)\n" % (elt['acronym'],elt['name'])
    return region_name
    
                
def data_print(cons,level):
    #print res_arr
    if cons == None:
        return None

    car,cdr = cons
    if cdr is None:
        return regiondetails(car, level)
    #print ' '*level + "%s,%d,%d" %(i['name'],i['id'],level)
    #print ' '*level + "%s,%d,%d" %(j['name'],j['id'],level)
    return regiondetails(cdr,level) + data_print(car,level+1)

def term_details(cons):
    if cons is None:
        return (None,0)
    car,cdr = cons
    lev = 1
    while cdr is not None:
        car,cdr = car
	lev = lev + 1
    return car,lev

    
#def dfs_leaf(node_1,ff,level = 0 ):
#    if len(node_1['children'])==0:
#        print node_1['id'],node_1['name']
#        ff.write(str(node_1['id'])+','+node_1['name']+'\t'+str(level)+'\n')
#    for ci in node_1['children']:
#        dfs_leaf(ci,ff,level+1)
#

#ontologyfile = '/home/balamurali/json_data/ABAadultMouseBrainOntology.json'
#ontologyfile = 'LDDMM/Atlas/ABAtlas/annotation_ccf_2015/ABAadultMouseBrainOntology.json' 
ontologyfile = '1.json' 

def ontology_find(id_no):
    with open(ontologyfile) as data_file:    
        data = json.load(data_file)
        df=data['msg'][0]
        
        res_arr = dfs_find(df, id_no)
        #print res_arr
        return term_details(res_arr)
    return (None,0)

import sys
if __name__ =="__main__":
    idno = sys.argv[1]

    with open(ontologyfile) as data_file:    
        data = json.load(data_file)
        df=data['msg'][0]
        
        res_arr = dfs_find(df, int(idno))
        #print res_arr
        print data_print(res_arr,0)
	print term_details(res_arr)


