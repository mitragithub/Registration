import MySQLdb
import glob
import sys
import os
#Please change the DB credential
mdevel = MySQLdb.connect(host="mitramba1.cshl.edu",user="xli",passwd="mba3515",db="MBAStorageDB")
cmdevel = mdevel.cursor()
qc = MySQLdb.connect(host="mitramba1.cshl.edu",user="xli",passwd="mba3515",db="qc_pipeline")
cqc = qc.cursor()
qc.autocommit(True)
sql = "INSERT INTO Navigator_brain_info (brain_id,tracer,plane,locationM) VALUES (%s,%s,%s,%s)"


#cmdevel.execute("SELECT DISTINCT injection.tracer ,brain.sectionPlane,section.path FROM Navigator_brain brain LEFT OUTER JOIN Navigator_injection injection on (injection.brain_id = brain.name) LEFT OUTER JOIN Navigator_section section on (section.brain_id = brain.name AND section.modeIndex=12 and label = 'N') where brain.name=%s")
#all_info = list(cmdevel.fetchall())
#cmdevel.execute("SELECT DISTINCT name from Navigator_brain")
brain = sys.argv[1]
print(brain)
if os.path.exists('/grid/mitra/home/xli/RegistrationPipelineV3/Data/' + brain + '/registration_done.txt'):

    cqc.execute("delete from Navigator_brain_info where brain_id = '%s'" % brain)
    cqc.execute("delete from Navigator_brain_qc_section where brain_info_id = '%s'" % brain)
    #cqc.execute(sql,val)
    i = "/nfs/data/main/M32/RegistrationData/Data/" + brain
    print(i)
    count = 0
    tracer_count = 0
        #print i
    cmdevel.execute("SELECT DISTINCT brain.name, injection.tracer ,brain.sectionPlane,section.path FROM Navigator_brain brain LEFT OUTER JOIN Navigator_injection injection on (injection.brain_id = brain.name) LEFT OUTER JOIN Navigator_section section on (section.brain_id = brain.name AND section.modeIndex=12 and label = 'N') where brain.name like '%s'" % brain)
    output = cmdevel.fetchall()
    all_info = []
    for j in output:
        all_info.append(list(j))
    if len(all_info)==1:
        all_info[0][3] = i
        cqc.execute(sql,all_info[0])
    elif len(all_info)>0:
        temp = all_info[0]
        temp[3] = i
            #print temp
            #print len(all_info)
        for index in range(1,len(all_info)):
            temp[1] = temp[1] + '/ ' + all_info[index][1]
        cqc.execute(sql,temp)

    sql = "INSERT INTO Navigator_brain_qc_section (internal_section_id,brain_info_id,section_name,thickness,sectionQC_fail,sectionRegQC_fail,sectionTranQC_BC_fail) VALUES (%s,%s,%s,20,0,0,0)"

    id = cqc.lastrowid
    print(id)


    sec = glob.glob("/grid/mitra/home/xli/RegistrationPipelineV3/Data/"+brain+"/Registration_OUTPUT/*1_straight.png")
    cmdevel.execute("SELECT id,filename from Navigator_section where brain_id = '%s' " % brain)
    MBA_section = cmdevel.fetchall()
    for section in sec:
        file = '_'.join(section.split('/')[-1].split('_')[:4])
        for item in MBA_section:
            item = list(item)
            if '_'.join(file.split('_')[1:]) in item[1] and file.split('-')[1] in item[1]:
                print('_'.join(file.split('_')[1:]))
                print(file.split('-')[1])
                section = "/nfs/data/main/M32/RegistrationData/Data/" + brain + '/Registration_OUTPUT/' + os.path.basename(section)
                # print section
                val = (item[0],brain,section)
                print(sql,val)
                cqc.execute(sql,val)
                break
