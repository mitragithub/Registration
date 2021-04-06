import numpy as np 
import MySQLdb
import glob
import os
import sys


pmdno = sys.argv[1]
if os.path.exists('/sonas-hs/mitra/hpc/home/xli/RegistrationPipelineV3/Data/' + pmdno + '/registration_done.txt'):
    db = MySQLdb.connect('mitramba1.cshl.edu', 'xli', 'mba3515', 'analysis_db')
    cursor = db.cursor()
    sql = "UPDATE brain_tracking SET IsRegistered=1 WHERE brain_id = '" + pmdno +"'"
    print sql
    try:
        cursor.execute(sql)
        db.commit()
    except:
        db.rollback()
    db.close()

    db = MySQLdb.connect('mitramba1.cshl.edu', 'xli', 'mba3515', 'analysis_db')
    cursor = db.cursor()
    sql = "UPDATE IsRegistered_T SET IsRegistered=1 WHERE brain_id = '" + pmdno +"'"
    try:
        cursor.execute(sql)
        db.commit()
    except:
        db.rollback()
    db.close()

