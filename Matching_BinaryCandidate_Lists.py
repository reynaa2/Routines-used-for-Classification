import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import corner

# Find the number of apogee IDs that are in both Kevin's gaia candidates, my candidates, and ML's list
# def matches(x,y,z,a,b,c,binpercent,ID1, ID2, ID3):
    # # x = biggest list of candidates apogee id
    # # y = 2nd largest list of candidates apogee id
    # # z = smallest list of candidates apogee id
    # # a = biggest location id
    # # b = 2nd biggest location id
    # # c = smallest location id 
    # match_locationID = []
    # match_apogeeID = []
    # match_identifier = []
    # match_id2 = []
    # percent = []
    # x = np.asarray(x) # array of strings
    # y = np.asarray(y)
    # z = np.asarray(z)
    # for i in range(len(z)): # only allow the smaller list length
    #     if z[i] in x:
    #         match_locationID.append(c[i])
    #         match_apogeeID.append(z[i])
    #         match_identifier.append(ID1[i])
    #         match_id2.append(ID2[i])
    #         percent.append(binpercent[i])
    #         # if ID1[i] != ID2[i] or ID1[i] != ID3[i]:
    #         #     match_identifier.append(ID1[i])
    #         # else:
    #         #     match_identifier.append(ID2[i])
    #     if z[i] in y:
    #         print('I pass Gaia')
    #         match_locationID.append(c[i])
    #         match_apogeeID.append(z[i])
    #         match_identifier.append(ID1[i])
    #         match_id2.append(ID3[i])
    #         percent.append(binpercent[i])
    # return match_locationID, match_apogeeID, match_identifier, match_id2, percent

def matches(x,y,z,a,b,c,binpercent,ID1, ID2, ID3, n):
    # x = biggest list of candidates apogee id
    # y = 2nd largest list of candidates apogee id
    # z = smallest list of candidates apogee id
    # a = biggest location id
    # b = 2nd biggest location id
    # c = smallest location id 
    # n = determines if you want specific matches or all the matches 
    match_locationID = []
    match_apogeeID = []
    match_identifier = []
    match_id2 = []
    # match_id3 = []
    percent = []
    x = np.asarray(x) # array of strings
    y = np.asarray(y)
    z = np.asarray(z)
    if n ==1:
        for i in range(len(z)):
            if z[i] in x and z[i] in y: # Check if a star is in ML, JMR, and Gaia lists
                match_locationID.append(c[i])
                match_apogeeID.append(z[i])
                match_identifier.append(ID1[i])
                match_id2.append(ID2[i])
                # match_id3.append(ID3[i])
                percent.append(binpercent[i])
    if n ==2:
        if z[i] in x and z[i] not in y: # check if Gaia is in JMR cuts
            match_locationID.append(a[i])
            match_apogeeID.append(x[i])
            match_identifier.append(ID1[i])
            match_id2.append(ID2[i])
            JG_percent.append(binpercent[i])
    if n ==3:
        if z[i] in y z[i] and not in x: # check if Gaia is in ML cuts
            match_locationID.append(b[i])
            match_apogeeID.append(y[i])
            match_identifier.append(ID1[i])
            match_id2.append(ID2[i])
            ML_percent.append(binpercent[i])
            
    return match_locationID, match_apogeeID, match_identifier, match_id2, percent

# Report candidate stars in a csv file
def csv_writer(filename,header1,header2,header3,header4,header5,x,y,z,u,delimiter_choice):
    cols = [header1, header2]
    # Construct dataframe to store variable data
    dataframe = pd.DataFrame(columns=cols)
    dataframe[header1] = x # Location ID
    dataframe[header2] = y # Apogee ID
    dataframe[header3] = z # Identified Method Tag (M, J, G)
    dataframe[header4] = 3 # Visual Identification Tag assumed to be 3 = binary ; 2 = possible binary ; 1 = Non-binary 
    dataframe[header5] = u # binary percentage from ML
    dataframe.to_csv(filename,sep=delimiter_choice,index_label=False)
    return

# Function to extract column data from csv file
def read_csvfile(filename,header1,header2,header3,header4,header5,header6,header7, header8):
    data = pd.read_csv(filename)
    colm_1 = np.asarray(data[header1])
    colm_2 = np.asarray(data[header2])
    colm_3 = np.asarray(data[header3])
    colm_4 = np.asarray(data[header4])
    colm_5 = np.asarray(data[header5])
    colm_6 = np.asarray(data[header6])
    colm_7 = np.asarray(data[header7])
    colm_8 = np.asarray(data[header8])
    return colm_1, colm_2, colm_3, colm_4, colm_5, colm_6, colm_7, colm_8

### --- Main Rountine ---- ##
# Read in training set binaries via read_csvfile function
tsb_locationid,tsb_apogeeid, tsb_minr51, tsb_minr101, tsb_minr151, tsb_minratio1, tsb_minratio2, tsb_maxXR = read_csvfile('TrainingSet_SmallR_LargeXR.csv','LocationID','ApogeeID','log(R51)',
'log(R101)','log(R151)','log(R151/R101)','log(R101/R51)','log(xr)')
# Give training set an identifier of TS
ts_id = ['TS' for d in range(len(tsb_apogeeid))]

# Read in JMR binary candidates via read_csvfile function
jmr_locationid,jmr_apogeeid, jmr_minr51, jmr_minr101, jmr_minr151, jmr_minratio1, jmr_minratio2, jmr_maxXR = read_csvfile('JMR_Cuts_RecentRevision.csv','LocationID','ApogeeID','min_log(R51)',
'min_log(R101)','min_log(R151)','min_log(R151/R101)','min_log(R101/R51)','max_log(xr)')
# Give JMR Cut Stars an Identifier of J
JMR_id = ['J' for j in range(len(jmr_apogeeid))]

# Read in ML binary canidates vai Pandas
ml_stars = pd.read_csv('KMEVP1-10R51R401NewXR51R401LogAvREDO.csv')
ml_apo = np.asarray(ml_stars['Star_Name'])
ml_loc = np.asarray(ml_stars['ID'])
percent = np.asarray(ml_stars['1'])

# Give ML found candidates an identifier of M
id_ml = ['M' for i in range(len(ml_apo))]
ml_locationid = []
ml_apogeeid = []
binpercent = []
# Only accept candidates that have a binary confidence % of > 90%
for i in range(len(ml_apo)):
    if percent[i] >= 0.90:
        ml_apogeeid.append(ml_apo[i])
        ml_locationid.append(ml_loc[i])
        binpercent.append(percent[i])

# Read in the new file that Kevin made of targets that match JMR cuts and Gaia method
gaia_apogeeid = pd.read_table('binaries_from_gaia_and_Jessica.tbl',header=None, delim_whitespace=True,usecols=[0],dtype=str)
gaia_locationid = pd.read_table('binaries_from_gaia_and_Jessica.tbl',header=None,delim_whitespace=True,usecols=[2],dtype=np.int32)
# gaia_flag = pd.read_csv('binaries_from_gaia_and_Jessica.tbl',header=None,delim_whitespace=True,usecols=[3])
gaia_apogeeid = gaia_apogeeid[0]
gaia_locationid = gaia_locationid[2]

# Give Gaia candidates an identifier of G
gaia_id = ['G' for k in range(len(gaia_apogeeid))]
# check the lengths of the lists
print(len(ml_apogeeid), len(jmr_apogeeid), len(gaia_apogeeid))

# Find the number of matched apogee IDs in Kevin's list, mine, and ML
match_locationid, match_apogeeid, match_ID1, match_ID2, BinaryConfidence = matches(jmr_apogeeid, ml_apogeeid, gaia_apogeeid, jmr_locationid, ml_locationid, gaia_locationid, binpercent, id_ml, JMR_id, gaia_id,1)
print(len(match_apogeeid))

# Find the matches for Gaia and JMR Cuts that are not in ML list
inJList_locationid, inJList_apogeeid, JG_ID1, JG_ID2, JG_binconf = matches(jmr_apogeeid, ml_apogeeid, gaia_apogeeid, jmr_locationid, ml_locationid, gaia_locationid, binpercent, id_ml, JMR_id, gaia_id,2)
print(len(inJList_apogeeid))

# Find the matches for Gaia and ML that are not in JMR Cuts list
inML_locationid, inML_apogeeid, ML_ID1, ML_ID2, MG_binconf = matches(jmr_apogeeid, ml_apogeeid, gaia_apogeeid, jmr_locationid, ml_locationid, gaia_locationid, binpercent, id_ml, JMR_id, gaia_id,3)
print(len(inML_apogeeid))

# Combine the matched ids into combo strings ie: JM if found by JMR Cut and ML method
match_ID = [i + str(match_ID2[0]) for i, n in zip(match_ID1,match_ID2)]
JG_match = 

# Output the results into a csv file via csv_writer function
matched_candidates = csv_writer('Matched_BinaryCandidates_Updated.csv', 'LocationID','ApogeeID','Method_ID','Class','Bin%',match_locationid, match_apogeeid,match_ID,BinaryConfidence,',')

# Output the list that has targets that are in Gaia and JMR Cuts but not ML 
Gaia_inJ_not_ML = csv_writer('GaiaCandidates_in_JRCuts.csv','LocationID','ApogeeID','Method_ID','Class','Bin%',inJList_locationid, inJList_apogeeid, JG_binconf, ',' )

# Output the list that has targets that are in Gaia and ML but not JMR Cuts
Gaia_inML_not_J = csv_writer('GaiaCandidates_in_ML.csv','LocationID','ApogeeID','Method_ID','Class','Bin%',inML_locationid, inML_apogeeid, MG_binconf, ',')