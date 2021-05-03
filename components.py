# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 12:55:23 2021

@author: luc

"""

'''
class list:
    sample: single sampels with information from run planner website
            important infor: id, pool id, run length, request name, recipe, reads num
            function: addinfor: can add infor to sample class by given two list values
            
    function ImportFromExcel: get grouped_sample as output from excel(run-planner)
    
    grouped_samples: sample list with information of totalreads, projectID, is_user_project, barcode list
            function: create from list of samples
                        group samples by isUser, infor(RunLength, RequestID)
                    
            
    ?lane: one lane of a run, have four types, with infor of list of samples
            with function is_valid
    ?run: one run

'''
import pandas as pd


# sample class will contain all the information we can get from run planner
class sample:
    infor_list = ['Pool','SampleID','Recipe','PoolConc.','RequestID','RequestName',
                  'SampleConc.','Units','Volume','BarcodeSequence','RunLength',
                  'ReadsRequested.','ReadsRemaining.','ReadsTotal']
    
    def __init__(self, id):
        self.id = id
        self.infor =  {key: [] for key in self.infor_list}
        self.isUser = ''
        self.isPool = ''
        
    def addInfor(self, key, values):
        for item in zip(key, values):
            if item[0] in self.infor_list:
                self.infor[item[0]] = str(item[1])
                
        if self.infor['Pool'] == 'nan':
            self.isPool = False
        else:
            self.isPool = True
            
        if self.infor['RequestName'] == 'Investigator Prepared Libraries':
            self.isUser = True
        else:
            self.isUser = False
        

# grouped_sample will be a list of samples 
# function cal_totalreads: sums up all the reads numebr from each sample in the list
class grouped_samples(sample):   
    
    def cal_totalreads(self):
        for i in self.obj:
            self.totalreads += float(i.infor['ReadsRequested.'])
    
    def __init__(self, sample_list):
        self.obj = sample_list
        self.isUser = ''
        self.runLength = ''
        self.totalreads = 0
        self.id = ''
        self.barcodeList = []    
        self.cal_totalreads()
 # function that seperate samples by isUser property, return dictionary with key user/IGO, value group_samples           
    def group_by_isUser(self):
        group_list_temp = {'IsUser':[], 'IsIGO':[]}
        for i in self.obj:            
            if i.isUser == True:
                group_list_temp['IsUser'].append(i)        
            else:
                group_list_temp['IsIGO'].append(i) 
        
        group_list = {}
        for key in group_list_temp:
            group_list[key] = grouped_samples(group_list_temp[key])
            if group_list_temp[key][0].isUser == True:
                group_list[key].isUser = True
            else:
                group_list[key].isUser = False
        return group_list

    # function that seperate samples by run length, return dictionary with key run length, value group_samples  
    def group_by_runLength(self):
        # get run length list
        runLengthList = []        
        for i in self.obj:
            if i.infor['RunLength'] not in runLengthList:
                runLengthList.append(i.infor['RunLength'])
                
        group_list_temp = {}
        for n in runLengthList:
            group_list_temp[n] = []
            
        for i in self.obj:            
            group_list_temp[i.infor['RunLength']].append(i)
        
        group_list = {}
        for key in group_list_temp:
            group_list[key] = grouped_samples(group_list_temp[key])
            group_list[key].isUser = self.isUser
            group_list[key].runLength = group_list_temp[key][0].infor['RunLength']
            
        return group_list
    
    # function that seperate samples by project ID, a list of grouped_samples 
    # if wes or capture group by it self         
    def group_by_project(self):
         # get project list
        projectList = []        
        for i in self.obj:
            if (i.infor['RequestID'] not in projectList) and i.isPool == False and ("Exome" not in i.infor['Recipe']):
                projectList.append(i.infor['RequestID'])
        
        group_list_temp = {}
        if projectList:
            for n in projectList:
                group_list_temp[n] = []
            
        for i in self.obj:            
            if i.infor['RequestID'] in projectList:
                group_list_temp[i.infor['RequestID']].append(i)
            else:
                group_list_temp[i.id] = [i]
        
        group_list = []
        for key in group_list_temp:
            group_list.append(grouped_samples(group_list_temp[key]))
            group_list[-1].isUser = self.isUser
            group_list[-1].runLength = self.runLength
            group_list[-1].id = key
            for i in group_list_temp[key]:
                group_list[-1].barcodeList.append(i.infor['BarcodeSequence'])
                
        return group_list
    
    
# get sample infor from excel file and create list of samples
def ImportFromExcel(file_address):
    sample_list = pd.read_excel(file_address, index_col = 1)
    sample_name = sample_list.index
    values = sample_list.values.tolist()
    marker = [x.replace(" ","") for x in list(sample_list.columns.values)]
    
    sample_list_1 = []
    for idx, val in enumerate(sample_name):
            val = sample(val)
            val.addInfor(marker, values[idx])
            sample_list_1.append(val)
            
    return grouped_samples(sample_list_1)
    
#class lane:   
#    
#class run: