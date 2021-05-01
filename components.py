# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 12:55:23 2021

@author: luc

"""

'''
class list:
    sample: single sampels with information from run planner website
            important infor: id, pool id, run length, request name, recipe, reads num
            with function: create from excel
    grouped_samples: sample list with information of id, totalreads, is_project, is_user_project, barcode list
            function: create from list of samples
            
    ?lane: one lane of a run, have four types, with infor of list of samples
            with function is_valid
    ?run: one run

'''

# sample class will contain all the information we can get from run planner
class sample:
    def __init__(self, id, runlength, reads, recipe):
        self.id = id
        self.runlength = runlength
        self.reads = int(reads)
        self.recipe = recipe

    
# grouped_sample will be a list of samples that belongs to one project.
# will have additional character totalreads that sums up all the reads numebr from each sample
class grouped_samples(sample):   
    def __init__(self, id, totalreads):
        self.id = id
        self.totalreads = int(totalreads)
    
    
#class lane:   
#    
#class run: