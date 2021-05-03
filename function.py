#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 12:55:38 2021

@author: luc
"""

'''
funciton list:
    index_collision: method to compare two input barcode by given length to compare and number of mismatch allowed
    index_collision_exact: method to get back of collision list from list of barcodes. -need improve later
    index_collision_mismatch: under construction
    group_by_barcode: return number of groups based on {sampleL:barcode} list given (under construction)
    optimize_runs: method to get back of optimized run type and sample in it by given reads number -need to add barcode into it later
    optimize_lanes: under construction

'''

import iteration_utilities 
from collections import OrderedDict

flow_cells = OrderedDict({"S4": [9000, 10000, 1000, 4],"S2":[3600, 3800, 200, 2],"S1":[1600, 1800, 200, 2],"SP":[700, 800, 100, 2]})
flow_cell_lanes = OrderedDict({"S4_lane": [2400, 2600, 200],"S2_lane":[1800, 1900, 100],"S1_lane":[800, 900, 100],"SP_lane":[350, 400, 50]})

def index_collision(seq1, seq2, length, num_mis):
    seq1_frag = seq1[0:length]
    seq2_frag = seq2[0:length]
    mismatch1 = 0
    mismatch2 = 0
    if '-' in seq1_frag:
        for i in range(seq1_frag.index('-')):
            if seq1_frag[i] != seq2_frag[i]:
                mismatch1 += 1
        for i in range(seq1_frag.index('-')+1, length):
            if seq1_frag[i] != seq2_frag[i]:
                mismatch2 += 1
    else:
        for i in range(length):
            if seq1_frag[i] != seq2_frag[i]:
                mismatch1 += 1
                
    if mismatch1 <= num_mis and mismatch2 <= num_mis:
        return True
    else:
        return False



# function for checking index collision with exact match
# Input: list of index(NNNNNNNN-NNNNNNNN)(list of sample class later, sample barcode is list)
# Return: a dictionary with dup barcode seq as key and barcode index in the list as value(may change to sample name later)

def index_collision_exact(index_list):
    length_list = []
    minlen = 30
    dup_list = []
    collision_list = {}
    
    # get minimum length of the index_list and list of index length
    for i in range(len(index_list)):
        if len(index_list[i]) not in length_list:
            length_list.append(len(index_list[i]))
        if len(index_list[i]) < minlen:
            minlen = len(index_list[i])
    
    # check barcode collision depending on whether they are same length        
    if len(length_list) == 1: # all same length
        dup_list = list(iteration_utilities.unique_everseen(iteration_utilities.duplicates(index_list)))
        for i in dup_list:
            collision_list[i] = []
            for j in range(len(index_list)):
                if i == index_list[j]:
                    collision_list[i].append(j)
                                    
    else: # different length
        index_list_spliced = []
        for i in index_list:
            index_list_spliced.append(i[0:minlen])
        dup_list = list(iteration_utilities.unique_everseen(iteration_utilities.duplicates(index_list_spliced)))
        for i in dup_list:
            collision_list[i] = []
            for j in range(len(index_list_spliced)):
                if i == index_list_spliced[j]:
                    collision_list[i].append(j)
             
    return collision_list


# function for optimizing runs based on reads number(algorithm comes from bin packing problem)
# will check barcode collision later
# Input: list of grouped_samples class
# Return: a list of list with run name (plus "remaining"), name of grouped sample

def optimize_runs(list_of_grouped_samples):
    run = [['remaining']]
    # get id and total reads infor from input and save as a list
    new_list = []
    for i in list_of_grouped_samples:
        new_list.append([i.id,i.totalreads])
    
    # function that can calculate total reads from the sorted_list
    def get_total_reads(a_list):
        total = 0
        for i in a_list:
            total += i[1]
            
        return total
    
    # sort the list by totalreads and calculater total reads
    sorted_list = sorted(new_list, key = lambda x: x[1], reverse = True)
    total_reads = get_total_reads(sorted_list)
    
    # fill in run and remaining list until sorted_list is empty
    temp = []
    while total_reads > 0:

        # check whether the toal reads fit in flow cell directly
        for key in flow_cells:
            if total_reads >= flow_cells[key][0] and total_reads <= flow_cells[key][1]:
                run.append([key])
                run[-1].extend(sorted_list)
                sorted_list = []
                if temp:
                        sorted_list = sorted_list + temp
                        sorted_list = sorted(sorted_list, key = lambda x: x[1], reverse = True)
                        temp = []
                        total_reads = get_total_reads(sorted_list)
                break

        # check if total_reads < 700, then put everthing into remaining
        if total_reads < 700:
            run[0].extend(sorted_list)
            break

        if sorted_list:  # if the total reads don't fit any flow cell directly
            for key in flow_cells: # decide the biggest flow cell can be use first
                if total_reads > flow_cells[key][1]:
                    run.append([key])
                    remaining = flow_cells[key][1]
                    gap = flow_cells[key][2]
                    break
            for i in sorted_list: # check samples one by one
                if remaining - i[1] > 0:
                    run[-1].append(i)
                    remaining -= i[1]
                                      
                elif remaining - i[1] == 0:
                    run[-1].append(i)
                    remaining -= i[1]
                    for i in run[-1]: # remove the item from sorted_list if already append to the run
                        if i in sorted_list:
                            sorted_list.remove(i)                    
                
                    if temp:
                        sorted_list = sorted_list + temp
                        sorted_list = sorted(sorted_list, key = lambda x: x[1], reverse = True)
                        temp = []
                        
                    break
                
                else:
                    if i == sorted_list[0]:
                        run[0].append(i)
                        sorted_list.remove(i)
                        break
                    else:
                        continue
            # after checking all the samples in the list
            
            if remaining != 0 and remaining <= gap: # if fits a run
                for i in run[-1]: # remove the item from sorted_list if already append to the run
                        if i in sorted_list:
                            sorted_list.remove(i)
                if temp:
                        sorted_list = sorted_list + temp
                        sorted_list = sorted(sorted_list, key = lambda x: x[1], reverse = True)
                        temp = []

            if remaining == flow_cells[run[-1][0]][1]: # if run is empty
                run.pop()
            elif remaining > gap: # if doesn't fit a run, remove the last sample in the potential run and try again with left
                if len(run[-1]) == 2:
                    run[0].append(run[-1][1])
                    sorted_list.remove(run[-1][1])
                    run.pop()
                    if temp:
                        sorted_list = sorted_list + temp
                        sorted_list = sorted(sorted_list, key = lambda x: x[1], reverse = True)
                        temp = []                        
                else:
                    if temp:
                        to_be_remove = []
                        for i in temp:                           
                            if run[-1][-1] > i:
                                sorted_list.append(i)
                                to_be_remove.append(i)
                        for j in to_be_remove:
                            temp.remove(j)
                        temp.append(run[-1][-1])                    
                        sorted_list.remove(run[-1][-1])
                        sorted_list = sorted(sorted_list, key = lambda x: x[1], reverse = True)
                        run.pop()
                    else:
                        temp.append(run[-1][-1])                    
                        sorted_list.remove(run[-1][-1])
                        run.pop()
                        
        total_reads = get_total_reads(sorted_list)
        
    return run




