#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 12:55:56 2021

@author: luc
"""
import function
import components
import pandas as pd



#test = components.ImportFromExcel("RunPlanner-04-22-2021.xlsx")
#test1 = test.group_by_isUser()
#
#
#test2 = test1['IsUser'].group_by_runLength()
#
#
#test3 = test2['PE50'].group_by_project()
#
##
##for i in test3:
##    print(i.id, i.barcodeCollision, i.barcodeList)  
#
#test4 = function.group_by_barcode(test3)
##for i in test4:
##    print(i.id, i.barcodeList) 
#print (test4)
#
#
test = components.ImportFromExcel("RunPlanner-05-10-2021.xlsx")
test1 = test.group_by_isUser()
for key in test1:
    if key != 'nan':
        test2 = test1[key].group_by_runLength()
    for key1 in test2:
        if key1 != 'nan':
            test3= test2[key1].group_by_project()
            test4 = function.optimize_runs(test3)
            for i in test4:
                if len(i) > 1:
                    function.WriteToExcel(i[1:],key + key1 + i[0]+'.xlsx')

#test_group = function.group_lanes(test4[1][1:],2)
print(function.group_by_barcode(test3))
#result = []
#for i in test_group:
#    temp = []
#    for l in i:
#        temp = temp + l.obj
#    result.append(components.grouped_samples(temp))
#
#function.WriteToExcel(result, 'result.xlsx')