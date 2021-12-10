#!/usr/bin/env python
import sys, os, subprocess, linecache

def intersection_point_fixy(line1, y):
    xdiff = line1[0][0] - line1[1][0]
    ydiff = line1[0][1] - line1[1][1]
    if xdiff == 0 :
        x = line1[0][0]
    else:
        slope = ydiff / xdiff
        intercept = line1[0][1] - slope*line1[0][0]
        x = (y-intercept)/slope
    return x

def intersection_point_fixx(line1, x):
    xdiff = line1[0][0] - line1[1][0]
    ydiff = line1[0][1] - line1[1][1]
    slope = ydiff / xdiff
    intercept = line1[0][1] - slope*line1[0][0]
    y = slope*x + intercept
    return y

def behavior(filname,flats,flatn,flonl,flonr):
    '''
    read total cyclone number in ff_250_1980-2020_2_3045-6060
    then judge whether the cyclone is NTN, STN or LYS in this region
    if not, it is considered to be able to pass over the TP

    20211014
    renql
    '''
    behv = ["NTN" ,"STN" ,"PAS" ,"LYS"]#,"DIF"]
    
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    number = [int(term[0]),]
   
    # create output file and tid lists for different cyclones
    outfile=[]
    for nb in range(0,len(behv),1):
        outfile.append(open(filname+"_"+behv[nb],"w"))
        outfile[nb].write(line1)
        outfile[nb].write(line2)
        outfile[nb].write(line3)
        outfile[nb].write(line4)
        locals()['tid_%s'%behv[nb]]=[] 

    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            data=[]
            value=[]
            for nl in range(0,num,1):
                line = ff.readline()
                value.append(line)
                data.append(list(map(float, line.strip().split(" "))))
           
            signal=-10
            for nl in range(0,num-1,1):
                if data[nl][2]<=flats[0] and data[nl+1][2]>flats[0]:
                    point = intersection_point_fixy([data[nl][1:3], data[nl+1][1:3]], flats[0])
                    if point<=flonr[0] and point>=flonl[0] :
                        signal = 0 # NTN
                        break 
                if signal==-10 and data[nl][2]>=flats[1] and data[nl+1][2]<flats[1]:
                    point = intersection_point_fixy([data[nl][1:3], data[nl+1][1:3]], flats[1])
                    if point<=flonr[1] and point>=flonl[1] :
                        signal = 1 # STN
                        break 
                #if signal==-10 and data[nl][1]>=flonl[-1] and data[nl+1][2]<flonl[-1]:
                #    point = intersection_point_fixx([data[nl][1:3], data[nl+1][1:3]], flonl[-1])
                #    if point<=flatn[-1] and point>=flats[-1] :
                #        signal = -1 # RET 
                #        break 
            if signal == -10 :
                #if data[-1][1]<=flonr[-1] and data[-1][1]>=flonl[-1] \
                #and data[-1][2]<=flatn[-1] and data[-1][2]>=flats[-1] :
                if data[-1][1]>=flonr[2]: 
                    signal = 2 # PAS
                else:
                    signal = -1 # LYS
            
            locals().get('tid_%s'%behv[signal]).append(term[2])
            outfile[signal].write(lineid)
            outfile[signal].write(linenum)
            for nll in value:
                outfile[signal].write(nll)

        line = ff.readline()

    ff.close()
    for nb in range(0,len(behv),1):
        number.append(len(locals().get('tid_%s'%behv[nb])))
        outfile[nb].seek(0,0) # Go back to the beginning of the file
        outfile[nb].write(line1)
        outfile[nb].write(line2)
        outfile[nb].write(line3)
        outfile[nb].write("TRACK_NUM %9d ADD_FLD    0   0 &" %number[nb+1])

        print("filter cyclone number in %s : %d" %(outfile[nb].name,number[nb+1]))
        outfile[nb].close()
    
    return number

def line_filt(filname,flats,flatn,flonl,flonr,option=2,time=24):
    '''
    read total cyclone number in ff_250_1980-2020
    then chose the cyclone that can pass through the line 

    20211014
    renql
    '''
    suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)

    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    outfile = open(filname+"_"+suffix,"w")
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write(line4)

    tid=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            data=[]
            value=[]
            for nl in range(0,num,1):
                line = ff.readline()
                value.append(line)
                data.append(list(map(float, line.strip().split(" "))))
           
            signal = 0
            if num >= time :
                for nl in range(0,num-1,1):
                    if flonl == flonr :
                        if data[nl][1] <= flonl and data[nl+1][1] > flonl and data[nl+1][1] <= 180 :
                            point = intersection_point_fixx([data[nl][1:3], data[nl+1][1:3]], flonl)
                            if point <= flatn and point >= flats :
                                signal = 1
                                break 
                    if flats == flatn :
                        #if (data[nl][2]-flats)*(data[nl+1][2]-flats) <= 0 :
                        if data[nl][2] >= flatn and data[nl+1][2] < flatn:
                            point = intersection_point_fixy([data[nl][1:3], data[nl+1][1:3]], flats)
                            if point <= flonr and point >= flonl :
                                signal = 1
                                break 
            
            if signal == 1:
                tid.append(term[2])
                outfile.write(lineid)
                outfile.write(linenum)
                for nll in value:
                    outfile.write(nll)

        line = ff.readline()

    ff.close()
    outfile.seek(0,0) # Go back to the beginning of the file
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write("TRACK_NUM %9d ADD_FLD    0   0 &" %len(tid))
    print("filter cyclone number in %s : %d" %(outfile.name,len(tid)))
    outfile.close()

    return len(tid)

def time_filt(filname,time):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
   
    fname = filname.split("_",1)
    outfile = open("%s%d_%s"%(fname[0],time,fname[1]),"w")
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write(line4)

    tid=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            if num >= time :
                tid.append(term[2])
                outfile.write(lineid)
                outfile.write(linenum)
                for nl in range(0,num,1):
                    line = ff.readline()
                    outfile.write(line)

        line = ff.readline()

    ff.close()
    outfile.seek(0,0) # Go back to the beginning of the file
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write("TRACK_NUM %9d ADD_FLD    0   0 &" %len(tid))
    print("%s liftime longer than %d hours: %d" %(outfile.name,time,len(tid)))
    outfile.close()

    return len(tid)

