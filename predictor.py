# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 22:25:40 2018

@author: hasee
"""
from collections import Counter
import datetime
import math 

def predict_vm(train, inputfile):
    temp = 0
    shitianpingjun=0
    def expsmooth3(timeseries, alph):
        temp = 0
        for i in range(3):
            temp = timeseries[i] + temp
        ori = temp / 3.0

        s1 = []
        s2 = []
        s3 = []
        s1.append(ori)
        for i in range(1, len(timeseries) + 1):
            ta = float((alph * timeseries[i - 1]) + (1 - alph) * s1[i - 1])
            s1.append(ta)

        s2.append(ori)
        for i in range(len(s1) - 1):
            s2.append(alph * s1[i + 1] + (1 - alph) * s2[i])

        s3.append(ori)
        for i in range(len(s2) - 1):
            s3.append(alph * s2[i + 1] + (1 - alph) * s3[i])

        a = 3 * s1[-1] - 3 * s2[-1] + s3[-1]
        b = (alph / (2 * (1 - alph) * (1 - alph))) * (
                (6 - 5 * alph) * s1[-1] - 2 * (5 - 4 * alph) * s2[-1] + (4 - 3 * alph) * s3[-1])
        c = (alph * alph / 2 / (1 - alph) / (1 - alph)) * (s1[-1] - 2 * s2[-1] + s3[-1])

        y1 = a + b + c
        if y1<=0:
            y1=0
        return (y1)

    def meanstd(y):
        mean = sum(y) / float(len(y))
        stdtemp = 0
        for i in range(len(y)):
            stdtemp = stdtemp + (y[i] - mean) * (y[i] - mean)
        std = math.sqrt(stdtemp / float(len(y)))
        return (mean, std)

    def expsmooth2(timeseries, alph):
        temp=0
        for i in range(3):
            temp = timeseries[i] + temp
        ori = temp / 3


        s1 = []
        s2 = []
        s1.append(ori)
        for i in range(1, len(timeseries) + 1):
            ta = float((alph * timeseries[i - 1]) + (1 - alph) * s1[i - 1])

            s1.append(ta)
        '''
        timeseries.append(s1[-1])
        ta=float((alph*timeseries[-1])+(1-alph)*s1[i-1])
        s1.append(ta)
        '''
        s2.append(s1[1])
        for i in range(1, len(s1) - 1):
            s2.append(alph * s1[i + 1] + (1 - alph) * s2[i - 1])

        a = 2 * s1[-1] - s2[-1]
        b = alph / (1 - alph) * (s1[-1] - s2[-1])
        y1 = a + b
        y2 = a + 2 * b
        y3 = a + 3 * b

        return (y1, y2, y3)








    def expsmooth2(timeseries, alph):
        temp=0
        for i in range(3):
            temp = timeseries[i] + temp
        ori = temp / 3


        s1 = []
        s2 = []
        s1.append(ori)
        for i in range(1, len(timeseries) + 1):
            ta = float((alph * timeseries[i - 1]) + (1 - alph) * s1[i - 1])

            s1.append(ta)
        '''
        timeseries.append(s1[-1])
        ta=float((alph*timeseries[-1])+(1-alph)*s1[i-1])
        s1.append(ta)
        '''
        s2.append(s1[1])
        for i in range(1, len(s1) - 1):
            s2.append(alph * s1[i + 1] + (1 - alph) * s2[i - 1])

        a = 2 * s1[-1] - s2[-1]
        b = alph / (1 - alph) * (s1[-1] - s2[-1])
        y1 = a + b
        y2 = a + 2 * b
        y3 = a + 3 * b

        return (y1, y2, y3)
    def split_fp_delet(phy, n, flag, xuqiu_vim, flavorlist, xuqiu):
        tp = sum(xuqiu_vim2.values())

        while (tp != 0):
            tp = sum(xuqiu_vim2.values())
            xuqiu_vim = dict(zip(flavorlist, xuqiu))
            for i in range(n):
                phylist[i] = []
            a = [int(phy[flag])] * n
            b = [int(phy[flag - 1])] * n
            i = 0
            enough = 1
            while ((tp != 0) & enough):

                dictzero(xuqiu_vim)
                biggest = dictmax(dictcom(xuqiu_vim, vim[2]))[0]
                biglist = dictmax(dictcom(xuqiu_vim, vim[2]))[1]
                biggest = dictmax(dictcom2(biglist, vim[flag]))[0]
                key = biggest[0]
                while (xuqiu_vim[key] != 0):
                    t = 0
                    if ((a[i] >= vim[flag][key]) & (b[i] >= vim[re(flag)][key])):
                        xuqiu_vim[key] = xuqiu_vim[key] - 1
                        a[i] = a[i] - vim[flag][key]
                        b[i] = b[i] - vim[re(flag)][key]
                        phylist[i].append(key)
                        tp = tp - 1
                        i = (i + 1) % n
                    else:
                        while ((t < n) & (tp != 0)):
                            t = t + 1
                            i = (i + 1) % n
                            if ((a[i] >= vim[flag][key]) & (b[i] >= vim[re(flag)][key])):
                                xuqiu_vim[key] = xuqiu_vim[key] - 1
                                a[i] = a[i] - vim[flag][key]
                                b[i] = b[i] - vim[re(flag)][key]
                                phylist[i].append(key)
                                i = (i + 1) % n
                                tp = tp - 1
                                dictzero(xuqiu_vim)
                                t = 0
                            if t==n:
                                xuqiu_vim[key]=xuqiu_vim[key]-1

                        if tp == 0:
                            return (a, n, phylist, b)

                        if t == n:
                            enough = 0
                            n = n + 1
                            break

        return (a, n, phylist, b)

    def split_fp(phy, n, flag, xuqiu_vim, flavorlist, xuqiu):
        tp = sum(xuqiu_vim2.values())
     ## banben1
        while (tp != 0):
            tp = sum(xuqiu_vim2.values())
            xuqiu_vim = dict(zip(flavorlist, xuqiu))
            for i in range(n):
                phylist[i] = []
            a = [int(phy[flag])] * n
            b = [int(phy[flag - 1])] * n

            enough = 1
            while ((tp != 0) & enough):
                i = 0
                dictzero(xuqiu_vim)
                biggest = dictmax(dictcom(xuqiu_vim, vim[2]))[0]
                biglist = dictmax(dictcom(xuqiu_vim, vim[2]))[1]
                biggest = dictmax(dictcom2(biglist, vim[flag]))[0]
                key = biggest[0]
                while (xuqiu_vim[key] != 0):
                    t = 0
                    if ((a[i] >= vim[flag][key]) & (b[i] >= vim[re(flag)][key])):
                        xuqiu_vim[key] = xuqiu_vim[key] - 1
                        a[i] = a[i] - vim[flag][key]
                        b[i] = b[i] - vim[re(flag)][key]
                        phylist[i].append(key)
                        tp = tp - 1
                        i = (i + 1) % n
                    else:
                        while ((t < n) & (tp != 0)&(xuqiu_vim[key] != 0)):
                            t = t + 1
                            i = (i + 1) % n
                            if ((a[i] >= vim[flag][key]) & (b[i] >= vim[re(flag)][key])):
                                xuqiu_vim[key] = xuqiu_vim[key] - 1
                                a[i] = a[i] - vim[flag][key]
                                b[i] = b[i] - vim[re(flag)][key]
                                phylist[i].append(key)
                                i = (i + 1) % n
                                tp = tp - 1
                                t = 0
                                break
                        if tp == 0:
                            return (a, n, phylist, b)

                        if t == n:
                            enough = 0
                            n = n + 1
                            break

        return (a, n, phylist, b)


    def split_fp(phy, n, flag, xuqiu_vim, flavorlist, xuqiu):
        tp = sum(xuqiu_vim2.values())
     ## banben1
        while (tp != 0):
            tp = sum(xuqiu_vim2.values())
            xuqiu_vim = dict(zip(flavorlist, xuqiu))
            for i in range(n):
                phylist[i] = []
            a = [int(phy[flag])] * n
            b = [int(phy[flag - 1])] * n

            enough = 1
            while ((tp != 0) & enough):
                i = 0
                dictzero(xuqiu_vim)
                biggest = dictmax(dictcom(xuqiu_vim, vim[2]))[0]
                biglist = dictmax(dictcom(xuqiu_vim, vim[2]))[1]
                biggest = dictmax(dictcom2(biglist, vim[flag]))[0]
                key = biggest[0]
                while (xuqiu_vim[key] != 0):
                    t = 0
                    if ((a[i] >= vim[flag][key]) & (b[i] >= vim[re(flag)][key])):
                        xuqiu_vim[key] = xuqiu_vim[key] - 1
                        a[i] = a[i] - vim[flag][key]
                        b[i] = b[i] - vim[re(flag)][key]
                        phylist[i].append(key)
                        tp = tp - 1
                        i = (i + 1) % n
                    else:
                        while ((t < n) & (tp != 0)&(xuqiu_vim[key] != 0)):
                            t = t + 1
                            i = (i + 1) % n
                            if ((a[i] >= vim[flag][key]) & (b[i] >= vim[re(flag)][key])):
                                xuqiu_vim[key] = xuqiu_vim[key] - 1
                                a[i] = a[i] - vim[flag][key]
                                b[i] = b[i] - vim[re(flag)][key]
                                phylist[i].append(key)
                                i = (i + 1) % n
                                tp = tp - 1
                                t = 0
                                break
                        if tp == 0:
                            return (a, n, phylist, b)

                        if t == n:
                            enough = 0
                            n = n + 1
                            break

        return (a, n, phylist, b)












    def get_keys(d, value):
        return [k for k,v in d.items() if v == value]
    
    def dictmin(dic):
        minnest=list(min(dic.items(), key=lambda x: x[1])) 
        
        return minnest
    
    def dictmax(dic):
        biggest=list(max(dic.items(), key=lambda x: x[1])) 
        k=get_keys(dic,biggest[1])
        return [biggest,k]
    def dictzero(dic):
        l=[]
        for key in dic:
            if dic[key]==0:
                l.append(key)
        for i in range(len(l)):
            del dic[l[i]]
        return dic
    
    def linearreg(train_x,train_y):
        x_mean=float(sum(train_x))/len(train_x)
        y_mean=float(sum(train_y))/len(train_y)
        m=0
        n=0
        for i in range(len(train_x)):
            m=m+(train_x[i]-x_mean)*(train_y[i]-y_mean)
            n=n+(train_x[i]-x_mean)*(train_x[i]-x_mean)
        a=m/n
        b=y_mean-a*x_mean
        return(a,b)
    def dictcom2(dic1,dic2):
        dic3={}
        for i in range(len(dic1)):
            dic3[dic1[i]]=dic2[dic1[i]]
            
        return dic3
    
    
    
    
    def datedif(a,b):
        start=datetime.datetime.strptime(str(a),"%Y-%m-%d")
        end=datetime.datetime.strptime(str(b),"%Y-%m-%d")
        cha=(end-start).days+1
        return cha
    
    
    def getresult(w,trl,tel):
        sumsult=0
        for i in range(trl,trl+tel):
            sumsult=sumsult+i*w[0]+w[1]
        return sumsult
    
    def dateRange(start, end, step=1, format="%Y-%m-%d"):
        strptime, strftime = datetime.datetime.strptime, datetime.datetime.strftime
        days = (strptime(end, format) - strptime(start, format)).days
        return [strftime(strptime(start, format) + datetime.timedelta(i), format) for i in xrange(0, days, step)]
    
    def dictcom(dic1,dic2):
        dic3={}
        for key in dic1.keys():
            dic3[key]=dic2[key]
        return dic3
    
    def re(a):
        if a==0:
            return 1
        if a==1:
            return 0
        
        
    
    def greedysplit(phy,vim,flag,xuqiu_vim,n):
        
        phylist={}
        fflag=0
        tp=sum(xuqiu_vim.values())
        while(fflag==0):
            
            for i in range(n):
                phylist[i]=[]
            a=[int(phy[flag])]*sum(xuqiu_vim.values())
            b=[int(phy[flag-1])]*sum(xuqiu_vim.values())
            i=0
            while (tp!=0):
                dictzero(xuqiu_vim)
                biggest=dictmax(dictcom(xuqiu_vim,vim[2]))[0]  
                biglist=dictmax(dictcom(xuqiu_vim,vim[2]))[1] 
                biggest=dictmax(dictcom2(biglist,vim[flag]))[0]
                key=biggest[0]
                while(xuqiu_vim[key]!=0):     
                    if((a[i]>vim[flag][key])&(b[i]>vim[re(flag)][key])):
                        xuqiu_vim[key]=xuqiu_vim[key]-1
                        a[i]=a[i]-vim[flag][key]
                        b[i]=b[i]-vim[re(flag)][key]
                        tp=tp-1
                        i=(i+1)%n
                    else:
                        t=1
                        while(t<n&(a[i+t]<vim[flag][key])&(b[i+t]<vim[re(flag)][key])):
                            t=t+1                     
                        if t<n:
                            xuqiu_vim[key]=xuqiu_vim[key]-1
                            a[i+t]=a[i+t]-vim[flag][key]
                            b[i+t]=b[i+t]-vim[re(flag)][key]
                            tp=tp-1
                            i=(i+1)%n
                        else:
                            n=n+1
                            fflag=1
                            
        return (a,b)
                   
                   
            
        
        
        
        
        
        
        
    def greedy2(phy,vim,flag,xuqiu_vim): 
        number=0
        phylist={}
        tp=sum(xuqiu_vim.values())
        for i in range(tp):
            phylist[i]=[]
        
        a=[int(phy[flag])]*sum(xuqiu_vim.values())
        b=[int(phy[flag-1])]*sum(xuqiu_vim.values())
        for i in range(sum(xuqiu_vim.values())):
            
            n=0
            while n<len(a):
                dictzero(xuqiu_vim)
                biggest=dictmax(dictcom(xuqiu_vim,vim[flag]))[0]  
                biglist=dictmax(dictcom(xuqiu_vim,vim[flag]))[1] 
                biggest=dictmax(dictcom2(biglist,vim[re(flag)]))[0]
                
                if xuqiu_vim[biggest[0]]>0: 
                    if (a[n]>=vim[flag][biggest[0]])&(b[n]>=vim[re(flag)][biggest[0]]):
                        a[n]=a[n]-vim[flag][biggest[0]]
                        b[n]=b[n]-vim[re(flag)][biggest[0]]
                        xuqiu_vim[biggest[0]]=xuqiu_vim[biggest[0]]-1
                        phylist[n].append(biggest[0])
                        break
                    else:
                        n=n+1
        for i in range(len(a)):
            if a[i]-int(phy[flag])<0:
                number=number+1
        
        return (a,number,phylist,b)    
    
    

    data={}
    temp=[]
    
    
    cpu=[]
    mem=[]
    flavorlist=[]
    i=3
    while inputfile[i]!='\r\n':
        flavorlist.append(inputfile[i].split()[0])
        cpu.append(int(inputfile[i].split()[1]))
        mem.append(int(inputfile[i].split()[2]))
        i=i+1
    
      
        
        
    name=[]
    time=[]
    
    
    traindata=[]
    while inputfile[-1]=='\r\n':
        del inputfile[-1]
    testdaterange=(datedif(inputfile[-2].split()[0],inputfile[-1].split()[0])-1 ) 
    start=train[0].split()[2]
    end=train[-1].split()[2]







    dayrange=dateRange(str(train[0].split()[2]),str(train[-1].split()[2]))
    dayrange.append(train[-1].split()[2])

    traindaterange = 42
    alph2 = 0.38
    
    for index, item in enumerate(train):
        values = item.split()
        uuid = values[0]
        flavorName = values[1]
        createTime = values[2]
        if flavorName in flavorlist:
            name.append(flavorName)
            time.append(createTime)
    for j in range(len(flavorlist)):
        for i in range(len(name)):
            if name[i]==flavorlist[j]:
                temp.append(time[i])
                data[flavorlist[j]]=temp
        temp=[]
    xuqiu=[]

    for i in range(len(flavorlist)):

        ccc = data[flavorlist[i]]
        teeet=flavorlist[i]
        rrr = {}
        bbb = Counter(ccc)
        for item in dayrange:
            if item in ccc:
                rrr[item] = bbb[item]
            else:
                rrr[item] = 0

        y = []
        tk = traindaterange / testdaterange
        trk=traindaterange/tk*tk
        #整周期个天数
        for i in range(len(dayrange) - trk, len(dayrange)):
            y.append(rrr[dayrange[i]])


        tendays=0
        for i in range(1,11):
            tendays=tendays+y[-i]
        tendays=(tendays/10.0)

        noise=meanstd(y)

        week = []



        for m in range(tk):
            temp = 0
            for j in range(testdaterange):
                temp = temp + y[j + m * testdaterange]
            week.append(temp)




        def doubleweek(week,a1,a2):
            return (a1*week[-1]+a2*week[-2])
        '''
        temp=0
        for i in range(1,11):
            temp=y[-i]+temp
        
        for v in week:
            if ((v<(noise[0]-4*noise[1]))|(v>noise[0]+4*noise[1])):
                week.remove(v)
        '''
        doubleweekf=0
        if(doubleweekf):
            if testdaterange==14:
                xuqiu.append(int(math.ceil(doubleweek(week,0.6,0.5))*2))
            else:
                xuqiu.append(int(math.ceil(doubleweek(week, 0.6, 0.4))))
        if(shitianpingjun):
            xuqiu.append(int(math.ceil(tendays*testdaterange)))
        else:
            xuqiu.append(int(expsmooth3(week, alph2)))
    '''
    for i in range(len(flavorlist)):    
        ccc=data[flavorlist[i]]
        rrr={}
        bbb=Counter(ccc)
        for item in dayrange:
            if item  in ccc:
                rrr[item]=bbb[item]
            else:
                rrr[item]=0
        x=range(traindaterange)
        y=rrr.values()
        w=linearreg(x,y)
        xuqiu.append(int(round(getresult(w,traindaterange,testdaterange))))
        
        #plt.plot(x,y)   
        #print getresult(w,cha,testdaterange)
    '''      
    if inputfile[-4]=='CPU\r\n':
        flag=0
    else:
        flag=1


    
    vim=[]
    phcom=[]
    phcom.append(inputfile[0].split()[0])  
    phcom.append(int((inputfile[0].split()[1]))*1024)   
    
    
    vim.append(dict(zip(flavorlist,cpu)))
    vim.append(dict(zip(flavorlist,mem)))
    bili={}
    rebili={}
    for value in flavorlist:
        bili[value]=(vim[1][value]/1024)/(vim[0][value])
        rebili[value]=(vim[0][value])/float((vim[1][value]/1024))
    if(flag==0):
        vim.append(bili)
    else:
        vim.append(rebili)
    xuqiu_vim=dict(zip(flavorlist,xuqiu))

    xuqiu_vim2=dict(zip(flavorlist,xuqiu))
    
    
    
    
    
    dictzero(xuqiu_vim)
    
    
    
    

    

    cpusum=0
    memsum=0
    for va in flavorlist:
        cpusum=cpusum+xuqiu_vim2[va]*vim[0][va]
        memsum=memsum+xuqiu_vim2[va]*vim[1][va]
    cpumaxnum=math.ceil(float(cpusum)/float(phcom[0]))
    memmaxnum=math.ceil(float(memsum)/float(phcom[1]))
    phynumber=int(max(cpumaxnum,memmaxnum))
    
    
    
    
    phy=phcom
    
    phylist={}
    '''
    这段是用来区分cpu mem
    if(flag==0):
        rrrr = split_fp(phcom, phynumber, flag, xuqiu_vim, flavorlist, xuqiu)
   

    if(flag==1):
        rrrr = greedy2(phcom, vim, flag, xuqiu_vim)
        #rrrr = split_fp(phcom, phynumber, flag, xuqiu_vim, flavorlist, xuqiu)
     '''

    #flag改为1的时候是测CPU
    '''
    if(flag==0):
        return 0
    '''

    if (testdaterange==8):
        rrrr = greedy2(phcom, vim, flag, xuqiu_vim)

    if(testdaterange==7):
        rrrr = split_fp(phcom, phynumber, flag, xuqiu_vim, flavorlist, xuqiu)
    else:
        return 0
    #rrrr=split_fp(phcom,phynumber,flag,xuqiu_vim,flavorlist,xuqiu)
    #rrrr = greedy2(phcom, vim, flag, xuqiu_vim)

    '''
    #这段是填满 使用最小的填满
    minre=dictmin(vim[flag])
    minreother=vim[re(flag)][minre[0]]
    n=0
    while n!=rrrr[1]:
        
        if (rrrr[0][n]>=minre[1])&(rrrr[-1][n]>=minreother):
            rrrr[0][n]=rrrr[0][n]-minre[1]
            rrrr[-1][n]=rrrr[-1][n]-minreother
            rrrr[2][n].append(minre[0])
            xuqiu_vim2[minre[0]]=xuqiu_vim2[minre[0]]+1
        else:
            n=n+1
    
    '''






    outputresult=[]
    outputresult.append(sum(xuqiu_vim2.values()))
    for item in flavorlist:
        temp=str(item+' '+str(xuqiu_vim2[item]))
        outputresult.append(temp)
    outputresult.append('') 
    outputresult.append(rrrr[1]) 
    for i in range(rrrr[1]):
        ct=Counter(rrrr[2][i])
        j=0
        x=str(ct.keys()[j])+' '+str(ct.values()[j])
        if len(ct)>1:
            for j in range(1,len(ct)):                
                x=x+' '+str(ct.keys()[j])+' '+str(ct.values()[j])
        outputresult.append(str(i+1)+' '+x) 
    return outputresult