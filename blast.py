#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 14:34:57 2018

@author: hannes
"""


import sys
import getopt


#instructions
form='blast.py\t\t\t\t\tVersion 1.2\t23/04/2019\n\t--f1 <fasta/flat>\n\t--f2 <fasta/flat> *opt*\n\t--first <number> *opt, only take first n sequences*\n\t--contrast *opt,reduces markersize*\n\t--equal *opt,axis ratios equal*\n\t-p image_name *opt,default=noimage*\n\t-w width *opt,sliding window size, default=20*\n\t-l *opt,create txt with axis legend*'

def slyce(line, width):
    d = {}
    
    for i in range(len(line)+1-width):
        seq_short = line[i:i+width:]
        if seq_short in d:
            d[seq_short].append(i)
        else:
            d[seq_short] = [i]
    return d



def match(line1, line2, width):
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    l1 = []
    l2 = []
    d1 = slyce(line1, width)
    d2 = slyce(line2, width)


    for i in d1:
        seq = Seq(i, generic_dna)
        
        #i in d1 and in d2??
        if i in d2:         
            for j in range(len(d1[i])):
                for k in range(len(d2[i])):
                    l1.append(d1[i][j])
                    l2.append(d2[i][k])

        #complementary?
        if str(seq.complement()) in d2:
            for j in range(len(d1[i])):
                for k in range(len(d2[str(seq.complement())])):
                    l1.append(d1[i][j])
                    l2.append(d2[str(seq.complement())][k])

        #it is in rev list?
        if i[::-1] in d2:
            for j in range(len(d1[i])):
                for k in range(len(d2[i[::-1]])):
                    l1.append(d1[i][j])
                    l2.append(d2[i[::-1]][k])

        #reverse complementary?
        if str(seq[::-1].complement()) in d2:
            for j in range(len(d1[i])):
                for k in range(len(d2[str(seq[::-1].complement())])):
                    l1.append(d1[i][j])
                    l2.append(d2[str(seq[::-1].complement())][k])

    return l1,l2




def plot(x, y, lablex, labley, width, marker_size, figure, contrast):
    #N's should not give hits
    x = x.replace('n' and 'N', '*')
    y = y.replace('n' and 'N', '#')
    x,y = match(x,y,width)
    figure.plot(x,y, 'k.', ms=marker_size, alpha = contrast)
    figure.set_xlabel(lablex)
    figure.set_ylabel(labley)
    figure.set_title('BLAST {x22} against {y22}'.format(x22 = lablex, y22 = labley))           

def plot_no_title(x, y, lablex, labley, width, marker_size, figure, contrast):
    x = x.replace('n' and 'N', '*')
    y = y.replace('n' and 'N', '#')    
    x,y = match(x,y,width)
    figure.plot(x,y, 'k.', ms=marker_size, alpha = contrast)
    #ax1.set_xlabel(lablex)
    #ax1.set_ylabel(labley)
    #ax1.set_title('BLAST {x22} against {y22}'.format(x22 = lablex, y22 = labley))       
    



def combinations_of_two(lst):
    
    c = []
    for i in range(len(lst)):
        for j in range(len(lst)-1-i):
            if lst[i] != lst[len(lst)-1-j]:
                c.append([lst[i], lst[len(lst)-1-j]])
    return c

def combinations_of_two_incl_self(lst):
    
    c = []
    for i in range(len(lst)):
        for j in range(len(lst)-1-i):
            if [lst[i], lst[len(lst)-1-j]] not in c:
                c.append([lst[i], lst[len(lst)-1-j]])
        c.append([lst[i], lst[i]])
    return sorted(c)



def plot_update(l1, l2, figure1, width, marker_size, contrast, pad, f1, f2, init, ax_arr):

    import numpy as np

    #two single fastas
    if len(l1) <= 1 and len(l2) <= 1:

        #screen info
        info = np.array([[str(0)]*3], dtype=object)
        info[:,1] = ['against']
        for i in range(len(l1)):     
            info[i,0] = l1[i][0]
            info[i,2] = l2[i][0]

        #print(info[(i+1)*(j+1)-1,:])
        print(info[0,0], '\t', info[0,1], '\t', info[0,2])
        seq1 = l1[0][1]
        seq2 = l2[0][1]
        for ax in ax_arr.flatten():
            plot(seq1, seq2, l1[0][0], l2[0][0], width, marker_size, ax, contrast)
            ax.set_xlim((0,len(seq1)))
            ax.set_ylim((0,len(seq2)))



    #flat and fasta
    elif (len(l1) >= 1 and len (l2) == 1) or (len(l1) == 1 and len(l2) >= 1):
        
        if len(l1) > 1:
            temp = l1
            l1 = l2
            l2 = temp
            temp = ''
        l2 = l2[::-1]
        #contig index
        print('contigs:')
        for i in range(len(l2)):
            print(i, l2[len(l2)-i-1][0])
        print()
        print()
        print(len(l2),'plots:')
        print()
        #screen info
        info = np.array(len(l1)*len(l2)*[[str(0)]*3], dtype=object)
        info[:,1] = ['against']
        #print(info)
        
        for i in range(len(l1)):
            for j in range(len(l2)):     
                info[(i+1)*(j+1)-1,0] = l1[i][0]
                info[(i+1)*(j+1)-1,2] = l2[j][0]

        #name the rows
        rows = ['{}'.format(row) for row in range(len(l2))]
        i=0
        j=-1
        temp = 0
        for ax in ax_arr.flatten():
            if j == len(l2)-1:
                i+=1
                j=-1
            temp+=1
            j+=1

            #print(info[(i+1)*(j+1)-1,:])
            print(info[(i+1)*(j+1)-1,0], '\t', info[(i+1)*(j+1)-1,1], '\t', info[(i+1)*(j+1)-1,2])
            seq1 = l1[i][1]
            seq2 = l2[j][1]
            print(str((i+1)*(j+1)))

            #label contig rows by numbers
            ax.annotate(rows[len(l2)-j-1], xy=(1, 0.5), xytext=(pad, 0),
                        xycoords='axes fraction', textcoords='offset points',
                        size='large', ha='left', va='center')  

            #which axes labels & ticks to display
            if j==len(l2)-1:
                ax.set_xlabel(l1[i][0])
                
            else:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected                        
                                labelbottom=False) # labels along the bottom edge are off    
                 
            plot_no_title(seq1, seq2, l1[i][0], l2[j][0], width, marker_size, ax, contrast)
            ax.set_xlim((0,len(seq1)))
            ax.set_ylim((0,len(seq2)))

        temp=''



    #two flat files or single flat file 
    elif len(l1)>1 and len(l2):
        d_row = {}
        d_col = {}
        l_combo = []
        
        #make dictionary
        for i in range(len(l1)):            
            d_col[i] = l1[i][0]
        
        for i in range(len(l2)):
            d_row[i] = l2[i][0]
        
        #combinations
        for j in range(len(l2)):
            for i in range(len(l1)):            
                    l_combo.append([i, j])

        
        cols = ['{}'.format(col) for col in range(len(l1))]
        rows = ['{}'.format(row) for row in range(len(l2))]
        
        #order plots
        l_ordered = []
        for i in range(len(l_combo)):
            l_ordered.append( [l_combo[i][0], len(l2)-1 - l_combo[i][1]] )
        temp=''
        axes_index = l_combo
        l_combo = l_ordered
        print(axes_index)
        print(l_ordered)
        #screen info
        info = np.array(len(l1)*len(l2)*[[str(0)]*3], dtype=object)
        info[:,1] = ['against']
        for i in range(len(l_combo)):
            info[i][0] = l1[l_combo[i][0]][0]
            info[i][2] = l2[l_combo[i][1]][0]

        #for multiplots
        for i,ax in zip(range(len(l_combo)),ax_arr.flatten()):
            print(info[i])

            seq1 = l1[l_combo[i][0]][1] #x-axis
            seq2 = l2[l_combo[i][1]][1] #y-axis
            
            y=l_combo[i][1]
            x=axes_index[i][0]
            
            ax.set_xlim((0,len(seq1)))
            ax.set_ylim((0,len(seq2)))

            #label columns
            if y == len(l2)-1:
                ax.annotate(cols[x], xy=(0.5, 1), xytext=(0, pad),
                            xycoords='axes fraction', textcoords='offset points',
                            size='large', ha='center', va='baseline') 
            
            if x == len(l1)-1:
                ax.annotate(rows[y], xy=(1, 0.5), xytext=(pad, 0),
                            xycoords='axes fraction', textcoords='offset points',
                            size='large', ha='left', va='center')
            
            #which axes labels & ticks to display            
            if x == 0 and y == 0:    
                ax.tick_params( axis='both',
                                which='both')                
            
            elif y != 0 and x !=0:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom=False,      # ticks along the bottom edge are off
                                left=False,
                                top=False,         # ticks along the top edge are off
                                labelbottom=False,
                                labelleft=False) # labels along the bottom edge are off     
            elif y == 0 and x !=0:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,
                                top=False,         # ticks along the top edge are off
                                labelleft=False) # labels along the bottom edge are off   
               
            elif y !=0 and x ==0:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom=False,      # ticks along the bottom edge are off
                                
                                top=False,         # ticks along the top edge are off
                                labelbottom=False) 

            plot_no_title(seq1, seq2, d_col[l_combo[i][0]], d_row[l_combo[i][1]], width, marker_size, ax, contrast)









































def main(argv):
    contrast = 1
    pad = 20 #space between plots and numbers
    width = 20 #default min alignment
    f1 = 0
    f2 = 0
    p  = 0
    eq = 0
    log = 0
    marker_size = 1
    first = 0
    try:
        opts, args = getopt.getopt(argv,"hlvp:w:",["f1=","f2=","first=","contrast","equal","version"])
    except getopt.GetoptError:
        print ('{}'.format(form))
        sys.exit()
    for opt, arg in opts:        
        if opt == '-h' or opt == '-help' or opt == '--help' or opt == '--version' or opt == '-v':
            print ('{}'.format(form))
            sys.exit()
        elif opt == '-w':
            width = int(arg)
        elif opt == '-p':
            image_name = arg
            p = 1
        elif opt == '-l':
            log = 1
        elif opt == '--f1':
            fasta1 = arg
            f1 = 1
        elif opt == '--f2':
            fasta2 = arg
            f2 = 1
        elif opt == '--contrast':
            marker_size = 0.005
        elif opt == '--equal':
            eq = 1
        elif opt == '--first':
            first = int(arg)


    print(argv)
    print()
    print()
    
      
    import numpy as np
    import re
    import matplotlib.pyplot as pl
    from matplotlib.widgets import Slider
    
    l1=[]
    l2=[]
    t=0
    
    # open the fasta file
    if f1 == 1:
        ass=open(fasta1)
        for line in ass:
            line = line.strip('\n')
            #check if .fasta header
            match= re.search("^(>)", line)
            #first contig  
            if match and t == 0:       
                l1.append([line])           
                temp=''
                t=1
            #if not first contig
            elif match and t==1:
                l1[-1].append(temp)
                l1.append([line])
                temp=''
            else:
                temp += line
        #last contig
        l1[-1].append(temp)
        temp=''
        t=0
        #print(l1)
        ass.close()
    
    if f2 == 1:
        # open second fasta file
        ass=open(fasta2)
        for line in ass:
            line = line.strip('\n')
            #check if .fasta header
            match= re.search("^(>)", line)
            #first contig  
            if match and t == 0:       
                l2.append([line])           
                temp=''
                t=1
            #if not first contig
            elif match and t==1:
                l2[-1].append(temp)
                l2.append([line])
                temp=''
            else:
                temp += line
        #last contig
        l2[-1].append(temp)
        temp=''
        t=0
        ass.close()

    #ax1.set_subplots(figsize=(15,10))    
    if f1 == 0 and f2 == 0:
        print()
        print ('{}'.format(form))
        sys.exit()

    #only take first n sequences of plots
    if first != 0:
        l1 = l1[:first]
        l2 = l2[:first]

    #prepare plots and sequences
    if len(l1) == 1 and len(l2) == 0 or len(l1) == 0 and len(l2) == 1:
        subplots = 1
        figure1, ax_arr = pl.subplots(nrows=1, ncols=1, figsize=(20,10), dpi=100)
        if len(l2) == 0:
            l2 = l1
        else:
            l1 = l2
    elif len(l1) == 1 and len(l2) == 1:
        subplots = 1
        figure1, ax_arr = pl.subplots(nrows=1, ncols=1, figsize=(20,10), dpi=100)
    elif len(l1)>1 and f2==0 or len(l2)>1 and f1==0:
        if f2==0:
            subplots = len(l1)**2
            figure1, ax_arr = pl.subplots(nrows=len(l1), ncols=len(l1), figsize=(20,10), dpi=100)
            l2,fasta2 = l1,fasta1
        else:
            subplots = len(l2)**2
            print(subplots, len(l2))
            figure1, ax_arr = pl.subplots(nrows=len(l2), ncols=len(l2), figsize=(20,10), dpi=100)
            l1,fasta1 = l2,fasta2
    elif (len(l1) >= 1 and len (l2) == 1):
        figure1, ax_arr = pl.subplots(nrows=len(l1), ncols=len(l2), figsize=(20,10), dpi=100)
        subplots =len(l1) *len(l2)
    else:
        subplots =len(l1) *len(l2)
        figure1, ax_arr = pl.subplots(nrows=len(l2), ncols=len(l1), figsize=(20,10), dpi=100) #, sharex=True, sharey=True)

    #creat subplot index
    if subplots == 1:
        ax_arr = np.array(ax_arr)

    ax_arr_dic = {}
    temp = 0
    for ax in ax_arr.flatten():
        temp += 1
        ax_arr_dic[ax]= temp
    temp=''



    #two single fastas
    if len(l1) <= 1 and len(l2) <= 1:
        #screen info
        info = np.array([[str(0)]*3], dtype=object)
        info[:,1] = ['against']
        for i in range(len(l1)):     
            info[i,0] = l1[i][0]
            info[i,2] = l2[i][0]
        #print(info[(i+1)*(j+1)-1,:])
        print(info[0,0], '\t', info[0,1], '\t', info[0,2])
        seq1 = l1[0][1]
        seq2 = l2[0][1]
        for ax in ax_arr.flatten():
            plot(seq1, seq2, l1[0][0], l2[0][0], width, marker_size, ax, contrast)
            ax.set_xlim((0,len(seq1)))
            ax.set_ylim((0,len(seq2)))



    #flat and fasta
    elif (len(l1) >= 1 and len (l2) == 1) or (len(l1) == 1 and len(l2) >= 1):
        
        if len(l1) > 1:
            temp = l1
            l1 = l2
            l2 = temp
            temp = ''
        l2 = l2[::-1]
        #contig index
        print('contigs:')
        for i in range(len(l2)):
            print(i, l2[len(l2)-i-1][0])
        #write legend.txt
        if log == 1:
            txt = open('{}_legend.txt'.format(image_name),'w')
            print('contigs:', file=txt)
            for i in range(len(l2)):
                print(i, l2[i][0], file=txt)     
            txt.close()

        print()
        print()
        print(len(l2),'plots:')
        print()
        #screen info
        info = np.array(len(l1)*len(l2)*[[str(0)]*3], dtype=object)
        info[:,1] = ['against']
        #print(info)
        
        for i in range(len(l1)):
            for j in range(len(l2)):     
                info[(i+1)*(j+1)-1,0] = l1[i][0]
                info[(i+1)*(j+1)-1,2] = l2[j][0]

        #name the rows
        rows = ['{}'.format(row) for row in range(len(l2))]
        i=0
        j=-1
        temp = 0
        for ax in ax_arr.flatten():
            if j == len(l2)-1:
                i+=1
                j=-1
            temp+=1
            j+=1
           

            
            #print info while calculating
            print(info[(i+1)*(j+1)-1,0], '\t', info[(i+1)*(j+1)-1,1], '\t', info[(i+1)*(j+1)-1,2])
            seq1 = l1[i][1]
            seq2 = l2[j][1]
            print(str((i+1)*(j+1)))

            #label contig rows by numbers
            ax.annotate(rows[len(l2)-j-1], xy=(1, 0.5), xytext=(pad, 0),
                        xycoords='axes fraction', textcoords='offset points',
                        size='large', ha='left', va='center')  

            #which axes labels & ticks to display
            if j==len(l2)-1:
                ax.set_xlabel(l1[i][0])

            else:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected                        
                                labelbottom=False) # labels along the bottom edge are off    

            plot_no_title(seq1, seq2, l1[i][0], l2[j][0], width, marker_size, ax, contrast)
            ax.set_xlim((0,len(seq1)))
            ax.set_ylim((0,len(seq2)))

            #ax1.set_show() #for seperate plots
        temp=''



    #two flat files or single flat file 
    elif len(l1)>1 and len(l2)>1:
        d_row = {}
        d_col = {}
        l_combo = []
        
        #make dictionary
        for i in range(len(l1)):            
            d_col[i] = l1[i][0]
        
        for i in range(len(l2)):
            d_row[i] = l2[i][0]
        
        #combinations
        for j in range(len(l2)):
            for i in range(len(l1)):            
                    l_combo.append([i, j])
        #l_combo = sorted(l_combo)
        print()


        print('contigs x-axis:', fasta1)
        for i in range(len(d_col)):
            print(i, d_col[i])
        print()
        print()
        print('contigs y-axis:', fasta2)
        for i in range(len(d_row)):
            print(i, d_row[i])        
        print()
        print()
        print('combinations:')
        print(l_combo)
        print()
        print()
        print(len(l_combo),'plots:')
        print()

        if log == 1:
            txt = open('{}_legend.txt'.format(image_name),'w')
            print('contigs x-axis:', fasta1, file=txt)
            for i in range(len(d_col)):
                print(i, d_col[i], file=txt)
            print(file=txt)
            print(file=txt)
            print('contigs y-axis:', fasta2, file=txt)
            for i in range(len(d_row)):
                print(i, d_row[i], file=txt)        
            txt.close()

        cols = ['{}'.format(col) for col in range(len(l1))]
        rows = ['{}'.format(row) for row in range(len(l2))]

        #order plots
        l_ordered = []
        for i in range(len(l_combo)):
            l_ordered.append( [l_combo[i][0], len(l2)-1 - l_combo[i][1]] )

        temp=''
        axes_index = l_combo
        l_combo = l_ordered
        print(axes_index)
        print(l_ordered)

        #screen info
        info = np.array(len(l1)*len(l2)*[[str(0)]*3], dtype=object)
        info[:,1] = ['against']
        for i in range(len(l_combo)):
            info[i][0] = l1[l_combo[i][0]][0]
            info[i][2] = l2[l_combo[i][1]][0]

        #for multiplots
        for i,ax in zip(range(len(l_combo)),ax_arr.flatten()):
            print(info[i])
            
            #print(info[i,0], '\t', info[i,1], '\t', info[i,2])
            #print(d_col[l_combo[i][0]], '\t', 'against', '\t', d_row[l_combo[i][1]])
            
            seq1 = l1[l_combo[i][0]][1] #x-axis
            seq2 = l2[l_combo[i][1]][1] #y-axis
            
            y=l_combo[i][1]
            x=axes_index[i][0]
            
            ax.set_xlim((0,len(seq1)))
            ax.set_ylim((0,len(seq2)))
            #label columns
            if y == len(l2)-1:
                ax.annotate(cols[x], xy=(0.5, 1), xytext=(0, pad),
                            xycoords='axes fraction', textcoords='offset points',
                            size='large', ha='center', va='baseline') 
            
            if x == len(l1)-1:
                ax.annotate(rows[y], xy=(1, 0.5), xytext=(pad, 0),
                            xycoords='axes fraction', textcoords='offset points',
                            size='large', ha='left', va='center')
            
            #which axes labels & ticks to display            
            if x == 0 and y == 0:    
                ax.tick_params( axis='both',
                                which='both')                
            
            elif y != 0 and x !=0:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom=False,      # ticks along the bottom edge are off
                                left=False,
                                top=False,         # ticks along the top edge are off
                                labelbottom=False,
                                labelleft=False) # labels along the bottom edge are off     
            elif y == 0 and x !=0:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                left=False,
                                top=False,         # ticks along the top edge are off
                                labelleft=False) # labels along the bottom edge are off   
               
            elif y !=0 and x ==0:
                ax.tick_params( axis='both',          # changes apply to the x-axis
                                which='both',      # both major and minor ticks are affected
                                bottom=False,      # ticks along the bottom edge are off
                                
                                top=False,         # ticks along the top edge are off
                                labelbottom=False) 

            plot_no_title(seq1, seq2, d_col[l_combo[i][0]], d_row[l_combo[i][1]], width, marker_size, ax, contrast)

    #set contrast by option
    if p != 1:
        figure1.subplots_adjust(bottom=0.25)
        ax1_contrast = figure1.add_axes([0.12, 0.1, 0.78, 0.03])
        s_contrast = Slider(ax1_contrast, 'Contrast', 0, 1, valinit=(marker_size-0.005)/.995)

        #update contrast
        def update(val):
            contrast = s_contrast.val*.995+.005
            for ax in ax_arr.flatten():
                ax.clear()
            plot_update(l1, l2, figure1, width, marker_size, contrast, pad, f1, f2, s_contrast.val, ax_arr)
        
        s_contrast.on_changed(update)

    #print sequences on hovering with mouse    
    def callback(event):
        if event.xdata != None and event.ydata != None:
            
            for i, ax in enumerate(ax_arr.flatten()):
        
                # For infomation, print which plot the mouse is in
                if ax == np.any(event.inaxes):
                    print ("mouse is in axes plot{}".format((len(l2)-1-(i)//len(l1))*len(l1)+(i)%len(l1)+1))
            # final print of seq
            if event.inaxes in ax_arr.flatten():
                
                print(int(event.xdata), '\t', int(event.ydata))
                #print(seq1)
                print(l1[(i)%len(l1)][1][int(event.xdata):int(event.xdata)+20:], '\t',  l1[(i)%len(l1)][1][int(event.xdata):int(event.xdata)+20:])
                print(l2[len(l2)-1-(i)//len(l1)][1][int(event.ydata):int(event.ydata)+20:], '\t',   l2[len(l2)-1-(i)//len(l1)][1][int(event.ydata):int(event.ydata)+20:][::-1])
                print()

    #connect mouse movement with callback event
    figure1.canvas.mpl_connect('motion_notify_event', callback)

    # option equal
    if eq == 1 and len(ax_arr.flatten()) == 1:
        for ax in ax_arr.flatten():
            ax.set_aspect('equal')

    #normal interactive plot
    if p == 0:
        pl.show()

    #option image save as png
    if p == 1:
        pl.savefig('{}.png'.format(image_name), dpi=300)
    sys.exit()

if __name__ == "__main__":
    main(sys.argv[1:])
