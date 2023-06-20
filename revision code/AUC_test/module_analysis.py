import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.stats import norm
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.neighbors import KernelDensity
import statsmodels.api as sm
import time
from matplotlib_venn import venn2

def get_data(fil_data,fil_annotation,if_normalize=True,output_folder=None):  
    ## read annotation file
    print('### read annotation file: ')
    fil=open(fil_annotation)
    ct=0
    meta_annotation=[]
    for line in fil: 
        if ct==0:
            header_annotation=line.strip().split('\t')
            print('header_annotation')
            print(header_annotation)
        else:
            meta_annotation.append(line.strip().split('\t'))
        ct+=1
    fil.close()
    print('\n')
    
    ## read data file 
    print('### read data file: ')
    fil=open(fil_data)
    ct=0
    tempMetaData={}
    RawCountData=[]
    FullGeneNames=[]
    for line in fil: 
        if ct==0:
            tempMetaData['ID']=line.strip().split('\t')[1:]
        elif ct==1:
            tempMetaData['Batch']=line.strip().split('\t')[1:]
        elif ct==2:
            tempMetaData['Date']=line.strip().split('\t')[1:]
        elif ct==3:
            tempMetaData['IRIS']=line.strip().split('\t')[1:]
        else:
            temp=line.strip().split('\t')
            FullGeneNames.append(temp[0])
            RawCountData.append(temp[1:])
        ct+=1
    fil.close()
    
    ## convert the RawCountData matrix to float matrix with nan values 
    RawCountData = np.asarray(RawCountData).T     
    temp_shape = RawCountData.shape
    RawCountData = RawCountData.astype('bytes').reshape(-1)
    RawCountData = np.genfromtxt(RawCountData).reshape(temp_shape)
    #print('test')
    #print(RawCountData[0:5,0:5])
    
    n_sample,n_gene=RawCountData.shape
    print('n_sample=%d, n_gene=%d\n'%(n_sample,n_gene))
    
    ## cross reference the annotation file and the data file
    print('### cross reference the annotation and the data: ')
    FullMetaData=[]
    for i in range(n_sample):
        sample_id = tempMetaData['ID'][i]
        if_appended = False
        for j in range(len(meta_annotation)):
            temp_sample_id=meta_annotation[j][0]
            if sample_id==temp_sample_id:
                FullMetaData.append(meta_annotation[j])
                if_appended=True
                break
        if if_appended==False:
            FullMetaData.append([sample_id,'non_info'])
            print('!!! non info in the annotation file: ',FullMetaData[i])
    print('\n')
    
    #for i in range(len(FullMetaData)):
    #    print(FullMetaData[i])
    
    ## filtering and normalization
    print('### filtering and normalization')
    print('if_normalize: %s'%str(if_normalize))
    if if_normalize:            
        # filter out low-count genes:
        filGeneMeanCtNum=1
        print('Gene filter threshold = %0.2f'%filGeneMeanCtNum)
        mean_gene_count=np.nanmean(RawCountData,axis=0)
        GeneNames=[]
        for i in range(n_gene):
            if mean_gene_count[i]>filGeneMeanCtNum:
                GeneNames.append(FullGeneNames[i])
        CountData=RawCountData[:,mean_gene_count>filGeneMeanCtNum]

        # filter out low-count samples 
        filSampMeanCtNum=np.max([1,np.nanmean(CountData)*0.2])
        #filSampMeanCtNum = 1
        print('Sample filter threshold = %0.2f'%filSampMeanCtNum)
        MetaData=[]
        mean_samp_count=np.nanmean(CountData,axis=1)
        for i in range(mean_samp_count.shape[0]):
            if mean_samp_count[i]>filSampMeanCtNum:
                MetaData.append(FullMetaData[i])

        CountData=CountData[mean_samp_count>filSampMeanCtNum,:]   
         
        ## normalize by the size factor
        n_sample,n_gene = CountData.shape
        ref_gene=np.zeros([n_gene],dtype=float)
        for i in range(n_gene):
            temp=CountData[:,i]+1
            ref_gene[i]=sp.stats.mstats.gmean(temp[np.isnan(temp)==0])
        
        ratio_gene=CountData/ref_gene        
        size_factor=np.zeros([n_sample],dtype=float)
        for i in range(n_sample):
            temp=ratio_gene[i,:]
            temp=temp[np.isnan(temp)==0]
            size_factor[i]=np.median(temp[temp>0])  
            
        #plt.figure()
        #plt.hist(size_factor,bins=20)
        #plt.show()
        
        CountData=(CountData.T/size_factor).T          
    else:
        CountData = RawCountData
        MetaData  = FullMetaData
        GeneNames = FullGeneNames
        
    print('### filtering and normalization completed')
    print('n_sample=%d, n_gene=%d\n'%(n_sample,n_gene))
    return CountData,MetaData,GeneNames,header_annotation


def data_write(count_data,meta_data,gene_name,header,output_file):
    f=open(output_file,'w')
    n_sample=len(meta_data)
    n_gene=len(gene_name)
    for i in range(len(header)):
        f.write(header[i])
        for j in range(n_sample):
            f.write('\t'+meta_data[j][i])
        f.write('\n')
    for i in range(n_gene):
        f.write(gene_name[i])
        for j in range(n_sample):
            f.write('\t'+str(count_data[j][i]))
        f.write('\n')
    f.close()
    return

def get_label(MetaData,TT_list,group_config,daytol=30):
    print('### extract the data label')
    name_list={}
    labels={}
    n_sample=len(MetaData)
    
    ## Extract the ID labels. 
    label_ID=np.zeros([n_sample],dtype=int)
    ID_list=[]
    for i in range(n_sample):
        find_ID=0
        for j in range(len(ID_list)):
            if MetaData[i][1]==ID_list[j]:
                label_ID[i]=j
                find_ID=1
        if find_ID==0:
            ID_list.append(MetaData[i][1])
            label_ID[i]=len(ID_list)-1  
    name_list['ID_list']=ID_list
    labels['label_ID']=label_ID
    
    ## Extract the healthy status (HS) labels:
    label_HS=np.zeros(n_sample,dtype=int)
    HS_list=['Healthy', 'Infection_Late', 'Infection_Recovery_Late', 'Infection_Early', 'Infection_Middle',\
             'Infection_Recovery_Early', 'Ant_Middle', 'Ant_Recovery_Early', 'Imz_Early', 'Imz_Late', \
             'Ant_Recovery_Late', 'Imz_Middle', 'Imz_Recovery_Early', 'Imz_Recovery_Late', 'Ant_Early', 'Ant_Late'] 
    for i in range(n_sample):
        for j in range(len(HS_list)):
            if MetaData[i][6] == HS_list[j]:
                label_HS[i]=j
    name_list['HS_list']=HS_list
    labels['label_HS']=label_HS 
    
    ## Extract IRIS label:   
    label_IRIS=np.zeros(n_sample,dtype=int)
    IRIS_list=['IR','IS','Unknown']   
    for i in range(n_sample):
        for j in range(len(IRIS_list)):
            if MetaData[i][4] == IRIS_list[j]:
                label_IRIS[i]=j
    name_list['IRIS_list']=IRIS_list
    labels['label_IRIS']=label_IRIS          

    ## Extract Diabetic class (DBC) label:   
    label_DBC=np.zeros(n_sample,dtype=int)
    DBC_list=['Diabetic', 'NA', 'Violator', 'Prediabetic', 'Control']   
    for i in range(n_sample):
        for j in range(len(DBC_list)):
            if MetaData[i][5] == DBC_list[j]:
                label_DBC[i]=j
    name_list['DBC_list']=DBC_list
    labels['label_DBC']=label_DBC 
    
    ## Extract batch label:   
    batch_list=[]
    for i in range(n_sample):
        if_find=0
        for j in range(len(batch_list)):
            if MetaData[i][2]==batch_list[j]:
                if_find=1
        if if_find==0:
            batch_list.append(MetaData[i][2])
    label_batch=np.zeros(n_sample,dtype=int)
    for i in range(n_sample):
        for j in range(len(batch_list)):
            if MetaData[i][2] == batch_list[j]:
                label_batch[i]=j
    name_list['batch_list']=batch_list
    labels['label_batch']=label_batch  
    
    ## Extract the time-series label
    # standardizing the date
    mon_list=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    label_date=[]
    for i in range(n_sample):
        s=MetaData[i][3]
        s=s.strip().split('-')
        if len(s)<3:
            label_date.append(None)
        else:
            yyyy='20'+s[2]
            dd=s[0]
            mm=0
            for i in range(len(mon_list)):
                if s[1]==mon_list[i]:
                    mm=str(i+1)
            if mm==0:
                label_date.append(None)
            else:
                label_date.append(yyyy+'-'+mm+'-'+dd)
    labels['label_date']=np.array(label_date)    
          
    label_TT=np.zeros([n_sample],dtype=int)+len(TT_list)-1
    dis_start=1
    dis_end=5
    
    # insert time points of the category of interest 
    for i in range(n_sample):
        temp_HS=MetaData[i][6]
        for j in range(dis_start,dis_end+1):
            if temp_HS==TT_list[j]:
                label_TT[i]=j           

    ## insert the pre/ post dicease time points 
    samp_idx=np.arange(n_sample)
    for idx_ID in range(len(ID_list)):
        ID_mask=(label_ID==idx_ID)
        temp_idx=samp_idx[ID_mask]
        temp_n_sample=temp_idx.shape[0]
        ## for every healthy sample, whose label is else and has a date
        for i in range(temp_n_sample): 
            if label_TT[temp_idx[i]]==len(TT_list)-1 \
            and (label_date[temp_idx[i]] is not None) and label_HS[temp_idx[i]]==0:   
                for j in range(temp_n_sample):      
                    if label_TT[temp_idx[j]]>=dis_start and label_TT[temp_idx[j]]<=dis_end \
                    and (label_date[temp_idx[j]] is not None):
                        datedif=date_dif(label_date[temp_idx[i]],label_date[temp_idx[j]])
                        #print(datedif)
                        if 0<datedif<daytol:
                            label_TT[temp_idx[i]]=dis_start-1
                        elif -daytol<datedif<=0:
                            label_TT[temp_idx[i]]=dis_end+1               
                              
    ## Eliminating samples taken within the same day
    for i in range(len(TT_list)-1):
        mask=(label_TT==i)
        mask_idx=np.arange(n_sample)[mask]
        n_mask=np.sum(mask)
        for j in range(n_mask):
            for k in range(j):
                if label_date[mask_idx[j]]!=None and label_date[mask_idx[k]]!=None:
                    if (date_dif(label_date[mask_idx[j]],label_date[mask_idx[k]])==0)\
                    and (label_ID[mask_idx[j]]==label_ID[mask_idx[k]]):
                        print('hit')
                        print(MetaData[mask_idx[j]])
                        print(MetaData[mask_idx[k]])
                        label_TT[mask_idx[j]]=len(TT_list)-1
                        
    ## Summary:
    for i in range(len(TT_list)):
        print('# ',TT_list[i],': ',np.sum(label_TT==i))
        
    name_list['TT_list']=TT_list
    labels['label_TT']=label_TT 
    
    ## converting the label_TT to label_group
    label_group=np.zeros([label_TT.shape[0]],dtype=int)
    for i in range(len(group_config)):
        for j in group_config[i]:
            label_group[label_TT==j]=i
    labels['label_group']=label_group
    return name_list,labels

def impute(CountData,label_ID,label_HS):
    temp=np.where(np.isnan(CountData))
    mean_count=np.nanmean(CountData)
    n_nan=temp[0].shape[0]
    for i in range(n_nan):
        x=temp[0][i]
        y=temp[1][i]
        cd_gene=CountData[:,y]
        ## searching for samples to impute the data
        if np.sum((label_ID==label_ID[x])*(label_HS==label_HS[x])*(np.isnan(cd_gene)==0))>1:
            CountData[x,y]=np.mean(cd_gene[(label_ID==label_ID[x])*(label_HS==label_HS[x])*(np.isnan(cd_gene)==0)])
        elif np.sum((label_HS==label_HS[x])*(np.isnan(cd_gene)==0))>1:
            CountData[x,y]=np.mean(cd_gene[(label_HS==label_HS[x])*(np.isnan(cd_gene)==0)])
        else:
            CountData[x,y]=mean_count
    return CountData

def date_dif(d1,d2):
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return (d2 - d1).days

## balance the sample proportion in each group
def sample_balance(label_group,group_list,label_ID,label_select,sample_balance_percentage=0.5): 
    #return np.ones([label_group.shape[0]],dtype=bool)
    for i in range(len(group_list)-1):
        thres=np.max([int(round(np.sum((label_group==i)*label_select)*sample_balance_percentage)),1])
        #print(group_list[i],np.sum((label_group==i)*label_select),thres,np.sum((label_group==i)*label_select*(label_ID==0)))
        idx_group=np.arange(label_group.shape[0])[(label_group==i)*label_select] # the indices of all samples for a group
        dump_idx=[]
        #print(label_ID[idx_group])
        for j in np.unique(label_ID[idx_group]):
            if np.sum(label_ID[idx_group]==j)>thres: ## the number of samples for one ID is larger than the threshold
                temp_idx=idx_group[label_ID[idx_group]==j][np.random.permutation(np.sum(label_ID[idx_group]==j))[thres:]]
                dump_idx+=list(temp_idx)
        #print(dump_idx)
        label_select[dump_idx]=0
    return label_select
            

def analysis(CountData,GeneNames,config,vis=0):
    ## load config parameters
    FDR=config['FDR']
    TT_list=config['TT_list']
    label_TT=config['label_TT']
    label_group=config['label_group']
    label_date=config['label_date']
    label_ID=config['label_ID']
    label_HS=config['label_HS']
    output_folder=config['output_folder']
    suffix=config['suffix']
    group_config=config['group_config']
    group_list=config['group_list']  
    win_size=config['win_size']
    sample_balance_percentage=config['sample_balance_percentage']
    debug_mode=config['debug_mode']
    
    ## write the basic information     
    f=open(output_folder+'/analysis_summary'+suffix+'.txt','w')
    f.write(config['Title']+'\n')
    for j in range(len(group_list)):
        f.write(group_list[j]+'\t'+str(np.sum(label_group==j))+'\n')
   
    ## analysis
    n_sample,n_gene=CountData.shape
    logCountData=np.log2(CountData+1)
    
    ## testing using paired_AUC statistic
    label_select=sample_balance(label_group,group_list,label_ID,np.ones([n_sample],dtype=bool),\
                                sample_balance_percentage=sample_balance_percentage)
    check_logCountData,G_mean,G_var,G_n=Stats_Cal(CountData[label_select,:],label_group[label_select],group_list)
    if check_logCountData is None:
        f.write('## not enough samples!!!\n')
        f.close()
        return None,None,None
    
    ## output the log2 fold change # fixit 
    #write_log2_fc(CountData,label_TT,TT_list,label_group,group_list,GeneNames,output_folder,suffix)
    p_cal_each_group(logCountData[label_select,:],label_group[label_select],\
                     label_ID[label_select],label_HS[label_select],\
                     label_date[label_select],group_list,output_folder,\
                     label_TT[label_select],TT_list,GeneNames,vis=vis,suffix=suffix,\
                     win_size=win_size,debug_mode=debug_mode)
    
    ## test by paired_AUC 
    analysis_label=np.zeros([n_sample],dtype=int)+len(group_list)-1
    p_TT,analysis_label[label_select],T_auc = p_cal_pair_auc(logCountData[label_select,:],label_group[label_select],\
                                                             label_ID[label_select],label_HS[label_select],\
                                                             label_date[label_select],group_list,output_folder,\
                                                             label_TT[label_select],TT_list,GeneNames,vis=vis,suffix=suffix,\
                                                             win_size=win_size) 
    sig_idx_pauc = post_p_analysis(p_TT,T_auc,GeneNames,group_list,analysis_label,FDR,\
                                   output_folder,suffix='_pAUC'+suffix,debug_mode=debug_mode)
       
    ## test by linear regression
    if debug_mode:
        analysis_label = np.zeros([n_sample],dtype=int)+len(group_list)-1
        p_lr,analysis_label[label_select],T_lr = p_cal_lr(logCountData[label_select,:],label_group[label_select],\
                                                          label_ID[label_select],label_HS[label_select],\
                                                          label_date[label_select],group_list,\
                                                          label_TT[label_select],TT_list,GeneNames)
        sig_idx_lr = post_p_analysis(p_lr,T_lr,GeneNames,group_list,\
                                     analysis_label,FDR,output_folder,suffix='_LR'+suffix)


        ## generate the venn graph
        n_common = np.intersect1d(sig_idx_pauc,sig_idx_lr).shape[0]
        plt.figure()
        venn2(subsets = (sig_idx_pauc.shape[0]-n_common, sig_idx_lr.shape[0]-n_common, n_common),set_labels = ('pauc','lr'))
        plt.savefig(output_folder+'/venn'+suffix+'.pdf')
        plt.savefig(output_folder+'/venn'+suffix+'.png')
        plt.close()   
    f.close()
    
    ## write down the cluster data information (z-scored data)
    write_cluster_data(GeneNames,G_mean,sig_idx_pauc,output_folder)
    
    return sig_idx_pauc,None,analysis_label
    
""" 
    testing methods for calculating the p-value
"""

def p_cal_each_group(logCountData,label_group,label_ID,label_HS,label_date,group_list,output_folder,\
                   label_TT,TT_list,GeneNames,vis=0,suffix='',win_size=183,debug_mode=True):     
    
    ## some calculation 
    n_sample,n_gene=logCountData.shape
    ## convert the date label 
    date_ref='2000-1-1'
    label_date_to_ref=np.zeros([n_sample],dtype=int)
    for i in range(n_sample):
        if label_date[i] is None:
            label_date_to_ref[i]=-1000
        else:
            label_date_to_ref[i]=date_dif(date_ref,label_date[i])
            
    ## normalizing the data: find out the reference list 
    temp_logCountData = logCountData.copy()
    label_select = np.zeros([n_sample],dtype=bool)
    for ID in np.unique(label_ID): # iterate by ID 
        for j in np.arange(n_sample)[(label_ID==ID)*(label_group>0)*(label_group<4)]:
            ref_idx = np.arange(n_sample)[(label_ID==ID)*(label_HS==0)\
                                          *(np.absolute(label_date_to_ref-label_date_to_ref[j])<win_size)]        
            if ref_idx.shape[0]>0:
                ref_data = logCountData[ref_idx].mean(axis=0)
                temp_logCountData[j] -= ref_data
                label_select[j]=True
          
    f_fc=open(output_folder+'/log2fc_personal'+suffix+'.txt','w')
    f_fc.write('Attribute')  
    for i in range(n_gene):
        f_fc.write('\t'+GeneNames[i]+'\t'+'p-value')
    f_fc.write('\n')
    
    ## write the fc and p-value by TT
    for i in range(1,len(TT_list)-1):        
        f_fc.write('TT_'+TT_list[i])
        if np.sum((label_TT==i)*label_select)<2:
            f_fc.write('not enough samples\n')
        else:            
            temp_mean= np.mean(temp_logCountData[(label_TT==i)*label_select],axis=0)
            temp_fc=temp_mean
            temp_var = np.var(temp_logCountData[(label_TT==i)*label_select],axis=0)
            temp_var[temp_var<0.001]=1e12
            temp_n = np.sum((label_TT==i)*label_select)

            temp_T = temp_mean/np.sqrt(temp_var/temp_n)
            temp_p_value=2*(1-norm.cdf(np.absolute(temp_T)))
            for j in range(n_gene):
                if np.sum(label_TT==i)>0:
                    f_fc.write('\t'+str(temp_fc[j])+'\t'+str(temp_p_value[j]))
                else:
                    f_fc.write('\tna')
            f_fc.write('\n')  
            
            n_rej,_,_ = BHq(temp_p_value, alpha = 0.1)
            
            if debug_mode:
                plt.figure()
                plt.hist(temp_p_value,bins=100)
                plt.savefig(output_folder+'/p_hist_%s%s.png'%(TT_list[i],suffix))
                plt.savefig(output_folder+'/p_hist_%s%s.pdf'%(TT_list[i],suffix))
                plt.close()
        
    ## write the fc and p-value by group             
    for i in range(1,len(group_list)-1):        
        f_fc.write('group'+group_list[i])
        if np.sum((label_group==i)*label_select)<2:
            f_fc.write('not enough samples\n')
        else:
            
            temp_mean= np.mean(temp_logCountData[(label_group==i)*label_select],axis=0)
            temp_fc=temp_mean
            temp_var = np.var(temp_logCountData[(label_group==i)*label_select],axis=0)
            temp_var[temp_var<0.001]=1e12
            temp_n = np.sum((label_group==i)*label_select)

            temp_T = temp_mean/np.sqrt(temp_var/temp_n)
            temp_p_value=2*(1-norm.cdf(np.absolute(temp_T)))
            
            for j in range(n_gene):
                if np.sum(label_group==i)>0:
                    f_fc.write('\t'+str(temp_fc[j])+'\t'+str(temp_p_value[j]))
                else:
                    f_fc.write('\tna')
            f_fc.write('\n')
            
            n_rej,_,_ = BHq(temp_p_value, alpha = 0.1)
            #print(group_list[i], temp_n,n_rej)
            
            if debug_mode:
                plt.figure()
                plt.hist(temp_p_value,bins=100)
                plt.savefig(output_folder+'/p_hist_%s%s.png'%(group_list[i],suffix))
                plt.savefig(output_folder+'/p_hist_%s%s.pdf'%(group_list[i],suffix))
                plt.close()
    f_fc.close()
        
def p_cal_pair_auc(logCountData,label_group,label_ID,label_HS,label_date,group_list,output_folder,\
                   label_TT,TT_list,GeneNames,vis=0,suffix='',win_size=183):                  
    n_sample,n_gene=logCountData.shape
  
    ## convert the date label 
    date_ref='2000-1-1'
    label_date_to_ref=np.zeros([n_sample],dtype=int)
    for i in range(n_sample):
        if label_date[i] is None:
            label_date_to_ref[i]=-1000
        else:
            label_date_to_ref[i]=date_dif(date_ref,label_date[i])
            
    ## normalizing the data: find out the reference list 
    ref_list = {}
    label_select = np.zeros([n_sample],dtype=bool)
    for ID in np.unique(label_ID): # iterate by ID 
        for j in np.arange(n_sample)[(label_ID==ID)*(label_group>0)*(label_group<4)]:
            ref_idx = np.arange(n_sample)[(label_ID==ID)*(label_HS==0)\
                                          *(np.absolute(label_date_to_ref-label_date_to_ref[j])<win_size)]           
            if ref_idx.shape[0]>0:
                ref_list[j] = ref_idx                 
                label_select[j]=True
                label_select[ref_idx] = True
    
    ## calculate the testing stastistics
    a = np.zeros([n_sample],dtype=float)
    for i in range(1,4):
        n_g = np.sum(label_select*(label_group==i)) 
        a[label_select*(label_group==i)] = 1/n_g
        for j in np.arange(n_sample)[label_select*(label_group==i)]:
            a[ref_list[j]] -= 1/n_g/(ref_list[j].shape[0])
    
     
    ## estimate the variance 
    temp = logCountData[(label_HS==0)*label_select]
    temp_label_ID = label_ID[(label_HS==0)*label_select]
    temp_select = np.ones([temp.shape[0]],dtype=bool)
    for ID in np.unique(temp_label_ID):
        if np.sum(temp_label_ID==ID)<=2 or ID==0:
            temp_select[temp_label_ID==ID] = False
        else:
            temp[temp_label_ID==ID] -= temp[temp_label_ID==ID].mean(axis=0)          
    sigma_h = np.var(temp[temp_select,:],axis=0)
    
    #print('n_sigma',temp_select.sum())
    #print('sigma_h',sigma_h[0:10]) 
    #print(np.var(logCountData[label_HS==1],axis=0)[0:10])
    #print('mean',a.dot(logCountData)[0:10])
    #plt.figure()
    #plt.hist(np.log10(sigma_h))
    #plt.show()   
    sigma_h[sigma_h<1e-5] = 1e10
              
    ## calculate the statistics 
    T_pair_auc= a.dot(logCountData) / np.sqrt(np.sum(a**2)*sigma_h)
    
    #plt.figure()
    #plt.hist(np.log10(np.absolute(a.dot(logCountData))))
    #plt.show()  
    p_TT=2*(1-norm.cdf(np.absolute(T_pair_auc)))
    
    ## generate the analysis label
    analysis_label=label_group
    analysis_label[label_select==0]=len(group_list)-1
    return p_TT,analysis_label,T_pair_auc

## maintain an estimate of the variance 
    #temp = logCountData[(label_HS==0)*label_select]
    #temp_label_ID = label_ID[(label_HS==0)*label_select]
    #temp_date = label_date_to_ref[(label_HS==0)*label_select]
    #temp_date = (temp_date - temp_date.mean())/np.std(temp_date)
    #
    #temp_select = np.ones([temp.shape[0]],dtype=bool)
    #for ID in np.unique(temp_label_ID):
    #    if np.sum(temp_label_ID==ID)<=1 or ID==0:
    #        temp_select[temp_label_ID==ID] = False
    #
    #temp = temp[temp_select,:]
    #temp_label_ID = temp_label_ID[temp_select]
    #temp_date = temp_date[temp_select]
    #
    #A = np.zeros([temp.shape[0],np.unique(temp_label_ID).shape[0]],dtype=float)
    #for i,ID in enumerate(np.unique(temp_label_ID)):
    #    A[temp_label_ID==ID,i] = 1
   #
    #sigma_h = np.zeros([n_gene],dtype=float)
    #for i in range(n_gene):
    #    if temp[:,i].mean()<1e-4:
    #        continue
    #    mod = sm.OLS(temp[:,i],A)
    #    res = mod.fit()
    #    sigma_h[i] = res.scale
   #
    #sigma_h[np.isnan(sigma_h)]=1000
    #sigma_h[sigma_h<1e-4]=1000   

def p_cal_lr(logCountData,label_group,label_ID,label_HS,label_date,group_list,label_TT,TT_list,GeneNames):
    ## create the design matrix    
    n_sample,n_gene = logCountData.shape
    n_ID = np.unique(label_ID).shape[0]
    A = np.zeros([n_sample,4+n_ID])
    samp_select = np.ones([n_sample],dtype=bool) 
    feature_select = np.ones([A.shape[1]],dtype=bool)
    
    ## update the group information 
    samp_select[label_group==5] = False
    A[label_group==1,0] = 1
    A[label_group==2,1] = 1
    A[label_group==3,2] = 1
    
    ## update the time information
    date_ref='2000-1-1'
    label_date_to_ref=np.zeros([n_sample],dtype=int)
    for i in range(n_sample):
        if label_date[i] is None:
            label_date_to_ref[i]=-1000
            samp_select[i] = False
        else:
            label_date_to_ref[i]=date_dif(date_ref,label_date[i])
    A[:,3] = (label_date_to_ref-label_date_to_ref.mean())/np.std(label_date_to_ref)
    #A[:,4] = (label_date_to_ref-label_date_to_ref.mean())/np.std(label_date_to_ref)
    
    ## update the sample ID 
    for i,ID in enumerate(np.unique(label_ID)):
        A[label_ID==ID,4+i] = 1
        if np.sum((label_ID==ID)*samp_select)<=1:
            samp_select[label_ID==ID] = False
            feature_select[i+4] = False
    
    logCountData = logCountData[samp_select,:]
    A = A[samp_select,:]
    A = A[:,feature_select]
    
    #for i,ID in enumerate(np.unique(label_ID)):
    #    if np.sum((label_ID==ID)*samp_select)>0:
    #        print(i,np.sum((label_ID==ID)*samp_select))
    F_val = np.zeros([n_gene],dtype=float)
    p_val = np.zeros([n_gene],dtype=float)
    
    start_time = time.time()
    for i in range(n_gene):
        mod = sm.OLS(logCountData[:,i],A)
        res = mod.fit()
        F = res.f_test("x1=x2=x3=0")
        F_val[i],p_val[i] = F.fvalue[0][0],F.pvalue

        #print(res.summary())
        #print(F.pvalue)
        
        ## plot the data 
        #temp_data = logCountData[:,i]
        #plt.figure(figsize=[18,5])
        #for i,ID in enumerate(np.unique(label_ID)):
        #    
        #    if np.sum((label_ID==ID)*samp_select)>0:
        #        plt.figure(figsize=[18,5])
        #        temp = (label_ID==ID)*samp_select
        #        plt.scatter(label_group[temp],temp_data[temp[samp_select]],label=ID,s=80,alpha=0.6)
        #        plt.show() 
        #
        #break
            
    analysis_label=label_group
    analysis_label[samp_select==0]=len(group_list)-1
    
    return p_val,analysis_label,F_val
   
    
""" 
    functions for outputing the results
"""

def write_cluster_data(GeneNames,G_mean,sig_idx,output_folder):
    
    temp_std=np.std(G_mean,axis=0)
    idx_bad=(temp_std==0)
    temp_G_mean=np.zeros(G_mean.shape,dtype=float)
    temp_G_mean[:,idx_bad==0]=sp.stats.mstats.zscore(G_mean[:,idx_bad==0],axis=0)
        
    f=open(output_folder+'/cluster_data.txt','w')
    f.write('gene_name\tPre\tEARLY\tLATE\tRECOVERY\tPost\n')
    for i in sig_idx:
        f.write(GeneNames[i]+'\t'+str(temp_G_mean[0,i])+'\t'+str(temp_G_mean[1,i])+'\t'+\
                str(temp_G_mean[2,i])+'\t'+str(temp_G_mean[3,i])+'\t'+str(temp_G_mean[4,i])+'\n')
    f.close()

def write_log2_fc(CountData,label_TT,TT_list,label_group,group_list,GeneNames,output_folder,suffix):
    f_fc=open(output_folder+'/log2fc'+suffix+'.txt','w')
    cds=CountData+1
    logCountData=np.log2(cds)
    n_sample,n_gene=CountData.shape
    f_fc.write('Attribute')  
    for i in range(n_gene):
        f_fc.write('\t'+GeneNames[i]+'\t'+'p-value')
    f_fc.write('\n')
    
    
    if np.sum(label_TT==0)>1:
        mean_pre=np.mean(cds[label_TT==0,:],axis=0)
        log_mean_pre=np.mean(logCountData[label_TT==0,:],axis=0)
        var_pre=np.var(logCountData[label_TT==0,:],axis=0)
        var_pre[var_pre<0.001]=1e12
        n_pre=np.sum(label_TT==0)
    else:
        f_fc.write('not enough samples')
        f_fc.close()
        return
    
    
    for i in range(1,len(TT_list)-1):        
        f_fc.write('TT_'+TT_list[i])
        if np.sum(label_TT==i)>1:
            temp_mean=np.mean(cds[label_TT==i,:],axis=0)
            temp_fc=np.log2(temp_mean/mean_pre)
            
            log_temp_mean=np.mean(logCountData[label_TT==i,:],axis=0)
            temp_var=np.var(logCountData[label_TT==i,:],axis=0)
            temp_var[temp_var<0.001]=1e12
            temp_T=(log_temp_mean-log_mean_pre)/np.sqrt(temp_var/np.sum(label_TT==i)+var_pre/n_pre)
            temp_p_value=2*(1-norm.cdf(np.absolute(temp_T)))
        for j in range(n_gene):
            if np.sum(label_TT==i)>0:
                f_fc.write('\t'+str(temp_fc[j])+'\t'+str(temp_p_value[j]))
            else:
                f_fc.write('\tna')
        f_fc.write('\n')                                
        
    if np.sum(label_group==0)>1:
        mean_pre=np.mean(cds[label_group==0,:],axis=0)
        log_mean_pre=np.mean(logCountData[label_group==0,:],axis=0)
        var_pre=np.var(logCountData[label_group==0,:],axis=0)
        var_pre[var_pre<0.001]=1e12
        n_pre=np.sum(label_group==0)
    else:
        f_fc.write('not enough samples')
        f_fc.close()
    
    for i in range(1,len(group_list)-1):        
        f_fc.write('group'+group_list[i])
        if np.sum(label_group==i)>1:
            temp_mean=np.mean(cds[label_group==i,:],axis=0)
            temp_fc=np.log2(temp_mean/mean_pre)
            
            log_temp_mean=np.mean(logCountData[label_group==i,:],axis=0)
            temp_var=np.var(logCountData[label_group==i,:],axis=0)
            temp_var[temp_var<0.001]=1e12
            
            temp_T=(log_temp_mean-log_mean_pre)/np.sqrt(temp_var/np.sum(label_group==i)+var_pre/n_pre)
            temp_p_value=2*(1-norm.cdf(np.absolute(temp_T)))
            
        for j in range(n_gene):
            if np.sum(label_group==i)>0:
                f_fc.write('\t'+str(temp_fc[j])+'\t'+str(temp_p_value[j]))
            else:
                f_fc.write('\tna')
        f_fc.write('\n')
        
    return

"""
    calculate the mean and variance of each groups
"""

def Stats_Cal(CountData,label_group,group_list): ## automatically leave out the last group, 'Else'
    n_sample,n_gene=CountData.shape
    logCountData=np.log2(CountData+1)
    
    G_mean=np.zeros([len(group_list)-1,n_gene])
    G_var=np.zeros([len(group_list)-1,n_gene])
    G_n=np.zeros([len(group_list)-1],dtype=int)

    for i in range(len(group_list)-1):
        if np.sum(label_group==i)>1:            
            G_mean[i,:]=np.mean(logCountData[label_group==i],axis=0)
            G_var[i,:]=np.var(logCountData[label_group==i],axis=0)
            G_n[i]=np.sum(label_group==i)  
        else:
            return None,None,None,None
    return logCountData,G_mean,G_var,G_n

"""
    writing information after p-value calculation
"""
def post_p_analysis(p_val,T_val,GeneNames,group_list,analysis_label,FDR,output_folder,suffix='',debug_mode=False):
    if debug_mode is False:
        
        ## generate the figures for p_val and for T_val 
        plt.figure(figsize=[10,8])
        plt.hist(p_val,bins=100)
        plt.xlabel('p-value')
        plt.savefig(output_folder+'/p_val_histogram'+suffix+'.pdf')
        plt.savefig(output_folder+'/p_val_histogram'+suffix+'.png')
        plt.close()
    else:
        ## generate the figures for p_val and for T_val 
        plt.figure(figsize=[18,8])
        plt.subplot(121)
        plt.hist(T_val,bins=100)
        #plt.ylim([0,2000])
        plt.title('test statistics')
        plt.subplot(122)
        plt.hist(p_val,bins=100)
        #plt.ylim([0,2000])
        plt.title('p-value')
        plt.suptitle(suffix)
        plt.savefig(output_folder+'/p_val_histogram'+suffix+'.pdf')
        plt.savefig(output_folder+'/p_val_histogram'+suffix+'.png')
        plt.close()
    
    ## write the AUC statistics as well as the p-values
    f_auc=open(output_folder+'/test_stats'+suffix+'.txt','w')    
    f_auc.write(' \t')
    for i in range(len(GeneNames)):
        f_auc.write(GeneNames[i]+'\t')
    f_auc.write('\np_values\t')
    for i in range(len(GeneNames)):
        f_auc.write(str(p_val[i])+'\t')
    f_auc.write('\nstatistics\t')
    for i in range(len(GeneNames)):
        f_auc.write(str(T_val[i])+'\t')
    f_auc.close()
   
    # write the analysis statistics
    n_gene = p_val.shape[0]
    f=open(output_folder+'/summary'+suffix+'.txt','w')    
    f.write('## sample number in each group in actual analysis##\n')
    for i in range(len(group_list)-1):
        f.write('# %s: %s, '%(group_list[i],str(np.sum(analysis_label==i))))
    f.write('\n')
    #n_rej,n_t,rej_list=BHq(p_val,alpha=FDR)
    n_rej,n_t,rej_list=sBHq(p_val,alpha=FDR)
    sig_idx=np.arange(n_gene)[rej_list==1]
    sort_idx=np.argsort(p_val[rej_list==1])
    sig_idx=sig_idx[sort_idx]
    sig_GeneNames=[]
    for i in sig_idx:
        sig_GeneNames.append(GeneNames[i])
    f.write('number of rejections: %s'%str(n_rej)+'\t FDR level: %s\n'%str(FDR))
    f.write('rejected gene list:\n')
    for i in range(len(sig_GeneNames)):
        f.write(sig_GeneNames[i])
        f.write('\t')
    f.write('\n')
    f.close()
    return sig_idx

""" 
    multiple control procedure
"""

def BHq(x, alpha = 0.05):
    x_s = sorted(x)
    #n = len(x_s)
    n = np.sum(x<0.99).astype(np.int)
    ic = 0
    for i in range(n):
        if x_s[i] < i*alpha/float(n):
            ic = i
    return ic, x_s[ic],(x<x_s[ic])

def sBHq(x, alpha = 0.05):
    x_s = sorted(x)
    #n = len(x_s)
    n = np.sum(x<0.99).astype(np.int)
    pi_0 = np.sum((x>0.5)*(x<0.9))/0.4/n
    pi_0 = np.min([pi_0,1])
    alpha = alpha/pi_0
    ic = 0
    for i in range(n):
        if x_s[i] < i*alpha/float(n):
            ic = i
    return ic, x_s[ic],(x<x_s[ic])