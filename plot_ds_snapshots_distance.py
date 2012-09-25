#!/usr/bin/env python

#plot all statistics obtained by parsing snapshots from detailed disksim simulations.

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

import sys
sys.path.append("/uf21/lj5bp/py_lib")
sys.path.append("/uf21/lj5bp/Documents/research/opencv/lib/python2.6/dist-packages")
sys.path.append("/uf21/lj5bp/Documents/research/opencv/lib")

import os,csv,re,math
#os.chdir("/uf21/lj5bp/Documents/research/opencv/lib")
import cv

import numpy as np
import scipy.stats as stats
import pdb, gc, pickle

from SnapshotFile import SnapshotFile
from SnapshotRecords import SnapshotRecords
from fileutils import *
from csvfileutils import *
from block_metadata import block_metadata

import argparse, colorama
from colorama import Fore
from get_label import *


year = 365 * 24 * 60 * 60
fmt = ''
# Disksim simulations snapshots are saved for every 15 day simulation period.
snapshot_frequency = 0.5 #frequency of disksim snapshots.
start_year = 1
num_years = 1
max_retention = 1
histogram = False
normalize = False
out_dir = ''
save_table=False
min_intensity = 1
max_intensity = 15
scale = 12
filename = []
def main(argv=None):

  colorama.init(autoreset=True)
  if argv is None:
    argv = sys.argv

  parser = argparse.ArgumentParser(description='Script to plot graphs from reliability snapshots')
  parser.add_argument('-i', '--in_root', required=True,help='Root input directory')
  parser.add_argument('-o', '--out_dir', default='.',help='Output directory to save figures')
  parser.add_argument('-w', '--workload', default=None,action='append',help='List of workloads to be analyzed')
  parser.add_argument('-f','--fmt',default='png',choices=['eps','png'],help='Output file in png/eps format')
  parser.add_argument('-s','--start_year',type=float,default=0.0,help='Starting snapshot (in years)')
  parser.add_argument('-n','--num_years',type=float,default=5.0,help='Number of years to plot (in years)')
  parser.add_argument('-m','--max_retention',type=float,default=35.0,help='Max retention period to consider')
  parser.add_argument('-I','--min_intensity',type=float,default=1,help='Min intensity to consider')
  parser.add_argument('--max_intensity',type=float,default=20,help='Max intensity to consider')
  parser.add_argument('--scale',type=float,default=12,help='Pick points every x*0.5 months')
  parser.add_argument('-H',action="store_true",default=False,help='Plots histogram instead of CDF')
  parser.add_argument('-N',action="store_true",default=False,help='Normalize results?')
  parser.add_argument('-t',action="store_true",default=False,help='Save data in table')
  parser.add_argument('-F','--file',default=None,action='append',help='Files to process')
  parser.add_argument('-r',action="store_true",default=True,help='Model retention?')
  parser.add_argument('-l',action="store_true",default=False,help='Model logical stress?')
  parser.add_argument('-e',action="store_true",default=False,help='Model equivalent cycles?')
  parser.add_argument('-V',action="store_true",default=False,help='Model equivalent cycles?')
  parser.add_argument('--hellinger_distance_discrete', action="store_true", default=False, help='Plot discrete Hellinger Distance')
  parser.add_argument('--hellinger_distance_normalize', action="store_true", default=False, help='Plot normalized Hellinger Distance and shift mean distance')
  parser.add_argument('--earth_mover_distance', action="store_true", default=False, help="Plot earth mover distance")
  parser.add_argument('--analyze_error', action="store_true", default=False, help="Plot Error Analysis")

  args = parser.parse_args()
  workloads = args.workload
  trace_root = os.path.abspath(args.in_root)
  global fmt
  fmt = args.fmt
  global start_year
  start_year = args.start_year
  global num_years
  num_years = args.num_years
  global max_retention
  max_retention = args.max_retention
  global histogram
  histogram = args.H
  global normalize
  normalize = args.N
  global out_dir
  out_dir = args.out_dir
  global save_table
  save_table = args.t
  global min_intensity
  min_intensity = args.min_intensity
  global max_intensity
  max_intensity = args.max_intensity
  global filename
  filename = args.file
  mkdir_p(out_dir)
  global plot_retention
  plot_retention = args.r
  global plot_stress
  plot_stress = args.l
  global plot_eqn
  plot_eqn = args.e
  global scale
  scale = args.scale  
  global plot_delvth
  plot_delvth = args.V
  global hellinger_discrete
  hellinger_discrete = args.hellinger_distance_discrete
  global hellinger_normalize
  hellinger_normalize = args.hellinger_distance_normalize
  global earth_mover
  earth_mover = args.earth_mover_distance
  global analyze_err
  analyze_err = args.analyze_error
  global earth_mover_dist
  earth_mover_dist = []
  global err_in_dist
  err_in_dist = []

  if plot_delvth == True:
    print 'Plotting delvth'
  if plot_retention == True:
    print 'Plotting retention'
  if plot_eqn == True:
    print 'Plotting equivalent cycles'
  if plot_stress == True:
    print 'Plotting logical stresses'
  if hellinger_discrete == True:
	print 'Plotting Hellinger Discrete'
  if hellinger_normalize == True:
	print 'Plotting Hellinger Normalized and shift mean distance'
  if earth_mover == True:
	print 'Plotting Earth Mover Distance'
  if args.workload != None:
    workloads = args.workload
  else:
    workloads=('exchange','dapps','msnfs','msncfs', 'radius')
    #workloads=('exchange','dapps','msncfs','msnfs')

  print 'Analyzing following workloads', workloads

  #writecsv('dump_stats',(('Workload','Intensity','Snapshot (months)','Mean', 'Median', 'Variance'),))
  for workload in workloads:
    #analyze_by_intensity(trace_root,workload)
    analyze_by_time(trace_root,workload)
    #get_first_snapshot_with_min_lifetime(trace_root,workload)

  colorama.deinit()


def analyze_by_intensity(trace_root,workload):

  prefix=None
  if filename!=None and len(filename)>0:
    prefix = []
    snapshot_series = []
    for fn in filename:
      bn = os.path.basename(fn)
      snapshot_id = get_snapshot_id(fn) * snapshot_frequency
      if len(bn.split('_')) >= 2:
        intensity = int(bn.split('_')[1].split('x')[0])
        #prefix.append('extrapolated_%.1f'%snapshot_id)
        prefix.append(os.path.basename(os.path.dirname(fn)))
      else:
        dirname = os.path.basename(os.path.dirname(fn))
        intensity = int(dirname.split('_')[2])
        prefix.append('detailed_%.1f' %snapshot_id)

      snapshot_records = get_all_records_in_file(fn,intensity)
      if snapshot_records != None:
        snapshot_series.append(snapshot_records)
    compare_and_plot_across_intensities(snapshot_series,workload,snapshot_id,prefix=prefix)
  else:
    snapshot_dir =  os.path.join(trace_root,workload)
    dirs = get_all_dirs(snapshot_dir,'*disksim_snapshot_*180')
    print dirs

    start_index = int(start_year * 12 * 1/snapshot_frequency)
    limit = start_index+int(num_years*12/snapshot_frequency)+1
    snapshot_files = []
    snapshot_files = [ ('snapshot%d' % i) for i in np.arange(start_index,limit,scale) ]
    #snapshot_files.append('snapshot1')
    print snapshot_files
    #snapshot_files = ('snapshot12','snapshot24','snapshot36','snapshot48','snapshot60','snapshot72','snapshot84','snapshot96','snapshot108','snapshot120')
    #snapshot_files = ('snapshot2','snapshot4','snapshot6','snapshot8','snapshot10') #intervals that we would like to analyze
    #snapshot_files = ('snapshot96',) #intervals that we would like to analyze

    for snapshots in snapshot_files:
      snapshot_series = []
      gc.collect()
      snapshot_id = get_snapshot_id(snapshots) * snapshot_frequency
      for ds_snapshot_directory in dirs:
        intensity = int(ds_snapshot_directory.split('_')[2])
        if intensity < min_intensity or intensity >max_intensity:
          continue;
        regex = snapshots+'.outv'
        #print 'Checking', ds_snapshot_directory, 'for', regex
        ds_files = get_all_files(ds_snapshot_directory,regex)
        if len(ds_files) != 0: #if the file exists
          assert len(ds_files) == 1, 'More than one %s exists' % snapshots
          ds_file = ds_files[0]
          snapshot_records = get_all_records_in_file(ds_file,intensity)
          if snapshot_records != None:
            snapshot_series.append(snapshot_records)
            print >>sys.stderr,Fore.GREEN+'Successfully read', ds_file

      compare_and_plot_across_intensities(snapshot_series,workload,snapshot_id,prefix=prefix)
    #end of for loop.

def compare_and_plot_across_intensities(snapshot_series,workload,snapshot_id,prefix=None):
  #intensity_pairs = ((1,2),(2,3),(3,4),(4,5))
  intensity_pairs = ((1,2,3,4,5,6,10,15,20,max_intensity),)
  for intensity_pair in intensity_pairs:
    compare_series = []
    for snapshots in snapshot_series:
      if snapshots.intensity in intensity_pair:
        compare_series.append(snapshots)

    if len(compare_series) >= 1:
      stats = None
      #print 'Plotting intensities', intensity_pair
      #stats = compute_stats(snapshot_series,get_logical_stress)
      if plot_retention == True:
        #stats = compute_rmse(snapshot_series)
        #stats+= compute_rsquared(snapshot_series)
        plot_intensities(compare_series,workload,snapshot_id,get_retention_period,stats,prefix=prefix)
      elif plot_stress == True:
        plot_intensities(compare_series,workload,snapshot_id,get_logical_stress,stats,prefix=prefix)
      elif plot_eqn == True:
        plot_intensities(compare_series,workload,snapshot_id,get_eqn_cycle,stats,prefix=prefix)
      else:
        plot_intensities(compare_series,workload,snapshot_id,get_delvth,stats,prefix=prefix)
    else:
      print 'Unable to find files for intensities', intensity_pair

def compute_rsquared(snapshot_series):
  stat_text = None
  if snapshot_series!=None and len(snapshot_series) == 2:
    stat_text = []
    retention1 = get_retention_period(snapshot_series[0])
    retention2 = get_retention_period(snapshot_series[1])
    sserr= sum( [math.pow((y1 -y2),2) for y1,y2 in zip(retention1,retention2)] )
    average = sum(retention2)/len(retention2) #second set of points corresponds to measured values ->detailed simulations.
    sstot = sum ([math.pow((r2-average),2) for r2 in retention2])
    r_square = 1 - (sserr/sstot)
    stat_text.append('R^2=%.2f'% r_square)

  return stat_text

def compute_rmse(snapshot_series):
  stat_text = None
  if snapshot_series!=None and len(snapshot_series) == 2:
    stat_text = []
    retention1 = get_retention_period(snapshot_series[0])
    retention2 = get_retention_period(snapshot_series[1])
    mse = 0
    for r1,r2 in zip(retention1,retention2):
      mse += math.pow((r1-r2),2)
    mse = mse/min(len(retention1),len(retention2))
    stat_text.append('RMSE=%.2f'% (math.sqrt(mse)))

  return stat_text
    
def compute_stats(snapshot_series,compute_func):
  statistics = []
  for snapshot in snapshot_series:
    sample = compute_func(snapshot)
    stat_text = mean_variance_median(sample) # normaltest(sample) #chisquaretest(sample) 
    statistics.append(stat_text)
  return statistics

def mean_variance_median(sample):
  return 'mean=%.2f,var=%.2E,med=%.2f,mode=%.2E'% (np.mean(sample),np.var(sample),np.median(sample),stats.mode(sample)[0])

def chisquaretest(sample):
  np_sample = np.array(sample)
  #bins = int(math.sqrt(len(sample)))
  bins = 100
  histos = stats.histogram(np_sample,bins)
  (chi_sq,p) = stats.chisquare(filter(lambda x:x>5,histos[0]))
  dof = (histos[0]>5).sum() - 3
  print filter(lambda x:x>5,histos[0])
  #pdb.set_trace()
  return 'chi_square=%.2f,p=%.4f,dof=%d' % (chi_sq,p,dof)

def normaltest(sample):
  return 'normaltest = %6.3f pvalue = %6.4f' % stats.normaltest(sample)

def chunks(l,n=1):
  num_chunks = len(l)/n
  return [l[i:i+num_chunks] for i in range(0,len(l),num_chunks)]

def plot_intensities(snapshot_series,workload,duration,compute_func,statistics=None,prefix=None):
 
  samples = []
  legends = []
  output_file = '_'.join([workload,str(int(duration))])
  if save_table == True:
    csv_files = []
    csv_output_file = os.path.join(out_dir,'_'.join([workload,str(duration)]))
  j = 0
  for snapshot in snapshot_series:
    gc.collect()
   
    sample = compute_func(snapshot)
    #split sample into 16 chunks corresponding to each chip.
    num_chunks = 1 #16 *8
    chunk_sample = chunks(sample,num_chunks)
    #chunk_sample = chunk_sample[11:16]
    for i in range(0,len(chunk_sample)): 
      samples.append(chunk_sample[i])
      legend_str = str(snapshot.intensity)+'x'
      if prefix!=None and len(prefix)>j:
        legend_str = legend_str + 'x_' + prefix[j]
      legends.append(legend_str)
      output_file = '_'.join([output_file,str(int(snapshot.intensity))])
    j+=1
    if save_table == True:
      csv_files.append('_'.join([csv_output_file,str(int(snapshot.intensity))])+'.csv')

  if save_table == True:
    write_multi_csvs(samples,csv_files)
    return

  annotations = []
  if statistics!=None:
    for stat_text,snapshot in zip(statistics,snapshot_series):
      annotations.append(stat_text)    


  #plot them first.
  
  plt.figure()
  sub_plot = plt.subplot(111)
  

  if normalize == True:
    sub_plot.set_ylabel("% of blocks")
  else:
    sub_plot.set_ylabel("# of blocks")

  if plot_retention == True:
    string = 'retention period'
    sub_plot.set_xlabel('Retention period (years)')
    output_file = output_file+'_r'
  elif plot_stress == True:
    string = 'P/E cycles'
    sub_plot.set_xlabel('P/E cycles')
    output_file = output_file+'_pe'
  elif plot_eqn == True:
    string = 'Equivalent cycles'
    output_file = output_file+'_eqn'
  elif plot_delvth == True:
    string = 'delvth'
    sub_plot.set_xlabel('delvth')
    output_file = output_file+'_vth'

  if normalize == True:
    output_file = output_file+'_norm'

  if histogram == False:
    output_file = output_file+'_cdf.%s' %fmt
  else:
    output_file = output_file+'.%s' %fmt

  if histogram == False:
    title = ('CDF of %s for various intensities for ' + get_label(workload) +' after ' + str(duration) +' months') % (string)
    sub_plot.set_title(title,fontproperties= FontProperties(size='small'))
    plot_cdf(samples,output_file,legends,sub_plot,normalize,annotations)
  else:
    title = ('Histogram of %s for various intensities for ' + get_label(workload) +' after ' + str(duration/12) +' years') % (string)
    print annotations
    sub_plot.set_title(title,fontproperties= FontProperties(size='small'))
    plot_histogram(samples,output_file,legends,sub_plot,normalize,annotations)


def write_multi_csvs(records,csv_files):
  if csv_files!=None and len(csv_files)>0:
    assert len(csv_files) == len(records), 'Error in dumping csv files'
    for csv_file,record in zip(csv_files,records):
      writecsv(csv_file,record)
      print 'csv filed saved in', os.path.abspath(csv_file)

def get_logical_stress(snapshotrecord):
  return snapshotrecord.get_logical_stress()

def get_eqn_cycle(snapshotrecord):
 return snapshotrecord.get_eqn_cycle()

def get_delvth(snapshotrecord):
 return snapshotrecord.get_delvth()

def get_retention_period(snapshotrecord):
  return snapshotrecord.get_retention_period()

def find_avg_cycling(si):
  avg = []
  for snapshot_records in si:
    avg.append(snapshot_records.compute_avg_cycling())
  return map(lambda y: y.logical_stress,avg)

#only return blocks with less than specified lifetime
def get_retention_period_under_limit(snapshotrecord):
  ret = get_retention_period(snapshotrecord)
  print 'before filtering', len(ret)
  filtered_ret = filter( lambda x: x<=max_retention, ret)
  print 'after filtering', len(filtered_ret)
  return filtered_ret

def get_least_retention_record(snapshotrecord):
  ret = get_retention_period(snapshotrecord)
  return min(ret)
  
def get_first_snapshot_with_min_lifetime(trace_root,workload):
  snapshot_dir =  os.path.join(trace_root,workload)
  dirs = get_all_dirs(snapshot_dir,'disksim_snapshot_*')
  dirs.sort(reverse=True)
  print dirs

  cache_file_name = 'cache_snapshot_'+workload
  snapshot_series = load_cache(cache_file_name)
  
  if snapshot_series == None:
    strs = []
    for ds_snapshot_directory in dirs: #for every intensity
      intensity = int(ds_snapshot_directory.split('_')[2])
      if intensity < min_intensity:
        continue;
      ds_files = get_all_files(ds_snapshot_directory,'snapshot*')
      ds_files = filter_files(ds_files)
      for fname in ds_files:
        gc.collect()
        records = get_all_records_in_file(fname,intensity)
        if records!=None:
          xvalues = get_retention_period_under_limit(records)
          if xvalues!=None and len(xvalues)>0:
            strs.append('First snapshot with retention period less than %f year is %s' % (max_retention,fname))
            break;
    if len(strs)>0:
      writecsv('start_snapshot.csv',strs,mode='a')

def analyze_by_time(trace_root,workload):

  snapshot_dir =  os.path.join(trace_root,workload)
  dirs = get_all_dirs(snapshot_dir,'disksim_snapshot_*')
  print dirs

  cache_file_name = 'cache_snapshot_'+workload
  snapshot_series = load_cache(cache_file_name)
  
  if snapshot_series == None:
    for ds_snapshot_directory in dirs: #for every intensity
      snapshot_series = []
      gc.collect()
      intensity = int(ds_snapshot_directory.split('_')[2])
      if intensity < min_intensity:
        continue;
      ds_files = get_all_files(ds_snapshot_directory,'snapshot*')
      ds_files = filter_files(ds_files)

      print 'Processing', ds_snapshot_directory
      snapshot_records = get_all_records_in_dir(ds_snapshot_directory,ds_files,intensity)
      snapshot_series.append(snapshot_records)
      print 'Finishing getting all the file'
		# snapshot_series: [[SnapshotRecords1_intensity1,SnapshotRecords2_intensity1...],SnapshotRecords1_intensity2,SnapshotRecords2_intensity2...]... ]
		# python has object-oriented features

      #output_file='lifetime_'+workload+'_'+str(intensity)+'x_average.'+fmt
      #plot_series(snapshot_series,workload,get_avg_lifetime_series,output_file)

      #output_file='lifetime_'+workload+'_'+str(intensity)+'x_min.'+fmt
      #plot_series(snapshot_series,workload,get_min_lifetime_series,output_file)

      if histogram == False:
        if normalize == True:
          output_file = workload+'_'+str(intensity)+'x_'+str(max_retention)+'_cdf.'+fmt
        else:
          output_file = workload+'_'+str(intensity)+'x_'+str(max_retention)+'_no_norm_cdf.'+fmt
      else:
        output_file = workload+'_'+str(intensity)+'x_'+str(max_retention)+'_hist.'+fmt
      if histogram  == True:
        compare_and_plot_across_time(snapshot_series,workload,intensity,output_file,get_retention_period_under_limit)
      if hellinger_discrete == True:
        output_file = workload+'_'+str(intensity)+'x_'+str(max_retention)+'_hellinger_discrete.'+fmt
      if hellinger_normalize == True:
	output_file = workload+'_'+str(intensity)+'x_'+str(max_retention)+'_hellinger_normalize.'+fmt
      if earth_mover == True:
	output_file = workload+'_'+str(intensity)+'x_'+str(max_retention)+'_earth_mover.'+fmt
        compute_distance(snapshot_series, workload, intensity, output_file, get_retention_period_under_limit)
  
      if analyze_err == True:
        output_file = workload+'_'+str(intensity)+'x_'+str(max_retention)+'_error_percentage.'+fmt
        compute_error_and_time_saving(snapshot_series, workload, intensity, output_file, get_retention_period_under_limit)
      
# snapshot_series contains all the information extracted from snapshots;

def filter_files(ds_files):
  new_ds_files = []
  scale = 2 #snapshot corresponding to every snapshot_frequency * scale month are analyzed
  scale = int(snapshot_frequency*scale) 
  start_index = int(start_year * 12 * 1/snapshot_frequency)
  num_files = int(num_years * 12 * 1/snapshot_frequency)
  limit = min(num_files,len(ds_files))
  ds_files = ds_files[start_index:start_index+limit] #limit the analysis to these files.
  new_ds_files = ds_files[::scale] #pick every scale file.
  if (len(new_ds_files)<1 or new_ds_files[-1] != ds_files[-1]):
    new_ds_files.append(ds_files[-1])
  return new_ds_files

def compare_and_plot_across_time(snapshot_series,workload,intensity,output_file,compute_func):

  for si in snapshot_series: #corresponds to every intensity
    xvalues = []
    legends = []
    for snapshots in si: #corresponds to different time periods
      print 'Processing ', snapshots.filename 
      print 'compute_func in compare_and_plot_across_time'
      xvalue = compute_func(snapshots)
      if xvalue!=None and len(xvalue)!=0:      
        xvalues.append(xvalue)
        snapshot_id = get_snapshot_id(snapshots.filename)*snapshot_frequency
        legends.append('%.1f' % snapshot_id);

    if len(legends)>0:
      title = 'plot of delvth at different months for ' + workload + ' - '+str(intensity) + 'x'
      plt.figure() # create a new figure
      sub_plot = plt.subplot(111)
      sub_plot.set_xlabel('delvth')

      if normalize == True:
        sub_plot.set_ylabel("% of blocks")
      else:
        sub_plot.set_ylabel("# of blocks")

      sub_plot.set_title(title,fontproperties= FontProperties(size='small'))
      print Fore.GREEN+'Plotting %s' % workload

      if histogram == False:
        plot_cdf(xvalues,output_file,legends,sub_plot,normalized=normalize)
      else:
        plot_histogram(xvalues,output_file,legends,sub_plot,normalized=normalize)

def plot_earth_mover(x, output_file, sub_plot):
#  earth_mover_dist = []
  if analyze_err == True:
    cost_last = []
    [count_last, bin_last] = np.histogram(x[len(x)-1], bins=3000, normed=True)
    for i in range(len(bin_last)-1):
      cost_last.append((bin_last[i]+bin_last[i+1])/2)   
    sig_last = cv.CreateMat(len(count_last), 2, cv.CV_32FC1)
    for j in range(len(count_last)):
      sig_last[j, 0] = count_last[j]
      sig_last[j, 1] = cost_last[j]

  for i in range(len(x)-1):
    count1 = []
    count2 = []
    bin1 = []
    bin2 = []
    [count1, bin1] = np.histogram(x[i], bins=3000, normed=True)
    [count2, bin2] = np.histogram(x[i+1], bins=3000, normed=True)
   
    # get cost1 and cost2
    cost1 = []
    cost2 = []
    gc.collect()
    for s in range(len(bin1)-1):
      cost1.append((bin1[s]+bin1[s+1])/2)
      cost2.append((bin2[s]+bin2[s+1])/2)

     # get sig1 and sig2
    sig1 = cv.CreateMat(len(count1), 2, cv.CV_32FC1)
    sig2 = cv.CreateMat(len(count2), 2, cv.CV_32FC1)
    gc.collect()
    # sig1 and sig2 has 2 columns
    for j in range(len(count1)):
      sig1[j,0] = count1[j]
      sig1[j,1] = cost1[j]
      sig2[j,0] = count2[j]
      sig2[j,1] = cost2[j]
    # compute cost matrix
 #   costs = cv.CreateMat(len(count1), len(count2), cv.CV_32FC1)
 #   for m in range(len(count1)):
    #  cost = []
  #    for n in range(len(count2)):
  #      costs[m, n] = abs(cost1[m] - cost2[n])
  #     cost.append(math.fabs(cost1[m] - cost2[n]))
  #     costs.append(cost)
    # compute EMD
    dist = cv.CalcEMD2(sig1, sig2, cv.CV_DIST_L1)
    earth_mover_dist.append(dist)
    if analyze_err == True:
      dist_to_last = cv.CalcEMD2(sig1, sig_last, cv.CV_DIST_L1)
      err_in_dist.append(dist_to_last)
  
  # plot the dist curve
  sub_plot.plot(range(1, len(earth_mover_dist)+1), earth_mover_dist, marker='o')
#  plt.ylim([0, 4])
  output_path = os.path.join(out_dir, output_file)
  print 'Saving the graph in ', os.path.abspath(output_path)
  plt.savefig(output_path) 


def compute_distance(snapshot_series,workload,intensity,output_file,compute_func):
  xvalues = []

  for si in snapshot_series:  
    for snapshots in si:
      print 'Processing ', snapshots.filename
      print 'compute distance'
      xvalue = compute_func(snapshots)
      if xvalue!= None and len(xvalue)>0:
        xvalues.append(xvalue)
    
  title = 'Distance Plot of delvth curves for' + workload +' - '+str(intensity) + 'x'
  plt.figure()
  sub_plot = plt.subplot(111)
  sub_plot.set_xlabel('snapshot ID')
  sub_plot.set_ylabel('Distance of Curve')

  sub_plot.set_title(title, fontproperties=FontProperties(size='small'))
  print Fore.GREEN+'Plotting %s' % workload    

  if hellinger_discrete == True:
    plot_hellinger_discrete(xvalues, output_file, sub_plot)
  if hellinger_normalize == True:
    plot_hellinger_normalize(xvalues, output_file, sub_plot)
  if earth_mover == True:
    plot_earth_mover(xvalues, output_file, sub_plot)

def compute_error_and_time_saving(snapshot_series, workload, intensity, output_file, compute_func):

  title = 'Error Percentage of Snapshots for ' + workload +' - '+str(intensity) + 'x'
  fig = plt.figure()
  sub_plot = fig.add_subplot(211)
  sub_plot_2 = fig.add_subplot(212)
  sub_plot.set_xlabel('Error Allowed in Distance between neighboring snapshots')
  sub_plot.set_ylabel('Total Error')
  sub_plot_2.set_xlabel('Time Saving Percentage')
  sub_plot_2.set_ylabel('Total Error percentage')
  sub_plot.set_title(title, fontproperties=FontProperties(size='small'))
  print Fore.GREEN+'Plotting %s' % workload
  plot_error_percentage(output_file, sub_plot, sub_plot_2)

def plot_error_percentage(output_file, sub_plot, sub_plot_2):

  legend = []
  legend.append('error')
  legend.append('time saving percentage')

  # compute time saving
#  snapshot_to_skip = [x1 for x1 in range(0, len(err_in_dist))]
#  time_saving = [(((len(err_in_dist))-x1)/float(len(err_in_dist)))  for x1 in snapshot_to_skip]

#  sub_plot.plot(earth_mover_dist, err_in_dist, label=legend[0], marker='o', color='r', linestyle=' ')
#  sub_plot.legend(loc='upper left', shadow=True, fancybox=True, ncol=5,prop= FontProperties(size='small'))
#  plt.ylim([0, 12])
#  sub_plot_2.plot(earth_mover_dist, time_saving, label=legend[1], marker='o', linestyle=' ')
#  sub_plot_2.legend(loc='upper left', shadow=True, fancybox=True, ncol=20,prop= FontProperties(size='small'))
#  plt.ylim([0, 1])

  max_err = []
  time_saving = []
  neighbor_dist_allowed = []
  max_err_percentage = []
  
  i = 0
  while i < len(earth_mover_dist):
    for j in range(i+1, len(earth_mover_dist)):
      if earth_mover_dist[j] < earth_mover_dist[i]:
        break
    neighbor_dist_allowed.append(earth_mover_dist[i])
    max_err.append(err_in_dist[j])
    max_err_percentage.append(err_in_dist[j]/float(err_in_dist[0]))
    time_saving.append((len(earth_mover_dist) - j)/float(len(earth_mover_dist) + 1))
    if j == (len(earth_mover_dist) - 1):
      break
    i = j

  neighbor_dist_allowed.append(0)
  max_err.append(0)
  time_saving.append(0)
  max_err_percentage.append(0)
  sub_plot.plot(neighbor_dist_allowed, max_err, label=legend[0], marker='o', color='g', linestyle='--')
  sub_plot.legend(loc='upper left', shadow=True, fancybox=True, ncol=5,prop= FontProperties(size='small'))
 # sub_plot.set_ylim(0, 10)  

  sub_plot2 = sub_plot.twinx()
  sub_plot2.plot(neighbor_dist_allowed, time_saving, label=legend[1], marker='o', linestyle='--' )  
  sub_plot2.legend(loc='upper right', shadow=True, fancybox=True, ncol=5,prop= FontProperties(size='small'))
  sub_plot2.set_ylim(0, 1)
  
 
 # title = 'Time Saving vs Total Error'
 # sub_plot_2.set_title(title, fontproperties=FontProperties(size='small'))
  sub_plot_2.plot(time_saving, max_err_percentage, marker='o', color='r', linestyle='--')
  sub_plot_2.set_xlim(0, 1)
  sub_plot_2.set_ylim(0, 1)

  output_path = os.path.join(out_dir, output_file)
  print 'Saving the graph in', os.path.abspath(output_path)
  plt.savefig(output_path)
  
  
	
def plot_hellinger_discrete(x, output_file, sub_plot):
  hellinger_distance = []
  count_all = []
  
  for i in range(len(x)):
    (counts, bins) = np.histogram(x[i], bins = 3000, range=(0.0, 30.0), normed = True)
    count_all.append(counts)
#	bins = bins[:len(counts)]
#    counts = []
#    for j in range(0, 300):
 #     count_interval = 0
 #     for s in range(len(x[i])):
  #      if x[i][s] >= 0.1 * j and x[i][s] < 0.1 * (j+1):
  #        count_interval = count_interval+1
  #    counts.append(count_interval)

 #   counts = [float(x1)*100/sum(counts) for x1 in counts]


  for i in range(0, len(count_all)-1):
    sum_of_all = 0
    for j in range(len(count_all[i])):
      sum_of_all = sum_of_all + (math.sqrt(count_all[i+1][j])-math.sqrt(count_all[i][j]))*(math.sqrt(count_all[i+1][j])-math.sqrt(count_all[i][j]))
    dist = math.sqrt(sum_of_all)
    hellinger_distance.append(dist)
  sub_plot.plot(range(1, len(count_all)), hellinger_distance)
  output_path = os.path.join(out_dir, output_file)
  print 'Saving the graph in', os.path.abspath(output_path)
  plt.savefig(output_path)

def plot_series(snapshot_series,workload,compute_func,output_file):

  points = []
  xvalues = []
  legends = []  

  for si in snapshot_series:
    xvalue = create_xvalues(map(lambda x: x.filename, si))  
    yvalue = compute_func(si)
    points.append(yvalue)
    xvalues.append(xvalue)
    print len(yvalue),len(xvalue),xvalue,yvalue
    legends.append('%dx' % si[0].intensity)

  title = 'Disksim detailed evaluation for ' + workload 
  plt.figure()
  sub_plot = plt.subplot(111)
  sub_plot.set_xlabel('Used Lifetime(months)')
  sub_plot.set_ylabel("Remaining lifetime(years)")
  sub_plot.set_title(title)

  print Fore.GREEN+'Plotting %s' % workload
  plot_lifetime_workloads(xvalues,points,output_file,legends,sub_plot)


def get_all_records_in_dir(directory,input_files,intensity):
  if input_files == None:
    input_files = get_all_files(directory,'snapshot*')

  records = []
  for fname in input_files:
    snapshot_records = get_all_records_in_file(fname,intensity)
    if snapshot_records != None:
      records.append(snapshot_records)
  return records

def get_all_records_in_file(fname,intensity):
  #print 'Processing ', fname
  if os.path.isfile(fname) == True and os.path.getsize(fname) > 0:
    ta_snapshot=True
    if fnmatch.fnmatch(fname,'*disksim_snapshot*') == True:
      ta_snapshot=False
    ss_file = SnapshotFile(fname,is_ta_snapshot=ta_snapshot)
    snapshot_records = SnapshotRecords(ss_file.readrecords(),fname,intensity)
    return snapshot_records
  else:
    print >>sys.stderr,Fore.RED+'Invalid file', fname
    return None

def get_avg_lifetime_series(si):
  avg = []
  for snapshot_records in si:
    avg.append(snapshot_records.compute_avg_retention_block())

  return map(lambda y: y.retention_period/year,avg)

def get_min_lifetime_series(si):
  min_blocks = []
  for snapshot_records in si:
    min_blocks.append(snapshot_records.find_min_retention_block())

  return map(lambda y: y.retention_period/year, min_blocks)

def create_xvalues(files):
  xvalue = []
  for eachfile in files:
    snapshot_id = get_snapshot_id(eachfile)
    xvalue.append(snapshot_id*snapshot_frequency)
  return xvalue

def get_snapshot_id(filename):
  snapshot_id = os.path.splitext(os.path.basename(filename))[0]
  snapshot_id = int(re.split('snapshot',snapshot_id)[1])
  return snapshot_id
  
def load_cache(filename):
  if os.path.isfile(filename):
    print 'Loading data from ',filename
    fp = open(filename,"rb")
    data = pickle.load(fp)
    fp.close()
    return data
  else:
    return None

def save_cache(data,filename):
  if data != None and len(data) >0:
    fp = open(filename,"wb")
    pickle.dump(data,fp)
    fp.close()

def plot_lifetime_workloads(x,y,output_file,legends,sub_plot=None):
  plot_objs = []
  legend_text = []

  for i in range(len(legends)):
    print 'Plotting', legends[i]
    legend_text.append(legends[i])
    if y[i]!= None and len(y[i])>0:
      plot_objs.append(sub_plot.plot(x[i],y[i]))

  plt.legend(plot_objs, legend_text, loc='upper right', shadow=True, fancybox=True,
             ncol=5,prop= FontProperties(size='small'))
  plt.legend(plot_objs)

  output_path = os.path.join(out_dir,output_file)
  print 'Saving the graph in', os.path.abspath(output_path)
  plt.savefig(output_path)


def plot_cdf(x,output_file,legends,sub_plot,normalized=False,annotations=None):

  if len(x)>0:
    for i in range(len(legends)):
      print 'Plotting', legends[i]
      (counts,edges) = np.histogram(x[i],bins=int(math.sqrt(len(x[i]))),normed=normalized)
      edges = edges[:len(counts)] #sometimes the size of edges and counts are not the same.
      yvalue = np.cumsum(counts) 

      if normalized == True:
        scale = 1/yvalue[-1]
        ncdf = (scale * yvalue * 100)
        yvalue = ncdf

      #print yvalue,edges
      sub_plot.plot(edges,yvalue,label=legends[i])

    x_start = 0.8
    y_start = 0.9
    x_off = 0.0
    y_off = 0.05
    x_o = x_start
    y_o = y_start
    if annotations != None and len(annotations)>0:
      for text_str in annotations:
        print text_str
        x_o += x_off
        y_o -= y_off
        sub_plot.text(x_o,y_o,text_str,transform = sub_plot.transAxes,fontproperties= FontProperties(size='small'))


    plt.legend(loc='upper left', shadow=True, fancybox=True,
             ncol=5,prop= FontProperties(size='small'))

    output_path = os.path.join(out_dir,output_file)
    print 'Saving the graph in', os.path.abspath(output_path)
    plt.savefig(output_path)

def plot_histogram(x,output_file,legends,sub_plot,normalized=False,annotations=None):

  if len(x)>0:
    for i in range(len(legends)):
      print 'Plotting', legends[i],
      #n,bins,rectangles = sub_plot.hist(x[i],bins=300,label=legends[i],normed=normalized)
      (counts,bins) = np.histogram(x[i],bins=int(math.sqrt(len(x[i]))),normed=False)
      bins = bins[:len(counts)]
      if normalized == True:
        counts = [float(x1)*100/sum(counts) for x1 in counts]
      sub_plot.plot(bins,counts,label=legends[i])
     
    x_start = 0.8
    y_start = 0.9
    x_off = 0.0
    y_off = 0.05
    x_o = x_start
    y_o = y_start
   # CHANGED: changed
    if annotations != None and len(annotations)>0:
      for text_str in annotations:
        print text_str
        x_o += x_off
        y_o -= y_off
        sub_plot.text(x_o,y_o,text_str,transform = sub_plot.transAxes,fontproperties= FontProperties(size='small'))
    
#    plt.legend(loc='upper left', shadow=True, fancybox=True,
 #            ncol=5,prop= FontProperties(size='small'))
    if plot_retention == True and normalized == False:
      plt.xlim(xmin=-0.5)
      #plt.ylim(ymin=-100)
    
    output_path = os.path.join(out_dir,output_file)
    print 'Saving the graph in', os.path.abspath(output_path)
    plt.savefig(output_path)

if __name__ == "__main__":
    sys.exit(main())

