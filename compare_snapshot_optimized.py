#!/usr/bin/env python

# compare the mean, std_dev and earth_mover_dist

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

import sys
sys.path.append("/uf21/lj5bp/py_lib")
sys.path.append("/uf21/lj5bp/Documents/research/opencv/lib/python2.6/dist-packages")
sys.path.append("/uf21/lj5bp/Documents/research/opencv/lib")

import os,csv,re,math
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
out_dir = ''

def main(argv=None):
  colorama.init(autoreset=True)
  
  # handle some args
  if argv is None:
    argv = sys.argv

  parser = argparse.ArgumentParser(description='Script to compare the difference between snapshots of two workloads')
  parser.add_argument('-d', '--detailed', required=True, help='detailed simulation snapshots directory')
  parser.add_argument('-f', '--fastforward', required=True, help='fastforward simulation snapshots directory')
  parser.add_argument('-o', '--out_dir', default='.',help='Output directory to save figures')
  parser.add_argument('-w', '--workload', default=None,help='Name of workload to be analyzed')
  parser.add_argument('-i', '--intensity',type=float,default=1,help='intensity to be considered')
  

  args = parser.parse_args()
  detailed_dir = os.path.abspath(args.detailed)
  fastforward_dir = os.path.abspath(args.fastforward)

  global fmt
  fmt = 'png'
  
  global out_dir
  out_dir = args.out_dir
  mkdir_p(out_dir)

  global max_retention
  max_retention = 35.0

  global workload
  workload = args.workload
  
  global intensity
  intensity = args.intensity
  
  print 'Processing ', detailed_dir
  # get the records of snapshots under dir1 and dir2 in 2 dim array
  detailed_files = get_all_files(detailed_dir, 'snapshot*')
#  detailed_snapshot_records = get_all_records_in_dir(detailed_dir, detailed_files, intensity)

  # optimize memory cost
  print 'Processing ', fastforward_dir
  fastforward_files = get_all_files(fastforward_dir, 'snapshot*')
#  fastforward_snapshot_records = get_all_records_in_dir(fastforward_dir, fastforward_files, intensity)

  # records is [[blk1_ss1, blk2_ss2, ..., blkn_ssn], ...]

  # get the retention time of the extracted records
#  detailed_retention = extract_retention_from_records(detailed_snapshot_records)
#  fastforward_retention = extract_retention_from_records(fastforward_snapshot_records)


  # compute mean and std-dev with the array
#  plot_mean(detailed_retention, fastforward_retention)
#  plot_stddev(detailed_retention, fastforward_retention)
  plot_mean(detailed_files, fastforward_files)
  plot_stddev(detailed_files, fastforward_files)
 
  # compute and plot the earth mover's distance diffs
  plot_earth_mover_distance(detailed_files, fastforward_files)

  colorama.deinit()
  
def get_all_records_in_dir(directory, input_files, intensity):
  if input_files == None:
    input_files = get_all_files(directory,'snapshot*')

  records = []
  for fname in input_files:
    snapshot_records = get_all_records_in_file(fname, intensity)
    if snapshot_records != None:
      records.append(snapshot_records)
  return records

def get_all_records_in_file(fname,intensity):
  #print 'Processing ', fname
  if os.path.isfile(fname) == True and os.path.getsize(fname) > 0:
    ta_snapshot=True
#    if fnmatch.fnmatch(fname,'*disksim_snapshot*') == True:
    ta_snapshot=False
    ss_file = SnapshotFile(fname, is_ta_snapshot=ta_snapshot)
    snapshot_records = SnapshotRecords(ss_file.readrecords(),fname,intensity)
    return snapshot_records
  else:
    print >>sys.stderr,Fore.RED+'Invalid file', fname
    return None

# return retention time array [[blk1_ss1, blk2_ss2, ..., blkn_ssn], ...]
def extract_retention_from_records(dir_record):
  
  dir_retention = []
  for snapshot_record in dir_record:
    snapshot_retention = get_retention_period_under_limit(snapshot_record)
    if snapshot_retention != None:
      dir_retention.append(snapshot_retention) 
  return dir_retention

#only return blocks with less than specified lifetime
def get_retention_period_under_limit(snapshotrecord):
  ret = get_retention_period(snapshotrecord)
  return filter( lambda x: x<=max_retention, ret)

def get_retention_period(snapshotrecord):
  return snapshotrecord.get_retention_period()


def plot_mean(detailed_files, fastforward_files):
  title = 'Plot for mean of the Retention time for ' + workload + ' - ' + str(intensity)
  fig = plt.figure()
  sub_plot = fig.add_subplot(111)
  sub_plot.set_xlabel('Snapshot ID')
  sub_plot.set_ylabel('Mean of retention time')
  sub_plot.set_title(title, fontproperties=FontProperties(size='small'))
  print 'Plotting %s' % workload
  
  legend = []
  legend.append('detailed_mode')
  legend.append('fastforward_mode')

  detailed_mean = []
  fastforward_mean = []

  # calculate the mean
  for detailed_file in detailed_files:
    detailed_record = []
    detailed_retention = []
    gc.collect()
    detailed_record = get_all_records_in_file(detailed_file, intensity)
    detailed_retention = get_retention_period_under_limit(detailed_record)
    detailed_mean.append(np.mean(detailed_retention))

  for fastforward_file in fastforward_files:
    fastforward_record = []
    fastforward_retention = []
    gc.collect()
    fastforward_record = get_all_records_in_file(fastforward_file, intensity)
    fastforward_retention = get_retention_period_under_limit(fastforward_record)
    fastforward_mean.append(np.mean(fastforward_retention))

  sub_plot.plot(range(1, len(detailed_mean)+1), detailed_mean, label=legend[0], marker='*',linestyle='--')
  sub_plot.plot(range(1, len(fastforward_mean)+1), fastforward_mean, label=legend[1], marker='*', linestyle='--')
  sub_plot.legend(loc='upper left', shadow=True, fancybox=True, ncol=5,prop= FontProperties(size='small'))
  output_file = 'mean_plot.' + fmt 
  output_path = os.path.join(out_dir, output_file)
  print 'Saving the graph in', os.path.abspath(output_path)
  plt.savefig(output_path)

def plot_stddev(detailed_files, fastforward_files):
  title = 'Plot for standard deviation of the Retention time for ' + workload + ' - ' + str(intensity)
  fig = plt.figure()
  sub_plot = fig.add_subplot(111)
  sub_plot.set_xlabel('Snapshot ID')
  sub_plot.set_ylabel('std-dev of retention time')
  sub_plot.set_title(title, fontproperties=FontProperties(size='small'))
  print 'Plotting %s' % workload

  legend = []
  legend.append('detailed_mode')
  legend.append('fastforward_mode')

  detailed_stddev = []
  fastforward_stddev = []

  # calculate the mean
  for detailed_file in detailed_files:
    detailed_record = []
    detailed_retention = []
    gc.collect()
    detailed_record = get_all_records_in_file(detailed_file, intensity)
    detailed_retention = get_retention_period_under_limit(detailed_record)
    detailed_stddev.append(np.std(detailed_retention))

  for fastforward_file in fastforward_files:
    fastforward_record = []
    fastforward_retention = []
    gc.collect()
    fastforward_record = get_all_records_in_file(fastforward_file, intensity)
    fastforward_retention = get_retention_period_under_limit(fastforward_record)
    fastforward_stddev.append(np.std(fastforward_retention))

  sub_plot.plot(range(1, len(detailed_stddev)+1), detailed_stddev, label=legend[0], marker='*', linestyle='--')
  sub_plot.plot(range(1, len(fastforward_stddev)+1), fastforward_stddev, label=legend[1], marker='*', linestyle='--')
  sub_plot.legend(loc='upper left', shadow=True, fancybox=True, ncol=5,prop= FontProperties(size='small'))
  output_file = 'stddev_plot.' + fmt 
  output_path = os.path.join(out_dir, output_file)
  print 'Saving the graph in', os.path.abspath(output_path)
  plt.savefig(output_path)

def plot_earth_mover_distance(detailed_files, fastforward_files):
  title = 'Plot for difference in earth mover distance of the Retention time for ' + workload + ' - ' + str(intensity)
  fig = plt.figure()
  sub_plot = fig.add_subplot(111)
  sub_plot.set_xlabel('Snapshot ID')
  sub_plot.set_ylabel('earth mover distance of retention time')
  sub_plot.set_title(title, fontproperties=FontProperties(size='small'))
  print 'Plotting %s' % workload

  earth_mover_dist = []
 
  for i in range(len(detailed_files)):
    detailed_file = detailed_files[i]
    fastforward_file = fastforward_files[i]
    detailed_record = []
    detailed_retention = []
#    gc.collect()
    detailed_record = get_all_records_in_file(detailed_file, intensity)
    detailed_retention = get_retention_period_under_limit(detailed_record)
    
    fastforward_record = []
    fastforward_retention = []
 #   gc.collect()
    fastforward_record = get_all_records_in_file(fastforward_file, intensity)
    fastforward_retention = get_retention_period_under_limit(fastforward_record)
    
    count1 = []
    count2 = []
    bin1 = []
    bin2 = []
#    gc.collect()
    [count1, bin1] = np.histogram(detailed_retention, bins=3000, normed=True)
    [count2, bin2] = np.histogram(fastforward_retention, bins=3000, normed=True)

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
 
    # compute EMD
    dist = cv.CalcEMD2(sig1, sig2, cv.CV_DIST_L1)
    earth_mover_dist.append(dist)


  # plot the dist curve
  sub_plot.plot(range(1, len(earth_mover_dist)+1), earth_mover_dist, marker='o')
  output_file = 'earthmoverdist_plot.' + fmt 
  output_path = os.path.join(out_dir, output_file)
  print 'Saving the graph in', os.path.abspath(output_path)
  plt.savefig(output_path)


if __name__ == "__main__":
    sys.exit(main())
  
   
  




 
  
 









