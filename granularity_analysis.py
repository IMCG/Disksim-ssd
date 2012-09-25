#!/usr/bin/env python
# By Luyao Jiang, University of Virginia

# This script is to analyze the granularity of spatial and temporal bucket as well as the smallest unit for detailed and fastforward simulation

# Usage: TBA

# initialization
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


def main(argv=None):
  colorama.init(autoreset=True)
  
  # handle some args
  if argv is None:
    argv = sys.argv

  parser = argparse.ArgumentParser(description='Script to analyze the granularity for different workloads')
  parser.add_argument('-i', '--in_dir', default='.',help='Input directory to read in all the stress distribution matrices')
  parser.add_argument('-o', '--out_dir', default='.',help='Output directory to save figures')
  parser.add_argument('-r', '--row', type=int, default=1, help='Spatial granularity')
  parser.add_argument('-c', '--column', type=int, default=1, help='Temporal granularity')
  parser.add_argument('-s', '--stride', type=int, default=1, help='Stride to combine the smallest unit for simulation')
  # might need the start index for the matrix
  

  args = parser.parse_args()
  
  global in_dir
  in_dir = os.path.abspath(args.in_dir)
  
  global out_dir
  out_dir = os.path.abspath(args.out_dir)

  global fmt
  fmt = 'png'
  
  global spatial_compress_rate
  spatial_compress_rate = args.row

  global temporal_compress_rate
  temporal_compress_rate = args.column

  global stride
  stride = args.stride

  global num_row_origin
  num_row_origin = 4096

  global num_col_origin
  num_col_origin = 48

# ########################################### beginning of granularity analysis code
  # 1 get all the file names under the in_dir
  # 2 get stress distribution array within one simulation unit (one stride)
  # 3 pre-filter each matrix getting from step 2 to the required spatial and temporal granularity
  # 4 plot spatial and temporal distribution for each simulation unit under specified granularity
  # 5 calc the distance between 2 adjacent matrices;
  print 'Getting all files under input directory'
  # get_all_files return an array containing absolute directory + filename
  matrix_files = get_all_files(in_dir, '*matrix*')
  print 'Test if the order is correct'
  print matrix_files

  distance = []
  
  last_unit = (len(matrix_files) - 1) / stride
  print 'Last unit is %d' %last_unit 
  for i in range(0, last_unit):
    prev_unit_matrix = []
    curr_unit_matrix = []
    # prev_init_matrix is the unit matrix for specified granularity including the prefiltering stage;
    # TODO:
    prev_unit_matrix = get_matrix_for_unit(matrix_files, i)
    curr_unit_matrix = get_matrix_for_unit(matrix_files, i+1)
    plot_spatial_and_temporal_distribution(prev_unit_matrix, i)
    plot_spatial_and_temporal_distribution(curr_unit_matrix, i+1)
    normalize_matrix(prev_unit_matrix)
    normalize_matrix(curr_unit_matrix)
    dist = compute_distance(prev_unit_matrix, curr_unit_matrix)
    distance.append(dist)

  # TODO: plot distance between adjacent matrices
  print 'Plot Distance between adjacent matrices'
  title = 'Plot Distance between adjacent matrices'
  fig = plt.figure()
  sub_plot = fig.add_subplot(111)
  sub_plot.set_xlabel('Matrix ID')
  sub_plot.set_ylabel('Distance')
  sub_plot.set_title(title, fontproperties=FontProperties(size='small'))
  sub_plot.plot(range(1, len(distance)+1),  distance, marker='*', linestyle='--')
  output_file = 'distance_plot' + '.' +fmt
  output_path = os.path.join(out_dir, output_file)
  print 'Saving distance graph in', os.path.abspath(output_path)
  plt.savefig(output_path)
    
def get_matrix_for_unit(matrix_files, index):
  # initialize unit matrix according to specified granularity
  num_of_row = num_row_origin / spatial_compress_rate
  num_of_col = num_col_origin / temporal_compress_rate
  unit_matrix = [([0] * num_of_col) for j in range(num_of_row)]

  for i in range(index * stride, (index+1) * stride):
    # read matrix from file
    curr_matrix = []
    # TODO:
    curr_matrix = read_matrix_from_file(matrix_files[i])
    # justify current matrix into specified granularity
    # curr_matrix_after_filtering = [([0] * num_of_col) for j in range(num_of_row)]
    for m in range(num_row_origin):
      for n in range(num_col_origin):
        unit_matrix[m/spatial_compress_rate][n/temporal_compress_rate] += curr_matrix[m][n]

  # for debugging, print out the matrix
 # for j in range(num_of_row):
 #   print unit_matrix[j]
  
  return unit_matrix

def read_matrix_from_file(file_name):
  f = open(file_name, 'r')
  # TODO: check if the last line has /n; for now, assume true
  curr_matrix = []
  for i in range(num_row_origin):
    line = f.readline()
    # excluding both ' ' and '\n' at the end of a line
    row = []
    row = line[0:-2].split(' ')
    # convert the type of elements in row to integer
    row = [int(row[i]) for i in range(len(row))]
    assert len(row) == num_col_origin
    # for debugging purpose
  #  print 'Read a row from the original matrix:'
  #  print row
    curr_matrix.append(row)

  return curr_matrix

def plot_spatial_and_temporal_distribution(unit_matrix, index):
  # plot spatial distribution
  num_row = len(unit_matrix)
  num_col = len(unit_matrix[0])
  assert num_row == num_row_origin/spatial_compress_rate
  assert num_col == num_col_origin/temporal_compress_rate
  spatial_dist = []
  for i in range(num_row):
    stress_per_row = 0
    for j in range(num_col):
      stress_per_row += unit_matrix[i][j]
    spatial_dist.append(stress_per_row)

  print 'Plot spatial distribution for index %d' %index
  title1 = 'Plot Spatial Distribution for stride ' + str(stride)
  fig1 = plt.figure()
  sub_plot1 = fig1.add_subplot(111)
  sub_plot1.set_xlabel('Block Range')
  sub_plot1.set_ylabel('Number of Stresses')
  sub_plot1.set_title(title1, fontproperties=FontProperties(size='small'))
  sub_plot1.plot([i*32*stride for i in range(num_row)],  spatial_dist, marker='*', linestyle='--')
  output_file1 = 'spatial_dist_plot_for_index_' + str(index) + '.' +fmt
  output_path1 = os.path.join(out_dir, output_file1)
  print 'Saving spatial distribution graph in', os.path.abspath(output_path1)
  plt.savefig(output_path1)

  # plot temporal distribution
  temporal_dist = []
  for i in range(num_col):
    stress_per_col = 0 
    for j in range(num_row):
      stress_per_col += unit_matrix[j][i]
    temporal_dist.append(stress_per_col)

  
  print 'Plot temporal distribution for index %d' %index
  title2 = 'Plot Temporal Distribution for stride ' + str(stride)
  fig2 = plt.figure()
  sub_plot2 = fig2.add_subplot(111)
  sub_plot2.set_xlabel('Time Range')
  sub_plot2.set_ylabel('Number of Stresses')
  sub_plot2.set_title(title1, fontproperties=FontProperties(size='small'))
  sub_plot2.plot(range(1, num_col + 1),  temporal_dist, marker='*', linestyle='--')
  output_file2 = 'temporal_dist_plot_for_index_' + str(index) + '.' +fmt
  output_path2 = os.path.join(out_dir, output_file2)
  print 'Saving temporal distribution graph in', os.path.abspath(output_path2)
  plt.savefig(output_path2)

def normalize_matrix(unit_matrix):
  print 'Normalize matrix ...'
  num_row = len(unit_matrix)
  num_col = len(unit_matrix[0])
  assert num_row == num_row_origin/spatial_compress_rate
  assert num_col == num_col_origin/temporal_compress_rate
  tot_stress = 0
  for i in range(num_row):
    for j in range(num_col):
      tot_stress += unit_matrix[i][j]
  print 'Total stresses is %d' %tot_stress
  
  for i in range(num_row):
    for j in range(num_col):
      unit_matrix[i][j] = float(unit_matrix[i][j]) / float(tot_stress)

# compute distance between two normalized matrices
def compute_distance(matrix1, matrix2):
  print 'Compute distance between two normalized matrices...'  
  num_row = len(matrix1)
  num_col = len(matrix1[0])
  assert num_row == num_row_origin/spatial_compress_rate
  assert num_col == num_col_origin/temporal_compress_rate
  
  dist = 0
  for i in range(num_row):
    for j in range(num_col):
      assert 0 <= matrix1[i][j] < 1
      assert 0 <= matrix2[i][j] < 1
      dist += abs(matrix1[i][j] - matrix2[i][j])

  return dist



# ########################################### end of granularity analysis code 




if __name__ == "__main__":
    sys.exit(main())
  
   
  




 
  
 









