#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 19:02:45 2022

Auxillary Functions:
 * Provides a rudimentary reader (boundary_search) to recover the plasma boundary from the used test-g-file (G-EQDSK format).
 * Compability with other g-files besides the provided test-case is NOT guaranteed!

@author: Georg Graßler
"""

import numpy as np

def boundary_search(filename,exp_symbol):

  # Open the g-file in read mode and check if valid first line
  fileID = open(filename, "r")
  desc = fileID.readline()
  if not desc:
      raise IOError("Cannot read from input file")

  # Check first line for matching format
  s = desc.split()
  if len(s) < 3:
      raise IOError("First line must contain at least 3 numbers")

  # Determine the structure of entries by looking at second line (first float data line)
  desc = fileID.readline()

  # Use the number of occurrences of exponent indentifier e.g. 'e' in 1.234567e-1 to determine number of entries in line
  entries_per_line = desc.count(exp_symbol)

  # Then get the number of characters per entry
  characters_per_line = len(desc) - 1
  characters_per_entry = characters_per_line/entries_per_line

  # Check if characters per entry is integer -> else the assumption of uniform entry lenght is not fullfilled
  if (abs(np.floor(characters_per_entry) - characters_per_entry) != 0):
    print(characters_per_entry)
    raise IOError("Number of characters per entry have to be integer! Can't read file!")
  characters_per_entry = int(characters_per_entry)

  # Continue down file until floating point data ends, header of boundary data has only integers (aka no exponent identifiers)
  while (True):
    desc = fileID.readline()
    length = len(desc) - 1
    if (desc.count(exp_symbol) == 0):
      break
    elif(length < 1):
      raise IOError('End of File. Could not find boundary data!')

  # Retrieve from boundary data header the number of points
  s = desc.split()
  n_boundary_points = int(s[0])

  # Determine the number of lines containing boundary data: to each point belong 2 entries (R,Z) and per line there are entries_per_line-many entries
  n_boundary_rows = int(np.ceil(n_boundary_points*2/entries_per_line))
  
  # Fill the boundary list with all the entries (l-loop) of all the relevant lines (k-loop)
  boundary_list = []
  for k in range(n_boundary_rows):
    desc = fileID.readline()
    length = len(desc) - 1
    n_entries = int(length/characters_per_entry) # Have to acount for the fact that the last line might only be partially filled
    for l in range(n_entries):
      index_start = l*characters_per_entry
      index_end = index_start + characters_per_entry
      boundary_list.append(float(desc[index_start:index_end]))

  # Assigning the points to R and Z coordinate
  rbbbs = np.array(boundary_list[0::2])
  zbbbs = np.array(boundary_list[1::2])

  # Close the g-file
  fileID.close()
  
  return rbbbs, zbbbs
