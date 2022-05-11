import numpy as np

def boundary_search(filename,exp_symbol):

  fileID = open(filename, "r")
  desc = fileID.readline()
  if not desc:
      raise IOError("Cannot read from input file")

  s = desc.split()
  if len(s) < 3:
      raise IOError("First line must contain at least 3 numbers")

  desc = fileID.readline()
  entries_per_line = desc.count(exp_symbol)
  characters_per_line = len(desc) - 1
  characters_per_entry = characters_per_line/entries_per_line

  if (abs(np.floor(characters_per_entry) - characters_per_entry) != 0):
    print(characters_per_entry)
    raise IOError("Number of characters per entry have to be integer! Can't read file!")

  characters_per_entry = int(characters_per_entry)

  while (True):
    desc = fileID.readline()
    length = len(desc) - 1

    if (desc.count(exp_symbol) == 0):
      break
    elif(length < 1):
      raise IOError('End of File. Could not find boundary data!')

  s = desc.split()
  n_boundary_points = int(s[0])
  n_boundary_rows = int(np.ceil(n_boundary_points*2/5))
  boundary_list = []

  for k in range(n_boundary_rows):
    desc = fileID.readline()
    length = len(desc) - 1
    n_entries = int(length/characters_per_entry)
    for l in range(n_entries):
      index_start = l*characters_per_entry
      index_end = index_start + characters_per_entry
      boundary_list.append(float(desc[index_start:index_end]))

  rbbbs = np.array(boundary_list[0::2])
  zbbbs = np.array(boundary_list[1::2])

  fileID.close()
  
  return rbbbs, zbbbs
