def score_match_mismatch(row_nucleotide, col_nuclotide, match_score, mismatch_score):
  """Function to calculate score between two single nucleotides

  Parameters
  ----------
  row_nucleotide: character
                  One single nucelotide from the first sequence
  col_nuclotide: character
                  One single nucelotide from the second sequence
  match_score: int
               Score for when row and column bases match
  mismatch_score: int
                  Score for when row and column bases mismatch

  Returns
  -------
  score: int
         Match/mismatch/indel score
  """
  if row_nucleotide == col_nuclotide:
    return match_score
  else:
    return mismatch_score

def initialize_score_matrix_nw(row_seq_length, col_seq_length, gap_score):
  """Function to initialize scoring matrix

  Parameters
  -----------
  row_seq_length: int
                 Length of the sequence that you stack along the rows
  col_seq_length: int
                 Length of the sequence that you stack along the columns
  gap_score: int
             Score (penalty) for a gap (indel)
  Returns
  --------
  s: list of list
     Scoring matrix
  """
  s = [[0 for j in range(1+col_seq_length)] for i in range(1+row_seq_length)]
  for i in range(1, row_seq_length+1):
    s[i][0] = s[i-1][0] + gap_score
  for j in range(1, col_seq_length+1):
    s[0][j] = s[0][j-1] + gap_score
  return s

def run_needleman_wunsch(row_seq, col_seq, score_function, gap_score, match_score, mismatch_score):
  """Function to perform needleman wunsch algorithm

  Parameters
  -----------
  row_seq: str
           Sequence that you stack along the rows
  col_seq: str
           Sequence that you stack along the columns
  score_function: function
                  Function to calculate match/mismatch/gapscores given two nucleotides
  gap_score: int
             Score (penalty) for a gap (indel)
  match_score: int
               Score for when row and column bases match
  mismatch_score: int
                  Score for when row and column bases mismatch

  Returns
  -------
  traceback_matrix: list of list
                    A matrix to traceback the highest score calculations

  s: list of list
     Score matrix
  """
  # intialize scoring matrix
  row_seq_length = len(row_seq)
  col_seq_length = len(col_seq)


  """" WRITE CODE HERE TO INITIALIZE THE SCORING MATRIX
  """

  s = initialize_score_matrix_nw(row_seq_length = len(row_seq),
                                 col_seq_length = len(col_seq),
                                 gap_score = gap_score)

  traceback_matrix = [['-' for j in range(col_seq_length+1)] for i in range(row_seq_length+1)]
  for i in range(1, row_seq_length+1):
    traceback_matrix[i][0] =  'go_up'
  for j in range(1, col_seq_length+1):
    traceback_matrix[0][j] = 'go_left'


  for i in range(1, row_seq_length+1):
    for j in range(1, col_seq_length+1):
      col_base = col_seq[j-1]
      row_base = row_seq[i-1]
      diagonal_score = s[i-1][j-1] + score_function(row_base, col_base, match_score, mismatch_score)
      deletion_in_col = s[i-1][j] + gap_score
      deletion_in_row = s[i][j-1] + gap_score
      s[i][j] = max(diagonal_score, deletion_in_row, deletion_in_col)
      if s[i][j] == diagonal_score:
        traceback_matrix[i][j] = 'go_diagonal'
      elif s[i][j] == deletion_in_col:
        traceback_matrix[i][j] = 'go_up'
      else:
        traceback_matrix[i][j] =  'go_left'
  return traceback_matrix, s

def get_aligned_seqs(traceback_matrix, row_seq, col_seq):
  """Get the aligned sequences

  Parameters
  -----------
  traceback_matrix: list of list
                    A matrix to traceback the highest score calculations
  row_seq: str
           Sequence that you stack along the rows
  col_seq: str
           Sequence that you stack along the columns

  Returns
  -------
  row_seq_aligned: str
                   Aligned sequence along the rows
  col_seq_aligned: str
                   Aligned sequence along the columns
  """

  i = len(row_seq)
  j = len(col_seq)
  col_seq_aligned = ''
  row_seq_aligned = ''

  while (i>0 or j>0):
    direction = traceback_matrix[i][j]
    if direction == 'go_diagonal':
      row_seq_aligned = row_seq[i-1] + row_seq_aligned
      col_seq_aligned = col_seq[j-1] + col_seq_aligned
      i -= 1
      j -= 1
    elif i>0 and direction == "go_up":
      row_seq_aligned = row_seq[i-1] + row_seq_aligned
      col_seq_aligned = '-' + col_seq_aligned
      i -= 1
    elif j>0 and direction == "go_left":
      row_seq_aligned = '-' + row_seq_aligned
      col_seq_aligned = col_seq[j-1] + col_seq_aligned
      j -= 1
    elif direction == "-":
      pass
  return row_seq_aligned, col_seq_aligned

row_seq = "WHY"
col_seq = "WHAT"
gap_score = -2
match_score = 1
mismatch_score = -1
traceback_matrix, s = run_needleman_wunsch(row_seq = row_seq,
                                           col_seq = col_seq,
                                           score_function = score_match_mismatch,
                                           gap_score = gap_score,
                                           match_score = match_score,
                                           mismatch_score = mismatch_score)
rowseq_aligned, colseq_aligned = get_aligned_seqs(traceback_matrix=traceback_matrix, row_seq=row_seq, col_seq=col_seq)

print(colseq_aligned)
print(rowseq_aligned)

# also possible to get
#WHAT
#WHY-

row_seq = "GCATGCT"
col_seq = "GATTACA"

gap_score = -2
match_score = 1
mismatch_score = -1
traceback_matrix, s = run_needleman_wunsch(row_seq = row_seq,
                                           col_seq = col_seq,
                                           score_function = score_match_mismatch,
                                           gap_score = gap_score,
                                           match_score = match_score,
                                           mismatch_score = mismatch_score)
rowseq_aligned, colseq_aligned = get_aligned_seqs(traceback_matrix=traceback_matrix, row_seq=row_seq, col_seq=col_seq)

print(colseq_aligned)
print(rowseq_aligned)

row_seq = "GCATGCT"
col_seq = "GATTACA"

gap_score = -1
match_score = 1
mismatch_score = -1
traceback_matrix, s = run_needleman_wunsch(row_seq = row_seq,
                                           col_seq = col_seq,
                                           score_function = score_match_mismatch,
                                           gap_score = gap_score,
                                           match_score = match_score,
                                           mismatch_score = mismatch_score)
rowseq_aligned, colseq_aligned = get_aligned_seqs(traceback_matrix=traceback_matrix, row_seq=row_seq, col_seq=col_seq)

print(colseq_aligned)
print(rowseq_aligned)

row_seq = "GTTTGACCAGCC"
col_seq = "CTGACCCACCGC"

gap_score = -1
match_score = 2
mismatch_score = -1
traceback_matrix, s = run_needleman_wunsch(row_seq = row_seq,
                                           col_seq = col_seq,
                                           score_function = score_match_mismatch,
                                           gap_score = gap_score,
                                           match_score = match_score,
                                           mismatch_score = mismatch_score)
rowseq_aligned, colseq_aligned = get_aligned_seqs(traceback_matrix=traceback_matrix, row_seq=row_seq, col_seq=col_seq)

print(colseq_aligned)
print(rowseq_aligned)

col_seq = 'GCATGCT'
row_seq = 'GATACCA'

gap_score = -1
match_score = 2
mismatch_score = -2
traceback_matrix, s = run_needleman_wunsch(row_seq = row_seq,
                                           col_seq = col_seq,
                                           score_function = score_match_mismatch,
                                           gap_score = gap_score,
                                           match_score = match_score,
                                           mismatch_score = mismatch_score)
rowseq_aligned, colseq_aligned = get_aligned_seqs(traceback_matrix=traceback_matrix, row_seq=row_seq, col_seq=col_seq)
print(colseq_aligned)
print(rowseq_aligned)

def score_match_mismatch2(row_nucleotide, col_nuclotide, match_score, transition_score, transversion_score):
  """
  """
  chemical_class = {"A":"purine", "G":"purine", "T":"pyrimidine","C":"pyrimidine"}
  if row_nucleotide == col_nuclotide:
      return match_score

  if chemical_class[row_nucleotide] == chemical_class[col_nuclotide]:
      return transition_score
  else:
    return transversion_score

def run_needleman_wunsch2(row_seq, col_seq, score_function, gap_score, match_score,  transition_score, transversion_score):
  """Function to perform needleman wunsch algorithm

  Parameters
  -----------
  row_seq: str
           Sequence that you stack along the rows
  col_seq: str
           Sequence that you stack along the columns
  score_function: function
                  Function to calculate match/mismatch/gapscores given two nucleotides
  gap_score: int
             Score (penalty) for a gap (indel)
  match_score: int
               Score for when row and column bases match
  mismatch_score: int
                  Score for when row and column bases mismatch

  Returns
  -------
  traceback_matrix: list of list
                    A matrix to traceback the highest score calculations

  s: list of list
     Score matrix
  """
  # intialize scoring matrix
  row_seq_length = len(row_seq)
  col_seq_length = len(col_seq)


  """" WRITE CODE HERE TO INITIALIZE THE SCORING MATRIX
  """

  s = initialize_score_matrix_nw(row_seq_length = len(row_seq),
                                 col_seq_length = len(col_seq),
                                 gap_score = gap_score)

  traceback_matrix = [['-' for j in range(col_seq_length+1)] for i in range(row_seq_length+1)]
  for i in range(1, row_seq_length+1):
    traceback_matrix[i][0] =  'go_up'
  for j in range(1, col_seq_length+1):
    traceback_matrix[0][j] = 'go_left'


  for i in range(1, row_seq_length+1):
    for j in range(1, col_seq_length+1):
      col_base = col_seq[j-1]
      row_base = row_seq[i-1]
      diagonal_score = s[i-1][j-1] + score_match_mismatch2(row_base, col_base, match_score,  transition_score, transversion_score)
      deletion_in_col = s[i-1][j] + gap_score
      deletion_in_row = s[i][j-1] + gap_score
      s[i][j] = max(diagonal_score, deletion_in_row, deletion_in_col)
      if s[i][j] == diagonal_score:
        traceback_matrix[i][j] = 'go_diagonal'
      elif s[i][j] == deletion_in_col:
        traceback_matrix[i][j] = 'go_up'
      else:
        traceback_matrix[i][j] =  'go_left'
  return traceback_matrix, s

row_seq = "GCATGCT"
col_seq = "GATTACA"

gap_score = -1
match_score = 2
transition_score = -1
transversion_score = -2


traceback_matrix, s = run_needleman_wunsch2(row_seq = row_seq,
                                           col_seq = col_seq,
                                           score_function = score_match_mismatch2,
                                           gap_score = gap_score,
                                           match_score = match_score,
                                            transition_score = transition_score,
                                            transversion_score = transition_score)
rowseq_aligned, colseq_aligned = get_aligned_seqs(traceback_matrix=traceback_matrix, row_seq=row_seq, col_seq=col_seq)
print(rowseq_aligned)
print(colseq_aligned)

arrow_up = "\u2191"
arrow_right = "\u2192"
arrow_down = "\u2193"
arrow_left = "\u2190"
arrow_diag_down = "\u2198"
arrow_diag_up = "\u2196"

print(arrow_up)
print(arrow_right)

print(arrow_down)
print(arrow_left)
print(arrow_diag_down)
print(arrow_diag_up)

def pretty_print_score_matrix(score_matrix):
  """Pretty prints the scoring matrix

  Parameters
  ----------
  score_matrix: list of list
                The s matrix from any of the above cells

  Returns
  -------
  The matrix s printed in a more readable way (so it is reconisable as a matrix)
  """

  for row in score_matrix:
    print(" ".join(f"{x:4}" for x in row))


pretty_print_score_matrix(s)

def pretty_print_traceback_matrix(traceback_matrix):
    """Pretty prints the scoring matrix

    Parameters
    ----------
    score_matrix: list of list
                  The s matrix from any of the above cells

    Returns
    -------
    The traceback matrix s printed in a more readable way (so it is reconisable as a matrix and one can trace back using the arrows)
    """


    for index, row in enumerate(traceback_matrix):
      str_to_print = ""
      if index == 0:
        str_to_print = str_to_print + "  "
      for x in row:
        if x == 'go_left':
          str_to_print = str_to_print + " " + arrow_left
        elif x == 'go_up':
            str_to_print = str_to_print + " " + arrow_up
        elif x == 'go_diagonal':
            str_to_print = str_to_print + " " + arrow_diag_up
      print(str_to_print)


pretty_print_traceback_matrix(traceback_matrix)

