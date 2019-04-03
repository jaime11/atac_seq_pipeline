import copy
import pysam
import collections
import os
import CGAT.Experiment as E
from CGATPipelines.Pipeline import cluster_runnable
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import tempfile
import re
import pandas as pd
from pandas.core.frame import DataFrame
import numpy as np
import matplotlib as mpl
mpl.use('Agg') # So that matplotlib.pyplot doesn't try to find the X-server and give an error
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, "/home/mbp15ja/dev/AuxiliaryPrograms/StringOperations/")
import StringOperations




# NOTE: Overwrites the infile
# Based on shortenBedHeaderStatement
# Generates a bash statement to be run which substitutes the header of the infile
# with another header where each Sample's text in header_add_text_file is added to header's sample
# Note that since sample names are unique, text will only be added to the first coincidence.
# 
# Inputs: 
#    -infile: dataframe_file: Uncompressed file with the information
#         chrom   start   end     15548_12_segments       22965_12_segments       23977_12_segments
#         chr1    0       200     E3      E3      E3
#         chr1    200     400     E3      E3      E3
#         chr1    400     600     E3      E3      E3
#
#    -header_add_text_file: A tab separated file without header and two columns. 
#      The first column containing sample names as they appear in the header of the 
#      dataframe_file (all prefixes and suffixes common to all samples are removed
#      when generating the tables with the data to cluster, so the sample names here must take that into account).
#      The second column, the text to add to each label. Only labels with text to be
#      added to the sample name should appear here, although if a complete sample name match is not found,
#      the label won't be added. Example:
#          Sample1 text_label_to_be_added
#          Sample2 text_label_to_be_added
#
#
# Outputs: The statement required to produce the output (uncompressed)
def addHeaderTextToSamples(infile, header_add_text_file):
    
    header_fields = []
    
    # First get the header fields
    with IOTools.openFile(infile, "r") as reader:
        
        for line in reader:
             
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
             
            header_fields = line.split("\t")
            
            break
    
    reader.close()
    
    
    # Go through each line in the header_add_text_file and add any corresponding text
    with IOTools.openFile(header_add_text_file, "r") as reader:
        
        for line in reader:
             
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
             
            adding_text_fields = line.split("\t")
            
            if(len(adding_text_fields) != 2):
                
                raise Exception("The file contains lines with number of fields different to 2")
            
            # If it gets here we have 2 fields. Go through each header_field and see if we have
            # the id to add text
            for position, header_field in enumerate(header_fields):
                
                if(header_field == adding_text_fields[0]):
                    
                    header_fields[position] += "_"+adding_text_fields[1]
                    
                    # There should only by 1 id with the text
                    break
    
    reader.close()
    
    
    
    # Create the statement to be filled
    statement = ""
    
    
    # We are going to modify the headers    
        
    # First delete the first line of the file
    statement += '''sed -i '1d' '''+infile+'''; '''
    
    # Then add the modified string with the header
    statement += '''sed -i '1i'''+("\\t".join(header_fields))+'''' '''+infile+''';'''
        
    return statement 


# Based on PipelineChromHMM.coverageHistogram
#  
# Gets the input file and outputs a file containing the most variable rows given
# the parameters.
#
# Inputs:
#     -infile: File containing the input data with rows as features and columns as samples.
#         Format (includes header):
#         
#         chrom   start   end     sample1       sample2
#         chr1    0       200     0.52      1.42
#         chr1    200     400     1.44      0.0012
#         chr1    400     600     12.45     5.122
#     -outfile: Output file containing at most the specified num_max_rows ordered by variability.
#     -start_data_col: The first column containing data, included (0 based position).
#     -end_data_col: The last column containing data, included (0 based position).
#        # Note that all the columns in between must contain numeric and non NaN data.
#     -num_max_rows: The maximum number of rows for the output to contain.
#
# Outputs:
#     -Writes the filtered file in outfile.
@cluster_runnable
def getTopVarianceRows(infile, 
                       outfile, 
                       start_data_col, 
                       end_data_col, 
                       num_max_rows):
    
    total_positions = 0
    
    full_df = None
    
    # Read the infile
    full_df = pd.read_table(infile,
                            header=0,
                            delimiter='\t')
        
    # Obtain the number of rows in the dataframe
    total_positions = len(full_df.index)
        
    # Get the variance into a new column
    # Take into account start_data_col:(end_data_col+1) since the last column specified is excluded    
    full_df['variance'] = (full_df.iloc[:,start_data_col:(end_data_col+1)]).var(axis=1)
    
    # Sort the values by variance
    full_df.sort_values('variance', axis=0, ascending=False, inplace=True, kind='quicksort')
    
    # Pick the top num_max_rows
    output_df = full_df.head(num_max_rows)
    
    # Delete the variance column
    del output_df['variance']
    
    
    
    
    # Output it
    
    chunksize = 10000000
    
    # To create the outfile we need to know when to append headers (first chunk)
    first_chunk=True
    
    
    for g, chunk in output_df.groupby(np.arange(len(output_df)) // chunksize):
        
        # Append headers
        if(first_chunk):
            
            # Open in write mode and output headers
            with open(outfile, 'w') as writer:
                chunk.to_csv(writer,
                             sep="\t",
                             header=True,
                             index_label=None,
                             index=False,
                             line_terminator='\n')
                
            
            writer.close()
            
            first_chunk = False
        
        # Don't append headers
        else:
            
            # Open in append mode and output headers
            with open(outfile, 'a') as writer:
                chunk.to_csv(writer,
                             sep="\t",
                             header=False,
                             index_label=None,
                             index=False,
                             line_terminator='\n')
                
            
            writer.close()
    
    
    
    
    



# Based on PipelineChromHMM.getRandomSampleRows
#  
# Gets the input file and outputs a file containing a random sample of rows given
# the parameters.
#
# Inputs:
#     -infile: File containing the input data with rows as features and columns as samples.
#         Format (includes header):
#         
#         chrom   start   end     sample1       sample2
#         chr1    0       200     0.52      1.42
#         chr1    200     400     1.44      0.0012
#         chr1    400     600     12.45     5.122
#     -outfile: Output file containing at most the specified num_max_rows ordered by variability.
#     -start_data_col: The first column containing data, included (0 based position).
#     -end_data_col: The last column containing data, included (0 based position).
#        # Note that all the columns in between must contain numeric and non NaN data.
#     -num_max_rows: The maximum number of rows for the output to contain.
#
# Outputs:
#     -Writes the filtered file in outfile.
@cluster_runnable
def getRandomSampleRows(infile, 
                       outfile, 
                       start_data_col, 
                       end_data_col, 
                       num_max_rows):
    
    total_positions = 0
    
    full_df = None
    
    # Read the infile
    full_df = pd.read_table(infile,
                            header=0,
                            delimiter='\t')
    
    # Pick a random sample of num_max_rows
    output_df = full_df.sample(n=num_max_rows)
    
    # Output it
    chunksize = 10000000
    
    # To create the outfile we need to know when to append headers (first chunk)
    first_chunk=True
    
    
    for g, chunk in output_df.groupby(np.arange(len(output_df)) // chunksize):
        
        # Append headers
        if(first_chunk):
            
            # Open in write mode and output headers
            with open(outfile, 'w') as writer:
                chunk.to_csv(writer,
                             sep="\t",
                             header=True,
                             index_label=None,
                             index=False,
                             line_terminator='\n')
                
            
            writer.close()
            
            first_chunk = False
        
        # Don't append headers
        else:
            
            # Open in append mode and output headers
            with open(outfile, 'a') as writer:
                chunk.to_csv(writer,
                             sep="\t",
                             header=False,
                             index_label=None,
                             index=False,
                             line_terminator='\n')
                
            
            writer.close()    
    
    
    
    
    
    

    
    
    


# Based on PipelineChromHMM.coverageHistogram
# Requirements: The values are equal or higher to 0
#  
# Computes a histogram taking the values from values_file and
# using the specified number_bins. Stores the output on the outfile.
#
# Inputs:
#     -values_file: File containing one value per line.
#     -number_bins: Number of bins to produce for the histogram.
#     -outfile: Outfile containing the histogram
#
# Outputs:
#     -A histogram with the characteristics specified
@cluster_runnable
def histogramFromOneDValues(values_file, number_bins, outfile):
    
    # We are going to do two passes:
    # 1) Get the maximum score of all the series (for the bins).
    # 2) Create the numpy array of instances to populate the histogram
    
    
    # We are going to read the values_file in chunks at a time, this indicates the size of each chunk
    chunk_size = 10000000

    
    # Pass 1
    # Go through the values_file, get:
    # -the maximum score
    maximum_score = float(0)
    
    
    # Read a chunk of lines in the values_file
    for chunk in pd.read_table(values_file,
                               names=['ends_distance'], 
                               header=None, 
                               chunksize=chunk_size, 
                               delimiter='\t'):
        
        # Update maximum distance if it has been surpassed
        max_score_chunk = chunk['ends_distance'].max()
        
        if(max_score_chunk > maximum_score):
            maximum_score = max_score_chunk
            
        
    # Now we are ready for pass 2
     
    # Get the bin edges for the histogram based on the max
    # min is always 0.0
    bin_edges = np.linspace(0.0, maximum_score, number_bins + 1)
    
    # These are going to be the actual bins, where in each
    # position we put the count
    # np.ndarray filled with np.uint32 zeros, CHANGED TO int64
    total_hist_data = np.zeros(number_bins, np.int64)
    
    
    # Read a chunk of lines in the values_file
    for chunk in pd.read_table(values_file,
                               names=['ends_distance'], 
                               header=None, 
                               chunksize=chunk_size, 
                               delimiter='\t'):
        
        
        # Iterates rows
        for row in chunk.itertuples():
            
            # Each row will contain score and the counts for that score
            # Bin the data for one score
            subtotal, edges = np.histogram(row[1], bins=bin_edges)
            
            # accumulate bin counts over chunks
            total_hist_data += subtotal
    
    
    # Turn interactive plotting off, we don't want to show it
    plt.ioff()
    
    # Plot the histogram
    plt.hist(bin_edges[:-1], bins=bin_edges, weights=total_hist_data)
    
    # Save the file, remove whitespace around the image    
    plt.savefig(outfile)
    
    # Clear the plot for the next one
    plt.clf() 



# Generates a report with all the grouped data
# 
# Inputs: 
#    -dataframe_file: Dataframe with the information:
#         chrom   start   end     15548_12_segments       22965_12_segments       23977_12_segments
#         chr1    0       200     E3      E3      E3
#         chr1    200     400     E3      E3      E3
#         chr1    400     600     E3      E3      E3
#    -output: Output file
#
# Outputs: Writting of the report to the output file
def generateGroupStatsReport(flagstats_stats,
                             library_complexity_stats,
                             mark_duplicates_stats,
                             number_peaks_stats,
                             number_reads_stats,
                             assigned_fraction_stats,
                             outfile):
    
    output_string = ""
    
    output_string += generateMarkDuplicatesSectionReport(mark_duplicates_stats)
    
    output_string += generateLibraryComplexitySectionReport(library_complexity_stats)
    
    output_string += generateFlagstatsSectionReport(flagstats_stats)
    
    output_string += "NUMBER OF SINGLE ENDS INPUTTED TO PEAK CALLER\t\n"
    
    output_string += ("Final number of single ends inputted to the peak caller\t"
                      +str(int(number_reads_stats))+
                      "\n")
    
    
    output_string += generateNumberPeaksSectionReport(number_peaks_stats)
    
    output_string += "ASSIGNED FRACTION (FILTERED READS IN PEAKS/TOTAL FILTERED READS)\t\n"
    
    output_string += ("Assigned fraction %\t"
                      +'%d' % int(round(float(assigned_fraction_stats)))+
                      "\n")
    
    
    # Get after_chr_sel_stats
    with IOTools.openFile(outfile, "w") as writer:
        
        writer.write(output_string)
        
    writer.close()
    
    



def generateFlagstatsSectionReport(flagstats_stats):
    
    output_string = ""
    
    # Output the section
    output_string += "FLAGSTATS\t\n"
    
    output_string += ("Total read pairs in proper pairs after deduplication\t"
                      +str(int(int(flagstats_stats['Total_reads_proper_pairs'])/2))+
                      "\n")
    
    return output_string



def generateLibraryComplexitySectionReport(library_complexity_stats):
    
    output_string = ""
    
    # Output the section
    output_string += "LIBRARY COMPLEXITY\t\n"
    
    
    output_string += ("Total read pairs\t"
                      +str(int(library_complexity_stats['TotalReadPairs']))+
                      "\n")
    
    output_string += ("Distinct read pairs\t"
                      +str(int(library_complexity_stats['DistinctReadPairs']))+
                      "\n")
    
    output_string += ("One read pair\t"
                      +str(int(library_complexity_stats['OneReadPair']))+
                      "\n")
    
    output_string += ("Two read pairs\t"
                      +str(int(library_complexity_stats['TwoReadPairs']))+
                      "\n")
    
    output_string += ("NRF=Distinct/Total\t"
                      '%.2f' % library_complexity_stats['NRF=Distinct/Total']+
                      "\n")
    
    output_string += ("PBC1=OnePair/Distinct\t"
                      '%.2f' % library_complexity_stats['PBC1=OnePair/Distinct']+
                      "\n")
    
    output_string += ("PBC2=OnePair/TwoPair\t"
                      '%.2f' % library_complexity_stats['PBC2=OnePair/TwoPair']+
                      "\n")
    
    return output_string
    
    


def generateMarkDuplicatesSectionReport(mark_duplicates_stats):
    
    output_string = ""
    
    # Output the section
    output_string += "MARK DUPLICATES\t\n"
    
    output_string += ("Reads pairs examined\t"
                      +str(int(mark_duplicates_stats['READ_PAIRS_EXAMINED']))+
                      "\n")
    
    output_string += ("Unpaired reads examined\t"
                      +str(int(mark_duplicates_stats['UNPAIRED_READS_EXAMINED']))+
                      "\n")
    
    output_string += ("Unpaired read duplicates\t"
                      +str(int(mark_duplicates_stats['UNPAIRED_READ_DUPLICATES']))+
                      "\n")
    
    output_string += ("Unmapped reads\t"
                      +str(int(mark_duplicates_stats['UNMAPPED_READS']))+
                      "\n")
    
    output_string += ("Read pairs marked as duplicates\t"
                      +str(int(mark_duplicates_stats['READ_PAIR_DUPLICATES']))+
                      "\n")
    
    output_string += ("Read pairs duplicates that were caused by optical (machine) duplication\t"
                      +str(int(mark_duplicates_stats['READ_PAIR_OPTICAL_DUPLICATES']))+
                      "\n")
    
    output_string += ("Fraction duplication\t"
                      +'%.2f' % mark_duplicates_stats['FRACTION_DUPLICATION']+
                      "\n")
    
    output_string += ("Estimated library size\t"
                      +str(int(mark_duplicates_stats['ESTIMATED_LIBRARY_SIZE']))+
                      "\n")
    
    return output_string




def generateNumberPeaksSectionReport(number_peaks_stats):
    
    
    output_string = ""
    
    # Output the section
    output_string += "NUMBER OF PEAKS GENERATED\t\n"
    
    output_string += ("Narrow pre-processing\t"
                      +str(int(number_peaks_stats['narrow pre-processing']))+
                      "\n")
    
    output_string += ("Broad pre-processing\t"
                      +str(int(number_peaks_stats['broad pre-processing']))+
                      "\n")
    
    output_string += ("Narrow post-processing\t"
                      +str(int(number_peaks_stats['narrow post-processing']))+
                      "\n")
    
    output_string += ("Broad post-processing\t"
                      +str(int(number_peaks_stats['broad post-processing']))+
                      "\n")
    
    return output_string
    


   


# Generates a bash statement to be run which filters out rows where
# the columns 4-end (all but the chr coordinates) are equal
# 
# Inputs: 
#    -dataframe_file: Compressed file with the information:
#         chrom   start   end     15548_12_segments       22965_12_segments       23977_12_segments
#         chr1    0       200     1      0.3      2
#         chr1    200     400     0.1      0.1      0.1
#         chr1    400     600     0.1      0.2      0.1
#    -output: Output file
#
# Outputs: The statement required to produce the output (uncompressed).
def deleteEqualBedRowsStatement(infile, output):
    
    number_of_fields = 0
    
    # First determine the number of fields
    with IOTools.openFile(infile, "r") as reader:
        
        for line in reader:
             
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
             
            fields = line.split("\t")
            
            number_of_fields = len(fields)
            
            break
    
    reader.close()
    
    # At least there has to be 4 fields: chr, start, end and label
    if number_of_fields < 4:
        
        raise Exception("The file: "+infile+" contains too few fields.")
    
    # If it gets here, the number of fields are correct
    # Create the statement to be filled
    statement = ""
    
    # If there are 4 fields, then we only have one cell line and there's nothing to do
    # except to copy the file
    if number_of_fields == 4:
        
        statement += '''zcat '''+infile+''' > '''+output+'''; '''
    
    else:
        
        # First we need to add the header of the file (which will not be involved in the processing with AWK):
        statement += '''zcat '''+infile+''' | head -n 1 > '''+output+'''; '''
        
        
        # Next we start outputting the file for processing from the second line
        # We want to create a statement of the type:
        # tail -n +2 infile | awk '!(($4==$5) && ($5==$6) && ...) {printf ("%s\n", $0)'} >> test
        statement += '''zcat '''+infile+''' | tail -n +2 | awk -F $'\\t' '!('''
        
        
        # Now we have to add the awk statement, we need to compare all the fields
        # after chr, start and end, make sure they are all the same and then inverse
        # the condition to choose only the different ones.
        # AWK uses 1 based numbering, we need the range to go from 4 to the
        # number of fields-1 
        # (we are always outputting the field and the next field, so 
        # we need to go to end-1). Eg. 10 fields, go to 9, last output 9 and 10
        for field in range(4, (number_of_fields)):
            
            # In the first adding we don't add &&
            if field == 4:
                statement += '''($'''+str(field)+'''==$'''+str((field+1))+''')'''
                
            else:
                statement += ''' && ($'''+str(field)+'''==$'''+str((field+1))+''')'''
                
        
        # Finish the statement        
        statement += ''') {printf ("%%s\\n", $0)'} >> '''+output+''';'''
                
    return statement
        
        

# NOTE: Overwrites the infile
# Generates a bash statement to be run which gets the shortest 
# different cell line strings from the header. Eliminating 
# common prefixes and suffixes to all cell lines.
# 
# Inputs: 
#    -dataframe_file: dataframe_file: Uncompressed file with the information
#         chrom   start   end     15548_12_segments       22965_12_segments       23977_12_segments
#         chr1    0       200     E3      E3      E3
#         chr1    200     400     E3      E3      E3
#         chr1    400     600     E3      E3      E3
#
# Outputs: The statement required to produce the output (uncompressed)
def shortenBedHeaderStatement(infile):
    
    # Get the number of fields and the header and the text of the header
    number_of_fields = 0
    
    header_fields = []
    
    # First determine the number of fields
    with IOTools.openFile(infile, "r") as reader:
        
        for line in reader:
             
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
             
            header_fields = line.split("\t")
            
            number_of_fields = len(header_fields)
            
            break
    
    reader.close()
    
    
    # At least there has to be 4 fields: chr, start, end and label
    if number_of_fields < 4:
        
        raise Exception("The file: "+infile+" contains too few fields.")
    
    # If it gets here, the number of fields are correct
    # Create the statement to be filled
    statement = ""
    
    
    # We are going to shorten the headers only for the cell lines
    # to avoid errors when plotting the figures
    # due to the cell line names being too long
    headers_to_modify = header_fields[3:]
    
    # Shorten the cell line headers, returns a list of shortened strings
    new_headers = StringOperations.shortenStringsByCommonPrefixAndSuffix(headers_to_modify)
    
    # Assign the new headers back (list of strings)
    all_new_headers = (header_fields[0:3] + new_headers)
    
    
    # First delete the first line of the file
    statement += '''sed -i '1d' '''+infile+'''; '''
    
    # Then add the modified string with the header
    statement += '''sed -i '1i'''+("\\t".join(all_new_headers))+'''' '''+infile+''';'''
        
    return statement    



# Creates a statement which gets all the narrow peaks provided in 
# "infiles", puts all the peaks into one file, coordinate sorts the file
# and merges the features taking into account merging_range
# Features separated by merging_range or less will be merged into 1.
# Outputs into the specified outfile
#
# Inputs:
#     -infiles: Narrow or broad peak files:
#                Format narrow peaks:
#                chr21   7926000 7926657 Peak_2  11675   .       22.95633        1174.55322      1167.59265      409
#                chr4    55328053        55328396        Peak_4  7221    .       41.97443        728.88098       722.18463       162
#
#                Format broad peaks:
#                chr21   7926000 7926657 Peak_2  5018    .       11.94926        506.90253       501.84232
#                chr4    55328053        55328396        Peak_4  3598    .       24.55315        365.07687       359.89346
#
#     -outfile: The outfile to write the output to
#                Format:
#                chromosome    start    end    strand
#
#     -merging_range: The maximum range used to merge the features separated
#         by this range.
def merge_all_narrow_peaks_statement(infiles, outfile, merging_range):
    
    # Begin the statement to put all the contents of the files into one
    statement = '''zcat '''
    
    for infile in infiles:
       
       statement += infile+ " "
       
    statement += ''' | sort -k1,1 -k2,2n | ''' # Coordinate sort the file
    
    # Get the chr, start, end and strand
    statement += ''' awk 'FS="\\t" {printf ("%%s\\t%%s\\t%%s\\t%%s\\n", $1, $2, $3, $6)'} | ''' 
    
    statement += ''' bedtools merge -i stdin -d ''' +merging_range # Merge the features
    
    statement += ''' | sort -k1,1 -k2,2n | gzip -c > '''+outfile # Zip the file
    
    return statement 
    
        
        
    


# Creates a statement to first get all the samples which correspond to
# each type (control or treatment), including the pooled sample.
# Then overlaps the pooled sample peaks with the peaks of the
# first sample of the type, then overlaps this result with the second
# sample, and so on...
# Does this for both "control" and "treatment" samples.
#
# Inputs:
#     
#     -infiles: All the peaks files (including control, treatment, pooled and others).
#     -outfile_control: The outfile of the process to contain the intersection of the peaks
#         corresponding to the samples of control condition.
#     -output_treatment: The outfile of the process to contain the intersection of the peaks
#         corresponding to the samples of control condition.
#     -sample_table: The tab separated table specifying the following 
#        (doesn't contain a header)
#        Sample_name    condition(control/treatment)    five_prime_quality_trimming(number of base pairs trimmed from both reads' 5' end in a pair)
#        nd2345    control    2
#
#     -file_suffix: The file suffix to add to each file to get from sample name to 
#            the target file basename
#     -sample_files_directory: The directory for the sample files
#
# Exception:
#     -If the pooled control and treatment is missing, can't do anything
def createNaiveOverlapThresholdingPeaksStatement(sample_table,
                                                 infiles,
                                                 output_control,
                                                 output_treatment,
                                                 file_suffix,
                                                 sample_files_directory):
     
    # To determine if at least one control sample and treatment sample are found
    control_sample_found = False
    treatment_sample_found = False
     
     
    samples_array = []
     
    # First get all the samples into an array
    with open(sample_table, 'r') as reader:
     
        for line in reader:
             
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
             
            fields = line.split("\t")
             
            samples_array.append(fields[0])
             
    reader.close()
     
     
    # Now we are going to extract the samples into 4:
    # -pooled_control
    # -pooled_treatment
    # -control_samples
    # -treatment_samples
    pooled_control = ""
    pooled_treatment = ""
    control_samples = []
    treatment_samples = []
    
    
    # Get the pooled files
    for infile in infiles:
         
        # Extract the pooled_control file
        match_pooled_control = re.match(sample_files_directory+"/pooled_control"+file_suffix, infile)
        
        if match_pooled_control:
           
            pooled_control = infile
            
            # To avoid double matching (in case)
            continue
            
        # Extract the pooled_control file
        match_pooled_treatment = re.match(sample_files_directory+"/pooled_treatment"+file_suffix, infile)
        
        if match_pooled_treatment:
           
            pooled_treatment = infile 
           
        
    # If after trying to match the pooled control and pooled treatment
    # are missing, Except
    if (pooled_treatment == "" and
        pooled_control == ""):
        
        raise Exception("Pooled control and Pooled treatment is missing")    
         
    
    # Create the statements, empty
    statement_control = ""
    statement_treatment = ""
      
     
    # Go through the array of samples, get the sample type and create the statements
    for sample in samples_array:
         
        # Get the condition: "control", "treatment" or "" if any other
        condition = getSampleCondition(sample, sample_table)
         
         
        if condition == "control":
             
            # Set the flag to true
            control_sample_found = True
             
            # Add the sample_name, the first sample_name added
            # creates a different statement to the rest
            if statement_control == "":
                statement_control += "intersectBed -wo -a "+pooled_control+" "
                statement_control += "-b " 
                statement_control += os.path.join(sample_files_directory, (sample+file_suffix)) # sample file name
                statement_control += ''' | awk 'BEGIN{FS="\\t";OFS="\\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq '''
            
            # Second and subsequent samples    
            else:
                statement_control += " | intersectBed -wo -a stdin "
                statement_control += "-b " 
                statement_control += os.path.join(sample_files_directory, (sample+file_suffix)) # sample file name
                statement_control += ''' | awk 'BEGIN{FS="\\t";OFS="\\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq '''
                            
             
        elif condition == "treatment":
             
            # Set the flag to true
            treatment_sample_found = True
             
            # Add the sample_name, the first sample_name added
            # creates a different statement to the rest
            if statement_treatment == "":
                statement_treatment += "intersectBed -wo -a "+pooled_treatment+" "
                statement_treatment += "-b " 
                statement_treatment += os.path.join(sample_files_directory, (sample+file_suffix)) # sample file name
                statement_treatment += ''' | awk 'BEGIN{FS="\\t";OFS="\\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq '''
            
            # Second and subsequent samples    
            else:
                statement_treatment += " | intersectBed -wo -a stdin "
                statement_treatment += "-b " 
                statement_treatment += os.path.join(sample_files_directory, (sample+file_suffix)) # sample file name
                statement_treatment += ''' | awk 'BEGIN{FS="\\t";OFS="\\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq '''
    
     
    # Add the final command to each statement
    statement_control += ''' | gzip -c > '''+output_control
    statement_treatment += ''' | gzip -c > '''+output_treatment
     
     
    # If control or treatment don't have any sample 
    # or if pooled_control or pooled_treatment are not found 
    # assign the empty string
    if ((not control_sample_found) or
        (pooled_control == "")):
         
        statement_control = ""
         
    if ((not treatment_sample_found) or
        (pooled_treatment == "")):
         
        statement_treatment = ""
         
     
    output_statement = statement_control
     
    # If there is a control statement and a treatment statement then add ; and checkpoint
    if ((statement_control != "") and
        statement_treatment != ""): 
         
        output_statement += ("; checkpoint; ")
     
     
    output_statement += statement_treatment
     
    return output_statement 
    
    


# Separates all the samples into control and treatment condition.
# Creates a statement to concatenate the contents of the files from each
# condition into the corresponding output file.
# Rows doesn't contain "treatment" or "control" (case insensitive) for the condition
# are ignored.
#
# Inputs:
#     -sample_table: The tab separated table specifying the following 
#        (doesn't contain a header):
#    
#        Sample_name    condition(control/treatment)    five_prime_quality_trimming(number of base pairs trimmed from both reads' 5' end in a pair)
#        nd2345    control    2
#     
#     -output_pooled_control: Output file containing a concatenation of all the 
#            single end control files.
#     -output_pooled_treatment: Output file containing a concatenation of all the 
#            single end treatment files.
#     -file_suffix: The file suffix to add to each file to get from sample name to 
#            the target file basename
#     -sample_files_directory: The directory for the sample files
#    
#
# Outputs:
#     -returns the statement to concatenate the contents of the files.
def virtualSEConditionPooling(sample_table,
                              output_pooled_control, 
                              output_pooled_treatment,
                              file_suffix,
                              sample_files_directory):
    
    # Begin the construct of the statements
    statement_control = '''zcat '''
    statement_treatment = '''zcat '''
    
    # To determine if at least one control and treatment are found
    control_found = False
    treatment_found = False
    
    
    samples_array = []
    
    # First get all the samples into an array
    with open(sample_table, 'r') as reader:
    
        for line in reader:
            
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
            
            fields = line.split("\t")
            
            samples_array.append(fields[0])
            
    reader.close()
    
    
    # Go through the array of samples, get the sample type and create the statements
    for sample in samples_array:
        
        # Get the condition: "control", "treatment" or "" if any other
        condition = getSampleCondition(sample, sample_table)
        
        
        if condition == "control":
            
            # Set the flag to true
            control_found = True
            
            # Add the sample_name 
            statement_control += (os.path.join(sample_files_directory, (sample+file_suffix)) + " ")
                
            
        elif condition == "treatment":
            
            # Set the flag to true
            treatment_found = True
            
            # Add the sample_name 
            statement_treatment += (os.path.join(sample_files_directory, (sample+file_suffix)) + " ")
            
    
    # Add the final command to each statement
    statement_control += ''' | gzip -c > '''+output_pooled_control
    statement_treatment += ''' | gzip -c > '''+output_pooled_treatment
    
    
    # If control or treatment don't have any sample assign the empty string
    if not control_found:
        
        statement_control = ""
        
    if not treatment_found:
        
        statement_treatment = ""
        
    
    output_statement = statement_control
    
    # If there is a control statement and a treatment statement then add ; and checkpoint
    if control_found and treatment_found:
        
        output_statement += ("; checkpoint; ")
    
    
    output_statement += statement_treatment
    
    return output_statement   
    



# Based on virtualSEConditionPooling
# Given a sample table with format specified and a condition
# returns True if the table doesn't contain at least one sample
# under that condition, False otherwise.
#
# Inputs:
#     -sample_table: The tab separated table specifying the following 
#        (doesn't contain a header):
#    
#        Sample_name    condition(control/treatment)    five_prime_quality_trimming(number of base pairs trimmed from both reads' 5' end in a pair)
#        nd2345    control    2
#
#     -condition: "control" or "treatment"
#
# Outputs:
#     -True if the table doesn't contain any sample with that condition,
#     False, otherwise.
#
# Exception:
#     -If the condition supplied is not either 'control' or 'treatment'
def emptyConditionSamplesTable(sample_table,
                               condition):
    
    condition_lower_case = condition.lower()
    
    if (condition_lower_case == "control" or 
         condition_lower_case == "treatment"):
        
        empty_table = True
    
    
        with open(sample_table, 'r') as reader:
    
            for line in reader:
                
                # Remove the new line character from the line (avoids taking this into account downstream)
                line = line.rstrip('\n')
                
                fields = line.split("\t")
                
                # Go through each sample, if the condition is found
                # return False
                if(fields[1].lower() == condition_lower_case):
                    
                    empty_table = False
                    
                    break
            
        reader.close()
        
        # Return the result after going scanning the file
        return empty_table
    
    
    # If the condition is not valid raise exception    
    else:
        
        raise Exception("Condition must be either 'control' or 'treatment'")
     






# Looks into sample_table for the sample specified by sample_name
# (it performs a case insensitive comparison).
# If it finds it, returns the shift specified for that sample.
# If it finds an empty string, it returns "0".
# If it doesn't find it, or if it is not an integer or empty returns an exception
# The shift is defined as the previously 5' constant (performed on all reads within a read pair) 
# trimming.
# For example, if 2bp have been trimmed from the 5' end in qc, 
# five_prime_quality_trimming=2
# Negative numbers are also allowed if the reads have been extended.
# Inputs:
#     -sample_name: The name of the sample
#     -sample_table: The tab separated table specifying the following 
#        (doesn't contain a header):
#    
#        Sample_name    condition(control/treatment)    five_prime_quality_trimming(number of base pairs trimmed from both reads' 5' end in a pair)
#        nd2345    control    2
#    
#
# Outputs:
#     -returns the shift specified for that sample or an exception
def getSampleQCShift(sample_name, sample_table):
    
    shift = ""
    
    # To make sure the sample is found
    found = False
    
    with open(sample_table, 'r') as reader:
    
        for line in reader:
            
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
            
            fields = line.split("\t")
            
            # If the name of the sample is equal to the name specified in the file
            if fields[0].lower() == sample_name.lower():
                
                shift = fields[2]
                
                # If it is not empty test whether it is an integer
                if shift != "":
                    
                    try:
                        
                        int(shift)
                        
                    except ValueError:
                        
                        raise Exception("The sample specified: "+sample_name+ " contains a non empty non number string")
                
                # If it is empty return 0
                else:
                    shift = "0"
                
                found = True
                
                break
    
    reader.close()
    
    if not found:
        
        raise Exception("The sample specified: "+sample_name+ " has not been found")
    
    else:
        
        return shift
    
    
    
    
    

# Looks into sample_table for the sample specified by sample_name
# (it performs a case insensitive comparison).
# If it finds it, returns the condition ("control" or "treatment") specified
# If it doesn't find it, or if it is not either "control" or "treatment" returns an empty string
# Inputs:
#     -sample_name: The name of the sample
#     -sample_table: The tab separated table specifying the following 
#        (doesn't contain a header):
#    
#        Sample_name    condition(control/treatment)    five_prime_quality_trimming(number of base pairs trimmed from both reads' 5' end in a pair)
#        nd2345    control    2
#    
#
# Outputs:
#     -returns the condition "control" or "treatment" for that sample,
#        If it doesn't find it, or if it is not either "control" or "treatment" returns an empty string
#        If the sample is not found returns an exception
def getSampleCondition(sample_name, sample_table):
    
    condition = ""
    
    # To make sure the sample is found
    found = False
    
    with open(sample_table, 'r') as reader:
    
        for line in reader:
            
            # Remove the new line character from the line (avoids taking this into account downstream)
            line = line.rstrip('\n')
            
            fields = line.split("\t")
            
            # If the name of the sample is equal to the name specified in the file
            if fields[0].lower() == sample_name.lower():
                
                # get in lowercase
                condition = fields[1].lower()
                
                # If it is not empty and its not "control" or "treatment" return an empty string    
                if (condition.lower() != "control"
                      and condition.lower() != "treatment"):
                            
                    condition = ""
                
                
                found = True
            
                break
    
    reader.close()
    
    if not found:
        
        raise Exception("The sample specified: "+sample_name+ " has not been found")
    
    else:
        
        return condition





# Creates a statement getting all the peak_files and counts the number of lines from each, 
# reporting the result to outfile.
#
# Inputs:
#     -peak_files: The peak files containing one peak per line
#     -outfile: The report file to store the results.
#     -tmp_dir: Temp dir to be used in the creation of files.
#
# Outputs:
#     -returns the statement to be run
def reportPeaks(peak_files, outfile, tmp_dir):
    
    # Create a temporal directory name for the run
    peaks_temp_dir = tempfile.mkdtemp(dir=tmp_dir)
    
    outfile_tmp = os.path.join(peaks_temp_dir, os.path.basename(outfile))
    
    statement = ""
    
    # We want to return the files in the order
    # 1) Preprocessed
    # 2) Within preprocessed, narrow and then broad peaks.
    # We rearrange the infiles (peak_files) to reflect this arrangement
    
    ordered_peak_files = len(peak_files)*[str]
    
    for peak_file in peak_files:
        
        if "filtered_peaks" in peak_file:
            
            if "narrow" in peak_file:
            
                ordered_peak_files[2] = peak_file
                
            else:
                
                ordered_peak_files[3] = peak_file
        else:
            
            if "narrow" in peak_file:
            
                ordered_peak_files[0] = peak_file
                
            else:
                
                ordered_peak_files[1] = peak_file

    
    for peak_file in ordered_peak_files:
        
        peak_type = ""
        
        processing_stage = ""
        
        # Decide from the name if its narrow or broad peaks
        if "narrow" in peak_file:
            
            peak_type = "narrow"
            
        else:
            
            peak_type = "broad"
            
            
        if "filtered_peaks" in peak_file:
            
            processing_stage = "post-processing"
            
        else:
            
            processing_stage = "pre-processing"
        
        statement += '''echo -n "'''+peak_type+''' '''+processing_stage+'''" >> '''+outfile_tmp+''';
                    echo -e -n '\\t' >> '''+outfile_tmp+''';
                    zcat '''+peak_file+''' | wc -l >> '''+outfile_tmp+''';
                    checkpoint;
                    '''
        
    # After everything has been outputed to the file, compress it to the outfile        
    statement += ''' cat '''+outfile_tmp+ ''' | gzip > '''+outfile+''';
    checkpoint;
    
    rm '''+outfile_tmp
        
    return statement


# Creates a statement to exclude from the infile the chrs
# indicated in excluded_chrs and saves the processed file to outfile
#
# Inputs:
#    -infile: peaks infile (compressed or uncompressed) to process.
#    -excluded_chrs: list of chrs to exclude (separated by ,)
#    -outfile: processed file without the chrs regions.
#
# Outputs:
#    -Writes the processed file to outfile.
@cluster_runnable
def createExcludingChrFromBedStatement(infile, excluded_chrs, outfile):
    
    statement = ""
    
    # If no excluded chrs provided just copy the file
    if excluded_chrs == "":
        
        statement = "cp "+infile+" "+outfile
    
    
    else:        
        
        # See the extension of the file to determine the cat/zcat command
        _, file_extension = os.path.splitext(infile)
        
        if file_extension == ".gz":
        
            statement = "zcat "
            
        else:
            
            statement = "cat "
            
        
        # Create the statement
        statement += infile
        statement += " | egrep -v '("
        statement += excluded_chrs
        statement += ")' | gzip -c > "
        statement += outfile
        statement += ";"
            
    return statement





# Creates a statement to exclude from the peaks infile the bed regions
# indicated in excluded_beds and saves the processed file to outfile
#
# Inputs:
#    -infile: bed infile (compressed or uncompressed) to process.
#    -excluded_beds: list of chrs to exclude (separated by |)
#    -outfile: processed file without the chrs.
#
# Outputs:
#    -Writes the processed file to outfile.
@cluster_runnable
def createExcludingBedsFromBedStatement(infile, excluded_beds, outfile):
    
    statement = ""
    
    # Get the string into a list, separating the elements
    list_excluded_beds = excluded_beds.split(",")
    
    # If no excluded beds provided just copy the file
    if len(list_excluded_beds) == 0:
        
        statement = "cp "+infile+" "+outfile
    
    
    else:                    
        
        # Create the statement
        # Because sometimes inputing a compressed file to bedtools intersect
        # can cause problems (https://github.com/arq5x/bedtools2/issues/513)
        # The file is zcat first and inputted to bedtools intersect -a stdin
        if infile.endswith(".gz"):
            
            statement += "zcat "+infile
            
        else:
            
            statement += "cat "+infile
        
        statement += " | bedtools intersect -v -a stdin"
        statement += " -b "
        
        for excluded_bed in list_excluded_beds:
            
            statement += excluded_bed
            statement += " "
        
        
        statement += " | gzip -c > "
        statement += outfile
        statement += ";"
            
    return statement
    
    
    


# inspired pipelineCaptCPerl.filterReadMappings
# Assuming a Sam/Bam file resulting from read pairs mapping
# and containing ONLY UNIQUE mappings per pair in each read
# and SORTED BY READNAME from any mapper output. 
# Because a pair R1 and R2 can only have one alignment each, 
# at most there would be 2 alignments, both alignments must 
# reflect either a proper pair or not a proper pair, same with
# pcr_duplicates. Therefore this information can be extracted 
# after both reads have been seen  
# Obtains statistics from the reads and read pairs.
# For individual reads:
#    -Unmapped reads
#    -Reads which are supplementary alignment (subread parts mapping to
#        multiple places.
#    -Exception if any non primary alignments found
# For read pairs:
#    -Read is not paired (will count 1 read per supposed read pair: not found)
#    -Read pair is not mapped in proper pair
#    -PCR duplicate
#    -Read pairs which are all of the below: 
#        -Read pair mapped in proper pair
#        -Doesn't contain supplementary alignments
#        -Doesn't contain secondary alignments
#        -Not PCR duplicate
#        -Contains R1 and R2
#    -Read pairs which are all of the below: 
#        -Read pair mapped in proper pair
#        -Contains supplementary alignments (maybe translocations, inversions, etc..)
#        -Doesn't contain secondary alignments
#        -Not PCR duplicate
#        -Contains R1 and R2
#
# Inputs:
#     -infile: Sam/Bam file containing only unique mappings per pair in each read
#         sorted by readname sorted by readname from any mapper output
#     -stats_outfile: Outfile with the statistics
#
# Output:
#     -returns printable output with the stats
#
# Exception:
#     -If a multimapping (read alignment not being primary is found).
@cluster_runnable
def getUniquelyMappedPairsNoMultimapping(infile, stats_outfile):
    
    c = collections.Counter()
    
    # Determine the extension of the file (BAM or SAM), done this for SAM/BAM compatibility
    # At the moment, all calls should only pass SAM, if it stops working just leave 
    # samfile = pysam.AlignmentFile(infile, "r")
    filename, file_extension = os.path.splitext(infile)
    
    if (file_extension.lower() == '.sam'):
        samfile = pysam.AlignmentFile(infile, "r")
    
    elif (file_extension.lower() == '.bam'):
        samfile = pysam.AlignmentFile(infile, "rb")
    
    
    # We need to know when the read name changes
    last_read_name = "first"
    
    # We need to know the last read when the read name changes to assign
    # the read_pair details from the last read when we have seen all the read pair
    # alignments
    last_read = ""
    
    # Create new dictionary template to store the details
    read_pair_dict_template = {'read_pair_is_paired':"", # Read pair contains R1 and R2    
            'read_pair_mapped_in_proper_pair':"", # Read pair contains proper mapped pair
            'read_pair_is_pcr_duplicate':"", # Read pair is a PCR duplicate
            'read_number_unmapped':0, # Number of R1 and R2 unmapped
            'read_number_supplementary':0, # Number of R1 and R2 supplementary (Eg. Translocations) 
            }
    
    # Create new copy of the template dictionary to store
    read_pair_dict = copy.deepcopy(read_pair_dict_template)
            
    
    # Process the reads
    for read in samfile:
        
        # If it is the first read, assign the name (it will consider
        # it as read name hasn't changed)
        if last_read_name == "first":
            last_read_name = read.query_name
            
            
        
        # The read has changed excepting the first read from the file
        if (last_read_name != read.query_name):
            
            # We have seen all the reads from a pair, fill the details
            # of the pair
            
            # reads contain a pair (R1 and R2)
            read_pair_dict['read_pair_is_paired'] = last_read.is_paired
            
            # Read pair contains proper mapped pair
            read_pair_dict['read_pair_mapped_in_proper_pair'] = last_read.is_proper_pair
            
            # Read pair is a PCR duplicate
            read_pair_dict['read_pair_is_pcr_duplicate'] = last_read.is_duplicate
            
            
            # Store everything updating the counters
            c = getUniquelyMappedPairsUpdateCounters(c, read_pair_dict)
            
            
            # Create new copy of the template dictionary to store
            read_pair_dict = copy.deepcopy(read_pair_dict_template)
            
        
        
        # Fill the dictionary with the details for the individual reads
        if read.is_unmapped:
            
            read_pair_dict['read_number_unmapped'] += 1
            
        if read.is_supplementary:
            
            read_pair_dict['read_number_supplementary'] += 1
            
        
        # If the read alignment is secondary, it means multimapping is enabled
        # exit with exception
        if read.is_secondary:
        
            raise Exception("Secondary read found: "+read.query_name+ " Multimapping was enabled")
        
        
        # Last command in the loop: update the last_read with
        # the read currently viewed
        last_read = read
        
        last_read_name = read.query_name
        
        
    # Process the last read
    # We have seen all the reads from a pair, fill the details
    # of the pair
    
    # reads contain a pair (R1 and R2)
    read_pair_dict['read_pair_is_paired'] = last_read.is_paired
    
    # Read pair contains proper mapped pair
    read_pair_dict['read_pair_mapped_in_proper_pair'] = last_read.is_proper_pair
    
    # Read pair is a PCR duplicate
    read_pair_dict['read_pair_is_pcr_duplicate'] = last_read.is_duplicate
    
    
    # Store everything updating the counters
    c = getUniquelyMappedPairsUpdateCounters(c, read_pair_dict)
    
    createOutputStringFromCounters(c, stats_outfile)
    




# Updates the counters of reads taking into account the guidelines from
# getUniquelyMappedPairs  
# Inputs:
#    -counter: Counter with the different keys and the counts
#    -dictionary: Dictionary containing the counts for the different keys
#        of the last read pair
#
# Outputs:
#    -updated_counter: Counter after adding the different keys and corresponding counts
def getUniquelyMappedPairsUpdateCounters(counter, dictionary):
    
    # Get proper mapped pairs
    if dictionary["read_pair_is_paired"] and dictionary["read_pair_mapped_in_proper_pair"]:
        
        counter["total_proper_mapped_pairs"] += 1
        
        # For supplementary, for now, all we want to know is whether a read
        # has these alignments or not
        if dictionary['read_number_supplementary'] != 0:
            
            counter["proper_mapped_pairs_with_supp_alignment"] += 1
        
        
        if dictionary["read_pair_is_pcr_duplicate"]:
            
            counter["proper_mapped_pairs_pcr_duplicate"] += 1
            
        # -Read pairs which are all of the below (group1): 
        #        -Read pair mapped in proper pair
        #        -Doesn't contain supplementary alignments
        #        -Doesn't contain secondary alignments
        #        -Not PCR duplicate
        #        -Contains R1 and R2
        if ((dictionary['read_number_supplementary'] == 0) and \
             not dictionary["read_pair_is_pcr_duplicate"]):
            
            counter["proper_mapped_pairs_with_no_supp_or_pcr_duplicate"] += 1
        
        #    -Read pairs which are all of the below (group2): 
        #        -Read pair mapped in proper pair
        #        -Contains supplementary alignments (maybe translocations, inversions, etc..)
        #        -Doesn't contain secondary alignments
        #        -Not PCR duplicate
        #        -Contains R1 and R2 
        elif ((dictionary['read_number_supplementary'] != 0) and \
             not dictionary["read_pair_is_pcr_duplicate"]):
    
            counter["proper_mapped_pairs_with_supp_or_pcr_duplicate"] += 1
        
    return counter    
        


# inspired pipelineCaptCPerl.filterReadMappings
# Assuming a Sam/Bam file resulting from read pairs mapping with BWA mem -M
# and containing unique or multiple mappings (1 pair - multiple alignments)
# per pair in each read and SORTED BY READNAME from any mapper
# output. 
# Obtains the number of correctly mapped read pairs and the number of those pairs
# correctly mapped read pairs which are marked as PCR duplicates
# Note: For each read pair, it only needs one proper read pair mapping in it's alignments
#     to be considered a correctly mapped read pair.
#     As long as one correctly mapped read pair is classified as a PCR duplicate,
#     the read pair will be considered a PCR duplicate.
# 
#
# Inputs:
#     -infile: Sam/Bam file containing only unique mappings per pair in each read
#         sorted by readname sorted by readname from BWA mem -M
#     -stats_outfile: Outfile with the statistics
#
# Output:
#     -returns printable output with the stats
@cluster_runnable
def getCorrectReadPairs(infile, stats_outfile):
    
    c = collections.Counter()
    
    # Determine the extension of the file (BAM or SAM), done this for SAM/BAM compatibility
    # At the moment, all calls should only pass SAM, if it stops working just leave 
    # samfile = pysam.AlignmentFile(infile, "r")
    filename, file_extension = os.path.splitext(infile)
    
    if (file_extension.lower() == '.sam'):
        samfile = pysam.AlignmentFile(infile, "r")
    
    elif (file_extension.lower() == '.bam'):
        samfile = pysam.AlignmentFile(infile, "rb")
    
    
    # We need to know when the read name changes
    last_read_name = "first"
    
    # We need to know the last read when the read name changes to assign
    # the read_pair details from the last read when we have seen all the read pair
    # alignments
    last_read = ""
    
    # Create new dictionary template to store the details
    read_pair_dict_template = {'read_pair_is_paired':False, # Read pair contains R1 and R2    
            'read_pair_mapped_in_proper_pair':False, # Read pair contains proper mapped pair
            'read_pair_is_pcr_duplicate':False, # Read pair is a PCR duplicate
            'proper_pair_additional_hit_equally_good':False, # If additional hits with equally good score are found
            'proper_pair_with_chimeric_read':False, # Chimeric read
            }
    
    # Create new copy of the template dictionary to store
    read_pair_dict = copy.deepcopy(read_pair_dict_template)
            
    
    # Process the reads
    for read in samfile:        
        
        # If it is the first read, assign the name (it will consider
        # it as read name hasn't changed
        if last_read_name == "first":
            last_read_name = read.query_name
            
        
        # The read has changed excepting the first read from the file
        if (last_read_name != read.query_name):
            
            # We have seen all the reads from a pair, update the counters
            
            # Store everything updating the counters
            c = getCorrectReadPairsUpdateCounters(c, read_pair_dict)
            
            
            # Create new copy of the template dictionary to store
            read_pair_dict = copy.deepcopy(read_pair_dict_template)
        
        
        # Checks for the read
        # For every read pair we check whether it is a proper pair, 
        # We only need one mapping from the read pair as a proper pair to
        # consider it a proper pair towards the output
        if(read.is_paired and read.is_proper_pair):
            
            read_pair_dict['read_pair_is_paired'] = read.is_paired
            
            read_pair_dict['read_pair_mapped_in_proper_pair'] = read.is_proper_pair
            
            # We only need one proper pair marked as duplicate for a read pair
            # to consider it a duplicate. 
            if read.is_duplicate:
                read_pair_dict['read_pair_is_pcr_duplicate'] = read.is_duplicate
            
            # Check additional hits when the read pair happens
            if read.has_tag("XA"):
                
                # If there are additional hits
                if (len(read.get_tag("XA").split(","))>1):
                    
                    # If the score on the additional hit (XS) is as good or better
                    # than the score of the original hit (AS)
                    if (read.get_tag("AS") <= read.get_tag("XS")):
                    
                        read_pair_dict['proper_pair_additional_hit_equally_good'] = True
                        
            # Check chimeric hits (split reads)
            if read.has_tag("SA"):
                 
                # If there are chimeric hits
                if (len(read.get_tag("SA").split(","))>1):
                    
                    read_pair_dict['proper_pair_with_chimeric_read'] = True
        
        # Last command in the loop: update the last_read with
        # the read currently viewed
        last_read = read
        
        last_read_name = read.query_name
        
        
    # Process the last read
    # We have seen all the reads from a pair, update the counters
            
    # Store everything updating the counters
    c = getCorrectReadPairsUpdateCounters(c, read_pair_dict)
    
    createOutputStringFromCounters(c, stats_outfile)




# Updates the counters of reads taking into account the guidelines from
# getCorrectReadPairs  
# Inputs:
#    -counter: Counter with the different keys and the counts
#    -dictionary: Dictionary containing the counts for the different keys
#        of the last read pair
#
# Outputs:
#    -updated_counter: Counter after adding the different keys and corresponding counts
def getCorrectReadPairsUpdateCounters(counter, dictionary):
    
    # Variable to control the pair isn't into any "incorrect category"
    correct_pair = True
    
    if (dictionary["read_pair_mapped_in_proper_pair"] and 
        dictionary["read_pair_is_paired"]):
        
        counter["read_pairs_mapped_in_proper_pair"] += 1
    
        if dictionary["read_pair_is_pcr_duplicate"]:
            
            counter["proper_read_pair_is_pcr_duplicate"] += 1
            
            correct_pair = False

    
        if dictionary['proper_pair_additional_hit_equally_good']:
            
            counter["proper_read_pair_additional_hit_equally_good"] += 1
            
            correct_pair = False
            
            
        if dictionary['proper_pair_with_chimeric_read']:
            
            counter["proper_read_pair_with_chimeric_read"] += 1
            
            correct_pair = False
        
        # After performing all the checks, see which are "correct" pairs    
        if correct_pair:
            
            counter["proper_read_pair_no_dup_no_additional_hit_equally_good_no_chimeric"] += 1
    
    # Not proper pairs    
    else:
        
        counter["read_pairs_mapped_not_in_proper_pair"] += 1
        
                    
    return counter 


# Writes the counters in the stats_outfile
# Inputs:
#     -counters: Counters with different keys and counts
#     -stats_outfile: File with the output of the stats
def createOutputStringFromCounters(counters, stats_outfile):
    
    with open(stats_outfile, 'w') as writer:
    
        for key in counters:
            
            writer.write(key + ":\t" +str(counters[key]) + "\n")
            
    writer.close()




###########################################################
# Heavily inspired on /pipeline_capt_c_perl_bwa_all_alignments/src/pipelineCaptCPerl.getUnmappedReads
# Only changed readname being regular read name
# Output adding the possibility to be compressed
###########################################################
# ASSUMING THE INFILE IS ORDERED BY READNAME!!!!
# If a read contains at least one mapping, it is considered a mapped read.
# If it doesn't, it is considered an unmapped read. 
#
#    Inputs:
#        -infile: SAM/BAM file ordered by readname
#        -outfile: Table with statistics containing total reads, mapped and unmapped reads.
#                 The output file must finish in ".tsv"
#
#    Outputs:
#        -Writing the table outfile with counters with mapped, unmapped and total reads.
#        -Writing output files containing unique lists of reads which are mapping and not mapping:
#            -.mapped.txt.gz: Mapped read names
#            -.unmapped.txt.gz: Unmapped read names
@cluster_runnable
def getMappedUnmappedReads(infile, outfile):
    
    # Determine the extension of the file (BAM or SAM), done this for SAM/BAM compatibility
    # At the moment, all calls should only pass SAM, if it stops working just leave 
    # samfile = pysam.AlignmentFile(infile, "r")
    filename, file_extension = os.path.splitext(infile)
    
    if (file_extension.lower() == '.sam'):
        samfile = pysam.AlignmentFile(infile, "r")
    
    elif (file_extension.lower() == '.bam'):
        samfile = pysam.AlignmentFile(infile, "rb")
    
    # Unmapped readnames
    unmapped = []
    
    # Mapped readnames
    mapped = []
    
    
    dictionary = {'read_name': 'first', 'mapped': 'no'}
    
    c = collections.Counter()
       
    # Go through the reads
    for read in samfile:
           
        # Get the read name
        read_name = read.query_name
           
           
        # First case, assign the read name to the dictionary
        if dictionary['read_name'] == 'first':
            dictionary['read_name'] = read_name
               
            # Add one read to the total
            c["total"] += 1
               
   
           
        # If a fragment from the read has been seen before
        if dictionary['read_name'] == read_name:  
               
            # If that fragment was mapped
            if dictionary['mapped'] == 'yes':
                   
                # Skip to the next data line
                continue
               
            # If all the prior fragments from the read haven't been mapped
            else:
               
                # If it is mapped, change the state to mapped (the default beginning state is unmapped)
                if not read.is_unmapped:
                   
                    # Put the read as mapped
                    dictionary['mapped'] = 'yes'
                       
                   
                
                
        # If a fragment from the same read hasn't been seen before (new read)        
        else:
                           
               
            # If the last read wasn't mapped, store it
            if dictionary['mapped'] == 'no':
                unmapped.append(dictionary['read_name'])
                   
                # Add one read to the unmapped
                c["unmapped"] += 1
                  
            else:
                  
                mapped.append(dictionary['read_name'])
                  
                # Add one read to the mapped
                c["mapped"] += 1
               
               
            # Add one read to the unmapped
            c["total"] += 1
               
            # Create a new dictionary with unmapped and the current read name
            dictionary = {'read_name': read_name, 'mapped': 'no'}
               
            # If it is mapped, change the state to mapped (the default beginning state is unmapped)
            if not read.is_unmapped:
               
                # Put the read as mapped
                dictionary['mapped'] = 'yes'
                   
       
    # The last read is hanging
    # If it is unmapped, store it
    if dictionary['mapped'] == 'no':
        unmapped.append(read_name)
           
        # Add one read to the unmapped
        c["unmapped"] += 1
          
    else:
                  
        mapped.append(dictionary['read_name'])
          
        # Add one read to the mapped
        c["mapped"] += 1
      
     
     
    # Create the name of the files for the unmapped and mapped reads
    out_unmapped_basename = P.snip(os.path.basename(outfile), '.tsv') + '.unmapped.txt.gz'
    out_unmapped = os.path.join(os.path.dirname(outfile), (out_unmapped_basename))
     
     
    out_mapped_basename = P.snip(os.path.basename(outfile), '.tsv') + '.mapped.txt.gz'
    out_mapped = os.path.join(os.path.dirname(outfile), (out_mapped_basename))
     
     
    # Store the unmapped read names in the outfile specified
    with IOTools.openFile(out_unmapped, "w") as output_file_write:
         
        for unmapped_read in unmapped:
            output_file_write.write(unmapped_read + '\n')
             
    output_file_write.close()
     
     
     
    # Store the mapped read names in the outfile specified
    with IOTools.openFile(out_mapped, "w") as output_file_write:
         
        for mapped_read in mapped:
            output_file_write.write(mapped_read + '\n')
             
    output_file_write.close()
             
     
     
    # Store the counters        
     
    # Create a name for the counters
    with IOTools.openFile(outfile, "w") as output_file_write: 
         
        for counter in c:
             
            output_file_write.write(counter + ':\t' + str(c[counter]) + '\n')
             
    output_file_write.close()


# Gets the statistics reports generated by the different states in the pipeline
# Puts them all formatted in a file
#
# Inputs:
#    -sample: The sample name used.
#    -initial_stats: Initial mapped, not mapped, total statistics.
#    -after_chr_sel_stats: After chr selection mapped, not mapped, total statistics.
#    -after_marking_dups_stats: After marking duplicates stats (various fields).
#    -outfile: The outfile containing the concatenation of tables.
#
# Outputs: writes the contents of the stats to a file.
def combineStatsReports(sample,
                        initial_stats,
                        after_chr_sel_stats,
                        after_marking_dups_stats,
                        outfile):
    
    # Create an empty stats dictionary to fill
    stats_dictionary = {}
    
    # Get initial stats
    with IOTools.openFile(initial_stats, "r") as reader:
         
        for line in reader:
            
            line = line.rstrip('\n')
            
            # Get the field and value
            field = line.split("\t")[0]
            
            value = line.split("\t")[1]
            
            # The field can contain : at the end, strip it
            if field[-1] == ":":
                
                field = field[0:-1]
            
            stats_dictionary[("initial_"+field)] = value
             
    reader.close()
    
    
    # Get after_chr_sel_stats
    with IOTools.openFile(after_chr_sel_stats, "r") as reader:
         
        for line in reader:
            
            line = line.rstrip('\n')
            
            # Get the field and value
            field = line.split("\t")[0]
            
            value = line.split("\t")[1]
            
            # The field can contain : at the end, strip it
            if field[-1] == ":":
                
                field = field[0:-1]
            
            stats_dictionary[("after_chr_sel_"+field)] = value
             
    reader.close()
    
    
    
    # Get after_marking_dups_stats
    with IOTools.openFile(after_marking_dups_stats, "r") as reader:
         
        for line in reader:
            
            line = line.rstrip('\n')
            
            # Get the field and value
            field = line.split("\t")[0]
            
            value = line.split("\t")[1]
            
            # The field can contain : at the end, strip it
            if field[-1] == ":":
                
                field = field[0:-1]
            
            stats_dictionary[("after_dedupping_"+field)] = value
             
    reader.close()
    
    # Get after_chr_sel_stats
    with IOTools.openFile(outfile, "w") as writer:
        
        # First output the sample
        writer.write("Sample\t")
        writer.write(sample+"\n")
        
        # Write the initial stats first
        if "initial_total" in stats_dictionary:
            writer.write("Initial_total\t")
            writer.write(stats_dictionary["initial_total"]+"\n")
        
        if "initial_unmapped" in stats_dictionary:    
            writer.write("Initial_unmapped\t")
            writer.write(stats_dictionary["initial_unmapped"]+"\n")
        
        if "initial_mapped" in stats_dictionary:
            writer.write("Initial_mapped\t")
            writer.write(stats_dictionary["initial_mapped"]+"\n")
        
        # Write the after_chr_sel
        if "after_chr_sel_total" in stats_dictionary:
            writer.write("After_chr_sel_total\t")
            writer.write(stats_dictionary["after_chr_sel_total"]+"\n")
            
        if "after_chr_sel_unmapped" in stats_dictionary:
            writer.write("After_chr_sel_unmapped\t")
            writer.write(stats_dictionary["after_chr_sel_unmapped"]+"\n")
        
        if "after_chr_sel_mapped" in stats_dictionary:
            writer.write("After_chr_sel_mapped\t")
            writer.write(stats_dictionary["after_chr_sel_mapped"]+"\n")
        
        # Write the deduplication (multiple categories)
        for key, value in stats_dictionary.iteritems():
            
            if key.startswith("after_dedupping_"):
                
                writer.write(key+"\t")
                writer.write(value+"\n")
        
         
    writer.close()
