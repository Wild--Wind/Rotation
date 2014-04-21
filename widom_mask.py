__author__ = 'greene'

import os, sys, getopt, operator

def all_blank(items):
    return all(x == [] for x in items)

def check_isdir(folder_path):
    try:
        assert os.path.isdir(folder_path)
        return True
    except AssertionError:
        print "%s folder path invalid" % folder_path
        return False

def check_empty(filename):
    if os.stat(filename).st_size > 0:
        empty= False
    elif os.path.isdir(filename):
        # this is a folder within the folder, skip it
        empty= True
    elif ''.join(line for line in open(filename).readlines() if not line.isspace()).isspace():
        empty = True
    else:
        print "%s\tFILE IS EMPTY" % filename
        empty= True
    return empty

def check_file_exists(filename):
    try:
        open(filename, 'r')
        return True
    except OSError:
        print "%s\tFILE DOES NOT EXIST" % filename
        return False

def valid_mfi(mfi):
    try:
        float(mfi)
        if float(mfi) > 0:
            return True
        else:
            print "\nERROR: min value (%s) must be positive, nonzero number\n" % mfi
            return False
    except ValueError:
        print "\nERROR: min value (%s) not valid number\n" % mfi
        return False

def valid_fraction(fraction):
    try:
        float(fraction)
        if float(fraction) > 0 and float(fraction) < 1:
            return True
        else:
            print "\nERROR: max_fraction value (%s) not valid fraction\n" % fraction
            return False
    except ValueError:
        print "\nERROR: max_fraction value (%s) not valid number\n" % fraction
        return False

def valid_cutoffs(mfi, upper, fraction):
    if valid_mfi(mfi) and valid_mfi(upper) and valid_fraction(fraction):
        if float(mfi) > float(upper):
            print "\nERROR: upper_fmi_min (%s) must be greater than mfi_min (%s)\n" % (upper, mfi)
            return False
        else:
            return True

def list_folder_files(folder_path):
    if check_isdir(folder_path):
        snp_files = []
        empty_files = []
        for snp_file in os.listdir(folder_path):
            snp_file= os.path.join(folder_path, snp_file)
            if check_empty(snp_file)== False and check_file_exists(snp_file) == True:
                snp_files.append(snp_file)
            elif check_empty(snp_file) and not check_isdir(snp_file):
                empty_files.append(snp_file)
            else:
                pass
        return snp_files, empty_files
    else:
        exit(1)

def replace_dbl_quotes(line, block_list):
    # for some annoying reason, a lot of the .csv files that were opened in excel
    # have lines of form '"item","item2","item3",...'
    # so if this is the case, remove the double quotes
    if '"' in line:
        block_list.append(line.replace('"','').strip().split(','))
    else:
        block_list.append(line.strip().split(','))
    return block_list

####################################################################################################################

def mfi_reader(mfi_file):
    # reads the luminex output files. requires that the line 'DataType:,Median' be in the file
    # Once DataType:,Median is found, collect the block of MFI value lines
    # After this block if there is a blank line or line of empty commas, stop collecting lines
    mfi_file= open(mfi_file)
    mfi_block = [] # collects the lines of relevant data (MFI's)
    program_info = [] # collects the lines preceding the MFI data block
    found = False

    for line in mfi_file:
        if found: # we're in the block of relevant data
            if line == '\n' or ",,,,,," in line: # many of the .csv's opened in excel have empty lines = ',,,,,'
                break
            elif 'Location' in line: # this is the first item in the header field
                program_info = replace_dbl_quotes(line, program_info)
            mfi_block = replace_dbl_quotes(line, mfi_block)
        else:
            if 'DataType' in line and 'Median' in line: # This is the line that starts the relevant info
                found = True
                program_info = replace_dbl_quotes(line, program_info)
            else:
                program_info = replace_dbl_quotes(line, program_info)
    if found == False: # the file doesn't contain info in the expected format
        print "\n\nERROR: File in unexpected format\nFile: %s\n" % mfi_file
        print "Expect to find 'DataType:,Median' \n"
        exit(1)

    mfi_file.close()
    return mfi_block, program_info


def map_sample_to_snps(mfi_block):
    # returns a list of dictionaries with keys = col headers vals = row value
    # also returns the headers
    headers = mfi_block[0]
    samples = []
    for sample in mfi_block[1:]:
        sample_dict = dict(zip(headers, sample))
        samples.append(sample_dict)

    return samples, headers


def map_snps_to_samples(samples, snp_names):
    # takes in a list of sample dictionaries (ouput from map_sample_to_snps)
    # returns dictionary with key = snp_allele_name (ex rs129770-C)
    # values = { (Location, Sample): MFI ... }
    # { rs12256-C : { (Location, Sample) : MFI... } ... }
    snps = {}
    for snp_name in snp_names:
        well_samps = []
        samp_mfis = []
        for sample in samples:
            well_samp = (sample['Location'], sample['Sample'])
            well_samps.append(well_samp)
            samp_mfis.append(sample[snp_name])
        snp_dict = dict(zip(well_samps, samp_mfis))
        snps[snp_name] = snp_dict

    return snps


def update_masked_samples(masked_snps, masked_samples):
    # masked_samples is a dict with keys = (Location, Sample), values = { snp:MFI, ...}
    # masked_snps is a dict with keys = snp-allele, values = { (Location, Sample):MFI ...}
    # return an updated version of masked_samples based on the masked_snps dictionary
    done = []
    for well_samp in masked_samples.keys():
        for snpname in masked_snps.keys():
            if (well_samp, snpname) not in done:
                done.append((well_samp, snpname))
                if masked_samples[well_samp][snpname] != masked_snps[snpname][well_samp]:
                    if masked_snps[snpname][well_samp] != '-1':
                        print 'ERROR: Misalignment in sample \ snp dictionaries!'
                        exit(1)
                    else:
                        masked_samples[well_samp][snpname] = '-1'
    return masked_samples #is now updated based on masked snps


def mask_sample(sample, allele_matches, rs_index, mfi_cutoff, upper_mfi_cutoff):
    cutoff = float(mfi_cutoff)
    upper_mfi = float(upper_mfi_cutoff)
    done = []
    below_upper = []
    masked_sample = sample
    well_samp = (sample['Location'], sample['Sample'])
    for key in masked_sample.keys():
        if '-' in key and key.split('-')[rs_index] not in done:
            match_set = allele_matches[key.split('-')[rs_index]]
            a1_mfi = float(masked_sample[match_set[0]])

            a2_mfi = float(masked_sample[match_set[1]])
            done.append(key.split('-')[rs_index])
            if a1_mfi < 0 or a2_mfi < 0:
                masked_sample[match_set[0]] = '-1'
                masked_sample[match_set[1]] = '-1'
            elif a1_mfi >= cutoff or a2_mfi >= cutoff:
                if a1_mfi <= upper_mfi and a1_mfi >= cutoff:
                    below_upper.append(match_set[0])
                if a2_mfi <= upper_mfi and a2_mfi >= cutoff:
                    below_upper.append(match_set[1])
            else:
                masked_sample[match_set[0]] = '-1'
                masked_sample[match_set[1]] = '-1'
    return masked_sample, well_samp, below_upper


def find_matching_alleles(headers):
    # headers is a list of header titles in a luminex file. includes rs numbers
    done= []
    matches = {}
    for key in headers:
        if '-' in key: # the snp names are expected to have '-base' in their names
            if key.count('-') == 1:
                rs_index = 0
            elif key.count('-') == 2:
                rs_index = 1 # Marco told me about some instances where the rs items form= num-rs_name-allele
            else: # otherwise I don't know what the format is, so fail
                print "Unknown SNP name format. %s" % key
                print "SNP Names should be in format <name>-<allele> OR"
                print "format <rand_num>-<name>-<allele>"
                exit(1)
            #scan the keys for its match
            for key2 in headers:
                if key != key2 and key not in done and key2 not in done and '-' in key2:
                    if key.split('-')[rs_index] == key2.split('-')[rs_index]:
                        matching_set = (key, key2)
                        done.append(key)
                        done.append(key2)
                        matches[key.split('-')[rs_index]] = matching_set
    return matches, rs_index


def should_mask(values, fraction_masked):
    # values will be a list of MFI values, with '-1' representing masked values.
    # if the fraction of masked values is greater than or equal to that indicated by the user,
    # all values should then be masked

    cutoff_number = int(len(values)*float(fraction_masked))
    count = 0
    for val in values:
        if val == '-1':
            count += 1
    if count >= cutoff_number:
        return True
    else:
        return False


def mask_values(dictionary, fraction_masked):
    # mask the all values for an entire snp or sample
    masked_dict = dictionary
    if should_mask(dictionary.values(), fraction_masked): #check if the fraction of -1 values exceeds given max_frac
       for key in dictionary.keys():
           if '-' in key or type(key) == tuple: #only the snp name fields should be masked
                masked_dict[key] = '-1'

    else:
        pass
    return dictionary


def update_below_upper(final_masked_dict, below_upper_dict):
    # this is a really ugly function
    # takes in the final masked dictionary and the below_upper_dict and if a sample at a snp = -1, remove
    # the snp from the sample list in below_upper_dict
    final_below_upper = below_upper_dict
    for sample in final_below_upper.keys():
        for snp_allele in final_below_upper[sample]:
            for masked_sample in final_masked_dict:
                if sample == masked_sample:
                    if final_masked_dict[sample][snp_allele] == '-1':
                        final_below_upper[sample].remove(snp_allele)


    return final_below_upper


def mask_program(input_file, upper_min, mfi_cutoff, fraction_cutoff):
    # three passes of masking occur:
    # 1. mask individual snp values in each sample (mask = -1) -- collect list of snps per sample
    #    that are above the mfi_cutoff but below the upper_min.
    # 2. mask entire snps if > fraction_cutoff of the snp are -1. this is done in a diff dictionary organized by snp
    #    instead of sample, so before moving on to next step update the sample based dict
    #    if an entire snp is masked, remove any instances of that snp from the list
    #    that are above the mfi_cutoff but below the upper_min.
    # 3. mask entire sample if > fraction_cutoff of the sample are -1. this is the final dictionary.
    #    if an entire sample is masked, remove any instances of that sample from the list
    #    that are above the mfi_cutoff but below the upper_min
    # write the output file based on this dictionary and program_info.

    if valid_cutoffs(mfi_cutoff, upper_min, fraction_cutoff): # make sure the input from the user is ok
        # mfi_block = list of lines of relevant data, starting with headers 'Location,Sample,snp_names...'
        # program_info = list of all lines in the file leading up to and including the headers
        mfi_block, progam_info = mfi_reader(input_file)
        samples, headers = map_sample_to_snps(mfi_block)
        matches, rs_index = find_matching_alleles(headers)
        below_upper_min = {}
        masked_samples = {}
        final_masked_samples = {}
        for sample in samples: # first pass of masking. per sample, mask individual snps that are below cutoff
            masked_sample, well_samp, below_upper = mask_sample(sample, matches, rs_index, mfi_cutoff, upper_min)
            below_upper_min[well_samp] = below_upper
            masked_samples[well_samp]= masked_sample
        snps = map_snps_to_samples(masked_samples.values(), headers[2:-2])
        for snp in snps.keys(): # second pass of masking. per snp, mask entire snp if > max_fraction have -1 MFI
            snps[snp]= mask_values(snps[snp], fraction_cutoff)
        updated_masked_samples = update_masked_samples(snps, masked_samples) # update the masked samples dict based on snps
        for sample in updated_masked_samples.keys(): # final pass of masking, mask entire sample if > max_fraction have -1 MFI
            final_masked_samples[sample] = mask_values(updated_masked_samples[sample], fraction_cutoff)

        final_below_upper = update_below_upper(final_masked_samples, below_upper_min)
        return final_masked_samples, progam_info, final_below_upper
    else:
        exit(1)


def write_masked_line(sample, headers):
    masked_line= ['']*len(headers)
    for key in sample.keys():
        index = headers.index(key)
        if key == 'Location':
            mfi = int(sample[key])
        else:
            mfi = sample[key]
        masked_line[index]= mfi

    return masked_line


def write_masked_file(final_masked_dict, orig_name, output_folder, program_info, min_mfi, max_fraction,
                      upper_mfi, below_upper):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_name = '%s\\%s_masked_%s_%s.csv' % (output_folder, orig_name.split('\\')[-1].split('.')[0],
                                               str(int(float(min_mfi))), max_fraction.split('.')[1] )

    headers = program_info[-1]
    masked_lines = []
    for sample in final_masked_dict:
        masked_line = write_masked_line(final_masked_dict[sample], headers)
        masked_lines.append(masked_line)
    sorted_list = sorted(masked_lines, key=operator.itemgetter(0))

    with open(output_name,"w") as f:
        for line in program_info:
            f.write(",".join(line) + '\n')
        for mfi_line in sorted_list:
            mfi_line[0] = str(mfi_line[0])
            f.write(",".join(mfi_line) + '\n')
        f.close()

    log_folder = '%s\\below_upper_limit' % output_folder
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    if not all_blank(below_upper.values()): # don't bother writing the file if all samples have empty lists of middling snps
        log_name = '%s\\%s_below_upper_cutoff.txt' % (log_folder, orig_name.split('\\')[-1].split('.')[0] )
        with open(log_name, 'w') as l:
            l.write('%s\nAlleles with MFI > %s but < %s\n' % (log_name, min_mfi, upper_mfi))
            for well_samp in sorted(below_upper.keys()):
                if below_upper[well_samp] != []:
                    l.write('\n' + well_samp[0] + '\t' + well_samp[1])
                    for allele in below_upper[well_samp]:
                        l.write('\t%s' % (str(allele)))
        l.close()

###################################################################################################################

def main(argv):
   inputfolder = ''
   outputfolder = ''
   min_MFI = ''
   max_fraction = ''
   uppermin = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:u:m:f:",["ifolder=","ofolder=","uppermin=","minmfi=","maxfraction="])
   except getopt.GetoptError:
      print '\n\nERROR: widom_mask.py requires 5 arguments: python widom_mask.py -i <inputfolder> -o <outputfolder> -u <uppermin_MFI> -m <min_MFI> -f <max_fraction>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print '\nUsage:\n'
         print 'python widom_mask.py -i <inputfolder>  -o <outputfolder> -u <uppermin_MFI> -m <min_MFI> -f <max_fraction>\n'
         print '\tIf spaces appear in input folder or output file path name'
         print '\t"(ie "...Brzustowicz Lab\Personal Folders...")'
         print '\tenclose entire path name in double quotes'
         sys.exit()
      elif opt in ("-i", "--ifolder"):
        inputfolder = arg
      elif opt in ("-o", "--ofolder"):
        outputfolder = arg
      elif opt in ("-u, --uppermin"):
          uppermin = arg
      elif opt in ("-m", "--minmfi"):
          min_MFI = arg
      elif opt in ("-f", "--maxfraction"):
          max_fraction = arg

   if inputfolder == '' or outputfolder == '' or min_MFI == '' or max_fraction == '' or uppermin == '':
       print '\n\nERROR: widom_mask.py requires 5 arguments: -i <inputfolder> -o <outputfolder> -u <uppermin_MFI> -m <min_MFI> -f <max_fraction>\n '
   else:
       print '\n'
       print 'masking...\n'
       block_files, empty_files = list_folder_files(inputfolder)
       if empty_files != []:
           print 'FYI: Some input files were empty:\n%s\n' % empty_files
       for lum_file in block_files:
            final_masked_dict, program_info, below_upper = mask_program(lum_file, uppermin, min_MFI, max_fraction)
            write_masked_file(final_masked_dict, lum_file, outputfolder, program_info, min_MFI, max_fraction,
                              uppermin, below_upper)
            print lum_file
       print 'done.'

if __name__ == "__main__":
   main(sys.argv[1:])