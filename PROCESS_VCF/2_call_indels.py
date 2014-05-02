# the purpose of this is to call the indels and store the data in a more useable format for mutation calls
# it savesthe data as SAMPLE_normal.txt and SAMPLE_tumor.txt in DATA/INDELS

# for if you lose a file again...
# gpg --no-tty --batch --cipher-algo AES256 --passphrase-file /share/WilkeLab/work/dzd58/SS_passphrase.txt --decrypt TCGA_GBM.tumor.indels.recal.vcf.gpg > TCGA_GBM.tumor.indels.recal.vcf

## FORMAT IS GT:AD:DP:GQ:PL
## FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
## FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
## FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
## FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
## FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">

def set_up( indel_dict, line ):

    sample_order = []

    line = line.strip( '\n\r' )
    fields = line.split( '\t' )
    for j in range( 9, len( fields ) ):
        sample_key = fields[ j ][ 0:17 ]
        indel_dict[ sample_key ] = {}
        sample_order.append( sample_key )

    return sample_order

def process_indel( indel_dict, line ):

    # take it if the quality is good enough
    line = line.strip( '\n\r' )
    fields = line.split( '\t' )
    qual = fields[ 5 ]
    if float( qual ) >= 50:
        chrom = fields[ 0 ]
        pos = fields[ 1 ]
        ref = fields[ 3 ]
        #print chrom, pos

        # since we've taken the indel, which samples have it?
        for j in range( 9, len( fields ) ):
            if fields[ j ] == "./.":
                continue
            else:
                sample_key = sample_order[ j-9 ]
                if sample_key not in indel_dict:
                    print "FAIL! Did not fine %s in indel_dict" % sample_key
                else:
                    if chrom in indel_dict[ sample_key ]:
                        indel_dict[ sample_key ][ chrom ].append( pos )
                    else:
                        indel_dict[ sample_key ][ chrom ] = []
                        indel_dict[ sample_key ][ chrom ].append( pos )

                # this is just input about the indels that we don't wind up needing to locate them...
                # format = fields[ j ]
                # f_fields = format.split( ':' )
                # genotype = f_fields[ 0 ]
                # allelic_depths =
                # read_depths 
                # print sample_key# , format
                
    print "done with %s" % i
    #print indel_dict

    return

def save_output( indel_dict, tissue ):

    for sample in indel_dict:
        new_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/INDELS/%s_%s.txt" % ( sample, tissue )
        new_fh = open( new_file, 'w' )
        for chrom in indel_dict[ sample ]:
            for loc in indel_dict[ sample ][ chrom ]:
                new_fh.write( "%s\t%s\n" % ( chrom, loc ) )
        new_fh.close()

    return

###################
## MAIN FUNCTION ##
###################

import os.path

gatk_indel_dir = "/share/WilkeLab/work/MattCC/TCGA_NextGen_Data/GATK_Indel_Calls_FINAL"

tumor_indel_file = os.path.join( gatk_indel_dir, "TCGA_GBM.tumor.indels.recal.vcf" )
blood_indel_file = os.path.join( gatk_indel_dir, "TCGA_GBM.blood.indels.recal.vcf" )

file_list = [ tumor_indel_file, blood_indel_file ]
# file_list = [ blood_indel_file ]

for file in file_list:
    tissue = file[ -22:-17 ]
    indel_dict = {}

    # open the file, and enumerate the lines, so that the whole thing doesn't go into memory (it's too big)
    file_fh = open( file, 'r' )
    for i, line in enumerate( file_fh ):

        if i == 64:
            sample_order = set_up( indel_dict, line )

        if i > 64:
        #if i == 65:
            process_indel( indel_dict, line )
            
    # save as individual sample/tissue files
    save_output( indel_dict, tissue )

## i forgot about the order; order matters of course, which is so, so, so stupid
## sleeping on this problem, but...
## maybe lists? make a list of samples that is in order...
## and use that to find the right sample for each line... 
