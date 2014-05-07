## these functions make a list of paths to be fed into some analysis pipeline; the four current functions make the following lists:
##
## existing_VCFs: VCFs that exist and are in Matt's folder of unencrypted VCFs on phylocluster
## existing_VCF_doubles: from the list of "existing_VCFs," those samples that appear twice
## analyzed_VCFs: VCFs that have been filtered, and whose data exist in my DATA/MUTATION_CALLS/ folder
## analyzed_VCF_doubles: 
##

def analyzed_VCFs():

    import os

    analyzed_VCFs = []
    
    mut_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS"
    one = os.walk( mut_dir )
    data = [x[0] for x in os.walk( mut_dir ) ]

    for dir in data:
        if "TCGA" in dir and "FILT" not in dir:
            if dir != "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS":
                analyzed_VCFs.append( dir )

    return analyzed_VCFs

def analyzed_doubles():

    # this part is the same as the existing VCFs function above
    import os

    mut_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS"
    one = os.walk( mut_dir )
    data = [x[0] for x in os.walk( mut_dir ) ]

    # this is the part that is different
    analyzed_doubles = {}
    doub_count = 0

    samples = {}
    for dir in data:
        if "TCGA" in dir and "FILT" not in dir:
            if dir != "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS":
                sample = dir[ -12: ]
                if sample not in samples:
                    samples[ sample ] = [ dir ]
                else:
                    samples[ sample ].append( dir )

    for sample in samples:
        if len( samples[ sample ] ) >1:
            doub_count += 1
            analyzed_doubles[ sample ] = sorted( samples[ sample ] )

    return analyzed_doubles, doub_count

def existing_VCFs():

    existing_VCFs = []

    return existing_VCFs

def existing_doubles():

    existing_doubles = []

    return existing_doubles


## MAIN TEST
# analyzed_doubles = analyzed_doubles()
# print analyzed_doubles
