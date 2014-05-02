## this is a function that just gets the names of all the VCFs that are living in Matt's VCF folder at this time...

def sample_dir_paths():

    import os

    sample_dir_paths = []
    
    mut_dir = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS"
    one = os.walk( mut_dir )
    data = [x[0] for x in os.walk( mut_dir ) ]

    for dir in data:
        if "TCGA" in dir and "FILT" not in dir:
            if dir != "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/MUTATION_CALLS":
                sample_dir_paths.append( dir )

    return sample_dir_paths
