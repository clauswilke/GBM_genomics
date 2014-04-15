# the purpose of this script is to parse the doubles_files list, to see where each pair of samples differs

def sort_doubles():

    doubles_dict = {}

    # look at the file in ../DATA with the double files
    doubles_file = "/share/WilkeLab/work/dzd58/TCGA_Reanalysis/DATA/doubles.txt"
    doubles_fh = open( doubles_file, 'r' )
    data = doubles_fh.readlines()

    # now the parsing begins
    for line in data:
        line = line.strip( '\n' )
        #print line
        fields = line.split( '\t' )
        sample = fields[ 0 ]
        # file one, info extraction
        file_one = fields[ 1 ]
        prefix_1 = file_one[ :4 ]
        suffix_1 = file_one[ -2: ]
        barcode_1 = file_one[ 5:-2 ]
        fields_1 = barcode_1.split( '-' )
        study_1 = fields_1[ 0 ]
        tissue_source_site_1 = fields_1[ 1 ]
        participant_1 = fields_1[ 2 ]
        sample_type_1 = fields_1[ 3 ][ 0:2 ]
        vial_1 = fields_1[ 3 ][ 2 ]
        portion_1 = fields_1[ 4 ][ 0:2 ]
        analyte_1 = fields_1[ 4 ][ 2 ]
        if len( fields_1 ) > 5:
            plate_1 = fields_1[ 5 ]
            center_1 = fields_1[ 6 ]
        else:
            plate_1 = "BLANK"
            center_1 = "BLANK"
        # file two, info extraction       
        file_two = fields[ 2 ]
        prefix_2 = file_two[ :4 ]
        suffix_2 = file_two[ -2: ]
        barcode_2 = file_two[ 5:-2 ]
        fields_2 = barcode_2.split('-' )
        study_2 = fields_2[ 0 ]
        tissue_source_site_2 = fields_2[ 1 ]
        participant_2 = fields_2[ 2 ]
        sample_type_2 = fields_2[ 3 ][ 0:2 ]
        vial_2 = fields_2[ 3 ][ 2 ]
        portion_2 = fields_2[ 4 ][ 0:2 ]
        analyte_2 = fields_2[ 4 ][ 2 ]
        if len( fields_2 ) > 5:
            plate_2 = fields_2[ 5 ]
            center_2 = fields_2[ 6 ]
        else:
            plate_2 = "BLANK"
            center_2 = "BLANK"
        # compare the two
        if study_1 == study_2 and tissue_source_site_1 == tissue_source_site_2 and participant_1 == participant_2 and sample_type_1 == sample_type_2 == "01": ## where '01' is the code for a tumor
            print "good to go"
            if vial_1 == vial_2:
                VIAL = 1
            else:
                VIAL = 0
            if portion_1 == portion_2:
                PORT = 1
            else:
                PORT = 0
            if analyte_1 == analyte_2:
                ANAL = 1
            else:
                ANAL = 0
                # print analyte_1, analyte_2
            if plate_1 == plate_2:
                PLATE = 1
            else:
                PLATE = 0
            if center_1 == center_2:
                CENT = 1
            else:
                CENT = 0
            if prefix_1 == prefix_2:
                PREF = 1
            else:
                PREF = 0
                #print prefix_1, prefix_2
            if suffix_1 == suffix_2:
                SUFF = 1
            else:
                SUFF = 0
                #print suffix_1, suffix_2
        # print sample, VIAL, PORT, ANAL, PLATE, CENT, PREF, SUFF


        # store the information
    
    
    #print "this function is running"

    return doubles_dict
