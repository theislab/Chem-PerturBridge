def ensembl_mapping_op3():
    """
    Manual mapping for OP3 gene symbols to Ensembl IDs where automatic lookup may fail.
    
    The annotation is based on https://www.genenames.org/ and 
    https://www.ncbi.nlm.nih.gov/gene/ sources for
    RPL9P9, C18ORF21, TTTY15 and LINC-PINT-1
    the rest of the genes are annotated using GSE279945_multiome_counts_processed.h5ad from
    GEO accession GSE279945. The mapping is constructed for genes which have None
    Ensembl IDs after the HGNC API lookup.

    Returns
    -------
    dict
        Dictionary with 'op3' key containing gene symbol to Ensembl ID mappings
    """
    gene2ensg = {
        'op3': {
            'RPL9P9': 'ENSG00000293320',
            'C18ORF21': 'ENSG00000141428',
            'TTTY15': 'ENSG00000114374',
            'LINC-PINT-1': 'ENSG00000231721',
            'AC004687.1': 'ENSG00000265206',
            'AC004865.2': 'ENSG00000239636',
            'AC005332.4': 'ENSG00000274712',
            'AC008105.3': 'ENSG00000267121',
            'AC008124.1': 'ENSG00000273015',
            'AC008969.1': 'ENSG00000176593',
            'AC010642.2': 'ENSG00000283103',
            'AC012146.1': 'ENSG00000234327',
            'AC020656.1': 'ENSG00000257764',
            'AC025164.1': 'ENSG00000245904',
            'AC027644.3': 'ENSG00000272831',
            'AC058791.1': 'ENSG00000273319',
            'AC068587.4': 'ENSG00000283674',
            'AC087190.1': 'ENSG00000260349',
            'AC093323.1': 'ENSG00000170846',
            'AC100810.1': 'ENSG00000253982',
            'AC114760.2': 'ENSG00000272211',
            'AC116407.2': 'ENSG00000277511',
            'AC127502.2': 'ENSG00000270055',
            'AC243829.5': 'ENSG00000278097',
            'AC243960.1': 'ENSG00000268027',
            'AC245060.5': 'ENSG00000274422',
            'AC245297.3': 'ENSG00000274265',
            'AL133342.1': 'ENSG00000278231',
            'AL133453.1': 'ENSG00000258757',
            'AL139246.5': 'ENSG00000272449',
            'AL662844.4': 'ENSG00000272501',
            'FO393401.1': 'ENSG00000230155',
            'MATR3-1': 'ENSG00000015479',
            
        }
    }
    return gene2ensg
