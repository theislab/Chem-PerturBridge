def pubchem_mapping_l1000():
    """
    Manual mapping for L1000 compounds where automatic PubChem lookup may fail.
    
    This provides fallback mappings for:
    - Ambiguous compound names
    - Non-standard nomenclature
    - Known problematic lookups
    
    Returns
    -------
    dict
        Dictionary with 'l1000' key containing compound name to PubChem CID mappings
    """
    sm2pubchem = {
            'l1000': {

                        'BRD-U19872303': 5702243, 
                        'GNF2': 5311510, 
                        'CMAP-ERG-AZD2281': 23725625
                        
            }
    }
    return sm2pubchem