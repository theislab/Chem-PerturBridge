"""Manual PubChem CID overrides for the CIGS dataset."""


def pubchem_mapping_cigs() -> dict:
    return {
        "cigs": {
            # CAS 33289-85-9 is shared by Vindolinine (correct) and Dipsacoside B
            # (wrong CAS in the TCM supplement). Pin Vindolinine's CID via its
            # TCM catalog_id so the CAS-first PubChem lookup never silently
            # mislabels it.
            "Cpd0165": 24148538,
        }
    }
