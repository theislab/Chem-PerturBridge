import re
from typing import Any, Dict, Optional

import requests

import pandas as pd

from src.utils.parsing_utils import *

def _fetch_cellosaurus_details(cvcl_id: str, n_retries: int = 3) -> Optional[Dict[str, Any]]:
    """Fetch detailed information for a specific Cellosaurus ID."""
    def _get_details(cvcl_id: str) -> Dict[str, Any]:
        url = f"https://api.cellosaurus.org/cell-line/{cvcl_id}"
        params = {"format": "json"}
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        return r.json()

    return _call_with_retry(_get_details, cvcl_id, n_retries=n_retries, sleep_s=2)

def _call_with_retry(func, *args, n_retries: int = 3, sleep_s: int = 2, **kwargs):
    last_exc = None
    for attempt in range(n_retries):
        try:
            return func(*args, **kwargs)
        except Exception as exc:
            last_exc = exc
            if attempt + 1 < n_retries:
                time.sleep(sleep_s)
    return None


def _parse_age_value(age_value: str) -> Any:
    """
    Parse age string and extract numeric value in years.

    Cellosaurus age format: "54Y", "30Y6M", "3M", "45D", "Fetus", "Age unspecified"
    where Y=years, M=months, D=days
    """
    if not age_value or age_value == "Age unspecified":
        return "-666"

    if age_value.lower() in ["fetus", "embryo", "newborn"]:
        return age_value

    try:
        year_match = re.search(r"(\d+)Y", age_value)
        if year_match:
            return year_match.group(1)

        if "Y" not in age_value and (
            re.search(r"\d+[MD]", age_value) or re.search(r"<.*?[MD]", age_value)
        ):
            return 0
    except (ValueError, AttributeError):
        pass

    return "-666"


def _map_gender_value(gender_value: Optional[str]) -> str:
    """Map gender value to single letter code."""
    if not gender_value or gender_value == "Sex unspecified":
        return "-666"

    gender_lower = gender_value.lower()
    if gender_lower == "female":
        return "F"
    if gender_lower == "male":
        return "M"
    return gender_value


def _extract_ethnicity_from_comments(comment_list: list) -> str:
    """Extract ethnicity/population from Cellosaurus comment list."""
    if not comment_list:
        return "-666"

    for comment in comment_list:
        if comment.get("category") == "Population":
            return comment.get("value", "-666")

    return "-666"


def _extract_donor_info(data: Dict[str, Any], progress: bool) -> Dict[str, Any]:
    """Extract donor age, sex, and ethnicity from Cellosaurus API response."""
    donor_info = {
        "donor_age": "-666",
        "donor_sex": "-666",
        "donor_ethnicity": "-666",
    }

    if not data:
        return donor_info

    try:
        cellosaurus = data.get("Cellosaurus", {})
        cell_line_list = cellosaurus.get("cell-line-list", [])

        if not cell_line_list:
            return donor_info

        cell_line = cell_line_list[0]
        donor_info["donor_age"] = _parse_age_value(cell_line.get("age", "Age unspecified"))
        donor_info["donor_sex"] = _map_gender_value(cell_line.get("sex", "Sex unspecified"))
        donor_info["donor_ethnicity"] = _extract_ethnicity_from_comments(
            cell_line.get("comment-list", [])
        )
    except Exception as e:
        if progress:
            logger.warning(f"Error extracting donor info: {e}")

    return donor_info


def fetch_donor_info_from_cellosaurus(
    df: pd.DataFrame,
    cell_id_col: str = "cell_id",
    base_cell_id_col: str = "base_cell_id",
    cellosaurus_id_col: str = "cellosaurus_id",
    n_retries: int = 3,
    progress: bool = True,
) -> pd.DataFrame:
    """
    Fetch donor information from Cellosaurus API and add to dataframe.

    Parameters:
    -----------
    df : pd.DataFrame
        Annotated cell line dataframe with cellosaurus_id column
    cell_id_col : str, default="cell_id"
        Column containing cell line identifiers for filtering
    base_cell_id_col : str, default="base_cell_id"
        Column containing base cell line names
    cellosaurus_id_col : str, default="cellosaurus_id"
        Column containing Cellosaurus IDs (CVCL_xxxx)
    n_retries : int, default=3
        Number of retries for API calls
    progress : bool, default=True
        Whether to log progress information

    Returns:
    --------
    pd.DataFrame
        Original dataframe with added columns: donor_age, donor_sex, donor_ethnicity
    """
    
    # Filter dataframe where cell_id == base_cell_id
    filtered_df = df[df[cell_id_col] == df[base_cell_id_col]].copy()
    
    # Get unique cellosaurus IDs that need donor info
    unique_cvcl_ids = filtered_df[cellosaurus_id_col].dropna().unique()
    n_ids = len(unique_cvcl_ids)
    
    if progress:
        logger.info(f"Fetching donor information for {n_ids} unique Cellosaurus IDs")
    
    # Fetch donor info for each unique CVCL ID
    cvcl_to_donor = {}
    iteration_count = 0
    
    for cvcl_id in unique_cvcl_ids:
        if pd.isna(cvcl_id) or not str(cvcl_id).strip():
            continue
        
        cvcl_id = str(cvcl_id).strip()
        
        # Fetch details from API
        details = _fetch_cellosaurus_details(cvcl_id, n_retries=n_retries)
        donor_info = _extract_donor_info(details, progress=progress)
        cvcl_to_donor[cvcl_id] = donor_info
        
        # Progress logging
        iteration_count += 1
        if progress and iteration_count % 10 == 0:
            n_fetched = sum(1 for info in cvcl_to_donor.values() if info["donor_age"] != "-666")
            logger.info(f"Processed {iteration_count}/{n_ids} Cellosaurus IDs ({n_fetched} with donor info)")
    
    # Final progress summary
    if progress:
        n_with_info = sum(1 for info in cvcl_to_donor.values() if info["donor_age"] != "-666")
        logger.info(f"Completed: {n_with_info}/{n_ids} IDs have donor information")
    
    # Add donor info to filtered dataframe with temporary column names
    filtered_df["donor_age_new"] = filtered_df[cellosaurus_id_col].map(
        lambda x: cvcl_to_donor.get(x, {}).get("donor_age", "-666") if pd.notna(x) else "-666"
    )
    filtered_df["donor_sex_new"] = filtered_df[cellosaurus_id_col].map(
        lambda x: cvcl_to_donor.get(x, {}).get("donor_sex", "-666") if pd.notna(x) else "-666"
    )
    filtered_df["donor_ethnicity_new"] = filtered_df[cellosaurus_id_col].map(
        lambda x: cvcl_to_donor.get(x, {}).get("donor_ethnicity", "-666") if pd.notna(x) else "-666"
    )
    
    # Create a mapping from the filtered df with new column names
    filtered_subset = filtered_df[[cell_id_col, "donor_age_new", "donor_sex_new", "donor_ethnicity_new"]].copy()
    
    # Merge the donor info back
    result_df = df.merge(
        filtered_subset,
        left_on=base_cell_id_col,
        right_on=cell_id_col,
        how="left",
        suffixes=('', '_donor_annotation')
    )
    
    # Fill NaN values in new columns with defaults
    result_df["donor_age_new"] = result_df["donor_age_new"].fillna("-666")
    result_df["donor_sex_new"] = result_df["donor_sex_new"].fillna("-666")
    result_df["donor_ethnicity_new"] = result_df["donor_ethnicity_new"].fillna("-666")
    
    # Update existing columns or create new ones, only filling empty values
    # For donor_age
    result_df['donor_age'] = result_df['donor_age'].astype(str)
    if "donor_age" in result_df.columns:
        mask = result_df["donor_age_new"] != "-666"
        result_df.loc[mask, "donor_age"] = result_df.loc[mask, "donor_age_new"]
    else:
        result_df["donor_age"] = result_df["donor_age_new"]
    
    # For donor_sex
    if "donor_sex" in result_df.columns:
        mask = result_df["donor_sex_new"] != "-666"
        result_df.loc[mask, "donor_sex"] = result_df.loc[mask, "donor_sex_new"]
    else:
        result_df["donor_sex"] = result_df["donor_sex_new"]
    
    # For donor_ethnicity
    if "donor_ethnicity" in result_df.columns:
        # Only update where existing value is "-666" (empty)
        mask = result_df["donor_ethnicity"] == "-666"
        result_df.loc[mask, "donor_ethnicity"] = result_df.loc[mask, "donor_ethnicity_new"]
    else:
        result_df["donor_ethnicity"] = result_df["donor_ethnicity_new"]
    
    # Drop temporary columns
    result_df = result_df.drop(columns=["donor_age_new", "donor_sex_new", "donor_ethnicity_new"])
    
    return result_df, filtered_subset
