from typing import List, Literal
import logging
import pandas as pd
from uniprot_meta_tool import UniprotData, get_pathways, get_molecular_functions, get_biological_processes

enrich_params = {
    'protein_name': lambda uniprot_id: UniprotData(uniprot_id).name,
    'gene': lambda uniprot_id: UniprotData(uniprot_id).genes[0],
    'organism': lambda uniprot_id: UniprotData(uniprot_id).organism.scientific_name,
    'sequence': lambda uniprot_id: UniprotData(uniprot_id).sequence,
    'molecular_functions': lambda uniprot_id: '; '.join(get_molecular_functions(uniprot_id)),
    'biological_processes': lambda uniprot_id: '; '.join(get_biological_processes(uniprot_id)),
    'pathways': lambda uniprot_id: '; '.join(get_pathways(uniprot_id))
}


def wrap(fn, value):
    try:
        return fn(value)
    except Exception as e:
        logging.error(value)
        logging.exception(e)
        return None


def enrich_df(
        df: pd.DataFrame,
        id_column: str,
        columns: List[Literal[
            'protein_name', 'gene', 'organism',
            'sequence', 'molecular_functions', 'biological_processes',
            'pathways'
        ]]
) -> pd.DataFrame:
    for col in columns:
        df[col] = df[id_column].apply(lambda x: wrap(enrich_params[col], x))
    return df
