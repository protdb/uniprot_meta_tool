from uniprot_meta_tool.data_parser import UniprotData
from typing import Union, List


def _transform_from_string(fn):
    def wrapper(uniprot_id: Union[str, UniprotData]):
        print(uniprot_id)
        if type(uniprot_id) is str:
            return fn(UniprotData(uniprot_id))
        if type(uniprot_id) is UniprotData:
            return fn(uniprot_id)
        raise TypeError('invalid parameter passed!')
    return wrapper


@_transform_from_string
def get_molecular_functions(uniprot_data: UniprotData) -> List[str]:
    print(type(uniprot_data))
    return uniprot_data.keywords.get('Molecular function', [])


@_transform_from_string
def get_biological_processes(uniprot_data: UniprotData) -> List[str]:
    return uniprot_data.keywords.get('Biological process', [])


@_transform_from_string
def get_pathways(uniprot_data: UniprotData) -> List[str]:
    pathways = []
    pathway_list = [x.get('PathwayName', '') for x in uniprot_data.cross_references.get('Reactome', [])]
    for pw in pathway_list:
        if pw not in pathways and pw != '':
            pathways.append(pw)
    return pathways
