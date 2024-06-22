# Uniprot meta tool

this tool is designed to be useful and simple interface for [UniProt.org](https://uniprot.org)
knowlege base witch can  be used to simply get protein metadata by its ID

**NOTE:** this package is in early stage of development so be careful!

## Examples

```
from uniprot_meta_tool import UniprotData

>>>up_data = UniprotData('P05067')
>>>up_data.name
'Amyloid-beta precursor protein'
>>>up_data.genes[0]
'APP'
>>>up_data.organism.common_name
'Human'
```

some util functions are also available:

```
get_molecular_functions(uniprot_data) -> List[str]:
returns list of molecular functions from keywords. accepts uniprot_id as string or UniprotData object

get_biological_processes(uniprot_data) -> List[str]:
returns list of biological processes from keywords.

get_pathways(uniprot_data: UniprotData) -> List[str]:
returns list of pathways (from Reactome DB using UniProt cross-reference)    
```

## Raw search

You can reach raw metadata from `UniprotData(u_id).raw_data` property.

You can use [JMESPath](https://jmespath.org/) query language to get something from raw data using 
`UniprotData(u_id).search()` method


## Configuration

You can set up basic parameters of how uniprot_meta_tool works with remote API
and local storage, but be careful with this

```python
from uniprot_meta_tool.data_parser import set_config
set_config(
    meta_cache_dir='/your/local/folder', # default is %tmp%/uniprot_meta
    use_cache=True, # if you set it to False, module will not store retrieved info to your disk 
    #                 and will perform new HTTP request every time you create UniprotData object
    uniprot_rest_url='https://rest.uniprot.org/uniprotkb/{}' # have no idea why you could want to change it, but you can
)
```

## MR and feature requests

Feel free to suggest features using GitHub issues and to implement something you
need if you want to!