# Uniprot meta tool

this tool is designed to be useful and simple interface for [UniProt.org](https://uniprot.org)
knowlege base witch can  be used to simply get protein metadata by its ID

**NOTE:** this package is in early stage of development so be careful!

## installation

as simple as possible:

```shell
pip install uniprot-meta-tool
```

only dependencies are `requests` and `jmespath`

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

utility to enrich dataframes (requires Pandas to be installed):

```python
from uniprot_meta_tool.dataframe import enrich_df

df = enrich_df(
    df, # pandas DataFrame
    'protein_id', # name of column contains UniProt ID
    ['gene', 'protein_name', 'organism'] # list of parameters that should be added
)
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
    gzip_cache=True, # Cache files in Gzip archives
    uniprot_rest_url='https://rest.uniprot.org/uniprotkb/{}', # have no idea why you could want to change it, but you can
    list_separator='; ' # affects dataframe tools; for list values (pathways, processes, etc) value in the df will be
    #                     joined with given string
)
```

### Configuration using environment variables

Designed specially to use in Docker containers. You can set this values in Dockerfile or with docker run command if 
your container uses this module.

```shell
UNIPROT_META_CACHE_DIR=/path/to/directory
UNIPROT_USE_CACHE=true
UNIPROT_GZIP_CACHE=true
UNIPROT_REST_URL=https://rest.uniprot.org/uniprotkb/{}
```

## MR and feature requests

Feel free to suggest features using GitHub issues and to implement something you
need if you want to!