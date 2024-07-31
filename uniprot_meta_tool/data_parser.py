import os.path
import json
import gzip
from typing import List, Optional, Dict
from dataclasses import dataclass
import tempfile

import requests
import jmespath

config = {
    'meta_cache_dir': os.path.join(tempfile.gettempdir(), 'uniprot_meta'),
    'use_cache': True,
    'uniprot_rest_url': 'https://rest.uniprot.org/uniprotkb/{}',
    'gzip_cache': False,
    'list_separator': '; '

}
os.makedirs(config['meta_cache_dir'], exist_ok=True)


class ObsoleteException(Exception):
    pass


def set_config(meta_cache_dir=None, uniprot_rest_url=None, use_cache=None, gzip_cache=None, list_separator=None):
    if meta_cache_dir is not None:
        config['meta_cache_dir'] = meta_cache_dir
        os.makedirs(meta_cache_dir, exist_ok=True)
    if uniprot_rest_url is not None:
        config['uniprot_rest_url'] = uniprot_rest_url
    if use_cache is not None:
        config['use_cache'] = use_cache
    if gzip_cache is not None:
        config['gzip_cache'] = bool(gzip_cache)
    if list_separator is not None:
        config['list_separator'] = str(list_separator)




@dataclass
class PDBEntry:
    pdb_id: str
    method: Optional[str] = None
    resolution: Optional[float] = None
    chains: Optional[str] = None


@dataclass
class Feature:
    feature_type: str
    feature_id: Optional[str]
    description: Optional[str]
    start: int
    end: int


@dataclass
class Organism:
    scientific_name: str
    common_name: str
    taxon_id: str


class UniprotData:
    uniprot_id: Optional[str] = None
    raw_data: Optional[Dict] = None
    name: Optional[str] = None
    alt_names: List[str] = []
    genes: List[str] = []
    gene: str = ''
    pdb_ids: List[str] = None
    pdb_entries: List[PDBEntry] = None
    cross_references: Optional[Dict] = None
    keywords: Optional[Dict] = None
    organism: Optional[Organism] = None
    features: List[Feature] = None
    sequence: Optional[str] = ''
    comments: Optional[List[Dict]] = None
    subcellular_locations: Optional[List[str]]

    def search(self, query):
        return jmespath.search(query, self.raw_data)

    def get_position_features(self, start: int, end: Optional[int] = None, inner_only=False) -> List[Feature]:
        """
        get features over given position in protein sequence
        :param start: start of the inspected sequence. if end is not given, looks for features in this specific point
        :param end: end of the inspected sequence
        :param inner_only: if true, returns only features contained between given start and end;
            otherwise returns overlapping features. Ignored if no end value provided
        :return: list of Feature found for sequence
        """
        if end is None:
            return [x for x in self.features if x.start <= start <= x.end]
        else:
            if inner_only:
                return [x for x in self.features if x.start >= start and x.end <= end]
            else:
                return [x for x in self.features if x.start <= start <= x.end and x.start <= end <= x.end]

    def __init__(self, uniprot_id, raw_meta=None):
        cache_file = CacheFile(uniprot_id)
        if raw_meta is not None:
            self.raw_data = raw_meta
        else:
            if not cache_file.exists:
                req = requests.get(config['uniprot_rest_url'].format(uniprot_id))
                if req.status_code == 200:
                    self.raw_data = req.json()
                else:
                    raise requests.exceptions.RequestException(
                        f'UniProt API Returned {req.status_code}',
                        response=req
                    )
                if self.raw_data is not None:
                    cache_file.write(self.raw_data)
            else:
                self.raw_data = cache_file.read()
        if self.search('entryType') == 'Inactive':
            raise ObsoleteException(f'Entry {uniprot_id} is obsolete!')
        self.uniprot_id = uniprot_id
        self.name = self.search('proteinDescription.recommendedName.fullName.value')
        self.alt_names = self.search('proteinDescription.alternativeNames[*].fullName.value')
        self.pdb_ids = self.search("uniProtKBCrossReferences[?database=='PDB'].id")
        self.genes = self.search("genes[*].geneName.value")
        if self.genes is None:
            self.genes = []
        else:
            if len(self.genes) > 0:
                self.gene = self.genes[0]
        self.sequence = self.search('sequence.value')
        # proces pdbentries
        pdb_entry_list = self.search(
            "uniProtKBCrossReferences[?database=='PDB']."
            "{pdb_id: id, resolution:properties[?key=='Resolution'].value|[0],"
            "method:properties[?key=='Method'].value|[0],"
            "chains:properties[?key=='Chains'].value|[0]}"
        )
        self.pdb_entries = []
        if pdb_entry_list is not None:
            for entry in pdb_entry_list:
                try:
                    entry['resolution'] = float(entry['resolution'].replace('A', ''))
                except ValueError:
                    entry['resolution'] = None
                self.pdb_entries.append(PDBEntry(**entry))
        # process features
        try:
            self.features = [Feature(**x) for x in self.search(
                "features[*].{feature_type: type, start: location.start.value, end: location.end.value,"
                "description: description, feature_id: featureId}"
            )]
        except TypeError:
            print(uniprot_id)
            print(self.raw_data)
            raise
        self.comments = self.search('comments')
        self.subcellular_locations = self.search(
            "comments[?commentType=='SUBCELLULAR LOCATION'].subcellularLocations[*].location.value|[0]"
        )
        if self.subcellular_locations is None:
            self.subcellular_locations = []
        # process organism
        self.organism = Organism(
            scientific_name=self.search('organism.scientificName'),
            common_name=self.search('organism.commonName'),
            taxon_id=self.search('organism.taxonId')
        )
        # process keywords
        self.keywords = dict()
        for kw in self.search('keywords'):
            category = kw.get('category', None)
            if not category:
                category = 'uncategorized'
            if category not in self.keywords.keys():
                self.keywords[category] = []
            self.keywords[category].append(kw['name'])

        # process crossReferences
        self.cross_references = dict()
        for reference in self.search('uniProtKBCrossReferences'):
            db = reference.get('database', None)
            if db is None:
                db = 'unknown'
            if db not in self.cross_references.keys():
                self.cross_references[db] = []
            refer = {
                'id': reference.get('id', ''),
            }
            refer.update({x['key']: x['value'] for x in reference.get('properties')})
            self.cross_references[db].append(refer)


class CacheFile:
    path: Optional[str] = None
    uniprot_id: Optional[str] = None

    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        if not config['gzip_cache']:
            self.path = os.path.join(config['meta_cache_dir'], f"{uniprot_id}.meta.json")
        else:
            self.path = os.path.join(config['meta_cache_dir'], f"{uniprot_id}.meta.json.gz")

    @property
    def exists(self) -> bool:
        return os.path.exists(self.path)

    def read(self):
        if not config['gzip_cache']:
            with open(self.path) as in_file:
                data = json.loads(in_file.read())
        else:
            with gzip.open(self.path, 'rb') as file:
                data = json.loads(file.read().decode('utf-8'))
        return data

    def write(self, data):
        if not config['gzip_cache']:
            with open(self.path, 'w') as out_file:
                out_file.write(json.dumps(data))
        else:
            with gzip.open(self.path, 'wb') as out_file:
                out_file.write(json.dumps(data).encode('utf-8'))