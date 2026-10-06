"""Conservative evidence projections. Intensities never choose annotations."""
import hashlib
import json
import math
import re
from dataclasses import dataclass, field

LEVELS = ('KO', 'OG', 'Pathway')
RANKS = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
PREFIXES = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')

def canonical(sequence):
    return sequence.upper().replace('I', 'L')

def parse_terms(value, level):
    if level == 'KO':
        return tuple(sorted(set(re.findall(r'K\d{5}', value or ''))))
    if level == 'OG':
        return tuple(sorted(set(re.findall(r'[^;,]+@(?:2\|Bacteria|2157\|Archaea)', value or ''))))
    if level == 'Pathway':
        # Include every map term. Overlap is retained, never equal-split.
        return tuple(sorted(set(re.findall(r'map\d{5}', value or ''))))
    raise ValueError(level)

def lineage(text):
    tokens = dict((v[:3], v[3:]) for v in (text or '').split(';') if len(v) > 3)
    result = []
    for prefix in PREFIXES:
        if not tokens.get(prefix):
            break
        result.append(prefix + tokens[prefix])
    return tuple(result)

def common_prefix(a, b):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    return a[:n]

def eggnog_domains(raw_og):
    """Only explicit domain-level OG tokens; no species inference from homologs."""
    known={'2':'Bacteria','2157':'Archaea','2759':'Eukaryota'}
    return frozenset(known[x] for x in re.findall(r'@(2|2157|2759)\|',raw_og or ''))

def taxonomy_check(tax_lineage, domains):
    if not tax_lineage:
        return 'no_resolved_lineage'
    if not domains:
        return 'no_eggnog_domain_annotation'
    if len(domains)>1:
        return 'domain_annotation_disagreement_withheld'
    return 'domain_compatible' if tax_lineage[0][3:] in domains else 'domain_conflict_withheld'

@dataclass
class FunctionEvidence:
    variants: set = field(default_factory=set)
    missing: bool = False

    def add(self, alternatives):
        known={tuple(sorted(set(terms))) for terms in alternatives if terms}
        # Different annotation sources describe the SAME parent. A missing field
        # in one source does not erase a positive annotation in another source.
        if not known:
            self.missing = True
        self.variants.update(known)

    def result(self):
        variants = sorted(self.variants)
        union = sorted({t for terms in variants for t in terms})
        if not variants:
            state = 'unannotated'
        elif self.missing:
            state = 'incomplete'
        elif len(variants) > 1:
            state = 'conflicting_annotation_sets'
        elif len(variants[0]) == 1:
            state = 'single_consistent'
        else:
            state = 'consistent_multiple_membership'
        # Only memberships present for every possible annotation are guaranteed
        # within this reference; these are not a disjoint signal partition.
        common = sorted(set.intersection(*(set(x) for x in variants))) if variants and not self.missing else []
        return dict(state=state, terms=union, annotation_sets=variants,
                    missing=self.missing, common_memberships=common)

@dataclass
class PeptideEvidence:
    parents: int = 0
    missing_taxonomy: int = 0
    lca: tuple | None = None
    functions: dict = field(default_factory=lambda: {level: FunctionEvidence() for level in LEVELS})

    def add(self, tax_lineage, annotations):
        self.parents += 1
        if not tax_lineage:
            self.missing_taxonomy += 1
        else:
            self.lca = tuple(tax_lineage) if self.lca is None else common_prefix(self.lca, tax_lineage)
        for level in LEVELS:
            self.functions[level].add(annotations.get(level, ()))

    def taxonomy(self):
        known = self.lca or ()
        if not self.parents:
            state, supported = 'no_parent', ()
        elif self.missing_taxonomy:
            state, supported = 'incomplete_taxonomy', ()
        elif not known:
            state, supported = 'root_only', ()
        else:
            state, supported = 'supported_lca', known
        return dict(state=state, supported_lineage=list(supported),
                    rank=RANKS[len(supported)-1] if supported else 'unresolved',
                    known_parent_lca=list(known), missing_parents=self.missing_taxonomy)

    def result(self):
        return dict(parents=self.parents, taxonomy=self.taxonomy(),
                    functions={level: self.functions[level].result() for level in LEVELS})

def signature(level, functional, taxonomic):
    record = dict(level=level, function=functional, taxonomy=taxonomic)
    # Parent counts are evidence metadata, not biological feature identities.
    record['taxonomy'] = {k: taxonomic[k] for k in ('state', 'supported_lineage', 'rank')}
    body = json.dumps(record, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(body.encode()).hexdigest(), record

def valid_intensity(value):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError('Accepted feature intensity must be positive and finite')
    return value
