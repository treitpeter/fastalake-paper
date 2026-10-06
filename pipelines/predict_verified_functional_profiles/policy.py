RANKS=('domain','phylum','class','order','family','genus','species')
def taxonomy(old,parents):
    missing=int(old['missing_parents']);known=list(old['known_parent_lca']);annotated=parents-missing
    assert parents>0 and 0<=missing<=parents and len(known)<=7
    assert not known or annotated>0
    if known:state='annotated_parent_lca' if missing else 'supported_lca'
    else:state='root_only' if annotated else 'unannotated_taxonomy'
    result=dict(state=state,supported_lineage=known,known_parent_lca=known,rank=RANKS[len(known)-1] if known else 'unresolved',missing_parents=missing,annotated_parents=annotated,annotation_complete=missing==0)
    if not missing:assert result['supported_lineage']==old['supported_lineage']
    return result
def tests():
    def old(known,missing):return dict(known_parent_lca=known,missing_parents=missing,supported_lineage=known if not missing else [])
    s=['d__Bacteria','p__P','c__C','o__O','f__F','g__G','s__G s']
    assert taxonomy(old(s,1),2)['rank']=='species'
    assert taxonomy(old(s[:5],1),3)['rank']=='family'
    assert taxonomy(old([],3),3)['state']=='unannotated_taxonomy'
    assert taxonomy(old([],1),3)['state']=='root_only'
    assert taxonomy(old(s,0),2)['state']=='supported_lca'
