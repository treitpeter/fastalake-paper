"""Source classification of protein accessions for the A20 to A28 chain: host (human) versus microbial membership, with accessions of the form PREFIX_<64 hex> treated as unclassifiable because a hash accession carries no source. classify_members returns the class of a set of parent accessions."""
import re
_HEX64=re.compile(r"^[0-9a-fA-F]{64}$")
def is_classifiable_accession(member: str) -> bool:
    """False for `<PREFIX>_<64 hex>` -- a uniform hash accession carries no source."""
    if "_" in member:
        rest = member.rsplit("_", 1)[1]
        if len(rest) == 64 and _HEX64.match(rest):
            return False
    return True
def is_host_member(member: str) -> bool:
    if (member.startswith("sp|") or member.startswith("tr|")) and "_HUMAN" in member:
        return True
    return member.startswith(("GENCODE|", "REFSEQ|", "CHESS|", "UNIPROT|"))
def classify_members(members) -> str:
    has_human = has_micro = any_readable = False
    for m in members:
        if not m or not is_classifiable_accession(m):
            continue
        any_readable = True
        if is_host_member(m):
            has_human = True
        else:
            has_micro = True
        if has_human and has_micro:
            return "Mixed"
    if not any_readable:
        return "Unknown"
    return "Human" if has_human else "Microbial"
