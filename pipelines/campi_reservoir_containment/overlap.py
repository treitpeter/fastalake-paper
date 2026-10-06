#!/usr/bin/env python3
"""R8 Task B: presence of the eight SIHUMIx reference proteomes in the blind CAMPI reservoir.

Modes (one per Slurm step, all streaming, single CPU):
  hash      sha256 of the large inputs, written into inputs/manifest.json (pilot step 0)
  protein   one pass over the reservoir FASTA: exact-protein presence (sha256 of the
            I/L-collapsed sequence) for the 8 proteomes + 8 shuffled controls, and
            tryptic-context peptide presence (KR/P, 0 missed cleavages, 9-50, I/L)
  sources   one pass over the joint-headers TSV: every source tag of every lake entry
            that matched a reference protein -> gut / uniprot / both
  digest    write the B. longum tryptic peptides (real + shuffled) for the production extractor
  extractor-summary
            summarize the extractor hits TSV (substring presence, the production lookup)
  selftest  synthetic lake + synthetic reference; asserts the expected counts

Hash contract (canonical/rust/fasta_extractor_v2/src/bin/lake_builder.rs, hash v5.3):
uppercase A-Z only, everything else dropped, sha256. The lake deduplicates on that raw
hash (no I/L). This analysis keys on the same bytes with I->L applied (production
normalize_il), so one reference protein can match several lake entries (I/L variants).
"""
import argparse, hashlib, json, os, random, re, resource, sys, time
from collections import defaultdict

SOURCE_TAGS = ["UHGG", "GMGC", "GMSC", "UNIPROT_SP", "UNIPROT_TR"]
GUT = {"UHGG", "GMGC", "GMSC"}
UNI = {"UNIPROT_SP", "UNIPROT_TR"}
# Header prefix -> lake source tag (verified on the reservoir headers, see SCOPE.md)
PREFIX = [(b"MGYG", "UHGG"), (b"GMGC10.", "GMGC"), (b"GMSC10.", "GMSC"),
          (b"sp|", "UNIPROT_SP"), (b"tr|", "UNIPROT_TR")]
SRCBIT = {t: 1 << i for i, t in enumerate(SOURCE_TAGS)}
MIN_LEN, MAX_LEN = 9, 50
SPLIT = re.compile(rb"(?<=[KR])(?!P)")  # Sage KR/P rule, no initiator-M clipping

# byte translation tables: uppercase letters kept, everything else deleted
_RAW = bytearray(range(256))
for c in range(ord("a"), ord("z") + 1):
    _RAW[c] = c - 32
RAW_TABLE = bytes(_RAW)
IL_TABLE = bytearray(RAW_TABLE)
IL_TABLE[ord("I")] = ord("L")
IL_TABLE[ord("i")] = ord("L")
IL_TABLE = bytes(IL_TABLE)
DELETE = bytes(b for b in range(256) if not (65 <= b <= 90 or 97 <= b <= 122))


def raw_clean(s):
    return s.translate(RAW_TABLE, DELETE)


def il_clean(s):
    return s.translate(IL_TABLE, DELETE)


def source_of(protein_id):
    for p, t in PREFIX:
        if protein_id.startswith(p):
            return t
    return "UNKNOWN"


def classify(tags):
    tags = set(tags) - {"UNKNOWN"}
    if not tags:
        return "unresolved"
    g, u = bool(tags & GUT), bool(tags & UNI)
    return "both" if (g and u) else ("gut_only" if g else "uniprot_only")


def sha256_file(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def maxrss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def dump(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def write_tsv(path, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def fasta_records(path):
    with open(path, "rb") as f:
        header, parts = None, []
        for line in f:
            if line.startswith(b">"):
                if header is not None:
                    yield header, b"".join(parts)
                header, parts = line.rstrip(b"\r\n"), []
            else:
                parts.append(line)
        if header is not None:
            yield header, b"".join(parts)


# --------------------------------------------------------------------------- reference
def load_members(path):
    members = {}
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            members[r["campi_os_string"]] = (r["member"], int(r["reference_proteins"]))
    return members


def shuffled(seq, acc):
    rnd = random.Random(int(hashlib.sha256(acc).hexdigest()[:16], 16))
    l = list(seq)
    rnd.shuffle(l)
    return bytes(l)


def load_reference(ref_fasta, members_tsv):
    """Return proteins [(label, acc, raw_seq, il_seq)], proteome order, nonmember count."""
    members = load_members(members_tsv)
    os_re = re.compile(rb" OS=(.+?) (?:GN|PE|OX)=")
    proteins, nonmember, per_os = [], 0, defaultdict(int)
    for header, seq in fasta_records(ref_fasta):
        m = os_re.search(header)
        osname = m.group(1).decode() if m else None
        if osname is None or osname not in members:
            nonmember += 1
            continue
        label = members[osname][0]
        per_os[osname] += 1
        acc = header[1:].split()[0]
        raw = raw_clean(seq)
        proteins.append((label, acc.decode(), raw, il_clean(raw)))
        sh = shuffled(raw, acc)
        proteins.append((label + "_SHUFFLED", acc.decode(), sh, il_clean(sh)))
    for osname, (label, n) in members.items():
        if per_os[osname] != n:
            raise SystemExit(f"member count mismatch for {osname}: {per_os[osname]} != {n}")
    labels = [members[o][0] for o in members]
    return proteins, labels, nonmember


def tryptic(il_seq):
    return [p for p in SPLIT.split(il_seq) if MIN_LEN <= len(p) <= MAX_LEN]


# --------------------------------------------------------------------------- modes
def mode_hash(a):
    man = json.load(open(a.manifest))
    for key in a.keys:
        e = man["inputs"][key]
        st = os.stat(e["path"])
        if st.st_size != e["bytes"] or st.st_mtime_ns != e["mtime_ns"]:
            raise SystemExit(f"{key}: size/mtime changed since manifest was written")
        t0 = time.time()
        e["sha256"] = sha256_file(e["path"])
        e["sha256_computed_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        e["sha256_job_id"] = os.environ.get("SLURM_JOB_ID")
        e["sha256_seconds"] = round(time.time() - t0, 1)
        print(f"{key}\t{e['sha256']}\t{e['sha256_seconds']}s", flush=True)
    dump(a.manifest, man)


def mode_protein(a):
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    proteins, labels, nonmember = load_reference(a.reference, a.members)
    all_labels = labels + [l + "_SHUFFLED" for l in labels]
    lbit = {l: 1 << i for i, l in enumerate(all_labels)}
    # exact-protein index: I/L-collapsed sha256 digest -> list of protein indices
    prot_index = defaultdict(list)
    for i, (label, acc, raw, il) in enumerate(proteins):
        prot_index[hashlib.sha256(il).digest()].append(i)
    # peptide index: I/L peptide -> label mask
    pep_mask = {}
    for label, acc, raw, il in proteins:
        b = lbit[label]
        for p in tryptic(il):
            pep_mask[p] = pep_mask.get(p, 0) | b
    print(f"reference: {len(proteins)//2} member proteins (+ shuffled), {nonmember} nonmember "
          f"entries skipped, {len(prot_index)} distinct I/L keys, {len(pep_mask)} distinct "
          f"tryptic peptides {MIN_LEN}-{MAX_LEN} (real+shuffled)", flush=True)

    # ---- streaming pass over the reservoir
    limit = a.max_bytes if a.max_bytes else None
    slice_out = open(a.slice_out, "wb") if a.slice_out else None
    fh = hashlib.sha256()
    n_entries, n_bytes, src_counts = 0, 0, defaultdict(int)
    prot_hits = []  # (il_hex, raw_hex, source, protein_id, protein_index)
    pep_found = {}  # peptide -> source bit mask
    header, parts = None, []
    stop = False

    def process(header, seq):
        nonlocal n_entries
        n_entries += 1
        if n_entries % 25_000_000 == 0:
            print(f"  {n_entries} entries, {n_bytes/1e9:.1f} GB, {time.time()-t0:.0f}s, "
                  f"{len(prot_hits)} protein hits, {len(pep_found)} peptides found", flush=True)
        pid = header[1:].split(None, 1)[0]
        src = source_of(pid)
        src_counts[src] += 1
        il = il_clean(seq)
        d = hashlib.sha256(il).digest()
        if d in prot_index:
            raw = raw_clean(seq)
            raw_hex = hashlib.sha256(raw).hexdigest()
            for i in prot_index[d]:
                prot_hits.append((d.hex(), raw_hex, src, pid.decode(), i))
        bit = SRCBIT.get(src, 0)
        for p in SPLIT.split(il):
            if MIN_LEN <= len(p) <= MAX_LEN and p in pep_mask:
                pep_found[p] = pep_found.get(p, 0) | bit

    with open(a.reservoir, "rb", buffering=1 << 24) as f:
        for line in f:
            if line.startswith(b">"):
                if limit is not None and n_bytes >= limit:
                    stop = True
                    break
                if header is not None:
                    process(header, b"".join(parts))
                header, parts = line.rstrip(b"\r\n"), []
            else:
                parts.append(line)
            fh.update(line)
            n_bytes += len(line)
            if slice_out:
                slice_out.write(line)
        if header is not None and not stop:
            process(header, b"".join(parts))
        elif header is not None and stop:
            process(header, b"".join(parts))
    if slice_out:
        slice_out.close()
    elapsed = time.time() - t0
    full_pass = limit is None
    file_sha = fh.hexdigest() if full_pass else None
    print(f"reservoir pass: {n_entries} entries, {n_bytes} bytes, {elapsed:.0f}s, "
          f"maxrss {maxrss_gb():.2f} GB, full_pass={full_pass}", flush=True)
    if full_pass:
        man = json.load(open(a.manifest))
        exp = man["inputs"]["reservoir_fasta"].get("sha256")
        if exp and exp != file_sha:
            raise SystemExit(f"reservoir sha256 mismatch: {file_sha} != manifest {exp}")
        if a.expect_entries and n_entries != a.expect_entries:
            raise SystemExit(f"entry count {n_entries} != expected {a.expect_entries}")

    # ---- protein-level outputs
    write_tsv(os.path.join(a.out, "protein_hits.tsv"),
              ["il_sha256", "lake_raw_sha256", "kept_source", "lake_protein_id", "label", "accession",
               "exact_raw_match"],
              [(ih, rh, s, pid, proteins[i][0], proteins[i][1],
                int(rh == hashlib.sha256(proteins[i][2]).hexdigest())) for ih, rh, s, pid, i in prot_hits])
    per_protein = defaultdict(lambda: {"n_lake_entries": 0, "kept_sources": set(), "exact_raw": 0})
    for ih, rh, s, pid, i in prot_hits:
        e = per_protein[i]
        e["n_lake_entries"] += 1
        e["kept_sources"].add(s)
        e["exact_raw"] |= int(rh == hashlib.sha256(proteins[i][2]).hexdigest())
    detail = []
    for i, (label, acc, raw, il) in enumerate(proteins):
        e = per_protein.get(i)
        detail.append((label, acc, len(raw), int(e is not None), e["n_lake_entries"] if e else 0,
                       ",".join(sorted(e["kept_sources"])) if e else "", e["exact_raw"] if e else 0))
    write_tsv(os.path.join(a.out, "protein_presence_detail.tsv"),
              ["label", "accession", "length", "present_exact_il", "n_lake_entries", "kept_sources",
               "exact_raw_match"], detail)
    rows = []
    for label in all_labels:
        idx = [i for i, p in enumerate(proteins) if p[0] == label]
        pres = [i for i in idx if i in per_protein]
        uniq = len({proteins[i][3] for i in idx})
        kept = defaultdict(int)
        for i in pres:
            for s in per_protein[i]["kept_sources"]:
                kept[s] += 1
        rows.append((label, len(idx), uniq, len(pres), f"{len(pres)/len(idx):.4f}",
                     sum(per_protein[i]["exact_raw"] for i in pres),
                     len(pres) - sum(per_protein[i]["exact_raw"] for i in pres),
                     *[kept[t] for t in SOURCE_TAGS]))
    write_tsv(os.path.join(a.out, "protein_presence_by_proteome.tsv"),
              ["label", "n_proteins", "n_unique_il_sequences", "n_present_exact", "fraction",
               "n_present_exact_raw", "n_present_il_variant_only",
               *[f"n_kept_{t}" for t in SOURCE_TAGS]], rows)

    # ---- peptide-level (tryptic context) outputs
    pep_rows, len_rows = [], []
    for label in all_labels:
        b = lbit[label]
        peps = [p for p, m in pep_mask.items() if m & b]
        found = [p for p in peps if p in pep_found]
        cls = defaultdict(int)
        for p in found:
            cls[classify([t for t in SOURCE_TAGS if pep_found[p] & SRCBIT[t]])] += 1
        pep_rows.append((label, len(peps), len(found), f"{len(found)/max(1,len(peps)):.4f}",
                         cls["gut_only"], cls["uniprot_only"], cls["both"]))
        by_len = defaultdict(lambda: [0, 0])
        for p in peps:
            by_len[len(p)][0] += 1
            by_len[len(p)][1] += p in pep_found
        for L in sorted(by_len):
            len_rows.append((label, L, by_len[L][0], by_len[L][1], f"{by_len[L][1]/by_len[L][0]:.4f}"))
    write_tsv(os.path.join(a.out, "peptide_tryptic_presence_by_proteome.tsv"),
              ["label", "n_peptides_9_50_il_unique", "n_present_tryptic_context", "fraction",
               "n_gut_only_kept", "n_uniprot_only_kept", "n_both_kept"], pep_rows)
    write_tsv(os.path.join(a.out, "peptide_tryptic_presence_by_length.tsv"),
              ["label", "length", "n_peptides", "n_present", "fraction"], len_rows)
    with open(os.path.join(a.out, "peptide_tryptic_found.tsv"), "w") as f:
        f.write("peptide\tlabels\tkept_sources\n")
        for p, m in pep_found.items():
            f.write(p.decode() + "\t" + ",".join(l for l in all_labels if pep_mask[p] & lbit[l]) + "\t"
                    + ",".join(t for t in SOURCE_TAGS if m & SRCBIT[t]) + "\n")
    dump(os.path.join(a.out, "protein_pass.json"), {
        "reservoir": a.reservoir, "full_pass": full_pass, "bytes_read": n_bytes, "entries_read": n_entries,
        "reservoir_sha256": file_sha, "entries_by_source": dict(src_counts), "elapsed_s": round(elapsed, 1),
        "maxrss_gb": round(maxrss_gb(), 3), "n_protein_hits": len(prot_hits), "n_peptides_found": len(pep_found),
        "n_reference_member_proteins": len(proteins) // 2, "n_nonmember_reference_entries": nonmember,
        "min_len": MIN_LEN, "max_len": MAX_LEN, "missed_cleavages": 0, "rule": "KR/P, no initiator-M clipping",
        "job_id": os.environ.get("SLURM_JOB_ID")})


def mode_sources(a):
    t0 = time.time()
    hits = []
    with open(os.path.join(a.out, "protein_hits.tsv")) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            hits.append(dict(zip(cols, line.rstrip("\n").split("\t"))))
    want = {h["lake_raw_sha256"].encode() for h in hits}
    print(f"sources: {len(hits)} protein hits, {len(want)} distinct lake entries to resolve", flush=True)
    found = {}
    fh = hashlib.sha256()
    n_lines, n_bytes = 0, 0
    limit = a.max_bytes if a.max_bytes else None
    with open(a.joint_headers, "rb", buffering=1 << 24) as f:
        for line in f:
            fh.update(line)
            n_bytes += len(line)
            n_lines += 1
            if limit is not None and n_bytes > limit:
                break
            h = line[:64]
            if h in want:
                fields = line.rstrip(b"\r\n").split(b"\t")
                n_src = int(fields[1])
                rest = b"\t".join(fields[3:])
                tags = [piece.split(b"\t", 1)[0].decode() for piece in rest.split(b"|||") if piece]
                if len(tags) != n_src:
                    raise SystemExit(f"joint header parse error at line {n_lines}: {line[:200]!r}")
                found[h.decode()] = tags
    elapsed = time.time() - t0
    full_pass = limit is None
    file_sha = fh.hexdigest() if full_pass else None
    print(f"joint headers pass: {n_lines-1} rows, {n_bytes} bytes, {elapsed:.0f}s, "
          f"resolved {len(found)}/{len(want)}", flush=True)
    if full_pass:
        man = json.load(open(a.manifest))
        exp = man["inputs"]["joint_headers_tsv"].get("sha256")
        if exp and exp != file_sha:
            raise SystemExit(f"joint headers sha256 mismatch: {file_sha} != manifest {exp}")
        if len(found) != len(want):
            raise SystemExit(f"unresolved lake entries in full pass: {len(want)-len(found)}")
    # per reference protein: union of all source tags over all matching lake entries
    per_prot = defaultdict(lambda: {"tags": set(), "n_lake_entries": 0, "unresolved": 0, "label": None})
    for h in hits:
        k = (h["label"], h["accession"])
        e = per_prot[k]
        e["label"] = h["label"]
        e["n_lake_entries"] += 1
        tags = found.get(h["lake_raw_sha256"])
        if tags is None:
            e["unresolved"] += 1
            e["tags"].add(h["kept_source"])  # fallback: the kept (highest-priority) source
        else:
            e["tags"].update(tags)
    write_tsv(os.path.join(a.out, "protein_sources_detail.tsv"),
              ["label", "accession", "n_lake_entries", "n_unresolved", "all_sources", "class"],
              [(k[0], k[1], e["n_lake_entries"], e["unresolved"], ",".join(sorted(e["tags"])),
               classify(e["tags"]) + ("*" if e["unresolved"] else ""))
               for k, e in sorted(per_prot.items())])
    # by proteome, joined with the protein-level table
    base = {}
    with open(os.path.join(a.out, "protein_presence_by_proteome.tsv")) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            base[r["label"]] = r
    rows = []
    for label, r in base.items():
        cls = defaultdict(int)
        tagc = defaultdict(int)
        unresolved = 0
        for k, e in per_prot.items():
            if k[0] != label:
                continue
            cls[classify(e["tags"])] += 1
            unresolved += bool(e["unresolved"])
            for t in e["tags"]:
                tagc[t] += 1
        rows.append((label, r["n_proteins"], r["n_present_exact"], r["fraction"],
                     cls["gut_only"], cls["uniprot_only"], cls["both"], unresolved,
                     *[tagc[t] for t in SOURCE_TAGS]))
    write_tsv(os.path.join(a.out, "protein_presence_by_source.tsv"),
              ["label", "n_proteins", "n_present_exact", "fraction", "n_gut_only", "n_uniprot_only", "n_both",
               "n_with_unresolved_entries", *[f"n_in_{t}" for t in SOURCE_TAGS]], rows)
    dump(os.path.join(a.out, "sources_pass.json"), {
        "joint_headers": a.joint_headers, "full_pass": full_pass, "bytes_read": n_bytes, "rows_read": n_lines - 1,
        "joint_headers_sha256": file_sha, "resolved": len(found), "wanted": len(want),
        "elapsed_s": round(elapsed, 1), "maxrss_gb": round(maxrss_gb(), 3), "job_id": os.environ.get("SLURM_JOB_ID")})


def mode_digest(a):
    proteins, labels, _ = load_reference(a.reference, a.members)
    os.makedirs(a.out, exist_ok=True)
    real, shuf = set(), set()
    for label, acc, raw, il in proteins:
        if label == a.label:
            real.update(tryptic(il))
        elif label == a.label + "_SHUFFLED":
            shuf.update(tryptic(il))
    peps = sorted(real | shuf)
    with open(os.path.join(a.out, "extractor_queries.txt"), "w") as f:
        for p in peps:
            f.write(p.decode() + "\n")
    write_tsv(os.path.join(a.out, "extractor_query_labels.tsv"), ["peptide", "real", "shuffled"],
              [(p.decode(), int(p in real), int(p in shuf)) for p in peps])
    print(f"digest {a.label}: {len(real)} real, {len(shuf)} shuffled, {len(real & shuf)} shared, "
          f"{len(peps)} query lines", flush=True)


def mode_extractor_summary(a):
    labels = {}
    with open(os.path.join(a.out, "extractor_query_labels.tsv")) as f:
        f.readline()
        for line in f:
            p, r, s = line.rstrip("\n").split("\t")
            labels[p.encode()] = (int(r), int(s))
    found = {}
    n_rows = 0
    with open(a.hits, "rb", buffering=1 << 24) as f:
        f.readline()
        for line in f:
            n_rows += 1
            pep, pid, _ = line.rstrip(b"\r\n").split(b"\t")
            found[pep] = found.get(pep, 0) | SRCBIT.get(source_of(pid), 0)
    unknown = [p for p in found if p not in labels]
    if unknown:
        raise SystemExit(f"{len(unknown)} hit peptides not in the query set, e.g. {unknown[:3]}")
    rows, len_rows = [], []
    for name, col in (("real", 0), ("shuffled", 1)):
        peps = [p for p, l in labels.items() if l[col]]
        fnd = [p for p in peps if p in found]
        cls = defaultdict(int)
        for p in fnd:
            cls[classify([t for t in SOURCE_TAGS if found[p] & SRCBIT[t]])] += 1
        rows.append((a.label + ("" if col == 0 else "_SHUFFLED"), len(peps), len(fnd),
                     f"{len(fnd)/max(1,len(peps)):.4f}", cls["gut_only"], cls["uniprot_only"], cls["both"]))
        by_len = defaultdict(lambda: [0, 0])
        for p in peps:
            by_len[len(p)][0] += 1
            by_len[len(p)][1] += p in found
        for L in sorted(by_len):
            len_rows.append((a.label + ("" if col == 0 else "_SHUFFLED"), L, by_len[L][0], by_len[L][1],
                             f"{by_len[L][1]/by_len[L][0]:.4f}"))
    write_tsv(os.path.join(a.out, "peptide_substring_presence_extractor.tsv"),
              ["label", "n_peptides_9_50_il_unique", "n_present_substring", "fraction",
               "n_gut_only_kept", "n_uniprot_only_kept", "n_both_kept"], rows)
    write_tsv(os.path.join(a.out, "peptide_substring_presence_by_length.tsv"),
              ["label", "length", "n_peptides", "n_present", "fraction"], len_rows)
    # calibration against the tryptic-context measure from the protein pass
    tryp = {}
    with open(os.path.join(a.out, "peptide_tryptic_found.tsv")) as f:
        f.readline()
        for line in f:
            p, ls, _ = line.rstrip("\n").split("\t")
            tryp[p.encode()] = ls.split(",")
    real = [p for p, l in labels.items() if l[0]]
    both = sum(1 for p in real if p in found and p in tryp)
    only_sub = sum(1 for p in real if p in found and p not in tryp)
    only_tryp = sum(1 for p in real if p not in found and p in tryp)
    dump(os.path.join(a.out, "extractor_summary.json"), {
        "label": a.label, "hits_rows": n_rows, "n_real_peptides": len(real),
        "present_substring_and_tryptic": both, "present_substring_only": only_sub,
        "present_tryptic_only_UNEXPECTED": only_tryp, "hits_file": a.hits,
        "job_id": os.environ.get("SLURM_JOB_ID"), "maxrss_gb": round(maxrss_gb(), 3)})
    if only_tryp:
        print(f"WARNING: {only_tryp} peptides present in tryptic context but not found by the extractor", flush=True)
    print(f"extractor summary: {n_rows} hit rows; real {len(real)} peptides, substring-present "
          f"{both+only_sub}, tryptic-context-present {both+only_tryp}", flush=True)


def mode_selftest(a):
    """Synthetic negative/positive control. Exits non-zero on any deviation."""
    d = a.out
    os.makedirs(d, exist_ok=True)
    A = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGEVQAGVWPAAVRESVPSLL"
    B = b"MSKIVKIIGREIIDSRGNPTVEAEVHLEGGFVGMAAAPSGASTGSREALELRDGDKSRFLGKGVTKAVAAVNGPIAQALIGKDAKDQAGIDKIMIDLDGTENKSKFGANAILAVSLANAKAAAAAKGMPLYEHIAELNGTPGKYSMPVPMMNIINGGEHADNNVDIQEFMIQPVGAKTVKEAIRMGSEVFHHLAKVLKAKGMNTAVGDEGGYAPNLGSNAEALAVIAEAVKAAGYELGKDITLAMDCAASEFYKDGKYVLAGEGNKAFTSEEFTHFLEELTKQYPIVSIEDGLDESDWDGFAYQTKVLGDKIQLVGDDLFVTNTKILKEGIEKGIANSILIKFNQIGSLTETLAAIKMAKDAGYTAVISHRSGETEDATIADLAVGTAAGQIKTGSMSRSDRVAKYNQLIRIEEALGEKAPYNGRKEIKGQA"
    C = b"MAKLTKRMRVIREKVDATKQYDINEAIALLKELATAKFVESVDVAVNLGIDARKSDQNVRGATVLPHGTGRSVRVAVFTQGANAEAAKAAGAELVGMEDLADQIKKGEMNFDVVIASPDAMRVVGQLGQVLGPRGLMPNPKVGTVTPNVAEAVKNAKAGQVRYRNDKNGIIHTTIGKVDFDADKLKENLEALLVALKKAKPTQAKGVYIKKVSISTTMGAGVAVDQAGLSASVN"
    ref = os.path.join(d, "ref.faa")
    with open(ref, "w") as f:
        f.write(">tr|T1|T1_TEST protein one OS=Anaerostipes caccae DSM 14662 GN=x PE=3 SV=1\n" + A.decode() + "\n")
        f.write(">tr|T2|T2_TEST protein two OS=Anaerostipes caccae DSM 14662 GN=y PE=3 SV=1\n" + B.decode() + "\n")
        f.write(">tr|T3|T3_TEST protein three OS=Anaerostipes caccae DSM 14662 GN=z PE=3 SV=1\n" + C.decode() + "\n")
        f.write(">sp|ALBU_BOVIN|\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE\n")
    members = os.path.join(d, "members.tsv")
    with open(members, "w") as f:
        f.write("member\tcampi_os_string\treference_proteins\taliases\n")
        f.write("Anaerostipes caccae\tAnaerostipes caccae DSM 14662\t3\t\n")
    # lake: A exact under a UHGG header, A again under TrEMBL header (same raw seq -> dedup in a
    # real lake, here kept to test union), B as I->L variant under GMSC, C absent, one shuffled
    # decoy that must NOT count, plus an unrelated protein.
    B_il = B.replace(b"I", b"L")
    lake = os.path.join(d, "lake.fasta")
    with open(lake, "w") as f:
        f.write(">MGYG000000001_00001 hypothetical protein\n" + A.decode() + "\n")
        f.write(">GMSC10.100AA.000_000_001\n" + B_il.decode() + "\n")
        f.write(">tr|X1|X1_TEST unrelated OS=Foo bar OX=1 PE=4 SV=1\nMSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA\n")
    joint = os.path.join(d, "joint.tsv")
    hA = hashlib.sha256(A).hexdigest()
    hB = hashlib.sha256(B_il).hexdigest()
    with open(joint, "w") as f:
        f.write("sha256_hash\tn_sources\tkept_header\tall_source_headers\n")
        f.write(f"{hA}\t2\tMGYG000000001_00001 hypothetical protein\tUHGG\tMGYG000000001_00001 hypothetical protein|||UNIPROT_TR\ttr|T1|T1_TEST protein one OS=Anaerostipes caccae DSM 14662 GN=x PE=3 SV=1\n")
        f.write(f"{hB}\t1\tGMSC10.100AA.000_000_001\tGMSC\tGMSC10.100AA.000_000_001\n")
    man = os.path.join(d, "manifest.json")
    dump(man, {"inputs": {"reservoir_fasta": {"path": lake}, "joint_headers_tsv": {"path": joint}}})
    ns = argparse.Namespace(reference=ref, members=members, reservoir=lake, manifest=man, out=d,
                            max_bytes=0, slice_out=None, expect_entries=3, joint_headers=joint)
    mode_protein(ns)
    mode_sources(ns)
    rows = {}
    with open(os.path.join(d, "protein_presence_by_source.tsv")) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            rows[r["label"]] = r
    real, shuf = rows["Anaerostipes caccae"], rows["Anaerostipes caccae_SHUFFLED"]
    checks = [
        ("n_proteins", real["n_proteins"], "3"), ("present", real["n_present_exact"], "2"),
        ("both (A)", real["n_both"], "1"), ("gut_only (B via I/L)", real["n_gut_only"], "1"),
        ("uniprot_only", real["n_uniprot_only"], "0"),
        ("shuffled present", shuf["n_present_exact"], "0"),
    ]
    det = {}
    with open(os.path.join(d, "protein_presence_detail.tsv")) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            det[(r["label"], r["accession"])] = r
    checks += [("B is I/L variant only", det[("Anaerostipes caccae", "tr|T2|T2_TEST")]["exact_raw_match"], "0"),
               ("A exact raw", det[("Anaerostipes caccae", "tr|T1|T1_TEST")]["exact_raw_match"], "1"),
               ("C absent", det[("Anaerostipes caccae", "tr|T3|T3_TEST")]["present_exact_il"], "0")]
    pep = {}
    with open(os.path.join(d, "peptide_tryptic_presence_by_proteome.tsv")) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            pep[r["label"]] = r
    nA = len(set(tryptic(A.replace(b"I", b"L"))) | set(tryptic(B_il)))
    checks += [("tryptic peptides of A+B all present, none of C", pep["Anaerostipes caccae"]["n_present_tryptic_context"], str(nA)),
               ("shuffled peptides absent", pep["Anaerostipes caccae_SHUFFLED"]["n_present_tryptic_context"], "0")]
    ok = True
    for name, got, exp in checks:
        flag = "PASS" if got == exp else "FAIL"
        ok &= got == exp
        print(f"  selftest {flag}: {name}: got {got!r} expected {exp!r}")
    # digest + extractor summary path with a fake hits file
    ns2 = argparse.Namespace(reference=ref, members=members, out=d, label="Anaerostipes caccae")
    mode_digest(ns2)
    fake_hits = os.path.join(d, "hits.tsv")
    with open(fake_hits, "w") as f:
        f.write("peptide\tprotein_id\tposition\n")
        for p in tryptic(A.replace(b"I", b"L")):
            f.write(f"{p.decode()}\tMGYG000000001_00001\t0\n")
    ns3 = argparse.Namespace(out=d, hits=fake_hits, label="Anaerostipes caccae")
    mode_extractor_summary(ns3)
    print("SELFTEST", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["hash", "protein", "sources", "digest", "extractor-summary", "selftest"])
    ap.add_argument("--manifest")
    ap.add_argument("--keys", nargs="*", default=["reservoir_fasta", "joint_headers_tsv"])
    ap.add_argument("--reference")
    ap.add_argument("--members")
    ap.add_argument("--reservoir")
    ap.add_argument("--joint-headers")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-bytes", type=int, default=0, help="stop after this many bytes (pilot slice)")
    ap.add_argument("--slice-out", help="write the consumed part of the reservoir here (pilot)")
    ap.add_argument("--expect-entries", type=int, default=0)
    ap.add_argument("--label", default="Bifidobacterium longum")
    ap.add_argument("--hits")
    a = ap.parse_args()
    {"hash": mode_hash, "protein": mode_protein, "sources": mode_sources, "digest": mode_digest,
     "extractor-summary": mode_extractor_summary, "selftest": mode_selftest}[a.mode](a)


if __name__ == "__main__":
    main()
