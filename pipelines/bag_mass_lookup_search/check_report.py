"""Check report arithmetic, references, embedded PDF fonts and page boundaries."""
from pathlib import Path
import hashlib
import json
import fitz
import pandas as pd

HERE = Path(__file__).resolve().parent
G = HERE.parent


def main():
    contrasts = pd.read_csv(G / "results/CONTRASTS.tsv", sep="\t")
    mechanisms = pd.read_csv(G / "results/PEPTIDE_MECHANISMS.tsv", sep="\t")
    for arm, expected in (("bags", (782,452,330)), ("mass", (4222,5161,-939)), ("combined", (4226,5181,-955))):
        d = contrasts[contrasts.contrast == arm + "_vs_exact"]
        assert tuple(d[["gained","lost","net"]].sum()) == expected
    assert (contrasts[contrasts.contrast == "combined_vs_bags"].net < 0).all()
    cm = contrasts[contrasts.contrast == "combined_vs_mass"]
    assert tuple(cm[["gained","lost","net"]].sum()) == (29,45,-16)
    mass = mechanisms[mechanisms.contrast == "mass_vs_exact"]
    removed = mass[(mass.change == "lost") & (mass.category == "sequence_removed")]
    new = mass[(mass.change == "gained") & (mass.category == "sequence_added")]
    retained = mass[(mass.change == "lost") & (mass.category == "sequence_retained_above_q_cutoff")]
    assert len(removed) == 1992 and removed.length.median() == 9 and (removed.own_min_q <= .001).sum() == 1526
    assert len(new) == 3342 and new.length.median() == 11 and (new.own_min_q <= .001).sum() == 2193
    assert len(retained) == 2142 and retained.length.median() == 8 and retained.shared_flanking_context.sum() == 2140
    assert abs(retained.own_min_q.median()-.00876) < .00001
    report_inputs = json.loads((HERE / "REPORT_INPUTS.json").read_text())
    for path, expected in report_inputs["sources"].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected["sha256"], path
    pdf = HERE / "Bags_Flips_Verification_Peter_Philipp_2026-09-24.pdf"
    diagnostics = []
    with fitz.open(str(pdf)) as doc:
        assert len(doc) == 17
        for i, page in enumerate(doc):
            spans = [s for block in page.get_text("dict")["blocks"] if "lines" in block for line in block["lines"] for s in line["spans"]]
            assert len(page.get_text()) > 500, i+1
            bad = [s["text"] for s in spans if s["bbox"][0] < 35 or s["bbox"][2] > page.rect.width-35 or s["bbox"][3] > page.rect.height-12]
            assert not bad, (i+1, bad)
            fonts = page.get_fonts()
            assert any("hne" in f[3] for f in fonts), (i+1, fonts)
            diagnostics.append(dict(page=i+1, text_characters=len(page.get_text()), text_bounds_pass=True))
        text = "\n".join(p.get_text() for p in doc)
        assert "116/127" in text and "0 / 24" in text and "1.84" in text
    result = dict(status="PASS_ARITHMETIC_SOURCE_HASHES_PDF_STRUCTURE", pages=diagnostics,
                  pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),
                  caveat="Page-boundary checks supplement, rather than replace, visual review of charts and typography.")
    (HERE / "REPORT_VALIDATION.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(status=result["status"], pages=len(diagnostics))))


if __name__ == "__main__":
    main()
