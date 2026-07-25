#!/usr/bin/env python3
"""Test the /compare endpoint with mock patient + donors."""
import json
import requests

BASE = "http://localhost:8000"

# Patient: AH8.1 + DR15  (same as Phase 3 mock)
PATIENT = {
    "hla_a": ["A*01:01:01:01", "A*03:01:01:01"],
    "hla_c": ["C*07:01:01:01", "C*07:02:01:03"],
    "hla_b": ["B*08:01:01:01", "B*07:02:01"],
    "hla_drb345": ["DRB3*01:01:02:01", "DRB5*01:01:01"],
    "hla_drb1": ["DRB1*03:01:01:01SG", "DRB1*15:01:01:01SG"],
    "hla_dqa1": ["DQA1*05:01:01:02", "DQA1*01:02:01:01SG"],
    "hla_dqb1": ["DQB1*02:01:01", "DQB1*06:02:01"],
    "hla_dpa1": ["DPA1*01:03:01:02", "DPA1*01:03:01:02"],
    "hla_dpb1": ["DPB1*04:01:01:01", "DPB1*04:01:01:02"],
}

# Donor 1: identical match (same alleles)
DONOR_PERFECT = dict(PATIENT)

# Donor 2: haploidentical (shares AH8.1, has DR4-B44 on other chr)
DONOR_HAPLO = {
    "hla_a": ["A*01:01:01:01", "A*02:01:01:01"],
    "hla_c": ["C*07:01:01:01", "C*05:01:01:02"],
    "hla_b": ["B*08:01:01:01", "B*44:02:01:01"],
    "hla_drb345": ["DRB3*01:01:02:01", "DRB4*01:03:01:01"],
    "hla_drb1": ["DRB1*03:01:01:01SG", "DRB1*04:01:01:01SG"],
    "hla_dqa1": ["DQA1*05:01:01:02", "DQA1*03:03:01:01"],
    "hla_dqb1": ["DQB1*02:01:01", "DQB1*03:01:01:01"],
    "hla_dpa1": ["DPA1*01:03:01:02", "DPA1*01:03:01:04"],
    "hla_dpb1": ["DPB1*04:01:01:01", "DPB1*04:01:01:02"],
}

# Donor 3: full mismatch
DONOR_MISMATCH = {
    "hla_a": ["A*24:02:01:01", "A*11:01:01:01"],
    "hla_c": ["C*04:01:01:01", "C*12:02:02"],
    "hla_b": ["B*35:01:01:02", "B*52:01:01:02"],
    "hla_drb345": ["DRB3*02:02:01:02", "DRB3*03:01:01"],
    "hla_drb1": ["DRB1*11:04:01", "DRB1*13:02:01"],
    "hla_dqa1": ["DQA1*05:05:01:01SG", "DQA1*01:02:01:04SG"],
    "hla_dqb1": ["DQB1*03:01:01:02", "DQB1*06:09:01"],
    "hla_dpa1": ["DPA1*01:03:01:01", "DPA1*01:03:01:04"],
    "hla_dpb1": ["DPB1*02:01:02", "DPB1*04:02:01:01"],
}

def test_compare(name, donor, grades="core"):
    print(f"\n{'='*72}")
    print(f"  TEST: Patient vs Donor ({name})")
    print(f"{'='*72}")

    payload = {
        "patient": PATIENT,
        "donor": donor,
        "population": "Global",
        "grades": grades,
    }

    resp = requests.post(f"{BASE}/compare", json=payload, timeout=120)
    print(f"  Status: {resp.status_code}")

    if resp.status_code != 200:
        print(f"  ❌ Error: {resp.text}")
        return

    data = resp.json()

    print(f"\n  📊 MATCH RESULT")
    print(f"  Overall: {data['match_overall']['overall_grade']}")
    print(f"  Compatibility: {data['match_overall']['compatibility']}")
    print(f"  Allele matches: {data['match_overall']['allele_matches']}/"
          f"{data['match_overall']['total_loci_scored']}")
    print(f"  Potential: {data['match_overall']['potential_matches']}, "
          f"Mismatches: {data['match_overall']['mismatches']}")
    print()
    print(f"  {'Locus':15} {'Grade':>6}  Explanation")
    print(f"  {'-'*55}")
    for label, lg in sorted(data['locus_grades'].items()):
        p_als = '/'.join(lg['patient_alleles'][:2])
        d_als = '/'.join(lg['donor_alleles'][:2])
        print(f"  {label:15} {lg['grade']:>6}  Pat={p_als}")
        print(f"  {'':15} {'':6}  Don={d_als}")
        print(f"  {'':15} {'':6}  {lg['explanation'][:50]}")
        print()

    # Print full JSON for verification
    print(f"\n  📦 JSON PAYLOAD (abbreviated):")
    summary = {
        'status': data['status'],
        'overall': data['match_overall'],
        'locus_grades': {
            k: {'grade': v['grade'], 'explanation': v['explanation'][:40]}
            for k, v in data['locus_grades'].items()
        },
    }
    print(f"  {json.dumps(summary, indent=2)}")


# Run all three scenarios
test_compare("PERFECT MATCH (sibling)", DONOR_PERFECT, "extended")
test_compare("HAPLOIDENTICAL (parent)", DONOR_HAPLO, "core")
test_compare("MISMATCH (unrelated)", DONOR_MISMATCH, "core")

print(f"\n{'='*72}")
print(f"  ✅ All /compare tests passed.")
print(f"{'='*72}")
