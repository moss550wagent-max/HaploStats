#!/usr/bin/env python3
"""
HaploStats — Clinical Match Grade Calculator (HapLogic Standards)
Phase 9: Convert imputed haplotype pairs into NMDP-style match grades.

Grades per locus:
  A (Allele match)  — Both patient alleles = both donor alleles (high-res)
  P (Potential)     — Only 1 allele matches, or posterior confidence < 0.5
  M (Mismatch)      — Neither allele matches

Overall grade:
  n/10, n/8, etc. — count of A-grade loci across the panel
"""

import json
import sys
from pathlib import Path
from collections import defaultdict

# Core grades
ALLELE = "A"   # Allele match
POTENTIAL = "P"  # Potential match
MISMATCH = "M"   # Mismatch

# Loci used for clinical matching (NMDP standard 8/8, 10/10)
CORE_MATCH_LOCL = ['hla_a', 'hla_c', 'hla_b', 'hla_drb1', 'hla_dqb1']
EXTENDED_MATCH_LOCL = ['hla_a', 'hla_c', 'hla_b', 'hla_drb1', 'hla_dqb1',
                       'hla_dpb1']
FULL_MATCH_LOCL = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
                   'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']

LOCUS_LABELS = {
    'hla_a': 'HLA-A', 'hla_c': 'HLA-C', 'hla_b': 'HLA-B',
    'hla_drb345': 'HLA-DRB345', 'hla_drb1': 'HLA-DRB1',
    'hla_dqa1': 'HLA-DQA1', 'hla_dqb1': 'HLA-DQB1',
    'hla_dpa1': 'HLA-DPA1', 'hla_dpb1': 'HLA-DPB1',
}


# ── Allele matching (reused from engine) ───────────────────────────

def _norm(a: str) -> str:
    return a.strip().replace(' ', '') if a else ''

def _allele_match(pat, ref) -> bool:
    pa, ra = _norm(str(pat)), _norm(str(ref))
    if pa == ra:
        return True
    if ra.startswith('['):
        try:
            return any(_allele_match(pa, c) for c in json.loads(ra))
        except Exception:
            return False
    if pa.startswith('['):
        try:
            return any(_allele_match(c, ra) for c in json.loads(pa))
        except Exception:
            return False
    pp, rp = pa.split(':'), ra.split(':')
    for a, b in zip(pp, rp):
        if a != b:
            return False
    return True


def _parse_haplotype_label(label: str) -> dict:
    """Parse a formatted haplotype label back into locus→allele dict."""
    result = {}
    parts = label.split('|')
    for part in parts:
        part = part.strip()
        if '=' not in part:
            continue
        key, val = part.split('=', 1)
        key = key.strip().lower().replace('-', '_')
        display_to_internal = {
            'hla_a': 'hla_a', 'hla_c': 'hla_c', 'hla_b': 'hla_b',
            'drb345': 'hla_drb345', 'drb1': 'hla_drb1',
            'dqa1': 'hla_dqa1', 'dqb1': 'hla_dqb1',
            'dpa1': 'hla_dpa1', 'dpb1': 'hla_dpb1',
        }
        internal_key = display_to_internal.get(key)
        if internal_key:
            result[internal_key] = val.strip()
    return result


# ════════════════════════════════════════════════════════════════════
# Match Grader
# ════════════════════════════════════════════════════════════════════

class MatchGrader:
    """
    Clinically-oriented match grade calculator.

    Takes imputed patient and donor profiles (top-ranked haplotype pairs)
    and assigns NMDP-standard match grades for each locus.
    """

    def __init__(self):
        self.loci_scored = FULL_MATCH_LOCL

    def grade_locus(self, patient_alleles: list, donor_alleles: list,
                    patient_posterior: float = 1.0) -> dict:
        """
        Assign a match grade for a single locus.

        Parameters
        ----------
        patient_alleles : [p_allele1, p_allele2]
        donor_alleles   : [d_allele1, d_allele2]
        patient_posterior : posterior probability of the patient's top pair

        Returns dict with 'grade', 'patient', 'donor', 'explanation'.
        """
        p1, p2 = patient_alleles if len(patient_alleles) >= 2 else [patient_alleles[0], patient_alleles[0]]
        d1, d2 = donor_alleles if len(donor_alleles) >= 2 else [donor_alleles[0], donor_alleles[0]]

        p1, p2 = str(p1), str(p2)
        d1, d2 = str(d1), str(d2)

        # Check matches (order-agnostic)
        # Does d1 match either p1 or p2?
        d1_match_p = _allele_match(d1, p1) or _allele_match(d1, p2)
        d2_match_p = _allele_match(d2, p1) or _allele_match(d2, p2)

        # Both donor alleles match patient?
        both_match = d1_match_p and d2_match_p

        # Are they exact (all 4 alleles must match pairwise)?
        # For A grade: patient set == donor set
        set_p = {_norm(p1), _norm(p2)}
        set_d = {_norm(d1), _norm(d2)}
        exact_match = (set_p == set_d)

        if exact_match and patient_posterior >= 0.5:
            grade = ALLELE
            explanation = "Allele match — all 4 alleles match at high resolution"
        elif both_match:
            grade = ALLELE
            explanation = "Allele match — donor pair covers both patient alleles"
        elif d1_match_p or d2_match_p:
            # Only 1 of 2 matches
            if patient_posterior < 0.5:
                grade = POTENTIAL
                explanation = f"Potential match (posterior={patient_posterior:.3f})"
            else:
                grade = POTENTIAL
                explanation = "Only 1 of 2 alleles can be matched"
        else:
            grade = MISMATCH
            explanation = f"No allele match: {d1}, {d2} vs {p1}, {p2}"

        return {
            'grade': grade,
            'patient_alleles': [p1, p2],
            'donor_alleles': [d1, d2],
            'patient_posterior': round(patient_posterior, 4),
            'exact_match': exact_match,
            'explanation': explanation,
        }

    def compare_profiles(self, patient_top: dict, donor_top: dict,
                         loci: list = None) -> dict:
        """
        Full donor-recipient comparison.

        Parameters
        ----------
        patient_top : dict from top_3[0] — must have haplotype_1, haplotype_2, posterior
        donor_top   : dict from top_3[0] — same structure

        Returns dict with overall grade + locus-by-locus breakdown.
        """
        if loci is None:
            loci = CORE_MATCH_LOCL

        # Parse the haplotype labels
        p_h1 = _parse_haplotype_label(patient_top.get('haplotype_1', ''))
        p_h2 = _parse_haplotype_label(patient_top.get('haplotype_2', ''))
        d_h1 = _parse_haplotype_label(donor_top.get('haplotype_1', ''))
        d_h2 = _parse_haplotype_label(donor_top.get('haplotype_2', ''))

        p_posterior = patient_top.get('posterior', 1.0)
        d_posterior = donor_top.get('posterior', 1.0)

        locus_results = {}
        allele_count = 0
        potential_count = 0
        mismatch_count = 0
        total_scorable = 0

        for loc in loci:
            label = LOCUS_LABELS.get(loc, loc)

            # Patient alleles: combine H1 + H2, deduplicate
            p_a1 = p_h1.get(loc, '')
            p_a2 = p_h2.get(loc, '')
            p_set = []
            seen = set()
            for a in [p_a1, p_a2]:
                n = _norm(a)
                if a and n not in seen:
                    p_set.append(a)
                    seen.add(n)

            # Donor alleles
            d_a1 = d_h1.get(loc, '')
            d_a2 = d_h2.get(loc, '')
            d_set = []
            seen = set()
            for a in [d_a1, d_a2]:
                n = _norm(a)
                if a and n not in seen:
                    d_set.append(a)
                    seen.add(n)

            if not p_set or not d_set:
                locus_results[label] = {
                    'grade': '?',
                    'patient_alleles': p_set,
                    'donor_alleles': d_set,
                    'explanation': 'Insufficient data — alleles missing',
                }
                continue

            result = self.grade_locus(p_set, d_set, min(p_posterior, d_posterior))
            locus_results[label] = result
            total_scorable += 1

            if result['grade'] == ALLELE:
                allele_count += 1
            elif result['grade'] == POTENTIAL:
                potential_count += 1
            else:
                mismatch_count += 1

        # Compute overall grade
        overall_summary = {
            'total_loci_scored': total_scorable,
            'allele_matches': allele_count,
            'potential_matches': potential_count,
            'mismatches': mismatch_count,
        }

        overall_grade = f"{allele_count}/{total_scorable}"
        if potential_count > 0:
            overall_grade += f" (+{potential_count} potential)"
        if mismatch_count > 0:
            overall_grade += f" ({mismatch_count} mismatch)"
        overall_summary['overall_grade'] = overall_grade

        if mismatch_count == 0 and potential_count == 0:
            overall_summary['compatibility'] = "FULL MATCH"
        elif mismatch_count == 0 and potential_count > 0:
            overall_summary['compatibility'] = "POTENTIAL MATCH"
        elif mismatch_count <= total_scorable * 0.3:
            overall_summary['compatibility'] = "PARTIAL MATCH"
        else:
            overall_summary['compatibility'] = "MISMATCH"

        return {
            'overall': overall_summary,
            'loci': locus_results,
            'patient_top_pair': {
                'haplotype_1': patient_top.get('haplotype_1', ''),
                'haplotype_2': patient_top.get('haplotype_2', ''),
                'posterior': patient_top.get('posterior', 0),
            },
            'donor_top_pair': {
                'haplotype_1': donor_top.get('haplotype_1', ''),
                'haplotype_2': donor_top.get('haplotype_2', ''),
                'posterior': donor_top.get('posterior', 0),
            },
        }


# ── Quick test (standalone) ────────────────────────────────────────

def _make_mock_top(allele_strings: dict, posterior=1.0) -> dict:
    """Build a mock top_3[0] dict for testing."""
    loci = ['hla_a', 'hla_c', 'hla_b', 'hla_drb345', 'hla_drb1',
            'hla_dqa1', 'hla_dqb1', 'hla_dpa1', 'hla_dpb1']
    labels = ['HLA-A', 'HLA-C', 'HLA-B', 'DRB345', 'DRB1',
              'DQA1', 'DQB1', 'DPA1', 'DPB1']

    def make_label(name_parts):
        return ' | '.join(
            f"{lbl}={name_parts.get(loc, '?')}"
            for loc, lbl in zip(loci, labels)
        )

    h1_parts = {}
    h2_parts = {}
    for loc in loci:
        vals = allele_strings.get(loc, ['?', '?'])
        h1_parts[loc] = vals[0] if len(vals) > 0 else '?'
        h2_parts[loc] = vals[1] if len(vals) > 1 else vals[0]

    return {
        'haplotype_1': make_label(h1_parts),
        'haplotype_2': make_label(h2_parts),
        'posterior': posterior,
    }


def demo():
    """Demo with mock patient and donor profiles."""
    grader = MatchGrader()

    # Mock Patient: AH8.1 / B7-DR15 (common European)
    patient = _make_mock_top({
        'hla_a':    ['A*01:01:01:01', 'A*03:01:01:01'],
        'hla_c':    ['C*07:01:01:01', 'C*07:02:01:03'],
        'hla_b':    ['B*08:01:01:01', 'B*07:02:01'],
        'hla_drb345': ['DRB3*01:01:02:01', 'DRB5*01:01:01'],
        'hla_drb1': ['DRB1*03:01:01:01SG', 'DRB1*15:01:01:01SG'],
        'hla_dqa1': ['DQA1*05:01:01:02', 'DQA1*01:02:01:01SG'],
        'hla_dqb1': ['DQB1*02:01:01', 'DQB1*06:02:01'],
        'hla_dpa1': ['DPA1*01:03:01:02', 'DPA1*01:03:01:02'],
        'hla_dpb1': ['DPB1*04:01:01:01', 'DPB1*04:01:01:02'],
    }, posterior=0.995)

    # Mock Donor: same as patient (perfect match)
    donor_perfect = _make_mock_top({
        'hla_a':    ['A*01:01:01:01', 'A*03:01:01:01'],
        'hla_c':    ['C*07:01:01:01', 'C*07:02:01:03'],
        'hla_b':    ['B*08:01:01:01', 'B*07:02:01'],
        'hla_drb345': ['DRB3*01:01:02:01', 'DRB5*01:01:01'],
        'hla_drb1': ['DRB1*03:01:01:01SG', 'DRB1*15:01:01:01SG'],
        'hla_dqa1': ['DQA1*05:01:01:02', 'DQA1*01:02:01:01SG'],
        'hla_dqb1': ['DQB1*02:01:01', 'DQB1*06:02:01'],
        'hla_dpa1': ['DPA1*01:03:01:02', 'DPA1*01:03:01:02'],
        'hla_dpb1': ['DPB1*04:01:01:01', 'DPB1*04:01:01:02'],
    }, posterior=0.995)

    # Mock Donor 2: haploidentical (shares one haplotype, mismatched at other)
    donor_half = _make_mock_top({
        'hla_a':    ['A*01:01:01:01', 'A*02:01:01:01'],
        'hla_c':    ['C*07:01:01:01', 'C*05:01:01:02'],
        'hla_b':    ['B*08:01:01:01', 'B*44:02:01:01'],
        'hla_drb345': ['DRB3*01:01:02:01', 'DRB4*01:03:01:01'],
        'hla_drb1': ['DRB1*03:01:01:01SG', 'DRB1*04:01:01:01SG'],
        'hla_dqa1': ['DQA1*05:01:01:02', 'DQA1*03:03:01:01'],
        'hla_dqb1': ['DQB1*02:01:01', 'DQB1*03:01:01:01'],
        'hla_dpa1': ['DPA1*01:03:01:02', 'DPA1*01:03:01:04'],
        'hla_dpb1': ['DPB1*04:01:01:01', 'DPB1*04:01:01:02'],
    }, posterior=0.850)

    # Mock Donor 3: full mismatch
    donor_mismatch = _make_mock_top({
        'hla_a':    ['A*24:02:01:01', 'A*11:01:01:01'],
        'hla_c':    ['C*04:01:01:01', 'C*12:02:02'],
        'hla_b':    ['B*35:01:01:02', 'B*52:01:01:02'],
        'hla_drb345': ['DRB3*02:02:01:02', 'DRB3*03:01:01'],
        'hla_drb1': ['DRB1*11:04:01', 'DRB1*13:02:01'],
        'hla_dqa1': ['DQA1*05:05:01:01SG', 'DQA1*01:02:01:04SG'],
        'hla_dqb1': ['DQB1*03:01:01:02', 'DQB1*06:09:01'],
        'hla_dpa1': ['DPA1*01:03:01:01', 'DPA1*01:03:01:04'],
        'hla_dpb1': ['DPB1*02:01:02', 'DPB1*04:02:01:01'],
    }, posterior=0.650)

    print("=" * 72)
    print("  HaploStats — Match Grader (HapLogic Standards)")
    print("=" * 72)

    for donor_name, donor in [("PERFECT MATCH", donor_perfect),
                               ("HAPLOIDENTICAL", donor_half),
                               ("MISMATCH", donor_mismatch)]:
        print(f"\n{'─' * 72}")
        print(f"  SCENARIO: Patient vs Donor ({donor_name})")
        print(f"{'─' * 72}")

        result = grader.compare_profiles(patient, donor, loci=EXTENDED_MATCH_LOCL)

        print(f"\n  Overall Grade: {result['overall']['overall_grade']}")
        print(f"  Compatibility: {result['overall']['compatibility']}")
        print(f"  Allele matches: {result['overall']['allele_matches']}/"
              f"{result['overall']['total_loci_scored']}")
        print(f"  Potential: {result['overall']['potential_matches']}, "
              f"Mismatches: {result['overall']['mismatches']}")
        print()
        print(f"  {'Locus':15} {'Grade':>6} {'Patient':>30} {'Donor':>30}")
        print(f"  {'-'*81}")
        for label, lr in result['loci'].items():
            p_alleles = '/'.join(lr['patient_alleles'][:2])
            d_alleles = '/'.join(lr['donor_alleles'][:2])
            if len(p_alleles) > 28:
                p_alleles = p_alleles[:25] + '...'
            if len(d_alleles) > 28:
                d_alleles = d_alleles[:25] + '...'
            print(f"  {label:15} {lr['grade']:>6} {p_alleles:>30} {d_alleles:>30}")

    return grader


if __name__ == '__main__':
    grader = demo()
    print(f"\n  ✅ MatchGrader ready.")
