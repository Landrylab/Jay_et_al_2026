"""Resolve conflicting signed FungAMR evidence before defining markers."""
import math
import re


SIGNED_EVIDENCE_RULE = (
    'Within each species, canonical reference gene, exact mutation and drug, '
    'retain the nonzero score with the smallest absolute value; positive '
    'resistance evidence wins equal-magnitude ties. Resolve signed evidence '
    'before applying the selected positive confidence-score threshold.'
)


def is_exact_single_substitution(record):
    """Only standalone single-AA records can define exact-substitution markers.

    A component of a multi-substitution record is a separate evidence context;
    its confidence score must never be promoted to the component mutation.
    """
    return (
        record.get('mutation_type', '').upper() in {'SNP', 'SNPS'}
        and record.get('mutation_composition', '').strip().casefold() == 'single'
        and re.fullmatch(r'[A-Z*]\d+[A-Z*]', record.get('mutation', '').strip().upper()) is not None
    )


def best_signed_score(scores):
    """Return the best signed score, or None when no finite nonzero score exists.

    The caller must group equivalent evidence before passing the score iterable.
    Exact substitutions and complete mutation combinations are separate groups.
    """
    eligible = []
    for value in scores:
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(value) and value:
            eligible.append(value)
    if not eligible:
        return None
    result = min(eligible, key=lambda value: (abs(value), value < 0))
    return int(result) if result.is_integer() else result


DEFAULT_CONFIDENCE_SCORES = (1, 2, 3, 4)
EXCLUDED_GENE_MUTATIONS = {('ERG11', 'D116E'), ('ERG11', 'V488I')}


def validate_confidence_scores(scores):
    """Accept original positive database scores; signed resolution precedes filtering."""
    selected = tuple(sorted(set(scores)))
    if not selected or any(not isinstance(v, int) or v < 1 or v > 8 for v in selected):
        raise ValueError('Confidence scores must be integers from 1 to 8')
    return selected


def is_excluded_substitution(gene, mutation):
    """Apply the marker definition documented in Methods."""
    return (str(gene).upper(), str(mutation).upper()) in EXCLUDED_GENE_MUTATIONS
