"""Resolve conflicting signed FungAMR evidence before defining markers."""
import math
import re


SIGNED_EVIDENCE_RULE = (
    'Within each species, canonical reference gene, exact mutation and drug, '
    'retain the nonzero score with the smallest absolute value; positive '
    'resistance evidence wins equal-magnitude ties. Resolve signed evidence '
    'before restricting high-confidence resistance markers to scores +1 or +2.'
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
