"""Presentation-only normalization; retain genuinely different interpretations."""
import re


def canonical_fragment_label(value: str) -> str:
    # 18:3,O / 18:3,O2 / 18:3;O are equivalent oxygen-count spellings.
    def oxygen(match):
        chain, before, after = match.group(1),match.group(2),match.group(3)
        return f'{chain}({before or after or "1"}O)'
    normalized=re.sub(r'(\d+:\d+)[,;](\d*)O(\d*)(?![A-Za-z0-9])',oxygen,str(value))
    labels=[]
    for label in normalized.split(' | '):
        label=label.strip()
        if label and label not in labels:labels.append(label)
    return ' | '.join(labels)

