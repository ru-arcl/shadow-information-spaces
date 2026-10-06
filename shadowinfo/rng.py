"""Tiny deterministic PRNG shared bit-for-bit with the JS demo (DESIGN §5).

Mulberry32: a 32-bit state, 32-bit output generator.  Every operation is done
modulo ``2**32`` so that :mod:`docs/js/rng.js` (which uses ``Math.imul`` and
``>>> 0``) produces exactly the same integer stream for the same seed.
"""

from __future__ import annotations

_MASK = 0xFFFFFFFF
_TWO32 = 4294967296


def _imul(a: int, b: int) -> int:
    """32-bit wrapping multiplication (JS ``Math.imul`` on unsigned operands)."""
    return (a * b) & _MASK


class Rng:
    """Mulberry32 generator.

    ``next_u32()`` returns an integer in ``[0, 2**32)``, ``random()`` a float in
    ``[0, 1)`` equal to ``u32 / 2**32`` and ``randint(n)`` an integer in
    ``[0, n)`` computed as ``floor(u32 * n / 2**32)`` (exact for ``n < 2**21``).
    """

    def __init__(self, seed: int = 1) -> None:
        self.state = int(seed) & _MASK

    def next_u32(self) -> int:
        self.state = (self.state + 0x6D2B79F5) & _MASK
        t = self.state
        t = _imul(t ^ (t >> 15), t | 1)
        t ^= (t + _imul(t ^ (t >> 7), t | 61)) & _MASK
        return (t ^ (t >> 14)) & _MASK

    def random(self) -> float:
        return self.next_u32() / _TWO32

    def randint(self, n: int) -> int:
        if not 0 < n < (1 << 21):
            raise ValueError("randint(n) requires 0 < n < 2**21")
        return (self.next_u32() * n) >> 32
