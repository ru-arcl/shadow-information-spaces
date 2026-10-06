/* Mulberry32 PRNG, bit-for-bit identical to shadowinfo/rng.py (DESIGN §5). */
(function (SI) {
  "use strict";

  function Rng(seed) {
    this.state = (seed === undefined ? 1 : seed) >>> 0;
  }

  Rng.prototype.nextU32 = function () {
    this.state = (this.state + 0x6d2b79f5) >>> 0;
    var t = this.state;
    t = Math.imul(t ^ (t >>> 15), t | 1) >>> 0;
    t = (t ^ ((t + (Math.imul(t ^ (t >>> 7), t | 61) >>> 0)) >>> 0)) >>> 0;
    return (t ^ (t >>> 14)) >>> 0;
  };

  Rng.prototype.next_u32 = Rng.prototype.nextU32;

  Rng.prototype.random = function () {
    return this.nextU32() / 4294967296;
  };

  /* Integer in [0, n): floor(u32 * n / 2^32), exact for n < 2^21. */
  Rng.prototype.randint = function (n) {
    if (!(n > 0 && n < 2097152)) throw new RangeError("randint(n) requires 0 < n < 2^21");
    return Math.floor((this.nextU32() * n) / 4294967296);
  };

  SI.Rng = Rng;
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
