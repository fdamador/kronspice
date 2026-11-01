# Licensing

kronspice is **dual-licensed**. You may use it under either of two licenses:

1. **GNU General Public License, version 3 only (GPL-3.0-only)** — free and
   open source. The full text is in [LICENSE](LICENSE).
2. **A commercial license** from the copyright holder — for uses the GPL does
   not allow or that you prefer not to make under GPL terms. See
   [COMMERCIAL-LICENSE.md](COMMERCIAL-LICENSE.md).

Each source file carries the SPDX identifier:

```
SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
```

## Which license do I need?

| Your situation | License |
|---|---|
| Research, study, personal use, internal evaluation | GPL-3.0 is enough |
| You distribute software containing kronspice and release your whole program's source under GPL-3.0 | GPL-3.0 is enough |
| You distribute kronspice, or a product that includes or links it, **without** releasing your source under GPL-3.0 (closed-source software, firmware, embedded devices, SaaS appliances shipped to customers) | Commercial license |
| You want warranty, indemnity, support or other terms the GPL does not provide | Commercial license |

Go programs are statically linked, so a program that imports any kronspice
package is generally a combined work under the GPL. When in doubt, ask.

This table is a plain-language guide, not legal advice; the license texts
govern.

## Contributions

To keep dual licensing possible, the copyright holder must be able to offer
every part of the code under both licenses. Contributions are therefore
accepted only under a Contributor License Agreement — see [CLA.md](CLA.md)
and [CONTRIBUTING.md](CONTRIBUTING.md).

## Third-party code

kronspice currently contains no third-party code. Any third-party code added
later must carry a license compatible with *both* GPL-3.0 and the commercial
license (for example MIT, BSD or Apache-2.0, with its notice kept). Code
under the GPL or other copyleft licenses must not be copied in, because it
could not be offered commercially (CLAUDE.md §6).
