# Contributing

Thank you for your interest in kronspice. Please read
[CLAUDE.md](CLAUDE.md) first: it describes the design, the coding rules and
the milestones, and every change must keep the tests passing
(`go test -race ./...`).

## Contributor License Agreement (required)

kronspice is dual-licensed under GPL-3.0-only and a commercial license
([LICENSING.md](LICENSING.md)). So that every part of the code can be offered
under both, contributions are accepted only from contributors who have
signed the project's Contributor License Agreement (CLA). In short, the CLA:

- confirms that you wrote the contribution, or otherwise have the right to
  submit it;
- lets the copyright holder license your contribution under GPL-3.0 **and**
  under commercial terms;
- does not take away your own right to use your contribution.

The agreement is [CLA.md](CLA.md). Sign it by email to fdamador@comcast.net as
described there before opening your first pull request. Pull requests
without a signed CLA cannot be merged.

## Rules for code

- Write the code yourself or take it from a source whose license is
  compatible with both GPL-3.0 and commercial use (MIT, BSD, Apache-2.0),
  keeping its notice. Never copy GPL or other copyleft code (CLAUDE.md §6).
- Start every new Go file with the license header:

  ```go
  // SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
  // Copyright (C) 2026 Franklin Amador
  ```
