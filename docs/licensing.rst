License recommendation and scope
================================

BSD-3-Clause is recommended for the new independent Python package. It permits
academic and commercial use and redistribution, requires preservation of the
copyright/license notices, and contains an explicit non-endorsement condition.
Keep the standard text intact. Citation requests for both primary papers belong
in the README and citation files as scholarly credit, separate from the standard
license terms.

MIT offers a similarly short permissive grant. Apache-2.0 includes an express
patent license and additional notice provisions. BSD-3-Clause fits the requested
broad distribution and is also the license used by Trekker. This recommendation
is not an institutional ownership or release-authorization determination.

Proposal status
---------------

The local licensing/ folder contains the proposed BSD text and the future
packaging metadata fields. The copyright-holder field is deliberately unresolved.
No active LICENSE file or License-Expression is asserted while that remains open.
The --require-license archive check should fail in this preparation state.
The proposal directory is excluded from Python distributions and normal Git
tracking. Move an approved, completed text to top-level LICENSE before activation.

The eventual package metadata would use license = "BSD-3-Clause" and
license-files = ["LICENSE"], with setuptools >=77.0.3. Preserve existing
third-party licenses when shipping any such code or bundled environments.
The independent source archive does not vendor the declared numerical libraries;
see THIRD_PARTY_NOTICES.rst for the direct dependency inventory.

Trekker comparison
------------------

Trekker's LICENSE.md declares BSD-3-Clause and a 2021--2025 copyright notice naming
Dogu Baran Aydogan. Its README requests paper citations in a separate Publications
section. This is a useful format precedent, not evidence of SurfTracker's owner
or of the approvals underlying either project. Its documentation-site footer
is also separate from the software's legal notice.

* https://github.com/dmritrekker/trekker/blob/main/LICENSE.md
* https://github.com/dmritrekker/trekker/blob/main/README.md
* https://spdx.org/licenses/BSD-3-Clause.html
* https://www.apache.org/licenses/LICENSE-2.0
* https://packaging.python.org/en/latest/guides/licensing-examples-and-user-scenarios/

Research attribution and paper rights
-------------------------------------

AUTHORS.rst records the human roles supplied by Xinyu Nie and acknowledges
OpenAI Codex/GPT development assistance separately. Neither paper coauthorship
nor AI-tool acknowledgment is used to assign legal software ownership. The
research group and institute are recorded as affiliations.

CITATION.cff represents the software with both primary papers in references;
there is no single preferred-citation override selecting only one paper. The
README explicitly requests both, and CITATIONS.bib provides both entries.
Journal years are 2008 and 2024; DOI years remain 2007 and 2023 respectively.
The paper PDFs are not redistributed or relicensed by this package.

See ownership_review.rst for the limited institutional question that remains.
