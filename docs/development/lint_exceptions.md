# nf-core lint exceptions

Every entry in the `lint` section of `.nf-core.yml` must be justified here. Add an entry when an exception is introduced and remove it when the exception is no longer needed.

## Exceptions set by the template

The nf-core/tools 4.2.2 template adds these entries for custom (non-nf-core) pipelines. They are kept as generated.

| Lint test         | Entry                                                                                                                                                                | Reason                                                                        |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| `files_unchanged` | `CODE_OF_CONDUCT.md`, `assets/email_template.txt`, `.github/ISSUE_TEMPLATE/bug_report.yml`, `.github/PULL_REQUEST_TEMPLATE.md`, `docs/README.md`, nf-core logo files | Custom pipelines do not carry nf-core branding or the nf-core code of conduct |
| `files_exist`     | nf-core logo files, `CODE_OF_CONDUCT.md`, `.github/ISSUE_TEMPLATE/config.yml`, AWS test workflows, iGenomes configs                                                  | Not generated for custom pipelines, or skipped features (`igenomes`)          |
| `nextflow_config` | `manifest.name`, `manifest.homePage`                                                                                                                                 | The name and home page are not under the `nf-core/` organisation              |
| `multiqc_config`  | `report_comment`                                                                                                                                                     | The report comment differs from the nf-core wording                           |

## Exceptions added for this pipeline

| Lint test         | Entry                        | Reason                                                                                                                                                              | Added      |
| ----------------- | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------- |
| `files_unchanged` | `assets/email_template.html` | The HTML email embeds the nf-core pipeline logo, which custom pipelines do not have. The template already exempts the plain-text twin, `assets/email_template.txt`. | 2026-11-14 |

## How to lint

Run `dev/nfdev lint` (developer machine). It lints a git snapshot of the staged files, which matches what CI checks and ignores untracked local files.
