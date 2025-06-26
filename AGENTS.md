# AGENTS

## build-environment
Role: Prepare and validate the PyMOL build environment.
Responsibilities:
- Install all system packages (apt-get).
- Initialize Git submodules.
- Upgrade pip and install Python build dependencies.
- Build and install the PyMOL core and extensions.
Failure Mode: Stop on any missing headers or pip errors and report which package is missing.

## code-assistant
Role: Implement, refactor, and troubleshoot PyMOL plugin scripts.
Responsibilities:
- Operate on files under `script_repo/`.
- Follow PyMOL plugin conventions (proper imports, naming).
- Use existing `COMMANDS.md` as guidance for usage documentation.
- Avoid touching CMakeLists.txt or core source code unless explicitly asked.
Failure Mode: Flag any import or API mismatches with the core.

## documentation-assistant
Role: Maintain in-repo docs for usage, installation, and extensions.
Responsibilities:
- Update `md_docs/` markdown with clear headings and links.
- Keep `README.md` in sync with setup script and AGENTS.md.
- Ensure code examples in docs match actual APIs.
Failure Mode: Report any code snippets that fail to render or run.
