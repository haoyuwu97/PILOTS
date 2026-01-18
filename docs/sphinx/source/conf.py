import os
import sys

# -- Path setup --------------------------------------------------------------

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, ROOT)

# -- Project information -----------------------------------------------------

project = "PILOTS"
copyright = ""
author = ""

# The short X.Y version.
release = os.environ.get("PILOTS_VERSION", "dev")

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = []

language = "en"

# -- Options for HTML output -------------------------------------------------

try:
    import sphinx_rtd_theme  # noqa: F401
    html_theme = "sphinx_rtd_theme"
except Exception:
    html_theme = "alabaster"

html_static_path = ["_static"]

html_title = "PILOTS Manual"

# Optional: show "Edit on GitHub" links when building in GitHub Actions.
# This requires sphinx-rtd-theme, but is harmless if the theme is not available.
repo = os.environ.get("GITHUB_REPOSITORY")  # e.g. "owner/repo"
ref_name = os.environ.get("GITHUB_REF_NAME", "main")
if repo and "/" in repo:
    owner, name = repo.split("/", 1)
    html_context = {
        "display_github": True,
        "github_user": owner,
        "github_repo": name,
        "github_version": ref_name,
        "conf_py_path": "/docs/sphinx/source/",
    }

# Intersphinx (optional; helps link to external docs).
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", "https://docs.python.org/3/objects.inv"),
}

python_intersphinx = intersphinx_mapping.get("python")
if (
    python_intersphinx
    and isinstance(python_intersphinx, tuple)
    and len(python_intersphinx) >= 2
):
    python_url, python_inventory = python_intersphinx[0], python_intersphinx[1]
    if not python_inventory or isinstance(python_inventory, dict):
        intersphinx_mapping["python"] = (
            python_url,
            f"{python_url.rstrip('/')}/objects.inv",
        )
