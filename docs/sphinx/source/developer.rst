Developer Notes
===============

Adding a new measure
--------------------

PILOTS is configured so that adding a new measure typically requires adding a single ``.cpp`` file under
``src/measures/``.

Start from:

* ``src/measures/template_measure.cpp``

Workflow
~~~~~~~~

1. Copy the template to a new file, for example ``RegisterISF.cpp``.
2. Implement ``caps`` and ``create``.
3. Enable registration at the bottom with a ``static MeasureRegistrar``.
4. Rebuild and confirm the new type appears in:

   .. code-block:: bash

      ./build/pilots --list-measures

Building the manual
-------------------

The long-term documentation workflow uses Sphinx.

.. code-block:: bash

   ./docs/build_docs.sh

This writes HTML to ``docs/manual/``.

Editing the manual on GitHub
----------------------------

The Sphinx sources live under ``docs/sphinx/source/``.

If GitHub Pages is enabled for the repository, you can edit any ``.rst`` file from the GitHub web UI,
commit to ``main``, and let the site rebuild from source.

Publishing the manual with GitHub Pages
---------------------------------------

This repository includes a workflow that builds the manual and deploys it to GitHub Pages:

* ``.github/workflows/docs-pages.yml``

One-time setup in repository settings:

1. Go to **Settings** → **Pages**.
2. Under **Build and deployment**, set **Source** to **GitHub Actions**.

After that, pushes to ``main`` that touch ``docs/sphinx/`` trigger a docs build and deployment.
The published manual is available at:

* ``https://haoyuwu97.github.io/PILOTS/``
