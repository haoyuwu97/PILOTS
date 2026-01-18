Developer notes
===============

Adding a new measure
--------------------

PILOTS is configured so that adding a new measure typically requires adding a
single ``.cpp`` file under ``src/measures/``.

Start from:

* ``src/measures/template_measure.cpp``

Workflow
~~~~~~~~

1. Copy the template to a new file (e.g. ``RegisterISF.cpp``)
2. Implement ``caps`` and ``create``
3. Enable registration at the bottom with a ``static MeasureRegistrar``
4. Build and run with:

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

If your repository uses GitHub Pages (see below), you can edit any ``.rst`` file
from the GitHub web UI ("Edit this file"), commit to ``main``, and the website
will be rebuilt and deployed automatically.

Publishing the manual with GitHub Pages
---------------------------------------

This repository includes a workflow that builds the manual and deploys it to
GitHub Pages:

* ``.github/workflows/docs-pages.yml``

One-time setup (repository settings):

1. Go to **Settings** -> **Pages**.
2. Under **Build and deployment**, set **Source** to **GitHub Actions**.

After that, every push to ``main`` that touches ``docs/sphinx/`` triggers a docs
build and deployment.

The published manual will be available at:

* ``https://<owner>.github.io/<repo>/``

