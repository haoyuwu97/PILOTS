Adding a measure
================

PILOTS is designed so that adding a new physical observable usually means adding a **single C++ file** under
``src/measures/``.
The runner, registry, selection plumbing, and output indexing are already provided by the framework.

The recommended workflow
------------------------

1. Copy ``src/measures/template_measure.cpp`` into a new file, for example ``src/measures/RegisterAlpha2.cpp``.
2. Implement the physics.
3. Keep registration in the same translation unit.
4. Reference the new measure in the INI config with a ``[measure.<name>]`` section.
5. Rebuild and validate the config.

Single translation unit registration
------------------------------------

PILOTS uses an explicit registry. For a new measure type, keep everything in a single translation unit:

* a ``caps(const MeasureBuildEnv&)`` function,
* a ``create(const MeasureBuildEnv&, const ParsedSection&)`` factory,
* the ``IMeasure`` implementation.

Then add one line at the bottom of the same ``.cpp`` file:

.. code-block:: cpp

   static MeasureRegistrar reg_alpha2(Alpha2Measure::kType, alpha2_create, alpha2_caps);

This style keeps registration simple and avoids macros.

Capabilities checklist
----------------------

Every measure should declare:

* ``selection_policy``: most time-correlation observables should require static selection.
* ``requires_fields``: core or extra fields required from the dump.
* ``requires_topology``: topology sections needed (bonds/angles/dihedrals/impropers).
* ``requires_identity_consistent``: true for most time-correlation functions.
* ``scale_compatibility``: AA/UA/CG.
* ``requires_mapping``: set true if the observable should run on bead-level data.

Respect --validate-config (dry-run)
-----------------------------------

When ``--validate-config`` is used, PILOTS instantiates measures and runs ``on_start()`` to validate requirements.
A measure must not create outputs or perform filesystem side effects in this mode.

Use the SDK ``TextWriter`` with ``dry_run`` enabled (it becomes a no-op writer), and avoid manual ``mkdir`` calls.

Example: Non-Gaussian parameter ``alpha2(t)``
---------------------------------------------

Definition
~~~~~~~~~~

In ``d`` dimensions, the self non-Gaussian parameter is

.. math::

   \alpha_2(t) = \frac{d}{d+2} \frac{\langle |\Delta r(t)|^4 \rangle}{\langle |\Delta r(t)|^2 \rangle^2} - 1.

In 3D (``d=3``):

.. math::

   \alpha_2(t) = \frac{3}{5} \frac{\langle r^4 \rangle}{\langle r^2 \rangle^2} - 1.

This observable depends on two moments as a function of lag time: ``<r2(t)>`` and ``<r4(t)>``.

Implementation outline in PILOTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to implement ``alpha2(t)`` in the current correlator infrastructure is to reuse the T6 correlator:

* push the selected particle positions into the correlator exactly like ``MSDMeasure`` does,
* provide a pair operator that computes two scalars per lag and stores them in two Tensor6 channels:

  * ``XX`` = ``<r2>`` (mean over particles at a given time origin)
  * ``YY`` = ``<r4>`` (mean over particles at a given time origin)

Then compute ``alpha2`` from the two stored channels during flush/output.

Pair operator sketch
~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   struct Alpha2PairOpT6 {
     Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
       const std::size_t n = cur.xx.size();
       double sum_r2 = 0.0;
       double sum_r4 = 0.0;
   #if PILOTS_HAS_OPENMP
   #pragma omp parallel for reduction(+:sum_r2,sum_r4)
   #endif
       for (std::size_t p = 0; p < n; ++p) {
         const double dx = cur.xx[p] - org.xx[p];
         const double dy = cur.yy[p] - org.yy[p];
         const double dz = cur.zz[p] - org.zz[p];
         const double r2 = dx*dx + dy*dy + dz*dz;
         sum_r2 += r2;
         sum_r4 += r2*r2;
       }
       Tensor6 out;
       out.v[T6_XX] = sum_r2 / static_cast<double>(n);
       out.v[T6_YY] = sum_r4 / static_cast<double>(n);
       return out;
     }
   };

Flush/output sketch
~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   const CorrelationSeriesT6 series = corr_->snapshot();
   for (std::size_t i = 0; i < series.lags.size(); ++i) {
     const double r2 = series.values[i].v[T6_XX];
     const double r4 = series.values[i].v[T6_YY];
     const double denom = r2 * r2;
     const double alpha2 = (denom > 0.0) ? (3.0/5.0 * r4 / denom - 1.0) : 0.0;
     // write: lag, time, r2, r4, alpha2, count_pairs, ...
   }

Config example
~~~~~~~~~~~~~~

.. code-block:: ini

   [measure.alpha2]
   type = alpha2
   enabled = true
   group = all
   topo_group = all
   combine = A & T
   output = alpha2.txt

   correlator.type = multitau
   correlator.lag_axis = timestep
   correlator.dt = 0.001

Notes
~~~~~

* ``alpha2(t)`` should almost always use static selection.
* For AA/UA, you will often want bead mapping first, and compute ``alpha2`` on beads rather than atoms.
* Error bars for ``alpha2(t)`` require more care than MSD because it is a ratio of moments. Start by outputting ``r2`` and ``r4``.
