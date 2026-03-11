Measures Reference
==================

This page is the **measure quick reference** for the current PILOTS source tree.
It is organized by measure family and should be updated in the same commit as the
corresponding C++ files under ``src/measures/``.

For each measure or family, the page records:

* source file(s)
* required input fields and topology dependencies
* important config keys
* output columns
* physical meaning
* typical use cases

Common conventions
------------------

Static selection and identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most time-correlation and displacement measures require:

* a **static** selection (configured by ``group``, ``topo_group`` and ``combine``)
* **identity-consistent** frames, i.e. the same canonical particle ordering across frames

Frame-window controls
~~~~~~~~~~~~~~~~~~~~~

Many measures support:

* ``frame_start`` / ``frame_end``: analysis window in frame index space
* ``remove_drift`` / ``drift_group``: subtract COM drift before evaluating displacements
* ``components``: choose ``xx``, ``yy``, ``zz``, ``xxyy``, ``xxzz``, ``yyzz``, or ``xxyyzz``

Common list-valued parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``q_values``: wave-vector magnitudes
* ``a_values``: overlap cutoff radii
* ``p_values``: Rouse mode indices
* ``waiting_frames`` / ``lag_frames``: explicit waiting-time and lag grids for two-time observables

MSD and single-particle glassy dynamics
---------------------------------------

``src/measures/RegisterMSD.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``msd``
^^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``

**Dependencies**
  static selection; identity-consistent frames; optional static ``drift_group``

**Outputs**
  lag, time, tensor components, count-pairs, SEM, block statistics

**Physical meaning**
  Mean-squared displacement (or selected Cartesian tensor block) averaged over time origins.

**Typical use**
  Diffusion, caging, subdiffusion, anisotropic transport.

.. math::

   \mathrm{MSD}(t)=\left\langle \frac{1}{N}\sum_i |\mathbf r_i(t_0+t)-\mathbf r_i(t_0)|^2\right\rangle_{t_0}

``src/measures/RegisterAlpha2.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``alpha2``
^^^^^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``

**Dependencies**
  static selection; identity-consistent frames; optional drift removal; 1D correlator infrastructure

**Outputs**
  lag, time, :math:`\langle r^2\rangle`, :math:`\langle r^4\rangle`, :math:`\alpha_2`, count-pairs, SEM

**Physical meaning**
  Non-Gaussian parameter for self displacements.

**Typical use**
  Detect departure from Gaussian diffusion, caging-to-hopping crossover, dynamic heterogeneity onset.

.. math::

   \alpha_2(t)=\frac{3}{5}\frac{\langle r^4(t)\rangle}{\langle r^2(t)\rangle^2}-1

``src/measures/RegisterScattering.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``self_isf``
^^^^^^^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``

**Dependencies**
  static selection; identity-consistent frames; optional drift removal; ``q_values``

**Outputs**
  lag, time, q, :math:`F_s(q,t)`, count-pairs, SEM

**Physical meaning**
  Isotropic self intermediate scattering function.

.. math::

   F_s(q,t)=\left\langle\frac{1}{N}\sum_i \frac{\sin(q|\Delta\mathbf r_i|)}{q|\Delta\mathbf r_i|}\right\rangle_{t_0}

``chi4_overlap``
^^^^^^^^^^^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``

**Dependencies**
  static selection; identity-consistent frames; optional drift removal; ``a_values``

**Outputs**
  lag, time, a, :math:`Q`, :math:`Q^2`, :math:`\chi_4`, count-pairs, SEM

**Physical meaning**
  Four-point susceptibility based on the overlap indicator.

.. math::

   Q_a(t)=\left\langle\frac{1}{N}\sum_i \Theta(a-|\Delta\mathbf r_i|)\right\rangle_{t_0},
   \qquad
   \chi_4(a,t)=N\left[\langle Q_a^2\rangle-\langle Q_a\rangle^2\right]

``chi4_isf``
^^^^^^^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``

**Dependencies**
  static selection; identity-consistent frames; optional drift removal; ``q_values``

**Outputs**
  lag, time, q, :math:`\overline F_s`, :math:`\overline{F_s^2}`, :math:`\chi_4^{F_s}`

**Physical meaning**
  Four-point susceptibility built from the self ISF amplitude.

.. math::

   \chi_4^{F_s}(q,t)=N\left[\overline{F_s^2}(q,t)-\overline F_s(q,t)^2\right]

``src/measures/RegisterTwoTime.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file implements **aging / two-time observables**. These are needed when the system is not
stationary and the observable depends on both the waiting time :math:`t_w` and the lag :math:`\tau`.

The family uses an explicit waiting-time grid and does **not** collapse over all time origins like a
stationary correlator. Common keys are:

* ``waiting_frames`` or ``waiting_start_frame`` / ``waiting_end_frame`` / ``waiting_stride_frames``
* ``waiting_window_frames``: number of consecutive origins averaged inside each waiting-time bin
* ``origin_stride_frames``: subsampling stride inside the waiting-time bin
* ``lag_frames`` or ``max_lag_frames`` / ``lag_stride_frames``
* finite ``frame_end`` is required; follow mode is intentionally rejected

``msd_tw``
^^^^^^^^^^

**Outputs**
  waiting frame, waiting timestep, waiting time, lag in frames and timesteps, lag time,
  two-time MSD, SEM across origins inside the waiting-time bin, number of origins

.. math::

   \mathrm{MSD}(t_w,\tau)=\left\langle \frac{1}{N}\sum_i
   |\mathbf r_i(t_w+\tau)-\mathbf r_i(t_w)|^2\right\rangle_{t_w\text{-bin}}

``fs_tw``
^^^^^^^^^

**Additional key**
  ``q_values``

**Outputs**
  waiting frame/time, lag frame/time, q, :math:`F_s(q;t_w,\tau)`, SEM, number of origins

.. math::

   F_s(q;t_w,\tau)=\left\langle \frac{1}{N}\sum_i
   \frac{\sin(q|\mathbf r_i(t_w+\tau)-\mathbf r_i(t_w)|)}{q|\mathbf r_i(t_w+\tau)-\mathbf r_i(t_w)|}
   \right\rangle_{t_w\text{-bin}}

``q_tw``
^^^^^^^^

**Additional key**
  ``a_values``

**Outputs**
  waiting frame/time, lag frame/time, a, :math:`Q(t_w,\tau)`, SEM, number of origins

.. math::

   Q_a(t_w,\tau)=\left\langle \frac{1}{N}\sum_i
   \Theta\!\left(a-|\mathbf r_i(t_w+\tau)-\mathbf r_i(t_w)|\right)
   \right\rangle_{t_w\text{-bin}}

``chi4_tw``
^^^^^^^^^^^

**Additional keys**
  ``a_values`` and ``waiting_window_frames >= 2``

**Outputs**
  waiting frame/time, lag frame/time, a, :math:`Q`, :math:`Q^2`, :math:`\chi_4(t_w,\tau)`,
  SEM of :math:`Q`, number of origins in the waiting-time bin

**Physical meaning**
  Waiting-time-resolved fluctuation amplitude of the overlap field.

.. math::

   \chi_4(t_w,\tau)=N\left[\langle Q_a(t_w,\tau)^2\rangle_{t_w\text{-bin}}-
   \langle Q_a(t_w,\tau)\rangle_{t_w\text{-bin}}^2\right]

``waiting_time_scan``
^^^^^^^^^^^^^^^^^^^^^

**Additional key**
  ``observable = msd|fs|q|chi4``

**Meaning**
  Convenience front-end that reuses the same kernels as ``msd_tw``, ``fs_tw``, ``q_tw`` or
  ``chi4_tw`` while keeping a single config style for waiting-time scans.

Polymer dynamics and static chain measures
------------------------------------------

``src/measures/RegisterPolymerDynamics.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``g1``
^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``

**Dependencies**
  static selection; identity-consistent frames; optional drift removal

.. math::

   g_1(t)=\left\langle\frac{1}{\sum_c N_c}\sum_c\sum_n
   |\mathbf r_{c,n}(t_0+t)-\mathbf r_{c,n}(t_0)|^2\right\rangle_{t_0}

``g2``
^^^^^^

**Additional dependencies**
  chain membership from ``chain_id_field`` or ``mol`` or topology bonds

.. math::

   g_2(t)=\left\langle\frac{1}{\sum_c N_c}\sum_c\sum_n
   |(\mathbf r_{c,n}-\mathbf R_{\mathrm{cm},c})(t_0+t)-
   (\mathbf r_{c,n}-\mathbf R_{\mathrm{cm},c})(t_0)|^2\right\rangle_{t_0}

``g3``
^^^^^^

**Additional dependencies**
  chain membership from ``chain_id_field`` or ``mol`` or topology bonds

.. math::

   g_3(t)=\left\langle\frac{1}{N_{\mathrm{chain}}}\sum_c
   |\mathbf R_{\mathrm{cm},c}(t_0+t)-\mathbf R_{\mathrm{cm},c}(t_0)|^2\right\rangle_{t_0}

``rouse``
^^^^^^^^^

**Additional dependencies**
  chain membership plus linear chain ordering from ``chain_pos_field`` or topology bonds

**Additional key**
  ``p_values``

**Outputs**
  lag, time, mode index, :math:`C_p(t)`, normalized :math:`C_p(t)/C_p(0)`

.. math::

   \mathbf X_{p,c}(t)=\frac{1}{N_c}\sum_{n=0}^{N_c-1}\mathbf r_{c,n}(t)
   \cos\!\left[\frac{p\pi(n+1/2)}{N_c}\right]

   C_p(t)=\left\langle\frac{1}{N_{\mathrm{chain}}}\sum_c
   \mathbf X_{p,c}(t_0+t)\cdot \mathbf X_{p,c}(t_0)\right\rangle_{t_0}

``src/measures/RegisterPolymerStatics.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``rg_ree``
^^^^^^^^^^

**Dependencies**
  linear chain ordering from ``chain_pos_field`` or topology bonds

**Outputs**
  frame, time, :math:`\langle R_g^2\rangle`, :math:`\langle R_g\rangle`,
  :math:`\langle R_{ee}^2\rangle`, :math:`\langle R_{ee}\rangle`

.. math::

   R_g^2=\frac{1}{N}\sum_n |\mathbf r_n-\mathbf R_{\mathrm{cm}}|^2,
   \qquad
   R_{ee}^2=|\mathbf r_N-\mathbf r_1|^2

Structure and van Hove family
-----------------------------

``src/measures/RegisterStructure.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``rdf``
^^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``

**Dependencies**
  static selection; optional second group ``group_b`` / ``topo_group_b`` / ``combine_b``

**Outputs**
  r-bin, shell volume, :math:`g(r)`, cumulative coordination number

.. math::

   g_{AB}(r)=\frac{\langle n_{AB}(r)\rangle}{N_A\rho_B\Delta V(r)},
   \qquad
   \Delta V(r)=\frac{4\pi}{3}(r_{hi}^3-r_{lo}^3)

``static_structure_factor`` / ``sq_static``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   S(q)=1+\frac{2}{N}\sum_{i<j}\frac{\sin(qr_{ij})}{qr_{ij}}

``dynamic_structure_factor`` / ``coherent_isf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Outputs**
  lag, time, q, coherent :math:`F(q,t)`, self part :math:`F_s(q,t)`, distinct part :math:`F_d(q,t)`

.. math::

   F(q,t)=\left\langle \frac{1}{N}\sum_i\sum_j j_0\!\left(q|\mathbf r_j(t_0+t)-\mathbf r_i(t_0)|\right)\right\rangle_{t_0}

``van_hove_self`` and ``van_hove_distinct``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Outputs**
  lag, time, radial bin, shell-averaged :math:`G_s(r,t)` or :math:`G_d(r,t)`

.. math::

   G_s(\mathbf r,t)=\frac{1}{N}\left\langle\sum_i
   \delta\!\left(\mathbf r-[\mathbf r_i(t_0+t)-\mathbf r_i(t_0)]\right)\right\rangle_{t_0}

   G_d(\mathbf r,t)=\frac{1}{N}\left\langle\sum_{i\neq j}
   \delta\!\left(\mathbf r-[\mathbf r_j(t_0+t)-\mathbf r_i(t_0)]\right)\right\rangle_{t_0}

Dielectric and orientation / relaxation family
----------------------------------------------

``src/measures/RegisterDielectric.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``bulk_dielectric``
^^^^^^^^^^^^^^^^^^^

**Required fields**
  ``xu``, ``yu``, ``zu``, charge field (default ``q``)

**Dependencies**
  entity grouping from ``entity_id_field`` or ``mol`` or topology bonds; temperature and unit metadata

**Meaning**
  Classic neutral-entity dipole-fluctuation dielectric. Use this when grouped entities are approximately neutral.

**Outputs**
  per-component permittivity and isotropic average

.. math::

   \varepsilon_\alpha = 1 + \frac{\mathrm{Var}(M_\alpha)}{\varepsilon_0 V k_B T},
   \qquad
   \mathbf M=\sum_e \boldsymbol\mu_e

``bulk_fragment_dielectric``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Additional keys**
  ``charge_model = fragment_internal`` (default for this alias), ``entity_reference = geometric|mass``

**Meaning**
  Conduction-blocked charged-mixture dielectric built from fragment-internal dipoles

.. math::

   \boldsymbol\mu_e^{\mathrm{int}} = \sum_{i\in e} q_i \left(\mathbf r_i-\mathbf R_e\right)

This yields a finite, origin-invariant entity-internal polarization even when grouped fragments carry net charge.
It is useful as a charged-mixture **surrogate/background-polarization dielectric**, but should not be confused
with the full zero-frequency static dielectric of a conducting electrolyte.

``slab_dielectric`` / ``layered_dielectric``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Additional keys**
  ``axis`` and ``n_bins``

**Outputs**
  slab coordinate, :math:`\varepsilon_\parallel(z)`, :math:`\varepsilon_\perp(z)` and related susceptibilities

.. math::

   \chi_{\parallel}(z)=\frac{C_{\parallel}(z)}{\varepsilon_0 V_{\mathrm{bin}} k_B T},
   \qquad
   \chi_{\perp}(z)=\frac{C_{\perp}(z)}{\varepsilon_0 V_{\mathrm{bin}} k_B T}

   \varepsilon_{\parallel}(z)=1+\chi_{\parallel}(z),
   \qquad
   \varepsilon_{\perp}(z)=\frac{1}{1-\chi_{\perp}(z)}

``slab_fragment_dielectric`` / ``layered_fragment_dielectric``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Additional keys**
  ``charge_model = fragment_internal`` (default for this alias), ``entity_reference = geometric|mass`` and the slab-workflow keys
  ``slab_align_recenter``, ``halfcell_fold``, ``target_center_frac``, ``anchor_group`` / ``anchor_topo_group`` / ``anchor_combine``.

**Meaning**
  Charged-mixture, fragment-internal layered dielectric. This is the finite, conduction-blocked local dielectric surrogate for systems where grouped entities may carry net charge.

``slab_bg_dielectric`` / ``layered_bg_dielectric``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Additional keys**
  ``species_groups`` (comma-separated static groups to resolve),
  ``background_groups`` (comma-separated subset whose sum defines the background polarization; defaults to ``species_groups``),
  plus the same ``entity_reference`` and slab-workflow keys used by ``slab_fragment_dielectric``.

**Meaning**
  Species-resolved **background-polarization** layered dielectric inspired by local polarization-density fluctuation methods.
  The measure constructs species-resolved fragment-internal polarization fields, sums the requested
  ``background_groups`` into a total background polarization, and outputs both the total background dielectric
  and per-species susceptibility contributions.

**Outputs**
  For each slab bin, one ``TOTAL_BG`` row plus one row per species group. The total row contains
  :math:`\varepsilon_\parallel^{\mathrm{bg}}(z)` and :math:`\varepsilon_\perp^{\mathrm{bg}}(z)`. Species rows contain
  their partial covariances and susceptibility contributions with respect to ``TOTAL_BG``.

.. math::

   \boldsymbol\mu_{e,g}^{\mathrm{int}} = \sum_{i\in e\cap g} q_i \left(\mathbf r_i-\mathbf R_e\right)

   \mathbf M_g = \sum_e \boldsymbol\mu_{e,g}^{\mathrm{int}},
   \qquad
   \mathbf M_{\mathrm{bg}} = \sum_{g\in \mathrm{background}} \mathbf M_g

   \chi_{g,\parallel}^{\mathrm{bg}}(z)=
   \frac{\mathrm{Cov}\!\left(M_{g,\parallel}^{\mathrm{bin}}(z),\,M_{\mathrm{bg},\parallel}\right)}
        {\varepsilon_0 V_{\mathrm{bin}} k_B T}

   \chi_{g,\perp}^{\mathrm{bg}}(z)=
   \frac{\mathrm{Cov}\!\left(M_{g,\perp}^{\mathrm{bin}}(z),\,M_{\mathrm{bg},\perp}\right)}
        {\varepsilon_0 V_{\mathrm{bin}} k_B T}

   \chi_{\mathrm{bg}}(z)=\sum_{g\in \mathrm{background}} \chi_g^{\mathrm{bg}}(z)

The ``TOTAL_BG`` row then reports :math:`\varepsilon_\parallel^{\mathrm{bg}}(z)=1+\chi_{\parallel}^{\mathrm{bg}}(z)`
and :math:`\varepsilon_\perp^{\mathrm{bg}}(z)=1/[1-\chi_{\perp}^{\mathrm{bg}}(z)]`.

``src/measures/RegisterOrientation.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``bond_vector_acf``
^^^^^^^^^^^^^^^^^^^

**Dependencies**
  topology bonds

``segment_reorientation``
^^^^^^^^^^^^^^^^^^^^^^^^^

**Dependencies**
  ordered chains from ``chain_pos_field`` or topology bonds

``ring_normal_acf``
^^^^^^^^^^^^^^^^^^^

**Dependencies**
  ring identifiers via ``ring_id_field``

``dipole_acf`` and ``dielectric_spectrum``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Dependencies**
  entity grouping from ``entity_id_field`` or ``mol`` or topology bonds, plus charge field

**Outputs**
  lag-resolved orientational or dipolar ACFs; for ``dielectric_spectrum`` also frequency-domain spectrum

**Physical meaning**
  Rotational relaxation, segmental reorientation, ring-plane relaxation, dielectric loss.

The common orientational kernel is

.. math::

   C_l(t)=\left\langle P_l[\hat{\mathbf u}(t_0+t)\cdot\hat{\mathbf u}(t_0)]\right\rangle_{t_0}

For the total dipole family,

.. math::

   C_\mu(t)=\langle \mathbf M(t_0+t)\cdot \mathbf M(t_0)\rangle_{t_0}

Transport family
----------------

``src/measures/RegisterTransport.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``current_acf`` / ``conductivity_gk``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Required fields**
  charge field and a vector field triplet (for example ``vx``, ``vy``, ``vz``)

**Outputs**
  lag, time, current ACF tensor block, and Green--Kubo running conductivity

.. math::

   \mathbf J(t)=\sum_i q_i\mathbf v_i(t),
   \qquad
   \sigma_{\alpha\alpha}^{\mathrm{GK}}(t)=\frac{1}{V k_B T}\int_0^t C_{\alpha\alpha}(\tau)\,d\tau

``conductivity_einstein`` / ``ionic_conductivity_einstein`` / ``charge_msd``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \mathbf M(t;t_0)=\sum_i q_i[\mathbf r_i(t_0+t)-\mathbf r_i(t_0)]

   \sigma_E(t)=\frac{1}{6Vk_B T t}\left\langle|\mathbf M(t;t_0)|^2\right\rangle_{t_0}

``site_hop_stats`` / ``hop_stats``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Required field**
  discrete ``site_field``

**Outputs**
  committed transition counts, dwell-time statistics, hop rates

**Physical meaning**
  State-to-state hopping statistics with persistence filtering.

``src/measures/RegisterTransportSpecies.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``species_msd``
^^^^^^^^^^^^^^^

**Required field**
  integer-like ``species_field``

.. math::

   MSD^{(a)}_{\alpha\beta}(t)=\left\langle\frac{1}{N_a}\sum_{i\in a}
   \Delta r_{i,\alpha}(t;t_0)\Delta r_{i,\beta}(t;t_0)\right\rangle_{t_0}

``diffusion_tensor``
^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Einstein-form diffusion tensor inferred from the species-resolved MSD tensor.

``vacf``
^^^^^^^^

**Required fields**
  species field plus vector field triplet

.. math::

   C^{(a)}_{\alpha\beta}(t)=\left\langle\frac{1}{N_a}\sum_{i\in a}
   v_{i,\alpha}(t_0+t)v_{i,\beta}(t_0)\right\rangle_{t_0}

``species_current_corr``
^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   J_a^\alpha(t)=\sum_{i\in a} q_i v_i^\alpha(t),
   \qquad
   C^{ab}_{\alpha\beta}(t)=\langle J_a^\alpha(t_0+t)J_b^\beta(t_0)\rangle_{t_0}

``onsager_transport``
^^^^^^^^^^^^^^^^^^^^^

.. math::

   M_a^\alpha(t;t_0)=\sum_{i\in a} q_i\Delta r_{i,\alpha}(t;t_0)

   L_{ab}^{\alpha\beta}(t)=\frac{\langle M_a^\alpha M_b^\beta\rangle_{t_0}}{2Vk_B T t}

``transference_number``
^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Transport-number-like quantity constructed from the Einstein-form Onsager matrix block sums.
  The value is reference-frame dependent and should be interpreted together with the selected
  drift-removal / frame convention.

Profiles, coordination, and box statistics
------------------------------------------

``src/measures/RegisterOMIEC.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Even though this file was introduced for OMIEC-oriented workflows, the retained measures are all
**general soft-matter statistics**.

``profile1d``
^^^^^^^^^^^^^

**Modes**
  ``number``, ``charge``, ``mass``

**Important config keys**
  ``axis``, ``n_bins``, ``mode``, and, for slab workflows,
  ``slab_align_recenter``, ``halfcell_fold``, ``target_center_frac``, and
  optional ``anchor_group`` / ``anchor_topo_group`` / ``anchor_combine``.

**Outputs**
  1D bin center, average weight, density in the analysis coordinate. When slab
  alignment or half-cell folding is enabled, the coordinate reported in the text
  output is the transformed slab coordinate rather than the raw box axis.

``coordination``
^^^^^^^^^^^^^^^^

.. math::

   CN_i=\sum_{j\in B}\Theta(r_{\mathrm{cut}}-r_{ij})

``residence_acf``
^^^^^^^^^^^^^^^^^

.. math::

   h_{ij}(t)=\Theta(r_{\mathrm{cut}}-r_{ij}(t)),
   \qquad
   C_{\mathrm{res}}(t)=\left\langle\frac{\sum_{ij}h_{ij}(t_0)h_{ij}(t_0+t)}{\sum_{ij}h_{ij}(t_0)}\right\rangle_{t_0}

``box_stats`` / ``swelling_stats``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Outputs**
  frame, time, box lengths, box volume, and ratios relative to the first analyzed frame


Insertion and slab-kernel calibration family
--------------------------------------------

``src/measures/RegisterInsertionProfiles.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These measures share the same slab-analysis transform used by ``profile1d``:

* ``slab_align_recenter`` recenters a slab-like anchor selection to ``target_center_frac``
  using a circular mean in fractional coordinate space.
* ``halfcell_fold`` folds a periodic electrolyte / OMIEC / electrolyte cell into a single
  reservoir-to-slab half-cell.

``accessibility_profile`` / ``probe_accessibility_profile``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Probe-radius-dependent accessible-volume fraction profile :math:`h(z; r_p)` from hard-core
  insertion geometry.

.. math::

   h(z; r_p)=\left\langle \frac{N_{\mathrm{accessible}}(z; r_p)}{N_{\mathrm{trial}}(z)}\right\rangle

**Important config keys**
  ``probe_radius``, ``radius_field`` or ``occluder_radius``, ``samples_plane``,
  ``samples_axis``, ``convergence_refine_factor``, slab transform keys, and
  reservoir-window keys ``reservoir_lo_frac`` / ``reservoir_hi_frac`` with optional
  right-window counterparts.

**Outputs**
  Absolute accessibility, SEM across analyzed frames, reservoir-normalized
  accessibility, :math:`-\ln[h(z)/h_{\mathrm{res}}]`, and coarse-versus-refined voxel
  convergence diagnostics.

``widom_profile``
^^^^^^^^^^^^^^^^^

**Meaning**
  Total charge-off Widom insertion profile using the configured soft backend.

.. math::

   \mu^{\mathrm{ex}}_{0}(z)=-\beta^{-1}\ln\left\langle e^{-\beta \Delta U_0}\right\rangle_z

**Important config keys**
  All accessibility keys plus ``beta`` and ``energy_model``. The current public
  backend supports ``none`` and full ``lj126`` pair energies, with optional
  ``sigma_field`` / ``epsilon_field`` or constant probe/occluder LJ parameters.

``conditional_widom_profile``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Mode-A-style decomposition into hard accessibility and conditional charge-off
  insertion over the accessible region.

.. math::

   \mu_{\mathrm{hard}}(z)=-\beta^{-1}\ln h(z)

.. math::

   \mu_{\mathrm{cond}}(z)=-\beta^{-1}\ln\left\langle e^{-\beta \Delta U_0}\right\rangle_{z,\mathrm{accessible}}

.. math::

   \mu_{\mathrm{full}}(z)=\mu_{\mathrm{hard}}(z)+\mu_{\mathrm{cond}}(z)

**Outputs**
  Reservoir-referenced hard, conditional, and full charge-off insertion factors;
  corresponding free-energy profiles; SEM estimates from the per-frame distribution;
  and a conditional/full identity residual for audit.

Anisotropy, geometry, dynamic heterogeneity, and local structure
----------------------------------------------------------------

``src/measures/RegisterAnisotropy.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``directional_msd``
^^^^^^^^^^^^^^^^^^^

.. math::

   MSD_{\alpha\beta}(t)=\left\langle\frac{1}{N}\sum_i
   \Delta r_{i,\alpha}(t;t_0)\Delta r_{i,\beta}(t;t_0)\right\rangle_{t_0}

``directional_isf``
^^^^^^^^^^^^^^^^^^^

.. math::

   F_{s,\alpha}(q,t)=\left\langle\frac{1}{N}\sum_i \cos[q\,\Delta r_{i,\alpha}(t;t_0)]\right\rangle_{t_0}

``slab_rdf``
^^^^^^^^^^^^

**Meaning**
  Slab-resolved radial structure using a global density normalization.

``cylindrical_profile`` / ``2d_map`` / ``interface_excess``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Direction-resolved density / weight maps and Gibbs-style interfacial excess estimates.

``src/measures/RegisterHeterogeneity.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``g4_rt``, ``s4_qt``, ``xi4_t``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Spatial amplitude and length-scale measures of dynamic heterogeneity.

.. math::

   S_4(q,t)\approx \frac{S_4(0,t)}{1+[q\xi_4(t)]^2}

``string_motion`` / ``string_length_dist`` / ``mobile_cluster``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Cooperative motion diagnostics based on mobile-particle connectivity in real space or exchange space.

``excitation_map`` / ``facilitation_acf``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Per-particle excitation propensity and facilitation diagnostics.

``src/measures/RegisterLocalStructure.cpp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``q4_q6_w6``
^^^^^^^^^^^^

**Outputs**
  summary statistics and optional per-particle values for bond-orientational order invariants

.. math::

   q_l(i)=\sqrt{\frac{4\pi}{2l+1}\sum_{m=-l}^l |q_{lm}(i)|^2}

``locally_favored_structure``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning**
  Threshold-based classifier over :math:`q_4`, :math:`q_6`, :math:`\hat w_6`, and coordination number.

``voronoi_volume``
^^^^^^^^^^^^^^^^^^

**Meaning**
  Statistics over an externally supplied Voronoi-volume field.

``free_volume``
^^^^^^^^^^^^^^^

.. math::

   v_i^{\mathrm{free}}=V_i^{\mathrm{Vor}}-v_i^{\mathrm{excl}}

``softness_proxy``
^^^^^^^^^^^^^^^^^^

**Meaning**
  Transparent linear proxy built from local structural descriptors; not a hidden ML model.

Implementation notes for maintainers
------------------------------------

* Each measure family lives in one ``src/measures/*.cpp`` translation unit and registers itself at the bottom.
* If a new measure shares geometry, chain-ordering, entity-grouping, or integer-field logic, extend
  ``src/measures/MeasureCommon.hpp`` rather than duplicating helpers.
* Two-time / aging observables should stay separate from the stationary correlator path, because their
  output is a waiting-time/lag grid rather than a one-dimensional lag-average.
