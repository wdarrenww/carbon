## **Phase 1 – Core Engine & Minimal Physics (2–3 weeks)**

**Goal:** A working simulation loop with atoms you can place manually, bonds defined, and basic forces.

**Deliverables:**

1. **Data structures**

   * `Atom` struct (Z, mass, charge, sigma, epsilon, r, v, f)

   * `Bond` struct (i, j, k, r0)

   * `System` class (vectors of atoms/bonds \+ box info)

   * Basic constants table (H, C, O, N, Na, Cl)

2. **Math utils**

   * Vec3 type (or Eigen)

   * Distance, normalize, periodic wrapping

3. **Physics implementation**

   * Harmonic bonds  
      Ub=k(r−r0)2U\_b \= k(r \- r\_0)^2Ub​=k(r−r0​)2

   * Simple Lennard-Jones nonbonded forces (no neighbor list yet)

   * Vacuum Coulomb forces (k q1 q2 / r²)

   * Velocity-Verlet integration

4. **I/O**

   * Read/write XYZ files

   * Simple text-based “build” mode to place atoms, define bonds manually

5. **Simulation loop**

   * Forces → Integrate → Output positions each step

**End of phase milestone:**  
 You can create a water molecule by typing coordinates, define bonds, run the simulation, and see the atoms vibrate in an XYZ trajectory.

---

## **Phase 2 – Performance & Stability (2–3 weeks)**

**Goal:** Handle dozens–hundreds of atoms without big slowdowns, add stability to integration.

**Deliverables:**

1. **Neighbor list**

   * Verlet list with skin distance, updated every N steps

   * Optional cell list for large systems

2. **Periodic boundary conditions**

   * Wrap positions after movement

   * Apply minimum image convention in distance calc

3. **Basic thermostats**

   * Langevin or Berendsen to control temperature

4. **Energy tracking**

   * Track kinetic, potential, and total energy

   * Output to log file

5. **Improved input**

   * Load prebuilt molecule from XYZ

   * Auto-generate bonds based on covalent radii

**End of phase milestone:**  
 You can simulate a box of \~500 water molecules with reasonable speed and stable temperature, saving both XYZ trajectory and energy log.

---

## **Phase 3 – Chemistry Realism Upgrade (3–4 weeks)**

**Goal:** Add more force field terms and chemistry features.

**Deliverables:**

1. **Angle & dihedral potentials**

   * Harmonic angles

   * Periodic torsions

2. **Force field parameterization**

   * Small JSON file for element parameters (mass, sigma, epsilon, covalent radius)

   * Default bond/angle/dihedral parameters for common molecules

3. **Cutoffs & switching**

   * Lennard-Jones switching function to avoid force discontinuities

   * Coulomb cutoff for MVP (PME later)

4. **Constraints**

   * SHAKE for rigid bonds involving hydrogens → allows 2 fs timestep

5. **Default molecule templates**

   * Water, methane, ethanol, benzene, NaCl lattice

**End of phase milestone:**  
 You can load methane, ethanol, or a benzene dimer into the simulator, run dynamics with bonded \+ nonbonded terms, and get physically sensible vibrations and structure.

---

## **Phase 4 – Interactivity & GUI (4–6 weeks)**

**Goal:** Make it user-friendly to construct and visualize molecules.

**Deliverables:**

1. **Live rendering**

   * Dear ImGui \+ OpenGL viewer

   * Ball-and-stick or space-filling model

   * Color by element

2. **Interactive builder**

   * Place atoms by clicking

   * Drag bonds between atoms

   * Autogenerate angles/dihedrals

3. **Control panel**

   * Start/stop simulation

   * Adjust timestep, temperature

   * Load/save scenes

4. **Highlight & inspect**

   * Display selected atom’s element, position, velocity, partial charge

   * Show bond lengths & angles live

**End of phase milestone:**  
 You can build a molecule in the GUI, press "Run," and watch it move realistically in real time.

---

## **Phase 5 – Educational/Advanced Features (optional, ongoing)**

**Goal:** Make it a learning tool or more realistic simulator.

**Possible deliverables:**

1. **Reactive mode**

   * Simple bond creation/breaking when distances cross thresholds

   * Morse bond potential

2. **QM/MM hybrid mode**

   * Small subset of atoms calculated with a basic tight-binding method

3. **Implicit solvent**

   * Constant dielectric screening for Coulomb

4. **Analysis tools**

   * Radial distribution function (g(r))

   * Mean squared displacement

   * Energy decomposition by force term

5. **GPU acceleration**

   * OpenCL/CUDA for nonbonded loops

---

## **Summary Timeline**

| Phase | Duration | Atoms Target | Key Outcome |
| ----- | ----- | ----- | ----- |
| 1 – Core Engine | 2–3 weeks | \<50 | First molecules moving |
| 2 – Performance/Stability | 2–3 weeks | \~500 | Stable, fast, with thermostat |
| 3 – Chemistry Realism | 3–4 weeks | \~1000 | Bonded \+ nonbonded \+ angles/dihedrals |
| 4 – GUI & Interactivity | 4–6 weeks | \~200 live | Real-time construction & viewing |
| 5 – Advanced/Edu | Ongoing | Flexible | Advanced physics & analysis |
