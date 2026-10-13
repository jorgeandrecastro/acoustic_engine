# acoustic-engine

Solveur acoustique 2D par différences finies (FDTD, grille décalée) en `f32`.
Licence : **GPL-2.0-or-later**.

## Fonctionnalités

- Propagation d'ondes dans un milieu hétérogène (`rho`, `c` par cellule)
- Conditions aux limites ultra-absorbantes **PML (Perfectly Matched Layers)** ou murs rigides
- Sources temporelles : ondelette de Ricker, sinus, impulsion
- Récepteurs : enregistrement de `p(t)` en un point
- Pas de temps stable calculé (condition CFL)
- Parallélisation optionnelle avec la feature `rayon`
- Aucune dépendance par défaut

## Exemple

```rust
use acoustic_engine::{AcousticSolver2D, Waveform};

let mut s = AcousticSolver2D::new(300, 300, 0.5, 0.5, 1.2, 340.0); // air
s.set_pml_boundary(25, 15.0);
s.set_medium_rect(200, 100, 210, 200, 1200.0, 340.0); // mur dense
s.add_source(100, 150, Waveform::Ricker { freq: 100.0 }, 1.0);
let mic = s.add_receiver(250, 150);

let dt = s.stable_dt(0.9);
for _ in 0..1000 { s.step(dt); }
let signal: &[f32] = s.receiver(mic);
```

## Exemples fournis

```
cargo run --release --example impulse
cargo run --release --example receiver                 # écrit receiver.csv
cargo run --release --example bench --features rayon
cargo run --release --example viewer --features rayon  # fenêtre temps réel
```

## Règles de précision

- Garder au moins ~10 cellules par longueur d'onde minimale (`dx <= c / (10·f_max)`),
  sinon la dispersion numérique déforme l'onde.
- Toujours utiliser `dt <= max_stable_dt()` (`stable_dt(0.9)` recommandé).
- Pour modifier le milieu, utiliser `set_medium` / `set_medium_rect` (ils tiennent
  `1/rho` à jour).

## Tests

```
cargo test
cargo test --features rayon
```

Les tests vérifient la formule CFL, la symétrie de l'onde, la causalité, la stabilité
sous/sur la limite CFL,l'efficacité des limites PML et l'effet d'un milieu hétérogène.


## Changelog

- v0.2.0 : Intégration des conditions aux limites PML (Perfectly Matched Layers) directionnelles et mise en cache de la vitesse maximale (cached_c_max) pour optimiser les performances.

- v0.1.0 : Première version stable (solveur FDTD 2D, gestion des milieux hétérogènes, sources/récepteurs, support de Rayon et visualisation minifb).

## Licence

GPL-2.0-or-later. Placer le texte officiel dans `LICENSE`
(https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt).
