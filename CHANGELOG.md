# Changelog

Toutes les modifications notables de ce projet seront documentées dans ce fichier.

Le format est basé sur [Keep a Changelog](https://keepachangelog.com/fr/1.0.0/),
et ce projet adhère au versionnage sémantique.

## [0.2.0] - 2026-02-23

### Ajouté
- Implémentation de conditions aux limites **PML (Perfectly Matched Layers)** de haute précision avec profils d'amortissement directionnels indépendants (`sigma_x` et `sigma_y`).
- Mise en cache de la vitesse maximale du son (`cached_c_max`) pour éviter des recalculs lourds et inutiles à chaque pas de temps ou test de stabilité.

### Modifié
- Remplacement de l'ancienne couche éponge basique (`set_absorbing_boundary`) par la méthode scientifique `set_pml_boundary`.
- Mise à jour de la boucle de calcul `step` pour appliquer l'atténuation PML semi-implicite sur les vitesses et la pression.
- Mise à jour de la documentation (`README.md`) pour refléter les nouveautés de la version 0.2.

## [0.1.0] - Première version stable
- Solveur acoustique 2D FDTD sur grille décalée en `f32`.
- Gestion des milieux hétérogènes (`rho`, `c` par cellule).
- Sources (Ricker, sinus, impulsion) et récepteurs ponctuels.
- Parallélisation optionnelle avec `rayon` et rendu temps réel avec `minifb`.